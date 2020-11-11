# This code allows the user to remove a portion of the network by cutting a
# slice of the data.

# 0. Load the network
# 1. Define the cut plane
# 2. Find if any of the neurons are above the cut plane, remove those
# 3. Remova any synapses above the cut plane
# 4. Write a new hdf5
#

import os
import h5py
import numpy as np
from .load import SnuddaLoad
import time


class SnuddaCut(object):
  
  def __init__(self,networkFile,cutEquation="z>0",
               outFileName=None,
               plotOnly=False,showPlot=True):

    #self.snuddaLoader = SnuddaLoad(networkFile)
    #self.data = self.snuddaLoader.data

    self.cutEquation = cutEquation

    self.h5libver = "latest"
    self.h5driver = "sec2"
    self.inFile = None
    self.outFile = None
    
    if(outFileName is None):
      basePath = os.path.dirname(networkFile)
      self.outFileName = basePath + "/network-cut-slice.hdf5"
    else:
      self.outFileName = outFileName
      
      
    # We create a lambda expression, to avoid calling eval each time
    self.cutEquationLambda = eval("lambda x,y,z :" + cutEquation)

    self.openInputFile(networkFile)

    if(plotOnly):
      self.plotCut(includeSynapses=False,includeGapJunctions=False)
      exit(0)
    
    self.setupOutputFile(self.outFileName)
    self.writeCutSlice(self.cutEquationLambda)

    self.plotCut(includeSynapses=True,includeGapJunctions=True,
                 showPlot=showPlot)
    
    if(False):
      self.inFile.close()
      self.out_file.close()

    
    
  ############################################################################

  def writeCutSlice(self,cutEquationLambda):

    # Remove the neurons from the data    
    somaKeepFlag = self.somasInside(cutEquationLambda)
    somaKeepID = np.where(somaKeepFlag)[0]
    somaRemoveID = np.where(somaKeepFlag == False)[0]
    nSomaKeep = np.sum(somaKeepFlag)

    if(nSomaKeep == 0):
      print("No somas left, aborting!")
      exit(-1)
    
    print("Keeping " + str(nSomaKeep) + " out of " + str(len(somaKeepFlag)) \
          + " neurons (the others have soma outside of cut plane)")
    
    # We need to remap neuronID in the synapses and gap junction matrix
    remapID = dict([])
    for newID,oldID in enumerate(somaKeepID):
      remapID[oldID] = newID
    
    networkGroup = self.outFile.create_group("network")
    neuronGroup = networkGroup.create_group("neurons")

    for varName in self.inFile["network/neurons"]:

      data = self.inFile["network/neurons/" + varName]
      
      if(len(data.shape) == 0):
        # Scalar data, just copy
        self.inFile.copy("network/neurons/" +varName,neuronGroup)
        continue
        
      elif(len(data.shape) == 1):
        # 1D data, we only keep nSomaKeep of them
        dataShape = (nSomaKeep,)
      elif(len(data.shape) == 2):
        # 2D data, need to make sure to maintain dimensions
        dataShape = (nSomaKeep,data.shape[1])
      else:
        print("writeCutSlice: Only handle 0D, 1D and 2D data, update code!")
        import pdb
        pdb.set_trace()

      if(varName == "neuronID"):
        # We need to remap
        neuronGroup.create_dataset(varName, \
                                   dataShape,
                                   data.dtype,
                                   [remapID[data[x]] for x in somaKeepID],
                                   compression=data.compression)

        # Double check that it is OK, should be in order after
        assert (np.diff(neuronGroup["neuronID"][()]) == 1).all(), \
          "Problem with neuron remapping!"
        
      else:
        try:
          neuronGroup.create_dataset(varName, \
                                     dataShape,
                                     data.dtype,
                                     [data[x] for x in somaKeepID],
                                     compression=data.compression)
        except:
          import traceback
          tstr = traceback.format_exc()
          print(tstr)
          import pdb
          pdb.set_trace()
          

    # Next deal with synapses
    keepSynFlag = self.synapsesInside(cutEquationLambda,dataType="synapses")
    keepSynFlag = self.filterNeuronsSynapses(somaRemoveID,
                                             keepFlag=keepSynFlag,
                                             dataType="synapses")
    
    # Lastly deal with gap junctions
    keepGJFlag = self.synapsesInside(cutEquationLambda,dataType="gapJunctions")
    keepGJFlag = self.filterNeuronsSynapses(somaRemoveID,
                                            keepFlag=keepGJFlag,
                                            dataType="gapJunctions")

    nSyn = np.sum(keepSynFlag)
    nSynapses = np.zeros((1,),dtype=np.uint64) + nSyn
    
    nGJ = np.sum(keepGJFlag)
    nGapJunctions = np.zeros((1,),dtype=np.uint64) + nGJ
    
    networkGroup.create_dataset("nSynapses",data=nSynapses,dtype=np.uint64)
    networkGroup.create_dataset("nGapJunctions",data=nGapJunctions,
                                dtype=np.uint64)

    # !!!! TO BE CONTINUED... might need to handle cunk size differently based
    # on size... see snudda_prune line 1000

    synMat = self.inFile["network/synapses"]
    GJMat = self.inFile["network/gapJunctions"]

    print("Copying synapses and gap junctions")
    
    networkGroup.create_dataset("synapses", \
                                dtype=np.int32, \
                                shape = (nSyn,13), \
                                chunks = synMat.chunks, \
                                maxshape=(None,13), \
                                compression=synMat.compression)   

    networkGroup.create_dataset("gapJunctions", \
                                dtype=np.int32, \
                                shape = (nGJ,11), \
                                chunks = GJMat.chunks, \
                                maxshape=(None,11), \
                                compression=GJMat.compression)   

    for idx,rowIdx in enumerate(np.where(keepSynFlag)[0]):
      # We need to remap the neuronID if some neurons have been removed!!
      row = synMat[rowIdx,:]
      row[0] = remapID[row[0]]
      row[1] = remapID[row[1]]
      networkGroup["synapses"][idx,:] = row
        
    print("Keeping " + str(nSyn) + " synapses (out of " \
          + str(synMat.shape[0]) + ")")

    for idx,rowIdx in enumerate(np.where(keepGJFlag)[0]):
      # We need to remap the neuronID if some neurons have been removed!!
      row = GJMat[rowIdx,:]
      row[0] = remapID[row[0]]
      row[1] = remapID[row[1]]      
      networkGroup["gapJunctions"][idx,:] = row

    print("Keeping " + str(nGJ) + " gap junctions (out of " \
          + str(GJMat.shape[0]) + ")")

  
  ############################################################################
  
  # Tells which somas are inside the cut equation
  
  def somasInside(self,cutEquationLambda):

    pos = self.inFile["network/neurons/position"][()]
    insideFlag = np.array([cutEquationLambda(x,y,z) for x,y,z in pos], \
                          dtype=bool)

    return insideFlag

  ############################################################################
  
  def synapsesInside(self,cutEquationLambda,dataType="synapses"):

    voxelSize = self.inFile["meta/voxelSize"][()]
    simOrigo = self.inFile["meta/simulationOrigo"][()]

    if(dataType == "synapses"):
      pos = self.inFile["network/synapses"][:,2:5] * voxelSize + simOrigo
    elif(dataType == "gapJunctions"):
      pos = self.inFile["network/gapJunctions"][:,6:9] * voxelSize + simOrigo
    else:
      print("filterNeuronsSynapses: Unknown dataType: " + str(dataType))
      import pdb
      pdb.set_trace()
      
    insideFlag = np.array([cutEquationLambda(x,y,z) for x,y,z in pos], \
                          dtype=bool)

    return insideFlag

  ############################################################################

  # Returns the row numbers that do not contain the neuronID, ie filters
  # the synapses belonging to neuronID out...

  # dataType = "synapses" or "gapJunctions"
  
  def filterNeuronsSynapses(self,neuronID,keepFlag=None,dataType="synapses"):
  
    if(dataType == "synapses"):
      dataStr = "network/synapses"
    elif(dataType == "gapJunctions"):
      dataStr = "network/gapJunctions"
    else:
      print("filterNeuronsSynapses: Unknown dataType: " + str(dataType))
      import pdb
      pdb.set_trace()

    if(keepFlag is None):                          
      keepFlag = np.ones((self.inFile[dataStr].shape[0],),dtype=bool)

    srcID = self.inFile[dataStr][:,0]
    destID = self.inFile[dataStr][:,1]    
    
    for nID in neuronID:
      
      keepFlag = np.logical_and(keepFlag,
                                np.logical_and(srcID != nID, destID != nID))

    return keepFlag

  ############################################################################

  def openInputFile(self,networkFile):

    self.inFile = h5py.File(networkFile,"r",
                            libver=self.h5libver,
                            driver=self.h5driver)
    

  ############################################################################

  # This sets up the output file, copies the config and meta data,
  # but does not copy over the neurons, synapses or gap junctions
  
  def setupOutputFile(self,outFileName):

    print("Writing to " + str(outFileName))
    
    self.outFile = h5py.File(outFileName, "w",
                             libver=self.h5libver,
                             driver=self.h5driver)

    print("Copying 'config' and 'meta'")
    self.inFile.copy("config",self.outFile)
    self.inFile.copy("meta",self.outFile)

    if("morphologies" in self.inFile):
      print("Copying morphologies")
      self.inFile.copy("morphologies",self.outFile)


  ############################################################################

  # This is just used to verify
  
  def plotCut(self,includeSynapses=True,includeGapJunctions=True,showPlot=True):

    print("Plotting verification figure")
    
    voxelSize = self.inFile["meta/voxelSize"][()]
    simOrigo = self.inFile["meta/simulationOrigo"][()]

    
    inPos = self.inFile["network/neurons/position"][()]
    inSyn = self.inFile["network/synapses"][:,2:5]*voxelSize + simOrigo
    inGJ = self.inFile["network/gapJunctions"][:,6:9]*voxelSize + simOrigo

    if(self.outFile is not None):

      # Just double check that they match
      assert self.outFile["meta/voxelSize"][()] == voxelSize
      assert (self.outFile["meta/simulationOrigo"][()] == simOrigo).all()

      outPos = self.outFile["network/neurons/position"][()]    
      outSyn = self.outFile["network/synapses"][:,2:5]*voxelSize + simOrigo    
      outGJ = self.outFile["network/gapJunctions"][:,6:9]*voxelSize + simOrigo

    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    aVal=0.5
    aVal2=0.2    
    
    ax.scatter(inPos[:,0],inPos[:,1],inPos[:,2],c='black',s=25,alpha=aVal2)
    if(self.outFile is not None):
      ax.scatter(outPos[:,0],outPos[:,1],outPos[:,2],c='red',s=15,alpha=aVal) 

    if(includeSynapses):
      ax.scatter(inSyn[:,0],inSyn[:,1],inSyn[:,2],c='black',s=9,alpha=aVal2)

      if(self.outFile is not None):
        ax.scatter(outSyn[:,0],outSyn[:,1],outSyn[:,2],c='red',s=4,alpha=aVal)

    if(includeGapJunctions):
      ax.scatter(inGJ[:,0],inGJ[:,1],inGJ[:,2],c='blue',s=6,alpha=aVal2)
      if(self.outFile is not None):
        ax.scatter(outGJ[:,0],outGJ[:,1],outGJ[:,2],c='green',s=3,alpha=aVal)

    ax.view_init(elev=-3, azim=-95)

    if(showPlot):
      plt.ion()
      plt.show()

    figName = self.outFileName + ".pdf"
    print("Writing to figure " + figName)
    plt.savefig(figName,dpi=300)

    if(showPlot):
      print("Inspect plot, then quit debug")
      print("The viewing angle might not be good for your try, so leave it interactive")
      import pdb
      pdb.set_trace()
    else:
      print("Waiting 5 seconds")
      time.sleep(5)
      
  ############################################################################
                             
if __name__ == "__main__":

  from argparse import ArgumentParser

  parser = ArgumentParser(description="Cut a slice from a network")
  parser.add_argument("networkFile",help="Network file (hdf5)",type=str)
  parser.add_argument("cutEquation",
                      help="Equation defining parts left after cut, " \
                          + "e.g. 'z>0' or 'x+y>100e-6' (SI units)",
                      type=str)
  parser.add_argument("--plotOnly", \
                      help="Plots the network without cutting",
                      action="store_true")
  parser.add_argument("--hidePlot", help="Hide plot.", action="store_true")
  
  args = parser.parse_args()

  if(args.hidePlot):
    showPlot = False
  else:
    showPlot = True
    
  sc = SnuddaCut(networkFile=args.networkFile,
                 cutEquation=args.cutEquation,
                 plotOnly=args.plotOnly,
                 showPlot=showPlot)
  
