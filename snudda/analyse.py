# This script is custom written to handle very large datasets. It does so by
# not keeping all the information in memory, instead parsing the HDF5
# piece by piece

import numpy as np
import scipy.sparse as sps
import h5py
import timeit
import time
from glob import glob
import os
import sys
import json
import pickle

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D



from .load import SnuddaLoad

# !!! We need to parallelise the analysis script also!

class SnuddaAnalyse(object):

  # saveCache = should we save a pickled file with connection matrix?
  # loadCache = should we load cache file if available
  # lowMemory = if false, uses dense matrix which is faster (assuming lots of memory)

  def __init__(self,
               hdf5File=None,
               loadCache=False,
               saveCache=True,
               lowMemory=False,
               sideLen=250e-6,
               volumeType="cube",
               nMaxAnalyse=None,
               showPlots=False,
               closePlots=True): # "cube" or "full"

    self.debug = False
    self.showPlots = showPlots

    print("Assuming volume type: " + str(volumeType) \
          + "[cube or full]")

    self.volumeType = volumeType
    self.closePlots = closePlots

    if(nMaxAnalyse is None):
      if(volumeType == "cube"):
        nMaxAnalyse = 20000
      elif(volumeType == "full"):
        nMaxAnalyse = 20000

    print("Only using " + str(nMaxAnalyse) + "neurons of the connection data")

    self.nMaxAnalyse = nMaxAnalyse

    if(hdf5File is None or hdf5File == "last"):
      hdf5File = self.findLatestFile()

    self.figDir = os.path.dirname(hdf5File) + "/figures"
    if(not os.path.exists(self.figDir)):
      os.makedirs(self.figDir)


    # First load all data but synapses
    self.networkLoad = SnuddaLoad(hdf5File,loadSynapses=False)
    self.network = self.networkLoad.data

    if("config" in self.network):
      self.config = json.loads(self.network["config"])
    self.sideLen = sideLen

    self.lowMemory = lowMemory

    self.neuronNameRemap = {"FSN" : "FS"}

    cacheLoaded = False
    if(loadCache):
      try:
        cacheLoaded = self.loadCacheData(hdf5File)
      except:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)
        assert not cacheLoaded, "Load failed, cacheLoaded flag should be False"

    self.data = h5py.File(hdf5File,'r')

    if(not cacheLoaded):
      self.nNeurons = self.network["nNeurons"]
      print("Number of neurons: " + str(self.nNeurons))

      # GABA connection matrix (synType = 1) (ignore AMPA/NMDA = 2, GJ = 3)

      #self.connectionMatrix = self.createConnectionMatrix(synType=1,
      #                                                    lowMemory=lowMemory)

      self.connectionMatrix = self.createConnectionMatrix(lowMemory=lowMemory)
      self.connectionMatrixGJ = self.createConnectionMatrixGJ()

      #self.connectionMatrix = self.createConnectionMatrixSLOW(synType=1)
      self.makePopDict()
      self.positions = self.network["neuronPositions"]

      self.synapseDist()


      if(saveCache):
        self.saveCacheData(hdf5File)

    self.workerData = []

    self.neuronColors = {"dSPN" : (77./255,151./255,1.0),
                         "iSPN" : (67./255,55./255,181./255),
                         "FSN"  : (6./255,31./255,85./255),
                         "ChIN" : (252./266,102./255,0.0),
                         "LTS"  : (150./255,63./255,212./255),
                         "default" : [0.4, 0.4, 0.4] }

  ############################################################################

  def neuronName(self,neuronType):

    if(neuronType in self.neuronNameRemap):
      return self.neuronNameRemap[neuronType]
    else:
      return neuronType

  ############################################################################

  def getNeuronColor(self,neuronType):

    if(neuronType in self.neuronColors):
      return self.neuronColors[neuronType]
    else:
      return self.neuronColors["default"]

  ############################################################################

  # Reading the HDF5 files takes a lot of time, this stores a cached copy
  # of the connection matrix, positions and populations

  def saveCacheData(self,hdf5File):

    import h5py

    cacheFile = hdf5File + "-cache"

    print("Saving cache to " + cacheFile    )
    outFile = h5py.File(cacheFile,'w',libver='latest')

    # Connection Matrix
    outFile.create_dataset("conMat_data",data=self.connectionMatrix.data, \
                           compression='gzip')
    outFile.create_dataset("conMat_indices", \
                           data=self.connectionMatrix.indices, \
                           compression='gzip')
    outFile.create_dataset("conMat_indptr", \
                           data=self.connectionMatrix.indptr, \
                           compression='gzip')
    outFile.create_dataset("conMat_shape", \
                           data=self.connectionMatrix.shape)

    # GJ connection matrix
    outFile.create_dataset("conMatGJ_data",data=self.connectionMatrixGJ.data, \
                           compression='gzip')
    outFile.create_dataset("conMatGJ_indices", \
                           data=self.connectionMatrixGJ.indices, \
                           compression='gzip')
    outFile.create_dataset("conMatGJ_indptr", \
                           data=self.connectionMatrixGJ.indptr, \
                           compression='gzip')
    outFile.create_dataset("conMatGJ_shape", \
                           data=self.connectionMatrixGJ.shape)


    popGroup = outFile.create_group("populations")
    for k in self.populations:
      v = self.populations[k]
      popGroup.create_dataset(k,data=v)

    outFile["nNeurons"] = self.nNeurons
    outFile.create_dataset("positions",data=self.positions)

    try:

      dendPosBin = dict([])
      for prePost in self.dendPositionBin:
        pp = self.allTypes[prePost[0]] + "_" + self.allTypes[prePost[1]]
        dendPosBin[pp] = list(self.dendPositionBin[prePost])

      outFile.create_dataset("dendPositionBin", \
                             data=json.dumps(dendPosBin))
      outFile.create_dataset("dendPositionEdges", \
                             data=self.dendPositionEdges)

      allTypes = [x.encode("ascii","ignore") for x in self.allTypes]
      outFile.create_dataset("allTypes",
                             data=allTypes)
      outFile.create_dataset("neuronTypeID",data=self.neuronTypeID)

      outFile.create_dataset("nMaxAnalyse",data=self.nMaxAnalyse)

    except Exception as e:

      import traceback
      tstr = traceback.format_exc()
      print(tstr)

      print("Problem with writing hdf5")
      import pdb
      pdb.set_trace()

    outFile.close()

  ############################################################################

  def loadCacheData(self,hdf5File):

    import os
    import h5py

    cacheFile = hdf5File + "-cache"
    dataLoaded = False

    if(os.path.exists(cacheFile)):
      tOrig = os.path.getmtime(hdf5File)
      tCache = os.path.getmtime(cacheFile)

      # Make sure cache file is newer than data file
      if(tCache > tOrig):
        print("Loading from " + cacheFile)

        try:
          with h5py.File(cacheFile,'r') as data:

            assert self.nMaxAnalyse == data["nMaxAnalyse"][()], \
              "nMaxAnalyse has changed, have to reload connection matrix"

            self.connectionMatrix = sps.csr_matrix((data["conMat_data"], \
                                                    data["conMat_indices"], \
                                                    data["conMat_indptr"]),
                                                    data["conMat_shape"])

            self.connectionMatrixGJ = sps.csr_matrix((data["conMatGJ_data"], \
                                                    data["conMatGJ_indices"], \
                                                    data["conMatGJ_indptr"]),
                                                    data["conMatGJ_shape"])

            self.populations = dict([])

            for k in data["populations"].keys():
              self.populations[k] = data["populations"][k][:]

              self.nNeurons = data["nNeurons"][()]
              self.positions = data["positions"][:]

            dendPosBin = json.loads(data["dendPositionBin"][()])
            self.dendPositionBin = dict([])

            allTypes = list(data["allTypes"][()])
            self.allTypes = [x.decode() for x in allTypes]
            self.neuronTypeID = data["neuronTypeID"][()]

            for pp in dendPosBin:
              str = pp.split("_")
              preType = self.allTypes.index(str[0])
              postType = self.allTypes.index(str[1])

              self.dendPositionBin[(preType,postType)] \
                = np.array(dendPosBin[pp])

            self.dendPositionEdges = data["dendPositionEdges"][()]

            #import pdb
            #pdb.set_trace()

          dataLoaded = True
          print("Loading done.")

        except Exception as e:
          import traceback
          tstr = traceback.format_exc()
          print(tstr)
          print("Failed to load cache file.")
          dataLoaded = False
          #import pdb
          #pdb.set_trace()

    return dataLoaded

  ############################################################################

  def createConnectionMatrix(self,chunkSize=1000000, \
                             synType=None,
                             minDendDist=None,
                             maxDendDist=None,
                             lowMemory=False):

    t0 = timeit.default_timer()

    if(lowMemory):
      print("Trying to conserve memory, this is slower.")
      connectionMatrix = sps.lil_matrix((self.nNeurons,self.nNeurons),\
                                        dtype=np.int16)
    else:
      try:
        connectionMatrix = np.zeros((self.nNeurons,self.nNeurons),\
                                    dtype=np.int16)
      except:
        print("Unable to allocate full matrix, using sparse matrix instead")
        connectionMatrix = sps.lil_matrix((self.nNeurons,self.nNeurons),\
                                          dtype=np.int16)
    lastSrcID = 0
    lastDestID = 0
    lastCount = 0

    rowCtr = 0
    nSynTotal = self.network["nSynapses"]

    for synapses in self.networkLoad.synapseIterator(chunkSize=chunkSize):

      print("Synapse row " + str(rowCtr) \
            + " - " + str(100*rowCtr/float(nSynTotal)) + " %" \
            + " time: " + str(timeit.default_timer()-t0) + " seconds")

      for synRow in synapses:

        rowCtr += 1

        srcID = synRow[0]
        destID = synRow[1] # New format!
        dendDist = synRow[8]

        if((minDendDist is not None and dendDist < minDendDist)
           or (maxDendDist is not None and dendDist > maxDendDist)):
           # Not correct distance to soma on dendrite
           continue

        if(synType is None or synType == synRow[6]):
          # Only include specific synapse type

          if(lastSrcID == srcID and lastDestID == destID):
            lastCount += 1
          else:
            # Write the previous set of synapses to matrix
            # For first iteration of loop, lastCount is zero
            connectionMatrix[lastSrcID,lastDestID] += lastCount

            lastSrcID = srcID
            lastDestID = destID
            lastCount = 1

          # connectionMatrix[srcID,destID] += 1

    # Update the last row also
    connectionMatrix[lastSrcID,lastDestID] += lastCount
    lastCount = 0

    t1 = timeit.default_timer()

    print("Created connection matrix " + str(t1-t0) + " seconds")

    return sps.csr_matrix(connectionMatrix,dtype=np.int16)


  ############################################################################

  def createConnectionMatrixGJ(self):

    t0 = timeit.default_timer()

    connectionMatrixGJ = sps.lil_matrix((self.nNeurons,self.nNeurons),\
                                        dtype=np.int16)

    lastSrcID = 0
    lastDestID = 0
    lastCount = 0

    rowCtr = 0
    nGJTotal = self.network["nGapJunctions"]

    for gjList in self.networkLoad.gapJunctionIterator(chunkSize=100000):
      print("GJ row : " + str(rowCtr) \
            + " - " + str(100*rowCtr/float(nGJTotal)) + " % " \
            + " time : " + str(timeit.default_timer()-t0) + " seconds")

      for gjRow in gjList:
        rowCtr += 1

        srcID = gjRow[0]
        destID = gjRow[1]

        if(lastSrcID == srcID and lastDestID == destID):
          lastCount += 1
        else:
          connectionMatrixGJ[lastSrcID,lastDestID] += lastCount

          lastSrcID = srcID
          lastDestID = destID
          lastCount = 1

      connectionMatrixGJ[lastSrcID,lastDestID] += lastCount
      lastCount = 0

    t1 = timeit.default_timer()
    print("Created gap junction connection matrix " + str(t1-t0) + " seconds")

    return sps.csr_matrix(connectionMatrixGJ,dtype=np.int16)

  ############################################################################

  def makePopDict(self):

    print("Creating population dictionary")

    self.populations = dict([])

    for nid,neuron in enumerate(self.network["neurons"]):

      assert(nid == neuron["neuronID"])
      name = neuron["name"].split("_")[0]

      if(name not in self.populations):
        self.populations[name] = []

      self.populations[name].append(neuron["neuronID"])

    print("Done.")

  ############################################################################

  def getSubPop(self,volumeType="cube",volumePart="centre",sideLen=None,
                neuronID=None,volumeID="Striatum",nMaxAnalyse=None):

    # print("volumeType=" + volumeType + ",volumePart=" + volumePart + ",sideLen=" +str(sideLen))

    if(volumeType == "full"):
      # return all neurons

      if(volumeID is not None):
        idx = np.where([x["volumeID"] == volumeID \
                        for x in self.network["neurons"]])[0]

      if(neuronID is None):
        neuronID = idx

    elif(volumeType == "cube"):

      if(volumePart == "centre"):
        try:
          neuronID = self.centreNeurons(sideLen=sideLen,
                                        neuronID=neuronID,
                                        volumeID=volumeID)
        except:
          import traceback
          tstr = traceback.format_exc()
          print(tstr)
          import pdb
          pdb.set_trace()

      elif(volumePart == "corner"):
        neuronID = self.cornerNeurons(sideLen=sideLen,
                                      neuronID=neuronID,
                                      volumeID=volumeID)
    else:

      print("Unknown volume type: " + str(volumeType))
      import pdb
      pdb.set_trace()

    if(nMaxAnalyse is None):
      nMaxAnalyse = self.nMaxAnalyse

    if(nMaxAnalyse is not None):

      if(nMaxAnalyse < len(neuronID)):

        try:
          keepIdx = np.linspace(0,len(neuronID),nMaxAnalyse,
                                endpoint=False,dtype=int)
          print("Returning subset of neurons to analyse:" \
                + str(len(keepIdx)) + "/" + str(len(neuronID)))
          neuronID = np.array([neuronID[x] for x in keepIdx])
        except:
          import traceback
          tstr = traceback.format_exc()
          print(tstr)

          print("no no no wrong")
          import pdb
          pdb.set_trace()


    return neuronID

  ############################################################################

  def centreNeurons(self,sideLen=None,neuronID=None,volumeID="Striatum"):

    if(sideLen is None):
      sideLen = self.sideLen

    if(volumeID is None):
      idx = np.arange(0,self.network["nNeurons"])
    else:
      idx = np.where([x["volumeID"] == volumeID \
                      for x in self.network["neurons"]])[0]

    try:
      minCoord = np.min(self.network["neuronPositions"][idx,:],axis=0)
      maxCoord = np.max(self.network["neuronPositions"][idx,:],axis=0)
    except:
      import traceback
      tstr = traceback.format_exc()
      print(tstr)
      import pdb
      pdb.set_trace()


    xMin = minCoord[0]
    yMin = minCoord[1]
    zMin = minCoord[2]

    xMax = maxCoord[0]
    yMax = maxCoord[1]
    zMax = maxCoord[2]

    xCentre = (xMax + xMin)/2.0
    yCentre = (yMax + yMin)/2.0
    zCentre = (zMax + zMin)/2.0

    if(neuronID is None):
      neuronID = idx

    if(sideLen < 0):
      print("We want to have all but the outher layer, assume a cube")
      bufLen = -sideLen
      sideLen = np.min([xMax-xCentre-bufLen,
                        yMax-yCentre-bufLen,
                        zMax-zCentre-bufLen])

      assert sideLen > 0, "Unable to autodetect a good side len"

    if(sideLen is None):
      return neuronID

    cID = []

    for nid in neuronID:
      # pos = self.network["neurons"][nid]["position"]
      pos = self.positions[nid,:]

      assert volumeID is None \
        or self.network["neurons"][nid]["volumeID"] == volumeID, \
        "Neuron " + str(nid) + " does not belong to volumeID " + str(volumeID)

      if(abs(pos[0]-xCentre) <= sideLen
         and abs(pos[1]-yCentre) <= sideLen
         and abs(pos[2]-zCentre) <= sideLen):
        cID.append(nid)

    print("Centering in " + str(volumeID) + " : Keeping " + str(len(cID)) + "/" + str(len(neuronID)))

    return cID

  ############################################################################

  # If we use a corner, and account for missing 1/8th of the surrounding
  # synapses then we can get larger distances.
  #
  # <--->

  def cornerNeurons(self,sideLen=None,neuronID=None, volumeID="Striatum"):

    if(sideLen is None):
      sideLen = self.sideLen

    if(volumeID is None):
      idx = np.arange(0,self.network["nNeurons"])
    else:
      idx = np.where([x["volumeID"] == volumeID \
                      for x in self.network["neurons"]])[0]

    if(len(idx) == 0):
      print("No neurons found in volume " + str(volumeID))

      import pdb
      pdb.set_trace()

    minCoord = np.min(self.network["neuronPositions"][idx,:],axis=0)
    maxCoord = np.max(self.network["neuronPositions"][idx,:],axis=0)

    xMin = minCoord[0]
    yMin = minCoord[1]
    zMin = minCoord[2]

    xMax = maxCoord[0]
    yMax = maxCoord[1]
    zMax = maxCoord[2]

    if(sideLen > xMax - xMin \
      or sideLen > yMax - yMin \
       or sideLen > zMax - zMin):
      print("Warning: the analysis cube specified by sideLen is too large.")
      print("!!! Setting sideLen to None")

      sideLen = None

    if(neuronID is None):
      neuronID = idx

    if(sideLen is None):
      return neuronID

    cID = []

    for nid in neuronID:
      # pos = self.network["neurons"][nid]["position"]
      pos = self.positions[nid,:]

      # We assume centre is at zero
      if((pos[0] <= xMin + sideLen or pos[0] >= xMax - sideLen)
         and (pos[1] <= yMin + sideLen or pos[1] >= yMax - sideLen)
         and (pos[2] <= zMin + sideLen or pos[2] >= zMax - sideLen)):
        cID.append(nid)

    print("Taking corner neurons: Keeping " + str(len(cID)) + "/" + str(len(neuronID)))


    if(False):

      # Debug plot
      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')

      pos = self.data["network"]["neurons"]["position"][()]

      ax.scatter(pos[:,0],pos[:,1],pos[:,2],'black')
      ax.scatter(pos[cID,0],pos[cID,1],pos[cID,2],'red')

      plt.ion()
      #plt.draw()

      if(self.showPlots):
        plt.show()

      plt.pause(0.001)
      import pdb
      pdb.set_trace()

    return cID

  ############################################################################

  def saveFigure(self, plt, figName, figType="pdf"):

    if(not os.path.isdir(self.figDir)):
      os.mkdir(self.figDir)

    fullFigName = self.figDir + "/" + figName + "." + figType

    # Remove part of the frame
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)

    plt.tight_layout()
    plt.pause(0.001)
    plt.savefig(fullFigName)
    plt.savefig(fullFigName.replace('.pdf','.eps'))

    print("Wrote " + fullFigName)

    if(self.closePlots):
      time.sleep(1)
      plt.close()


  ############################################################################

  def plotNumSynapsesPerPair(self,preType,postType,sideLen=None, \
                             nameStr="",volumeID="Striatum",
                             connectionType="synapses"):

    if(sideLen is None):
      sideLen = self.sideLen

    print("Plotting number of connections")

    if(preType not in self.populations):
      print("plotNumSynapsesPerPair: " + str(preType) \
            + " is not in the simulation")
      return

    if(postType not in self.populations):
      print("plotNumSynapsesPerPair: " + str(postType) \
            + " is not in the simulation")
      return

    prePop = self.populations[preType]
    postPop = self.populations[postType]

    if(connectionType == "synapses"):
      conMat = self.connectionMatrix
    elif(connectionType == "gapjunctions"):
      conMat = self.connectionMatrixGJ
    else:
      print("Unknown connectionType: " + str(connectionType))
      print("Please use 'synapses' or 'gapjunctions'")
      import pdb
      pdb.set_trace()


    if(sideLen is not None):
      # We are only looking at post synaptic neurons at the centre,
      # to avoid edge effects
      print("Only analysing centre post synaptic neurons, sideLen = " \
            + str(sideLen) )
      # postPop = self.centreNeurons(neuronID=postPop,sideLen=sideLen)
      postPop = self.getSubPop(volumeType=self.volumeType,
                               volumePart="centre",
                               sideLen=sideLen,
                               neuronID=postPop,
                               volumeID=volumeID)

    print("Calculating max synapses")
    maxSynapses = conMat[prePop,:][:,postPop].max()

    print("Calculating mean synapses")

    # The prune tuning func might set data to 0, we want to exclude those
    meanSynapses=float(np.sum(conMat[prePop,:][:,postPop].data))\
                   / np.sum(conMat[prePop,:][:,postPop].data!=0)


    con = conMat[prePop,:][:,postPop]

    #con = con.toarray()
    #existingCon = con[con != 0]

    existingCon = con[np.nonzero(con)].transpose()

    # Any connections? Otherwise skip plot
    if((type(existingCon) == np.matrixlib.defmatrix.matrix \
        and len(existingCon) == 0)
       or (type(existingCon) != np.matrixlib.defmatrix.matrix \
           and existingCon.getnnz() == 0)):
      return

    print("Plotting " + str(existingCon.shape[0]) + " connections")

    #import pdb
    #pdb.set_trace()

    plt.figure()
    matplotlib.rcParams.update({'font.size': 22})

    plt.hist(existingCon,
             range(0,1+maxSynapses),
             density=True,
             align="left",
             color=self.getNeuronColor(preType))

    plt.xlabel("Number of " + connectionType)
    plt.ylabel('Probability density')
    #plt.title(preType + " to " + postType \
    #          + "(M=" + str(maxSynapses) \
    #          + ",m=" + '%.1f' % meanSynapses \
    #          + ",sl=" + '%.0f' % (sideLen*1e6) + ")")
    #plt.title(preType + " to " + postType \
    #          + " (total: " + str(np.sum(existingCon)) + ")")
    plt.title(self.neuronName(preType) + " to " + self.neuronName(postType))


    plt.tight_layout()

    plt.ion()
    plt.draw()
    if(self.showPlots):
      plt.show()

    plt.pause(0.001)
    figName = "Network-number-of-" + connectionType + "-from-" \
              + preType + "-to-" + postType + "-per-cell"

    self.saveFigure(plt,figName)


  ############################################################################

  def plotConnectionProbabilityParallel(self,
                                        preType=None,
                                        postType=None,
                                        nBins=86,dist3D=True):

    assert preType is not None
    assert postType is not None

    print("Plotting connection probability " + preType + " to " + postType)

    if(preType not in self.populations):
      print("plotConnectionProbability: " + str(preType) \
            + " is not in the simulation")
      return

    if(postType not in self.populations):
      print("plotConnectionProbability: " + str(postType) \
            + " is not in the simulation")
      return


    preID = self.populations[preType]
    postID = self.populations[postType]

    # We need to split the work between multiple workers
    import threading
    nThreads = 4
    threads = []

    self.workerData = []

    for iThread in range(0,nThreads):
      workerPreID = preID[range(iThread,len(preID),nThreads)]
      print("Worker " + str(iThread) + " PreID: " + str(workerPreID))
      t = threading.Thread(target=self.connectionProbabilityWrapper,
                           args=(workerPreID,postID,nBins,1000000.0,dist3D))
      threads.append(t)
      print("Starting " + t.getName())
      t.start()

    for t in threads:
      print("Joining " + t.getName())
      t.join()

    # Gather all the data
    dist = self.workerData[0][0]
    countCon = np.zeros((nBins,1))
    countAll = np.zeros((nBins,1))

    for data in self.workerData:
      countCon += data[2]
      countAll += data[3]

    Pcon = np.divide(countCon,countAll)

    # Now let's plot it
    matplotlib.rcParams.update({'font.size': 22})
    plt.figure()
    plt.plot(dist*1e6,Pcon)
    plt.xlabel("Distance ($\mu$m)")
    plt.ylabel("Connection probability")

    plt.title(self.neuronName(preType) + " to " \
              + self.neuronName(postType) + " connections")
    plt.tight_layout()

    plt.xlim([0, 250])
    plt.ylim([0, 1])

    plt.ion()
    plt.draw()

    if(self.showPlots):
      plt.show()

    plt.pause(0.001)
    figName = 'Network-distance-dependent-connection-probability-' \
              + str(preType) + "-to-" + str(postType)

    self.saveFigure(plt,figName)


  ############################################################################

  def plotConnectionProbability(self,
                                preType=None,
                                postType=None,
                                nBins=86,
                                nameStr="",
                                sideLen=None,
                                expMaxDist=[],
                                expData=[],
                                expDataDetailed = [],
                                dist3D=True,
                                volumeID=None,
                                xMax=250,
                                yMax=None,
                                connectionType="synapses",
                                drawStep=False):

    assert preType is not None
    assert postType is not None

    if(sideLen is None):
      sideLen = self.sideLen

    if(preType not in self.populations or postType not in self.populations):
      print("Missing " + preType + " or " + postType + " in network, " \
            + "skipping plot with their connectivity")
      return

    print("Plotting connection probability " + preType + " to " + postType \
          + " (" + str(connectionType) + ")")

    preID = self.populations[preType]
    postID = self.populations[postType]

    # We can in principle use all pairs, but here we restrict to just the
    # pairs who's post partner are in the centre
    # postID = self.centreNeurons(neuronID=postID,sideLen=sideLen)
    postID = self.getSubPop(volumeType=self.volumeType,
                            volumePart="centre",
                            sideLen=sideLen,
                            neuronID=postID,
                            volumeID=volumeID)

    if(preType not in self.populations):
      print("plotConnectionProbabilityChannels: " + str(preType) \
            + " is not in the simulation")
      return

    if(postType not in self.populations):
      print("plotConnectionProbabilityChannels: " + str(postType) \
            + " is not in the simulation")
      return

    if(len(expData) == 0 and len(expDataDetailed) > 0):
      expData = []
      for x in expDataDetailed:
        expData.append(x[0]/float(x[1]))

    if(len(expDataDetailed) == 0 and len(expData) > 0):
      expDataDetailed = [None for x in expData]

    (dist,Pcon,countCon,countAll) = \
      self.connectionProbability(preID,postID,nBins,dist3D=dist3D,
                                 connectionType=connectionType)

    # Now let's plot it

    #fig = plt.figure()
    fig,ax = plt.subplots(1)

    matplotlib.rcParams.update({'font.size': 24})


    pltCtr = 0

    # Add lines for experimental data and matching data for model
    for (dLimit,Pexp,expNum) in zip(expMaxDist,expData,expDataDetailed):
      cnt = 0
      cntAll = 0

      for (d,c,ca) in zip(dist,countCon,countAll):
        if(d <= dLimit):
          cnt += c
          cntAll += ca

      # Hack to avoid divide by zero
      cntAll[cntAll == 0] = 1

      Pmodel = float(cnt) / float(cntAll)

      print("P(d<" + str(dLimit) + ")=" + str(Pmodel))
      #ax = fig.get_axes()

      # Also add errorbars
      if(expNum is not None):
        P = expNum[0]/float(expNum[1])

        # Exp data specified should match
        assert Pexp is None or Pexp == P

        # stdExp = np.sqrt(P*(1-P)/expNum[1])

        # https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval
        # https://www.tandfonline.com/doi/abs/10.1080/01621459.1927.10502953
        # Wilson score
        # Wilson, Edwin B. "Probable inference, the law of succession, and statistical inference." Journal of the American Statistical Association 22.158 (1927): 209-212.
        ns = expNum[0]
        n = expNum[1]
        z = 1.96 # This gives us 95% confidence intervall
        barCentre = (ns + (z**2)/2) / (n + z *2)
        barHeight = z / (n + z**2) * np.sqrt((ns * (n - ns) / n + (z ** 2) / 4))

        #plt.errorbar(dLimit*1e6/2,P,stdExp,color="gray",
        #             elinewidth=1,capsize=5)
        plt.errorbar(dLimit*1e6/2,barCentre,barHeight,color="gray",
                     elinewidth=1,capsize=5)

        #import pdb
        #pdb.set_trace()

        #rectExpStd = patches.Rectangle((0,P-stdExp),
        #                               width=dLimit*1e6,height=2*stdExp,
        #                               alpha=0.2,linewidth=0,
        #                               color="red",fill=True)
        #ax.add_patch(rectExpStd)
      else:
        stdExp = 0

      if(Pexp is not None):
        plt.plot([0,dLimit*1e6],[Pexp, Pexp],
                 color=(0.8,0.3*pltCtr,0.3*pltCtr),linewidth=2)

        # Add a star also
        plt.plot(dLimit*1e6/2,Pexp,
                 color=(0.8,0.3*pltCtr,0.3*pltCtr),
                 marker="D",
                 markersize=10)

        pltCtr += 1
        plt.ion()
        plt.draw()

        if(self.showPlots):
          plt.show()

        #rectExp = patches.Rectangle((0,0), dLimit*1e6, Pexp, \
        #                            linewidth=2,color='red',fill=False)
        #ax.add_patch(rectExp)

      if(False): #
        # Show the binning from the model data
        plt.plot([0,dLimit*1e6],[Pmodel,Pmodel],color="blue",linewidth=2)

      #rect = patches.Rectangle((0,0), dLimit*1e6, P, \
      #                         linewidth=1,color='blue',fill=False)
      #ax.add_patch(rect)

      if(False):
        plt.text(dLimit*1e6,Pmodel,
                 "P(d<" + str(dLimit*1e6) + ")=%.3f" % Pmodel,size=9)


    # Draw the curve itself
    if(drawStep):
      plt.step(dist*1e6,Pcon,color='black',linewidth=2,where="post")
    else:
      dHalfStep = (dist[1]-dist[0])/2
      plt.plot((dist+dHalfStep)*1e6 ,Pcon,color='black',linewidth=2)

    plt.xticks(fontsize=14, rotation=0)
    plt.yticks(fontsize=14, rotation=0)

    labelSize = 22

    #if(dist3D):
    #  plt.xlabel("Distance ($\mu$m)",fontsize=labelSize)
    #else:
    #  plt.xlabel("2D Distance ($\mu$m)",fontsize=labelSize)

    # Hack to avoid divide by zero
    countAllB = countAll.copy()
    countAllB[countAllB == 0] = 1.0


    # This gives us 95% confidence intervall
    z = 1.96

    pCentre = np.array([(ns + (z**2)/2) / (n + z *2) \
                        for (ns,n) in zip(countCon,countAllB)]).flatten()
    pHeight = np.array([z / (n + z**2) \
                        * np.sqrt((ns * (n - ns) / n + (z ** 2) / 4)) \
                        for (ns,n) in zip(countCon,countAllB)]).flatten()

    # Use the last bin larger than xMax as the end
    dIdx = np.where(dist*1e6 > xMax)[0][0]

    pMin = pCentre - pHeight
    pMax = pCentre + pHeight

    if(drawStep):
      plt.fill_between(dist[:dIdx]*1e6, pMin[:dIdx], pMax[:dIdx],
                       color = 'grey', step="post",
                       alpha = 0.4)
    else:
      plt.fill_between((dist[:dIdx]+dHalfStep)*1e6, pMin[:dIdx], pMax[:dIdx],
                       color = 'grey', step=None,
                       alpha = 0.4)


    plt.xlabel("Distance ($\mu$m)",fontsize=labelSize)
    plt.ylabel("Con Prob (%)",fontsize=labelSize)


    # Adjust axis if needed
    # import pdb
    # pdb.set_trace()

    # Set plot limits
    #if(any(x is not None for x in expData)):
    #  if(max(plt.ylim()) < max(expData+stdExp)):
    #    plt.ylim([0, np.ceil(max(expData+stdExp)*10)/10])

    if(xMax is not None):
      plt.xlim([0, xMax])

    if(yMax is not None):
      plt.ylim([0, yMax])


    locs,labels = plt.yticks()
    newLabels = ["{0:g}".format(yy*100) for yy in locs]

    if(locs[0] < 0):
      locs = locs[1:]
      newLabels = newLabels[1:]

    if(locs[-1] > plt.ylim()[1]):
      locs = locs[:-1]
      newLabels = newLabels[:-1]

    plt.yticks(locs,newLabels)

    plt.title(self.neuronName(preType) + " to " + self.neuronName(postType))
    plt.tight_layout()
    plt.ion()
    plt.draw()

    if(self.showPlots):
      plt.show()

    plt.pause(0.001)

    if(dist3D):
      projText = '-3D-dist'
    else:
      projText = '-2D-dist'

    figName = 'Network-distance-dependent-connection-probability-' \
                + str(preType) + "-to-" + str(postType) \
                + "-" + str(connectionType) \
                + projText

    self.saveFigure(plt,figName)



  ############################################################################

  # expMaxDist=[50e-6,100e-6]
  # expData=[None, None] (none if no exp data, or values)
  # expDataDetailed=[(conA,totA),(conB,totB)]

  def plotConnectionProbabilityChannels(self,
                                        preType=None,
                                        postType=None,
                                        nBins=86,
                                        nameStr="",
                                        sideLen=None, # buffer around edges if neg
                                        expMaxDist=[],
                                        expData=[],
                                        expDataDetailed = [],
                                        dist3D=True,
                                        volumeID="Striatum"):

    if(preType not in self.populations):
      print("plotConnectionProbabilityChannels: " + str(preType) \
            + " is not in the simulation")
      return

    if(postType not in self.populations):
      print("plotConnectionProbabilityChannels: " + str(postType) \
            + " is not in the simulation")
      return


    if(sideLen is None):
      sideLen = self.sideLen

    if(len(expData) == 0 and len(expDataDetailed) > 0):
      expData = []
      for x in expDataDetailed:
        expData.append(x[0]/float(x[1]))

    if(len(expDataDetailed) == 0):
      expDataDetailed = []
      for x in expData:
        expDataDetailed.append(None)

    assert preType is not None
    assert postType is not None

    print("Plotting connection probability " + preType + " to " + postType)

    preID = self.populations[preType]
    postID = self.populations[postType]

    # We can in principle use all pairs, but here we restrict to just the
    # pairs who's post partner are in the centre
    # postID = self.centreNeurons(neuronID=postID,sideLen=sideLen)
    postID = self.getSubPop(volumeType=self.volumeType,
                            volumePart="centre",
                            sideLen=sideLen,
                            neuronID=postID,
                            volumeID=volumeID)

    (dist,PconWithin,PconBetween, \
     countConWithin, countConBetween, \
     countAllWithin, countAllBetween) = \
      self.connectionProbabilityChannels(preID,postID,nBins,dist3D=dist3D)

    # Now let's plot it

    #fig = plt.figure()
    fig,ax = plt.subplots(1)
    matplotlib.rcParams.update({'font.size': 22})

    # Draw the curve itself
    plt.plot(dist*1e6,PconWithin,color='black',linewidth=2)
    plt.plot(dist*1e6,PconBetween,color='grey',linewidth=2)

    labelSize=14
    if(dist3D):
      plt.xlabel("Distance ($\mu$m)",fontsize=labelSize)
    else:
      plt.xlabel("2D Distance ($\mu$m)",fontsize=labelSize)
    plt.ylabel("Connection probability",fontsize=labelSize)

    plt.xticks(fontsize=12, rotation=0)
    plt.yticks(fontsize=12, rotation=0)

    #plt.tick_params(axis='both', which='major', labelsize=10)
    #plt.tick_params(axis='both', which='minor', labelsize=8)

    #import pdb
    #pdb.set_trace()

    # Add lines for experimental data and matching data for model
    for (dLimit,Pexp,expNum) in zip(expMaxDist,expData,expDataDetailed):
      cntWithin = 0
      cntAllWithin = 0
      cntBetween = 0
      cntAllBetween = 0

      for (d,c,ca) in zip(dist,countConWithin,countAllWithin):
        if(d <= dLimit):
          cntWithin += c
          cntAllWithin += ca

      for (d,c,ca) in zip(dist,countConBetween,countAllBetween):
        if(d <= dLimit):
          cntBetween += c
          cntAllBetween += ca

      # Hack to avoid divide by zero, since the corresponding cntWithin
      # bin will also be zero, the resulting division is zero
      cntAllWithin[cntAllWithin == 0] = 1
      cntAllBetween[cntAllBetween == 0] = 1

      Pwithin = float(cntWithin) / float(cntAllWithin)
      Pbetween = float(cntBetween) / float(cntAllBetween)

      Ptotal = float(cntWithin+cntBetween) / float(cntAllWithin+cntAllBetween)

      print("Pwithin(d<" + str(dLimit) + ")=" + str(Pwithin))
      print("Pbetween(d<" + str(dLimit) + ")=" + str(Pbetween))
      print("Ptotal(d<" + str(dLimit) + ")=" + str(Ptotal))
      #ax = fig.get_axes()

      if(Pexp is not None):
        plt.plot([0,dLimit*1e6],[Pexp,Pexp],color="red",linewidth=2)

        oldRectFlag = False
        if(oldRectFlag):
          # Old rectangle
          rectExp = patches.Rectangle((0,0), dLimit*1e6, Pexp, \
                                      linewidth=2,color='red',fill=False,
                                      linestyle="--")
        # Also add errorbars
        if(expNum is not None):
          assert False, "Should change the error bars to wilson score!"
          P = expNum[0]/float(expNum[1])
          stdExp = np.sqrt(P*(1-P)/expNum[1])
          rectExpStd = patches.Rectangle((0,Pexp-stdExp),
                                         width=dLimit*1e6,height=2*stdExp,
                                         alpha=0.2,linewidth=0,
                                         color="red",fill=True)

          ax.add_patch(rectExpStd)

        if(oldRectFlag):
          # Add P-line on top
          ax.add_patch(rectExp)


      if(False):
        rectWithin = patches.Rectangle((0,0), dLimit*1e6, Pwithin, \
                                       linewidth=1,color='blue',fill=False)
        rectBetween = patches.Rectangle((0,0), dLimit*1e6, Pbetween, \
                                        linewidth=1,color='lightblue',
                                        fill=False)
        rectTotal = patches.Rectangle((0,0), dLimit*1e6, Ptotal, \
                                      linewidth=2,color='blue',fill=False)

        ax.add_patch(rectWithin)
        ax.add_patch(rectBetween)
        ax.add_patch(rectTotal)

      if(True):
        # Draw binned data as lines instead
        plt.plot([0,dLimit*1e6],[Pwithin, Pwithin],
                 color="blue",linewidth=1)
        plt.plot([0,dLimit*1e6],[Pbetween, Pbetween],
                 color="lightblue", linewidth=1)
        plt.plot([0,dLimit*1e6],[Ptotal,Ptotal],
                 color="blue", linewidth=2)

      if(False):
        Plt(dLimit*1e6,Pwithin,"Pw(d<" + str(dLimit*1e6) + ")=%.3f" % Pwithin,size=9)
        plt.text(dLimit*1e6,Pbetween,"Pb(d<" + str(dLimit*1e6) + ")=%.3f" % Pbetween,size=9)


    # Adjust axis if needed
    #import pdb
    #pdb.set_trace()

    if(any(x is not None for x in expData)):
      if(max(plt.ylim()) < max(expData)):
        plt.ylim([0, np.ceil(max(expData)*10)/10])

    plt.xlim([0, 250])
    #plt.xlim([0, 1000])

    plt.title(self.neuronName(preType) + " to " \
              + self.neuronName(postType) + " connections")

    plt.tight_layout()
    plt.ion()
    plt.draw()

    if(self.showPlots):
      plt.show()

    plt.pause(0.001)

    if(dist3D):
      projText = '-3D-dist'
    else:
      projText = '-2D-dist'

    figName = 'Network-distance-dependent-connection-probability-channels-' \
                + str(preType) + "-to-" + str(postType) \
                + projText

    self.saveFigure(plt,figName)


  ############################################################################

  def connectionProbabilityWrapper(self,preID,postID,nBins=86,dist3D=True):

    print("Worker started, preID: " + str(preID))

    (dist,Pcon,countCon,countAll) = self.connectionProbability(preID,
                                                               postID,
                                                               nBins,
                                                               dist3D=dist3D)

    self.workerData.append((dist,Pcon,countCon,countAll))

  ############################################################################

  # connectionType: "synapses" or "gapjunctions"

  def connectionProbability(self,
                            preID,
                            postID,
                            nBins=86,
                            nPoints=10000000.0,
                            dist3D=True,
                            connectionType="synapses"):

    # Count the connected neurons
    print("Counting connections")
    dist = np.linspace(0.0,1700.0e-6,num=nBins)
    deltaDist = dist[1]-dist[0]

    countCon = np.zeros((nBins,1))
    countAll = np.zeros((nBins,1))
    countRejected = 0

    nPerPair = nPoints/len(preID)

    if(connectionType == "synapses"):
      conMat = self.connectionMatrix
    elif(connectionType == "gapjunctions"):
      conMat = self.connectionMatrixGJ
    else:
      print("Unknown connectionType: " + str(connectionType))
      print("Please use 'synapses' or 'gapjunctions'")
      import pdb
      pdb.set_trace()


    # Make this loop use threads, to speed it up

    for xi,x in enumerate(preID):
      tA = timeit.default_timer()

      if(nPerPair - np.floor(nPerPair) > np.random.rand()):
        nPts = int(np.ceil(nPerPair))
      else:
        nPts = int(np.floor(nPerPair))

      try:
        poID = np.random.permutation(postID)[0:min(nPts,len(postID))]
      except:
        print("Something fucked up")
        import pdb
        pdb.set_trace()

      for y in poID: #postID:

        if(x == y):
          # Do not count self connections in statistics!!
          # This can lead to what looks like an artificial drop in
          # connectivity proximally
          continue

        if(dist3D):
          d = np.sqrt(np.sum((self.positions[x,:]-self.positions[y,:]) ** 2))
        else:
          d = np.sqrt(np.sum((self.positions[x,0:2]-self.positions[y,0:2]) ** 2))
          # We also need to check that z-distance is not too large

          # Gilad email 2017-11-21:
          # The 100 um is the lateral (XY) distance between somata of
          # recorded cells. In terms of the Z axis, the slice
          # thickness is 250 um but the recordings are all done in the
          # upper half of the slice due to visibility of the cells. So
          # lets say in a depth of roughly 40-110 um from the upper
          # surface of the slice.
          dz = np.abs(self.positions[x,2]-self.positions[y,2])

          # Using dzMax = 70, see comment above from Gilad.
          if(dz > 70e-6):
            # Skip this pair, too far away in z-depth
            countRejected += 1
            continue

        idx = int(np.floor(d/deltaDist))
        assert (idx < nBins), "Idx too large " + str(idx)

        if(conMat[x,y] > 0):
          countCon[idx] += 1

        countAll[idx] += 1

      tB = timeit.default_timer()

      if(self.debug):
        print(str(xi+1) + "/" + str(len(preID)) + " " + str(tB-tA) + "s")

    Pcon = np.divide(countCon,countAll)

    print("Requested: " + str(nPoints) + " calculated " + str(sum(countAll)))

    if(not dist3D):
      print("Rejected (too large z-depth): " + str(countRejected))

    #import pdb
    #pdb.set_trace()

    return (dist,Pcon,countCon,countAll)

  ############################################################################

  def connectionProbabilityChannels(self,
                                    preID,
                                    postID,
                                    nBins=86,
                                    nPoints=5000000.0,
                                    dist3D=True):

    neuronChannel = self.network["neuronChannel"]

    # Count the connected neurons
    print("Counting connections")
    dist = np.linspace(0.0,1700.0e-6,num=nBins)
    deltaDist = dist[1]-dist[0]

    countConWithinChannel = np.zeros((nBins,1))
    countAllWithinChannel = np.zeros((nBins,1))

    countConBetweenChannels = np.zeros((nBins,1))
    countAllBetweenChannels = np.zeros((nBins,1))

    countRejected = 0

    nPerPair = nPoints/len(preID)

    # Make this loop use threads, to speed it up

    for xi,x in enumerate(preID):
      tA = timeit.default_timer()

      if(nPerPair - np.floor(nPerPair) > np.random.rand()):
        nPts = int(np.ceil(nPerPair))
      else:
        nPts = int(np.floor(nPerPair))

      try:
        poID = np.random.permutation(postID)[0:min(nPts,len(postID))]
      except:
        print("Something fucked up")
        import pdb
        pdb.set_trace()

      for y in poID: #postID:

        if(x == y):
          # Dont include self-self
          continue

        if(dist3D):
          d = np.sqrt(np.sum((self.positions[x,:]-self.positions[y,:]) ** 2))
        else:
          d = np.sqrt(np.sum((self.positions[x,0:2]-self.positions[y,0:2]) ** 2))
          # We also need to check that z-distance is not too large

          # Gilad email 2017-11-21:
          # The 100 um is the lateral (XY) distance between somata of
          # recorded cells. In terms of the Z axis, the slice
          # thickness is 250 um but the recordings are all done in the
          # upper half of the slice due to visibility of the cells. So
          # lets say in a depth of roughly 40-110 um from the upper
          # surface of the slice.
          dz = np.abs(self.positions[x,2]-self.positions[y,2])

          # Using dzMax = 70, see comment above from Gilad.
          if(dz > 70e-6):
            # Skip this pair, too far away in z-depth
            countRejected += 1
            continue

        idx = int(np.floor(d/deltaDist))

        if(idx >= nBins):
          countRejected += 1
          continue

        # assert (idx < nBins), "Idx too large " + str(idx)

        sameChannelFlag = (neuronChannel[x] == neuronChannel[y])

        if(self.connectionMatrix[x,y] > 0):
          if(sameChannelFlag):
            countConWithinChannel[idx] += 1
          else:
            countConBetweenChannels[idx] += 1

        if(sameChannelFlag):
          countAllWithinChannel[idx] += 1
        else:
          countAllBetweenChannels[idx] += 1

      tB = timeit.default_timer()

      if(self.debug):
        print(str(xi+1) + "/" + str(len(preID)) + " " + str(tB-tA) + "s")

    PconWithinChannel = np.divide(countConWithinChannel,countAllWithinChannel)
    PconBetweenChannels = np.divide(countConBetweenChannels, \
                                    countAllBetweenChannels)

    print("Requested: " + str(nPoints) + " calculated " + str(sum(countAllWithinChannel) + sum(countAllBetweenChannels)))

    if(not dist3D):
      print("Rejected (too large z-depth): " + str(countRejected))

    return (dist,PconWithinChannel,PconBetweenChannels, \
            countConWithinChannel,countConBetweenChannels, \
            countAllWithinChannel,countAllBetweenChannels)



  ############################################################################

  def numIncomingConnections(self,neuronType,preType,sideLen=None,
                             volumeID="Striatum",
                             connectionType="synapses"):

    if(sideLen is None):
      sideLen=100e-6

    print("Calculating number of incoming connections " + preType \
      + " -> " + neuronType)

    # Only use post synaptic cell in central part of structure,
    # to minimize edge effect
    #neuronID = self.centreNeurons(neuronID=self.populations[neuronType],
    #                              sideLen=sideLen)

    neuronID = self.getSubPop(volumeType=self.volumeType,
                              volumePart="centre",
                              sideLen=sideLen,
                              neuronID=self.populations[neuronType],
                              volumeID=volumeID)


    preID = self.populations[preType]

    print("#pre = " + str(len(preID)) + ", #post = " + str(len(neuronID)))

    if(connectionType == "synapses"):
      conMat = self.connectionMatrix
    elif(connectionType == "gapjunctions"):
      conMat = self.connectionMatrixGJ
    else:
      print("Unknown connectionType: " + str(connectionType))
      print("Please use 'synapses' or 'gapjunctions'")
      import pdb
      pdb.set_trace()

    nCon = np.zeros((len(neuronID),1))
    nSyn = np.zeros((len(neuronID),1))

    for i,nID in enumerate(neuronID):
      cons = conMat[:,nID][preID,:]
      nCon[i] = np.sum(cons > 0)
      nSyn[i] = np.sum(cons)

    return (nCon,nSyn)

  ############################################################################

  def plotIncomingConnections(self,preType,neuronType,sideLen=None,
                              nameStr="",
                              connectionType="synapses",
                              fig=None,colour=None,histRange=None):

    if(preType not in self.populations):
      print("plotIncomingConnections: " + str(preType) \
            + " is not in the simulation")
      return

    if(neuronType not in self.populations):
      print("plotIncomingConnections: " + str(neuronType) \
            + " is not in the simulation")
      return


    (nCon,nSyn) = self.numIncomingConnections(neuronType=neuronType,
                                              preType=preType,
                                              sideLen=sideLen,
                                              connectionType=connectionType)

    # Plotting number of connected neighbours
    if(fig is None):
      fig = plt.figure()
      histType = 'bar'
    else:
      plt.sca(fig.axes[0])
      histType = 'step'
      
    matplotlib.rcParams.update({'font.size':22})

    if(sum(nCon) > 0):
      if(max(nCon) > 100):
        binSize=10
      else:
        binSize=1

      if(colour is None):
        colour = self.getNeuronColor(preType)

      if(histRange is None):
        histRange = range(0,int(np.max(nCon)),binSize)
        
      plt.hist(nCon,histRange,
               align="left",density=True,
               color=colour,histtype=histType)

    plt.xlabel("Number of connected neighbours")
    plt.ylabel("Probability density")
    plt.title(self.neuronName(preType) + " connecting to " \
              + self.neuronName(neuronType))
    plt.tight_layout()
    plt.ion()
    plt.draw()

    if(self.showPlots):
      plt.show()

    plt.pause(0.001)

    figName = "Network-" + connectionType + "-input-to-" \
      + neuronType + "-from-" + preType

    self.saveFigure(plt,figName)


    # Plotting number of input synapses
    plt.figure()
    matplotlib.rcParams.update({'font.size':22})

    if(sum(nSyn) > 0):
      if(max(nSyn) > 700):
        binSize = 100
      elif(max(nSyn) < 10):
        binSize = 1
      elif(max(nSyn) < 40):
        binSize = 5
      else:
        binSize = 10
      plt.hist(nSyn,range(0,int(np.max(nSyn)),binSize),
               align="left",density=True,
               color=self.getNeuronColor(preType))

    plt.xlabel("Number of incoming " + connectionType)
    plt.ylabel("Probability density")
    plt.title(self.neuronName(preType) + " " + connectionType \
              + " on " + self.neuronName(neuronType))
    plt.tight_layout()
    plt.ion()
    plt.draw()

    if(self.showPlots):
      plt.show()

    plt.pause(0.001)
    figName = "Network-" + connectionType + "-to-" + neuronType + "-from-" + preType

    self.saveFigure(plt,figName)

    return fig


  ############################################################################

  def findSynapses(self, neuronID, maxSynapses=10000, chunkSize=1000000):

    print("Finding synapses between: " + str(neuronID))

    synapses = np.zeros((maxSynapses,8))
    synapseCtr = 0

    nRows = self.data["network/synapses"].shape[0]
    chunkSize = min(nRows,chunkSize)

    nSteps = int(np.ceil(nRows/chunkSize))
    rowStart = 0

    for rowEnd in np.linspace(chunkSize,nRows,nSteps,dtype=int):

      print("Rows: " + str(rowStart) + ":" + str(rowEnd) \
        + " total: " + str(nRows))

      syn = self.data["network/synapses"][rowStart:rowEnd,:]

      for synRow in syn:
        if(synRow[0] in neuronID and synRow[1] in neuronID):
          synapses[synapseCtr,:] = synRow
          synapseCtr += 1

    print("Found " + str(synapseCtr) + " synapses")

    print(str(synapses[0:synapseCtr,:]))

    return synapses[0:synapseCtr,:]

  ############################################################################

  # Plots neuronID neuron, and all presynaptic partners

  def plotNeurons(self,neuronID,showSynapses=True,plotPreNeurons=True):

    axis = None

    # Finding synapses, this might take time
    (synapses,synapseCoords) = self.networkLoad.findSynapses(postID=neuronID)

    assert (synapses[:,1] == neuronID).all(), "!!!! Wrong synapses extracted"

    neurons = dict([])

    postNeuron = self.networkLoad.loadNeuron(neuronID)
    axis = postNeuron.plotNeuron(axis=axis,plotAxon=False, \
                                  plotDendrite=True)

    plottedNeurons = [neuronID]

    if(plotPreNeurons):
      for synRow in synapses:

        srcID = int(synRow[0])

        if(srcID not in plottedNeurons):
          neurons[srcID] = self.networkLoad.loadNeuron(srcID)
          axis = neurons[srcID].plotNeuron(axis=axis,plotAxon=True, \
                                           plotDendrite=False)
          plottedNeurons.append(srcID)

    if(showSynapses):
      axis.scatter(synapseCoords[:,0],
                   synapseCoords[:,1],
                   synapseCoords[:,2],c='red')

    plt.ion()
    plt.draw()

    if(self.showPlots):
      plt.show()

    plt.pause(0.001)

    return axis


  ############################################################################

  def findLatestFile(self):

    files = glob('save/network-connect-voxel-pruned-synapse-file-*.hdf5')
#    files = glob('save/network-connect-synapse-file-*.hdf5')

    modTime = [os.path.getmtime(f) for f in files]
    idx = np.argsort(modTime)

    print("Using the newest file: " + files[idx[-1]])

    return files[idx[-1]]

  ############################################################################

  # Loop through all synapses, create distance histogram

  def synapseDist(self,sideLen=None,volumeID="Striatum"):

    if(sideLen is None):
      sideLen = self.sideLen

    tA = timeit.default_timer()

    # cornerID = self.cornerNeurons(sideLen=sideLen)
    cornerID = self.getSubPop(volumeType=self.volumeType,
                              volumePart="corner",
                              sideLen=sideLen,
                              neuronID=None,
                              volumeID=volumeID)

    isCorner = np.zeros((self.nNeurons,),dtype=bool)
    isCorner[cornerID] = 1

    print("Calculating synapse distance histogram")

    nBins = 200
    maxDist = 1000e-6
    binWidth = maxDist/nBins
    binScaling = 1e-6 / binWidth # To avoid a division

    neuronType = [n["name"].split("_")[0] for n in self.network["neurons"]]
    self.allTypes = np.unique(neuronType).tolist()

    # To be quicker, we temporary define neuronTypeIDs
    self.neuronTypeID = [self.allTypes.index(x) for x in neuronType]

    self.dendPositionBin = dict([])
    self.dendPositionEdges = np.linspace(0,maxDist,nBins+1)

    for preTypeID in range(0,len(self.allTypes)):
      for postTypeID in range(0,len(self.allTypes)):
        self.dendPositionBin[(preTypeID,postTypeID)] = np.zeros((nBins,))

    print("Creating dist histogram")

    if(self.volumeType == "full"):
      synIncrement = 1
    elif(self.volumeType == "cube"):
      synIncrement = 8    # Becaues we use 1/8th of volume,
                          # count all synapses 8 times
    else:
      print("Unknown volume type: " + str(self.volumeType))
      import pdb
      pdb.set_trace()

    nSynapses = self.data["network"]["synapses"].shape[0]
    chunkSize = 1000000
    startIdx = 0
    endIdx = min(chunkSize,nSynapses)

    # Outer while loop is to split the data processing into smaller chunks
    # otherwise we use too much memory.
    while(startIdx < nSynapses):

      assert startIdx <= endIdx, "startIdx = " + str(startIdx) + ", endIdx = " + str(endIdx)

      allPreID = self.data["network"]["synapses"][startIdx:endIdx,0]
      allPostID = self.data["network"]["synapses"][startIdx:endIdx,1]

      # Dist to soma on dendrite
      allDistIdx = np.array(np.floor(self.data["network"]["synapses"][startIdx:endIdx,8] \
                                     * binScaling),
                            dtype=np.int)

      lastPre = None
      lastPost = None

      print("nSynapses = " + str(nSynapses) + ", at " + str(startIdx))

      for i in range(0,allPreID.shape[0]):

        if(not isCorner[allPostID[i]]):
          # We only want to include the corner neurons
          continue

        preTypeID = self.neuronTypeID[allPreID[i]]
        postTypeID = self.neuronTypeID[allPostID[i]]
        idx = allDistIdx[i]

        if(preTypeID != lastPre or postTypeID != lastPost):
          lastPre = preTypeID
          lastPost = postTypeID
          binCount = self.dendPositionBin[(preTypeID,postTypeID)]

        # +8 for cube corners, and +1 for full
        binCount[idx] += synIncrement


      startIdx += chunkSize
      endIdx += chunkSize
      endIdx = min(endIdx,nSynapses)


    tB = timeit.default_timer()

    print("Created distance histogram (optimised) in " + str(tB-tA) + " seconds")


  ############################################################################

  def plotSynapseCumDistSummary(self,pairList):

    matplotlib.rcParams.update({'font.size': 22})
    plt.figure()

    figName = "SynapseCumDistSummary"
    legendText = []

    for pair in pairList:

      try:
        pairID = (self.allTypes.index(pair[0]),
                  self.allTypes.index(pair[1]))
      except:
        print("Missing pair: " + str(pair))
        continue

      if(pairID not in self.dendPositionBin):
        print("Missing cum dist information for " + str(pair))
        continue

      if(sum(self.dendPositionBin[pairID]) == 0):
        print("Empty cum dist data for " + str(pair))
        continue

      cumDist = np.cumsum(self.dendPositionBin[pairID])  \
                 /np.sum(self.dendPositionBin[pairID])

      # Select range to plot
      endIdx = np.where(self.dendPositionEdges <= 400e-6)[0][-1]

      try:
        preType = pair[0]
        postType = pair[1]

        if(preType in self.neuronColors):
          plotCol = self.neuronColors[preType]
        else:
          plotCol = self.neuronColors["default"]

        # Hack: dSPN and iSPN are overlapping, need to make dSPN visible
        if(preType == "dSPN"):
          linewidth=7
        else:
          linewidth=3

        plt.plot(self.dendPositionEdges[:endIdx]*1e6,
                 cumDist[:endIdx],
                 linewidth=linewidth,color=plotCol)


        figName += "_" + preType + "-" + postType

        legendText.append(self.neuronName(preType) + "-" \
                          + self.neuronName(postType))

      except:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)
        import pdb
        pdb.set_trace()

    fontP = matplotlib.font_manager.FontProperties()
    fontP.set_size('small')
    plt.legend(legendText,prop=fontP)
    plt.xlabel('Distance from soma ($\mu$m)')
    plt.ylabel('Cumulative distrib.')

    figName += ".pdf"

    self.saveFigure(plt,figName)

  ############################################################################

  def plotSynapseCumDist(self):

    #import pdb
    #pdb.set_trace()

    for pair in self.dendPositionBin:

      if(sum(self.dendPositionBin[pair]) == 0):
        # Nothing in the bins, skip this plot
        continue

      #import pdb
      #pdb.set_trace()

      cumDist = np.cumsum(self.dendPositionBin[pair])  \
                 /np.sum(self.dendPositionBin[pair])

      idx = np.where(cumDist < 0.5)[0][-1]
      print(self.allTypes[pair[0]] + " to " + self.allTypes[pair[1]] \
            + " 50% of synapses are within " \
            + str(self.dendPositionEdges[idx]*1e6) + " micrometer")

      # Dont plot the full range
      endIdx = np.where(self.dendPositionEdges <= 300e-6)[0][-1]


      try:
        preType = pair[0]
        postType = pair[1]

        matplotlib.rcParams.update({'font.size': 22})

        plt.figure()
        plt.plot(self.dendPositionEdges[:endIdx]*1e6,
                 cumDist[:endIdx],
                 linewidth=3)
        plt.xlabel('Distance from soma ($\mu$m)')
        plt.ylabel('Cumulative distrib.')
        plt.title('Synapses ' + self.neuronName(self.allTypes[preType]) \
                  + " to " + self.neuronName(self.allTypes[postType]))
        plt.tight_layout()

        plt.ion()

        if(self.showPlots):
          plt.show()

        plt.draw()
        plt.pause(0.0001)
        figName = "SynapseCumulativeDistribution-" \
                    + self.neuronName(self.allTypes[preType]) + "-to-" \
                    + self.neuronName(self.allTypes[postType])

        self.saveFigure(plt,figName)


      except Exception as e:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)

        print("fucked up")
        import pdb
        pdb.set_trace()



  ############################################################################

  def plotSynapseDist(self,densityFlag=False,sideLen=None):

    if(sideLen is None):
      sideLen = self.sideLen

    #import pdb
    #pdb.set_trace()

    dendHistTot = self.dendriteDensity(nBins=len(self.dendPositionEdges),
                                       binWidth=self.dendPositionEdges[1],
                                       sideLen=sideLen)


    # self.dendPositionBin calculated by synapseDist(), done at time of init

    for pair in self.dendPositionBin:

      if(sum(self.dendPositionBin[pair]) == 0):
        # Nothing in the bins, skip this plot
        continue

      endIdx = np.where(self.dendPositionEdges <= 400e-6)[0][-1]

      try:
        preType = self.allTypes[pair[0]]
        postType = self.allTypes[pair[1]]

        plt.rcParams.update({'font.size': 16})

        plt.figure()
        if(densityFlag):

          # Skip first bin if we plot densities, since synapses from soma
          # are in first bin, but contribute no dendritic length

          plt.plot(self.dendPositionEdges[1:endIdx]*1e6,
                   np.divide(self.dendPositionBin[pair][1:endIdx],
                             dendHistTot[postType][1:endIdx]*1e6))
          plt.ylabel('Synapse/micrometer')
          plt.xlabel('Distance from soma ($\mu$m)')

          plt.title('Synapse density ' + self.neuronName(preType)\
                    + " to " + self.neuronName(postType))


        else:
          plt.plot(self.dendPositionEdges[:endIdx]*1e6,
                   self.dendPositionBin[pair][:endIdx])
          plt.ylabel('Synapse count')
          plt.xlabel('Distance from soma ($\mu$m)')
          plt.ylim([0,np.ceil(np.max(self.dendPositionBin[pair][:endIdx]))])

          plt.title('Synapses ' + self.neuronName(preType) \
                    + " to " + self.neuronName(postType))

        plt.tight_layout()

        plt.ion()

        if(self.showPlots):
          plt.show()

        plt.draw()
        plt.pause(0.0001)

        if(densityFlag):
          figName = "SynapseDistribution-density-dend-" \
                      + preType + "-to-" + postType
        else:
          figName = "SynapseDistribution-dend-" \
                      + preType + "-to-" + postType


        self.saveFigure(plt,figName)


      except Exception as e:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)

        print("fucked up")
        import pdb
        pdb.set_trace()

  ############################################################################


  ############################################################################

  # We want to get the number of synapses per unit length
  # We have the number of synapses at different distances from some (arc length)
  # So we want to calculate dendritic length, and then divide the synapses by
  # that

  def dendriteDensity(self,nBins,binWidth,sideLen=None,volumeID="Striatum"):

    if(sideLen is None):
      sideLen = self.sideLen

    print("Using nBins = " + str(nBins) + ", binWidth = " + str(binWidth))

    # 1. Loop through all morphologies
    morphFiles = set([self.data["morphologies"][name]["location"][()] \
                      for name in self.data["morphologies"]])

    morphHist = dict([])
    for mFile in morphFiles:
      morphHist[mFile] = self._dendDensity(mFile,nBins,binWidth)


    neuronTypes = set([m.split("_")[0] \
                       for m in self.data["morphologies"]])


    # 2. Sum all the histograms together, with the right multiplicity

    dendHistTot = dict([])
    for nt in neuronTypes:
      dendHistTot[nt] = np.zeros((nBins,))


    # cornerID = self.cornerNeurons(sideLen=sideLen)
    cornerID = self.getSubPop(volumeType=self.volumeType,
                              volumePart="corner",
                              sideLen=sideLen,
                              neuronID=None,
                              volumeID=volumeID)


    isCorner = np.zeros((self.nNeurons,),dtype=bool)
    isCorner[cornerID] = 1

    # OBS, neuronID unique, but many neurons can have same
    # morphology and properties and thus same neuronName
    # So the same name can appear twice in this list
    for nID, name in zip(self.data["network"]["neurons"]["neuronID"],
                         self.data["network"]["neurons"]["name"]):

      if(not isCorner[nID]):
        # Only post synaptic corner neurons included
        continue

      mFile = self.data["morphologies"][name]["location"][()]
      nType = name.decode().split("_")[0]
      dendHistTot[nType] += morphHist[mFile]

    return dendHistTot


  ############################################################################

  def _dendDensity(self,swcFile,nBins,binWidth):

    from .Neuron_morphology import NeuronMorphology

    dendHist = np.zeros((nBins,))

    nMorph = NeuronMorphology(swc_filename=swcFile)

    print("Parsing dendrite histogram : " + swcFile)

    # 0,1,2: x,y,z  3: radie, 4: dist to soma
    distToSoma = nMorph.dend[:,4]
    compartmentLength = np.zeros((nMorph.dend.shape[0],))

    for idx,link in enumerate(nMorph.dendLinks):
      compartmentLength[idx] = np.linalg.norm(nMorph.dend[link[0],:3] \
                                              - nMorph.dend[link[1],:3])
      distToSoma[idx] = (nMorph.dend[link[0],4] + nMorph.dend[link[1],4])/2

    for d,cl in zip(distToSoma,compartmentLength):
      idx = int(d/binWidth)
      dendHist[idx] += cl

    return dendHist


  ############################################################################

  # This

  def virtualAxonSynapses(self,postNeuronType):

    virtIdx = np.where([n["virtualNeuron"] for n in self.network["neurons"]])[0]

    postIdx = self.populations[postNeuronType]

    virtSyn = dict([])

    for vIdx in virtIdx:
      nType = self.network["neurons"][vIdx]["name"].split("_")[0]

      synMat = self.connectionMatrix[vIdx,:][:,postIdx]

      if(synMat.count_nonzero() == 0):
        # No synapses here, skip plot
        continue

      synMat = synMat.todense()
      synMat = synMat[np.where(synMat > 0)]

      if(nType not in virtSyn):
        virtSyn[nType] = [synMat]
      else:
        virtSyn[nType].append(synMat)

    for axonType in virtSyn:
      plt.figure()

      vSyn = np.concatenate(virtSyn[axonType]).flatten().transpose()

      # We dont want to show the count for zeros
      nBins = np.max(vSyn) + 1
      virtSynBins = np.array(range(1,nBins))

      plt.hist(vSyn,bins=virtSynBins,align="left")
      plt.xlabel("Synapses (only connected pairs)")
      plt.ylabel("Count")
      plt.title("Synapses from " + self.neuronName(axonType) \
                + " to " + self.neuronName(postNeuronType) \
                + " (nSyn=" + str(np.sum(vSyn)) + ")" )

      plt.tight_layout()

      plt.ion()
      plt.draw()

      if(self.showPlots):
        plt.show()

      figName = "VirtuaAxon-synapses-" + axonType + "-to-" \
                  + postNeuronType


      self.saveFigure(plt,figName)


      #import pdb
      #pdb.set_trace()

  ############################################################################

  def countMotifs(self,typeA,typeB,typeC,nRepeats=1000000):

    print("Counting motivs between " + str(typeA) + ", " + str(typeB) \
          + " " + str(typeC) + ". " + str(nRepeats) + " repeats.")

    IDA = self.populations[typeA]
    IDB = self.populations[typeB]
    IDC = self.populations[typeC]

    # Init counter
    motifCtr = np.zeros((64,),dtype=int)

    # bit 1 : A->B (1)
    # bit 2 : A<-B (2)
    # bit 3 : A->C (4)
    # bit 4 : A<-C (8)
    # bit 5 : B->C (16)
    # bit 6 : B<-C (32)

    # Faster to generate all int at once
    iAall = np.random.randint(len(IDA),size=(nRepeats,))
    iBall = np.random.randint(len(IDB),size=(nRepeats,))
    iCall = np.random.randint(len(IDC),size=(nRepeats,))

    for iRep in range(0,nRepeats):

      if(iRep % 100000 == 0 and iRep > 0):
        print("rep: " + str(iRep))

      iA = iAall[iRep]
      iB = iBall[iRep]
      iC = iCall[iRep]

      # In case the same neuron was picked twice, redo sampling
      while(iA == iB or iB == iC or iC == iA):
        iA = np.random.randint(len(IDA))
        iB = np.random.randint(len(IDB))
        iC = np.random.randint(len(IDC))

      idx =   int(self.connectionMatrix[iA,iB] > 0)*1 \
            + int(self.connectionMatrix[iB,iA] > 0)*2 \
            + int(self.connectionMatrix[iA,iC] > 0)*4 \
            + int(self.connectionMatrix[iC,iA] > 0)*8 \
            + int(self.connectionMatrix[iB,iC] > 0)*16 \
            + int(self.connectionMatrix[iC,iB] > 0)*32

      motifCtr[idx] += 1

    return (motifCtr,typeA,typeB,typeC)

  ############################################################################

  def analyseSingleMotifs(self,neuronType,nRepeats=10000000):

    (motifCtr,tA,tB,tC) = self.countMotifs(typeA=neuronType,
                                           typeB=neuronType,
                                           typeC=neuronType,
                                           nRepeats=nRepeats)


    # !!! For debug
    # motifCtr = [bin(x).count('1') for x in np.arange(0,64)]

    # No connections
    noCon = motifCtr[0]

    # One connection
    oneCon = motifCtr[1] + motifCtr[2] + motifCtr[4] + motifCtr[8] \
           + motifCtr[16] + motifCtr[32]

    # Two connections
    twoConDiverge = motifCtr[1+4] + motifCtr[2+16] + motifCtr[8+32]
    twoConConverge = motifCtr[2+8] + motifCtr[1+32] + motifCtr[4+16]
    twoConLine = motifCtr[1+16] + motifCtr[2+4] + motifCtr[4+32] \
               + motifCtr[8+1] + motifCtr[16+8] + motifCtr[32+2]
    twoConBi = motifCtr[1+2] + motifCtr[4+8] + motifCtr[16+32]

    # Three connections
    threeCircle = motifCtr[1+16+8] + motifCtr[2+4+32]
    threeCircleFlip = motifCtr[2+16+8] + motifCtr[1+32+8] + motifCtr[1+16+4] \
                    + motifCtr[1+4+32] + motifCtr[2+8+32] + motifCtr[2+4+16]
    threeBiDiv = motifCtr[1+2+4]   + motifCtr[1+2+16] \
               + motifCtr[4+8+1]   + motifCtr[4+8+32] \
               + motifCtr[16+32+2] + motifCtr[16+32+8]
    threeBiConv = motifCtr[1+2+8] + motifCtr[1+2+32] \
                + motifCtr[4+8+2] + motifCtr[4+8+16] \
                + motifCtr[16+32+1] + motifCtr[16+32+4]

    # Four connections
    fourDoubleBi = motifCtr[1+2+4+8] + motifCtr[1+2+16+32] + motifCtr[4+8+16+32]
    fourBiDiverge = motifCtr[1+2+4+16] + motifCtr[4+8+1+32] \
                    + motifCtr[16+32+2+8]
    fourBiConverge = motifCtr[1+2+8+32] + motifCtr[4+8+2+16] \
                    + motifCtr[16+32+1+4]
    fourBiCycle = motifCtr[1+2+4+32] + motifCtr[1+2+16+8] \
                + motifCtr[4+8+1+16] + motifCtr[4+8+32+2] \
                + motifCtr[16+32+2+4] + motifCtr[16+32+8+1]

    # Five connections
    fiveCon = motifCtr[2+4+8+16+32] \
            + motifCtr[1+4+8+16+32] \
            + motifCtr[1+2+8+16+32] \
            + motifCtr[1+2+4+16+32] \
            + motifCtr[1+2+4+8+32] \
            + motifCtr[1+2+4+8+16]

    # Six connections
    sixCon = motifCtr[1+2+4+8+16+32]


    conData = [("No connection", noCon),
               ("One connection", oneCon),
               ("Two connections, diverge", twoConDiverge),
               ("Two connections, converge", twoConConverge),
               ("Two connections, line", twoConLine),
               ("Two connections, bidirectional", twoConBi),
               ("Three connections, circle", threeCircle),
               ("Three connections, circular, one flipped", threeCircleFlip),
               ("Three connections, one bidirectional, one away", threeBiDiv),
               ("Three connections, one bidirectional, one towards",
                threeBiConv),
               ("Four connections, two bidirectional", fourDoubleBi),
               ("Four connections, one bi, two away", fourBiDiverge),
               ("Four connections, one bi, two towards", fourBiConverge),
               ("Four connections, cycle with one bidirectional", fourBiCycle),
               ("Five connections", fiveCon),
               ("Six connections", sixCon)]

    print("Motif analysis:")
    for (name,data) in conData:
      print(name + " " + str(100*data/nRepeats) + "% (" + str(data) + ")")


  ############################################################################


  # If A is connected to B and C, what is probability that
  # B and C are connected?

  def simpleMotif(self,typeA,typeB,typeC,nRep=1000):

    IDA = self.populations[typeA]
    IDB = self.populations[typeB]
    IDC = self.populations[typeC]

    # We need to reduce the IDA population, otherwise it will be too slow
    if(nRep < len(IDA)):
      try:
        rPermIdx = np.random.permutation(len(IDA))[:nRep]
        IDA = [IDA[x] for x in rPermIdx]
      except:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)
        import pdb
        pdb.set_trace()


    conBC = 0
    conCB = 0
    conBi = 0
    allCtr = 0

    nA = len(IDA)
    ctrA = 0


    for iA in IDA:

      ctrA += 1

      if(ctrA % 100 == 0):
        print(str(ctrA) + "/" + str(nA))

      matRow = self.connectionMatrix[iA,:]

      try:
        iBall = IDB[matRow[0,IDB].nonzero()[1]]
        iCall = IDC[matRow[0,IDC].nonzero()[1]]

      except:
        print("Ugh")
        import pdb
        pdb.set_trace()

      for iB in iBall:
        for iC in iCall:
          BC = self.connectionMatrix[iB,iC]
          CB = self.connectionMatrix[iC,iB]

          if(BC > 0):
            conBC += 1
            if(CB > 0):
              conCB += 1
              conBi += 1
          elif(CB > 0):
            conCB += 1

          allCtr += 1


    print("If " + str(typeA) + " connected to " + str(typeB) + " " \
          + str(typeC) + ":")
    print("P(" + typeB + "->" + typeC + ") = " + str((100.0*conBC)/allCtr)+"%")
    print("P(" + typeC + "->" + typeB + ") = " + str((100.0*conCB)/allCtr)+"%")
    print("P(" + typeB + "<->" + typeC + ") = " + str((100.0*conBi)/allCtr)+"%")


    return (conBC,conCB,conBi,allCtr)

  ############################################################################

  # Pick a post synaptic neuron, find out the distance to its closest
  # presynaptic neighbour

  def nearestPreNeighbourDistance(self,preType,postType,rabiesRate=1.0,
                                  nameStr=""):

    if(preType not in self.populations):
      print("nearestPreNeighbourDistance: " + str(preType) \
            + " is not in the simulation")
      return

    if(postType not in self.populations):
      print("nearestPreNeighbourDistance: " + str(postType) \
            + " is not in the simulation")
      return


    #postPop = self.getSubPop(neuronID=self.populations[postType],
    #                         volumePart="centre")

    postPop = self.populations[postType]
    prePop = self.populations[preType]


    # Assume 300 micrometer slice
    maxZ = np.max(self.positions[postPop,2])
    minZ = np.max(self.positions[postPop,2])

    sliceMinZ = (maxZ+minZ)/2 - 150e-6
    sliceMaxZ = (maxZ+minZ)/2 + 150e-6

    slicePop = np.where(np.bitwise_and(sliceMinZ <= self.positions[:,2],
                                       self.positions[:,2] <= sliceMaxZ))[0]

    # Only look at neurons in slice
    prePop = np.intersect1d(prePop,slicePop)
    postPop = np.intersect1d(postPop,slicePop)

    # conIdx = np.sum(self.connectionMatrix[prePop,postPop],axis=1)

    nearestDist = np.nan * np.zeros((len(postPop),))

    for c,postID in enumerate(postPop):

      # Calculate distance to all connected neighbours
      idx = np.where(self.connectionMatrix[prePop,postID].todense() > 0)[0]
      conIdx = prePop[idx]

      if(rabiesRate is not None and rabiesRate < 1):
        keepIdx = np.where(np.random.rand(len(conIdx)) < rabiesRate)[0]
        conIdx = conIdx[keepIdx]

      d = np.sqrt(np.sum((self.positions[conIdx,:3] \
                         - self.positions[postID,:3])**2,
                        axis=1))

      if(len(d) > 0):
        nearestDist[c] = np.min(d) * 1e6 # micrometers

    maxDist = np.nanmax(nearestDist)

    plt.figure()
    matplotlib.rcParams.update({"font.size":22})

    plt.hist(nearestDist,np.arange(0,maxDist+25,25))

    plt.xlabel("Distance")
    if(rabiesRate is not None and rabiesRate < 1):
      plt.ylabel("Count (Rabies rate: " + str(rabiesRate*100) + "%)")
    else:
      plt.ylabel("Count")
    plt.title("Nearest presynaptic neighbour " \
              + self.neuronName(preType) + " to " + self.neuronName(postType))

    # Data from Sabatini 2016
    if(preType == "LTS" and (postType == "dSPN" or postType == "iSPN")):
      SabatiniLTS = np.genfromtxt("DATA/LTS-nearest-neighbour-points-Sabatini2016.csv", delimiter = ",")
      LTSpoints = SabatiniLTS[:,1] * 1e3 # Get in micrometer
      plt.hist(LTSpoints,color='r',histtype="step")

    plt.ion()
    plt.draw()

    if(self.showPlots):
      plt.show()

    plt.pause(0.001)


    figName = "Nearest-presynaptic-slice-neighbour-to-" \
      + str(postType) + "-from-" + str(preType) + "-ID-" \
      +str(self.network["SlurmID"]) + nameStr + ".pdf"

    self.saveFigure(plt,figName)

    if(self.closePlots):
      time.sleep(1)
      plt.close()



  ############################################################################

  # !!! Doh, just realised using the connection matrix might have been smarter

  def nearestPreNeighbourDistanceOLD(self,preType,postType,nMax=1000,
                                  nameStr=""):

    postPop = self.getSubPop(neuronID=self.populations[postType],
                             volumePart="centre")

    if(len(postPop) > nMax):
      postPop = np.random.permutation(postPop)[1:nMax]

    nearestNeighbourDist = dict([])
    for idx in postPop:
      nearestNeighbourDist[idx] = np.inf

    potentialPrePop = self.populations[preType]

    lastSynPair = (np.nan, np.nan)

    for synapses in self.networkLoad.synapseIterator():
      for synRow in synapses:
        preIdx = synRow[0]
        postIdx = synRow[1]
        if((preIdx,postIdx) == lastSynPair):
          # Already did this pair
          continue

        if(postIdx in postPop and preIdx in potentialPrePop):
          d = np.sqrt(np.sum((self.positions[preIdx,:]
                              -self.positions[postIdx,:]) ** 2))
          if(d < nearestNeighbourDist[postIdx]):
            nearestNeighbourDist[postIdx] = d

        lastSynPair = (preIdx,postIdx)


    nnDist = np.array(list(nearestNeighbourDist.values()))
    nnDist = nnDist[nnDist < np.inf] * 1e6
    maxDist = max(nnDist)

    plt.figure()
    matplotlib.rcParams.update({"font.size":22})

    plt.hist(nnDist,np.arange(0,maxDist+25,25))

    plt.xlabel("Distance")
    plt.ylabel("Count")
    plt.title("Nearest presynaptic neighbour " \
              + self.neuronName(preType) + " to " + self.neuronName(postType))

    # Data from Sabatini 2016
    if(preType == "LTS" and (postType == "dSPN" or postType == "iSPN")):
      SabatiniLTS = np.genfromtxt("../DATA/LTS-nearest-neighbour-points-Sabatini2016.csv", delimiter = ",")
      LTSpoints = SabatiniLTS[:,1] * 1e3 # Get in micrometer
      plt.hist(LTSpoints,color='r',histtype="step")

    plt.ion()
    plt.draw()

    if(self.showPlots):
      plt.show()

    plt.pause(0.001)

    figName = "figures/Nearest-presynaptic-neighbour-to-" \
      + str(postType) + "-from-" + str(preType)

    self.saveFigure(plt,figName)


  ############################################################################

  # Inspired by:
  # Nao Chuhma, Kenji F. Tanaka, Rene Hen and Stephen Rayport 2011
  #
  # 10% of a neuron type are marked, fraction of presynaptic neurons
  # out of total population

  def ChuhmaVirtualExperiment(self,taggedType=["dSPN","iSPN"],tagFraction=0.1):

    print("Doing Chuma experiments: " + str(taggedType))

    idx = np.concatenate([self.populations[x] for x in taggedType])
    nTagged = np.round(len(idx)*tagFraction).astype(int)

    taggedNeurons = np.sort(np.random.permutation(idx)[:nTagged])

    # Find all presynaptic neurons
    preIdx = np.where(np.sum(self.connectionMatrix[:,taggedNeurons],
                             axis=1)>0)[0]

    # Ooops, no fun... pretty much all neurons got tagged.

    import pdb
    pdb.set_trace()


  ############################################################################

  # Number of ChINs connected to each MS

  ############################################################################

if __name__ == "__main__":

  assert False, "Do you want to run Network_analyse_striatum.py instead?"

  if(len(sys.argv) > 1):
    hdf5File = sys.argv[1]
    print("Using user supplied HDF5 file: " + hdf5File)

  else:
    hdf5File = None # Auto detect, use latest file
    # hdf5File = "save/network-connect-synapse-file-0.hdf5"

  #import cProfile
  #cProfile.run('na = NetworkAnalyse(hdf5File,loadCache=False,lowMemory=False)')
  #import pdb
  #pdb.set_trace()


  na = NetworkAnalyse(hdf5File,loadCache=True,lowMemory=False,sideLen=250e-6,
                      volumeType="full") # "cube"

  #na = NetworkAnalyse(hdf5File,loadCache=False)

  #na.plotNeurons(neuronID=5,plotPreNeurons=False)
  #na.plotNeurons(neuronID=5)


  enableAnalysis = True #True #False

  # No exp data for this
#  dist3D = False
#  na.plotConnectionProbabilityChannels("FSN","FSN", \
#                                       dist3D=dist3D, \
#                                       expMaxDist=[],\
#                                       expData=[])
#
#  import pdb
#  pdb.set_trace()

  # na.cornerNeurons(sideLen = 100e-6)

  #na.plotSynapseDist(densityFlag=True)
  #na.plotSynapseCumDist()

  #na.virtualAxonSynapses("dSPN")
  #na.virtualAxonSynapses("iSPN")
  #na.virtualAxonSynapses("FSN")

  #na.simpleMotif("dSPN","dSPN","dSPN")
  #na.simpleMotif("iSPN","iSPN","iSPN")
  #na.simpleMotif("dSPN","dSPN","iSPN")
  # na.analyseSingleMotifs("dSPN")


  na.nearestPreNeighbourDistance("LTS","dSPN")
  na.nearestPreNeighbourDistance("LTS","iSPN")


  # na.ChuhmaVirtualExperiment(taggedType=["dSPN","iSPN"])

  dist3D = False

  # 3/21 LTS->MS, Basal Ganglia book --- distance??
  # Ibanez-Sandoval, ..., Tepper  2011 3/21 -- if patching around visual axon
  # but 2/60 when patching blind
  # !!! Use the focused 3/21 statistics for validation!! --- please :)
  na.plotConnectionProbability("LTS","dSPN", \
                                       dist3D=dist3D,
                                       expMaxDist=[250e-6],
                                       expData=[2/60.0],
                                       expDataDetailed=[(2,60)])

  na.plotConnectionProbability("LTS","iSPN", \
                                       dist3D=dist3D,
                                       expMaxDist=[250e-6],
                                       expData=[2/60.0],
                                       expDataDetailed=[(2,60)])


  # Silberberg et al 2013, 2/12 FS-> LTS connected --- distance??
  na.plotConnectionProbability("FSN","LTS", \
                                       dist3D=dist3D,
                                       expMaxDist=[250e-6],
                                       expData=[2.0/12],
                                       expDataDetailed=[(2,12)])

  #import pdb
  #pdb.set_trace()




  na.plotConnectionProbability("dSPN","ChIN", \
                                       dist3D=dist3D)
  na.plotConnectionProbability("iSPN","ChIN", \
                                       dist3D=dist3D)

  na.plotConnectionProbability("ChIN","LTS", \
                                       dist3D=dist3D)
  na.plotConnectionProbability("ChIN","iSPN", \
                                       dist3D=dist3D)
  na.plotConnectionProbability("ChIN","dSPN", \
                                       dist3D=dist3D)

  print("Check the ChIN stuff")
  #import pdb
  #pdb.set_trace()



  #import pdb
  #pdb.set_trace()

  na.nearestPreNeighbourDistance("FSN","dSPN")
  na.nearestPreNeighbourDistance("FSN","iSPN")

  na.plotNumSynapsesPerPair("ChIN","dSPN")
  na.plotNumSynapsesPerPair("ChIN","iSPN")

  na.plotNumSynapsesPerPair("dSPN", "ChIN")
  na.plotNumSynapsesPerPair("iSPN", "ChIN")


  # 2-5 ChIN should connect to each MS (approx)
  na.plotIncomingConnections(neuronType="dSPN",preType="ChIN")
  na.plotIncomingConnections(neuronType="iSPN",preType="ChIN")

  na.plotIncomingConnections(neuronType="ChIN",preType="dSPN")
  na.plotIncomingConnections(neuronType="ChIN",preType="iSPN")

  # LTS plots
  na.plotNumSynapsesPerPair("LTS","dSPN")
  na.plotNumSynapsesPerPair("LTS","iSPN")
  na.plotNumSynapsesPerPair("LTS","ChIN")

  na.plotNumSynapsesPerPair("ChIN","LTS")
  na.plotNumSynapsesPerPair("FSN","LTS")

  na.plotIncomingConnections(neuronType="dSPN",preType="LTS")
  na.plotIncomingConnections(neuronType="iSPN",preType="LTS")
  na.plotIncomingConnections(neuronType="ChIN",preType="LTS")

  na.plotIncomingConnections(neuronType="LTS",preType="ChIN")
  na.plotIncomingConnections(neuronType="LTS",preType="FSN")


  if(True or enableAnalysis):
    dist3D = False

    # 100e-6 from Planert 2010, and 250e-6 data from Gittis 2010
    # 150e-6 from Gittis 2011 (actually 100 +/- 50 micrometers)
    na.plotConnectionProbability("FSN","iSPN", \
                                         dist3D=dist3D, \
                                         expMaxDist=[100e-6, 150e-6, 250e-6],
                                         expData=[6/9.0, 21/54.0, 27/77.0],
                                         expDataDetailed=[(6,9),(21,54),(27,77)])
    na.plotConnectionProbability("FSN","dSPN", \
                                         dist3D=dist3D, \
                                         expMaxDist=[100e-6, 150e-6, 250e-6],
                                         expData=[8/9.0, 29/48.0, 48/90.0],
                                         expDataDetailed=[(8,9),(29,48),(48,90)])
    na.plotConnectionProbability("dSPN","iSPN", \
                                         dist3D=dist3D, \
                                         expMaxDist=[50e-6,100e-6],\
                                         expData=[3/47.0,3/66.0],
                                         expDataDetailed=[(3,47),(3,66)])
    na.plotConnectionProbability("dSPN","dSPN", \
                                         dist3D=dist3D, \
                                         expMaxDist=[50e-6,100e-6],\
                                         expData=[5/19.0,3/43.0],
                                         expDataDetailed=[(5,19),(3,43)])
    na.plotConnectionProbability("iSPN","dSPN", \
                                         dist3D=dist3D, \
                                         expMaxDist=[50e-6,100e-6],\
                                         expData=[13/47.0,10/80.0],
                                         expDataDetailed=[(13,47),(10,80)])
    na.plotConnectionProbability("iSPN","iSPN", \
                                         dist3D=dist3D, \
                                         expMaxDist=[50e-6,100e-6],\
                                         expData=[14/39.0,7/31.0],
                                         expDataDetailed=[(14,39),(7,31)])

    # No exp data for this -- Gittis,...,Kreitzer 2010 (p2228) -- 7/12 (and 3/4 reciprocal) -- distance?
    # FS->FS synapses weaker, 1.1 +/- 1.5nS
    na.plotConnectionProbability("FSN","FSN", \
                                         dist3D=dist3D, \
                                         expMaxDist=[250e-6],\
                                         expData=[7/12.0],
                                         expDataDetailed=[(7,12)])


    # Do we have ChINs?
    # ChIN data, Johanna had ref. ????
    # Janickova

    if(True):
      # REF???!?!?!?!
      #na.plotConnectionProbability("ChIN","iSPN", \
      #                                     dist3D=dist3D,
      #                                     expMaxDist=[200e-6],
      #                                     expData=[62/89.0],
      #                                     expDataDetailed=[(62,89)])
      #na.plotConnectionProbability("ChIN","dSPN", \
      #                                     dist3D=dist3D,
      #                                     expMaxDist=[200e-6],
      #                                     expData=[62/89.0],
      #                                     expDataDetailed=[(62,89)])

      # Derived from Janickova H, ..., Bernard V 2017
      na.plotConnectionProbability("ChIN","iSPN", \
                                           dist3D=dist3D,
                                           expMaxDist=[200e-6],
                                           expData=[0.05])
      na.plotConnectionProbability("ChIN","dSPN", \
                                           dist3D=dist3D,
                                           expMaxDist=[200e-6],
                                           expData=[0.05])


      na.plotConnectionProbability("ChIN","FSN", \
                                           dist3D=dist3D)

      na.plotIncomingConnections(neuronType="dSPN",preType="ChIN")
      na.plotIncomingConnections(neuronType="iSPN",preType="ChIN")


    if(True):

      na.plotConnectionProbability("LTS","ChIN", \
                                           dist3D=dist3D)

      na.plotConnectionProbability("ChIN","LTS", \
                                           dist3D=dist3D)





  if(True):
    print("The synapse dist function needs a density func, which currently not working since we no longer include compartment length in the dend data, so need to calculate it")
    na.plotSynapseDist()
    na.plotSynapseCumDist()
    na.plotSynapseDist(densityFlag=True)

    #import pdb
    #pdb.set_trace()


  if(True and enableAnalysis):
    na.plotNumSynapsesPerPair("FSN","dSPN")
    na.plotNumSynapsesPerPair("FSN","iSPN")
    na.plotNumSynapsesPerPair("dSPN","dSPN")
    na.plotNumSynapsesPerPair("dSPN","iSPN")
    na.plotNumSynapsesPerPair("iSPN","dSPN")
    na.plotNumSynapsesPerPair("iSPN","iSPN")



  if(True and enableAnalysis):
    na.plotIncomingConnections(neuronType="dSPN",preType="iSPN")
    na.plotIncomingConnections(neuronType="dSPN",preType="dSPN")
    na.plotIncomingConnections(neuronType="dSPN",preType="FSN")



  if(True and enableAnalysis):
    na.plotIncomingConnections(neuronType="iSPN",preType="iSPN")
    na.plotIncomingConnections(neuronType="iSPN",preType="dSPN")
    na.plotIncomingConnections(neuronType="iSPN",preType="FSN")


  print("All done, exiting")
  #import pdb
  #pdb.set_trace()
