# Example:
#
# python3 Network_plot_spike_raster.py save/traces/network-output-spikes-0.txt save/network-connect-synapse-file-0.hdf5


import sys
import os
import numpy as np
import re
import ntpath
from snudda.load import SnuddaLoad
import time

class NetworkPlotSpikeRaster(object):

  def __init__(self,fileName,networkFile=None,skipTime=0.0,typeOrder=None,endTime=2.0):
    
    self.fileName = fileName

    self.time = []
    self.spikeID = []
    self.endTime = endTime

    try:
      self.ID = int(re.findall('\d+', ntpath.basename(fileName))[0])
    except:
      self.ID = 0
      
    self.neuronNameRemap = {"FSN" : "FS"}
    
    self.readCSV()

    if(networkFile is not None):
      self.networkInfo = SnuddaLoad(networkFile)
      self.networkFile = networkFile

      #assert(int(self.ID) == int(self.networkInfo.data["SlurmID"]))
    else:
      self.networkInfo = None
      self.networkFile = None

    if(self.networkInfo is None):
      print("If you also give network file, then the plot shows neuron types")
      self.plotRaster(skipTime=skipTime)
      time.sleep(1)
    else:

      self.sortTraces()
      self.plotColourRaster(skipTime=skipTime,typeOrder=typeOrder)
      time.sleep(1)
      
  ############################################################################

  def neuronName(self,neuronType):

    if(neuronType in self.neuronNameRemap):
      return self.neuronNameRemap[neuronType]
    else:
      return neuronType
    
  ############################################################################  
  
  def readCSV(self):

    data = np.genfromtxt(self.fileName, delimiter='\t')
    self.time = data[:,0]*1e-3
    self.spikeID = data[:,1].astype(int)
    

  ############################################################################
     
  def plotRaster(self,skipTime=0):
    
    import matplotlib.pyplot as plt
    plt.figure()
    plt.scatter(self.time-skipTime,self.spikeID,color='black',s=1)
    plt.xlabel('Time (s)')
    plt.ylabel('Neurons')
    xl = plt.xlim()
    newx = (0,np.max(self.time)-skipTime)
    plt.xlim(newx)
    plt.ion()
    plt.show()
    plt.draw()
    plt.pause(0.001)
    # plt.savefig('figures/Network-spike-raster-' + str(self.ID) + ".pdf")
    plt.savefig('figures/Network-spike-raster-' + str(self.ID) + ".png",dpi=600)    

    print("Figure done")

  ############################################################################

  def plotColourRaster(self,skipTime,plotIdx="sort",typeOrder=None):

    if(plotIdx == "sort"):
      plotIdx,tickPos,tickText,typedict = self.sortTraces(typeOrder)
    else:
      tickPos,tickText = None,None

    plotLookup = self.makePlotLookup(plotIdx)
      
    import matplotlib.pyplot as plt

    cellTypes = [n["name"].split("_")[0].lower() for n in self.networkInfo.data["neurons"]]

    colours = {"dSPN".lower(): (77./255,151./255,1.0),
               "iSPN".lower(): (67./255,55./255,181./255),
               "FS".lower(): (6./255,31./255,85./255),
               "ChIN".lower(): (252./255,102./255,0.0),
               "LTS".lower(): (150./255,63./255,212./255)}

    cols = [colours[c] for c in cellTypes]

    fig  = plt.figure(figsize=(6,4))
    r    = 4
    grid = plt.GridSpec(r, r, hspace=0, wspace=0)
    ax   = fig.add_subplot(grid[1:,:])
    atop = fig.add_subplot(grid[0,:])
    tIdx = np.where(self.time >= skipTime)[0]

    cols2 = [colours[cellTypes[int(s)]] for s in self.spikeID]

    
    ax.scatter(self.time[tIdx]-skipTime,
               plotLookup[self.spikeID[tIdx]],
               color=[cols2[t] for t in tIdx],s=1,
               linewidths=0.1)
    
    # histogram
    for t in typeOrder:
        pruned_spikes = [   self.time[int(i)]-skipTime for i in tIdx if i in typedict[t] ]

        nOfType = len([x["type"] for x in self.networkInfo.data["neurons"] \
                       if x["type"] == t])

        atop.hist(  pruned_spikes, 
                    bins=int(self.time[-1]*100), 
                    range=(0,self.time[-1]), 
                    density=0,
                    color=colours[t], 
                    alpha=1.0, 
                    histtype='step')
    
    ax.invert_yaxis()
    
    ax.set_xlabel('Time (s)')
    if(tickPos is not None):
      ax.set_yticks(tickPos)
      ax.set_yticklabels(tickText)
    else:
      ax.ylabel('Neurons')
    
    # set axes ---------------------------------------------------
    atop.axis('off')
    # UPDATE here to set specific range for plot window!!!!
    
    endTime = np.max([self.endTime,np.max(self.time)]) - skipTime
    
    atop.set_xlim([-0.01,endTime+0.01])
    ax.set_xlim(  [-0.01,endTime+0.01])
    
    m = len(self.networkInfo.data["neurons"])
    offset = m*0.05 # 5%
    ax.set_ylim([-offset,m+offset])
    # -----------------------------------------------------------
    plt.ion()
    plt.show()
    plt.draw()
    plt.pause(0.001)
    # plt.savefig('figures/Network-spike-raster-' + str(self.ID) + "-colour.pdf")
    
    # TODO: this is not working if run from the same folder as the networkFile
    # if so -> figPath = "/figs"
    figPath = os.path.dirname(self.networkFile) + "/figures"
    if(not os.path.exists(figPath)):
      os.makedirs(figPath)
    
    # have updated the name of the saved file to be the same as the fileName
    fn = os.path.basename(self.fileName)
    figName = '{}/{}{}'.format(figPath, fn.split('.')[0], '-colour.pdf')
    print("Saving " + figName)
    plt.savefig(figName,dpi=600)    

   ############################################################################

  def sortTraces(self,typeOrder=None):

    print("Sort the traces")

    allTypes = [x["type"].lower() for x in self.networkInfo.data["neurons"]]

    if(typeOrder is None):
      typeOrder = np.unique(allTypes)

    # This one was fun to write. For every type t in typeOrder we want
    # to find all the indexes associated with it (here using enumerate and if)
    idx = [ i for t in typeOrder for i,x in enumerate(allTypes) if x == t]
    
    # new dict with cell type specific spikes
    typedict = {'order':typeOrder}
    for t in typeOrder:
        typedict[t] = [ i for i,x in enumerate(self.spikeID) if allTypes[x] == t]
    
    prevPos = 0
    tickPos = []
    tickText = []
    
    for t in typeOrder:
      nElement = np.sum([x == t for x in allTypes])
      
      if(nElement == 0): # No labels for missing types
        continue
      
      tickPos.append(prevPos + nElement/2)
      tickText.append(self.neuronName(t))
      prevPos += nElement
      
    return idx,tickPos,tickText,typedict

  ############################################################################

  def makePlotLookup(self,plotIdx):

    plotLookup = np.nan*np.zeros(len(plotIdx))

    for i,p in enumerate(plotIdx):
      plotLookup[p] = i

    return plotLookup

  ############################################################################


    
  
############################################################################
    
if __name__ == "__main__":

  print("Usage: " + sys.argv[0] + " network-output-spikes-XXX.txt")
    
  if(len(sys.argv) > 1):
    fileName = sys.argv[1]
  else:
    fileName = None
  
  if(len(sys.argv) > 2):
    networkFile = sys.argv[2]
  else:
    networkFile = None

  if(len(sys.argv) > 3):
    endTime = float(sys.argv[3])
  else:
    endTime = 2.0
    
  if(fileName is not None):
    #type_order = ["FS", "dSPN", "LTS", "iSPN", "ChIN"]
    type_order = ["fs", "dspn", "lts", "ispn", "chin"]

    npsr = NetworkPlotSpikeRaster(fileName,networkFile,skipTime=0.0,
                                  endTime=endTime,
                                  typeOrder=type_order)


  #import pdb
  #pdb.set_trace()
