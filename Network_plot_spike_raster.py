# Example:
#
# python3 Network_plot_spike_raster.py save/traces/network-output-spikes-0.txt save/network-connect-synapse-file-0.hdf5


import sys
import os
import numpy as np
import re
import ntpath
from snudda_load import SnuddaLoad

class NetworkPlotSpikeRaster(object):

  def __init__(self,fileName,networkFile=None,skipTime=0.0,typeOrder=None):
    
    self.fileName = fileName

    self.time = []
    self.spikeID = []

    self.ID = int(re.findall('\d+', ntpath.basename(fileName))[0])

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
    else:

      self.sortTraces()
      self.plotColourRaster(skipTime=skipTime,typeOrder=typeOrder)

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
    newx = (0,xl[1])
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
      plotIdx,tickPos,tickText = self.sortTraces(typeOrder)
    else:
      tickPos,tickText = None,None

    plotLookup = self.makePlotLookup(plotIdx)
      
    import matplotlib.pyplot as plt

    cellTypes = [n["name"].split("_")[0] \
                 for n in self.networkInfo.data["neurons"]]

    colours = {"dSPN" : (77./255,151./255,1.0),
               "iSPN" : (67./255,55./255,181./255),
               "FSN" : (6./255,31./255,85./255),
               "ChIN" : (252./255,102./255,0.0),
               "LTS" : (150./255,63./255,212./255)}

    cols = [colours[c] for c in cellTypes]

    #import pdb
    #pdb.set_trace()


    
    fig=plt.figure()
    ax=fig.add_subplot(111)
#    for (t,sid) in zip(self.time,self.spikeID):
#      plt.scatter(t,sid,color=cols[int(sid)],s=1)
    tIdx = np.where(self.time >= skipTime)[0]


    cols2 = [colours[cellTypes[int(s)]] for s in self.spikeID]

    ax.scatter(self.time[tIdx]-skipTime,
               plotLookup[self.spikeID[tIdx]],
               color=[cols2[t] for t in tIdx],s=1)
    #ax.xlim(xmin=0)
    
    ax.set_xlabel('Time (s)')
    if(tickPos is not None):
      ax.set_yticks(tickPos)
      ax.set_yticklabels(tickText)
    else:
      ax.ylabel('Neurons')
      
    plt.ion()
    plt.show()
    plt.draw()
    plt.pause(0.001)
    # plt.savefig('figures/Network-spike-raster-' + str(self.ID) + "-colour.pdf")

    figPath = os.path.dirname(self.networkFile) + "/figs/"
    if(not os.path.exists(figPath)):
      os.makedirs(figPath)
    figName = figPath + 'Network-spike-raster-' + str(self.ID) + "-colour.png"
    print("Saving " + figName)
    plt.savefig(figName,dpi=600)    

   ############################################################################

  def sortTraces(self,typeOrder=None):

    print("Sort the traces")

    allTypes = [x["type"] for x in self.networkInfo.data["neurons"]]

    if(typeOrder is None):
      typeOrder = np.unique(allTypes)

    # This one was fun to write. For every type t in typeOrder we want
    # to find all the indexes associated with it (here using enumerate and if)
    idx = [ i for t in typeOrder for i,x in enumerate(allTypes) if x == t]

    prevPos = 0
    tickPos = []
    tickText = []
    
    for t in typeOrder:
      nElement = np.sum([x == t for x in allTypes])
      
      if(nElement == 0): # No labels for missing types
        continue
      
      tickPos.append(prevPos + nElement/2)
      tickText.append(t)
      prevPos += nElement
      
    return idx,tickPos,tickText

  ############################################################################

  def makePlotLookup(self,plotIdx):

    plotLookup = np.nan*np.zeros((np.max(self.spikeID)+1,))

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

  if(fileName is not None):
    npsr = NetworkPlotSpikeRaster(fileName,networkFile,skipTime=0.5,
                                  typeOrder=["FSN","dSPN","LTS","iSPN","ChIN"])


  import pdb
  pdb.set_trace()
