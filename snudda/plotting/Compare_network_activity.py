# Example:
# Modified code from Johannes Hjorth, by Johanna Frost Nylen 2020-07-10
# python3 Network_plot_spike_raster.py save/traces/network-output-spikes-0.txt save/network-connect-synapse-file-0.hdf5


import sys
import os
import numpy as np
import re
import ntpath
from snudda.load import SnuddaLoad
import time
import elephant
import neo
import pandas as pd
import matplotlib.pyplot as plt
import neo
import quantities as pq

class CompareNetwork(object):

  def __init__(self,fileNames,CompareNeuronType= None,networkFiles=None,skipTime=0.0,typeOrder=None,endTime=2.0):

    self.fileNames = dict()
    self.time = dict()
    self.spikeID = dict()
    self.CompareNeuronType = CompareNeuronType
    
    simulations_i = 0
    
    for fileName in fileNames:
      self.fileNames.update({simulations_i : fileName})
      self.time.update({simulations_i : list()})
      self.spikeID.update({simulations_i : list()})
      simulations_i+=1
      
    self.endTime = endTime

    try:
      
      self.ID = int(re.findall('\d+', ntpath.basename(fileName))[0])
    except:
      self.ID = 0
      
    self.neuronNameRemap = {"FSN" : "FS"}
    
    self.readCSV()

    self.networkInfos = dict()
    self.networkFiles = dict()

    simulations_i = 0
    if(type(networkFiles) is not list):
      networkFile = networkFiles
      for r in range(len(fileNames)):
        self.networkInfos.update({simulations_i: SnuddaLoad(networkFile)})
        self.networkFiles.update({simulations_i: networkFile})
        simulations_i = simulations_i + 1
    elif(type(networkFiles) is list):
      
      for networkFile in networkFiles:
        self.networkInfos.update({simulations_i: SnuddaLoad(networkFile)})
        self.networkFiles.update({simulations_i: networkFile})
        simulations_i = simulations_i + 1
        
      #assert(int(self.ID) == int(self.networkInfo.data["SlurmID"]))
    else:
      self.networkInfos = None
      self.networkFiles = None

    self.plotColourRaster(skipTime=skipTime,typeOrder=typeOrder)

    '''  
    if(self.networkInfos is None):
      print("If you also give network file, then the plot shows neuron types")
      self.plotRaster(skipTime=skipTime)
      time.sleep(1)
    else:

      self.sortTraces()
      self.plotColourRaster(skipTime=skipTime,typeOrder=typeOrder)
      time.sleep(1)
    '''  
  ############################################################################

  def neuronName(self,neuronType):

    if(neuronType in self.neuronNameRemap):
      return self.neuronNameRemap[neuronType]
    else:
      return neuronType
    
  ############################################################################  
  
  def readCSV(self):

    i=0
    
    for fileName in self.fileNames.values():
  
      data = np.genfromtxt(fileName, delimiter='\t')
      self.time[i] = data[:,0]*1e-3
      self.spikeID[i] = data[:,1].astype(int)
      i=i+1
  ############################################################################

  def plotColourRaster(self,skipTime,plotIdx="sort",typeOrder=None):

    colours = {"dSPN" : (77./255,151./255,1.0),\
               "iSPN" : (67./255,55./255,181./255),\
               "FSN" : (6./255,31./255,85./255),\
               "ChIN" : (252./255,102./255,0.0),\
               "LTS" : (150./255,63./255,212./255)}


    
    fig  = plt.figure(figsize=(10,10))
    axis = dict()
    binsize = 10*pq.ms

    k = 1
    
    for neuronNames in self.CompareNeuronType:

      positions = '1' + str(len(self.CompareNeuronType)) + str(k)
      axis.update({neuronNames : fig.add_subplot(eval(positions))})
      k=k+1             

    for ctr,networkInfo in self.networkInfos.items():
      
      time = self.time[ctr]
      spikeID = self.spikeID[ctr]              
      if(plotIdx == "sort"):
        plotIdx,tickPos,tickText,typedict = self.sortTraces(typeOrder=typeOrder,networkInfo=networkInfo,spikeID=spikeID) # add Networkinfo extra
      else:
        tickPos,tickText = None,None
      
      plotLookup = self.makePlotLookup(plotIdx)

      dicttype = typedict

      cellTypes = [n["name"].split("_")[0] \
                     for n in networkInfo.data["neurons"]]

      #cols = [colours[c] for c in cellTypes]

    
      tIdx = np.where(time >= skipTime)[0]

              
      #cols2 = [colours[cellTypes[int(s)]] for s in spikeID]

      cellTypeactivity = dict.fromkeys(self.CompareNeuronType,list())

        
      endTime = np.max([self.endTime,np.max(time)])
      for t in typeOrder:

        if t in cellTypeactivity.keys():
          
          spike_train = list()
          #for i in tIdx:
          for i in range(10):
            if i in dicttype[t]:
              print(i)
              spike_train.append(time[int(i)]-skipTime)
          print(t)
          cellTypeactivity[t].append(neo.SpikeTrain(spike_train*pq.s,t_stop=endTime)) 
      print('added')
      for neuronName, activity in cellTypeactivity.items():
        
        populationCount = elephant.statistics.time_histogram(activity, binsize,output='rate')
        number_series = pd.Series(np.transpose(populationCount)[0])
        windows = number_series.rolling(5)
        moving_averages = windows.mean().dropna()
        
        axis[neuronName].step(np.arange(len(moving_averages)),np.array(moving_averages),label=self.fileNames[ctr])
        axis[neuronName].set_ylim([0,np.max(moving_averages)*1.25])
        axis[neuronName].legend()


    figPath = os.path.dirname(self.networkFiles[ctr]) + "/figs"
    if(not os.path.exists(figPath)):
      os.makedirs(figPath)
    
    # have updated the name of the saved file to be the same as the fileName
    fn = os.path.basename(self.fileNames[ctr])
    figName = '{}/{}{}'.format(figPath, fn.split('.')[0], '-colour.svg')
    print("Saving " + figName)
    plt.savefig(figName,dpi=600)    

   ############################################################################

  def sortTraces(self,networkInfo = None, spikeID = None,typeOrder=None):

    
    print("Sort the traces")

    allTypes = [x["type"] for x in networkInfo.data["neurons"]]
    
    if(typeOrder is None):
      typeOrder = np.unique(allTypes)

    # This one was fun to write. For every type t in typeOrder we want
    # to find all the indexes associated with it (here using enumerate and if)
    idx = [ i for t in typeOrder for i,x in enumerate(allTypes) if x == t]
    
    # new dict with cell type specific spikes
    typedict = {'order':typeOrder}
    for t in typeOrder:
        typedict[t] = [ i for i,x in enumerate(spikeID) if allTypes[x] == t]
   
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
    npsr = CompareNetwork(fileName,networkFile,skipTime=0.0,
                                  endTime=endTime,
                                  typeOrder=["FSN","dSPN","LTS","iSPN","ChIN"])

