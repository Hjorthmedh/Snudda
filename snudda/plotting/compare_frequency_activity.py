

# python3 Network_plot_traces.py save/traces/network-voltage-0.csv save/network-connect-synapse-file-0.hdf5

# Modified code from Johannes Hjorth, by Johanna Frost Nylen 2020-07-12


import sys
import os
import numpy as np
from snudda.load import SnuddaLoad
import re
import ntpath
import time
import elephant
import neo
import quantities as pq
import pandas as pd

class ComparePlotTraces():

  ############################################################################
  
  def __init__(self,fileNames,networkFiles=None,labels=None,colours=None):

    self.fileNames = dict()
    self.time = dict()
    self.spikeID = dict()
    self.labels = labels
    self.networkFiles = networkFiles
    self.time = dict()
    self.voltage = dict()
    self.colours = colours
    self.neuronNameRemap = {"FSN" : "FS"}

    simulations_i = 0
    
    for fileName in fileNames:
      print(simulations_i)
      self.fileNames.update({simulations_i : fileName})
      self.time.update({simulations_i : list()})
      self.spikeID.update({simulations_i : list()})
      simulations_i = simulations_i + 1
      
    self.readCSV()
    print('read')
    self.networkInfos = dict()
    try:
      self.ID = int(re.findall('\d+', ntpath.basename(fileName))[0])
    except:
      print("Unable to guess ID, using 666.")
      self.ID = 666
    self.networkInfos = dict()
    simulations_i = 0
    if(type(networkFiles) is not list):
      networkFile = networkFiles
      for r in range(len(fileNames)):
        self.networkInfos.update({simulations_i: SnuddaLoad(networkFile)})
        simulations_i = simulations_i + 1
    elif(type(networkFiles) is list):
      
      for networkFile in networkFiles:
        self.networkInfos.update({simulations_i: SnuddaLoad(networkFile)})
        simulations_i = simulations_i + 1

    else:
      self.networkInfo = None
  ############################################################################
    
  def readCSV(self):

    i = 0

    for fileName in self.fileNames.values():
      print(fileName)
      
      data = np.genfromtxt(fileName, delimiter=',')
      assert(data[0,0] == -1) # First column should be time

      self.time[i] = data[0,1:] / 1e3

      self.voltage[i] = dict()
      print(len(data[1:,:]))
      for rows in data[1:,:]:
        cID = int(rows[0])
        self.voltage[i][cID] = rows[1:] * 1e-3
      i=i+1      
  ############################################################################

  def neuronName(self,neuronType):

    if(neuronType in self.neuronNameRemap):
      return self.neuronNameRemap[neuronType]
    else:
      return neuronType

  ############################################################################
  
  
  def plotTraces(self,traceID=None,offset=150e-3,colours=None,skipTime=None,
                 title=None):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    chosenColor = self.colours
    if(skipTime is not None):
      print("!!! Excluding first " + str(skipTime) + "s from the plot")
    
    if(colours is None):
      colours = {"dSPN" : (77./255,151./255,1.0),
                 "iSPN" : (67./255,55./255,181./255),
                 "FSN" : (6./255,31./255,85./255),
                 "ChIN" : [252./255,102./255,0],
                 "LTS" : [150./255,63./255,212./255]}
    for ctr,networkInfo in self.networkInfos.items():
      print("Plotting traces: " + str(traceID))
      
      print("Plotted " + str(len(traceID)) + " traces (total " + str(len(self.voltage[ctr])) + ")")
      
      

      typesInPlot = set()
    
      if(networkInfo is not None):
        cellTypes = [n["type"] for n in networkInfo.data["neurons"]]
        cellIDcheck = [n["neuronID"] for n in networkInfo.data["neurons"]]
        try:
          assert (np.array([cellIDcheck[x] == x for x in traceID])).all(), \
            "Internal error, assume IDs ordered"
        except:
          import traceback
          tstr = traceback.format_exc()
          print(tstr)
          print("This is strange...")
          import pdb
          pdb.set_trace()
      
        cols = [colours[c] if c in colours else [0,0,0] for c in cellTypes]
    
        #import pdb
        #pdb.set_trace()
    
        
    
        ofs = 0

        if(skipTime is not None):
          timeIdx = np.where(self.time[ctr] >= skipTime)[0]



        else:
          skipTime = 0.0
          timeIdx = range(0,len(self.time))

        plotCount = 0

        traceChosen = list()
        
        for r in traceID:

          if(r not in self.voltage[ctr]):
            print("Missing data for trace " + str(r))
            continue


          plotCount += 1
          typesInPlot.add(networkInfo.data["neurons"][r]["type"])
      
          if(colours is None or networkInfo is None):
            colour = "black"
          else:
            try:
              colour = cols[r]
            except:
              import traceback
              tstr = traceback.format_exc()
              print(tstr)
              import pdb
              pdb.set_trace()
          
          plt.figure(0)
          plt.plot(self.time[ctr][timeIdx]-skipTime,self.voltage[ctr][r][timeIdx] + ofs,color=chosenColor[ctr],label=self.labels[ctr])
          
          traceChosen.append(elephant.spike_train_generation.threshold_detection(neo.AnalogSignal(self.voltage[ctr][r][timeIdx] + ofs, units='mV',sampling_period = 0.5 * pq.ms),threshold = 0 *pq.mV))

          ofs += offset
        import pdb
        pdb.set_trace()
        binsize = 10 * pq.ms 
        populationCount = elephant.statistics.time_histogram(traceChosen, binsize,output='rate')
        number_series = pd.Series(np.transpose(populationCount)[0])
        windows = number_series.rolling(5)
        moving_averages = windows.mean().dropna()
        plt.figure(1)
        plt.step(np.arange(len(moving_averages)),np.array(moving_averages),label=self.fileNames[ctr])
        #plt.ylim([0,np.max(moving_averages)*1.25])
        plt.legend()

    import pdb
    pdb.set_trace()
      
    plt.xlabel('Time')
    plt.ylabel('Voltage')

    if(title is not None):
      plt.title(title)
    
    if(offset != 0):
      ax = fig.axes[0]
      ax.set_yticklabels([])
    #plt.legend()
    plt.tight_layout()
    #plt.show()
   

    #plt.savefig('figures/Network-spikes-' + str(self.ID) + "-colour.pdf")

    figPath = os.path.dirname(self.networkFiles) + "/figs"
    if(not os.path.exists(figPath)):
      os.makedirs(figPath)
 
    
    
    if(len(typesInPlot) > 1):
      figName = figPath + '/Network-spikes-' + str(self.ID) \
        + "-".join(typesInPlot) + "-colour.svg"
    else:
      figName = figPath + '/Network-spikes-' + str(self.ID) \
        + "-" + typesInPlot.pop() + "-colour.svg"
      
    plt.savefig(figName,
                dpi=300)
    print("Saving to figure " + str(figName))


  ############################################################################

  def plotTraceNeuronType(self,neuronType,nTraces=10,offset=0,skipTime=0.0):
    
    assert self.networkInfos is not None, "You need to specify networkInfo file"

    for dictIndx,networkInfo in self.networkInfos.items():
      neuronTypes = [x["type"] for x in networkInfo.data["neurons"]]

      # Find numbers of the relevant neurons
    
      traceID = [x[0] for x in enumerate(neuronTypes) if x[1] == neuronType]
    
      nTraces = min(len(traceID),nTraces)

      if(nTraces <= 0):
        print("No traces of " + str(neuronType) + " to show")
        return
    
      self.plotTraces(offset=offset,traceID=traceID[:nTraces],skipTime=skipTime,
                    title=self.neuronName(neuronType))
                                   
    time.sleep(1)
    
  ############################################################################
    
if __name__ == "__main__":

  if(len(sys.argv) > 1):
    fileName = sys.argv[1]

  if(len(sys.argv) > 2):
    networkFile = sys.argv[2]
  else:
    networkFile = None

  if(fileName is not None):
    npt = NetworkPlotTraces(fileName,networkFile)
    #npt.plotTraces(offset=0.2,traceID=[17,40,11,18])
    # npt.plotTraces(offset=0.2,traceID=[0,1,2,3,4,5,6,7,8,9])
    # npt.plotTraces(offset=0,traceID=[0,1,2,3,4,5,6,7,8,9])
    # npt.plotTraces(offset=0,traceID=[0,5,55]) #,5,55])
    #npt.plotTraces(offset=0,traceID=np.arange(0,100)) #,5,55])
    #npt.plotTraces(offset=-0.2,traceID=np.arange(0,20),skipTime=0.5)
    #npt.plotTraces(offset=-0.2,traceID=[5,54],skipTime=0.5)    
    #npt.plotTraces(offset=0.2,traceID=[1,5,7,15,16],skipTime=0.2)

    plotOffset = 0 # -0.2
    skipTime = 0 #0.5
    nTracesMax = 5
    
    npt.plotTraceNeuronType(neuronType="dSPN",nTraces=nTracesMax,offset=plotOffset,skipTime=skipTime)
    npt.plotTraceNeuronType(neuronType="iSPN",nTraces=nTracesMax,offset=plotOffset,skipTime=skipTime)
    npt.plotTraceNeuronType(neuronType="FSN",nTraces=nTracesMax,offset=plotOffset,skipTime=skipTime)
    npt.plotTraceNeuronType(neuronType="LTS",nTraces=nTracesMax,offset=plotOffset,skipTime=skipTime)
    npt.plotTraceNeuronType(neuronType="ChIN",nTraces=nTracesMax,offset=plotOffset,skipTime=skipTime)
    
    
    
  else:
    print("Usage: " + sys.argv[0] + " network-voltage-XXX.csv")
    
