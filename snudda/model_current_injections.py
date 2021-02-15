# In Chuhma et al 2011 they stimulate 10% of MS population and record the
# response in MS, ChIN and FS.
#
# Let us do the same thing but in our model network.
#
# Before you run this network you need to create a network with reduced MS
# density, no point in simulating the silent 90%.
#
# Suggestion: create a slice that is 0.6 x 0.6 mm in length, and 0.15 mm deep.
#
#
# Example usage:
#
# python3 snudda_model_current_injections.py setup Chuhma2011 networks/Chuhma2011-v15 
#
# mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_model_current_injections.py run Chuhma2011 networks/Chuhma2011-v15/ 
#
# python3 snudda_model_current_injections.py analyse Chuhma2011 networks/Chuhma2011-v15/ 
#
# OR
#
#
# python3 snudda_model_current_injections.py setup Straub2016FS networks/Straub2016FS-v9
# mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_model_current_injections.py run Straub2016FS networks/Straub2016FS-v9
# python3 snudda_model_current_injections.py analyse Straub2016FS networks/Straub2016FS-v9
#
#

import os

from snudda.simulate import SnuddaSimulate
from snudda.load import SnuddaLoad
from snudda.init import SnuddaInit
from snudda import Snudda

import numpy as np

import matplotlib
import matplotlib.pyplot as plt

import neuron

class SnuddaModelCurrentInjections(object):

  def __init__(self,
               simName,
               simType,
               curInj = 10e-9):

    self.simName = simName
    self.simType = simType
    self.snuddaSim = None
    self.expDataDict = dict()
    self.expTraceDict = dict()
    
    if(simType == "Chuhma2011"):
      self.tInj = 0.3
      self.injDuration = 1e-3
      self.curInj = 10e-9
      self.tWindow = 0.03
      self.simEnd = self.tInj + self.tWindow*2
      self.holdV = -60e-3
      self.GABArev = 2e-3 # 144mM inside, 133.5mM outside, 32-33C --NERNST--> 2mV
      self.nNrns = 100 # How many we measure from of each type
      
    elif(simType == "Straub2016FS" or simType == "Straub2016LTS"):
      self.tInj = 0.5
      self.injDuration = 1e-3
      self.curInj = 10e-9
      self.tWindow = 0.03
      self.simEnd = self.tInj + self.tWindow*2
      self.holdV = -70e-3
      self.GABArev = -0.3e-3 # Out: 133.5 mM chloride, In 131.5 mM, Temperature 33-34 C
      self.nNrns = 30
    elif(simType == "Szydlowski2013"):

      assert False, "Szydlowski did not stimulate all FS or did they when recording LTS?"
      self.tInj = 0.5
      self.injDuration = 1e-3
      self.curInj = 10e-9
      self.tWindow = 0.03
      self.simEnd = self.tInj + self.tWindow*2
      self.holdV = -70e-3
      self.GABArev = -30e-3 
      self.nNrns = 30
      
    else:
      print("Unknown simType: " + str(simType))
      
    self.plotExpTrace = False
    self.neuronNameRemap = {"FSN" : "FS"}

  ############################################################################

  def neuronName(self,neuronType):

    if(neuronType in self.neuronNameRemap):
      return self.neuronNameRemap[neuronType]
    else:
      return neuronType

  ############################################################################
  
  def defineNetwork(self,simName,simType=None):

    if(simType is None):
      simType = self.simType
            
    configName= simName + "/network-config.json"
    cnc = SnuddaInit(struct_def={}, config_file=configName, nChannels=1)

    # In a 1x1x0.15 mm slice there are 12000 neurons normally
    # We want 10% of MS population only, since those are the ones being
    # stimulated (Chuhma et al 2011)
    #
    # 47.5% dSPN normally, now 4.75% of normal density = 570 dSPN, 570 iSPN
    # We assume we measure only monosynaptic connections.
    #
    # Adding 10 FSN, 10 ChIN, 10 dSPN, 10 iSPN to measure from
    #


    if(False):
      #Small debug version
      #cnc.defineStriatum(nMSD1=20,nMSD2=20,nFS=0,nLTS=0,nChIN=0,
      #                   volumeType="slice",sideLen=200e-6)
      #cnc.defineStriatum(nMSD1=20,nMSD2=20,nFS=10,nLTS=0,nChIN=10,
      #                   volumeType="slice",sideLen=200e-6)
      cnc.define_striatum(num_dSPN=153, num_iSPN=153, num_FS=10, num_LTS=0, num_ChIN=10,
                          volume_type="slice", side_len=500e-6)



    if(simType == "Chuhma2011"):
      cnc.define_striatum(num_dSPN=1140 + self.nNrns,
                          num_iSPN=1140 + self.nNrns,
                          num_FS=5,
                          num_LTS=0,
                          num_ChIN=self.nNrns,
                          volume_type="slice",
                          side_len=1000e-6,
                          slice_depth=300e-6) # 400mum, assume 100 mum dead

    elif(simType == "Straub2016FS"):
      # nFS must be correct density, but readout neurons can be any density
      cnc.define_striatum(num_dSPN=self.nNrns,
                          num_iSPN=self.nNrns,
                          num_FS=182, num_LTS=0,
                          num_ChIN=self.nNrns,
                          volume_type="slice",
                          side_len=1000e-6,
                          slice_depth=175e-6)  #275e-6 m slice, assume 100e-6 dead

    elif(simType == "Straub2016LTS"):
      cnc.define_striatum(num_dSPN=self.nNrns,
                          num_iSPN=self.nNrns,
                          num_FS=0, num_LTS=98,
                          num_ChIN=self.nNrns,
                          volume_type="slice",
                          side_len=1000e-6,
                          slice_depth=175e-6)
    elif(simType == "Szydlowski2013"):
      cnc.define_striatum(num_dSPN=0,
                          num_iSPN=0,
                          num_FS=156,
                          num_LTS=self.nNrns,
                          num_ChIN=0,
                          volume_type="slice",
                          side_len=1000e-6,
                          slice_depth=150e-6)
    else:
      print("setup : Unkown simType: " + str(simType))
      exit(-1)
      
    dirName = os.path.dirname(configName)
  
    if not os.path.exists(dirName):
      os.makedirs(dirName)

    cnc.write_json(configName)

  ############################################################################

  def simulateNetwork(self,simName,simType=None):

    if(simType is None):
      simType = self.simType

    if(simType == "Chuhma2011"):
      self.simulateNetworkChuhma2011(simName)

    elif(simType == "Straub2016FS" or simType == "Straub2016LTS"):
      self.simulateNetworkStraub2016(simName,simType)

    elif(simType == "Szydlowski2013"):
      self.simulateNetworkSzydlowski2013(simName)
      
    else:
      print("simulateNetwork: unknown simType = " + str(simType))
      exit(-1)
      
  ############################################################################

  def simulateNetworkChuhma2011(self,simName):

    simType = "Chuhma2011"

    if(self.snuddaSim is None):
      logFile = simName + "/log/simlog.txt"

      cutFile = simName + "/network-cut-slice.hdf5"
      if(os.path.exists(cutFile)):
        self.networkFile = cutFile
      else:
        self.networkFile = simName + "/network-pruned-synapses.hdf5"

      print("Using network file: " + str(self.networkFile))
      
      self.snuddaSim = SnuddaSimulate(network_file=self.networkFile,
                                      input_file=None,
                                      log_file=logFile,
                                      disable_gap_junctions=True)
    

    # Get neuronID of neurons that will get artificial stimulation
    stimID = [x["neuronID"] \
                   for x in self.snuddaSim.network_info["neurons"] \
                   if "SPN" in x["type"]]


    # Pick neurons that we will measure from
    measureFSN = [x["neuronID"] \
                  for x in self.snuddaSim.network_info["neurons"] \
                  if x["type"] == "FSN"]
    
    measureChIN = [x["neuronID"] \
                   for x in self.snuddaSim.network_info["neurons"] \
                   if x["type"] == "ChIN"]

    # For future: maybe pick the ones more centrally
    measuredSPN = [x["neuronID"] \
                   for x in self.snuddaSim.network_info["neurons"] \
                   if x["type"] == "dSPN"][:self.nNrns]
    measureiSPN = [x["neuronID"] \
                   for x in self.snuddaSim.network_info["neurons"] \
                   if x["type"] == "iSPN"][:self.nNrns]

    
    # Remove the overlap, ie dont stimulate the neurons we measure from    
    stimID = np.setdiff1d(stimID,
                          np.union1d(measuredSPN,measureiSPN))

    measureID = np.union1d(np.union1d(measuredSPN,measureiSPN),
                           np.union1d(measureFSN,measureChIN))

    self._simulateNetworkHelper(simName,simType,stimID,measureID)
    

  ############################################################################

  def simulateNetworkStraub2016(self,simName,simType):

    if(self.snuddaSim is None):
      logFile = simName + "/log/simlog.txt"
      self.networkFile = simName + "/network-pruned-synapses.hdf5"
      
      self.snuddaSim = SnuddaSimulate(network_file=self.networkFile,
                                      input_file=None,
                                      log_file=logFile,
                                      disable_gap_junctions=True)
    
    # Get neuronID of neurons that will get artificial stimulation
    if(simType == "Straub2016FS"):
      stimID = [x["neuronID"] \
                for x in self.snuddaSim.network_info["neurons"] \
                if "FSN" in x["type"]]
    elif(simType == "Straub2016LTS"):
      stimID = [x["neuronID"] \
                for x in self.snuddaSim.network_info["neurons"] \
                if "LTS" in x["type"]]
    else:
      print("simulateNetworkStraub2016: Unknown simType : " + simType)
      exit(-1)

    measureID = [x["neuronID"] \
                 for x in self.snuddaSim.network_info["neurons"] \
                 if x["type"] == "ChIN" \
                 or x["type"] == "dSPN" \
                 or x["type"] == "iSPN"]
      
    self._simulateNetworkHelper(simName,simType,stimID,measureID)

  ############################################################################

  def simulateNetworkSzydlowski2013(self,simName):

    if(self.snuddaSim is None):
      logFile = simName + "/log/simlog.txt"
      self.networkFile = simName + "/network-pruned-synapses.hdf5"
      
      self.snuddaSim = SnuddaSimulate(network_file=self.networkFile,
                                      input_file=None,
                                      log_file=logFile,
                                      disable_gap_junctions=True)
    
    stimID = [x["neuronID"] \
              for x in self.snuddaSim.network_info["neurons"] \
              if "FSN" in x["type"]]

    measureID = [x["neuronID"] \
                 for x in self.snuddaSim.network_info["neurons"] \
                 if x["type"] == "LTS"]
      
    self._simulateNetworkHelper(simName,simType,stimID,measureID)
    
  ############################################################################
  
  def _simulateNetworkHelper(self,simName,simType,stimID,measureID):
    

    if(self.snuddaSim is None):
      logFile = simName + "/log/simlog.txt"
      self.networkFile = simName + "/network-pruned-synapses.hdf5"
          
      self.snuddaSim = SnuddaSimulate(network_file=self.networkFile,
                                      input_file=None,
                                      log_file=logFile,
                                      disable_gap_junctions=True)
    
    # Set up stimulation protocol
    for nID in stimID:
      self.snuddaSim.add_current_injection(neuron_id=nID,
                                           start_time=self.tInj,
                                           end_time=self.tInj + self.injDuration,
                                           amplitude=self.curInj)

    # Add recordings
    self.snuddaSim.add_voltage_clamp(cell_id= measureID,
                                     voltage = self.holdV,
                                     duration=self.simEnd,
                                     save_i_flag=True)
    
    # Also add voltage recording for debugging reasons
    saveVoltage = True # False #True
    if(saveVoltage):
      self.snuddaSim.add_recording(cell_id=stimID)

    self.setGABArev(self.GABArev)
    
    self.snuddaSim.run(self.simEnd*1e3)

    self.currentFile = simName + "/" + simType \
      + "-network-stimulation-current.txt"
    
    self.snuddaSim.write_current(self.currentFile)

    if(saveVoltage):
      voltageFile = simName  + "/" + simType \
        + "-network-stimulation-voltage.txt"
      
      self.snuddaSim.write_voltage(voltageFile)
      
  ############################################################################

  def createNetwork(self,simName):

    sn = Snudda(simName)

    class FakeArgs(object):
      def __init__(self):
        setattr(self,"h5legacy","latest")
        setattr(self, "volumeID", None)
        setattr(self, "hvsize", None) # Use default value
        setattr(self, "cont", None)
        setattr(self, "mergeonly", False)

      
    args = FakeArgs()
    
    
    sn.place_neurons(args)
    sn.touch_detection(args)
    sn.prune_synapses(args)

  ############################################################################

  def setupExpDataDict(self):

    self.expDataDict = dict()

    # Straub et al 2016
    LTS2SPN = np.array([0.0316, 0.0433, 0.0474, 0.1253, 0.1839,
                        0.1860, 0.1946, 0.1968, 0.2082, 0.2203,
                        0.2384, 0.2439, 0.2793, 0.3091, 0.3234,
                        0.3271, 0.3383, 0.3500, 0.3540, 0.3662,
                        0.3662, 0.3831, 0.4053, 0.4099, 0.4288,
                        0.4337, 0.4966, 0.5023, 0.5196, 0.5314,
                        0.5436, 0.5560, 0.5817, 0.6017, 0.6736,
                        0.6968, 0.7047, 0.7047, 0.7127, 0.7979,
                        0.9034, 1.0461])

    LTS2ChIN = np.array([0.2466, 0.5080, 0.5196, 0.6017, 0.6660,
                         0.7541, 0.7713, 0.8442, 1.1069, 1.2391,
                         1.2818, 1.4030, 2.3315])

    FSN2SPN = np.array([0.3091, 0.5137, 0.5255, 0.5687, 0.6890,
                        0.8161, 0.8832, 0.8932, 0.9667, 1.0228,
                        1.0228, 1.0822, 1.1844, 1.2391, 1.2964,
                        1.3111, 1.4189, 1.4350, 1.5530, 1.6247,
                        1.7385, 1.7984, 1.9028, 2.1063, 2.2539,
                        2.3580, 2.4669, 2.4949, 2.5232, 2.7307,
                        2.7930, 2.8247, 3.2711, 3.3458, 3.4222,
                        4.2648, 4.4617, 4.9668, 5.3148])

    FSN2ChIN = np.array([0.0233, 0.0378, 0.0419, 0.0428, 0.0666,
                         0.0762])
    
    self.expDataDict[("Straub2016LTS","dSPN")] = LTS2SPN
    self.expDataDict[("Straub2016LTS","iSPN")] = LTS2SPN
    self.expDataDict[("Straub2016LTS","ChIN")] = LTS2ChIN        
    self.expDataDict[("Straub2016FS","dSPN")]  = FSN2SPN
    self.expDataDict[("Straub2016FS","iSPN")]  = FSN2SPN
    self.expDataDict[("Straub2016FS","ChIN")]  = FSN2ChIN

    if(self.plotExpTrace):
      self.expTraceDict = dict()

      LTS2SPN = np.genfromtxt("DATA/Straub2016/LTSItoSPN_Straub2.txt")
      LTS2ChIN = np.genfromtxt("DATA/Straub2016/LTSItoChIN_Straub2.txt")
      FS2SPN = np.genfromtxt("DATA/Straub2016/FSItoSPN_Straub2_shorter.txt")

      # Convert current from pA to nA
      LTS2SPN[:,1:] = 1e-3*LTS2SPN[:,1:]
      LTS2ChIN[:,1:] = 1e-3*LTS2ChIN[:,1:]
      FS2SPN[:,1:] = 1e-3*FS2SPN[:,1:] 
    
      self.expTraceDict[("Straub2016LTS","dSPN")] = LTS2SPN
      self.expTraceDict[("Straub2016LTS","iSPN")] = LTS2SPN    
      self.expTraceDict[("Straub2016LTS","ChIN")] = LTS2ChIN
      self.expTraceDict[("Straub2016FS","dSPN")] = FS2SPN
      self.expTraceDict[("Straub2016FS","iSPN")] = FS2SPN    

      SPN2SPN = np.genfromtxt("DATA/Chuhma2011/SPNtoSPN_Chuhma.txt")
      SPN2ChIN = np.genfromtxt("DATA/Chuhma2011/SPNtoChIN_Chuhma.txt")

      # Convert current from pA to nA
      SPN2SPN[:,1:] = 1e-3 * SPN2SPN[:,1:]
      SPN2ChIN[:,1:] = 1e-3 * SPN2ChIN[:,1:]    
        
      self.expTraceDict[("Chuhma2011","dSPN")] = SPN2SPN
      self.expTraceDict[("Chuhma2011","iSPN")] = SPN2SPN    
      self.expTraceDict[("Chuhma2011","ChIN")] = SPN2ChIN
    
    
  ############################################################################
  
  def analyseNetwork(self,simName,simType=None,nPlotMax=10):

    figDir = simName + "/figures/"
    if(not os.path.exists(figDir)):
      os.makedirs(figDir)

    if(simType is None):
      simType = self.simType

    if(simType == "Straub2016LTS"):
      preType = "LTS"
      self.setupExpDataDict()
    elif(simType == "Straub2016FS"):
      preType = "FSN"
      self.setupExpDataDict()
    elif(simType == "Chuhma2011"):
      preType = "SPN"
      self.setupExpDataDict()
    elif(simType == "Szydlowski2013"):
      preType = "FSN"
    else:
      print("Unknown simType : " + simType)
      exit(-1)
      
    print("Analysing data in " + simName)
    voltFile = simName + "/" + simType + "-network-stimulation-current.txt"

    # Read data from file
    data = np.genfromtxt(voltFile,delimiter=",")

    assert(data[0,0] == -1) # First column should be time
    time = data[0,1:] * 1e-3
    
    current = dict()
    
    for rows in data[1:,:]:
      cID = int(rows[0])
      current[cID] = rows[1:] * 1e-9

    # Data in time, current now

    # Read the network info
    networkFile = simName + "/network-pruned-synapses.hdf5"    
    self.snuddaLoad = SnuddaLoad(networkFile)
    self.data = self.snuddaLoad.data
    
    recordedNeurons = [x for x in current]

    # Group the neurons by type

    if(simType == "Chuhma2011"):
      neuronTypeList = ["dSPN","iSPN","FSN","ChIN"]
    elif(simType == "Straub2016FS" or simType == "Straub2016LTS"):
      neuronTypeList = ["dSPN","iSPN","ChIN"]
    elif(simType == "Szydlowski2013"):
      neuronTypeList = ["LTS"]
    else:
      print("simulate: Unknown simType: " + simType)
      exit(-1)
      
    neuronPlotList = []

    minTimeIdx = np.where(time > self.tInj)[0][0]
    maxTimeIdx = np.where(time > self.tInj + self.tWindow)[0][0]
    
    for nt in neuronTypeList:
      IDList = [x for x in current if self.data["neurons"][x]["type"] == nt]
      maxIdx = [np.argmax(np.abs(current[x][minTimeIdx:maxTimeIdx] \
                                 -current[x][minTimeIdx])) + minTimeIdx \
                for x in IDList]

      neuronPlotList.append((IDList,maxIdx))
      

    matplotlib.rcParams.update({'font.size': 22})    
    
    for plotID,maxIdx in neuronPlotList:

      if(len(plotID) == 0):
        continue
      
      plotType = self.data["neurons"][plotID[0]]["type"]
      figName = figDir + "/" + simType + "-" + plotType + "-current-traces.pdf"
      figNameHist = figDir + "/" + simType + "-" + plotType + "-current-histogram.pdf"

      goodMax = []
      
      plt.figure()

      peakAmp = []
      peakTime = []
      voltCurve = []
      
      for pID,mIdx in zip(plotID,maxIdx):
        tIdx = np.where(np.logical_and(time > self.tInj,
                                       time < self.tInj+self.tWindow))[0]

        curAmp = current[pID][tIdx]-current[pID][tIdx[0]-1]
        maxAmp = current[pID][mIdx]-current[pID][tIdx[0]-1]
        
        if(mIdx < minTimeIdx or
           mIdx > maxTimeIdx or
           abs(maxAmp) < 1e-12):
          # No peaks
          continue

        goodMax.append(maxAmp*1e9)
        peakAmp.append(maxAmp*1e9)
        peakTime.append((time[mIdx]-time[tIdx[0]])*1e3)
        voltCurve.append(((time[tIdx]-time[tIdx[0]])*1e3,curAmp*1e9))
        
      # Pick which curves to plot
      sortIdx = np.argsort(peakAmp)
      if(len(sortIdx) < nPlotMax):
        keepIdx = sortIdx
      else:
        keepIdx = [sortIdx[int(np.round(x))] for x in \
                   np.linspace(0,len(sortIdx)-1,nPlotMax)]

      for x in keepIdx:
        plt.plot(voltCurve[x][0],voltCurve[x][1],'k-')
      
      plt.scatter(peakTime,peakAmp,marker=".",c="blue",s=100)


      nType = self.data["neurons"][plotID[0]]["type"]
      if((simType,nType) in self.expDataDict):
        expData = self.expDataDict[(simType,nType)]
        t = self.tWindow*1e3*(1+0.03*np.random.rand(expData.shape[0]))
        plt.scatter(t, -expData, marker=".", c="red",s=100)

      if(self.plotExpTrace and (simType,nType) in self.expTraceDict):
        data = self.expTraceDict[(simType,nType)]
        tExp = data[:,0]
        vExp = data[:,1:]
        tIdx = np.where(tExp < self.tWindow*1e3)[0]
        plt.plot(tExp[tIdx],vExp[tIdx,:],c="red")
        
      plt.title(self.neuronName(preType) + " to " + self.neuronName(plotType))
      plt.xlabel("Time (ms)")
      plt.ylabel("Current (nA)")

      # Remove part of the frame
      plt.gca().spines["right"].set_visible(False)
      plt.gca().spines["top"].set_visible(False)
      
      plt.tight_layout()
      plt.ion()
      plt.show()
      plt.savefig(figName,dpi=300)

      # Also plot histogram
      plt.figure()
      plt.hist(goodMax)
      plt.xlabel("Current (nA)")
      plt.title(self.neuronName(preType) + " to " + self.neuronName(plotType))

      # Remove part of the frame
      plt.gca().spines["right"].set_visible(False)
      plt.gca().spines["top"].set_visible(False)
      
      plt.tight_layout()
      plt.ion()
      plt.show()
      plt.savefig(figNameHist,dpi=300)
      
    import pdb
    pdb.set_trace()

    
  ############################################################################

  def setGABArev(self,vRevCl):

    print("Setting GABA reversal potential to " + str(vRevCl*1e3) + " mV")
    
    for s in self.snuddaSim.synapse_list:
      assert s.e == -65, "It should be GABA synapses only that we modify!"
      s.e = vRevCl * 1e3
    
  ############################################################################

    
if __name__ == "__main__":

  import argparse
  parser = argparse.ArgumentParser("Simulate Chuhma 2011 and Straub2016 experiments")
  parser.add_argument("action",choices=["setup","run","analyse"])
  parser.add_argument("simType",help="Experiment we want to perform",
                      choices=["Chuhma2011",
                               "Straub2016FS",
                               "Straub2016LTS",
                               "Szydlowski2013"])

  parser.add_argument("simName",
                      help="Simulation name, eg. networks/Chuhma2011-v1",
                      type=str)

  args = parser.parse_args()

  simName = args.simName
  simType = args.simType

  sm = SnuddaModelCurrentInjections(simName,simType)
  
  if(args.action == "setup"):
    print("Setup " + simName)
    # simName = "networks/Chuhma2011-v1"

    sm.defineNetwork(simName)
    sm.createNetwork(simName)

  if(args.action == "run"):
    print("Running " + simName)
    sm.simulateNetwork(simName)

  if(args.action == "analyse"):
    print("Analyse " + simName)
    sm.analyseNetwork(simType=simType,simName=simName)

    
