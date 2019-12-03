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
# python3 snudda_model_current_injections.py setup networks/Chuhma2011-v14 Chuhma2011
#
# mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_model_current_injections.py run networks/Chuhma2011-v14/ Chuhma2011
#
# python3 snudda_model_current_injections.py analyse networks/Chuhma2011-v14/ Chuhma2011
#
#

import os

from snudda_simulate import SnuddaSimulate
from snudda_load import SnuddaLoad
from snudda_init import SnuddaInit
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
      


  ############################################################################
  
  def defineNetwork(self,simName,simType=None):

    if(simType is None):
      simType = self.simType
            
    configName= simName + "/network-config.json"
    cnc = SnuddaInit(structDef={},configName=configName,nChannels=1)

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
      cnc.defineStriatum(nMSD1=153,nMSD2=153,nFS=10,nLTS=0,nChIN=10,
                         volumeType="slice",sideLen=500e-6)    



    if(simType == "Chuhma2011"):
      cnc.defineStriatum(nMSD1=570+self.nNrns,
                         nMSD2=570+self.nNrns,
                         nFS=5,
                         nLTS=0,
                         nChIN=self.nNrns,
                         volumeType="slice",sideLen=1000e-6)    

    elif(simType == "Straub2016FS"):
      # nFS must be correct density, but readout neurons can be any density
      cnc.defineStriatum(nMSD1=self.nNrns,
                         nMSD2=self.nNrns,
                         nFS=156, nLTS=0,
                         nChIN=self.nNrns,
                         volumeType="slice",sideLen=1000e-6)    

    elif(simType == "Straub2016LTS"):
      cnc.defineStriatum(nMSD1=self.nNrns,
                         nMSD2=self.nNrns,
                         nFS=0,nLTS=84,
                         nChIN=self.nNrns,
                         volumeType="slice",sideLen=1000e-6)    
    elif(simType == "Szydlowski2013"):
      cnc.defineStriatum(nMSD1=0,
                         nMSD2=0,
                         nFS=156,
                         nLTS=self.nNrns,
                         nChIN=0,
                         volumeType="slice",sideLen=1000e-6)          
    else:
      print("setup : Unkown simType: " + str(simType))
      exit(-1)
      
    dirName = os.path.dirname(configName)
  
    if not os.path.exists(dirName):
      os.makedirs(dirName)

    cnc.writeJSON(configName)

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
      self.networkFile = simName + "/network-pruned-synapses.hdf5"
      
      self.snuddaSim = SnuddaSimulate(networkFile=self.networkFile,
                                      inputFile=None,
                                      logFile=logFile,
                                      disableGapJunctions=True)
    

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
      
      self.snuddaSim = SnuddaSimulate(networkFile=self.networkFile,
                                      inputFile=None,
                                      logFile=logFile,
                                      disableGapJunctions=True)
    
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
      
      self.snuddaSim = SnuddaSimulate(networkFile=self.networkFile,
                                      inputFile=None,
                                      logFile=logFile,
                                      disableGapJunctions=True)
    
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
          
      self.snuddaSim = SnuddaSimulate(networkFile=self.networkFile,
                                      inputFile=None,
                                      logFile=logFile,
                                      disableGapJunctions=True)
    
    # Set up stimulation protocol
    for nID in stimID:
      self.snuddaSim.addCurrentInjection(neuronID=nID,
                                         startTime=self.tInj,
                                         endTime=self.tInj+self.injDuration,
                                         amplitude=self.curInj)

    # Add recordings
    self.snuddaSim.addVoltageClamp(cellID = measureID,
                                   voltage = self.holdV,
                                   duration=self.simEnd,
                                   saveIflag=True)
    
    # Also add voltage recording for debugging reasons
    saveVoltage = True # False #True
    if(saveVoltage):
      self.snuddaSim.addRecording(cellID=stimID)

    self.setGABArev(self.GABArev)
    
    self.snuddaSim.run(self.simEnd*1e3)

    self.currentFile = simName + "/" + simType \
      + "-network-stimulation-current.txt"
    
    self.snuddaSim.writeCurrent(self.currentFile)

    if(saveVoltage):
      voltageFile = simName  + "/" + simType \
        + "-network-stimulation-voltage.txt"
      
      self.snuddaSim.writeVoltage(voltageFile)
      
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
    
    
    sn.placeNeurons(args)
    sn.touchDetection(args)
    sn.pruneSynapses(args)

  ############################################################################

  def analyseNetwork(self,simName,simType=None):

    figDir = simName + "/figures/"
    if(not os.path.exists(figDir)):
      os.makedirs(figDir)

    if(simType is None):
      simType = self.simType

    if(simType == "Straub2016LTS"):
      preType = "LTS"
    elif(simType == "Straub2016FS"):
      preType = "FSN"
    elif(simType == "Chuhma2011"):
      preType = "SPN"
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

        
        plt.plot(time[tIdx]*1e3,
                 curAmp*1e9,
                 c="black")

        plt.plot(time[mIdx]*1e3,
                 maxAmp*1e9,
                 marker=".",c="blue")

        goodMax.append(maxAmp*1e9)

        
      plt.title(preType + " to " + plotType)
      plt.xlabel("Time (ms)")
      plt.ylabel("Current (nA)")
      plt.tight_layout()
      plt.ion()
      plt.show()
      plt.savefig(figName,dpi=300)

      # Also plot histogram
      plt.figure()
      plt.hist(goodMax)
      plt.xlabel("Current (nA)")
      plt.title(preType + " to " + plotType)      
      plt.tight_layout()
      plt.ion()
      plt.show()
      plt.savefig(figNameHist,dpi=300)
      
    import pdb
    pdb.set_trace()

    
  ############################################################################

  def setGABArev(self,vRevCl):

    print("Setting GABA reversal potential to " + str(vRevCl*1e3) + " mV")
    
    for s in self.snuddaSim.synapseList:
      assert s.e == -65, "It should be GABA synapses only that we modify!"
      s.e = vRevCl * 1e3
    
  ############################################################################

    
if __name__ == "__main__":

  import argparse
  parser = argparse.ArgumentParser("Simulate Chuhma 2011 and Straub2016 experiments")
  parser.add_argument("action",choices=["setup","run","analyse"])
  parser.add_argument("simName",
                      help="Simulation name, eg. networks/Chuhma2011-v1",
                      type=str)
  parser.add_argument("simType",help="Experiment we want to perform",
                      choices=["Chuhma2011",
                               "Straub2016FS",
                               "Straub2016LTS",
                               "Szydlowski2013"])

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
    sm.analyseNetwork(simName)

    
