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
# python3 snudda_model_current_injections.py setup networks/Chuhma2011-v14
#
# mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_model_current_injections.py run networks/Chuhma2011-v14/
#
# python3 snudda_model_current_injections.py analyse networks/Chuhma2011-v14/
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
               curInj = 10e-9):

    self.simName = simName

    self.tInj = 0.3
    self.injDuration = 1e-3
    self.curInj = 10e-9
    self.tWindow = 0.03
    self.simEnd = self.tInj + self.tWindow*2
    self.holdV = -60e-3
    self.GABArev = 2e-3 # 144mM inside, 133.5mM outside, 32-33C --NERNST--> 2mV


  ############################################################################
  
  def defineNetwork(self,simName):

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



    else:
      cnc.defineStriatum(nMSD1=590,nMSD2=590,nFS=5,nLTS=0,nChIN=20,
                         volumeType="slice",sideLen=1000e-6)    

    dirName = os.path.dirname(configName)
  
    if not os.path.exists(dirName):
      os.makedirs(dirName)

    cnc.writeJSON(configName)

  ############################################################################

  def simulateNetwork(self,simName):

    logFile = simName + "/log/simlog.txt"
    
    self.networkFile = simName + "/network-pruned-synapses.hdf5"
    
    self.snuddaSim = SnuddaSimulate(networkFile=self.networkFile,
                                    inputFile=None,
                                    logFile=logFile,
                                    disableGapJunctions=True)
    
    # Get neuronID of neurons that will get artificial stimulation
    self.stimID = [x["neuronID"] \
                   for x in self.snuddaSim.network_info["neurons"] \
                   if "SPN" in x["type"]]

    # Pick neurons that we will measure from
    self.measureFSN = [x["neuronID"] \
                       for x in self.snuddaSim.network_info["neurons"] \
                       if x["type"] == "FSN"]
    
    self.measureChIN = [x["neuronID"] \
                       for x in self.snuddaSim.network_info["neurons"] \
                       if x["type"] == "ChIN"]

    # For SPN we want to make sure we find ones in the centre
    self.measuredSPN = [x["neuronID"] \
                       for x in self.snuddaSim.network_info["neurons"] \
                       if x["type"] == "dSPN"][:20]
    self.measureiSPN = [x["neuronID"] \
                       for x in self.snuddaSim.network_info["neurons"] \
                       if x["type"] == "iSPN"][:20]

    # Remove the overlap, ie dont stimulate the neurons we measure from    
    self.stimID = np.setdiff1d(self.stimID,
                               np.union1d(self.measuredSPN,self.measureiSPN))
    
    # Set up stimulation protocol
    for nID in self.stimID:
      self.snuddaSim.addCurrentInjection(neuronID=nID,
                                         startTime=self.tInj,
                                         endTime=self.tInj+self.injDuration,
                                         amplitude=self.curInj)

    # Add recordings
    self.snuddaSim.addVoltageClamp(cellID = self.measureFSN,
                                   voltage = self.holdV,
                                   duration=self.simEnd,
                                   saveIflag=True)
    self.snuddaSim.addVoltageClamp(cellID = self.measureChIN,
                                   voltage = self.holdV,
                                   duration=self.simEnd,
                                   res=1e-6,
                                   saveIflag=True)
    self.snuddaSim.addVoltageClamp(cellID = self.measuredSPN,
                                   voltage = self.holdV,
                                   duration=self.simEnd,
                                   saveIflag=True)
    self.snuddaSim.addVoltageClamp(cellID = self.measureiSPN,
                                   voltage = self.holdV,
                                   duration=self.simEnd,
                                   saveIflag=True)
    
    # Also add voltage recording for debugging reasons
    saveVoltage = True # False #True
    if(saveVoltage):
      self.snuddaSim.addRecording(cellID=self.stimID)

    self.setGABArev(self.GABArev)
    
    self.snuddaSim.run(self.simEnd*1e3)

    self.currentFile = simName + "/Chuhma2011-network-stimulation-current.txt"
    self.snuddaSim.writeCurrent(self.currentFile)

    if(saveVoltage):
      voltageFile = simName  + "/Chuhma2011-network-stimulation-voltage.txt"
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

  def analyseNetwork(self,simName):

    figDir = simName + "/figures/"
    if(not os.path.exists(figDir)):
      os.makedirs(figDir)

    
    print("Analysing data in " + simName)
    voltFile = simName + "/Chuhma2011-network-stimulation-current.txt"

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

    dSPNID = [x for x in current if self.data["neurons"][x]["type"] == "dSPN"]
    iSPNID = [x for x in current if self.data["neurons"][x]["type"] == "iSPN"]
    FSNID = [x for x in current if self.data["neurons"][x]["type"] == "FSN"]
    ChINID = [x for x in current if self.data["neurons"][x]["type"] == "ChIN"]

    dSPNmaxIdx = [np.argmax(np.abs(current[x]-current[x][0])) \
                  for x in dSPNID]
    iSPNmaxIdx = [np.argmax(np.abs(current[x]-current[x][0])) \
                  for x in iSPNID]
    FSNmaxIdx = [np.argmax(np.abs(current[x]-current[x][0])) \
                  for x in FSNID]
    ChINmaxIdx = [np.argmax(np.abs(current[x]-current[x][0])) \
                  for x in ChINID]

    minTimeIdx = np.where(time > self.tInj)[0][0]
    maxTimeIdx = np.where(time > self.tInj + self.tWindow)[0][0]

    matplotlib.rcParams.update({'font.size': 22})    
    
    for plotID,maxIdx in [(dSPNID,dSPNmaxIdx),
                          (iSPNID,iSPNmaxIdx),
                          (FSNID,FSNmaxIdx),
                          (ChINID,ChINmaxIdx)]:

      if(len(plotID) == 0):
        continue
      
      plotType = self.data["neurons"][plotID[0]]["type"]
      figName = figDir + plotType + "-current-traces.pdf"
      figNameHist = figDir + plotType + "-current-histogram.pdf"      

      goodMax = []
      
      plt.figure()
      for pID,mIdx in zip(plotID,maxIdx):
        if(mIdx < minTimeIdx or
           mIdx > maxTimeIdx):
          # No peaks
          continue

        
        tIdx = np.where(np.logical_and(time > self.tInj,
                                       time < self.tInj+self.tWindow))[0]
        plt.plot(time[tIdx]*1e3,
                 (current[pID][tIdx]-current[pID][tIdx[0]-1])*1e9,
                 c="black")

        plt.plot(time[mIdx]*1e3,
                 (current[pID][mIdx]-current[pID][tIdx[0]-1])*1e9,
                 marker=".",c="red")

        goodMax.append((current[pID][mIdx]-current[pID][tIdx[0]-1])*1e9)

        
      plt.title(plotType)
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
      plt.title(plotType)      
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
  parser = argparse.ArgumentParser("Simulate Chuhma 2011 experiment")
  parser.add_argument("action",choices=["setup","run","analyse"])
  parser.add_argument("simName",
                      help="Simulation name, eg. networks/Chuhma2011-v1",
                      type=str)

  args = parser.parse_args()

  simName = args.simName

  sm = SnuddaModelCurrentInjections(simName)
  
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

    
