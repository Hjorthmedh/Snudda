# This code assumes you have created a small network of neurons, it will
# setup current injections
#
# How to:
#
# * Edit snudda_init_custom.py to have the neurons you want.
#
# * Generate network:
#
#  python3 snudda_init_custom.py
#  python3 snudda.py place networks/SynTest-v2
#  python3 snudda.py detect networks/SynTest-v2
#  python3 snudda.py prune networks/SynTest-v2

# * Figure out where to put the slcie cut (with plotOnly equation is ignored)
# 
#  python3 snudda_cut.py networks/SynTest-v2/network-pruned-synapses.hdf5 "z>0" --plotOnly
#
# * Look at networks/SynTest-v2/network-cut-slice.hdf5.pdf to decide cut plane
#
# * Cut the slice, so z > 0.00504 is kept
#
#  python3 snudda_cut.py networks/SynTest-v2/network-pruned-synapses.hdf5 "z>0.00504"
#
# * Look at networks/SynTest-v2/network-cut-slice.hdf5.pdf to verify cut plane
#
# * Run dSPN -> iSPN calibration (you get dSPN -> dSPN data for free then)
#
#  python3 snudda_calibrate_synapses.py run networks/SynTest-v2/network-cut-slice.hdf5 dSPN iSPN
#
# *  Analyse
#
#  python3 snudda_calibrate_synapses.py analyse networks/SynTest-v2/network-cut-slice.hdf5 dSPN iSPN
#
# * Look at plot with traces overlayed and histogram of voltage amplitudes

import os
import glob
import numpy as np
from snudda_simulate import SnuddaSimulate
from snudda_load import SnuddaLoad
import matplotlib.pyplot as plt

class SnuddaCalibrateSynapses():

  def __init__(self,networkFile,preType,postType,curInj=10e-9,logFile=None):

    if(os.path.isdir(networkFile)):
      self.networkFile = networkFile + "/network-pruned-synapses.hdf5"
    else:
      self.networkFile = networkFile
      
    self.preType = preType
    self.postType = postType
    self.curInj = curInj
    self.logFile = logFile
    
    print("Checking depolarisation/hyperpolarisation of " + preType \
          + " to " + postType + "synapses")

    self.injSpacing = 0.5
    self.injDuration = 1e-3

    # Voltage file
    self.voltFile = os.path.dirname(networkFile) \
      + "/synapse-calibration-volt-" \
      + self.preType + "-" + self.postType + ".txt"
    self.voltFileAltMask = os.path.dirname(networkFile) \
      + "/synapse-calibration-volt-" \
      + self.preType + "-*.txt"
    
  ############################################################################
    
  def runSim(self):
    
    self.snuddaSim = SnuddaSimulate(networkFile=self.networkFile,
                                    inputFile=None,
                                    logFile=self.logFile,
                                    disableGapJunctions=True)
    
    # A current pulse to all pre synaptic neurons, one at a time
    self.preID = [x["neuronID"] \
                  for x in self.snuddaSim.network_info["neurons"] \
                  if x["type"] == self.preType]

    # injInfo contains (preID,injStartTime)
    self.injInfo = list(zip(self.preID, \
                            self.injSpacing\
                            +self.injSpacing*np.arange(0,len(self.preID))))
    
    simEnd = self.injInfo[-1][1] + self.injSpacing
    
    # Add current injections defined in init
    for (nid,t) in self.injInfo:
      print("Current injection to " + str(nid) + " at " + str(t) + " s")
      self.snuddaSim.addCurrentInjection(neuronID=nid,
                                         startTime=t,
                                         endTime=t+self.injDuration,
                                         amplitude=self.curInj)
      
    # Record from all the potential post synaptic neurons
    self.snuddaSim.addRecordingOfType(self.postType)

    # Also save the presynaptic traces for debugging, to make sure they spike
    self.snuddaSim.addRecordingOfType(self.preType)

    
    # Run simulation
    self.snuddaSim.run(simEnd*1e3)
    
    # Write results to disk
    self.snuddaSim.writeVoltage(self.voltFile)


  ############################################################################

  def readVoltage(self,voltFile):

    if(not os.path.exists(voltFile)):
      print("Missing " + voltFile)
      
      if(self.preType == self.postType):
        fileList = glob.glob(self.voltFileAltMask)
        if(len(fileList) > 0):
          voltFile = fileList[0]
          print("Using " + voltFile + " instead, since pre and post are same")
      else:
        print("Aborting")
        exit(-1)
    
    data = np.genfromtxt(voltFile, delimiter=',')
    assert(data[0,0] == -1) # First column should be time
    time = data[0,1:] / 1e3
    
    voltage = dict()
    
    for rows in data[1:,:]:
      cID = int(rows[0])
      voltage[cID] = rows[1:] * 1e-3

    return (time,voltage)
      
  ############################################################################
  
  # This extracts all the voltage deflections, to see how strong they are
  
  def analyse(self):

    # Read the data
    self.snuddaLoad = SnuddaLoad(self.networkFile)
    self.data = self.snuddaLoad.data

    time,voltage = self.readVoltage(self.voltFile) # sets self.voltage
    checkWidth = 0.4

    # Generate current info structure
    # A current pulse to all pre synaptic neurons, one at a time
    self.preID = [x["neuronID"] \
                  for x in self.data["neurons"] \
                  if x["type"] == self.preType]

    self.possiblePostID = [x["neuronID"] \
                           for x in self.data["neurons"] \
                           if x["type"] == self.postType]
    
    # injInfo contains (preID,injStartTime)
    self.injInfo = zip(self.preID, \
                       self.injSpacing\
                       +self.injSpacing*np.arange(0,len(self.preID)))

    
    # For each pre synaptic neuron, find the voltage deflection in each
    # of its post synaptic neurons

    synapseData = []
    
    for (preID,t) in self.injInfo:
      # Post synaptic neurons to preID
      synapses,coords = self.snuddaLoad.findSynapses(preID=preID)

      postIDset = set(synapses[:,1]).intersection(self.possiblePostID)
      
      for postID in postIDset:

        # There is a bit of synaptic delay, so we can take voltage
        # at first timestep as baseline
        tIdx = np.where(np.logical_and(t <= time, time <= t + checkWidth))[0]
        synapseData.append((time[tIdx],voltage[postID][tIdx]))

    # Fig names:
    traceFig = os.path.dirname(self.networkFile) \
      + "/figures/synapse-calibration-volt-traces-" \
      + self.preType + "-" + self.postType + ".pdf"

    histFig = os.path.dirname(self.networkFile) \
      + "/figures/synapse-calibration-volt-histogram-" \
      + self.preType + "-" + self.postType + ".pdf"

    figDir = os.path.dirname(self.networkFile) + "/figures"
    
    if(not os.path.exists(figDir)):
      os.makedirs(figDir)
    
    # Now we have all synapse deflections in synapseData
    plt.figure()
    for t,v in synapseData:
      plt.plot((t-t[0])*1e3,(v-v[0])*1e3,color="black")
    plt.xlabel("Time (ms)")
    plt.ylabel("Voltage (mV)")
    plt.title(str(len(synapseData)) + " traces")
    plt.ion()
    plt.show()
    plt.savefig(traceFig)

      
    # Extract the amplitude of all voltage pulses
    amp = np.zeros((len(synapseData),))

    for i,(t,v) in enumerate(synapseData):
      # Save the largest deflection -- with sign
      amp[i] = v[np.argmax(np.abs(v-v[0]))]-v[0]

    plt.figure()
    plt.hist(amp*1e3)
    plt.xlabel("Voltage deflection (mV)")
    plt.show()
    plt.savefig(histFig)
    
    import pdb
    pdb.set_trace()

if __name__ == "__main__":

  from argparse import ArgumentParser

  parser = ArgumentParser(description="Calibrate synapse conductances")
  parser.add_argument("task", choices=["run","analyse"])
  parser.add_argument("networkFile",help="Network file (hdf5)")
  parser.add_argument("preType",help="Pre synaptic neuron type")
  parser.add_argument("postType",help="Post synaptic neuron type")
  args = parser.parse_args()
  
  scs = SnuddaCalibrateSynapses(networkFile=args.networkFile,
                                preType=args.preType,
                                postType=args.postType)

  if(args.task == "run"):
    scs.runSim()

  elif(args.task == "analyse"):
    scs.analyse()
  
  
  
