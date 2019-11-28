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
# * Compile mod files (we now have failure rates for GABA)
#
# nrnivmodl cellspecs/mechanisms/
#
# * Cut the slice, so z > 0.00504 is kept
#
#  python3 snudda_cut.py networks/SynTest-v2/network-pruned-synapses.hdf5 "abs(z)<100e-6"
#
# * Look at networks/SynTest-v2/network-cut-slice.hdf5.pdf to verify cut plane
#
# * Run dSPN -> iSPN calibration (you get dSPN -> dSPN data for free then)
#
#  mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run networks/SynTest-v2/network-cut-slice.hdf5 dSPN iSPN
#
# *  Analyse
#
#  python3 snudda_calibrate_synapses.py analyse networks/SynTest-v2/network-cut-slice.hdf5 dSPN iSPN
# python3 snudda_calibrate_synapses.py analyse networks/SynTest-v2/network-cut-slice.hdf5 dSPN dSPN
#
# * Look at plot with traces overlayed and histogram of voltage amplitudes
# (When you do preType to postType, you also get preType to preType for free
# since both voltages are recorded

import os
import glob
import numpy as np
from snudda_simulate import SnuddaSimulate
from snudda_load import SnuddaLoad
import matplotlib
import matplotlib.pyplot as plt

# We want to match Taverna 2008 data:

# The slices were 300 μm thick.  MSNs sampled were 25–100 μm from the
# surface of the slice to facilitate paired recordings.  Pairs always
# had cell bodies that were similar in depth (near the same focal
# plane) in an attempt to minimize the probability that the local axon
# collateral of one or the other cell was cut. After choosing a given
# pair of cells (lateral distance between somata, ≤50 μm)
#
# Assume loss of 100 micrometer in depth, at KI they start with 250 micrometers
# and get 150 micrometers after.

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

    self.injSpacing = 0.2 # 0.5
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

    # !!! We could maybe update code so that for postType == "ALL" we
    # record voltage from all neurons

    if(self.postType == "ALL"):
      self.snudda.addRecording()
    else:
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

      allFile = self.voltFileAltMask.replace("*","ALL")
      
      if(os.path.exists(allFile)):
        print("Using " + allFile + " instead")
        voltFile = allFile
      elif(self.preType == self.postType):
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
  
  def analyse(self,maxDist=50e-6):

    # Read the data
    self.snuddaLoad = SnuddaLoad(self.networkFile)
    self.data = self.snuddaLoad.data

    time,voltage = self.readVoltage(self.voltFile) # sets self.voltage
    checkWidth = 0.05

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
    tooFarAway = 0
    
    for (preID,t) in self.injInfo:
      # Post synaptic neurons to preID
      synapses,coords = self.snuddaLoad.findSynapses(preID=preID)

      postIDset = set(synapses[:,1]).intersection(self.possiblePostID)
      prePos = self.snuddaLoad.data["neuronPositions"][preID,:]
      
      for postID in postIDset:

        if(maxDist is not None):
          postPos = self.snuddaLoad.data["neuronPositions"][postID,:]
          if(np.linalg.norm(prePos-postPos) > maxDist):
            tooFarAway += 1
            continue
        
        # There is a bit of synaptic delay, so we can take voltage
        # at first timestep as baseline
        tIdx = np.where(np.logical_and(t <= time, time <= t + checkWidth))[0]
        synapseData.append((time[tIdx],voltage[postID][tIdx]))

    print("Number of pairs excluded, distance > " \
          + str(maxDist*1e6) + "mum : " + str(tooFarAway))
        
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

    # Extract the amplitude of all voltage pulses
    amp = np.zeros((len(synapseData),))
    idxMax = np.zeros((len(synapseData),),dtype=int)    
    tMax = np.zeros((len(synapseData),))
    
    for i,(t,v) in enumerate(synapseData):
      # Save the largest deflection -- with sign
      idxMax[i] = np.argmax(np.abs(v-v[0]))
      tMax[i] = t[idxMax[i]] - t[0]
      amp[i] = v[idxMax[i]]-v[0]

    assert len(amp) > 0, "No responses... too short distance!"
      
    print("Min amp: " + str(np.min(amp)))
    print("Max amp: " + str(np.max(amp)))
    print("Mean amp: " + str(np.mean(amp)) + " +/- " + str(np.std(amp)))
    print("Amps: " + str(amp))
      
      
    # Now we have all synapse deflections in synapseData
    matplotlib.rcParams.update({'font.size': 22})    
    plt.figure()
    for t,v in synapseData:
      plt.plot((t-t[0])*1e3,(v-v[0])*1e3,color="black")

    plt.scatter(tMax*1e3,amp*1e3,color="red",marker=".")
    plt.xlabel("Time (ms)")
    plt.ylabel("Voltage (mV)")
    plt.title(str(len(synapseData)) + " traces")
    plt.tight_layout()
    plt.ion()
    plt.show()
    plt.savefig(traceFig)

      

    plt.figure()
    plt.hist(amp*1e3,bins=30)
    plt.xlabel("Voltage deflection (mV)")
    plt.tight_layout()
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
  parser.add_argument("postType",help="Post synaptic neuron type (for run task, postType can be 'ALL' to record from all neurons)")
  args = parser.parse_args()
  
  scs = SnuddaCalibrateSynapses(networkFile=args.networkFile,
                                preType=args.preType,
                                postType=args.postType)

  if(args.task == "run"):
    scs.runSim()

  elif(args.task == "analyse"):
    scs.analyse()
  
  
  
