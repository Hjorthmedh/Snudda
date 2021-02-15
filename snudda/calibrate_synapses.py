# This code assumes you have created a small network of neurons, it will
# setup current injections
#
#
# OBS! We set a holding current to keep neuron at around -80mV
# We also change GABA reversal potential to -40mV, since internal Cl is 30mM
# and external Cl is 135 mM, temperature is 306K
#
#
# How to:
#
# * Edit snudda_init_custom.py to have the neurons you want.
#
# * Generate network 
#
#  python3 snudda_calibrate_synapses.py setup Planert2010 networks/Planert2010-v1
#  python3 snudda.py place networks/Planert2010-v1
#  python3 snudda.py detect networks/Planert2010-v1
#  python3 snudda.py prune networks/Planert2010-v1

# * Figure out where to put the slcie cut (with plotOnly equation is ignored)
# 
#  python3 snudda_cut.py networks/Planert2010-v1/network-pruned-synapses.hdf5 "z>0" --plotOnly
#
# * Look at networks/Planert2010-v1/network-cut-slice.hdf5.pdf to decide cut plane
#
# * Compile mod files (we now have failure rates for GABA)
#
# nrnivmodl cellspecs/mechanisms/
#
# * Cut the slice, so z > 0.00504 is kept
#
#  python3 snudda_cut.py networks/Planert2010-v1/network-pruned-synapses.hdf5 "abs(z)<100e-6"
#
# * Look at networks/Planert2010-v1/network-cut-slice.hdf5.pdf to verify cut plane
#
# * Run dSPN -> iSPN calibration (you get dSPN -> dSPN data for free then)
#
#  mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run Planert2010 networks/Planert2010-v1/network-cut-slice.hdf5 --pre dSPN --post iSPN
#
# *  Analyse
#
#  python3 snudda_calibrate_synapses.py analyse networks/Planert2010-v1/network-cut-slice.hdf5 dSPN iSPN
# python3 snudda_calibrate_synapses.py analyse Planert2010 networks/Planert2010-v1/network-cut-slice.hdf5 --pre dSPN --post dSPN
#
# * Look at plot with traces overlayed and histogram of voltage amplitudes
# (When you do preType to postType, you also get preType to preType for free
# since both voltages are recorded

import os
import glob
import numpy as np
from snudda.simulate import SnuddaSimulate
from snudda.load import SnuddaLoad
import matplotlib
import matplotlib.pyplot as plt
import neuron

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

class SnuddaCalibrateSynapses(object):

  def __init__(self,networkFile,
               preType,postType,
               curInj = 10e-9,
               holdV = -80e-3,
               maxDist = 50e-6,
               logFile=None):

    if(os.path.isdir(networkFile)):
      self.networkFile = networkFile + "/network-pruned-synapses.hdf5"
    else:
      self.networkFile = networkFile
      
    self.preType = preType
    self.postType = postType
    self.curInj = curInj
    self.holdV = holdV
    self.logFile = logFile
    self.maxDist = maxDist
    
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

    self.neuronNameRemap = {"FSN" : "FS"}

  ############################################################################
    
  def neuronName(self,neuronType):

    if(neuronType in self.neuronNameRemap):
      return self.neuronNameRemap[neuronType]
    else:
      return neuronType    
    
  ############################################################################

  def setup(self,simName,expType,nMSD1=120,nMSD2=120,nFS=20,nLTS=0,nChIN=0):

    from .init import SnuddaInit

    configName= simName + "/network-config.json"
    cnc = SnuddaInit(struct_def={}, config_file=configName, nChannels=1)
    cnc.define_striatum(num_dSPN=nMSD1, num_iSPN=nMSD2, num_FS=nFS, num_LTS=nLTS, num_ChIN=nChIN,
                        volume_type="slice", side_len=200e-6, slice_depth=150e-6)

    dirName = os.path.dirname(configName)
  
    if not os.path.exists(dirName):
      os.makedirs(dirName)

    cnc.write_json(configName)

    
    print("\n\npython3 snudda.py place " + str(simName))
    print("python3 snudda.py detect " + str(simName))
    print("python3 snudda.py prune " + str(simName))
    print("python3 snudda_cut.py " + str(simName) \
          + '/network-pruned-synapses.hdf5 "abs(z)<100e-6"')

    print("\nThe last command will pop up a figure and enter debug mode, press ctrl+D in the terminal window after inspecting the plot to continue")

    print("\n!!! Remember to compile the mod files: nrnivmodl cellspecs/mechanisms")

    print("\nTo run for example dSPN -> iSPN (and dSPN->dSPN) calibration:")
    print("mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run " + str(expType) + " " + str(simName) + "/network-cut-slice.hdf5 dSPN iSPN")

    print("\npython3 snudda_calibrate_synapses.py analyse " + str(expType) + " " + str(simName) + "/network-cut-slice.hdf5 --pre dSPN --post iSPN\npython3 snudda_calibrate_synapses.py analyse " + str(simName) + "/network-cut-slice.hdf5 --pre iSPN --post dSPN")
    
  ############################################################################

  def setupHoldingVolt(self,holdV=None,simEnd=None):

    assert simEnd is not None, \
      "setupHoldingVolt: Please set simEnd, for holding current"
    
    if(holdV is None):
      holdV = self.holdV

    if(holdV is None):
      print("Not using holding voltage, skipping.")
      return

    # Setup vClamps to calculate what holding current will be needed
    somaVClamp = []

    somaList = [self.snuddaSim.neurons[x].icell.soma[0] \
                for x in self.snuddaSim.neurons]
    
    for s in somaList:
      vc = neuron.h.SEClamp(s(0.5))
      vc.rs = 1e-9
      vc.amp1 = holdV*1e3
      vc.dur1 = 100

      somaVClamp.append((s,vc))

    neuron.h.finitialize(holdV*1e3)
    neuron.h.tstop = 100
    neuron.h.run()

    self.holdingIClampList = []

    # Setup iClamps    
    for s,vc in somaVClamp:
      cur = float(vc.i)
      ic = neuron.h.IClamp(s(0.5))
      ic.amp = cur
      ic.dur = 2*simEnd*1e3
      self.holdingIClampList.append(ic)
      
    # Remove vClamps
    vClamps = None
    vc = None  
    
  ############################################################################

  def setGABArev(self,vRevCl):

    print("Setting GABA reversal potential to " + str(vRevCl*1e3) + " mV")
    
    for s in self.snuddaSim.synapse_list:
      assert s.e == -65, "It should be GABA synapses only that we modify!"
      s.e = vRevCl * 1e3
          
  ############################################################################
  
  def runSim(self,GABArev):
    
    self.snuddaSim = SnuddaSimulate(network_file=self.networkFile,
                                    input_file=None,
                                    log_file=self.logFile,
                                    disable_gap_junctions=True)

    
    # A current pulse to all pre synaptic neurons, one at a time
    self.preID = [x["neuronID"] \
                  for x in self.snuddaSim.network_info["neurons"] \
                  if x["type"] == self.preType]

    # injInfo contains (preID,injStartTime)
    self.injInfo = list(zip(self.preID, \
                            self.injSpacing\
                            +self.injSpacing*np.arange(0,len(self.preID))))
    
    simEnd = self.injInfo[-1][1] + self.injSpacing
    
    # Set the holding voltage
    self.setupHoldingVolt(holdV=self.holdV,simEnd=simEnd)

    self.setGABArev(GABArev)
    
    
    # Add current injections defined in init
    for (nid,t) in self.injInfo:
      print("Current injection to " + str(nid) + " at " + str(t) + " s")
      self.snuddaSim.add_current_injection(neuron_id=nid,
                                           start_time=t,
                                           end_time=t + self.injDuration,
                                           amplitude=self.curInj)

    # !!! We could maybe update code so that for postType == "ALL" we
    # record voltage from all neurons

    if(self.postType == "ALL"):
      self.snuddaSim.add_recording()
    else:
      # Record from all the potential post synaptic neurons
      self.snuddaSim.add_recording_of_type(self.postType)

      # Also save the presynaptic traces for debugging, to make sure they spike
      self.snuddaSim.add_recording_of_type(self.preType)

    
    # Run simulation
    self.snuddaSim.run(simEnd * 1e3, hold_v=self.holdV)
    
    # Write results to disk
    self.snuddaSim.write_voltage(self.voltFile)


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
  
  def analyse(self,expType,maxDist=None,nMaxShow=10):

    self.setupExpData()
    
    if(maxDist is None):
      maxDist = self.maxDist
    
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
      synapses,coords = self.snuddaLoad.find_synapses(pre_id=preID)

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

    if(maxDist is not None):
      print("Number of pairs excluded, distance > " \
            + str(maxDist*1e6) + "mum : " + str(tooFarAway))
        
    # Fig names:
    traceFig = os.path.dirname(self.networkFile) \
      + "/figures/" + expType +"synapse-calibration-volt-traces-" \
      + self.preType + "-" + self.postType + ".pdf"

    histFig = os.path.dirname(self.networkFile) \
      + "/figures/" + expType + "synapse-calibration-volt-histogram-" \
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

    sortIdx = np.argsort(amp)
    if(len(sortIdx) > nMaxShow):
      keepIdx = [sortIdx[int(np.round(x))] \
                 for x in np.linspace(0,len(sortIdx)-1,nMaxShow)]
    else:
      keepIdx = sortIdx
      
    plt.figure()
    for x in keepIdx:

      t,v = synapseData[x]
      
      plt.plot((t-t[0])*1e3,(v-v[0])*1e3,color="black")
      
    plt.scatter(tMax*1e3,amp*1e3,color="blue",marker=".",s=100)

    if((expType,self.preType,self.postType) in self.expData):
      expMean,expStd = self.expData[(expType,self.preType,self.postType)]

      tEnd = (t[-1]-t[0])*1e3

      axes = plt.gca()
      ay = axes.get_ylim()
      # Plot SD or 1.96 SD?
      plt.errorbar(tEnd,expMean,expStd,ecolor="red",
                   marker='o',color="red")

      modelMean = np.mean(amp)*1e3
      modelStd = np.std(amp)*1e3
      plt.errorbar(tEnd-2,modelMean,modelStd,ecolor="blue",
                   marker="o",color="blue")
      
      axes.set_ylim(ay)
      
      
    plt.xlabel("Time (ms)")
    plt.ylabel("Voltage (mV)")
    #plt.title(str(len(synapseData)) + " traces")
    plt.title(self.neuronName(self.preType) \
              + " to " + self.neuronName(self.postType))

    # Remove part of the frame
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)
    
    plt.tight_layout()
    plt.ion()
    plt.show()
    plt.savefig(traceFig,dpi=300)

      

    plt.figure()
    plt.hist(amp*1e3,bins=20)
    plt.title(self.neuronName(self.preType) \
              + " to " + self.neuronName(self.postType))    
    plt.xlabel("Voltage deflection (mV)")

    # Remove part of the frame
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["top"].set_visible(False)    
    
    plt.tight_layout()
    plt.show()
    plt.savefig(histFig,dpi=300)
    
    import pdb
    pdb.set_trace()

############################################################################

  def setupExpData(self):

    self.expData = dict()

    planertD1D1 = (0.24,0.15)
    planertD1D2 = (0.33,0.15)
    planertD2D1 = (0.27,0.09)
    planertD2D2 = (0.45,0.44)
    planertFSD1 = (4.8,4.9)
    planertFSD2 = (3.1,4.1)

    self.expData[("Planert2010","dSPN","dSPN")] = planertD1D1
    self.expData[("Planert2010","dSPN","iSPN")] = planertD1D2
    self.expData[("Planert2010","iSPN","dSPN")] = planertD2D1
    self.expData[("Planert2010","iSPN","iSPN")] = planertD2D2    
    self.expData[("Planert2010","FSN","dSPN")]  = planertFSD1
    self.expData[("Planert2010","FSN","iSPN")]  = planertFSD2    


    
############################################################################
    
if __name__ == "__main__":

  from argparse import ArgumentParser

  parser = ArgumentParser(description="Calibrate synapse conductances")
  parser.add_argument("task", choices=["setup","run","analyse"])
  parser.add_argument("expType",help="Experiment we replicate",
                      choices=["Planert2010","Szydlowski2013"])
  parser.add_argument("networkFile", \
                      help="Network file (hdf5) or network directory")
  parser.add_argument("--preType","--pre",
                      help="Pre synaptic neuron type",
                      default="dSPN")
  parser.add_argument("--postType","--post",
                      help="Post synaptic neuron type (for run task, postType can be 'ALL' to record from all neurons)",
                      default="ALL")
  parser.add_argument("--maxDist",help="Only check neuron pairs within (mum)",
                      type=float,default=None)
  args = parser.parse_args()
  
  if(args.maxDist is None):
    maxDist = 50e-6
  elif(args.maxDist == "None"):
    maxDist = None
  else:
    maxDist = float(args.maxDist)

  print("Using maxDist = " + str(maxDist))
    

  if(args.expType == "Planert2010"):
    nMSD1 = 120
    nMSD2 = 120
    nFS   = 20
    nLTS  = 0
    nChIN = 0

    holdV = -80e-3
    maxDist = 100e-6 if args.maxDist is None else args.maxDist
    GABArev = -40e-3
    
  elif(args.expType == "Szydlowski2013"):
    nMSD1 = 10
    nMSD2 = 10
    nFS   = 20
    nLTS  = 20
    nChIN = 0

    holdV = -76e-3
    maxDist = 200e-6 if args.maxDist is None else args.maxDist
    GABArev = -39e-3
    
  else:
    print("Unknown expType = " + str(args.expType))
    exit(-1)

  scs = SnuddaCalibrateSynapses(networkFile=args.networkFile,
                                preType=args.preType,
                                postType=args.postType,
                                maxDist=maxDist,
                                holdV=holdV)
    
  if(args.task == "setup"):
    scs.setup(args.networkFile,
              expType=args.expType,
              nMSD1=nMSD1,nMSD2=nMSD2,
              nFS=nFS,nLTS=nLTS,nChIN=nChIN)
    
  elif(args.task == "run"):
    scs.runSim(GABArev=GABArev)

  elif(args.task == "analyse"):
    scs.analyse(args.expType)
  
  
  
