# snudda_detect.py
#
# Johannes Hjorth, Royal Institute of Technology (KTH)
# Human Brain Project 2019
#
# Requires Python version 3+
#
# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union’s Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).
#


import numpy as np
import math
import os
import sys
import itertools

import time
import timeit

import h5py
import json
import pickle

import re

from .Neuron_morphology import NeuronMorphology
from .place import SnuddaPlace
from .load import SnuddaLoad

status = None
hyperVoxelData = None

class SnuddaDetect(object):

  def __init__(self,
               configFile = None,
               positionFile = None,
               voxelSize=3e-6, #2e-6,
               hyperVoxelSize=100, #250, #100,
               verbose=True,
               logFileName=None,
               logFile=None,
               saveFile=None,
               workHistoryFile=None,
               restartDetectionFlag=True, # False = continue old detection
               SlurmID=0,
               volumeID=None,
               role="master",
               dView=None,
               lbView=None,
               rc=None,
               axonStumpIDFlag=False,
               h5libver="latest",
               debugFlag=False):


    if(rc is not None):
      dView = rc.direct_view(targets='all')
      lbView = rc.load_balanced_view(targets='all')
    
    assert role in ["master","worker"], \
      "SnuddaDetect: Role must be master or worker"
    self.role = role

    self.verbose = verbose
    self.h5libver = h5libver
    self.debugFlag = debugFlag

    self.logFile = logFile
    self.configFile = configFile
    self.positionFile = positionFile
    self.saveFile = saveFile

    self.workHistoryFile = workHistoryFile # Name of work history file
    self.workHistory = None # File pointer for actual file
    
    if(logFileName is None and logFile is not None):
      self.logFileName = logFile.name
    else:
      self.logFileName = logFileName
      
    self.setupLog()

    self.writeLog("Using hdf5 driver version: " + str(self.h5libver))
    
    mem = self.memory()
    self.writeLog(str(mem))

    
    self.SlurmID = int(SlurmID) # Make sure integer
    self.workersInitialised = False

    self.voxelSize=voxelSize
    self.hyperVoxelSize=hyperVoxelSize # = N,  N x N x N voxels in a hyper voxel
    self.hyperVoxelOrigo=np.zeros((3,))
    self.voxelOverflowCounter = 0
    
    self.nBins = hyperVoxelSize*np.ones((3,),dtype=int)
    self.writeLog("Each hyper voxel has %d x %d x %d voxels" \
                  % tuple(self.nBins))
    
    # These are voxels that axons/dend occupy
    self.axonVoxels = None
    self.dendVoxels = None

    self.volumeID = volumeID
    if(volumeID is not None):
      self.writeLog("Touch detection only " + str(volumeID))
    else:
      self.writeLog("Touch detecting all volumes")
    
    # These are counters, how many different axons/dend in the voxel
    self.axonVoxelCtr = None
    self.dendVoxelCtr = None

    self.maxAxonVoxelCtr = None
    self.maxDendVoxelCtr = None
    
    # Columns in hyperVoxelSynapses:
    # 0: sourceCellID, 1: destCellID, 2: voxelX, 3: voxelY, 4: voxelZ,
    # 5: hyperVoxelID, 6: channelModelID,
    # 7: sourceAxonSomaDist (not SI scaled 1e6, micrometers),
    # 8: destDendSomaDist (not SI scalled 1e6, micrometers)
    # 9: destSegID, 10: destSegX (int 0 - 1000, SONATA wants float 0.0-1.0)
    # 11: conductance (int, not SI scaled 1e12, in pS)
    # 12: parameterID
    #
    # Note on parameterID:
    # If there are n parameter sets for the particular synapse type, then
    # the ID to use is parameterID % n, this way we can reuse connectivity
    # if we add more synapse parameter sets later.

    self.hyperVoxelSynapses = None

    # Columns in hyperVoxelGapJunctions
    # 0: sourceCellID, 1: destCellID, 2: sourceSegID, 3: destSegID,
    # 4: sourceSegX, 5: destSegX, 6: voxelX, 7: voxelY, 8: voxelZ,
    # 9: hyperVoxelID, 10: conductance (integer, in pS)
    self.hyperVoxelGapJunctions = None

    self.hyperVoxelSynapseCtr = 0
    self.hyperVoxelGapJunctionCtr = 0

    self.hyperVoxelCoords = dict([])

    # This is used by the heap sort, when merging hdf5 files for the different
    # hyper voxels
    self.hyperVoxelSynapseLookup = None
    self.hyperVoxelGapJunctionLookup = None

    # Parameters for the HDF5 writing, this affects write speed
    self.synapseChunkSize = 10000
    self.gapJunctionChunkSize = 10000
    self.h5compression = "lzf"
    
    # This is an upper limit how many axon/dend we allow in each voxel max
    # 10 overflowed
    self.maxAxon = 45
    self.maxDend = 20

    self.maxNeurons = 10000 
    self.maxSynapses = 2000000
    self.maxGapJunctions = 100000

    # We have to dynamically create this lookup
    #self.synapseTypeLookup = { 1 : "GABA",
    #                           2 : "AMPA_NMDA",
    #                           3 : "GapJunction",
    #                           4 : "ACh",
    #                           5 : "NO"}
    #
    #self.synapseTypeReverseLookup = \
    #    {v: k for k, v in self.synapseTypeLookup.items()}
    
    self.connectivityDistributions = dict([])
    #self.connectivityDistributionsGJ = dict([])
    self.nextChannelModelID = 10
    
    self.prototypeNeurons = dict([])

    self.axonCumDensityCache = dict([])

    self.deleteOldMerge()
    
    # Rather than load all neuron morphologies, we only load prototypes
    self.readPrototypes(configFile=configFile,
                        axonStumpIDFlag=axonStumpIDFlag)
    
    # Read positions
    self.readNeuronPositions(positionFile)
    
    # Then we need to setup the workers

    if(self.role == "master"):
      
      self.setupParallel(dView=dView)

      if(self.workHistoryFile is None):        
        workDir = os.path.dirname(self.saveFile)
        workDir = workDir.replace("/voxels","/")
        logDir = workDir + "/log/"
        #self.workHistoryFile = self.saveFile.replace(".hdf5","-worklog.hdf5")
        #self.workHistoryFile = self.workHistoryFile.replace("/voxels/","/")
        self.workHistoryFile = logDir + "network-detect-worklog.hdf5"
        
        #workHistoryFile = re.sub("/voxels-\d+/","/",workHistoryFile)
      
        if(self.workHistoryFile == self.saveFile):
          self.writeLog("Unable to set a good worklog name")
          self.workHistoryFile = "worklog.hdf5"

      if(restartDetectionFlag):
        if(os.path.isfile(self.workHistoryFile)):
          self.writeLog("Removing old work history file")
          os.remove(self.workHistoryFile)

        # Setup new work history
        self.setupWorkHistory(self.workHistoryFile)
      else:
        # Open old file with work history
        print("Reusing old work history file " + str(self.workHistoryFile))
        self.workHistory = h5py.File(self.workHistoryFile,"r+",
                                     libver=self.h5libver)

        
      # For each neuron we need to find which hyper voxel it belongs to
      # (can be more than one)
      self.distributeNeuronsParallel(dView=dView)
      
      if(dView is not None):
        self.parallelProcessHyperVoxels(rc=rc,dView=dView)
        
      else:
        # We are running it in serial
        
        (allHyperIDs,nCompleted,remaining,self.voxelOverflowCounter) = \
          self.setupProcessHyperVoxelStateHistory()
        
        for hyperID in remaining: #self.hyperVoxels:
          (hyperID,nSyn,nGJ,execTime,voxelOverflowCtr) =\
            self.processHyperVoxel(hyperID)

          if(voxelOverflowCtr > 0):
            self.writeLog("!!! HyperID " + str(hyperID) + " OVERFLOWED " \
                          + str(voxelOverflowCtr) + " TIMES (" +\
                          str(execTime) + "s)")
            self.voxelOverflowCounter += voxelOverflowCtr
          else:
            self.writeLog("HyperID " + str(hyperID) + " completed - " \
                          + str(nSyn) + " synapses and " \
                          + str(nGJ) + " gap junctions found (" \
                          + str(execTime) + "s)")

          self.updateProcessHyperVoxelState(hyperID=hyperID,nSyn=nSyn,nGJ=nGJ,
                                            execTime=execTime,
                                         voxelOverflowCounter=voxelOverflowCtr)
          
          
    # We need to gather data from all the HDF5 files

  ############################################################################
    
  ############################################################################

  def __del__(self):

    if(self.workHistory is not None):
      try:
        self.workHistory.close()
      except:
        print("Work history already closed")

    if(self.logFile is not None):
      try:
        self.logFile.close()
      except:
        print("Log file already closed")
    

  ############################################################################

  # The original code had problems with not being able to access the nc
  # object on the worker from the map call.
  
  def parallelProcessHyperVoxels(self,rc=None,dView=None):

    self.writeLog("Starting parallelProcessHyperVoxels")
    
    startTime = timeit.default_timer()

    # Loads state if previously existed, otherwise creates new fresh history
    (allHyperIDs,nCompleted,remaining,self.voxelOverflowCounter) = \
        self.setupProcessHyperVoxelStateHistory()

    nWorkers = len(rc.ids)
    workerStatus = [None for x in rc.ids]
    workerIdx = 0
    jobIdx = 0
    busyCtr = 0
    noChangeCtr = 0
    nSyn = 1 # If nSyn is zero delay in loop is shorter

    self.writeLog("parallelProcessHyperVoxels: Using " \
                  + str(nWorkers) + " worker")
    
    while(jobIdx < len(remaining) or busyCtr > 0):

      if(workerStatus[workerIdx] is not None):
        
        # We have an async result, check status of it
        if(workerStatus[workerIdx].ready()):
          # Result is ready, get it

          hyperVoxelData = rc[workerIdx]["result"]

          hyperID          = hyperVoxelData[0]
          nSyn             = hyperVoxelData[1]
          nGJ              = hyperVoxelData[2]
          execTime         = hyperVoxelData[3]
          voxelOverflowCtr = hyperVoxelData[4]
          
          self.updateProcessHyperVoxelState(hyperID=hyperID,
                                            nSyn=nSyn,
                                            nGJ=nGJ,
                                            execTime=execTime,
                                            voxelOverflowCounter=\
                                              voxelOverflowCtr)
          workerStatus[workerIdx] = None
          rc[workerIdx]["result"] = None # Clear to be safe
          busyCtr -= 1

          if(voxelOverflowCtr > 0):
            self.writeLog("!!! HyperID " + str(hyperID) + " OVERFLOWED " \
                          + str(voxelOverflowCtr) + " TIMES (" +\
                          str(execTime) + "s)")
            self.voxelOverflowCounter += voxelOverflowCtr
          else:
            self.writeLog("HyperID " + str(hyperID) + " completed - " \
                          + str(nSyn) + " synapses found (" \
                          + str(execTime) + "s)")

      # Check that there are neurons in the hyper voxel, otherwise skip it.
      if(workerStatus[workerIdx] is None and jobIdx < len(remaining)):
        try:
          self.writeLog(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) \
                        + " Starting hyper voxel " + str(remaining[jobIdx]) \
                        + " on worker " + str(workerIdx))
          cmdStr = "result = nc.processHyperVoxel(" + str(remaining[jobIdx])+")"

          workerStatus[workerIdx] = rc[workerIdx].execute(cmdStr,block=False)

          jobIdx += 1
          busyCtr += 1
          noChangeCtr = 0
          
        except:
          print("!!! Problem with worker : " + str(workerIdx) + " (" \
                + str(nWorkers) + " total workers)")
          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)
          import pdb
          pdb.set_trace()
          
      else:
        noChangeCtr += 1

      workerIdx = (workerIdx + 1) % nWorkers

      if(noChangeCtr >= nWorkers):
        # If all workers were working, sleep for a bit, then reset counter
        if(nSyn > 0):
          time.sleep(0.010)
        else:
          # If previous hyper voxel had no synapses, be faster
          time.sleep(1e-6)
          
        noChangeCtr = 0
    
    endTime = timeit.default_timer()
    
    self.writeLog("Voxel overflows: " + str(self.voxelOverflowCounter))
    self.writeLog("Total number of synapses: " \
                  + str(np.sum(self.workHistory["nSynapses"][:])))
    self.writeLog("parallelProcessHyperVoxels: " \
                  +  str(endTime-startTime) + "s")

    self.workHistory.close()
    
    
  ############################################################################

  def setupWorkHistory(self,workHistoryFile=None):

    if(self.role != "master"):
      return
    
    if(workHistoryFile is None):
      workHistoryFile = self.workHistoryFile
    else:
      # Update internal state
      self.workHistoryFile = workHistoryFile

    self.writeLog("Work history file: " + str(self.workHistoryFile))
      
    self.workHistory = h5py.File(workHistoryFile,
                                 libver=self.h5libver)

    # We need to encode the connectivityDistributions tuple as a string
    # before saving it in json
    # we also need to parse this string and recreate a tuple afterwards

    
    tmpConDist = dict([])
    tmpConDistGJ = dict([])

    for keys in self.connectivityDistributions:
      tmpConDist["$$".join(keys)] = self.connectivityDistributions[keys]

    #for keys in self.connectivityDistributionsGJ:
    #  tmpConDistGJ["$$".join(keys)] = self.connectivityDistributionsGJ[keys]
      
    
    saveMetaData = [(self.SlurmID, "SlurmID"),
                    (self.configFile,"configFile"),
                    (self.positionFile,"positionFile"),
                    (self.voxelSize,"voxelSize"),
                    (self.hyperVoxelSize,"hyperVoxelSize"),
                    (self.axonStumpIDFlag,"axonStumpIDFlag"),
                    (json.dumps(self.config),"config"),
                    (json.dumps(tmpConDist),
                     "connectivityDistributions")]
#    ,
#                    (json.dumps(tmpConDistGJ),
#                     "connectivityDistributionsGJ")]

    
    if("meta" not in self.workHistory):
      self.writeLog("Writing metadata to work history file")
      meta = self.workHistory.create_group("meta")

      for data,dataName in saveMetaData:
        meta.create_dataset(dataName,data=data)
      
    else:
      self.writeLog("Work history file found, checking meta data")
      # Make sure config file etc match

      for data,dataName in saveMetaData:
        assert data == self.workHistory["meta/" + dataName][()], \
          dataName + " mismatch " + str(data) + " vs " \
          + str(self.workHistory["meta/" + dataName][()])              



    print("Write neuron data to file")

    networkGroup = self.workHistory.create_group("network")
    
    # Finally the neuron information
    neuronGroup = networkGroup.create_group("neurons")
    
    # If the name list is longer than 20 chars, increase S20
    nameList = [n["name"].encode("ascii","ignore") for n in self.neurons]
    strType = 'S'+str(max(1,max([len(x) for x in nameList])))
    neuronGroup.create_dataset("name", (len(nameList),), strType, nameList,
                               compression=self.h5compression)

    neuronIDlist = [n["neuronID"] for n in self.neurons]
    neuronGroup.create_dataset("neuronID",(len(neuronIDlist),), \
                               'int',neuronIDlist)

    # Just make sure there is at least one neuron in volumeIDlist
    # that i inside volumeID
    
    volumeSet = set([n["volumeID"] for n in self.neurons])
    assert self.volumeID is None or self.volumeID in volumeSet, "VolumeID contains no neurons: " + str(self.volumeID)

    
    volumeIDlist = [n["volumeID"].encode("ascii","ignore") \
                    for n in self.neurons]
    strTypeVID = 'S' + str(max(1,max([len(x) for x in volumeIDlist])))

    neuronGroup.create_dataset("volumeID", \
                               (len(volumeIDlist),),strTypeVID,volumeIDlist,
                               compression=self.h5compression)
      
    hocList = [n["hoc"].encode("ascii","ignore") for n in self.neurons]
    neuronGroup.create_dataset("hoc",(len(hocList),),'S100', hocList,
                               compression=self.h5compression)

      
    virtualNeuronList = np.array([n["virtualNeuron"] for n in self.neurons],
                                 dtype=bool)
    virtualNeuron = neuronGroup.create_dataset("virtualNeuron",
                                               data=virtualNeuronList,
                                               compression=self.h5compression)

    swcList = [n["morphology"].encode("ascii","ignore") for n in self.neurons]
    maxSwcLen = max([len(x) for x in swcList])
    neuronGroup.create_dataset("morphology",(len(swcList),),
                               'S'+str(maxSwcLen),swcList,
                               compression=self.h5compression)
    
    neuronPosition = neuronGroup.create_dataset("position",\
                                                (len(self.neurons),3),\
                                                "float",
                                                compression=self.h5compression)
    
    neuronRotation = neuronGroup.create_dataset("rotation",\
                                                (len(self.neurons),9),\
                                                "float",
                                                compression=self.h5compression)

    neuronDendRadius = neuronGroup.create_dataset("maxDendRadius", \
                                                  (len(self.neurons),),\
                                                  "float",
                                                 compression=self.h5compression)

    neuronAxonRadius = neuronGroup.create_dataset("maxAxonRadius", \
                                                  (len(self.neurons),),\
                                                  "float",
                                                 compression=self.h5compression)

    neuronParamID = neuronGroup.create_dataset("parameterID",
                                               (len(self.neurons),),\
                                               "int",
                                               compression=self.h5compression)

    neuronModulationID = neuronGroup.create_dataset("modulationID",
                                               (len(self.neurons),),\
                                               "int",
                                               compression=self.h5compression)

    
    for (i,n) in enumerate(self.neurons):
      neuronPosition[i]     = n["position"]
      neuronRotation[i]     = n["rotation"].reshape(1,9)
      neuronDendRadius[i]   = n["maxDendRadius"]
      neuronAxonRadius[i]   = n["maxAxonRadius"]
      neuronParamID[i]      = n["parameterID"]
      neuronModulationID[i] = n["modulationID"]
      

    # Store input information
    neuronGroup.create_dataset("populationUnitID", data=self.populationUnit,
                               compression=self.h5compression,dtype=int)
    
    neuronGroup.create_dataset("nPopulationUnits", data=self.nPopulationUnits)
    neuronGroup.create_dataset("populationUnitPlacementMethod", data=self.populationUnitPlacementMethod)

    # Variable for axon density "r", "xyz" or "" (No axon density)
    axonDensityType = [n["axonDensityType"].encode("ascii","ignore") \
                       if n["axonDensityType"] is not None \
                       else b"" \
                       for n in self.neurons]
    
    adStrType2 = "S"+str(max(1,max([len(x) if x is not None else 1 \
                                    for x in axonDensityType])))
    neuronGroup.create_dataset("axonDensityType", (len(axonDensityType),),
                               adStrType2,data=axonDensityType,
                               compression=self.h5compression)
    
    axonDensity = [n["axonDensity"].encode("ascii","ignore") \
                   if n["axonDensity"] is not None \
                   else b"" \
                   for n in self.neurons]
    adStrType = "S"+str(max(1,max([len(x) if x is not None else 1 \
                                   for x in axonDensity])))
    
    neuronGroup.create_dataset("axonDensity", (len(axonDensity),),
                               adStrType,data=axonDensity,
                               compression=self.h5compression)

    axonDensityRadius = [n["axonDensityRadius"] \
                         if n["axonDensity"] is not None \
                         and n["axonDensityType"] == "r" \
                         else np.nan for n in self.neurons]

    neuronGroup.create_dataset("axonDensityRadius",data=axonDensityRadius)

    # This is for the density function where it uses x,y,z
    axonDensityBoundsXYZ = np.nan*np.zeros((len(self.neurons),6))

    for ni, n in enumerate(self.neurons):
      
      if(n["axonDensity"] is None):
        # No axon density specified, skip
        continue
      
      if(n["axonDensityType"] == "xyz"):
        axonDensityBoundsXYZ[ni,:] = n["axonDensityBoundsXYZ"]
        
    neuronGroup.create_dataset("axonDensityBoundsXYZ",data=axonDensityBoundsXYZ)
    

    
  ############################################################################

  # Reading work
  
  def setupProcessHyperVoxelStateHistory(self):

    if("completed" in self.workHistory):
      self.writeLog("setupProcessHyperVoxelStateHistory: "\
                    + "Resuming from old state")
      # We already have a run in progress, load the state
      allHyperIDs = set(self.workHistory["allHyperIDs"])
      nCompleted = int(self.workHistory["nCompleted"][0])
      completed = set(self.workHistory["completed"][:nCompleted])
      remaining = self.sortRemainingBySize(allHyperIDs - completed)
      nSynapses = int(self.workHistory["nSynapses"][0])
      nGapJunctions = int(self.workHistory["nGapJunctions"][0])
      voxelOverflowCounter = self.workHistory["voxelOverflowCounter"][0]

    else:
      self.writeLog("setupProcessHyperVoxelStateHistory: " \
                    + "Creating new work history.")
      # No history, add it to work history file
      nHyperVoxels = len(self.hyperVoxels)
      minusone = -1*np.ones((nHyperVoxels,),dtype=np.int32)
      self.workHistory.create_dataset("completed",data=minusone)

      nCompleted = 0
      voxelOverflowCounter = 0

      # Could not rewrite scalars, so saving nCompleted as a vector of length 1
      self.workHistory.create_dataset("nCompleted",data=np.zeros(1,))
      allHyperIDs = np.array([x for x in self.hyperVoxels.keys()],dtype=np.int)

      # Remove the empty hyper IDs
      (validHyperIDs,emptyHyperIDs) = self.removeEmpty(allHyperIDs)
      allHyperIDs = validHyperIDs
      remaining = self.sortRemainingBySize(allHyperIDs)

      # This should never happen,
      try:
        assert (np.array([self.hyperVoxels[x]["neuronCtr"] for x in emptyHyperIDs]) == 0).all(), "All hyperIDs marked as empty are not empty!"
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)
        import pdb
        pdb.set_trace()
      
      self.writeLog("Skipping " + str(len(emptyHyperIDs)) + " empty hyper voxels")
      
      self.workHistory.create_dataset("allHyperIDs", data=allHyperIDs)
      self.workHistory.create_dataset("nSynapses",
                              data=np.zeros(nHyperVoxels,),dtype=np.int)
      self.workHistory.create_dataset("nGapJunctions",
                              data=np.zeros(nHyperVoxels,),dtype=np.int)
      self.workHistory.create_dataset("voxelOverflowCounter",
                              data=np.zeros(nHyperVoxels,),dtype=np.int)

      
    return (allHyperIDs,nCompleted,remaining,voxelOverflowCounter)

  ############################################################################

  # We want to do the hyper voxels with most neurons first, to minimize
  # the time waiting for lone cpu worker stragglers at the end.

  def sortRemainingBySize(self,remaining):

    remaining = np.array(list(remaining),dtype=int)
    
    # Minus since we want them in descending order
    numNeurons = [-self.hyperVoxels[x]["neuronCtr"] for x in remaining]
    sortIdx = np.argsort(numNeurons)
    
    return remaining[sortIdx]

  ############################################################################

  def removeEmpty(self,hyperID):

    nNeurons = np.array([self.hyperVoxels[x]["neuronCtr"] for x in hyperID])
    keepIdx = np.where(nNeurons > 0)[0]
    removeIdx = np.where(nNeurons == 0)[0]

    return (hyperID[keepIdx],hyperID[removeIdx])
  
  
  ############################################################################

  def getNeuronDistributionHistory(self):

    if("hyperVoxels" in self.workHistory):
      self.writeLog("Using neuron distribution from work history.")
      
      # We have hyper voxel information, load it
      hyperVoxels = dict([])

      for hIDstr in self.workHistory["hyperVoxels"]:
        hID = int(hIDstr)
        
        hyperVoxels[hID] = dict([])
        
        hyperVoxels[hID]["neurons"] = \
          self.workHistory["hyperVoxels"][hIDstr]["neurons"][()]
        
        hyperVoxels[hID]["neuronCtr"] = \
          self.workHistory["hyperVoxels"][hIDstr]["neuronCtr"][()]
        
        hyperVoxels[hID]["origo"] = \
          self.workHistory["hyperVoxels"][hIDstr]["origo"][()]

      hyperVoxelIDs = self.workHistory["meta/hyperVoxelIDs"][()]
      nHyperVoxels = self.workHistory["meta/nHyperVoxels"][()]
      simulationOrigo = self.workHistory["meta/simulationOrigo"][()]
      
      return (hyperVoxels,hyperVoxelIDs,nHyperVoxels,simulationOrigo)
    else:
      # No information stored
      return (None,None,None,None)

  ############################################################################

  def saveNeuronDistributionHistory(self,hyperVoxels,minCoord,maxCoord):

    self.writeLog("Writing neuron distribution history to file")
    
    assert "hyperVoxels" not in self.workHistory, \
      "saveNeuronDistributionHistory should only be called once"

    self.workHistory.create_dataset("meta/hyperVoxelIDs",
                                    data=self.hyperVoxelIDs)
    self.workHistory.create_dataset("meta/nHyperVoxels",
                                    data=self.nHyperVoxels)
    self.workHistory.create_dataset("meta/simulationOrigo",
                                    data=self.simulationOrigo)
    
    hv = self.workHistory.create_group("hyperVoxels")
    
    for hID in hyperVoxels:
      hData = hv.create_group(str(hID))
      neurons = hyperVoxels[hID]["neurons"]
      neuronCtr = hyperVoxels[hID]["neuronCtr"]
      origo = hyperVoxels[hID]["origo"]
      hData.create_dataset("neurons",data=neurons[:neuronCtr])
      hData.create_dataset("neuronCtr",data=neuronCtr)
      hData.create_dataset("origo",data=origo)
    
  ############################################################################

  def updateProcessHyperVoxelState(self,hyperID,nSyn,nGJ,execTime,
                                   voxelOverflowCounter):

    nCompleted = self.workHistory["nCompleted"][0]

    self.workHistory["completed"][nCompleted] = hyperID
    self.workHistory["nSynapses"][nCompleted] = nSyn
    self.workHistory["nGapJunctions"][nCompleted] = nGJ    
    self.workHistory["voxelOverflowCounter"][nCompleted] = voxelOverflowCounter

    nCompleted += 1
    self.workHistory["nCompleted"][0] = nCompleted
  
  ############################################################################
  
  def setupHyperVoxel(self,hyperVoxelOrigo, hyperVoxelID):

    # hypervoxel = a set of NxNxN voxels
    # hyperVoxelSynapses = list of all synapses detected in the hypervoxel

    self.hyperVoxelCoords[hyperVoxelID] = hyperVoxelOrigo # Used???
    
    self.hyperVoxelOrigo = hyperVoxelOrigo
    self.hyperVoxelID = hyperVoxelID

    if(self.hyperVoxelSynapses is None):
      self.hyperVoxelSynapses = np.zeros((self.maxSynapses,13),dtype=np.int32)
      self.hyperVoxelSynapseCtr = 0
    else:
      self.hyperVoxelSynapses[:] = 0
      self.hyperVoxelSynapseCtr = 0

    if(self.hyperVoxelGapJunctions is None):
      self.hyperVoxelGapJunctions = np.zeros((self.maxSynapses,11),
                                             dtype=np.int32)
      self.hyperVoxelGapJunctionCtr = 0
    else:
      self.hyperVoxelGapJunctions[:] = 0
      self.hyperVoxelGapJunctionCtr = 0

    # Clear lookup tables, just to be safe
    self.hyperVoxelSynapseLookup = None
    self.hyperVoxelGapJunctionLookup = None

    # Used by plotHyperVoxel to make sure synapses are displayed correctly
    self.hyperVoxelOffset = None
    
    # Which axons populate the different voxels  
    if(self.axonVoxels is None):
      self.axonVoxels = np.zeros((self.nBins[0],
                                  self.nBins[1],
                                  self.nBins[2],
                                  self.maxAxon),
                                 dtype=np.int32)
      self.axonVoxelCtr = np.zeros(self.nBins,dtype=np.int32)

      # How far from the soma is this point
      self.axonSomaDist = np.zeros((self.nBins[0],
                                    self.nBins[1],
                                    self.nBins[2],
                                    self.maxAxon),
                                   dtype=np.int16)
    else:
      # Already allocated, just clear it
      self.axonVoxels[:] = 0
      self.axonVoxelCtr[:] = 0
      self.axonSomaDist[:] = 0

    # Which dendrites populate the different voxels
    if(self.dendVoxels is None):
      self.dendVoxels = np.zeros((self.nBins[0],
                                  self.nBins[1],
                                  self.nBins[2],
                                  self.maxDend),
                                 dtype=np.int32)
      self.dendVoxelCtr = np.zeros(self.nBins,dtype=np.int32)

      # Which segment ID does the point belong to, and what segX
      self.dendSecID = np.zeros((self.nBins[0],
                                 self.nBins[1],
                                 self.nBins[2],
                                 self.maxDend),
                                dtype=np.int16)
      self.dendSecX = np.zeros((self.nBins[0],
                                self.nBins[1],
                                self.nBins[2],
                                self.maxDend),
                               dtype=np.float16) # 0 - 1.0, low pres

      # How far from the soma is this point
      self.dendSomaDist = np.zeros((self.nBins[0],
                                    self.nBins[1],
                                    self.nBins[2],
                                    self.maxDend),
                                   dtype=np.int16)
      
    else:
      # Already allocated, just clear it
      self.dendVoxels[:] = 0
      self.dendVoxelCtr[:] = 0
      self.dendSecID[:] = 0
      self.dendSecX[:] = 0
      self.dendSomaDist[:] = 0

    self.voxelOverflowCounter = 0
      
  ############################################################################

  # hyperID is only needed if we have neurons without axons, ie we use
  # axon density
  
  def detectSynapses(self):

    startTime = timeit.default_timer()
    
    # assert self.hyperVoxelSynapseCtr == 0 \
    #   and self.hyperVoxelSynapses is not None, \
    #   "setupHyperVoxel must be called before detecting synapses"

    # Find all voxels that contain axon and dendrites
    [xSyn,ySyn,zSyn] = np.where(np.bitwise_and(self.dendVoxelCtr>0,
                                               self.axonVoxelCtr>0))

    if(True):
      # This gives us some statistics, turn off later for speed
      self.maxAxonVoxelCtr = np.amax(self.axonVoxelCtr)
      self.maxDendVoxelCtr = np.amax(self.dendVoxelCtr)
    
    for x,y,z in zip(xSyn,ySyn,zSyn):
      axonIDs = self.axonVoxels[x,y,z,:self.axonVoxelCtr[x,y,z]]
      dendIDs = self.dendVoxels[x,y,z,:self.dendVoxelCtr[x,y,z]]

      axonDist = self.axonSomaDist[x,y,z,:self.axonVoxelCtr[x,y,z]]
      dendDist = self.dendSomaDist[x,y,z,:self.dendVoxelCtr[x,y,z]]
      
      dendSecID = self.dendSecID[x,y,z,:self.dendVoxelCtr[x,y,z]]
      dendSecX = self.dendSecX[x,y,z,:self.dendVoxelCtr[x,y,z]]

      # Maybe make dendrite loop outer, since it has more variables?
      # speedup??
      for (axID,axDist) in zip(axonIDs,axonDist):
        for (dID,dSegID,dSegX,dDist) \
            in zip(dendIDs,dendSecID,dendSecX,dendDist):
          
          if(axID == dID):
            # Avoid self connections
            continue
          
          preType = self.neurons[axID]["type"]
          postType = self.neurons[dID]["type"]

          if (preType,postType) in self.connectivityDistributions:
            
            conDict = self.connectivityDistributions[preType,postType]

            # We need to loop over conDict in case there are multiple
            # types of synapses from this neuron
            for conType in conDict:
              if conType == "GapJunction":
                # This part detects only axon-dend synapses, skip gap junctions
                continue
              
              meanSynapseCond,stdSynapseCond = conDict[conType]["conductance"]
              channelModelID = conDict[conType]["channelModelID"]
              
              # We can not do pruning at this stage, since we only see
              # synapses within hyper voxel, and pruning depends on
              # all synapses between two connected cells.

              # Do we have enough space allocated?
              if(self.hyperVoxelSynapseCtr >= self.maxSynapses):
                self.resizeHyperVoxelSynapsesMatrix()

              try:

                # Synapse conductance varies between synapses
                cond = np.random.normal(meanSynapseCond,
                                        stdSynapseCond)

                # Need to make sure the conductance is not negative,
                # set lower cap at 10% of mean value
                cond = np.maximum(cond,meanSynapseCond*0.1)
              
                paramID = np.random.randint(1000000)
              
                # Add synapse
                self.hyperVoxelSynapses[self.hyperVoxelSynapseCtr,:] = \
                    [axID,dID,x,y,z,self.hyperVoxelID,channelModelID,
                     axDist,dDist,dSegID,dSegX*1000,cond*1e12,paramID]
              except:
                import traceback
                tstr = traceback.format_exc()
                self.writeLog(tstr)
                import pdb
                pdb.set_trace()

              # !!! OBS, dSegX is a value between 0 and 1, multiplied by 1000
              # need to divide by 1000 later
            
              self.hyperVoxelSynapseCtr += 1

    # Sort the synapses (note sortIdx will not contain the empty rows
    # at the end.

    self.sortSynapses()

    # Convert from hyper voxel local coordinates to simulation coordinates
    # basically how many voxel steps do we need to take to go from
    # simulationOrigo to hyperVoxelOrigo (those were not included, so add them)
    try:
      hyperVoxelOffset = np.round((self.hyperVoxelOrigo - self.simulationOrigo)\
                                  /self.hyperVoxelWidth).astype(int) \
                                  * self.hyperVoxelSize

      # Just a double check...
      assert self.hyperVoxelIDs[int(np.round(hyperVoxelOffset[0]/self.hyperVoxelSize))][int(np.round(hyperVoxelOffset[1]/self.hyperVoxelSize))][int(np.round(hyperVoxelOffset[2]/self.hyperVoxelSize))] == self.hyperVoxelID, \
        "Internal inconsistency, have hyper voxel numbering or coordinates been changed?"
      
      self.hyperVoxelSynapses[:self.hyperVoxelSynapseCtr,:][:,range(2,5)]\
        += hyperVoxelOffset

      # We need this in case plotHyperVoxel is called
      self.hyperVoxelOffset = hyperVoxelOffset
      
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)

      print("ARgh, what happened...")
      import pdb
      pdb.set_trace()
      
      
    # These are used when doing the heap sort of the hyper voxels    
    self.hyperVoxelSynapseLookup \
      = self.createLookupTable(data=self.hyperVoxelSynapses,
                               nRows=self.hyperVoxelSynapseCtr)
    
    #if(self.hyperVoxelSynapseCtr > 0 and self.hyperVoxelSynapseCtr < 10):
    #  self.plotHyperVoxel()
    #  import pdb
    #  pdb.set_trace()

    endTime = timeit.default_timer()
    
    self.writeLog("detectSynapses: " + str(self.hyperVoxelSynapseCtr) \
                  + " took " + str(endTime-startTime) + "s")

    if(False and self.hyperVoxelSynapseCtr > 0):
      print("First plot shows dendrites, and the voxels that were marked")
      print("Second plot same, but for axons")
      self.plotHyperVoxel(plotNeurons=True,drawAxons=False)
      self.plotHyperVoxel(plotNeurons=True,drawDendrites=False)
      
      import pdb
      pdb.set_trace()
    
    
    return self.hyperVoxelSynapses[:self.hyperVoxelSynapseCtr,:]

  ############################################################################
  
  def placeSynapsesNoAxon(self,hyperID,voxelSpace,voxelSpaceCtr,
                                voxelAxonDist):
    
    startTime = timeit.default_timer()    

    # 1. Find neurons within hyper voxel that have no axon

    nNeurons = self.hyperVoxels[hyperID]["neuronCtr"]
    hNeurons = self.hyperVoxels[hyperID]["neurons"][:nNeurons]

    noAxonNeurons = [self.neurons[x] for x in hNeurons
                     if self.neurons[x]["axonDensity"] is not None]
    
    if(len(noAxonNeurons) == 0):
      # No neurons without axons
      return

    for naNeuron in noAxonNeurons:

      # There are two types of axon density specified
      # - Spherically symmetric
      # - f(x,y,z) in SWC coordinates

      if(naNeuron["axonDensityType"] == "r"):
      
        # 2. Check that we have cumulative probability distribution for
        #    radial distance, if not compute and cache
      
        if(naNeuron["type"] in self.axonCumDensityCache):
          (naCumDensity,naPoints) = self.axonCumDensityCache[naNeuron["type"]]

          self.writeLog("Placing " + str(naPoints) + " random axon points for "\
                        + str(naNeuron["name"]) + "(cached)")
        
        else:
          radius = np.arange(0,
                               naNeuron["axonDensityRadius"]+self.voxelSize,
                               self.voxelSize)
        
          density_as_func = eval('lambda r: ' + naNeuron["axonDensity"])
          naPDensity = np.array([density_as_func(r) for r in radius])

          # We need to scale by distance squared, since in the shell at distance
          # d further from the soma has more voxels in it than a shell closer
          # This cumulative distribution is only used to determine how far
          # from the soma a synapse is located (not the direction)

          # !!! Plot and verify this !!!
          naCumDensity = np.cumsum(np.multiply(naPDensity,radius**2))
          naCumDensity /= naCumDensity[-1] # Normalise

          # 3. Calculate how many points there should be within volume
          #    based on (unscaled raw) probability density
          # Volume at each distance is 4*pi*(r**2) * voxelSize
          naPoints = int(np.round(np.sum(4*np.pi*self.voxelSize \
                                         *np.multiply(radius**2,naPDensity))))

          self.writeLog("Placing " + str(naPoints) + " random axon points for "\
                        + str(naNeuron["name"]))
        
          self.axonCumDensityCache[naNeuron["type"]] = (naCumDensity,naPoints)

        #print("Check naPoints")
        #import pdb
        #pdb.set_trace()
        
        # 4. Randomize the points
        (naVoxelCoords,naAxonDist) = \
          self.noAxonPointsSphere(naNeuron["position"],
                                  naCumDensity,
                                  naPoints)

      elif(naNeuron["axonDensityType"] == "xyz"):

        axonDensityFunc = eval("lambda x,y,z: " + naNeuron["axonDensity"])
        
        (naVoxelCoords,naAxonDist) = \
          self.noAxonPointsXYZ(naNeuron["position"],
                               naNeuron["rotation"],
                               axonDensityFunc,
                               naNeuron["axonDensityBoundsXYZ"])
      else:
        self.writeLog("Unknown axonDensityType: " \
                      + str(naNeuron["axonDensityType"]) \
                      + "\n" + str(naNeuron))
        
        
      neuronID = naNeuron["neuronID"]
       
      for idx in range(0,naVoxelCoords.shape[0]):
        try:
          xIdx = naVoxelCoords[idx,0]
          yIdx = naVoxelCoords[idx,1]
          zIdx = naVoxelCoords[idx,2]
          axonDist = naAxonDist[idx]
        except:
          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)
         
          import pdb
          pdb.set_trace()
         
        vCtr = voxelSpaceCtr[xIdx,yIdx,zIdx]
        if(vCtr > 0 and voxelSpace[xIdx,yIdx,zIdx,vCtr-1] == neuronID):
          # Voxel already has neuronID, skip
          continue

        try:
          voxelSpace[xIdx,yIdx,zIdx,vCtr] = neuronID
          voxelAxonDist[xIdx,yIdx,zIdx,vCtr] = axonDist
          voxelSpaceCtr[xIdx,yIdx,zIdx] += 1
        except:
          self.voxelOverflowCounter += 1
          self.writeLog("!!! Axon voxel space overflow: " \
                        + str(voxelSpaceCtr[xIdx,yIdx,zIdx]))

      #if(True):
      #  # Debug plot
      #  self.plotHyperVoxel()
      #  import pdb
      #  pdb.set_trace()
    
    endTime = timeit.default_timer()

    self.writeLog("placeSynapsesNoAxonSphere: " + str(endTime-startTime) + "s" \
                  + ", hyperID: " + str(hyperID))


  
  ############################################################################
  
  # This picks points around soma centre. nPoints are randomized, points
  # outside the hyper sphere are rejected, so fewer than nPoints might be
  # returned.
  
  def noAxonPointsSphere(self,somaCentre,rCumDistribution,nPoints):
    
    uvr = np.random.rand(nPoints,3)
    theta = 2*np.pi*uvr[:,0]
    phi = np.arccos(2*uvr[:,1]-1)

    # Double check these are sorted
    # We want to sample from the supplied distance distribution
    rP = np.sort(uvr[:,2]*rCumDistribution[-1],axis=0)
    nextIdx = 0

    print("nPoints = " + str(nPoints))
    
    r = np.zeros((nPoints,))
    
    for posIdx,rP1 in enumerate(rP):
      try:
        while(rP1 > rCumDistribution[nextIdx+1]):
          nextIdx += 1

        r[posIdx] = nextIdx

      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)
        import pdb
        pdb.set_trace()
        
    # Rescale index to radie
    r = r * self.voxelSize
    
    xyz = np.array([r*np.sin(theta)*np.cos(phi),
                    r*np.sin(theta)*np.sin(phi),
                    r*np.cos(theta)]).transpose() + somaCentre

    # Check which points are inside this hyper voxel
    voxIdx = np.floor((xyz - self.hyperVoxelOrigo)/self.voxelSize).astype(int)

    #print("Verify that the points are on a sphere, and that inside check is ok")
    #import pdb
    #pdb.set_trace()
    
    insideIdx = np.where(np.sum(np.bitwise_and(0 <= voxIdx,
                                               voxIdx < self.hyperVoxelSize),
                                axis=1) == 3)[0]

    return (voxIdx[insideIdx,:],r[insideIdx])
  

  ############################################################################

  # Helper function to give points inside axon bounding box, that are
  # inside hyper voxel

  def getHVaxonPoints(self,
                      neuronPosition,
                      rotation,
                      axonDensityBoundsXYZ,
                      nPoints=1000):

    # Randomly place nPoints inside bounding box (SWC coordintes, soma (0,0,0))
    xMin = axonDensityBoundsXYZ[0]
    xWidth = axonDensityBoundsXYZ[1] - axonDensityBoundsXYZ[0]
    yMin = axonDensityBoundsXYZ[2]
    yWidth = axonDensityBoundsXYZ[3] - axonDensityBoundsXYZ[2]
    zMin = axonDensityBoundsXYZ[4]
    zWidth = axonDensityBoundsXYZ[5] - axonDensityBoundsXYZ[4]
    
    xyz = np.random.rand(nPoints,3)
    xyz[:,0] = xMin + xWidth * xyz[:,0]
    xyz[:,1] = yMin + yWidth * xyz[:,1]
    xyz[:,2] = zMin + zWidth * xyz[:,2]    
    
    # Check which of the points are inside hyper voxel (rotate+translate)
    voxIdx = ((np.matmul(rotation,xyz.transpose()).transpose() \
               + neuronPosition - self.hyperVoxelOrigo) \
              / self.voxelSize).astype(int)

    insideIdx = np.where(np.sum(np.bitwise_and(0 <= voxIdx,
                                               voxIdx < self.hyperVoxelSize),
                                axis=1) == 3)[0]

    return (xyz[insideIdx,:],voxIdx[insideIdx,:])
  
  ############################################################################

  # somaCentre and rotation of neuron
  # axonDensityFunc should be written so that it can handle x,y,z (SWC coords)
  # as vectors
  # axonDensityBoundsXYZ = [xmin,xmax,ymin,ymax,zmin,zmax] in SWC coordinates
  
  # axonDensityFunc = eval("lambda x,y,z: " + axonPstr)
  
  def noAxonPointsXYZ(self, neuronPosition,rotation,
                      axonDensityFunc,axonDensityBoundsXYZ):

    
    # Points for initial sample
    nPoints = 5000

    (xyzInside,voxIdx) = self.getHVaxonPoints(neuronPosition,
                                              rotation,
                                              axonDensityBoundsXYZ,
                                              nPoints)

    xWidth = axonDensityBoundsXYZ[1] - axonDensityBoundsXYZ[0]
    yWidth = axonDensityBoundsXYZ[3] - axonDensityBoundsXYZ[2]
    zWidth = axonDensityBoundsXYZ[5] - axonDensityBoundsXYZ[4]

    pointVolume = xWidth*yWidth*zWidth / nPoints
    voxelVolume = self.voxelSize ** 3
    
    # If no points were inside HV, then the intersection must be small
    # so we assume no voxels should be filled
    if(xyzInside.shape[0] == 0):
      # Double check correct data-types
      self.writeLog("Bounding box appears to be outside hyper voxel")
      return (np.zeros((0,3),dtype=int),np.zeros((0,1)))      
    
    # Calculate density at each of the points inside HV
    densityInside = axonDensityFunc(xyzInside[:,0],
                                    xyzInside[:,1],
                                    xyzInside[:,2])

    # Estimate number of synapses from density, in this case we use a volume
    # equal to bounding box volume / nPoints for each point.
    # OBS that we only want to know how many synapses to place inside HV
    nPointsToPlace = np.round(np.sum(densityInside*pointVolume))
    
    if(nPointsToPlace <= 0):
      # To little probability mass inside
      self.writeLog("Too little of axon appears to be inside hyper voxel")
      return (np.zeros((0,3),dtype=int),np.zeros((0,1)))
    
    # Calculate max density inside HV, divide all densities by that to
    # get Pkeep.
    maxDensity = np.max(densityInside)

    # We know that n out of N points placed were inside volume, so volume
    # acceptance rate is n/N.
    # In order to get about nPointsToPlace points placed we need to account
    # for how many outside volume, and also account for how many of those
    # inside volume gets accepted (Pkeep = density / maxDensity)
    nTries = np.round(nPointsToPlace * nPoints \
                      / np.sum(densityInside/maxDensity)).astype(int)

    if(nTries > 1e6):
      self.writeLog("!!! noAxonPointsXYZ: Warning trying to place " \
                    + str(nTries) + " points. Bounds: " \
                    + str(axonDensityBoundsXYZ))

    # Only print this in verbose mode
    if(self.verbose):
      self.writeLog("Trying to place " + str(nTries) + " points to get " \
                    + str(nPointsToPlace) + " axon voxels filled")

    if(nPoints >= nTries):
      # We already have enough points, use a subset
      useNum = np.round(voxIdx.shape[0]*nTries/nPoints).astype(int)
      pickedIdx = np.where(np.random.rand(useNum) \
                           < densityInside[:useNum] / maxDensity)[0]
      axonDist = np.sqrt(np.sum((xyzInside[pickedIdx,:])**2,axis=1))

      return (voxIdx[pickedIdx,:],axonDist)
    else:
      # Not enough points, use the ones we have, then get more
      pickedIdx = np.where(np.random.rand(voxIdx.shape[0])
                           < densityInside/maxDensity)[0]
      axonDist = np.sqrt(np.sum((xyzInside[pickedIdx,:])**2,axis=1))

      # Get more points
      
      (xyzInsideB,voxIdxB) = \
        self.getHVaxonPoints(neuronPosition,
                             rotation,
                             axonDensityBoundsXYZ,
                             nTries-nPoints)
      
      densityInsideB = axonDensityFunc(xyzInsideB[:,0],
                                       xyzInsideB[:,1],
                                       xyzInsideB[:,2])

      pickedIdxB = np.where(np.random.rand(voxIdxB.shape[0]) \
                            < densityInsideB/maxDensity)[0]
      
      axonDistB = np.sqrt(np.sum((xyzInsideB[pickedIdxB,:])**2,axis=1))

      return (np.concatenate([voxIdx[pickedIdx,:],
                              voxIdxB[pickedIdxB,:]]),
              np.concatenate([axonDist,axonDistB]))
    

  ############################################################################


  def resizeHyperVoxelSynapsesMatrix(self,newSize=None):                 

    if(newSize is None):
      newSize = int(np.ceil(1.5*self.maxSynapses))
    
    assert newSize >= self.hyperVoxelSynapseCtr, \
      " Can not shrink below existing number of synapses"
                 
    # We need to increase matrix size
    old = self.hyperVoxelSynapses
    self.maxSynapses = newSize
    self.writeLog("!!! Increasing max synapses to " + str(self.maxSynapses))
    self.hyperVoxelSynapses = np.zeros((self.maxSynapses,13),
                                       dtype=np.int32)
    self.hyperVoxelSynapses[:old.shape[0],:] = old
    del old
                     
                 
  ############################################################################

  # This truncates and sorts the hyper voxel synapse matrix
  
  def sortSynapses(self):

    sortIdx = np.lexsort(self.hyperVoxelSynapses[:self.hyperVoxelSynapseCtr,
                                                 [0, 1]].transpose())

    self.hyperVoxelSynapses[:self.hyperVoxelSynapseCtr,:] = \
      self.hyperVoxelSynapses[sortIdx,:]

  ############################################################################

  def sortGapJunctions(self):

    sortIdx = \
      np.lexsort(self.hyperVoxelGapJunctions[:self.hyperVoxelGapJunctionCtr,
                                             [0, 1]].transpose())

    self.hyperVoxelGapJunctions[:self.hyperVoxelGapJunctionCtr,:] = \
      self.hyperVoxelGapJunctions[sortIdx,:]

  ############################################################################

  # First and second column specify the source and destination ID of a synapse
  # or a gap junction.
  #
  # This creates a lookup table where all synapses in the hyper voxel
  # between the same pair of neurons are grouped together
  # returns a matrix where first column is a UID = srcID*nNeurons + destID
  # and the following two columns are start row and end row (-1) in matrix
  
  def createLookupTable(self,data,nRows):

    self.writeLog("Create lookup table")
    # nRows = data.shape[0] -- zero padded, cant use shape
    lookupTable = np.zeros((data.shape[0],3),dtype=int)
    
    nextIdx = 0
    startIdx = 0

    lookupIdx = 0
    nNeurons = len(self.neurons)
    
    while(nextIdx < nRows):
      srcID = data[nextIdx,0]
      destID = data[nextIdx,1]
      nextIdx += 1

      while(nextIdx < nRows \
            and data[nextIdx,0] == srcID \
            and data[nextIdx,1] == destID):
        nextIdx += 1

      lookupTable[lookupIdx,:] = [destID*nNeurons+srcID,startIdx,nextIdx]
                                  
      startIdx = nextIdx
      lookupIdx += 1
      
    return lookupTable[:lookupIdx,:]
    
  ############################################################################

  def includesGapJunctions(self):

    hasGapJunctions = False
    
    for key in self.connectivityDistributions:
      if("GapJunction" in self.connectivityDistributions[key]):
        hasGapJunctions = True
        
    return hasGapJunctions
  
  ############################################################################
  
  # Gap junctions are stored in self.hyperVoxelGapJunctions
    
  def detectGapJunctions(self):

    if(not self.includesGapJunctions()):
      self.writeLog("detectGapJunctions: No gap junctions defined in connectivity rules")
      return
    
    startTime = timeit.default_timer()
    
    assert self.hyperVoxelGapJunctionCtr == 0 \
      and self.hyperVoxelGapJunctions is not None, \
      "setupHyperVoxel must be called before detecting gap junctions"
    
    [xDV,yDV,zDV] = np.where(self.dendVoxelCtr > 0)
    
    for x,y,z in zip(xDV,yDV,zDV):

      neuronIDs = self.dendVoxels[x,y,z,:self.dendVoxelCtr[x,y,z]]
      segID = self.dendSecID[x,y,z,:self.dendVoxelCtr[x,y,z]]
      segX = self.dendSecX[x,y,z,:self.dendVoxelCtr[x,y,z]]

      # All possible pairs
      for pairs in itertools.combinations(np.arange(0,self.dendVoxelCtr[x,y,z]),2):
        neuronID1 = self.dendVoxels[x,y,z,pairs[0]]
        neuronID2 = self.dendVoxels[x,y,z,pairs[1]]

        # !!! Check no self connections??
        
        # Add type field, derived from name field MSD1_45 --> MSD1
        preType = self.neurons[neuronID1]["type"]
        postType = self.neurons[neuronID2]["type"]
        
        if (preType,postType) in self.connectivityDistributions:
          
          if("GapJunction" \
             in self.connectivityDistributions[preType,postType]):

            conInfo \
              = self.connectivityDistributions[preType,postType]["GapJunction"]
            
            segID1 = self.dendSecID[x,y,z,pairs[0]]
            segID2 = self.dendSecID[x,y,z,pairs[1]]        
        
            segX1 = self.dendSecX[x,y,z,pairs[0]]
            segX2 = self.dendSecX[x,y,z,pairs[1]]        

            meanGJcond,stdGJcond = conInfo["conductance"]

            # !!! Currently not using channelParamDict for GJ
          
            GJcond = np.random.normal(meanGJcond,stdGJcond)
            GJcond = np.maximum(GJcond,meanGJcond*0.1) # Avoid negative cond
          
            self.hyperVoxelGapJunctions[self.hyperVoxelGapJunctionCtr,:] = \
              [neuronID1,neuronID2,segID1,segID2,segX1*1e3,segX2*1e3,
               x,y,z,self.hyperVoxelID,GJcond*1e12]
            self.hyperVoxelGapJunctionCtr += 1

    self.sortGapJunctions()

    # We also translate from local hyper voxel coordinates to simulation
    # voxel coordinates

    hyperVoxelOffset = np.round((self.hyperVoxelOrigo \
                                 -self.simulationOrigo) \
                                /self.hyperVoxelWidth).astype(int) \
                                * self.hyperVoxelSize
    
    self.hyperVoxelGapJunctions[:self.hyperVoxelGapJunctionCtr,:][:,range(6,9)]\
      += hyperVoxelOffset

    
    self.hyperVoxelGapJunctionLookup \
      = self.createLookupTable(data=self.hyperVoxelGapJunctions,
                               nRows=self.hyperVoxelGapJunctionCtr)
        
    endTime = timeit.default_timer()

    self.writeLog("detectGapJunctions: " + str(endTime-startTime) + "s")
    
    return self.hyperVoxelGapJunctions[:self.hyperVoxelGapJunctionCtr,:]
          
  ############################################################################

  def setupLog(self,logFileName=None):

    if(logFileName is None):
      logFileName = self.logFileName

    if(logFileName is None or len(logFileName) == 0):
      # Not a valid log file name
      return
      
    if(self.logFile is not None):
      self.writeLog("Already have a log file setup, ignoring")
      return

    self.logFile = open(logFileName,'wt')    

  ############################################################################
  
  def writeLog(self,text,flush=True): # Change flush to False in future, debug
    if(self.logFile is not None):
      self.logFile.write(text + "\n")
      print(text)
      if(flush):
        self.logFile.flush()
    else:
      if(self.verbose):
        print(text)
    
  ############################################################################

  def readPrototypes(self, configFile=None, axonStumpIDFlag=False):

    if(configFile is None):
      configFile = self.configFile

    configFile = self.getPath(configFile)
      
    self.axonStumpIDFlag = axonStumpIDFlag
      
    print("Loading from " + configFile)
    
    cfgFile = open(str(configFile),'r')
    
    try:
      self.config = json.load(cfgFile)
    finally:
      cfgFile.close()

    self.prototypeNeurons = dict()
      
    for name, definition in self.config["Neurons"].items():

      self.writeLog("Reading prototype for: " + name)
      
      morph = self.getPath(definition["morphology"])
      param = self.getPath(definition["parameters"])
      mech = self.getPath(definition["mechanisms"])

      if("neuronType" in definition):
        neuronType = definition["neuronType"]
      else:
        neuronType = "neuron"
        
      if(neuronType == "virtual"):
        virtualNeuron = True
      else:
        virtualNeuron = False
      
      
      if 'hoc' in definition:
        hoc = definition["hoc"]
      else:
        hoc = None
        
      self.prototypeNeurons[name] \
        = NeuronMorphology(name=name,
                           swc_filename=morph,
                           param_filename=param,
                           mech_filename=mech,
                           hoc=hoc,
                           virtualNeuron=virtualNeuron,
                           axonStumpIDFlag=axonStumpIDFlag)
      
      if( "axonDensity" in definition):
        self.writeLog("Setting axon density for neuron without axon")
        axonDensityType = definition["axonDensity"][0]

        if(axonDensityType == "r"):
          density = definition["axonDensity"][1]
          maxRadius = definition["axonDensity"][2]

          self.prototypeNeurons[name].setAxonVoxelRadialDensity(density,
                                                                maxRadius)
        elif(axonDensityType == "xyz"):
          density = definition["axonDensity"][1]
          axonDensityBoundsXYZ = np.array(definition["axonDensity"][2])

          self.prototypeNeurons[name].setAxonVoxelXYZDensity(density,
                                                          axonDensityBoundsXYZ)
          
        else:
          self.writeLog(str(name) + ": Unknown axon density type : " \
                        + str(axonDensityType) \
                        + "\n" + str(definition["axonDensity"]))
        
      else:
        # If no axon density specified, then axon must be present in morphology
        assert (len(self.prototypeNeurons[name].axon) > 0), \
          "File: " + morph + " does not have an axon"
                
      assert len(self.prototypeNeurons[name].dend) > 0 \
        or self.prototypeNeurons[name].virtualNeuron,\
        "File: " + morph + " does not have a dendrite"
        
      # Since we already have the config file open, let's read connectivity
      # distributions also

    self.writeLog("Loading connectivity information")
      
    for name,definition in self.config["Connectivity"].items():
     
      preType,postType = name.split(",")

      conDef = definition.copy()
     
      for key in conDef:
        if(key == "GapJunction"):
          conDef[key]["channelModelID"] = 3
        else:
          conDef[key]["channelModelID"] = self.nextChannelModelID
          self.nextChannelModelID += 1

        # Also if conductance is just a number, add std 0
        if(type(conDef[key]["conductance"]) not in [list,tuple]):
          conDef[key]["conductance"] = [conDef[key]["conductance"],0]
        
      self.connectivityDistributions[preType,postType] = conDef

  ############################################################################

  def readNeuronPositions(self,positionFile):

    if(positionFile is None):
      positionFile = self.positionFile

    mem = self.memory()
    self.writeLog(str(mem))
      
    self.writeLog("Reading positions from file: " + positionFile)

    posInfo = SnuddaLoad(positionFile).data
      
    mem = self.memory()
    self.writeLog(str(mem))

    # TEMP DISABLE CHECK, TURN BACK ON LATER
    if(False):
      # Make sure we do not change config file unintentionally
      assert posInfo["configFile"] == self.configFile, \
        "Not using original config file: " \
        + str(posInfo["configFile"]) + " vs " + self.configFile

      
    self.neurons = posInfo["neurons"]
    nNeurons = len(self.neurons)
    
    self.neuronPositions = np.zeros((nNeurons,3))
    
    for ni,neuron in enumerate(posInfo["neurons"]):
      self.neuronPositions[ni,:] = neuron["position"]
      
      # Add a few sanity checks
      assert ni == neuron["neuronID"], \
          "NeuronID=" + str(neuron["neuronID"]) + "and ni=" + str(ni) \
           + " not equal, corruption?"
      assert neuron["name"] in self.prototypeNeurons, \
         "Neuron type " + neuron["name"] + " not in prototypeNeurons: " \
          + str(self.prototypeNeurons)

    # Also load the channel data
    self.nPopulationUnits = posInfo["nPopulationUnits"]
    self.populationUnit = posInfo["populationUnit"]
    self.populationUnitPlacementMethod = posInfo["populationUnitPlacementMethod"]
      
    self.populationUnits = dict([])
    for i in range(0,self.nPopulationUnits):
      self.populationUnits[i] = np.where(self.populationUnit == i)[0]

    self.writeLog("Position file read.")
    del posInfo

  ############################################################################

  # If the detect is rerun we need to make sure there are not old MERGE
  # files left that might remember old run accidentally
  
  def deleteOldMerge(self):

    if(self.role == "master"):

      workDir = os.path.dirname(self.saveFile)
      workDir = (workDir + "/").replace("/voxels/","/")

      delFiles =  [workDir + "network-putative-synapses-MERGED.hdf5",
                   workDir + "network-putative-synapses-MERGED.hdf5-cache",
                   workDir + "network-pruned-synapses.hdf5",
                   workDir + "network-pruned-synapses.hdf5-cache"]

      for f in delFiles:
        if os.path.exists(f):
          self.writeLog("Removing old files " + str(f))
          os.remove(f)
    
  ############################################################################

  def writeHyperVoxelToHDF5(self):

    startTime = timeit.default_timer()
    
    outputName = self.saveFile.replace(".hdf5","-" + str(self.hyperVoxelID) +".hdf5")

    with h5py.File(outputName, "w", libver=self.h5libver) as outFile:

      configData = outFile.create_dataset("config", \
                                          data=json.dumps(self.config))

      metaData = outFile.create_group("meta")
      
      metaData.create_dataset("hyperVoxelID",data=self.hyperVoxelID)
      metaData.create_dataset("hyperVoxelOrigo",data=self.hyperVoxelOrigo)
      metaData.create_dataset("simulationOrigo",data=self.simulationOrigo)
     
      metaData.create_dataset("SlurmID",data=self.SlurmID)
      metaData.create_dataset("voxelSize",data=self.voxelSize)
      metaData.create_dataset("hyperVoxelSize",data=self.hyperVoxelSize)
      metaData.create_dataset("nBins",data=self.nBins)
      metaData.create_dataset("voxelOverflowCounter",
                              data=self.voxelOverflowCounter)
      
      metaData.create_dataset("configFile",data=self.configFile)
      metaData.create_dataset("positionFile",data=self.positionFile)

      metaData.create_dataset("axonStumpIDFlag",data=self.axonStumpIDFlag)
      
      # These may or may not exist, if they do, write them to file
      if(self.maxAxonVoxelCtr is not None):
        metaData.create_dataset("maxAxonVoxelCtr",data=self.maxAxonVoxelCtr)

      if(self.maxDendVoxelCtr is not None):
        metaData.create_dataset("maxDendVoxelCtr",data=self.maxDendVoxelCtr)

        
      if(self.voxelOverflowCounter > 0):
        self.writeLog("!!! Voxel overflow detected, please increase maxAxon and maxDend")
      
      networkGroup = outFile.create_group("network")
      networkGroup.create_dataset("synapses", \
                                  data=self.hyperVoxelSynapses[:self.hyperVoxelSynapseCtr,:], \
                                  dtype=np.int32, \
                                  chunks = (self.synapseChunkSize,13), \
                                  maxshape=(None,13), \
                                  compression=self.h5compression)
      networkGroup.create_dataset("gapJunctions", \
                                  data=self.hyperVoxelGapJunctions[:self.hyperVoxelGapJunctionCtr,:], \
                                  dtype=np.int32, \
                                  chunks = (self.gapJunctionChunkSize,11), \
                                  maxshape=(None,11), \
                                  compression=self.h5compression)

      networkGroup.create_dataset("synapseLookup",
                                  data=self.hyperVoxelSynapseLookup,
                                  dtype=int)

      networkGroup.create_dataset("gapJunctionLookup",
                                  data=self.hyperVoxelGapJunctionLookup,
                                  dtype=int)

      # Additional information useful for debugging
      if(self.debugFlag):
        debugGroup = outFile.create_group("debug")
        
        debugGroup.create_dataset("dendVoxels",data=self.dendVoxels)
        debugGroup.create_dataset("axonVoxels",data=self.axonVoxels)
        
        debugGroup.create_dataset("dendVoxelCtr",data=self.dendVoxelCtr)
        debugGroup.create_dataset("axonVoxelCtr",data=self.axonVoxelCtr)
        
        
      endTime = timeit.default_timer()
      
      outFile.close()

    self.writeLog("Wrote hyper voxel " + str(self.hyperVoxelID) \
                  + " (" + str(self.hyperVoxelSynapseCtr) + " synapses, "
                  + str(self.hyperVoxelGapJunctionCtr) + " gap junctions)")

  ############################################################################
  
  def loadNeuron(self, neuronInfo):

    # Clone prototype neuron (it is centred, and not rotated)
    neuron = self.prototypeNeurons[neuronInfo["name"]].clone()

    # Rotate and place neuron in correct location
    neuron.place(rotation=neuronInfo["rotation"],
                 position=neuronInfo["position"])
          
    return neuron

  ############################################################################

  def distributeNeuronsParallel(self,dView=None):

    if(self.role != "master"):
      # Only run this as master
      return

    (hyperVoxels,hyperVoxelIDs,nHyperVoxels,simulationOrigo) = \
      self.getNeuronDistributionHistory()

    # Do we have old data that we can reuse?
    if(hyperVoxels is not None):
      self.writeLog("distributeNeuronsParallel: Reusing old neuron allocation")

      self.hyperVoxels = hyperVoxels
      self.hyperVoxelIDs = hyperVoxelIDs
      self.nHyperVoxels = nHyperVoxels
      self.simulationOrigo = simulationOrigo
      
      # We need to push the data to the workers also
      dView.push({"simulationOrigo":simulationOrigo,
                  "nc.hyperVoxels":hyperVoxels,
                  "nc.hyperVoxelIDs":hyperVoxelIDs,
                  "nc.nHyperVoxels":nHyperVoxels},block=True)
      return

    # No old data, we need to calculate it
    
    if(dView is None):
      self.writeLog("No dView specified, running distribute neurons in serial")
      (minCoord,maxCoord) = self.distributeNeurons()
      
      self.saveNeuronDistributionHistory(hyperVoxels=self.hyperVoxels,
                                         minCoord=minCoord,
                                         maxCoord=maxCoord)
      
      return

    (minCoord,maxCoord) = self.findMinMaxCoordParallell(dView=dView,
                                                        volumeID=self.volumeID)
    
    neuronIdx = np.random.permutation(np.arange(0,len(self.neurons),
                                                dtype=np.int32))

    # Split the neuronIdx between the workers
    dView.scatter("neuronIdx",neuronIdx,block=True)
    dView.push({"minCoord":minCoord,
                "maxCoord":maxCoord},block=True)

    self.writeLog("Distributing neurons, parallell.")

    # For the master node, run with empty list
    # This sets up internal state of master
    self.distributeNeurons(neuronIdx = [], minCoord=minCoord,maxCoord=maxCoord)
    
    cmdStr = "nc.distributeNeurons(neuronIdx=neuronIdx," \
                                + "minCoord=minCoord," \
                                + "maxCoord=maxCoord)"
    dView.execute(cmdStr,block=True)

    self.writeLog("Gathering neuron distribution from workers")

    # Collect all the neurons in the list from the workers
    # For each neuron we found out which hyper voxels it occupies,
    # now we want for each hyper voxel to know which neurons are in there
    hyperVoxelList = dView.gather("nc.hyperVoxels",block=True)

    self.writeLog("Distributions received.")
    
    for hv in hyperVoxelList:
      for hID in hv:
        
        assert (hv[hID]["origo"] == self.hyperVoxels[hID]["origo"]).all(), \
          "Origo for hyper voxels do not match --- should never happen"
        
        nNeurons = int(hv[hID]["neuronCtr"])
        startIdx = int(self.hyperVoxels[hID]["neuronCtr"] )
        endIdx = startIdx + nNeurons
        
        if(endIdx >= len(self.hyperVoxels[hID]["neurons"])):
          # Not enough space, reallocating

          old = self.hyperVoxels[hID]["neurons"]
          newMax = endIdx + self.maxNeurons
          
          self.hyperVoxels[hID]["neurons"] = np.zeros((newMax,),dtype=np.int32)
          
          # Copying back the old data to new vector
          if(len(old) > 0):
            self.hyperVoxels[hID]["neurons"][:len(old)] = old
            
          del old
            
        # Adding the new neurons
        self.hyperVoxels[hID]["neurons"][startIdx:endIdx] = \
          hv[hID]["neurons"][:nNeurons]
            
        # Increment counter
        self.hyperVoxels[hID]["neuronCtr"] += nNeurons

    # Sorting the list of neurons.
    # -- check why this order matters to number of synapses detected,
    #    it should not matter (except in case of voxel overflows).
    if(False):
      for hID in self.hyperVoxels:
        nCtr = self.hyperVoxels[hID]["neuronCtr"]
      
        self.hyperVoxels[hID]["neurons"] = \
          np.sort(self.hyperVoxels[hID]["neurons"][:nCtr])
        
    # Distribute the new list to all neurons
    dView.push({"nc.hyperVoxels":self.hyperVoxels},block=True)

    self.saveNeuronDistributionHistory(hyperVoxels=self.hyperVoxels,
                                       minCoord=minCoord,
                                       maxCoord=maxCoord)
    
  ############################################################################
  
  # This creates a list for each hyper voxel for the neurons that
  # has any neurites within its border (here defined as vertices inside region)
  
  def distributeNeurons(self, neuronIdx=None, minCoord=None, maxCoord=None):

    if(neuronIdx is None):
      neuronIdx = np.arange(0,len(self.neurons),dtype=np.int32)
    
    self.writeLog("distributeNeurons: neuronIdx = " + str(neuronIdx) \
                  + " (n=" + str(len(neuronIdx)) + ")")
                
    startTime = timeit.default_timer()
    
    if(maxCoord is None or minCoord is None):
      self.writeLog("distributeNeurons: calculating min and max coords")
      (minCoord,maxCoord) = self.findMinMaxCoord()

    # Simulation origo is in meters
    self.simulationOrigo = minCoord
    
    assert ((self.nBins - self.nBins[0]) == 0).all(), \
      "Hyper voxels should be cubes"
    
    self.hyperVoxelWidth = self.nBins[0] * self.voxelSize

    self.nHyperVoxels = np.ceil((maxCoord - minCoord )
                                / self.hyperVoxelWidth).astype(int) + 1

    self.hyperVoxelIDs = np.zeros(self.nHyperVoxels,dtype=int)
    
    self.hyperVoxelIDs[:] = \
      np.arange(0,self.hyperVoxelIDs.size).reshape(self.hyperVoxelIDs.shape)

    self.writeLog(str(self.hyperVoxelIDs.size) + " hyper voxels in total")
    
    # First assign hyperVoxelID to the space
    
    self.hyperVoxels = dict([])
    
    for ix in range(0,self.nHyperVoxels[0]):
      for iy in range(0,self.nHyperVoxels[1]):
        for iz in range(0,self.nHyperVoxels[2]):

          hID = self.hyperVoxelIDs[ix,iy,iz]

          self.hyperVoxels[hID] = dict([])
          self.hyperVoxels[hID]["origo"] = self.simulationOrigo \
                                          + self.hyperVoxelWidth\
                                            *np.array([ix,iy,iz])
          
          # Changed so we preallocate only empty, to preserve memory
          self.hyperVoxels[hID]["neurons"] = np.zeros((0,),dtype=np.int32)
          self.hyperVoxels[hID]["neuronCtr"] = 0

    self.writeLog("Pre allocation done.")
          
    ctr = 0

    if(neuronIdx is None):
      neurons = self.neurons
    elif(len(neuronIdx) == 0):
      neurons = []
    else:
      neurons = [self.neurons[idx] for idx in neuronIdx]

    for n in neurons:
      
      ctr = ctr + 1
      if(ctr % 10000 == 0):
        print("Assignment counter: " + str(ctr))
      
      neuron = self.loadNeuron(n)
      neuronID = n["neuronID"]
      
      if(neuron.dend.shape[0] > 0):
        dendLoc = np.floor((neuron.dend[:,:3] - self.simulationOrigo) \
                           / self.hyperVoxelWidth).astype(int)
      else:
        dendLoc = np.zeros((0,3))
        
      if(neuron.axon.shape[0] > 0):
        # We have an axon, use it
        axonLoc = np.floor((neuron.axon[:,:3] - self.simulationOrigo) \
                           / self.hyperVoxelWidth).astype(int)
        
      elif(neuron.axonDensityType == "r"):
        # axonLoc = np.zeros((0,3))

        # We create a set of points corresponding approximately to the
        # extent of the axonal density, and check which hyper voxels
        # they occupy

        # Radius of sphere in hyper voxels, rounded up
        rad = np.ceil(neuron.maxAxonRadius \
                          / (self.hyperVoxelSize*self.voxelSize))
        
        # Approximately how many hyper voxels will the dendritic tree occupy
        nHV = (2 * rad ) ** 3

        # Over sample
        nPoints = int(30*nHV)

        # Randomly place these many points within a sphere of the given radius
        # and then check which hyper voxels these points belong to

        theta = 2*np.pi*np.random.rand(nPoints)
        phi = np.arccos(2*np.random.rand(nPoints)-1)
        r = neuron.maxAxonRadius * (np.random.rand(nPoints) ** (1/3))

        x = np.multiply(r, np.multiply(np.sin(phi), np.cos(theta)))
        y = np.multiply(r, np.multiply(np.sin(phi), np.sin(theta)))
        z = np.multiply(r, np.cos(phi))

        axonCloud = np.zeros((len(x),3))
        axonCloud[:,0] = x + neuron.soma[0,0]
        axonCloud[:,1] = y + neuron.soma[0,1]
        axonCloud[:,2] = z + neuron.soma[0,2]

        axonLoc = np.floor((axonCloud[:,:3] - self.simulationOrigo) \
                           / self.hyperVoxelWidth).astype(int)

        axonInsideFlag = [xa >= 0 and xa < self.hyperVoxelIDs.shape[0] \
                          and ya >= 0 and ya < self.hyperVoxelIDs.shape[1] \
                          and za >= 0 and za < self.hyperVoxelIDs.shape[2] \
                          for xa,ya,za in axonLoc]

        axonLoc = axonLoc[axonInsideFlag,:]
        
        # We need to remove the axon volumes outside the modelled volume

        if(False):
          # Verify
          import matplotlib.pyplot as plt
          from mpl_toolkits.mplot3d import Axes3D
          
          fig = plt.figure()
          ax = fig.gca(projection='3d')
          ax.scatter(axonCloud[:,0],axonCloud[:,1],axonCloud[:,2])
          plt.ion()
          plt.show()
          
          import pdb
          pdb.set_trace()
        
        # !!! If there is no axon, and we use probability density,
        #     then we need to include the neuron in hyper voxels
        #     within the axon volume specified

      elif(neuron.axonDensityType == "xyz"):
        
        # Estimate how many points we need to randomly place
        nPoints = 100 *np.prod(neuron.axonDensityBoundsXYZ[1:6:2] \
                              - neuron.axonDensityBoundsXYZ[0:6:2]) \
                              / ((self.hyperVoxelSize*self.voxelSize)**3)
        nPoints = int(np.ceil(nPoints))
        
        if(nPoints > 1e4):
          self.writeLog("!!! Many many points placed for axon density of" \
                        + str(neuron.name) + ": " + str(nPoints))

        xmin = neuron.axonDensityBoundsXYZ[0]
        xwidth = neuron.axonDensityBoundsXYZ[1]-neuron.axonDensityBoundsXYZ[0]
        ymin = neuron.axonDensityBoundsXYZ[2]
        ywidth = neuron.axonDensityBoundsXYZ[3]-neuron.axonDensityBoundsXYZ[2]
        zmin = neuron.axonDensityBoundsXYZ[4]
        zwidth = neuron.axonDensityBoundsXYZ[5]-neuron.axonDensityBoundsXYZ[4]

        # The purpose of this is to find out the range of the axon bounding box
        axonCloud = np.random.rand(nPoints,3)
        axonCloud[:,0] = axonCloud[:,0]*xwidth + xmin
        axonCloud[:,1] = axonCloud[:,1]*ywidth + ymin
        axonCloud[:,2] = axonCloud[:,2]*zwidth + zmin        

        # Dont forget to rotate        
        axonCloud = np.matmul(neuron.rotation,
                              axonCloud.transpose()).transpose() \
                              + neuron.position

        axonLoc = np.floor((axonCloud[:,:3] - self.simulationOrigo) \
                           / self.hyperVoxelWidth).astype(int)

        axonInsideFlag = [x >= 0 and x < self.hyperVoxelIDs.shape[0] \
                          and y >= 0 and y < self.hyperVoxelIDs.shape[1] \
                          and z >= 0 and z < self.hyperVoxelIDs.shape[2] \
                          for x,y,z in axonLoc]

        axonLoc = axonLoc[axonInsideFlag,:]        

        if(False):
           # Verify
          import matplotlib.pyplot as plt
          from mpl_toolkits.mplot3d import Axes3D
          
          fig = plt.figure()
          ax = fig.gca(projection='3d')
          ax.scatter(axonCloud[:,0],axonCloud[:,1],axonCloud[:,2])
          plt.ion()
          plt.show()
          
          import pdb
          pdb.set_trace()
          
        
      else:
        self.writeLog(str(neuron.name)\
                      + ": No axon and unknown axon density type: " \
                      + str(neuron.axonDensityType))
        assert False, "No axon for " + str(neuron.name)
        
      # Find unique hyper voxel coordinates
      hLoc = np.unique(np.concatenate([axonLoc,dendLoc]),axis=0).astype(int)
      
      if(n["virtualNeuron"]):        
        # Range check since we have neurons coming in from outside the volume
        # the parts outside should be ignored
        try:
          hyperID = [self.hyperVoxelIDs[x,y,z] for x,y,z in hLoc \
                     if x >= 0 and x < self.hyperVoxelIDs.shape[0] \
                     and y >= 0 and y < self.hyperVoxelIDs.shape[1] \
                     and z >= 0 and z < self.hyperVoxelIDs.shape[2]]
        except:
          self.writeLog("Hyper ID problem")
          assert False, "Hyper ID problem. x=" + str(x) + " y=" \
            + str(y)+ " z=" + str(z)
          import pdb
          pdb.set_trace()
      else:
        # Not a virtual neuron, should all be inside volume
        try:
          hyperID = [self.hyperVoxelIDs[x,y,z] for x,y,z in hLoc]
        except Exception as e:
          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)
          self.writeLog("Affected neuron: " + str(n))
          self.writeLog("Range check fuckup : x = " + str(x) \
                        + " y = " + str(y)+ " z = " + str(z))
          assert False, "Range check fuckup : x = " + str(x) + " y=" \
            + str(y)+ " z=" + str(z)
          import pdb
          pdb.set_trace()
        
      # Add the neuron to the hyper voxel's list over neurons
      for hID in hyperID:
        
        nextPos =  self.hyperVoxels[hID]["neuronCtr"]

        if(nextPos >= len(self.hyperVoxels[hID]["neurons"])):
          old = self.hyperVoxels[hID]["neurons"]
          newMax = nextPos + self.maxNeurons
          self.hyperVoxels[hID]["neurons"] = np.zeros((newMax,),dtype=np.int32)

          if(nextPos > 0):
            self.hyperVoxels[hID]["neurons"][:len(old)] = old
            
          del old
          
        self.hyperVoxels[hID]["neurons"][nextPos] = neuronID
        self.hyperVoxels[hID]["neuronCtr"] += 1

    endTime = timeit.default_timer()

    if(len(neurons) > 0):
      self.writeLog("Calculated distribution of neurons: " \
                    + str(endTime-startTime) + " seconds")

    # For serial version of code, we need to return this, so we
    # can save work history
    return (minCoord,maxCoord)

        
  ############################################################################

  def setupParallel(self,dView=None):

    assert self.role == "master", \
      "setupParallel: Should only be called by master node"

    if(dView is None):
      self.writeLog("setupParallel called without dView, aborting.")
      return
   
    if(self.workersInitialised):
      self.writeLog("Workers already initialised.")
      return
    
    with dView.sync_imports():
      from snudda.detect import SnuddaDetect

    self.writeLog("Setting up workers: " \
                  + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    
    # Create unique log file names for the workers
    if(self.logFileName is not None):
      engineLogFile = [self.logFileName + "-" \
                       + str(x) for x in range(0,len(dView))]
    else:
      engineLogFile = [[] for x in range(0,len(dView))]

    self.writeLog("Scattering " + str(engineLogFile))
      
    dView.scatter('logFileName',engineLogFile,block=True)

    self.writeLog("Scatter done.")
    
    dView.push({"positionFile":self.positionFile,
                "configFile":self.configFile,
                "voxelSize":self.voxelSize,
                "hyperVoxelSize":self.hyperVoxelSize,
                "verbose":self.verbose,
                "SlurmID":self.SlurmID,
                "saveFile":self.saveFile},
               block=True)

    self.writeLog("Init values pushed to workers")
    
    cmdStr = "nc = SnuddaDetect(configFile=configFile, positionFile=positionFile,voxelSize=voxelSize,hyperVoxelSize=hyperVoxelSize,verbose=verbose,logFileName=logFileName[0],saveFile=saveFile,SlurmID=SlurmID,role='worker')"

    #import pdb
    #pdb.set_trace()
    
    dView.execute(cmdStr,block=True)

    self.writeLog("Workers setup: " \
                  + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    self.workersInitialised = True
      
  ############################################################################

  def findMinMaxCoordParallell(self,volumeID=None,dView=None):

    if(dView is None):
      self.writeLog("findMinMaxCoordParallell: dView is None")
      return self.findMinMaxCoord(volumeID=volumeID)
    
    self.writeLog("Finding min/max coords parallel")
    
    neuronIdx = np.random.permutation(np.arange(0,len(self.neurons),
                                                dtype=np.int32))    

    dView.scatter("neuronIdxFind",neuronIdx,block=True)
    dView.push({"volumeID":volumeID},block=True)

    cmdStr = "minMax = nc.findMinMaxCoord(volumeID=volumeID," \
                                       + "neuronIdx=neuronIdxFind)"
    
    dView.execute(cmdStr,block=True)
    
    self.writeLog("Execution of min/max complete")
    #allMinMax = dView.gather("minMax",block=True)
    allMinMax = dView["minMax"]
    self.writeLog("Gathered min/max - complete.")
    
    maxCoord = -1e6*np.ones((3,))
    minCoord = 1e6*np.ones((3,))

    try:
      for (minC,maxC) in allMinMax:
        maxCoord = np.maximum(maxCoord,maxC)
        minCoord = np.minimum(minCoord,minC)
    except:
      import traceback
      tstr = traceback.format_exc()
      print(tstr)
      import pdb
      pdb.set_trace()
        
    return (minCoord,maxCoord)
   
    
  ############################################################################
  
  def findMinMaxCoord(self,volumeID=None,neuronIdx=None):

    if(volumeID is None):
      volumeID = self.volumeID

    print("Finding minMax coord in volumeID = " + str(volumeID))
    
    maxCoord = -1e6*np.ones((3,))
    minCoord = 1e6*np.ones((3,))
    
    if(neuronIdx is None):
      neurons = self.neurons
    else:
      neurons = [self.neurons[idx] for idx in neuronIdx]
    
    for n in neurons:
      
      # By using "in" for comparison, we can pass a list of volumeID also
      if volumeID is not None and n["volumeID"] not in volumeID:
        self.writeLog("Skipping " + n["name"] \
                      + " when calculating hyper voxel size")
        # Only include neurons belonging to the volume ID
        # we are looking at now
        #import pdb
        #pdb.set_trace()
        continue
      
      neuron = self.loadNeuron(n)
      
      if(len(neuron.dend) > 0):
        maxCoord = np.maximum(maxCoord,np.max(neuron.dend[:,:3],axis=0))
        minCoord = np.minimum(minCoord,np.min(neuron.dend[:,:3],axis=0))
        
      if(len(neuron.axon) > 0):
        maxCoord = np.maximum(maxCoord,np.max(neuron.axon[:,:3],axis=0))
        minCoord = np.minimum(minCoord,np.min(neuron.axon[:,:3],axis=0))

      maxCoord = np.maximum(maxCoord,np.max(neuron.soma[:,:3],axis=0))
      minCoord = np.minimum(minCoord,np.min(neuron.soma[:,:3],axis=0))
        

    return (minCoord,maxCoord)

  ############################################################################

  
  def fillVoxelsSoma(self,voxelSpace,voxelSpaceCtr,
                     voxelSecID, voxelSecX,
                     somaCoord,neuronID,verbose=False):

    vCoords = np.floor((somaCoord[0,:3] - self.hyperVoxelOrigo) \
                       / self.voxelSize).astype(int)
    radius2 = somaCoord[0,3] ** 2
    vRadius = np.ceil(somaCoord[0,3]/self.voxelSize).astype(int)

    assert vRadius < 1000, "fillVoxelsSoma: vRadius=" + str(vRadius) \
      + " soma coords = " + str(somaCoord) + " (BIG SOMA, not SI units?)"

    # Range check, so we stay within hypervoxel
    vxMin = max(0,vCoords[0]-vRadius)
    vxMax = min(self.hyperVoxelSize,vCoords[0]+vRadius+1)
    
    vyMin = max(0,vCoords[1]-vRadius)
    vyMax = min(self.hyperVoxelSize,vCoords[1]+vRadius+1)
    
    vzMin = max(0,vCoords[2]-vRadius)
    vzMax = min(self.hyperVoxelSize,vCoords[2]+vRadius+1)

    if(verbose):
      print("Soma check x: " + str(vxMin) + " - " + str(vxMax) \
            + " y: " + str(vyMin) + " - " + str(vyMax) \
            + " z: " + str(vzMin) + " - " + str(vzMax))

    for vx in range(vxMin,vxMax):
      for vy in range(vyMin,vyMax):
        for vz in range(vzMin,vzMax):

          d2 = ((vx+0.5)*self.voxelSize \
                 + self.hyperVoxelOrigo[0]-somaCoord[0,0])**2 \
                + ((vy+0.5)*self.voxelSize\
                   + self.hyperVoxelOrigo[1]-somaCoord[0,1])**2 \
                + ((vz+0.5)*self.voxelSize\
                   +self.hyperVoxelOrigo[2]-somaCoord[0,2])**2
          
          if(d2 < radius2):
            # Mark the point
            try:
              vCtr = voxelSpaceCtr[vx,vy,vz]

              if(vCtr > 0 \
                 and voxelSpace[vx,vy,vz,vCtr-1] == neuronID):
                # Voxel already has neuronID, skip
                continue
              
              voxelSpace[vx,vy,vz,vCtr] = neuronID
              voxelSecID[vx,vy,vz,vCtr] = 0 # Soma is 0
              voxelSecX[vx,vy,vz,vCtr] = 0.5

              voxelSpaceCtr[vx,vy,vz] += 1
            except:
              self.voxelOverflowCounter += 1
              self.writeLog("!!! If you see this you need to increase " \
                            + "maxDend above " \
                            + str(voxelSpaceCtr[vx,vy,vz]))
              continue


  
  ############################################################################ 
  
  # This uses self.hyperVoxelOrigo, self.voxelSize, self.nBins

  # !!! OBS segX must be an integer here, so to get true segX divide by 10000
  
  def fillVoxelsDend(self,voxelSpace,voxelSpaceCtr,
                     voxelSecID, voxelSecX,
                     voxelSomaDist,
                     coords,links,
                     segID,segX,neuronID):
    
    # segID gives segment ID for each link
    # segX gives segmentX for each link
    
    for line,segmentID,segmentX in zip(links,segID,segX):
      p1 = coords[line[0],:3]
      p2 = coords[line[1],:3]
      p1Dist = coords[line[0],4] * 1e6 # Dist to soma
      p2Dist = coords[line[1],4] * 1e6
      
      vp1 = np.floor((p1 - self.hyperVoxelOrigo) / self.voxelSize).astype(int)
      vp2 = np.floor((p2 - self.hyperVoxelOrigo) / self.voxelSize).astype(int)
      
      vp1Inside = ((vp1 >= 0).all() and (vp1 < self.nBins).all())
      vp2Inside = ((vp2 >= 0).all() and (vp2 < self.nBins).all())      

      # Four cases, if neither inside, skip line
      # If one inside but not the other, start at inside point and
      # continue until outside
      # If both inside, add all points without checking if they are inside
      
      if(not vp1Inside and not vp2Inside):
        # No points inside, skip
        continue

      if((vp1 == vp2).all()):
        # Line is only one voxel, steps will be 0, so treat it separately
        # We know it is inside, since they are same and both not outside

        vCtr = voxelSpaceCtr[vp1[0],vp1[1],vp1[2]]
        if(vCtr > 0 and voxelSpace[vp1[0],vp1[1],vp1[2],vCtr-1] == neuronID):
            # Voxel already has neuronID, skip
            continue
        
        try:
          voxelSpace[vp1[0],vp1[1],vp1[2],vCtr] = neuronID
          voxelSecID[vp1[0],vp1[1],vp1[2],vCtr] = segmentID
          voxelSecX[vp1[0],vp1[1],vp1[2],vCtr] = segmentX[0]
          voxelSomaDist[vp1[0],vp1[1],vp1[2],vCtr] = p1Dist
          
          voxelSpaceCtr[vp1[0],vp1[1],vp1[2]] += 1

        except:
          self.voxelOverflowCounter += 1
          self.writeLog("!!! If you see this you need to increase " \
                        + "maxDend above " \
                        + str(voxelSpaceCtr[vp1[0],vp1[1],vp1[2]]))
          continue
        
        # Done, next voxel
        continue
      
      if(not vp1Inside):
        if(not vp2Inside):
          # No point inside, skip
          continue
        else:          
          # Start with vp2 continue until outside cube
          steps = max(np.abs(vp2-vp1))
          dv = (vp1 - vp2)/steps
          ds = (segmentX[0] - segmentX[1])/steps
          dd = (p1Dist - p2Dist)/steps

          # We want the end element "steps" also, hence +1          
          for i in range(0,steps+1):
            vp = (vp2 + dv * i).astype(int)
            sX = segmentX[1] + ds*i # float
            somaDist = (p2Dist + dd*i).astype(int)
            
            if((vp < 0).any() or (vp >= self.nBins).any()):
              # Rest of line outside
              break

            try:
              vCtr = voxelSpaceCtr[vp[0],vp[1],vp[2]]
              if(vCtr > 0 and voxelSpace[vp[0],vp[1],vp[2],vCtr-1] == neuronID):
                  # Voxel already contains neuronID, skip
                  continue
              
              voxelSpace[vp[0],vp[1],vp[2],vCtr] = neuronID
              voxelSecID[vp[0],vp[1],vp[2],vCtr] = segmentID
              voxelSecX[vp[0],vp[1],vp[2],vCtr] = sX
              voxelSomaDist[vp[0],vp[1],vp[2],vCtr] = somaDist
              
              voxelSpaceCtr[vp[0],vp[1],vp[2]] += 1
            except:
              # Increase maxAxon and maxDend
              self.writeLog("!!! If you see this you need to increase " \
                            + "maxDend above " \
                            + str(voxelSpaceCtr[vp[0],vp[1],vp[2]]))
              self.voxelOverflowCounter += 1
              continue
            
      elif(not vp2Inside):
        # Start with vp1 continue until outside cube
        steps = max(np.abs(vp2-vp1))
        dv = (vp2 - vp1)/steps
        ds = (segmentX[1] - segmentX[0])/steps
        dd = (p2Dist - p1Dist)/steps        

        # We want the end element "steps" also, hence +1
        for i in range(0,steps+1):
          vp = (vp1 + dv * i).astype(int)
          sX = segmentX[0] + ds*i # float
          somaDist = (p1Dist + dd*i).astype(int)
          
          if((vp < 0).any() or (vp >= self.nBins).any()):
            # Rest of line outside
            break

          try:
            vCtr = voxelSpaceCtr[vp[0],vp[1],vp[2]]

            if(vCtr > 0 and voxelSpace[vp[0],vp[1],vp[2],vCtr-1] == neuronID):
                # Voxel already contains neuronID, skip
                continue
            
            voxelSpace[vp[0],vp[1],vp[2],vCtr] = neuronID
            voxelSecID[vp[0],vp[1],vp[2],vCtr] = segmentID
            voxelSecX[vp[0],vp[1],vp[2],vCtr] = sX
            voxelSomaDist[vp[0],vp[1],vp[2],vCtr] = somaDist
            
            voxelSpaceCtr[vp[0],vp[1],vp[2]] += 1
          except:
            self.writeLog("!!! If you see this you need to increase " \
                          + "maxDend above " \
                          + str(voxelSpaceCtr[vp[0],vp[1],vp[2]]))
            self.voxelOverflowCounter += 1
            continue
            
      else:
        # Entire line inside        
        steps = max(np.abs(vp2-vp1))
        dv = (vp2 - vp1)/steps
        ds = (segmentX[1] - segmentX[0])/steps
        dd = (p2Dist - p1Dist)/steps
      
        for i in range(0,steps+1):
          vp = (vp1 + dv * i).astype(int)
          sX = segmentX[0] + ds*i # float
          somaDist = (p1Dist + dd*i).astype(int)

          try:
            vCtr = voxelSpaceCtr[vp[0],vp[1],vp[2]]
            
            if(vCtr>0 and voxelSpace[vp[0],vp[1],vp[2],vCtr-1] == neuronID):
              # Voxel already has neuronID, skip
              continue
            
            voxelSpace[vp[0],vp[1],vp[2],vCtr] = neuronID
            voxelSecID[vp[0],vp[1],vp[2],vCtr] = segmentID
            voxelSecX[vp[0],vp[1],vp[2],vCtr] = sX
            voxelSomaDist[vp[0],vp[1],vp[2],vCtr] = somaDist
            
            voxelSpaceCtr[vp[0],vp[1],vp[2]] += 1
          except:
            self.writeLog("!!! If you see this you need to increase " \
                          + "maxDend above " \
                          + str(voxelSpaceCtr[vp[0],vp[1],vp[2]]))
            self.voxelOverflowCounter += 1
            continue
            
        
      # Potentially faster?    
      # http://code.activestate.com/recipes/578112-bresenhams-line-algorithm-in-n-dimensions/

  ############################################################################
  
  def fillVoxelsAxon(self,voxelSpace,voxelSpaceCtr,
                     voxelAxonDist,
                     coords,links,
                     neuronID):
    
    # segID gives segment ID for each link
    # segX gives segmentX for each link
    
    for line in links:
      p1 = coords[line[0],:3]
      p2 = coords[line[1],:3]
      p1Dist = coords[line[0],4] * 1e6 # Dist to soma
      p2Dist = coords[line[1],4] * 1e6
      
      vp1 = np.floor((p1 - self.hyperVoxelOrigo) / self.voxelSize).astype(int)
      vp2 = np.floor((p2 - self.hyperVoxelOrigo) / self.voxelSize).astype(int)
      
      vp1Inside = ((vp1 >= 0).all() and (vp1 < self.nBins).all())
      vp2Inside = ((vp2 >= 0).all() and (vp2 < self.nBins).all())      

      # Four cases, if neither inside, skip line
      # If one inside but not the other, start at inside point and
      # continue until outside
      # If both inside, add all points without checking if they are inside
      
      if(not vp1Inside and not vp2Inside):
        # No points inside, skip
        continue

      if((vp1 == vp2).all()):
        # Line is only one voxel, steps will be 0, so treat it separately
        # We know it is inside, since they are same and both not outside
        try:
          vCtr = voxelSpaceCtr[vp1[0],vp1[1],vp1[2]]
          if(vCtr > 0 and voxelSpace[vp1[0],vp1[1],vp1[2],vCtr-1] == neuronID):
            # Voxel already has neuronID, skip
            continue
          
          voxelSpace[vp1[0],vp1[1],vp1[2],vCtr] = neuronID
          voxelAxonDist[vp1[0],vp1[1],vp1[2],vCtr] = p1Dist
          
          voxelSpaceCtr[vp1[0],vp1[1],vp1[2]] += 1

        except Exception as e:

          import traceback
          tstr = traceback.format_exc()
          print(tstr)
          
          self.voxelOverflowCounter += 1
          self.writeLog("!!! If you see this you need to increase " \
                        + "maxAxon above " \
                        + str(voxelSpaceCtr[vp1[0],vp1[1],vp1[2]]))
          continue
        
        # Done, next voxel
        continue
      
      if(not vp1Inside):
        if(not vp2Inside):
          # No point inside, skip
          continue
        else:          
          # Start with vp2 continue until outside cube
          steps = max(np.abs(vp2-vp1))
          dv = (vp1 - vp2)/steps
          dd = (p1Dist - p2Dist)/steps

          # We want the end element "steps" also, hence +1          
          for i in range(0,steps+1):
            vp = (vp2 + dv * i).astype(int)
            axDist = (p2Dist + dd*i).astype(int)
        
            if((vp < 0).any() or (vp >= self.nBins).any()):
              # Rest of line outside
              break

            try:
              vCtr = voxelSpaceCtr[vp[0],vp[1],vp[2]]
              if(vCtr>0 and voxelSpace[vp[0],vp[1],vp[2],vCtr-1]==neuronID):
                # Voxel already has neuronID, skip
                continue
              
              voxelSpace[vp[0],vp[1],vp[2],vCtr] = neuronID
              voxelAxonDist[vp[0],vp[1],vp[2],vCtr] = axDist

              voxelSpaceCtr[vp[0],vp[1],vp[2]] += 1
            except Exception as e:
              import traceback
              tstr = traceback.format_exc()
              print(tstr)
              
              # Increase maxAxon and maxDend
              self.writeLog("!!! If you see this you need to increase " \
                            + "maxAxon above " \
                            + str(voxelSpaceCtr[vp[0],vp[1],vp[2]]))
              self.voxelOverflowCounter += 1
              continue
            
      elif(not vp2Inside):
        # Start with vp1 continue until outside cube
        steps = max(np.abs(vp2-vp1))
        dv = (vp2 - vp1)/steps
        dd = (p2Dist - p1Dist)/steps

        # We want the end element "steps" also, hence +1
        for i in range(0,steps+1):
          vp = (vp1 + dv * i).astype(int)
          axDist = (p1Dist + dd*i).astype(int)
        
          if((vp < 0).any() or (vp >= self.nBins).any()):
            # Rest of line outside
            break

          try:
            vCtr = voxelSpaceCtr[vp[0],vp[1],vp[2]]
            if(vCtr>0 and voxelSpace[vp[0],vp[1],vp[2],vCtr-1] == neuronID):
              # Voxel already has neuronID, skip
              continue
            
            voxelSpace[vp[0],vp[1],vp[2],vCtr] = neuronID
            voxelAxonDist[vp[0],vp[1],vp[2],vCtr] = axDist
            
            voxelSpaceCtr[vp[0],vp[1],vp[2]] += 1
          except Exception as e:

            import traceback
            tstr = traceback.format_exc()
            print(tstr)

            self.writeLog("!!! If you see this you need to increase " \
                          + "maxAxon above " \
                          + str(voxelSpaceCtr[vp[0],vp[1],vp[2]]))
            self.voxelOverflowCounter += 1
            continue
            
      else:
        # Entire line inside        
        steps = max(np.abs(vp2-vp1))
        dv = (vp2 - vp1)/steps
        dd = (p2Dist - p1Dist)/steps
      
        for i in range(0,steps+1):
          vp = (vp1 + dv * i).astype(int)
          axDist = (p1Dist + dd*i).astype(int)

          try:
            vCtr = voxelSpaceCtr[vp[0],vp[1],vp[2]]
            if(vCtr>0 and voxelSpace[vp[0],vp[1],vp[2],vCtr-1] == neuronID):
              # Voxel already has neuronID, skip
              continue
            
            voxelSpace[vp[0],vp[1],vp[2],vCtr] = neuronID
            voxelAxonDist[vp[0],vp[1],vp[2],vCtr] = axDist
            
            voxelSpaceCtr[vp[0],vp[1],vp[2]] += 1
          except Exception as e:
            import traceback
            tstr = traceback.format_exc()
            print(tstr)
            
            self.writeLog("!!! If you see this you need to increase " \
                          + "maxAxon above " \
                          + str(voxelSpaceCtr[vp[0],vp[1],vp[2]]))
            self.voxelOverflowCounter += 1
            continue
            
        
      # Potentially faster?    
      # http://code.activestate.com/recipes/578112-bresenhams-line-algorithm-in-n-dimensions/


  ############################################################################

  def getPath(self,pathStr):

    return pathStr.replace("$DATA", os.path.dirname(__file__) + "/data")
          
  ############################################################################

  
  def processHyperVoxel(self,hyperID):

    startTime = timeit.default_timer()
    
    if(self.hyperVoxels[hyperID]["neuronCtr"] == 0):
      # No neurons, return quickly - do not write hdf5 file
      endTime = timeit.default_timer()
      return (hyperID,0,0,endTime-startTime,0)
    
    hOrigo = self.hyperVoxels[hyperID]["origo"]
    self.setupHyperVoxel(hOrigo,hyperID)
    
    nNeurons = self.hyperVoxels[hyperID]["neuronCtr"]

    self.writeLog("Processing hyper voxel : " + str(hyperID) \
                  + "/" + str(self.hyperVoxelIDs.size) \
                  + " (" + str(nNeurons) + " neurons)")

    # !!! Suggestion for optimisation. Place neurons with GJ first, then do
    # GJ touch detection, after that add rest of neurons (to get comlete set)
    # and then do axon-dend synapse touch detection
    
    for neuronID in self.hyperVoxels[hyperID]["neurons"][:nNeurons]:

      neuron = self.loadNeuron(self.neurons[neuronID])
      
      try:
        self.fillVoxelsAxon(self.axonVoxels,
                            self.axonVoxelCtr,
                            self.axonSomaDist,
                            neuron.axon,
                            neuron.axonLinks,
                            neuronID)

        self.fillVoxelsSoma(self.dendVoxels,
                            self.dendVoxelCtr,
                            self.dendSecID,
                            self.dendSecX,
                            neuron.soma,
                            neuronID)

        self.fillVoxelsDend(self.dendVoxels,
                            self.dendVoxelCtr,
                            self.dendSecID,
                            self.dendSecX,
                            self.dendSomaDist,
                            neuron.dend,
                            neuron.dendLinks,
                            neuron.dendSecID,
                            neuron.dendSecX, 
                            neuronID)
        
      except Exception as e:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)

        self.writeLog("Something went wrong: " +str(tstr))
        import pdb
        pdb.set_trace()

    # This should be outside the neuron loop
    try:
      # This places axon voxels for neurons without axon morphologies
      self.placeSynapsesNoAxon(hyperID,
                               self.axonVoxels,
                               self.axonVoxelCtr,
                               self.axonSomaDist)
    except Exception as e:
      import traceback
      tstr = traceback.format_exc()
      print(tstr)

      self.writeLog("Something wrong: " +str(tstr))
      import pdb
      pdb.set_trace()

        
    # This detects the synapses where we use a density distribution for axons
    # self.detectSynapsesNoAxonSLOW (hyperID) # --replaced by placeSynapseNoAxon

    # The normal voxel synapse detection
    self.detectSynapses()

    self.detectGapJunctions()

    
    self.writeHyperVoxelToHDF5()

    endTime = timeit.default_timer()

    #if(nNeurons < 3):
    #  print("DEBUGGING...")
    #  self.plotHyperVoxel(hyperID=hyperID,plotNeurons=True)
    #  import pdb
    #  pdb.set_trace()

    return (hyperID,self.hyperVoxelSynapseCtr,
            self.hyperVoxelGapJunctionCtr,endTime-startTime,
            self.voxelOverflowCounter)
    
  ############################################################################

  # hyperID is just needed if we want to plotNeurons also
  
  def plotHyperVoxel(self,plotNeurons=False,drawAxons=True,drawDendrites=True,
                     detectDone=True):
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  

    colors = np.zeros((self.dendVoxelCtr.shape[0], 
                       self.dendVoxelCtr.shape[1], 
                       self.dendVoxelCtr.shape[2],4))
    colors[:,:,:,3] = 0.3

    voxelData = np.zeros((self.dendVoxelCtr.shape[0], 
                          self.dendVoxelCtr.shape[1], 
                          self.dendVoxelCtr.shape[2]))

    if(drawAxons):
      colors[:,:,:,0] = self.axonVoxelCtr / max(np.max(self.axonVoxelCtr),1)
      voxelData += self.axonVoxelCtr
      
    if(drawDendrites):
      colors[:,:,:,2] = self.dendVoxelCtr / max(np.max(self.dendVoxelCtr),1)
      voxelData += self.dendVoxelCtr
      
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.voxels(voxelData > 0,
              facecolors=colors, edgecolor=None)

    if(self.hyperVoxelSynapseCtr > 0):
      sCoord = self.hyperVoxelSynapses[:self.hyperVoxelSynapseCtr,2:5]

      # In case hyperVoxelOffset has been applied, we need to subtract it
      # to draw within the hyper voxel
      if(self.hyperVoxelOffset is not None):
        sCoord - self.hyperVoxelOffset
      ax.scatter(sCoord[:,0],sCoord[:,1],sCoord[:,2],c="green")
    
    plt.ion()
    plt.show()
    
    if(plotNeurons):
      # Also plot the neurons overlayed, to verify

      nNeurons = self.hyperVoxels[self.hyperVoxelID]["neuronCtr"]

      for neuronID in self.hyperVoxels[self.hyperVoxelID]["neurons"][:nNeurons]:
        
        neuron = self.loadNeuron(self.neurons[neuronID])

        neuron.plotNeuron(axis=ax,
                          plotAxon=drawAxons,
                          plotDendrite=drawDendrites,
                          plotOrigo=self.hyperVoxelOrigo,plotScale=1/self.voxelSize)
      

    plt.show()
    plt.ion()

    plt.pause(0.001)
    figName = "figures/Hypervoxel-" + str(self.SlurmID) \
              + "-" + str(self.hyperVoxelID) + ".png"
    plt.savefig(figName)

  ############################################################################

  def exportVoxelVisualisationCSV(self, neuronID):

    # x,y,z = coords
    # shape = "cube" or "sphere"
    # type = "axon", "dendrite", "synapse"
    # id = neuronID
    # x,y,z,shape,type,id

    headerStr = "# x,y,z,shape,type,id\n"
    axonStr = ""
    dendStr = ""
    synapseStr = ""

    for x in range(0,self.axonVoxelCtr.shape[0]):
      for y in range(0,self.axonVoxelCtr.shape[1]):
        for z in range(0,self.axonVoxelCtr.shape[2]):
          for c in range(0,self.axonVoxelCtr[x,y,z]):
            nID = self.axonVoxels[x,y,z,c]
            if(nID in neuronID):
              axonStr += str(x) + "," + str(y) + "," + str(z) \
                + ",cube,axon," + str(nID) + "\n"

    for x in range(0,self.dendVoxelCtr.shape[0]):
      for y in range(0,self.dendVoxelCtr.shape[1]):
        for z in range(0,self.dendVoxelCtr.shape[2]):
          for c in range(0,self.dendVoxelCtr[x,y,z]):
            nID = self.dendVoxels[x,y,z,c]
            if(nID in neuronID):
              dendStr += str(x) + "," + str(y) + "," + str(z) \
                + ",cube,dend," + str(nID) + "\n"

    synList = []
    for ir, row in enumerate(self.hyperVoxelSynapses):
      if(row[0] in neuronID and row[1] in neuronID):
        synList.append(ir)

    for i in synList:
      xyz = self.hyperVoxelSynapses[i,2:5]
      synapseStr += str(xyz[0]) + "," + str(xyz[1]) + "," + str(xyz[2]) \
        + ",sphere,synapse," + str(self.hyperVoxelSynapses[i,1]) +"\n"

    fName = "hypervoxel-" + str(self.SlurmID) + ".csv"
    with open(fName,'w') as f:
      f.write(headerStr)
      f.write(axonStr)
      f.write(dendStr)
      f.write(synapseStr)
    
  ############################################################################

  # Example usage:
  # nc.plotNeuronsInHyperVoxel(neuronID=[1,20],
  #                            neuronColour=np.array([[0,0,1],[0,1,0]]),
  #                            axonAlpha=[1,0.3],dendAlpha=[0.3,1])
 
  
  # each row in neuronColour is a colour for a neuron
  
  def plotNeuronsInHyperVoxel(self,neuronID,neuronColour,
                              axonAlpha=None,dendAlpha=None):

    if(axonAlpha is None):
      axonAlpha = ones((len(neuronID),))

    if(dendAlpha is None):
      dendAlpha = ones((len(neuronID),))

    alphaAxonLookup = dict([])
    alphaDendLookup = dict([])
    neuronColourLookup = dict([])

    for ni,aa,da,nc in zip(neuronID,axonAlpha,dendAlpha,neuronColour):
      alphaAxonLookup[ni] = aa
      alphaDendLookup[ni] = da
      neuronColourLookup[ni] = nc
                      
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  

    colors = np.zeros((self.dendVoxelCtr.shape[0], 
                       self.dendVoxelCtr.shape[1], 
                       self.dendVoxelCtr.shape[2],4))

    voxelData = np.zeros((self.dendVoxelCtr.shape[0], 
                          self.dendVoxelCtr.shape[1], 
                          self.dendVoxelCtr.shape[2]))

    
    for ix in range(0,self.axonVoxelCtr.shape[0]):
      for iy in range(0,self.axonVoxelCtr.shape[1]):
        for iz in range(0,self.axonVoxelCtr.shape[2]):
          for ic in range(0,self.axonVoxelCtr[ix,iy,iz]):
            nID = self.axonVoxels[ix,iy,iz,ic]
            if(nID in neuronID):
              colors[ix,iy,iz,0:3] = neuronColourLookup[nID]
              colors[ix,iy,iz,3] = alphaAxonLookup[nID]
              voxelData[ix,iy,iz] = 1
    
    for ix in range(0,self.dendVoxelCtr.shape[0]):
      for iy in range(0,self.dendVoxelCtr.shape[1]):
        for iz in range(0,self.dendVoxelCtr.shape[2]):
          for ic in range(0,self.dendVoxelCtr[ix,iy,iz]):
            nID = self.dendVoxels[ix,iy,iz,ic]
            if(nID in neuronID):
              colors[ix,iy,iz,0:3] = neuronColourLookup[nID]
              colors[ix,iy,iz,3] = alphaDendLookup[nID]
              voxelData[ix,iy,iz] = 1
              
              
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.voxels(voxelData > 0,
              facecolors=colors, edgecolor=None)

    printList = []
    for ir,row in enumerate(self.hyperVoxelSynapses):
      if(row[0] in neuronID and row[1] in neuronID):
        printList.append(ir)
        
    
     # This should really only plot those between the neurons indicated
    if(len(printList) > 0):
      sCoord = self.hyperVoxelSynapses[printList,2:5]
      ax.scatter(sCoord[:,0]+0.5,sCoord[:,1]+0.5,sCoord[:,2]+0.5,c="red",s=100)

    plt.axis("off")
    plt.ion()
    plt.show()

    #import pdb
    #pdb.set_trace()
    
    plt.pause(0.001)
    figName = "figures/Hypervoxel-" + str(self.SlurmID) \
              + "-" + str(self.hyperVoxelID) + "-someNeurons.png"
    plt.savefig(figName,dpi=900)

      
    
  ############################################################################

  def trivialExample(self):

    # This places two neurons in a hyper voxel and tests it
    # --- It destroys the current state of the program.

    self.hyperVoxels[0]["neuronCtr"] = 2
    self.hyperVoxels[0]["origo"] = np.array([0,0,0])
    self.neurons[50]["position"] = np.array([-10e-6,120e-6,-10e-6])
    self.neurons[51]["position"] = np.array([120e-6,120e-6,-10e-6])
    self.hyperVoxels[0]["neurons"] = np.array([50,51])

    self.processHyperVoxel(0)

    #self.plotHyperVoxel(plotNeurons=True)
    self.plotHyperVoxel(plotNeurons=False)
    
    print("Synapses: " + str(self.hyperVoxelSynapseCtr))
    print(str(self.hyperVoxelSynapses[:self.hyperVoxelSynapseCtr,:]))

    import pdb
    pdb.set_trace()

    
    
  
  ############################################################################

  def testVoxelDraw(self):

    print("This changes internal state of the object, restart after test run.")

    self.hyperVoxelID = -1
    self.hyperVoxelOrigo = np.zeros((3,))
    self.voxelSize = 2
    self.nBins = np.ones((3,1))*10
    
    voxels = np.zeros((10,10,10,10),dtype=int)
    voxelCtr = np.zeros((10,10,10),dtype=int)
    voxelSecID = np.zeros((10,10,10,10),dtype=int)
    voxelSecX = np.zeros((10,10,10,10),dtype=float)
    voxelSomaDist = np.zeros((10,10,10,10),dtype=int)

    coords = np.array([[2,2,2,1.1,40],[8,10,8,1.2,50],[0,23,22,1.3,60]])
    links = np.array([[0,1],[0,2],[2,1]],dtype=int)
    segID = np.array([1,2,3],dtype=int)
    segX = np.array([[0.1,0.2],[0.3,0.4],[0.5,0.6]],dtype=float)

    if(False):
      self.fillVoxelsDend(voxelSpace=voxels,
                          voxelSpaceCtr=voxelCtr,
                          voxelSecID=voxelSecID,
                          voxelSecX=voxelSecX,
                          voxelSomaDist=voxelSomaDist,
                          coords=coords,
                          links=links,
                          segID=segID,
                          segX=segX,
                          neuronID=13)

    if(False):
      self.fillVoxelsSoma(voxelSpace=voxels,
                          voxelSpaceCtr=voxelCtr,
                          voxelSecID=voxelSecID,
                          voxelSecX=voxelSecX,
                          somaCoord=np.array([[10,10,10,8]]),
                          neuronID=14)
      

    #import pdb
    #pdb.set_trace()

    # We also need to check axon filling

    voxels[:] = 0
    voxelCtr[:] = 0
    voxelSomaDist[:] = 0
    
    self.fillVoxelsAxon(voxelSpace=voxels,
                        voxelSpaceCtr=voxelCtr,
                        voxelAxonDist=voxelSomaDist,
                        coords=coords,
                        links=links,
                        neuronID=13)

    import pdb
    pdb.set_trace()
    
############################################################################

  # Memory check code taken from
  # https://stackoverflow.com/questions/17718449/determine-free-ram-in-python/17718729#17718729
  #
  
  def memory(self):
      """
      Get node total memory and memory usage
      """
      try:
        with open('/proc/meminfo', 'r') as mem:
          ret = {}
          tmp = 0
          for i in mem:
            sline = i.split()
            if str(sline[0]) == 'MemTotal:':
              ret['total'] = int(sline[1])
            elif str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
              tmp += int(sline[1])
          ret['free'] = tmp
          ret['used'] = int(ret['total']) - int(ret['free'])
        return ret
      except:
        return "Non-linux system, /proc/meminfo unavailable"


    
############################################################################

  @staticmethod
  def _processHyperVoxelHelper(hyperID):

    mem = nc.memory()
    nc.writeLog("Memory status, before processing " + str(hyperID) \
                + ": "+ str(mem))
   
    return nc.processHyperVoxel(hyperID)
    
############################################################################

def nextRunID():

  runIDfile = ".runID.pickle"

  try:
    if(os.path.isfile(runIDfile)):
      
      with open(runIDfile,'rb') as f:
        runID = pickle.load(f)
        nextID = int(np.ceil(np.max(runID)) + 1)

      runID.append(nextID)
        
      with open(runIDfile,'wb') as f:
        pickle.dump(runID,f,-1)

    else:
      
      with open(runIDfile,'wb') as f:
        nextID = 1
        runID = [1]
        pickle.dump(runID,f,-1)
        
  except Exception as e:
    import traceback
    tstr = traceback.format_exc()
    print(tstr)

    print("Problem reading .runID.pickle file, setting runID to 0")
    import pdb
    pdb.set_trace()
    return 0

  print("Using runID = " + str(nextID))
  
  return nextID



############################################################################

if __name__ == "__main__":

  print("Please do not call this file directly, use snudda.py")
  exit(-1)

  
