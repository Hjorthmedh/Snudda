# snudda_prune.py
#
# Johannes Hjorth, Royal Institute of Technology (KTH)
# Human Brain Project 2019
#
# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union’s Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).
#
#
# Currently a n-way merge sort is done in serial using a heap
#
# New idea for how to speed up the merge:
#
# Divide neurons into groups depending on spatial location.
# We have a list with which hyper voxel each neuron is in.
# Each worker opens their files, and does an n-way merge
# sort ignoring all other neurons but their own
# Finally, the resulting files from the workers are merged together

import numpy as np
import scipy
# from scipy import sparse
import math
import os
import sys
import itertools
from glob import glob
import collections

import heapq # Priority queue
import queue

import time
import timeit

import h5py
import json
import pickle

from .Neuron_morphology import NeuronMorphology

# The x,y,z coordinates of the synapses in the original file are integer
# indexes refering to within the hyper voxel. Here they are transformed
# to refer to within the entire simulation, with the simulation origo as
# the origo (bin size unchanged)

# This code combines multiple voxel files
# !!! This also performs the pruning, since it requires knowledge of all
#     synapses between the neurons, something not possible within a hyper voxel
#     if the neuron crosses borders

class SnuddaPrune(object):

  ############################################################################
  
  def __init__(self,workHistoryFile,
               logFile=None,logFileName=None,
               dView=None,lbView=None,role="master",verbose=True,
               scratchPath=None,
               preMergeOnly=False,
               h5libver="latest",
               randomSeed=None,
               cleanVoxelFiles=True):

    # Help with parallel debugging, when workers cant print to screen:
    #self.writeToRandomFile("WH = " + str(workHistoryFile) \
    #                       + "\nlog = " + str(logFileName) \
    #                       + "\nrole = " + role)
    
    startTime = timeit.default_timer()
    self.workHistoryFile = workHistoryFile
    self.basePath = os.path.dirname(self.workHistoryFile)+"/"
    self.basePath = self.basePath.replace("/log/","/")
    
    self.logFile = logFile
    self.logFileName = logFileName
    self.verbose = verbose
    self.h5libver = h5libver
    
    if(logFile is None and logFileName is not None):
      self.logFile = open(logFileName,'w')
      self.writeLog("Log file created.")

    np.random.seed(randomSeed)
    self.writeLog("Setting random seed: " + str(randomSeed))
      
    self.h5driver = "sec2" #"core" # "stdio", "sec2"

    self.writeLog("Using hdf5 driver " + str(self.h5driver) \
                  + ", " + str(self.h5libver) + " version")
    
    self.dView = dView
    self.lbView = lbView
    self.role = role
    self.workersInitialised = False

    self.synapseTypeLookup = { 1 : "GABA",
                               2 : "AMPA_NMDA",
                               3 : "GapJunction",
                               4 : "ACh",
                               5 : "NO"}

    self.synapseTypeReverseLookup = \
        {v: k for k, v in self.synapseTypeLookup.items()}

    # Parameters for the HDF5 writing, this affects write speed
    self.synapseChunkSize = 10000
    self.gapJunctionChunkSize = 10000
    self.h5compression = "lzf"
    
    # These are for the merge code
    self.synapseWriteBuffer = None
    self.synapseBufferSize = 100000
    self.nextBufferWritePos = 0 
    self.nextFileWritePos = 0
    self.bufferOutFile = None
    self.fileList = []
    self.fileBuffers = []

    # These are for the pruning code
    self.mergeDataType = None # "synapses" or "gapJunctions"
    self.synapseReadBuffer = None
    self.nextSynapseBufferReadPos = 0
    self.nextSynapseFileReadPos = 0
    self.endOfSynapseReadRange = 0 
    
    self.histFile = None
    self.outFile = None
    self.tempFileList = [] # List of all temp files created
    
    self.nextMergeFileID = 0

    self.voxelOverflowCounter = 0
    self.overflowFiles = []
    
    self.openWorkHistoryFile(workHistoryFile=workHistoryFile)

    self.setScratchPath(scratchPath)
    self.loadPruningInformation()

    # (locationOfMatrix,locationOfN,locationOfCoords)
    self.dataLoc = {"synapses" : ("network/synapses",
                                  "nSynapses",
                                  "network/synapseLookup"), #range(2,5)),
                    "gapJunctions" : ("network/gapJunctions",
                                      "nGapJunctions",
                                      "network/gapJunctionLookup") }
                                      #range(6,9))}

    if(self.role == "master"):

      # This bit is only done by the master, which then delegates job to
      # the worker nodes

      # For the larger networks the merge can take hours, so we
      # check if the file already exists
      (mergeDone,synapses,synapseFile) = self.mergeFileExists()

      if(not mergeDone):

        if(self.dView):
          self.writeLog("Running parallel merge")

          synapseFile = self.bigMergeParallel()

          #print("We do not currently get all the synapses, where are the missing ones going? About 1.4% of synapses are missed")
          #import pdb
          #pdb.set_trace()

          
          assert not (synapseFile["network/synapses"][-1,:] == 0).all(), \
            "We are missing some synapses in the merge!"
            
          
        else:
          self.writeLog("Running merge in serial")
        
          (synapses,synapseFile) = \
            self.bigMergeLookup(mergeDataType="synapses",
                                cleanVoxelFiles=True)

          (gapJunctions,gapJunctionFile) = \
            self.bigMergeLookup(mergeDataType="gapJunctions",
                                cleanVoxelFiles=True)


      # When running on a cluster, we might want to do the serial parts
      # in a separte run, hence this option that allows us to prepare
      # the data for parallel execution.
      if(preMergeOnly):
        self.writeLog("Pre-merge of synapses done. preMergeOnly = " \
                      + str(preMergeOnly) + ", exiting.")

        self.histFile.close()
        self.histFile = None

        endTime = timeit.default_timer()

        self.writeLog("Total duration: " + str(endTime-startTime))
        return

      # We need to close and reopen file with other driver, to allow swmr mode
      fName = synapseFile.filename
      synapseFile.flush()
      synapseFile.close()
      synapseFile = h5py.File(fName,"r",libver=self.h5libver,swmr=True)

      
      # Delegate pruning
      if(self.dView is not None):
        self.pruneSynapsesParallel(synapseFile=synapseFile,
                                   mergeDataType="synapses",
                                   closeInputFile=False,
                                   setupOutFile=True)
        
        self.pruneSynapsesParallel(synapseFile=synapseFile,
                                   mergeDataType="gapJunctions",
                                   setupOutFile=False)
      else:
        # OBS, dont close the input file after first call
        self.pruneSynapses(synapseFile=synapseFile,
                           outputFileName=None,rowRange=None,
                           closeOutFile=False,
                           closeInputFile=False,
                           mergeDataType="synapses") 
        self.pruneSynapses(synapseFile=synapseFile,
                           outputFileName=None,rowRange=None,
                           closeOutFile=False,
                           closeInputFile=True,
                           mergeDataType="gapJunctions") 

      if(self.outFile is None):
        self.writeLog("No output file created, no synapses exist?")
        self.writeLog("Creating symbolic link to MERGE file instead")

        fDest=(os.path.dirname(self.workHistoryFile)+"/").replace("/log/","/")\
                                + "/network-pruned-synapses.hdf5"
        fSrc = os.path.basename(fName)
        
        print(str(fDest) + "->" + str(fSrc))
        os.symlink(fSrc,fDest)
        
        return
        
      nSynBefore = np.sum(self.histFile["nSynapses"][()])
      nSynAfter = self.outFile["network/nSynapses"][0]
        
      nOverflow = np.sum(self.histFile["voxelOverflowCounter"][()])
      nGJBefore = np.sum(self.histFile["nGapJunctions"][()])
      nGJAfter = self.outFile["network/nGapJunctions"][0]

      self.writeLog("Voxel overflows: " + str(nOverflow) \
                    + " (should be zero)")

      if(nSynBefore > 0):
        self.writeLog("Synapses before pruning: " + str(nSynBefore))
        self.writeLog("Synapses after pruning: " + str(nSynAfter) \
                      + " (" +str(round(100.0*nSynAfter/nSynBefore,2)) \
                      + " % kept)")
      else:
        self.writeLog("No synapses to prune")

      if(nGJBefore > 0):
        self.writeLog("Gap junctions before pruning " + str(nGJBefore))
        self.writeLog("Gap junctions after pruning " + str(nGJAfter) \
                      + " (" + str(round(100.0*nGJAfter/nGJBefore,2)) \
                      + " % kept)")
      else:
        self.writeLog("No gap junctions to prune.")
        
      self.cleanUpMergeFiles()
      
      self.histFile.close()
      self.histFile = None

      endTime = timeit.default_timer()

      self.writeLog("Total duration: " + str(endTime-startTime))

      if(self.voxelOverflowCounter > 0):
        print("Voxel overflows: " + str(self.voxelOverflowCounter) \
              + "\nIn files: ")
        for f in self.overflowFiles:
          print("Overflow in " + str(f))

  ############################################################################

  def setScratchPath(self,scratchPath=None):

    assert self.workHistoryFile is not None \
      and self.workHistoryFile is not "last", \
      "Need to call openWorkHistoryFile before setScratchPath"
    
    if(scratchPath is None):
      self.scratchPath = self.basePath + "/temp/"

      if not os.path.exists(self.scratchPath):
        os.makedirs(self.scratchPath)
      
      self.writeLog("Using default scratch path: " + self.scratchPath)
    else:
      self.scratchPath = scratchPath
      self.writeLog("User selected scratch path: " + self.scratchPath)
          
  ############################################################################

  def writeToRandomFile(self,text):
    
    import uuid
    tmp = open("save/tmp-log-file-" + str(uuid.uuid4()),'w')
    tmp.write(text)
    tmp.close()
    print(text)
    
  ############################################################################

  def __del__(self):

    try:
      if(self.histFile is not None):
        self.histFile.close()
    except:
      print("Hist file already closed?")

    try:
      if(self.outFile is not None):
        self.outFile.close()
        self.outFile = None
    except:
      print("Out file already closed?")

    self.cleanUpMergeFiles()
    
  ############################################################################

  def openWorkHistoryFile(self,workHistoryFile=None):

    if(workHistoryFile is None):
      workHistoryFile = self.workHistoryFile
    
    if(workHistoryFile == "last"):
      workHistoryFile = self.findLatestFile()

    self.writeLog("Opening work history file: " + workHistoryFile)

    self.histFile = h5py.File(workHistoryFile,'r')      
    self.workHistoryFile = workHistoryFile

    self.SlurmID = self.histFile["meta/SlurmID"][()]    
    self.hyperVoxelIDs = self.histFile["meta/hyperVoxelIDs"][()]
    self.allHyperIDs = self.histFile["allHyperIDs"][()]
    self.voxelSize = self.histFile["meta/voxelSize"][()]
    self.hyperVoxelSize = self.histFile["meta/hyperVoxelSize"][()] # num bins
    self.simulationOrigo = self.histFile["meta/simulationOrigo"][()]
    self.hyperVoxelWidth = self.voxelSize * self.hyperVoxelSize

    # Network_simulate.py uses axonStumpIDFlag = True
    # Neurodamus uses axonStumpIDFlag = False
    self.axonStumpIDFlag = self.histFile["meta/axonStumpIDFlag"]

    # We need a lookup table for offsets of hypervoxel origos
    self.hyperVoxelOffset = np.zeros((self.hyperVoxelIDs.size,3),dtype=int)
    for ix in range(0,self.hyperVoxelIDs.shape[0]):
      for iy in range(0,self.hyperVoxelIDs.shape[1]):
        for iz in range(0,self.hyperVoxelIDs.shape[2]):
          self.hyperVoxelOffset[self.hyperVoxelIDs[ix,iy,iz],:] \
            = np.array([ix,iy,iz]) * self.hyperVoxelSize
          
    # OBS the synapse and gap junction numbers are listed not in order of
    # hyperID but in order they were completed, to find out hyperID for
    # a particular one check self.histFile["completed"]
    self.nSynapsesTotal = np.sum(self.histFile["nSynapses"][()])
    self.nGapJunctionsTotal = np.sum(self.histFile["nGapJunctions"][()])

    self.configFile = self.histFile["meta/configFile"][()]
    self.positionFile = self.histFile["meta/positionFile"][()]


    self.detectConfig = json.loads(self.histFile["meta/config"][()])
    with open(self.histFile["meta/configFile"][()],"r") as f:
      self.config = json.load(f)    
    
    
  ############################################################################

  def checkHyperVoxelIntegrity(self,hFile,hFileName,verbose=False):

    if(verbose):
      self.writeLog("Checking that " + hFileName + " matches circuit settings")
    
    checkList = ["voxelSize", "hyperVoxelSize", "simulationOrigo",
                 "configFile", "positionFile", "SlurmID", "axonStumpIDFlag"]

    # Just some sanity checks
    for c in checkList:
      test = self.histFile["meta/"+c][()] == hFile["meta/"+c][()]
      if(type(test) == bool):
        assert test, \
          "Mismatch of " + c + " in file " + hFileName
      else:
        assert test.all(), \
          "Mismatch of " + c + " in file " + hFileName              
        
    # Get xyz coordinates of hyper voxel
    xyz = np.where(self.hyperVoxelIDs == hFile["meta/hyperVoxelID"][()])
    xyz = np.array([x[0] for x in xyz])

    # Just do a sanity check that the hypervoxel origo matches stored value
    hOrigo = self.simulationOrigo + self.hyperVoxelWidth * xyz
    assert (hOrigo == hFile["meta/hyperVoxelOrigo"][()]).all(), \
      "Hyper voxel origo mismatch in file " + hFileName

    ofc = hFile["meta/voxelOverflowCounter"][()] 

    if(ofc > 0):
      self.voxelOverflowCounter += ofc
      self.overflowFiles.append(hFileName)
      self.writeLog("Overflow of " + str(ofc) + " in " + hFileName)
    
        
  ############################################################################

  def openHyperVoxel(self,hyperVoxelID,verbose=False,verify=True):

    if(verbose):
      self.writeLog("Reading hypervoxel " + str(hyperVoxelID))
      
    hFileName = self.basePath + "/voxels/network-putative-synapse-" \
                + str(hyperVoxelID) + ".hdf5"

    
    hFile = h5py.File(hFileName)

    # Just make sure the data we open is OK and match the other data
    if(verify):
      self.checkHyperVoxelIntegrity(hFile,hFileName,verbose=verbose)

    return hFile
  
  ############################################################################

  # !!! Cache the results after creation, in case we want to rerun pruning 
  
  def createConnectionMatrix(self):

    if(self.histFile is None):
      self.openWorkHistoryFile()

    # Do we have the matrix cached?
    cMat = self.loadConnectionMatrixCache()
    if(cMat is not None):
      return cMat
      
    nNeurons = self.histFile["network/neurons/neuronID"].shape[0]

    self.writeLog("Creating the " +  str(nNeurons) +"x" + str(nNeurons) \
          + " sparse connection matrix") 
    
    # After creation, convert to csr or csc matrix
    conMat = scipy.sparse.lil_matrix((nNeurons,nNeurons),dtype=np.int16)


    self.writeLog("Parsing " + str(len(self.allHyperIDs)) + " hyper voxel files")
    
    for ctr, hID in enumerate(self.allHyperIDs):
      if(ctr % 1000 == 0 and ctr > 0):
        self.writeLog(str(ctr) + " files done")

      hFile = self.openHyperVoxel(hID)

      synapses = hFile["network/synapses"][()]

      for row in synapses:
        SrcID = row[0]
        DestID = row[1]

        conMat[SrcID,DestID] += 1

      hFile.close()

    self.writeLog("Converting to CSR sparse matrix format")
    cMat = conMat.tocsr()

    self.saveConnectionMatrixCache(cMat)
    
    return cMat

  ############################################################################

  def getConMatCacheFileName(self):

    cacheFile = self.basePath + "/connection-matrix-cache.pickle"

    return cacheFile
  
  ############################################################################
  
  def loadConnectionMatrixCache(self):

    cacheFile = self.getConMatCacheFileName()
    
    if(os.path.isfile(cacheFile)):
      self.writeLog("Loading connection matrix cache from " + cacheFile)
      
      with open(cacheFile,'rb') as f:
        data = pickle.load(f)
        conMat = data["conMat"]

        # Verify that the matrix matches the simulation
        checkList = [("meta/SlurmID","SlurmID"),
                     ("meta/simulationOrigo", "simulationOrigo"),
                     ("meta/voxelSize", "voxelSize"),
                     ("meta/hyperVoxelSize", "hyperVoxelSize"),
                     ("meta/configFile", "configFile"),
                     ("meta/positionFile", "positionFile"),
                     ("meta/axonStumpIDFlag", "axonStumpIDFlag")]

        for checkNames in checkList:
          test = self.histFile[checkNames[0]][()] == data[checkNames[1]]

          if(type(test) == bool):
            assert test, "Connection matrix cache - mismatch for " \
              + checkNames[1]
          else:
            assert test.all(), \
              "Connection matrix cache - mismatch for " + checkNames[1]

            
        assert self.nSynapsesTotal == data["nSynapsesTotal"] \
          and self.nGapJunctionsTotal == data["nGapJunctionsTotal"], \
          " Synapse or gap junction count mismatch -- corrupt file?"
      
    else:
      conMat = None
      
    return conMat

  ############################################################################

  def saveConnectionMatrixCache(self,conMat):
    
    cacheFile = self.getConMatCacheFileName()
    self.writeLog("Saving connection matrix cache to " + cacheFile)

    data = dict([])
    data["conMat"] = conMat
    data["SlurmID"] = self.SlurmID
    data["simulationOrigo"] = self.simulationOrigo
    data["voxelSize"] = self.voxelSize
    data["hyperVoxelSize"] = self.hyperVoxelSize
    data["nSynapsesTotal"] = self.nSynapsesTotal
    data["nGapJunctionsTotal"] = self.nGapJunctionsTotal
    data["configFile"] = self.configFile
    data["positionFile"] = self.positionFile
    data["axonStumpIDFlag"] = self.axonStumpIDFlag
    
    with open(cacheFile,'wb') as f:
      pickle.dump(data,f,-1) # -1 latest version

  ############################################################################

  # This checks that all connections included in the pruning, were present
  # in the detection. If some are missing in the detection phase, then they
  # would incorrectly be missing after pruning.
  
  def checkNetworkConfigIntegrity(self):

    detectConfig = json.loads(self.histFile["meta/config"][()])
    with open(self.histFile["meta/configFile"][()],"r") as f:
      pruneConfig = json.load(f)

    allPresent = True
    
    for con in pruneConfig["Connectivity"]:
      if(con not in detectConfig["Connectivity"]):
        self.writeLog("!!! Connection " + con + " has been added to " \
                      + self.histFile["meta/configFile"][()] \
                      + " after detection, please rerun snudda detect")
        allPresent = False

    assert allPresent, "Please rerun snudda detect."
    
    
  ############################################################################

  # Parse the connection information in the config file, stored in the
  # the work history file

  def loadPruningInformation(self):

    # self.config = json.loads(self.histFile["meta/config"][()])
    
    self.checkNetworkConfigIntegrity()
    with open(self.histFile["meta/configFile"][()],"r") as f:
      self.config = json.load(f)    
    
    self.populationUnitID = self.histFile["network/neurons/populationUnitID"][()]

    # Normally we use type names as lookups, but since we will do this
    # many millions of times, we create an temporary typeID number
    self.makeTypeNumbering()
    
    origConnectivityDistributions = \
      json.loads(self.histFile["meta/connectivityDistributions"][()])

    #origConnectivityDistributionsGJ = \
    #  json.loads(self.histFile["meta/connectivityDistributionsGJ"][()])

    self.connectivityDistributions = dict([])

    # For the pruning we merge the two into one
    for keys in origConnectivityDistributions:
      (preType,postType) = keys.split("$$")

      # Need to handle if preType or postType dont exist, then skip this
      if(preType not in self.typeIDLookup \
         or postType not in self.typeIDLookup):
        print("Skipping " + preType + " to " + postType + " connection")
        continue
      
      preTypeID = self.typeIDLookup[preType]
      postTypeID = self.typeIDLookup[postType]

      for conType in origConnectivityDistributions[keys]:
        conData = origConnectivityDistributions[keys][conType]
        
        pruning = self.completePruningInfo(conData["pruning"])
        
        if("pruningOther" in conData):
          pruningOther = self.completePruningInfo(conData["pruningOther"])
        else:
          pruningOther = None

        synapseTypeID = conData["channelModelID"]

        self.connectivityDistributions[preTypeID,postTypeID,synapseTypeID] \
          = (pruning,pruningOther)

  ############################################################################

  # This makes sure all the variables exist, that way pruneSynapseHelper
  # does not have to check, but can assume that they will exist
  
  def completePruningInfo(self,pruneInfo):

    if("distPruning" not in pruneInfo):
      pruneInfo["distPruning"] = None

    if("f1" not in pruneInfo or pruneInfo["f1"] is None):
      pruneInfo["f1"] = 1.0

    if("softMax" not in pruneInfo):
      pruneInfo["softMax"] = None

    if("mu2" not in pruneInfo):
      pruneInfo["mu2"] = None

    if("a3" not in pruneInfo):
      pruneInfo["a3"] = None

    return pruneInfo

  ############################################################################
  
  def loadPruningInformationOLD(self):

    self.config = json.loads(self.histFile["meta/config"][()])

    self.populationUnitID = self.histFile["network/neurons/populationUnitID"][()]
    
    # Normally we use type names as lookups, but since we will do this
    # many millions of times, we create an temporary typeID number
    self.makeTypeNumbering()
    
    self.connectivityDistributions = dict([])

    for name, definition in self.config.items():

      if(name in ["Volume", "PopulationUnits"]):
        # We are just loading the pruning information
        continue

      if "targets" in definition:
        conDist = definition["targets"]
      else:
        conDist = [] 
      
      for row in conDist:
        target = row[0]

        try:
          synapseType = self.synapseTypeReverseLookup[row[1][0]]
          # Skipping loading conductance mean and std,
          # and channelParamDictionary
          distrib = row[2]
        except:
          self.writeLog("Something wrong with your network json file. row: " \
                        + str(row))
          import pdb
          pdb.set_trace()
          
        if(len(row) > 3):
          # If distrib2 is specified, then distrib is within population Unit
          # distribution and distrib2 is between population unit distribution
          distrib2 = row[3]
        else:
          distrib2 = None
        
        typeName = name.split("_")[0]

        try:
          if(target in self.typeIDLookup):
            self.connectivityDistributions[self.typeIDLookup[typeName],
                                           self.typeIDLookup[target],
                                           synapseType] \
                                           = (distrib,distrib2)
          else:
            self.writeLog("WARNING: Target neuron " + str(target) + " not included in simulation")
        except:

          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)
          
          import pdb
          pdb.set_trace()
    
  ############################################################################

  def createConnectionMatrixOLD(self):

    if(self.histFile is None):
      self.openWorkHistoryFile()
    
    nNeurons = self.histFile["network/neurons/neuronID"].shape[0]
    
    self.writeLog("Creating the " +  str(nNeurons) +"x" + str(nNeurons) \
          + " sparse connection matrix") 
    
    # After creation, convert to csr or csc matrix
    conMat = scipy.sparse.lil_matrix((nNeurons,nNeurons),dtype=np.int16)

    assert self.synapses.shape[0] > 0, "No synapses in synapse matrix"
    
    prevSrcID = self.synapses[0,0]
    prevDestID = self.synapses[0,1]
    prevCtr = 1

    # The reason we dont write the synapse directly is that each access
    # to the sparse matrix takes time, so better group synapses
    # between same pairs together to one insert.
    # This requires self.synapess to be sorted.
    
    for row in self.synapses[1:,:]:
      SrcID = row[0]
      DestID = row[1]

      if(SrcID == prevSrcID and DestID == prevDestID):
        prevCtr += 1
      else:
        conMat[prevSrcID,prevDestID] += prevCtr
        prevSrcID = SrcID
        prevDestID = DestID
        prevCtr = 1

    # Dont forget last one
    conMat[prevSrcID,prevDestID] += prevCtr

    self.writeLog("Converting to CSR sparse matrix format")
    return conMat.tocsr()

  ############################################################################

  def makeTypeNumbering(self):

    nameList = [x if type(x) not in [bytes,np.bytes_] else x.decode() \
                for x in self.histFile["network/neurons/name"] ]

    typeNameList = [x.split("_")[0] for x in nameList]

    typeCtr = 1
    self.typeIDLookup = dict([])

    for x in typeNameList:
      if(x not in self.typeIDLookup):
        self.typeIDLookup[x] = typeCtr
        typeCtr += 1

    self.typeIDList = [self.typeIDLookup[x] for x in typeNameList]

    
  ############################################################################

  def setupOutputFile(self,outputFile=None):

    if(self.outFile is not None):
      self.writeLog("Output file already set: " + str(self.outFile.filename))
      return
    
    if(outputFile is None):
      outputFile = self.basePath + "/network-pruned-synapses.hdf5"
    
    self.writeLog("Writing to " + outputFile)

    outFile = h5py.File(outputFile, "w",
                        libver=self.h5libver,
                        driver=self.h5driver)

    
    outFile.create_dataset("config",data=json.dumps(self.config))

    # Copy over meta data
    self.histFile.copy("meta",outFile)
    
    # Morphologies not included at this stage -- maybe add them?
    # Same morphologies
    cfg = json.loads(self.histFile["meta/config"][()])

    morphGroup = outFile.create_group("morphologies")

    for name, definition in cfg["Neurons"].items():
      try:
        morphFile = definition["morphology"]
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)
        import pdb
        pdb.set_trace()
           
      with open(morphFile,"r") as f:
        swcData = f.read()

      self.writeLog("Saving morphology in HDF5 file: " + morphFile)
      swcGroup = morphGroup.create_group(name)
      swcGroup.create_dataset("swc",data=swcData)
      swcGroup.create_dataset("location",data=morphFile)

    
    networkGroup = outFile.create_group("network")

    # Copy over neuron data
    # self.histFile.copy("neurons",outFile)
    self.histFile.copy("network/neurons",networkGroup)
    
    networkGroup.create_dataset("synapses", \
                                dtype=np.int32, \
                                shape = (self.synapseChunkSize,13), \
                                chunks = (self.synapseChunkSize,13), \
                                maxshape=(None,13), \
                                compression=self.h5compression)   

    networkGroup.create_dataset("gapJunctions", \
                                dtype=np.int32, \
                                shape = (self.gapJunctionChunkSize,11), \
                                chunks = (self.gapJunctionChunkSize,11), \
                                maxshape=(None,11), \
                                compression=self.h5compression)   

    nSynapses = np.zeros((1,),dtype=np.uint64)
    nGapJunctions = np.zeros((1,),dtype=np.uint64)

    networkGroup.create_dataset("nSynapses",data=nSynapses,dtype=np.uint64)
    networkGroup.create_dataset("nGapJunctions",data=nGapJunctions,
                                dtype=np.uint64)
    

    self.outFile = outFile
    #self.outFileSynapseCtr = np.int64(0)
    #self.outFileGapJunctionCtr = np.int64(0)
    
  ############################################################################

  def findLatestFile(self):

    #files = glob('save/network_connect_voxel_log-*-worklog.hdf5')
    files = glob('save/network-connect-synapse-voxel-file-*-worklog.hdf5')

    modTime = [os.path.getmtime(f) for f in files]
    idx = np.argsort(modTime)

    self.writeLog("Using the newest file: " + files[idx[-1]])
    
    return files[idx[-1]]


  ############################################################################

  def mergeFileExists(self):

    # check if merge file exists
    mergeFileName = self.basePath + "/network-putative-synapses-MERGED.hdf5"
    
    mergeFileOK = False

    self.writeLog("Checking for merge file " + str(mergeFileName))

    try:
      if(os.path.isfile(mergeFileName)):

        mergeFileOK = True
      
        f = h5py.File(mergeFileName)

        # Check that SlurmID matches
        if(f["meta/SlurmID"][()] != self.SlurmID):
          mergeFileOK = False

        # Check that synapse matrix has right size
        if(f["network/synapses"].shape[0] != self.nSynapsesTotal):
          mergeFileOK = False

        # Check that last row is not empty
        if((f["network/synapses"][-1,:] == 0).all()):
          mergeFileOK = False

        if("gapJunctions" not in f["network"]):
          mergeFileOK = False
        elif(f["network/gapJunctions"].shape[0] != self.nGapJunctionsTotal):
          mergeFileOK = False
        elif((f["network/gapJunctions"][-1,:] == 0).all()):
          # Last row not set, not complete
          mergeFileOK = False
                
        # !!! Add tests for gap junctions also
        #import pdb
        #pdb.set_trace()
      else:
        f = None
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)

      self.writeLog("Something went wrong with reading old merge file, ignoring it")
      # Something went wrong
      mergeFileOK = False
      
    if(mergeFileOK):
      self.writeLog("Found old merge file " + str(mergeFileName))
      return (True,f["network/synapses"],f)
    else:
      # Bad merge file, close it
      if(f is not None):
        f.close()
        
      return (False,None,None)
    
  ############################################################################
  
  def setupMergeFile(self,verbose=False,bigCache=False,
                     outFileName=None,saveMorphologies=True,
                     nSynapses = None,
                     nGapJunctions = None,
                     deleteAfter=True):
    
    
    if(outFileName is None):
      outFileName = self.basePath + "/network-putative-synapses-MERGED.hdf5"

    #  Make a list of all temporary files so we can remove them
    if(deleteAfter):
      self.tempFileList.append(outFileName)
      
    self.writeLog("Setting up out file " + str(outFileName))
    if(bigCache):
      # !!! Special code to increase h5py cache size
      propfaid = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
      settings = list(propfaid.get_cache())
      # print(settings)
      # [0, 521, 1048576, 0.75]

      settings[2] *= 20
      # propfaid.set_cache(*settings)
      # settings = propfaid.get_cache()
      # print(settings)
      # (0, 521, 5242880, 0.75)


      # Low level opening hdf5 file, to have greater cache size
      fid = h5py.h5f.create(outFileName.encode(), \
                            flags=h5py.h5f.ACC_TRUNC, \
                            fapl=propfaid)
      outFile = h5py.File(fid, libver=self.h5libver,driver=self.h5driver)      
    else:
      outFile = h5py.File(outFileName, "w", libver=self.h5libver,
                          driver=self.h5driver)

    networkGroup = outFile.create_group("network")
    
    # Copy over meta data
    self.histFile.copy("meta",outFile)

    # Copy over neuron data
    #self.histFile.copy("neurons",outFile)
    self.histFile.copy("network/neurons",networkGroup)


    cfg = json.loads(self.histFile["meta/config"][()])


    # Save morphologies
    if(saveMorphologies):
      morphGroup = outFile.create_group("morphologies")

      for name, definition in cfg["Neurons"].items():
        try:
          morphFile = definition["morphology"]
        except:
          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)
          import pdb
          pdb.set_trace()
           
        with open(morphFile,"r") as f:
          swcData = f.read()

        self.writeLog("Saving morphology in HDF5 file: " + morphFile)
        swcGroup = morphGroup.create_group(name)
        swcGroup.create_dataset("swc",data=swcData)
        swcGroup.create_dataset("location",data=morphFile)
      
    chunkSize = self.synapseChunkSize      

    if(nSynapses is None):
      nSynapses = self.nSynapsesTotal

    if(nGapJunctions is None):
      nGapJunctions = self.nGapJunctionsTotal


    # !!! We are allocating a lot more space than we might be needing
    # The matrices are resized when the file write buffer is flushed
    # in bufferMergeWrite (flush=True)
      
    #import pdb
    #pdb.set_trace()
    if(nSynapses > chunkSize):
      networkGroup.create_dataset("synapses", \
                                  dtype=np.int32, \
                                  shape = (nSynapses,13), \
                                  chunks = (chunkSize,13), \
                                  maxshape=(nSynapses,13), \
                                  compression=self.h5compression)
    elif(nSynapses > 0):
      networkGroup.create_dataset("synapses", \
                                  dtype=np.int32, \
                                  shape = (nSynapses,13),
                                  chunks = (nSynapses,13))
    else:
      networkGroup.create_dataset("synapses", \
                                  dtype=np.int32, \
                                  shape = (nSynapses,13))

    if(nGapJunctions > chunkSize):
      networkGroup.create_dataset("gapJunctions", \
                                  dtype=np.int32, \
                                  shape = (nGapJunctions,11), \
                                  chunks = (chunkSize,11), \
                                  maxshape=(nGapJunctions,11), \
                                  compression=self.h5compression)
    elif(nGapJunctions > 0):
      networkGroup.create_dataset("gapJunctions", \
                                  dtype=np.int32, \
                                  shape = (nGapJunctions,11),
                                  chunks = (nGapJunctions,11))
    else:       
      networkGroup.create_dataset("gapJunctions", \
                                  dtype=np.int32, \
                                  shape = (nGapJunctions,11))
    
    self.nextMergeFileID += 1

    # Reset the position for the merge file
    self.nextFileWritePos = 0
    
    return(outFile,outFileName)
          
  ############################################################################

  def pruneSynapsesParallel(self,synapseFile,
                            outputFile=None,
                            mergeDataType="synapses",
                            setupOutFile=True,
                            closeInputFile=True):

    h5SynMat, h5SynN, h5SynLoc = self.dataLoc[mergeDataType]

    if(synapseFile[h5SynMat].shape[0] == 0):
      self.writeLog("pruneSynapsesParallel: No " + mergeDataType \
                    + " skipping pruning")
      return
    else:
      self.writeLog("pruneSynapsesParallel, before pruning : " \
                    + str(synapseFile[h5SynMat].shape[0]) + " "\
                    + mergeDataType)
    
    startTime = timeit.default_timer()
    
    if(self.dView is not None and self.role == "master"):
      self.setupParallel(dView=self.dView)

      try:
        # Make sure all synapse writes are on the disk
        synapseFile.flush()
        if(not synapseFile.swmr_mode):
          synapseFile.swmr_mode = True # Allow multiple readers from the file
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)    
        import pdb
        pdb.set_trace() 
               

    # 1. Pick names for the workers
    tempOutputFileName = [self.scratchPath + "worker-temp-" \
                          + mergeDataType + "-file-" + str(x) \
                          for x in range(0,len(self.dView))]
    self.dView.scatter("outputFileName",tempOutputFileName,block=True)

    # Add the files to a delete list, so we remove them after
    for f in tempOutputFileName:
      self.tempFileList.append(f)
    
    # 2. Define what ranges each worker should do
    synapseRanges = self.findRanges(synapseFile[h5SynMat],
                                    len(self.dView))

    if(synapseRanges is None or synapseRanges[-1][-1] is None):
      self.writeLog("There are few synapses, we will run it in serial instead")
      return self.pruneSynapses(synapseFile=synapseFile,
                                outputFileName=None,rowRange=None,
                                closeOutFile=False,
                                closeInputFile=closeInputFile,
                                mergeDataType=mergeDataType)

    # We got valid synapseRanges, continue
    
    self.dView.scatter("synapseRange",synapseRanges,block=True)
    self.writeLog("synapseRanges: " + str(synapseRanges))
    
    self.dView.push({"synapseFileName":synapseFile.filename},block=True)
    self.dView.push({"mergeDataType":mergeDataType},block=True)    

    fn = synapseFile.filename
    
    # Close the synapse file on the master node
    if(closeInputFile):
      synapseFile.close()
      synapseFile = None
    
    # 3. Let the workers prune
    self.writeLog("Sending pruning job to workers")    

    cmdStr = "nw.pruneSynapses(synapseFile=synapseFileName,outputFileName=outputFileName[0],rowRange=synapseRange[0],mergeDataType=mergeDataType)"

    self.dView.execute(cmdStr,block=True)

    endTime = timeit.default_timer()
    self.writeLog("Parallel pruning duration: " + str(endTime-startTime))

    
    # 4. Merge the resulting files -- this is easier, since synapses
    #    are already sorted in the file
    self.writeLog("Merging parallel pruning results")

    if(setupOutFile):
      assert self.outFile is None, \
        "pruneSynapsesParallel: Output file already setup"
      self.setupOutputFile(outputFile=outputFile)
    else:
      # If there were no synapses, but there are gap junctions
      # then this assert could theoretically happen
      assert self.outFile is not None, \
        "pruneSynapsesParallel: Out file not set up"

    tmpFiles = [h5py.File(f,'r') for f in tempOutputFileName]

    nSyn = np.sum(f[h5SynMat].shape[0] for f in tmpFiles)
    matWidthAll = [f[h5SynMat].shape[1] for f in tmpFiles]

    try:
    
      assert (np.array(matWidthAll) == matWidthAll[0]).all(), \
        "Internal error, width does not match"
      matWidth = matWidthAll[0]
    
      self.outFile[h5SynMat].resize((nSyn,matWidth))
      
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)     
      import pdb
      pdb.set_trace()
      
    nextSyn = 0
    
    for f in tmpFiles:
      
      n = f[h5SynMat].shape[0]

      if(n > 0):
        self.outFile[h5SynMat][nextSyn:(nextSyn+n),:] = \
          f[h5SynMat][()]
      
        nextSyn += n

      f.close()

    self.outFile["network/" + h5SynN][0] = nextSyn

    endTime2 = timeit.default_timer()
    self.writeLog("Parallel pruning + merging: " + str(endTime2-startTime))

    

  ############################################################################

  # Find which ranges of the synapse matrix that each worker should take care of
  
  def findRanges(self,synapses,nWorkers,startPos=0,nSyn=None):
    
    if(nSyn is None):
      nSyn = synapses.shape[0] - startPos

    blockSize = max(1,int(math.floor(float(nSyn)/nWorkers)))

    self.writeLog("Find block ranges. From " + str(startPos) \
                  + " to " + str(nSyn+startPos) \
                  + " block size " + str(blockSize))

    # We keep blockStart in a list, so we can remove if they overlap
    blockStart = [x+startPos for x \
                  in range(0,blockSize*(nWorkers+1),blockSize)]

    rangeEnd = startPos + nSyn
    
    # There can be multiple synapses between a cell pair, those must belong
    # to the same block since pruning is dependent on the number of synapses
    # between a pair of neurons
    for idx in range(1,len(blockStart)):
      startRow = blockStart[idx]
      
      while(startRow < rangeEnd
            and (synapses[blockStart[idx]-1,0:2] \
               == synapses[startRow,0:2]).all()):
        startRow += 1

      blockStart[idx] = startRow

    # If we have few synapses, and many workers, there might be some workers
    # assigned to the same range. Remove those, and pad with None

    blockStart = [x for x in blockStart if x <= rangeEnd]
    blockStart = list(collections.OrderedDict.fromkeys(blockStart))

    while(len(blockStart) < nWorkers+1):
      blockStart.append(None)

    synapseRange = [x for x in zip(blockStart[:-1],blockStart[1:])]

    self.writeLog("synapseRange=" + str(synapseRange))

    
    return synapseRange
  
  
  ############################################################################
   
  def findRangesOLD(self,synapses,nWorkers,startPos=0,nSyn=None):

    if(nSyn is None):
      nSyn = synapses.shape[0] - startPos

    blockSize = max(1,int(math.floor(float(nSyn)/nWorkers)))

    self.writeLog("Find block ranges. From " + str(startPos) \
                  + " to " + str(nSyn+startPos) \
                  + " block size " + str(blockSize))

    try:
      blockStart = np.array([x+startPos for x \
                             in range(0,blockSize*nWorkers,blockSize)])
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()
      

    self.writeLog("blockStart=" + str(blockStart))
    
    for idx in range(1,len(blockStart)):

      assert blockStart[idx-1] < blockStart[idx], \
        "Blocks seem to be quite small, run in serial instead. " \
        + "idx = " + str(idx) + " blockStart = " + str(blockStart)
      
      startRow = blockStart[idx]

      # Need to make sure that the synapses between pairs of neurons are
      # not split between workers
      try:
        while((synapses[blockStart[idx]-1,0:2] \
               == synapses[startRow,0:2]).all()):
          startRow += 1
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)

        self.writeLog("Problem with setting ranges, fall back to serial")
        return None
          
      blockStart[idx] = startRow

    assert (np.diff(blockStart) > 0).all(), \
      "All workers should have different block starts. Try running in serial."
      
    # Calculate the corresponding block ends, special treatment for last pos
    blockEnd = np.zeros(blockStart.shape,dtype=int)
    blockEnd[0:-1] = blockStart[1:]
    blockEnd[-1] = startPos + nSyn + 1 # +1, python dont include element at end 

    synapseRanges = [range(x,y) for x,y in zip(blockStart,blockEnd)]

    self.writeLog("synapseRanges = " + str(synapseRanges))
    
    return synapseRanges

    
  ############################################################################

  def cleanUpMergeFiles(self):

    if self.tempFileList is None:
      # Nothing to do
      return
    
    self.writeLog("Cleaning up old merge files")
    
    for f in self.tempFileList:
      # self.writeLog("Removing old merge file: " + str(f))
      try:
        if(os.path.exists(f)):
          os.remove(f)
      except:
        print("Closing of file failed: " + str(f))

    self.tempFileList = None

  ############################################################################
      
  def writeLog(self,text,flush=True): # Change flush to False in future, debug
    try:
      if(self.logFile is not None):
        self.logFile.write(text + "\n")
        print(text)
        if(flush):
          self.logFile.flush()
      else:
        if(self.verbose):
          print(text)
    except:
      print(text)
      print("Unable to write to log file. Is log file closed?")

  ############################################################################

  def setSeed(self,randomSeed):

    self.writeLog("Setting random seed: " + str(randomSeed))
    np.random.seed(randomSeed)

  ############################################################################
    
  def newWorkerSeeds(self,dView):

    nWorkers = len(self.dView)
    workerSeeds = np.random.randint(0,np.iinfo(np.uint32).max,
                                    dtype=np.uint32,
                                    size=(nWorkers,))
    self.dView.scatter("workerSeed",workerSeeds,block=True)
    self.dView.execute("nw.setSeed(workerSeed[0])",block=True)

    self.writeLog("New worker seeds: " + str(workerSeeds))
  
  ############################################################################
  
  def setupParallel(self, dView):

    assert self.role == "master", \
      "setupParallel: Should only be called by master node"

    if(dView is None):
      self.writeLog("setupParallel called without dView, aborting.")
      return
   
    if(self.workersInitialised):
      self.writeLog("Workers already initialised.")
      return
   
    with dView.sync_imports():
      from snudda.prune import SnuddaPrune

    self.writeLog("Setting up workers: " \
                  + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    
    # Create unique log file names for the workers
    if(self.logFileName is not None):
      engineLogFile = [self.logFileName + "-" \
                       + str(x) for x in range(0,len(dView))]
    else:
      engineLogFile = [[] for x in range(0,len(dView))]
      
    dView.scatter('logFileName',engineLogFile,block=True)

    dView.push({"workHistoryFile":self.workHistoryFile},
               block=True)

    cmdStr = "nw = SnuddaPrune(workHistoryFile=workHistoryFile, logFileName=logFileName[0],role='worker')"

    dView.execute(cmdStr,block=True)

    # Make sure we have different seeds for workers
    self.newWorkerSeeds(dView)
    
    self.writeLog("Workers setup: " \
                  + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))


    self.workersInitialised = True

############################################################################

  def lowMemory(self,threshold=0.1):
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

    #import pdb
    #pdb.set_trace()

    memoryRatio = ret['free'] / ret['total']

    # self.writeLog("Memory status: " + str(int(memoryRatio * 100)) + "% free")
    
    return memoryRatio < threshold

  ############################################################################

  # This code does a N-way merge for synapses

  # !!! We also need to merge gap junctions
  
  def bigMerge(self,mergeDataType="synapses",cleanVoxelFiles=True):

    assert False, "bigMerge is replaced by bigMergeLookup"
    
    # Since we want the code to work for both synapses and gap junction
    # we need to know location of synapse matrix, eg "network/synapses",
    # number of synapses, eg "nSynapses", and in some functions which
    # columns containe x,y,z voxel coords, eg. range(2,5)
    h5SynMat, h5SynN, h5SynLoc = self.dataLoc[mergeDataType]
    
    self.mergeDataType = mergeDataType
    self.writeLog("Doing bigMerge for " + mergeDataType)
    
    synapseHeap = []
    
    self.writeLog("Starting big merge of all data")

    nNeurons = len(self.histFile["network/neurons/neuronID"])
    assert np.max(self.histFile["network/neurons/neuronID"])+1 == nNeurons, \
      "bigMerge: There are neuron IDs missing"
    
    maxHyperID = np.max(self.allHyperIDs)+1
    fileList = [None] * maxHyperID
    numSynapses = np.zeros((maxHyperID,),dtype=np.int)
    
    # Open all files for reading
    hFileNameMask = self.basePath + "/voxels/network-putative-synapses-%s.hdf5"

    maxAxonVoxelCtr = 0
    maxDendVoxelCtr = 0

    nSynHist = self.histFile[h5SynN]
    nSynTotal = np.sum(nSynHist)

    for hID,nSyn,nOverflow in zip(self.histFile["completed"], nSynHist,
                                  self.histFile["voxelOverflowCounter"]):

        
      if(nSyn > 0):
        # Open file, and add first row to the heap
        hFileName = hFileNameMask % str(hID)
        fileList[hID] = h5py.File(hFileName,'r')
        numSynapses[hID] = nSyn

        srcID = fileList[hID][h5SynMat][0,0]
        destID = fileList[hID][h5SynMat][0,1]
        uniqueID = destID*nNeurons + srcID # This is used for heap priority

        # Create a heap containing the first element of all files
        heapq.heappush(synapseHeap,(uniqueID,hID))
            
        # This is so we can optimize the axon/dend voxelCtr and size
        if("maxAxonVoxelCtr" in fileList[hID]["meta"]):
          maxAxonVoxelCtr = max(maxAxonVoxelCtr,
                                fileList[hID]["meta/maxAxonVoxelCtr"][()])
        if("maxDendVoxelCtr" in fileList[hID]["meta"]):
          maxDendVoxelCtr = max(maxDendVoxelCtr,
                                fileList[hID]["meta/maxDendVoxelCtr"][()])

        if(cleanVoxelFiles):
          # This makes sure we remove the old voxel files afterwards
          self.tempFileList.append(hFileName)

          
    # Setup read buffer
    self.setupBufferedMergeRead(fileList,h5SynMat)

    if(self.bufferOutFile is None):
      # Create output file
      (self.bufferOutFile,outFileName) \
        = self.setupMergeFile(deleteAfter=False)
    else:
      # We need to reset the write pointer (GJ and synapses should start from 0)
      self.nextFileWritePos = 0

    # Only save this meta data if doing the synapses call
    if(maxAxonVoxelCtr > 0 and self.mergeDataType == "synapses"):
      self.bufferOutFile["meta"].create_dataset("maxAxonVoxelCtr",
                                                data=maxAxonVoxelCtr)
      self.writeLog("maxAxonVoxelCtr = " + str(maxAxonVoxelCtr))

    if(maxDendVoxelCtr > 0 and self.mergeDataType == "synapses"):
      self.bufferOutFile["meta"].create_dataset("maxDendVoxelCtr",
                                                data=maxDendVoxelCtr)
      self.writeLog("maxDendVoxelCtr = " + str(maxDendVoxelCtr))
    
    # 2. Pop the smallest element, and add it to the final file
    # -- check same file if there are more with the same source and dest
    # -- buffer the writes

    synCtr = 0

    if(len(synapseHeap) > 0):
      # Get the first file to read synapses from
      (uniqueID,hID) = heapq.heappop(synapseHeap)
    else:
      # No synapses at all, return
      self.cleanUpMergeReadBuffers()

      return (self.bufferOutFile[h5SynMat],self.bufferOutFile)
      
    # synapses coordinates are translated to simulation wide coordinates
    # from hyper voxel local coordinates
    (synapses,synapsesRemaining,nextSynapsePair) \
      = self.findPairSynapsesBufferedMerge(hID,h5SynMat,h5SynLoc)
    
    self.bufferMergeWrite(h5SynMat,synapses)
    synCtr += synapses.shape[0]

    # !!! DEBUG
    prevPair = synapses[0,0:2].copy()
    prevUID = uniqueID
    
    loopCtr = 0
    while(synapsesRemaining or len(synapseHeap) > 0):

      if(loopCtr % 1000000 == 0):
        self.writeLog("Synapses: " + str(synCtr) \
                      + "/" + str(nSynTotal) \
                      + " (heap size: " + str(len(synapseHeap)) + ")")
      
      if(synapsesRemaining):
        # Need to push the next synapse in line from that file onto the heap
        [srcID,destID] = nextSynapsePair

        uniqueID = destID*nNeurons + srcID # This is used for heap priority
        (uniqueID,hID) = heapq.heappushpop(synapseHeap,(uniqueID,hID))
      else:
        (uniqueID,hID) = heapq.heappop(synapseHeap)
                
      (synapses,synapsesRemaining,nextSynapsePair) \
        = self.findPairSynapsesBufferedMerge(hID,h5SynMat,h5SynLoc)

      try:
        assert nextSynapsePair is None or \
          synapses[0,1] < nextSynapsePair[1] \
          or (synapses[0,1] == nextSynapsePair[1] \
              and synapses[0,0] < nextSynapsePair[0]), \
              "They are not in order in file"
        
      except:
        print("This is realyl weird...")
        import pdb
        pdb.set_trace()
      
      assert uniqueID == synapses[0,1]*nNeurons+synapses[0,0], \
        "Oh no! Internal inconsistency"

      # !!! Just a check while writing code, to make sure it is ok
      assert (synapses[:,0] == synapses[0,0]).all() \
        and (synapses[:,1] == synapses[0,1]).all(), \
        "Source and dest should all be the same for pair synapses"
      
      try:
        assert synapses[0,1] > prevPair[1] \
          or (synapses[0,1] == prevPair[1] \
              and synapses[0,0] >= prevPair[0]), \
              "Synapses not in order"
      except:
        print("We have a big problem! Synapses not in order!")
        import pdb
        pdb.set_trace() 

      try:
        self.bufferMergeWrite(h5SynMat, synapses)
        synCtr += synapses.shape[0]
        loopCtr += 1
        
      except:
        print("Why is this wrong?")
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)    
        import pdb
        pdb.set_trace() 
       

      # DEBUG!!!
      prevPair = synapses[0,0:2].copy()
      prevUID = uniqueID

      
    # Make sure all of the buffer is written to file
    self.bufferMergeWrite(h5SynMat,flush=True)
    
    # The buffers should be cleared once the files are read from, 
    # but just as a precaution
    self.cleanUpMergeReadBuffers()

    return (self.bufferOutFile[h5SynMat],self.bufferOutFile)

  ############################################################################

  def findPairSynapses(self,synapses,startIdx):
    
    endIdx = startIdx + 1
    
    while(endIdx < synapses.shape[0] and
          synapses[startIdx,0] == synapses[endIdx,0] \
          and synapses[startIdx,1] == synapses[endIdx,1]):
      endIdx += 1

    if(endIdx < synapses.shape[0]):
      synapsesRemaining = True
    else:
      synapsesRemaining = False
      
    return (synapses[startIdx:endIdx,:],endIdx,synapsesRemaining)

  ############################################################################

  def setupBufferedMergeRead(self,fileList,synapseMatrixPath):

    self.fileList = fileList
    self.fileBuffers = [None] * len(fileList)
    self.nextBufferReadPos = np.zeros((len(fileList),),dtype=int)
    self.nextFileReadPos = np.zeros((len(fileList),),dtype=int)
    
    for idx, f in enumerate(fileList):
      if(f is not None):
        nSyn = f[synapseMatrixPath].shape[0]
        bufLen = min(nSyn,self.synapseBufferSize)
        self.fileBuffers[idx] = f[synapseMatrixPath][0:bufLen,:]
        self.nextFileReadPos[idx] = bufLen

  ############################################################################
  
  def updateMergeReadBuffer(self,hID, synapseMatrixPath):
    
    if(self.fileBuffers[hID] is None):
      # Nothing to do
      return

    if(self.fileList[hID] is None):
      self.fileBuffers[hID] = None
      return

    syn = self.fileList[hID][synapseMatrixPath]
    nSynToRead = min(syn.shape[0] - self.nextFileReadPos[hID],
                   self.synapseBufferSize)

    if(nSynToRead > 0):
      startPos = self.nextFileReadPos[hID]
      endPos = startPos + nSynToRead
      
      if(endPos - startPos != self.fileBuffers[hID].shape[0]):
        self.fileBuffers[hID] = syn[startPos:endPos,:].copy()
      else:
        # Reuse memory
        self.fileBuffers[hID][:,:] = syn[startPos:endPos,:]

      self.nextFileReadPos[hID] = endPos        
      self.nextBufferReadPos[hID] = 0

    else:
      self.fileBuffers[hID] = None
      self.nextBufferReadPos[hID] = 0        

  ############################################################################

  def cleanUpMergeReadBuffers(self):

    self.mergeDataType = None
    
    for hid, f in enumerate(self.fileList):
      
      if f is not None:
        try:
          f.close()
        except:
          self.writeLog("Problem closing file for HID: " + str(hid))

    if(self.fileBuffers is not None):
      for idx in range(0,len(self.fileBuffers)):
        self.fileBuffers[idx] = None

    self.fileBuffers = None
    self.nextFileReadPos = None
    self.nextBufferReadPos = None
      
  ############################################################################

  # Input: Hyper voxel ID, columns which contains location
  # (columns are different for synapses: range(2,5) and GJ range(6,9)
  # so we take them as parameter to reuse code for both cases

  # Returns: synapses between pair,
  #          synapses remaining flag,
  #          (srcID,destID) of next synapse
  #
  # OBS, voxel coordinates for synapses are translated from hyper voxel
  # coordiantes, to simulation wide coordinates
  
  def findPairSynapsesBufferedMerge(self,hID,synapseMatrixPath,locationColumns):
    
    startIdx = self.nextBufferReadPos[hID]
    endIdx = startIdx + 1

    synapses = self.fileBuffers[hID]
    bufLen = synapses.shape[0]

    oldSynapses = None
        
    while(True):

      if(endIdx >= bufLen):

        # End of buffer, we need to read more synapses from file
        assert oldSynapses is None, \
          "Buffer too small, should not overflow twice." #self.synapseBufferSize

        oldSynapses = synapses[startIdx:,:].copy()

        # This resets self.nextBufferReadPos[hID]
        self.updateMergeReadBuffer(hID,synapseMatrixPath) 

        startIdx = 0
        endIdx = 1

        if(self.fileBuffers[hID] is None):
          # This was the end of the file
          return (oldSynapses,False,None)        

        synapses = self.fileBuffers[hID]        
        bufLen = synapses.shape[0]
        
        assert oldSynapses.shape[0] > 0, "This should never be empty"          

        if(oldSynapses[0,0] != synapses[0,0] \
           or oldSynapses[0,1] != synapses[0,1]):
          # New pair in new buffer, send old pair synapses to caller
          return (oldSynapses,True,synapses[0,:2].copy())

      # Check if the next synapse belongs to pair  
      if(synapses[startIdx,0] != synapses[endIdx,0] \
          or synapses[startIdx,1] != synapses[endIdx,1]):
        # No more similar, stop
        break
        
      endIdx += 1

          
    # Update the buffer read pos
    self.nextBufferReadPos[hID] = endIdx
          
    if(endIdx < bufLen):
      synapsesRemaining = True
      nextSynapsePair = synapses[endIdx,0:2].copy()
    else:
      synapsesRemaining = False
      nextSynapsePair = None

          
    if(oldSynapses is not None):
      syn = np.concatenate((oldSynapses,synapses[startIdx:endIdx,:]))
      # Need to add first half also from previous buffer read
    else:
      syn = synapses[startIdx:endIdx,:]

    # OBS, we need to convert to global voxel coordinates from the local
    # hyper voxel coordinates

    # !!! We need to convert from hyper voxel local coordinates for synapse
    # to simulation wide coordinates
    # !!! THIS IS NOW DONE AFTER TOUCH DETECTION
    # syn[:,locationColumns] += self.hyperVoxelOffset[hID]

    return (syn,synapsesRemaining,nextSynapsePair)
       
  ############################################################################
  
  def bufferMergeWrite(self,synapseMatrixLoc,synapses=None,flush=False):

#    if(synapseMatrixLoc == "network/gapJunctions" \
#       and flush):
#      import pdb # !!DEBUG
#      pdb.set_trace()
    
    if(self.synapseWriteBuffer is None):
      # First time run
      
      if(synapses is not None):
        bufferWidth = synapses.shape[1]
      else:
        self.writeLog("bufferMergeWrite: Unable to guess width of matrix")
        import pdb
        pdb.set_trace()
      
      self.synapseWriteBuffer = np.zeros((self.synapseBufferSize,
                                          bufferWidth),
                                         dtype=np.int32)
      self.nextBufferWritePos = 0
      #self.nextFileWritePos = 0

      
    if(synapses is not None):

      # Is buffer almost full?
      if(synapses.shape[0] + self.nextBufferWritePos >= self.synapseBufferSize):
        # Flush buffer first
        endFileIdx = self.nextFileWritePos + self.nextBufferWritePos
        self.bufferOutFile[synapseMatrixLoc][self.nextFileWritePos:endFileIdx,:] \
          = self.synapseWriteBuffer[0:self.nextBufferWritePos,:]

        self.nextFileWritePos += self.nextBufferWritePos
        self.nextBufferWritePos = 0
        # Ok, now we have clean buffer, back to normal program

      endBufIdx = self.nextBufferWritePos + synapses.shape[0]
      self.synapseWriteBuffer[self.nextBufferWritePos:endBufIdx,:] = synapses
      self.nextBufferWritePos += synapses.shape[0]

    if(flush):
      self.writeLog("Flushing " + str(self.bufferOutFile.filename) \
                    + " data: " + synapseMatrixLoc )
      endFileIdx = self.nextFileWritePos + self.nextBufferWritePos
      self.bufferOutFile[synapseMatrixLoc][self.nextFileWritePos:endFileIdx,:] \
        = self.synapseWriteBuffer[0:self.nextBufferWritePos,:]

      self.nextFileWritePos += self.nextBufferWritePos
      self.nextBufferWritePos = 0

      # Resize matrix to fit data
      w = self.bufferOutFile[synapseMatrixLoc].shape[1]
      self.bufferOutFile[synapseMatrixLoc].resize((endFileIdx,w))
      self.writeLog(synapseMatrixLoc + " new size " + str((endFileIdx,w)))
      
      # Remove buffer
      self.synapseWriteBuffer = None

############################################################################

  def bigMergeParallel(self):

    if(self.role != "master"):
      self.writeLog("bigMergeParallel is only run on master node, aborting")
      return
  
    self.writeLog("bigMergeParallel, starting " + str(self.role))

    if(self.dView):
      self.setupParallel(dView=self.dView)

    # Split neurons between nodes, we need the neurons to be in order
    nNeurons = self.histFile["network/neurons/neuronID"].shape[0]
    assert nNeurons-1 == self.histFile["network/neurons/neuronID"][-1], \
      "neuronID should start from 0 and the end should be n-1"
    
    nWorkers = len(self.dView)

    neuronRanges = []
    rangeBorders = np.linspace(0,nNeurons,nWorkers+1).astype(int)
    
    for idx in range(0,nWorkers):
      neuronRanges.append((rangeBorders[idx],rangeBorders[idx+1]))

    assert neuronRanges[-1][-1] == nNeurons, \
      "bigMergeParallel: Problem with neuronRanges, last element incorrect"
    assert len(neuronRanges) == nWorkers, \
      "bigMergeParallel: Problem with neuronRanges, bad length"
    
    # Send list of neurons to workers
    self.dView.scatter("neuronRange",neuronRanges,block=True)

    # Each worker sorts a subset of the neurons and write it to separate files
    cmdStrSyn = "mergeResultSyn = nw.bigMergeHelper(neuronRange=neuronRange[0],mergeDataType='synapses')"
    
    self.dView.execute(cmdStrSyn,block=True)
    mergeResultsSyn = self.dView["mergeResultSyn"]

    # When we do scatter, it embeds the result in a list
    cmdStrGJ = "mergeResultGJ = nw.bigMergeHelper(neuronRange=neuronRange[0],mergeDataType='gapJunctions')"
    self.dView.execute(cmdStrGJ,block=True)
    mergeResultsGJ = self.dView["mergeResultGJ"]

    # We need to sort the files in order, so we know how to add them
    print("Processing merge...")
    
    mergeStartSyn = [x[1][0] for x in mergeResultsSyn]
    mergeStartGJ  = [x[1][0] for x in mergeResultsGJ]

    mergeOrderSyn = np.argsort(mergeStartSyn)
    mergeOrderGJ = np.argsort(mergeStartGJ)
        
    # We then need the file to merge them in
    (self.bufferOutFile,outFileName) \
      = self.setupMergeFile(bigCache=False,deleteAfter=False)

    #import pdb
    #pdb.set_trace()
    
    # Copy the data to the file

    for order, results, location \
        in zip([mergeOrderSyn, mergeOrderGJ],
               [mergeResultsSyn, mergeResultsGJ],
               ["network/synapses", "network/gapJunctions"]):

      startPos = 0
      
      for idx in order:
        
        dataFile = results[idx][0]
        nSynapses = results[idx][2]
        endPos = startPos + nSynapses

        if(dataFile is None):
          assert nSynapses == 0, "!!! Missing merge file " + str(dataFile) \
            + ", internal problem"
          continue
          
        self.writeLog("Extracting " + location + " from " + dataFile)
      
        fIn = h5py.File(dataFile,'r')

        # Some idiot checks...
        assert fIn[location].shape[0] == nSynapses, \
        "Internal inconsistency in number of rows stored"
        assert not (fIn[location][-1,:] == 0).all(), \
        "Last row all zero, that should not happen"
        
        self.bufferOutFile[location][startPos:endPos,:] = \
          fIn[location][:,:]

        fIn.close()
        startPos = endPos

    # We should check that the number of synapses in output file is correct
    
    assert self.bufferOutFile["network/synapses"].shape[0] \
      == np.sum(self.histFile["nSynapses"][:]), \
      "Not all synapses kept in merge, internal problem."
    assert self.bufferOutFile["network/gapJunctions"].shape[0] \
      == np.sum(self.histFile["nGapJunctions"][:]), \
      "Not all gap junctions kept in merge, internal problem."

    self.writeLog("bigMergeParallel: done")

    return self.bufferOutFile
      
    
############################################################################

  def getHyperVoxelList(self,neuronRange):
    
    # We need to find all the hypervoxels that might contain synapses.
    # During touch detection we generated a list for each hypervoxel with
    # all the neurons that had any part of it within that hypervoxel.
    
    hvList = []
    neuronSet = set(range(neuronRange[0],neuronRange[1]))
    
    for hid in self.histFile["hyperVoxels"]:
      hvNeurons = self.histFile["hyperVoxels"][hid]["neurons"]
      if(len(neuronSet.intersection(hvNeurons)) > 0):
        hvList.append(int(hid))

    return hvList

    
############################################################################

  # Needs to handle both gap junctions and synapses

  def bigMergeHelper(self,neuronRange,mergeDataType):

   # !!! TEMP
   #import uuid
   #tmpFile = open("delme-log-" + str(uuid.uuid4()), "w")
   #tmpFile.write("TEST!!")
   #tmpFile.close()
   # !!! TMP
   
   try:
    self.writeLog("Neuronrange = " + str(neuronRange))
     
    outputFileName = self.scratchPath + mergeDataType + "-for-neurons-" \
      + str(neuronRange[0]) + "-to-" + str(neuronRange[1]) \
      + "-MERGE-ME.hdf5"
    
    self.mergeDataType = mergeDataType
    
    # Which hyper voxels are the neurons located in?
    hvList = self.getHyperVoxelList(neuronRange)
    
    # Locations of data within the file
    h5SynMat, h5SynN, h5SynLookup = self.dataLoc[mergeDataType]

    synapseHeap = []
    fileList = dict([])
    fileMatIterator = dict([])
    
    # Next we need to open all the relevant files
    hFileNameMask = self.basePath + "/voxels/network-putative-synapses-%s.hdf5"

    maxAxonVoxelCtr = 0
    maxDendVoxelCtr = 0

    # !!! Special code to increase h5py cache size
    propfaid = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
    settings = list(propfaid.get_cache())
    settings[2] *= 20
    propfaid.set_cache(*settings)
    #settings = propfaid.get_cache()
    # !!! End special cache code

    nHV = int(self.histFile["nCompleted"][0])
    nTotal = 0

    #import pdb
    #pdb.set_trace()
    
    for hID,nSyn,nOverflow in zip(self.histFile["completed"][:nHV],
                                  self.histFile[h5SynN][:nHV],
                                  self.histFile["voxelOverflowCounter"][:nHV]):

      nTotal += nSyn
      
      if(hID not in hvList):
        # We only need a subset of the files
        continue

      if(nSyn > 0):
        hFileName = hFileNameMask % str(hID)
        self.writeLog("Opening voxel file: " + hFileName)
        
        # Low level opening hdf5 file, to have greater cache size #ACC_RDWR
        fid = h5py.h5f.open(hFileName.encode(), \
                            flags=h5py.h5f.ACC_RDONLY, \
                            fapl=propfaid)

        # !!! Temp print to check cache size
        settings = list(fid.get_access_plist().get_cache())
        print(settings)
        
        # fileList[hID] = h5py.File(hFileName,'r')
        fileList[hID] = h5py.File(fid,drive=self.h5driver)

        chunkSize = 10000
        lookupIterator = \
          self.fileRowLookupIteratorSubset(h5matLookup=fileList[hID][h5SynLookup],
                                           minDestID=neuronRange[0],
                                           maxDestID=neuronRange[1],
                                           chunkSize=chunkSize)
                                           
        
        fileMatIterator[hID] \
          = self.synapseSetIterator(h5matLookup=fileList[hID][h5SynLookup],
                                    h5mat=fileList[hID][h5SynMat],
                                    chunkSize=chunkSize,
                                    lookupIterator=lookupIterator)
        
        synSet, uniqueID = next(fileMatIterator[hID],(None,None))

        if(synSet is None):
          # There were synapses in the hyper voxel, but none relevant to our
          # selected files. Clear file List and fileMatIterator for this worker
          del fileList[hID]
          del fileMatIterator[hID]
          continue
        
        # Create a heap containing the first subset of all files
        heapq.heappush(synapseHeap,(uniqueID,hID,synSet))

       # This is so we can optimize the axon/dend voxelCtr and size
        if("maxAxonVoxelCtr" in fileList[hID]["meta"]):
          maxAxonVoxelCtr = max(maxAxonVoxelCtr,
                                fileList[hID]["meta/maxAxonVoxelCtr"][()])
        if("maxDendVoxelCtr" in fileList[hID]["meta"]):
          maxDendVoxelCtr = max(maxDendVoxelCtr,
                                fileList[hID]["meta/maxDendVoxelCtr"][()])

    if(mergeDataType == "synapses"):
      nSynapses = nTotal
      nGapJunctions = 0
    elif(mergeDataType == "gapJunctions"):
      nSynapses = 0
      nGapJunctions = nTotal
    else:
      assert False, "Unknown mergeDataType " + str(mergeDataType)
      
    # Setup output file
    (self.bufferOutFile,outFileName) \
      = self.setupMergeFile(bigCache=True,
                            outFileName=outputFileName,
                            saveMorphologies=False,
                            nSynapses = nSynapses,
                            nGapJunctions = nGapJunctions)

    # Here we store the sorted connection matrix
    sortedMat = self.bufferOutFile[h5SynMat]

    # Only save this meta data if doing the synapses call
    if(maxAxonVoxelCtr > 0 and self.mergeDataType == "synapses"):
      self.bufferOutFile["meta"].create_dataset("maxAxonVoxelCtr",
                                                data=maxAxonVoxelCtr)
      self.writeLog("maxAxonVoxelCtr = " + str(maxAxonVoxelCtr))

    if(maxDendVoxelCtr > 0 and self.mergeDataType == "synapses"):
      self.bufferOutFile["meta"].create_dataset("maxDendVoxelCtr",
                                                data=maxDendVoxelCtr)
      self.writeLog("maxDendVoxelCtr = " + str(maxDendVoxelCtr))

    # Take the first (smallest uniqueID) element from the heap
    
    if(len(synapseHeap) > 0):
      # Get the first file to read synapses from
      (uniqueID,hID,synSet) = heapq.heappop(synapseHeap)
    else:
      # No synapses at all, return
      self.cleanUpMergeReadBuffers()
      return (None,neuronRange,0)

    # Store synapse
    self.bufferMergeWrite(h5SynMat,synSet)
    synCtr = synSet.shape[0]    

    loopCtr = 0
    done = False

    while(not done):

      oldUniqueID = uniqueID
      
      if(loopCtr % 1000000 == 0):
        self.writeLog("Worker synapses: " + str(synCtr) \
                      + "/" + str(nTotal) \
                      + " (heap size: " + str(len(synapseHeap)) + ")")

      # Get the next set of synapses from this file from the iterator
      nextRowSet = next(fileMatIterator[hID],None)
        
      if(nextRowSet is not None):
        # More synapses in file, push next pair to heap, and pop top pair
        synSet,uniqueID = nextRowSet         
        (uniqueID,hID,synSet) = heapq.heappushpop(synapseHeap,
                                                  (uniqueID, hID,
                                                   synSet))
      elif(len(synapseHeap) > 0):
        (uniqueID,hID,synSet) = heapq.heappop(synapseHeap)

      else:
        done = True
        continue

      try:
        assert uniqueID >= oldUniqueID, "uniqueID should be increasing in file"
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)

        import pdb
        pdb.set_trace()
      
        
      self.bufferMergeWrite(h5SynMat,synSet)
      synCtr += synSet.shape[0]
      loopCtr += 1


    self.writeLog("Worker synapses: " + str(synCtr) \
                  + "/" + str(nTotal) \
                  + " (heap size: " + str(len(synapseHeap)) + ")")

    self.writeLog("Read " + str(synCtr) + " out of total " + str(nTotal) \
                  + " synapses")
      
    self.bufferMergeWrite(h5SynMat,flush=True)

    
    self.writeLog("bigMergeHelper: done")


    # Close the hyper voxel files
    for f in fileList:
      try:
        fileList[f].close()
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)

        print("Problems closing files")
        import pdb
        pdb.set_trace()
        
    
    self.bufferOutFile.close()
    
    return (outputFileName, neuronRange, synCtr)
   except:
    import traceback
    tstr = traceback.format_exc()
    self.writeLog(tstr)
    self.writeToRandomFile(tstr)
    import pdb
    pdb.set_trace() 
    
     
############################################################################

  def bigMergeLookup(self,mergeDataType="synapses",cleanVoxelFiles=True):

    # Since we want the code to work for both synapses and gap junction
    # we need to know location of synapse matrix, eg "network/synapses",
    # number of synapses, eg "nSynapses", the lookup table to quickly
    # find which synapse rows belongs to each pair of connected neurons
    # eg "network/synapseLookup"
    h5SynMat, h5SynN, h5SynLookup = self.dataLoc[mergeDataType]
    
    self.mergeDataType = mergeDataType
    self.writeLog("Doing bigMerge (lookup) for " + mergeDataType)
    
    synapseHeap = []

    nNeurons = len(self.histFile["network/neurons/neuronID"])
    assert np.max(self.histFile["network/neurons/neuronID"])+1 == nNeurons, \
      "bigMerge (lookup): There are neuron IDs missing"

    maxHyperID = np.max(self.allHyperIDs)+1
    fileList = [None] * maxHyperID
    fileMatIterator = [None] * maxHyperID
    
    #fileMat = [None] * maxHyperID # points to the synapse matrix in each file
    #fileMatLookup = [None] * maxHyperID # points to the matrix lookup in file
    
    numSynapses = np.zeros((maxHyperID,),dtype=np.int)
    
    # Open all files for reading
    hFileNameMask = self.basePath + "/voxels/network-putative-synapses-%s.hdf5"

    maxAxonVoxelCtr = 0
    maxDendVoxelCtr = 0

    nSynHist = self.histFile[h5SynN]
    nSynTotal = np.sum(nSynHist)

    # !!! Special code to increase h5py cache size
    propfaid = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
    settings = list(propfaid.get_cache())
    print(settings)
    # [0, 521, 1048576, 0.75]

    settings[2] *= 20
    propfaid.set_cache(*settings)
    settings = propfaid.get_cache()
    print(settings)
    # (0, 521, 5242880, 0.75)

    
    for hID,nSyn,nOverflow in zip(self.histFile["completed"], nSynHist,
                                  self.histFile["voxelOverflowCounter"]):

      hFileName = hFileNameMask % str(hID)
      
      if(cleanVoxelFiles):
        # This makes sure we remove the old voxel files afterwards
        self.tempFileList.append(hFileName)

      
      if(nSyn > 0):
        # Open file, and add info about first pairs synapses to the heap
        self.writeLog("Opening voxel file: " + hFileName)
        
        # Low level opening hdf5 file, to have greater cache size #ACC_RDWR
        fid = h5py.h5f.open(hFileName.encode(), \
                            flags=h5py.h5f.ACC_RDONLY, \
                            fapl=propfaid)

        # !!! Temp print to check cache size
        settings = list(fid.get_access_plist().get_cache())
        print(settings)
        
        # fileList[hID] = h5py.File(hFileName,'r')
        try:
          fileList[hID] = h5py.File(fid,drive=self.h5driver)
          fileMatIterator[hID] \
            = self.synapseSetIterator(h5matLookup=fileList[hID][h5SynLookup],
                                      h5mat=fileList[hID][h5SynMat],
                                      chunkSize=10000)
        except:
          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)

          print("This should not happen...")
          import pdb
          pdb.set_trace()
          
        numSynapses[hID] = nSyn

        # There should be at least the first row, otherwise nSyn = 0
        synSet, uniqueID = next(fileMatIterator[hID],None)
          
        # Create a heap containing the first subset of all files
        heapq.heappush(synapseHeap,(uniqueID,hID,synSet))

        # This is so we can optimize the axon/dend voxelCtr and size
        if("maxAxonVoxelCtr" in fileList[hID]["meta"]):
          maxAxonVoxelCtr = max(maxAxonVoxelCtr,
                                fileList[hID]["meta/maxAxonVoxelCtr"][()])
        if("maxDendVoxelCtr" in fileList[hID]["meta"]):
          maxDendVoxelCtr = max(maxDendVoxelCtr,
                                fileList[hID]["meta/maxDendVoxelCtr"][()])

          
    assert np.sum(numSynapses) == nSynTotal, \
      "Mismatch between work log file and data files: " \
      + str(nSynTotal) + " vs " + str(np.sum(numSynapses)) + " synapses"
          
    if(self.bufferOutFile is None):
      # Create output file
      (self.bufferOutFile,outFileName) \
        = self.setupMergeFile(bigCache=True,deleteAfter=False)
    else:
      # We need to reset the write pointer (GJ and synapses should start from 0)
      self.nextFileWritePos = 0

    # Here we store the sorted connection matrix
    sortedMat = self.bufferOutFile[h5SynMat]
      
    # Only save this meta data if doing the synapses call
    if(maxAxonVoxelCtr > 0 and self.mergeDataType == "synapses"):
      self.bufferOutFile["meta"].create_dataset("maxAxonVoxelCtr",
                                                data=maxAxonVoxelCtr)
      self.writeLog("maxAxonVoxelCtr = " + str(maxAxonVoxelCtr))

    if(maxDendVoxelCtr > 0 and self.mergeDataType == "synapses"):
      self.bufferOutFile["meta"].create_dataset("maxDendVoxelCtr",
                                                data=maxDendVoxelCtr)
      self.writeLog("maxDendVoxelCtr = " + str(maxDendVoxelCtr))

    # 2. Pop the smallest element, and add it to the final file
    # -- check same file if there are more with the same source and dest
    # -- buffer the writes

    synCtr = 0

    if(len(synapseHeap) > 0):
      # Get the first file to read synapses from
      (uniqueID,hID,synSet) = heapq.heappop(synapseHeap)
    else:
      # No synapses at all, return
      self.cleanUpMergeReadBuffers()
      return (sortedMat,self.bufferOutFile)

    # Store synapses
    # synEndIdx = synCtr + synSet.shape[0]
    # sortedMat[synCtr:synEndIdx,:] = synSet
    # synCtr = synEndIdx
    synCtr = synSet.shape[0]
    self.bufferMergeWrite(h5SynMat,synSet)

    loopCtr = 0
    done = False
    
    while(not done):
      
      if(loopCtr % 1000000 == 0):
        self.writeLog("Synapses: " + str(synCtr) \
                      + "/" + str(nSynTotal) \
                      + " (heap size: " + str(len(synapseHeap)) + ")")

      # Get the next set of synapses from this file from the iterator
      nextRowSet = next(fileMatIterator[hID],None)
      
      if(nextRowSet is not None):
        # More synapses in file, push next pair to heap, and pop top pair
        synSet,uniqueID = nextRowSet         
        (uniqueID,hID,synSet) = heapq.heappushpop(synapseHeap,
                                                  (uniqueID, hID,
                                                   synSet))
      elif(len(synapseHeap) > 0):
        (uniqueID,hID,synSet) = heapq.heappop(synapseHeap)

      else:
        done = True
        continue
      
      # Write synapses to file
      # synEndIdx = synCtr + synSet.shape[0]
      # sortedMat[synCtr:synEndIdx,:] = synSet
      # synCtr = synEndIdx
      self.bufferMergeWrite(h5SynMat,synSet)
      synCtr += synSet.shape[0]
      loopCtr += 1
      
      #assert uniqueID == synapses[0,1]*nNeurons+synapses[0,0], \
      #  "bigMergeLookup: Oh no! Internal inconsistency"

    # Flush the buffers to file
    #self.bufferOutFile.flush()
    self.bufferMergeWrite(h5SynMat,flush=True)
    
    self.writeLog("bigMergeLookup: done")
    
    return (sortedMat,self.bufferOutFile)

############################################################################

  def bigMergeLookupNOCACHE(self,mergeDataType="synapses"):

    # Since we want the code to work for both synapses and gap junction
    # we need to know location of synapse matrix, eg "network/synapses",
    # number of synapses, eg "nSynapses", the lookup table to quickly
    # find which synapse rows belongs to each pair of connected neurons
    # eg "network/synapseLookup"
    h5SynMat, h5SynN, h5SynLookup = self.dataLoc[mergeDataType]
    
    self.mergeDataType = mergeDataType
    self.writeLog("Doing bigMerge (lookup) for " + mergeDataType)
    
    synapseHeap = []

    nNeurons = len(self.histFile["network/neurons/neuronID"])
    assert np.max(self.histFile["network/neurons/neuronID"])+1 == nNeurons, \
      "bigMerge (lookup): There are neuron IDs missing"

    maxHyperID = np.max(self.allHyperIDs)+1
    fileList = [None] * maxHyperID
    fileMat = [None] * maxHyperID # points to the synapse matrix in each file
    fileMatLookup = [None] * maxHyperID # points to the matrix lookup in file
    
    nextLookupIdx = np.zeros((maxHyperID,),dtype=np.int)
    lookupLength = np.zeros((maxHyperID,),dtype=np.int)
    numSynapses = np.zeros((maxHyperID,),dtype=np.int)
    
    # Open all files for reading
    hFileNameMask = self.basePath + "/voxels/network-putative-synapses-%s.hdf5"

    maxAxonVoxelCtr = 0
    maxDendVoxelCtr = 0

    nSynHist = self.histFile[h5SynN]
    nSynTotal = np.sum(nSynHist)

    # !!! Special code to increase h5py cache size
    propfaid = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
    settings = list(propfaid.get_cache())
    print(settings)
    # [0, 521, 1048576, 0.75]

    settings[2] *= 20
    propfaid.set_cache(*settings)
    settings = propfaid.get_cache()
    print(settings)
    # (0, 521, 5242880, 0.75)

    
    for hID,nSyn,nOverflow in zip(self.histFile["completed"], nSynHist,
                                  self.histFile["voxelOverflowCounter"]):
      if(nSyn > 0):
        # Open file, and add info about first pairs synapses to the heap
        hFileName = hFileNameMask % str(hID)
        self.writeLog("Opening voxel file: " + hFileName)
        
        # Low level opening hdf5 file, to have greater cache size #ACC_RDWR
        fid = h5py.h5f.open(hFileName.encode(), \
                            flags=h5py.h5f.ACC_RDONLY, \
                            fapl=propfaid)

        # !!! Temp print to check cache size
        settings = list(fid.get_access_plist().get_cache())
        print(settings)
        
        # fileList[hID] = h5py.File(hFileName,'r')
        fileList[hID] = h5py.File(fid,drive=self.h5driver)
        fileMat[hID] = fileList[hID][h5SynMat]
        fileMatLookup[hID] = fileList[hID][h5SynLookup]
        numSynapses[hID] = nSyn

        try:
          uniqueID, startIdx, endIdx = fileMatLookup[hID][0,:]
          nextLookupIdx[hID] = 1
          lookupLength[hID] = fileMatLookup[hID].shape[0]
        except:
          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)

          print("No what happened...")
          import pdb
          pdb.set_trace()
          
        # Create a heap containing the first element of all files
        heapq.heappush(synapseHeap,(uniqueID,hID,startIdx,endIdx))

        # This is so we can optimize the axon/dend voxelCtr and size
        if("maxAxonVoxelCtr" in fileList[hID]["meta"]):
          maxAxonVoxelCtr = max(maxAxonVoxelCtr,
                                fileList[hID]["meta/maxAxonVoxelCtr"][()])
        if("maxDendVoxelCtr" in fileList[hID]["meta"]):
          maxDendVoxelCtr = max(maxDendVoxelCtr,
                                fileList[hID]["meta/maxDendVoxelCtr"][()])

    assert np.sum(numSynapses) == nSynTotal, \
      "Mismatch between work log file and data files: " \
      + str(nSynTotal) + " vs " + str(np.sum(numSynapses)) + " synapses"
          
    if(self.bufferOutFile is None):
      # Create output file
      (self.bufferOutFile,outFileName) \
        = self.setupMergeFile(bigCache=True)
    else:
      # We need to reset the write pointer (GJ and synapses should start from 0)
      self.nextFileWritePos = 0

    # Here we store the sorted connection matrix
    sortedMat = self.bufferOutFile[h5SynMat]
      
    # Only save this meta data if doing the synapses call
    if(maxAxonVoxelCtr > 0 and self.mergeDataType == "synapses"):
      self.bufferOutFile["meta"].create_dataset("maxAxonVoxelCtr",
                                                data=maxAxonVoxelCtr)
      self.writeLog("maxAxonVoxelCtr = " + str(maxAxonVoxelCtr))

    if(maxDendVoxelCtr > 0 and self.mergeDataType == "synapses"):
      self.bufferOutFile["meta"].create_dataset("maxDendVoxelCtr",
                                                data=maxDendVoxelCtr)
      self.writeLog("maxDendVoxelCtr = " + str(maxDendVoxelCtr))

    # 2. Pop the smallest element, and add it to the final file
    # -- check same file if there are more with the same source and dest
    # -- buffer the writes

    synCtr = 0

    if(len(synapseHeap) > 0):
      # Get the first file to read synapses from
      (uniqueID,hID,startIdx,endIdx) = heapq.heappop(synapseHeap)
    else:
      # No synapses at all, return
      self.cleanUpMergeReadBuffers()
      return (sortedMat,self.bufferOutFile)

    # Maybe store ref direct to synapses, avoiding one level of lookups
    synapses = fileMat[hID][startIdx:endIdx,:]
    synEndIdx = synCtr + synapses.shape[0]
    
    sortedMat[synCtr:synEndIdx,:] = synapses
    synCtr = synEndIdx

    loopCtr = 0
    done = False
    
    while(not done):
      
      if(loopCtr % 10000 == 0):
        self.writeLog("Synapses: " + str(synCtr) \
                      + "/" + str(nSynTotal) \
                      + " (heap size: " + str(len(synapseHeap)) + ")")

      if(nextLookupIdx[hID] < lookupLength[hID]):
        # More synapses in file, push next pair to heap, and pop top pair
        [uniqueID,startIdx,endIdx] = \
          fileMatLookup[hID][nextLookupIdx[hID],:]
        nextLookupIdx[hID] += 1
        
        (uniqueID,hID,startIdx,endIdx) = heapq.heappushpop(synapseHeap,
                                                           (uniqueID, hID,
                                                            startIdx, endIdx))
      elif(len(synapseHeap) > 0):
        (uniqueID,hID,startIdx,endIdx) = heapq.heappop(synapseHeap)

      else:
        done = True
        continue
      
      # Write synapses to file
      synapses = fileMat[hID][startIdx:endIdx,:]
      synEndIdx = synCtr + (endIdx-startIdx) #synapses.shape[0]
    
      sortedMat[synCtr:synEndIdx,:] = synapses
      synCtr = synEndIdx
      loopCtr += 1
      
      #assert uniqueID == synapses[0,1]*nNeurons+synapses[0,0], \
      #  "bigMergeLookup: Oh no! Internal inconsistency"

    # Flush the buffers to file
    self.bufferOutFile.flush()

    self.writeLog("bigMergeLookup NOCACHE: done")
    
    return (sortedMat,self.bufferOutFile)


############################################################################

  def pruneSynapses(self,synapseFile,outputFileName,rowRange,
                    mergeDataType,
                    closeInputFile=True,
                    closeOutFile=True):      

    try:
      h5SynMat, h5SynN, h5SynLoc = self.dataLoc[mergeDataType]
    
      if(rowRange is None):
        rowStart = 0
        rowEnd = synapseFile[h5SynMat].shape[0]
      else:
        rowStart = rowRange[0]
        rowEnd = rowRange[-1] 
      
      if(rowStart is None or rowEnd is None):
        self.writeLog("Nothing to do, empty row range")
        return
    
      self.writeLog("pruneSynapses called.")  
    
      if(type(synapseFile) != str and synapseFile[h5SynMat].shape[0] == 0):
        self.writeLog("pruneSynapses: No " + mergeDataType \
                      + " skipping pruning") 
        return

    
      self.writeLog("pruneSynapses: synapseFile=" + str(synapseFile) \
                    + ", outputFileName=" + str(outputFileName) \
                    + ", rowRange=" + str(rowRange) \
                    + " (" + mergeDataType + ")")

      if(type(synapseFile) == str):
        self.writeLog("Opening synapse file: " + synapseFile)
        if(self.role != "master"):
          # SWMR = one writer, multiple readers
          synapseFile = h5py.File(synapseFile,'r',swmr=True)        
        else:
          synapseFile = h5py.File(synapseFile,'r')

      nSyn = rowEnd - rowStart

      # We need to split the rowRange into smaller pieces that fit in memory    
      chunkSize = 1000000

      # To avoid division by zero
      nBlocks = max(1,int(np.ceil(float(rowEnd-rowStart)/chunkSize)))

      self.writeLog("About to calculate block ranges (" \
                    + str(nBlocks) + " blocks)")
    
      blockRanges = self.findRanges(synapses=synapseFile[h5SynMat],
                                    nWorkers=nBlocks,
                                    startPos=rowStart,
                                    nSyn=nSyn)

      self.writeLog("blockRanges="+str(blockRanges))
    
      self.setupOutputFile(outputFileName) # Sets self.outFile
    
      for synRange in blockRanges:
        self.writeLog("Pruning range : " + str(synRange))
      
        synapses = synapseFile[h5SynMat][synRange[0]:synRange[-1]]
        self.pruneSynapsesHelper(synapses=synapses,
                                 outputFile=self.outFile,
                                 mergeDataType=mergeDataType)

    except:
      print("pruneSynapses: Something went wrong... :/")
      
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      
      import pdb
      pdb.set_trace()

        
    # Close synapse input file
    if(closeInputFile):
      synapseFile.close()
      
    if(closeOutFile):
      self.outFile.close()
      self.outFile = None

  ############################################################################

  # This code prunes synapses, but you need to send it small chunks
  # so synapses are all in memory, also it buffers the write until the end

  # synapses -- subset of synapse matrix that fits in memory
  # rowRange -- which rows to read
  # outputFile -- where to write synapses, assumed to already exist
  # outFilePos -- which position to start writing from
  
  def pruneSynapsesHelper(self,synapses,outputFile,mergeDataType):

    h5SynMat, h5SynN, h5SynLoc = self.dataLoc[mergeDataType]

    keepRowFlag = np.zeros((synapses.shape[0],),dtype=bool)

    nextReadPos = 0
    readEndOfRange = synapses.shape[0]

    # Init some stats
    nAllRemoved = 0
    nSomeRemoved = 0
    nTooFewRemoved = 0
    nDistDepPruning = 0
    nTooManyRemoved = 0
    nNotConnected = 0

    oldPos = -1
    
    while(nextReadPos < readEndOfRange):

      if(oldPos == nextReadPos):
        print("pruneSynapsesHelper: Same again")
        import pdb
        pdb.set_trace()
        
      oldPos = nextReadPos
      
      # How many lines contain synapses between this pair of neurons
      readEndIdx = nextReadPos+1
      while(readEndIdx < readEndOfRange and \
            (synapses[nextReadPos,0:2] == synapses[readEndIdx,0:2]).all()):
        readEndIdx += 1

      # Temp check
      assert (synapses[nextReadPos:readEndIdx,0] \
              == synapses[nextReadPos,0]).all() and \
              (synapses[nextReadPos:readEndIdx,1] \
               == synapses[nextReadPos,1]).all(), \
               "pruneSynapsesHelper: Internal error, more than one neuron pair"

      #import pdb
      #pdb.set_trace()

      
      # Stats
      nPairSynapses = readEndIdx - nextReadPos

      srcID = synapses[nextReadPos,0]
      destID = synapses[nextReadPos,1]

      if(mergeDataType == "gapJunctions"):
        # All are gap junctions
        synapseType = 3
      else:
        synapseType = synapses[nextReadPos,6]
      
      conID = (self.typeIDList[srcID],self.typeIDList[destID],synapseType)
      
      if(conID in self.connectivityDistributions):

        # We have the option to separate between connections within a
        # population unit or not. If conInfo[1] != None then first
        # tuple is connection info within a population unit, and second item
        # is connection info between different population units
        conInfo = self.connectivityDistributions[conID]

        #
        if(conInfo[1] is None \
           or self.populationUnitID[srcID] == self.populationUnitID[destID]):
          # All or within population unit pruning parameters
          cInfo = conInfo[0]
        else:
          # Between population unit pruning parameters
          cInfo = conInfo[1]

        # These will always exist thanks to completePruningInfo function
        
        distP = cInfo["distPruning"] # Dist dep pruning
        f1 = cInfo["f1"]
        softMax = cInfo["softMax"]
        mu2 = cInfo["mu2"]
        a3 = cInfo["a3"]
          
      else:
        # Not listed in connectivityDistribution, skip neuron pair
        nextReadPos = readEndIdx
        # No need to update keepRowFlag since default set to 0
        
        # Stats
        nNotConnected += nPairSynapses
        continue
      
      # 3. This is the last step of pruning, but we move it to the top
      # since there is no point doing the other steps if we going to
      # throw them away anyway
      if(a3 is not None and np.random.random() > a3):
        # Prune all synapses between pair, do not add to synapse file
        nextReadPos = readEndIdx
        # No need to update keepRowFlag since default set to 0

        # Stats
        nAllRemoved += nPairSynapses
        continue
      
      if(distP is not None):
        # Distance dependent pruning, used for FS->MS connections

        # distP contains d (variable for distance to soma)
        d = synapses[nextReadPos:readEndIdx,8] * 1e-6 # dendrite distance
        P = eval(distP)
        fracFlag = np.random.random(nPairSynapses) < f1
        distFlag = np.random.random(nPairSynapses) < P
        
        keepRowFlag[nextReadPos:readEndIdx] = np.logical_and(fracFlag,distFlag)
     
        nFrac = sum(fracFlag)
        nSomeRemoved += nPairSynapses - nFrac
        nDistDepPruning += nFrac - sum(keepRowFlag[nextReadPos:readEndIdx])

      else:
        keepRowFlag[nextReadPos:readEndIdx] \
          = np.random.random(nPairSynapses) < f1
        nSomeRemoved += nPairSynapses - sum(keepRowFlag[nextReadPos:readEndIdx])
       
      # Check if too many synapses, trim it down a bit
      nKeep = np.sum(keepRowFlag[nextReadPos:readEndIdx])

      if(softMax is not None and nKeep > softMax):

        # pKeep = float(softMax)/nKeep # OLD implementation
        softMax = float(softMax)
        #pKeep = 2*softMax*np.divide(1-np.exp(-nKeep/softMax),1+np.exp(-nKeep/softMax))/nKeep
        pKeep = np.divide(2*softMax,(1+np.exp(-(nKeep-softMax)/5))*nKeep)
        
        keepRowFlag[nextReadPos:readEndIdx] = \
          np.logical_and(pKeep > np.random.random(nPairSynapses),
                         keepRowFlag[nextReadPos:readEndIdx])

        # Stats
        nTooManyRemoved += nKeep - sum(keepRowFlag[nextReadPos:readEndIdx])
        
        # Update count
        nKeep = np.sum(keepRowFlag[nextReadPos:readEndIdx])
       

      # If too few synapses, remove all synapses
      if(mu2 is not None):
        Pmu = 1.0/(1.0 + np.exp(-8.0/mu2 *(nKeep - mu2)))

        if(Pmu < np.random.random()):
          # Too few synapses, remove all -- need to update keepRowFlag
          keepRowFlag[nextReadPos:readEndIdx] = 0
          nextReadPos = readEndIdx
          
          # Stats
          nTooFewRemoved += nKeep
          continue

      nextReadPos = readEndIdx      

    # Time to write synapses to file
    nKeepTot = sum(keepRowFlag)
    writeStartPos = int(outputFile["network/" + h5SynN][0])
    writeEndPos = writeStartPos + nKeepTot
        
    if(nKeepTot > 0):
      outputFile[h5SynMat].resize((writeEndPos,outputFile[h5SynMat].shape[1]))
      outputFile[h5SynMat][writeStartPos:writeEndPos] = \
        synapses[keepRowFlag,:]

      # Update counters
      outputFile["network/" + h5SynN][0] = writeEndPos

    else:
      self.writeLog("No synapses kept, resizing")
      outputFile[h5SynMat].resize((writeEndPos,outputFile[h5SynMat].shape[1]))
      
    self.writeLog("Number of synapses removed where synapse connection not allowed: " + str(nNotConnected) \
                  + "\nNumber of synapses removed due to distance dependent pruning: " + str(nDistDepPruning) \
                  + "\nNumber of synapses removed randomly: " + str(nSomeRemoved) \
                  + "\nNumber of synapses removed due to too many synapses between connected pair: " + str(nTooManyRemoved) \
                  + "\nNumber of synapses removed due to too few synapses between connected pairs: " + str(nTooFewRemoved) \
                  + "\nNumber of synapses removed where all synapses between pairs are removed: " \
                  + str(nAllRemoved))

############################################################################
      
  def fileRowLookupIterator(self,h5mat,chunkSize=10000):

    matSize = h5mat.shape[0]

    # If matrix is smaller than chunk, buffer can be smaller than requested
    if(matSize < chunkSize):
      chunkSize = h5mat.shape[0]

    matBuf = np.zeros((chunkSize,h5mat.shape[1]),dtype=h5mat.dtype)
    endIdx = 0
    
    while(endIdx < matSize):
      startIdx = endIdx
      endIdx = startIdx+chunkSize
      
      if(endIdx < matSize or matSize == chunkSize):
        # Copy to existing buffer, to avoid memory allocation
        matBuf[:,:] = h5mat[startIdx:endIdx,:]
      else:
        # Create a new buffer
        matBuf = h5mat[startIdx:endIdx,:]
      
      for row in matBuf:
        yield row

############################################################################

  # minDestID and maxDestID are inclusive, only synapses with destID in that
  # range are iterated over

  def fileRowLookupIteratorSubset(self,h5matLookup,
                                  minDestID,maxDestID,
                                  chunkSize=10000):

    nNeurons = self.histFile["network/neurons/neuronID"].shape[0]
    minUniqueID = minDestID * nNeurons
    maxUniqueID = maxDestID * nNeurons 
    # minUniqueID <= destID < maxUniqueID

    self.writeLog("minUniqueID: " + str(minUniqueID) \
                  + ", maxUniqueID: " + str(maxUniqueID))
    
    matSize = h5matLookup.shape[0]

    # If matrix is smaller than chunk, buffer can be smaller than requested
    if(matSize < chunkSize):
      chunkSize = matSize

    matBuf = np.zeros((chunkSize, h5matLookup.shape[1]),
                      dtype=h5matLookup.dtype)
    endIdx = 0
    
    while(endIdx < matSize):
      startIdx = endIdx
      endIdx = startIdx+chunkSize
      
      if(endIdx < matSize or matSize == chunkSize):
        # Copy to existing buffer, to avoid memory allocation
        matBuf[:,:] = h5matLookup[startIdx:endIdx,:]
      else:
        # Create a new buffer
        matBuf = h5matLookup[startIdx:endIdx,:]
      
      for row in matBuf:
        if(minUniqueID <= row[0] and row[0] < maxUniqueID):
          # Only return synapses that terminate on the neurons we are
          # interested in here
          
          yield row

        
############################################################################

  # Returns (subset of rows, uniqueID)
  

  def synapseSetIterator(self,h5matLookup,h5mat,
                         chunkSize=10000,
                         lookupIterator=None):

    # Allow the user to set an alternative lookupIterator if we only
    # want to iterate over a subset of the synapses
    if(not lookupIterator):
      lookupIterator = self.fileRowLookupIterator(h5matLookup,
                                                  chunkSize=chunkSize)

    matSize = h5mat.shape[0]
    if(matSize < chunkSize):
      chunkSize = matSize

    oldSynapses = None
    # readBuffer = np.zeros((chunkSize,h5mat.shape[1]),dtype=h5mat.dtype)
    readBuffer = h5mat[:chunkSize,:].copy()
    bufferStart = 0 # What file pos does start of buffer correspond to
    bufferEnd = chunkSize # What file pos does end of buffer correspond to (+1)
    
    nextRowSet = next(lookupIterator,None)
    
    while(nextRowSet is not None):

      # startIdx and endIdx are the rows in the matrix we want to read between
      [uniqueID,startIdx,endIdx] = nextRowSet

      if(startIdx >= bufferEnd):
        # We need to jump forward...
        bufferStart = startIdx

        if(startIdx + chunkSize <= h5mat.shape[0]):
          bufferEnd = startIdx+chunkSize
          readBuffer[:,:] = h5mat[bufferStart:bufferEnd,:]
        else:
          bufferEnd = h5mat.shape[0]
          readBuffer = h5mat[bufferStart:bufferEnd,:].copy()
          
        oldSynapses = None
      
      if(endIdx > bufferEnd):
        assert oldSynapses is None, "getNextSynapseSet: chunkSize too small"

        # Part of the synapse range requested is outside buffer
        oldSynapses = readBuffer[(startIdx-bufferStart):,:].copy()

        bufferStart = bufferEnd
        bufferEnd = bufferStart + chunkSize
        
        if(bufferEnd > matSize):
          bufferEnd = matSize
          readBuffer = h5mat[bufferStart:bufferEnd,:].copy()
        else:
          # Reuse old buffer storage
          readBuffer[:,:] = h5mat[bufferStart:bufferEnd,:]

        # Need to concatenate with old synapses
        synMat = np.concatenate([oldSynapses,
                                 readBuffer[:(endIdx-bufferStart),:]],
                                axis=0)

        try:
          assert endIdx == bufferEnd \
            or (readBuffer[startIdx-bufferStart,:2] \
                != readBuffer[endIdx-bufferStart,:2]).any(), \
            "We missed one synpase! (2)"
          
          assert (synMat[:,0] == synMat[0,0]).all() \
            and (synMat[:,1] == synMat[0,1]).all(), \
            "Synapse matrix (2) contains more than one pair:\n" + str(synMat)

          assert synMat.shape[0] == endIdx - startIdx, \
            "Synapse matrix has wrong size"
        except:
          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)

          import pdb
          pdb.set_trace()
          
          
        yield (synMat,
               uniqueID)

        oldSynapses = None
          
      else:

        synMat = readBuffer[(startIdx-bufferStart):(endIdx-bufferStart),:]

        try:
          assert endIdx == bufferEnd \
            or (readBuffer[startIdx-bufferStart,:2] \
                != readBuffer[endIdx-bufferStart,:2]).any(), \
            "We missed one synpase! (1)"

          
          assert (synMat[:,0] == synMat[0,0]).all() \
            and (synMat[:,1] == synMat[0,1]).all(), \
            "Synapse matrix (1) contains more than one pair:\n" + str(synMat)

          assert synMat.shape[0] == endIdx - startIdx, \
            "Synapse matrix has wrong size"
        except:
          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)

          import pdb
          pdb.set_trace()
          
        
        yield (synMat,
               uniqueID)

      nextRowSet = next(lookupIterator,None)
      
##############################################################################

if __name__ == "__main__":

  print("Please do not call this file directly, use snudda.py")
  exit(-1)
