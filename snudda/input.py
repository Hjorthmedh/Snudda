# This code writes the input spikes for the neuron simulation --
#
#
# If nInputs is given then synapseDensity is scaled to give approximately
# that total number of synapses, otherwise it is used without scaling.
# see config/input-tinytest-v2.json for example config.
#

#
# !!!! Change how data is stored, many small datasets is inefficient
#

# Smith, Galvan, ..., Bolam 2014 -- Bra info om thalamic inputs, CM/PF
#


import numpy as np
import h5py
import json
import os
from glob import glob
import itertools

import matplotlib.pyplot as plt

from .Neuron_morphology import NeuronMorphology
from .load import SnuddaLoad

nl = None

class SnuddaInput(object):
  
  def __init__(self, spikeDataFileName, inputConfigFile,
               networkConfigFile=None,
               positionFile=None,
               HDF5networkFile=None,
               time=10.0,
               isMaster=True,
               h5libver="latest",
               randomSeed=None,
               logFile=None,
               verbose=True):
    
    if(type(logFile) == str):
      self.logFile = open(logFile,"w")
    else:
      self.logFile = logFile
      
    self.verbose = verbose
    
    self.writeLog("Time = " + str(time))

    # We need to set the seed, to avoid same seed on workers
    np.random.seed(randomSeed)
    self.writeLog("Setting random seed: " + str(randomSeed))
    
    self.h5libver = h5libver
    self.writeLog("Using hdf5 version " + str(h5libver))
    
    if(HDF5networkFile is not None):
      assert networkConfigFile is None and positionFile is None, \
        "If HDF5networkFile specified then positionFile " + \
        "and networkConfigFile should be left empty."

      if(HDF5networkFile == "last"):
        HDF5networkFile = self.findLatestFile()
      
      self.HDF5networkFile = HDF5networkFile
      self.readHDF5info(HDF5networkFile)
    else:
      self.networkConfigFile = networkConfigFile
      self.positionFile = positionFile

      #self.writeLog("Assuming axonStumpIDFlag is True (Running Network_simulate.py)")
      self.axonStumpIDFlag = False
    
    self.inputConfigFile = inputConfigFile

    if(spikeDataFileName is None):
      spikeDataFileName = "save/input-spikes-" + str(self.networkSlurmID) + ".hdf5"
    self.spikeDataFileName = spikeDataFileName

    self.time = time # How long time to generate inputs for

    self.neuronCache = dict([])
    
    # Read in the input configuration information from JSON file
    self.readInputConfigFile()

    # Read the network position file
    self.readNeuronPositions()

    # Read the network config file
    self.readNetworkConfigFile()

    # Only the master node should start the work
    if(isMaster):
      # Initialises lbView and dView (load balance, and direct view)
      self.setupParallell()

    
      # Make the "master input" for each channel
      self.makePopulationUnitSpikeTrains()
      
      # Generate the actual input spikes, and the locations
      # stored in self.neuronInput dictionary

      self.makeNeuronInputParallell()      

      # Consolidate code, so same code runs for serial and parallel case
      #if(self.lbView is None):
      #  self.makeNeuronInput()
      #else:
      #  self.makeNeuronInputParallell()      

      # Write spikes to disk, HDF5 format
      self.writeHDF5()

      # Verify correlation --- THIS IS VERY VERY SLOW
      #self.verifyCorrelation()

      self.checkSorted()

      
        
    


    
    
    # !!! TODO

    # 1. Define what the within correlation, and between correlation should be
    #    for each neuron type. Also what input frequency should we have for each
    #    neuron. --- Does it depend on size of dendritic tree?
    #    Store the info in an internal dict.
    
    # 2. Read the position file, so we know what neurons are in the network

    # 3. Create the "master input" for each population unit.

    # 4. Mix the master input with random input, for each neuron, to create
    #    the appropriate correlations

    # 5. Randomize which compartments each synaptic input should be on

    # 6. Verify correlation of input

    # 7. Write to disk
    
    # If more than one worker node, then we need to split the data
    # into multiple files
    # self.nWorkers=nWorkers


  ############################################################################

  def writeHDF5(self):

    import timeit
    import h5py

    self.writeLog("Writing spikes to " + self.spikeDataFileName)
    
    startTime = timeit.default_timer()

    outFile = h5py.File(self.spikeDataFileName,'w',libver=self.h5libver)

    configData = outFile.create_dataset("config",
                                        data=json.dumps(self.inputInfo,
                                                        indent=4))

    inputGroup = outFile.create_group("input")
    
    for neuronID in self.neuronInput:

      NIDGroup = inputGroup.create_group(str(neuronID))

      neuronType = self.neuronType[neuronID]
      # nName = self.neuronName[neuronID]
      
      for inputType in self.neuronInput[neuronID]:
        if(inputType.lower() != "VirtualNeuron".lower()):
          itGroup = NIDGroup.create_group(inputType)

          neuronIn = self.neuronInput[neuronID][inputType]
          spikeMat,nSpikes = self.createSpikeMatrix(neuronIn["spikes"])

          itGroup.create_dataset("spikes",data=spikeMat)
          itGroup.create_dataset("nSpikes",data=nSpikes)
                        
          itGroup.create_dataset("sectionID", data=neuronIn["location"][1])
          itGroup.create_dataset("sectionX", data=neuronIn["location"][2])

          itGroup.create_dataset("freq",data=neuronIn["freq"])
          itGroup.create_dataset("correlation",data=neuronIn["correlation"])
          itGroup.create_dataset("jitter",data=neuronIn["jitter"])
          itGroup.create_dataset("synapseDensity",
                                 data=neuronIn["synapseDensity"])
          itGroup.create_dataset("start",data=neuronIn["start"])
          itGroup.create_dataset("end",data=neuronIn["end"])
          itGroup.create_dataset("conductance",data=neuronIn["conductance"])
        
          populationUnitID = neuronIn["populationUnitID"]
          itGroup.create_dataset("populationUnitID",data=populationUnitID)

          chanSpikes = self.populationUnitSpikes[neuronType][inputType][populationUnitID]
          itGroup.create_dataset("populationUnitSpikes",
                                 data=chanSpikes)

          itGroup.create_dataset("generator",data=neuronIn["generator"])

          try:
            itGroup.create_dataset("modFile",data=neuronIn["modFile"])
            if(neuronIn["parameterFile"]):
              itGroup.create_dataset("parameterFile",data=neuronIn["parameterFile"])
            # We need to convert this to string to be able to save it
            itGroup.create_dataset("parameterList",
                                   data=json.dumps(neuronIn["parameterList"]))
            itGroup.create_dataset("parameterID",data=neuronIn["parameterID"])
          except:
            import traceback
            tstr = traceback.format_exc()
            self.writeLog(tstr)
       
            import pdb
            pdb.set_trace()
 
            
            
        else:          
          
          # Input is activity of a virtual neuron
          aGroup = NIDGroup.create_group("activity")
          spikes = self.neuronInput[neuronID][inputType]["spikes"]
            
          aGroup.create_dataset("spikes",data=spikes)
          generator = self.neuronInput[neuronID][inputType]["generator"]
          aGroup.create_dataset("generator", data=generator)

          
          
    outFile.close()
    

  ############################################################################

  def createSpikeMatrix(self,spikes):

    if(len(spikes) == 0):
      return np.zeros((0,0)), 0
    
    nInputTrains = len(spikes)
    nSpikes = np.array([len(x) for x in spikes])
    maxLen = max(nSpikes)
    
    spikeMat = -1*np.ones((nInputTrains,maxLen))
    for idx,st in enumerate(spikes):
      n = st.shape[0]
      spikeMat[idx,:n] = st

    return spikeMat, nSpikes
      
  ############################################################################

  # Reads from self.inputConfigFile
  
  def readInputConfigFile(self):

    self.writeLog("Loading input configuration from " + str(self.inputConfigFile))
    
    with open(self.inputConfigFile,'rt') as f:
      self.inputInfo = json.load(f)

    for neuronType in self.inputInfo:
      for inputType in self.inputInfo[neuronType]:
        if("parameterFile" in self.inputInfo[neuronType][inputType]):
          parFile = self.inputInfo[neuronType][inputType]["parameterFile"]

          # Allow user to use $DATA to refer to snudda data directory
          parFile = parFile.replace("$DATA",
                                    os.path.dirname(__file__) + "/data")
          
          parDataDict = json.load(open(parFile,'r'))
          
          # Read in parameters into a list
          parData = []
          for pd in parDataDict:
            parData.append(parDataDict[pd])
        else:
          parData = None

        self.inputInfo[neuronType][inputType]["parameterList"] = parData
         
    
  ############################################################################

  # Each synaptic input will contain a fraction of population unit spikes, which are
  # taken from a stream of spikes unique to that particular population unit
  # This function generates these correlated spikes
  
  def makePopulationUnitSpikeTrains(self,nPopulationUnits=None,timeRange=None):

    self.writeLog("Running makePopulationUnitSpikeTrains")
    
    if(nPopulationUnits is None):
      nPopulationUnits = self.nPopulationUnits

    if(timeRange is None):
      timeRange = (0,self.time)
      
    self.populationUnitSpikes = dict([])

    for cellType in self.inputInfo:
      
      self.populationUnitSpikes[cellType] = dict([])
      
      for inputType in self.inputInfo[cellType]:

        if(self.inputInfo[cellType][inputType]["generator"] == "poisson"):
        
          freq = self.inputInfo[cellType][inputType]["frequency"]
          self.populationUnitSpikes[cellType][inputType] = dict([])

          if("populationUnitID" in self.inputInfo[cellType][inputType]):
            popUnitList = \
              self.inputInfo[cellType][inputType]["populationUnitID"]

            if(type(popUnitList) != list):
              popUnitList = [popUnitList]
          else:
            popUnitList = range(0,self.nPopulationUnits)
              
          for idxPopUnit in popUnitList:
            self.populationUnitSpikes[cellType][inputType][idxPopUnit] = \
              self.generateSpikes(freq=freq,timeRange=timeRange)

    return self.populationUnitSpikes

  ############################################################################

  
  def makeNeuronInputParallell(self):

    self.writeLog("Running makeNeuronInputParallell")
    
    self.neuronInput = dict([])

    neuronIDList = []
    inputTypeList = []
    freqList = []
    startList = []
    endList = []
    synapseDensityList = []
    nInputsList = []
    PkeepList = []
    populationUnitSpikesList = []
    jitterDtList = []
    locationList = []
    populationUnitIDList = []
    conductanceList = []
    correlationList = []

    modFileList = []
    parameterFileList = []
    parameterListList = []
    
    for (neuronID,neuronType,populationUnitID) \
        in zip(self.neuronID, self.neuronType,self.populationUnitID):
      
      self.neuronInput[neuronID] = dict([])

      if(neuronType not in self.inputInfo):
        self.writeLog("!!! Warning, synaptic input to " + str(neuronType) \
                      + " missing in " + str(self.inputConfigFile) )
        continue
      
      for inputType in self.inputInfo[neuronType]:
        
        inputInf = self.inputInfo[neuronType][inputType]

        if("populationUnitID" in inputInf):
          popUnitID = inputInf["populationUnitID"]
          
          if(type(popUnitID) == list \
             and populationUnitID not in popUnitID):
            # We have a list of functional channels, but this neuron
            # does not belong to a functional channel in that list
            continue
          elif(populationUnitID != popUnitID):
            # We have a single functional channel, but this neuron is not
            # in that functional channel
            continue

        self.neuronInput[neuronID][inputType] = dict([])
          
        if(inputInf["generator"] == "poisson"):
          neuronIDList.append(neuronID)
          inputTypeList.append(inputType)
          freqList.append(inputInf["frequency"])
          PkeepList.append(np.sqrt(inputInf["populationUnitCorrelation"]))
          jitterDtList.append(inputInf["jitter"])

          if("start" in inputInf):
            startList.append(inputInf["start"])
          else:
            startList.append(0.0) # Default start at beginning

          if("end" in inputInf):
            endList.append(inputInf["end"])
          else:
            endList.append(self.time)

          if(inputType.lower() == "VirtualNeuron".lower()):
            # Virtual neurons spikes specify their activity, location and conductance not used
            cond = None
            nInp = 1
            
            modFile = None
            parameterFile = None
            parameterList = None
          else:
            assert "location" not in inputInf, \
              "Location in input config has been replaced with synapseDensity"
            cond = inputInf["conductance"]

            if("nInputs" in inputInf):
              nInp = inputInf["nInputs"]
            else:
              nInp = None

            modFile = inputInf["modFile"]
            if("parameterFile" in inputInf):
              parameterFile = inputInf["parameterFile"]
            else:
              parameterFile = None

            if("parameterList" in inputInf):
              parameterList = inputInf["parameterList"]
            else:
              parameterList = None

          if("synapseDensity" in inputInf):
            synapseDensity = inputInf["synapseDensity"]
          else:
            synapseDensity = "1"
            
          synapseDensityList.append(synapseDensity)
          nInputsList.append(nInp)
            
          populationUnitIDList.append(populationUnitID)
          conductanceList.append(cond)
          correlationList.append(inputInf["populationUnitCorrelation"])

          cSpikes = self.populationUnitSpikes[neuronType][inputType][populationUnitID]
          populationUnitSpikesList.append(cSpikes)

          modFileList.append(modFile)
          parameterFileList.append(parameterFile)
          parameterListList.append(parameterList)
          
        elif(inputInf["generator"] == "csv"):
          csvFile = inputInf["csvFile"] % neuronID
          
          self.neuronInput[neuronID][inputType]["spikes"] \
            = np.genfromtxt(csvFile, delimiter=',')
          self.neuronInput[neuronID][inputType]["generator"] = "csv"

        else:
          self.writeLog("Unknown input generator: " + inputInf["generator"]\
                        + " for " + str(neuronID))

    # The old code had so that all neurons within a population unit shared the same
    # mother process, which caused them all to activate at the same time
    # with high probability. By setting channelSpikeList to None we disable it
    self.writeLog("Clearing populationUnitSpikesList, thus all neurons will have their own mother process for each input")
    populationUnitSpikesList = [None for x in populationUnitSpikesList]
    amr = None
    
    #Lets try and swap self.lbView for self.dView
    if(self.dView is not None):
      
      #self.writeLog("Sending jobs to workers, using lbView")
      self.writeLog("Sending jobs to workers, using dView")

      # Changed the logic, the old input helper needed a global
      # variable to be visible, but it was not always so in its scope

      inputList = list(zip(neuronIDList,
                           inputTypeList,
                           freqList,
                           startList,
                           endList,
                           synapseDensityList,
                           nInputsList,
                           PkeepList,
                           populationUnitSpikesList,
                           jitterDtList,
                           populationUnitIDList,
                           conductanceList,
                           correlationList,
                           modFileList,
                           parameterFileList,
                           parameterListList))

      self.dView.scatter("inputList",inputList,block=True)
      cmdStr = "inpt = list(map(nl.makeInputHelperParallel,inputList))"
      self.dView.execute(cmdStr,block=True)

      inpt = self.dView["inpt"]
      amr = list(itertools.chain.from_iterable(inpt))

    else:
      # If no lbView then we run it in serial
      self.writeLog("Running input generation in serial")
      amr = map(self.makeInputHelperSerial,
                neuronIDList,
                inputTypeList,
                freqList,
                startList,
                endList,
                synapseDensityList,
                nInputsList,
                PkeepList,
                populationUnitSpikesList,
                jitterDtList,
                populationUnitIDList,
                conductanceList,
                correlationList,
                modFileList,
                parameterFileList,
                parameterListList)
          
    # Gather the spikes that were generated in parallell
    for neuronID, inputType, spikes, loc, synapseDensity, frq, \
        jdt, pUID,cond,corr,timeRange, \
        modFile,paramFile,paramList,paramID in amr:
      self.writeLog("Gathering " + str(neuronID) + " - " + str(inputType) )
      self.neuronInput[neuronID][inputType]["spikes"] = spikes
      
      if(inputType.lower() != "VirtualNeuron".lower()):
        # Virtual neurons have no location of their input, as the "input"
        # specifies the spike times of the virtual neuron itself
        self.neuronInput[neuronID][inputType]["location"] = loc
        self.neuronInput[neuronID][inputType]["synapseDensity"] \
          = synapseDensity
        self.neuronInput[neuronID][inputType]["conductance"] = cond
        
      self.neuronInput[neuronID][inputType]["freq"] = frq
      self.neuronInput[neuronID][inputType]["correlation"] = corr
      self.neuronInput[neuronID][inputType]["jitter"] = jdt
      self.neuronInput[neuronID][inputType]["start"] = timeRange[0]
      self.neuronInput[neuronID][inputType]["end"] = timeRange[1]
      self.neuronInput[neuronID][inputType]["populationUnitID"] = pUID

      assert pUID == self.populationUnitID[neuronID], \
        "Internal error: Neuron should belong to the functional channel "\
        + "that input is generated for" 
      
      self.neuronInput[neuronID][inputType]["generator"] = "poisson"
      self.neuronInput[neuronID][inputType]["modFile"] = modFile
      self.neuronInput[neuronID][inputType]["parameterFile"] = paramFile
      self.neuronInput[neuronID][inputType]["parameterList"] = paramList
      self.neuronInput[neuronID][inputType]["parameterID"] = paramID      

    return self.neuronInput
     
  
  ############################################################################

  # This generates poisson spikes with frequency freq, for a given time range
  
  def generateSpikes(self,freq,timeRange):

    # https://stackoverflow.com/questions/5148635/how-to-simulate-poisson-arrival
    start = timeRange[0]
    end = timeRange[1]
    duration = end-start
    
    tDiff = -np.log(1.0 - np.random.random(int(np.ceil(max(1,freq*duration)))))/freq

    tSpikes = []
    tSpikes.append(start+np.cumsum(tDiff))

    # Is last spike after end of duration
    while(tSpikes[-1][-1] <= end):
      tDiff = -np.log(1.0 - np.random.random(int(np.ceil(freq*duration*0.1))))/freq
      tSpikes.append(tSpikes[-1][-1] + np.cumsum(tDiff))

    # Prune away any spikes after end 
    if(len(tSpikes[-1]) > 0):
      tSpikes[-1] = tSpikes[-1][tSpikes[-1] <= end]

    # Return spike times
    return np.concatenate(tSpikes)


  ############################################################################

  # This takes a list of spike trains and returns a single spike train
  # including all spikes
  
  def mixSpikes(self,spikes):

    return np.sort(np.concatenate(spikes))
    
  ############################################################################

  def cullSpikes(self,spikes,Pkeep):

    return spikes[np.random.random(spikes.shape) < Pkeep]


  ############################################################################

  # timeRange --- (start,end time) of spike train
  # freq -- frequency of spike train
  # nSpikeTrains -- number of spike trains to generate
  # Pkeep -- fraction of channel spikes to include in spike train
  # retpopUnitSpikes -- if true, returns tuple with second item population unit spikes
  #                  if false, just return spikes
  # populationUnitSpikes --- if None, new population unit spikes will be generated
  #                  (population unit Spikes are the spikes shared between correlated
  #                   spike trains)
  
  def makeCorrelatedSpikes(self,freq,timeRange,nSpikeTrains,Pkeep, \
                           populationUnitSpikes=None, \
                           retPopUnitSpikes=False,jitterDt=None):

    assert(Pkeep >= 0 and Pkeep <= 1)

    if(populationUnitSpikes is None):
      populationUnitSpikes = self.generateSpikes(freq,timeRange)

    uniqueFreq = freq * (1-Pkeep)
    spikeTrains = []

    for i in range(0,nSpikeTrains):
      tUnique = self.generateSpikes(uniqueFreq,timeRange)
      tPopulationUnit = self.cullSpikes(populationUnitSpikes,Pkeep)
      
      spikeTrains.append(self.mixSpikes([tUnique,tPopulationUnit]))


    #if(False):
    #self.verifyCorrelation(spikeTrains=spikeTrains) # THIS STEP IS VERY VERY SLOW
      
    if(jitterDt is not None):
      spikeTrains = self.jitterSpikes(spikeTrains,jitterDt,timeRange=timeRange)
      
    if(retPopUnitSpikes):
      return (spikeTrains,populationUnitSpikes)
    else:
      return spikeTrains

  ############################################################################

  def makeUncorrelatedSpikes(self,freq,start,end,nSpikeTrains):

    spikeTrains = []

    for i in range(0,nSpikeTrains):
      spikeTrains.append(self.generateSpikes(freq,start,end))
    
    return spikeTrains


  ############################################################################

  # If a timeRange (start,endtime) is given then all spike times will
  # be modulo duration, so if we jitter and they go to before start time,
  # they wrap around and appear at end of the timeline
  
  def jitterSpikes(self,spikeTrains,dt,timeRange=None):

    jitteredSpikes = []
    
    for i in range(0,len(spikeTrains)):
      spikes = spikeTrains[i] + np.random.normal(0,dt,spikeTrains[i].shape)
      
      if(timeRange is not None):
        start = timeRange[0]
        end = timeRange[1]
        spikes = np.mod(spikes-start,end-start) + start
        
      s = np.sort(spikes)
      # Remove any spikes that happened to go negative
      s = s[np.where(s >= 0)]
      jitteredSpikes.append(s)

    return jitteredSpikes
  
  ############################################################################

  # Plot spikes as a raster plot, for debugging and visualisation purposes
  
  def rasterPlot(self,spikeTimes, \
                 markSpikes=None,markIdx=None, \
                 title=None,figFile=None,fig=None):

    if(fig is None):
      fig = plt.figure()
    # ax = plt.gca()

    for i, spikes in enumerate(spikeTimes):
      plt.vlines(spikes,i+1.5,i+0.5,color="black")

    plt.ylim(0.5,len(spikeTimes)+0.5)
      
    if(markSpikes is not None and markIdx is not None):
      for i, spikes in zip(markIdx,markSpikes):
        plt.vlines(spikes,i+1.5,i+0.5,color="red")

      plt.ylim(min(0.5,min(markIdx)-0.5), \
               max(max(markIdx)+0.5,len(spikeTimes))+0.5)
      
    plt.xlabel("Time")
    plt.ylabel("Inputs")

    plt.ion()
    plt.show()
    
    if(title is not None):
      plt.title(title)

    fig.show()

    if(figFile is not None):
      plt.savefig(figFile)

    return fig
    
  ############################################################################

  def readNeuronPositions(self):

    self.writeLog("Reading neuron postions")
    
    posInfo = SnuddaLoad(self.positionFile).data
    self.networkInfo = posInfo
    self.neuronInfo = posInfo["neurons"]
    
#    import pdb
#    pdb.set_trace()
    
    # Make sure the position file matches the network config file
    assert(posInfo["configFile"] == self.networkConfigFile)

    self.nPopulationUnits = posInfo["nPopulationUnits"]
    self.populationUnitID = posInfo["populationUnit"]

    self.neuronName = [n["name"] for n in self.neuronInfo]
    
    self.neuronID = [n["neuronID"] for n in self.neuronInfo]    
    self.neuronType = [n["type"] for n in self.neuronInfo]
    # self.nInputs =  [n["nInputs"] for n in self.neuronInfo]

  ############################################################################
    
  def readNetworkConfigFile(self):

    self.writeLog("Reading config file " + str(self.networkConfigFile))
    
    import json
    
    with open(self.networkConfigFile,'r') as f:
      self.networkConfig = json.load(f)

      
  ############################################################################

  def verifyCorrelation(self,spikeTrains,expectedCorr=None,dt=0):

    # THIS FUNCTION IS VERY VERY SLOW
    
    corrVec = []

    for si,s in enumerate(spikeTrains):
      for s2i,s2 in enumerate(spikeTrains):      
        if(si == s2i):
          # No self comparison
          continue

        corrVec.append(self.estimateCorrelation(s,s2,dt=dt))
        
    # print("corr = " + str(corrVec))
    self.writeLog("meanCorr = " + str(np.mean(corrVec)))
    
  ############################################################################

  def estimateCorrelation(self,spikesA,spikesB,dt=0):

    nSpikesA = len(spikesA)
    corrSpikes = 0

    for t in spikesA:
      if(np.min(abs(spikesB-t)) <= dt):
        corrSpikes += 1
        
    return (corrSpikes / float(nSpikesA))
    
  ############################################################################

  # inputDensity = f(d) where d is micrometers from soma,
  #                unit of f is synapses/micrometer

  # !!! Returns input locations only on dendrites, not on soma
  
  def dendriteInputLocations(self,
                             neuronID,
                             synapseDensity="1",
                             nSpikeTrains=None):

    neuronName = self.neuronName[neuronID]
    swcFile = self.networkConfig["Neurons"][neuronName]["morphology"]

    if(swcFile in self.neuronCache):
      morphology = self.neuronCache[swcFile]
    else:     
      morphology = NeuronMorphology(name=neuronID,
                                    swc_filename=swcFile,
                                    axonStumpIDFlag=self.axonStumpIDFlag)
      self.neuronCache[swcFile] = morphology

    return morphology.dendriteInputLocations(synapseDensity=synapseDensity,
                                             nLocations=nSpikeTrains)
      
  ############################################################################
  
  # Returns random dendrite compartments
  # Pdist = function of d (micrometers), e.g. "1", "2*d", ... etc
  
  def randomCompartment(self,neuronID,nLocations,Pdist="1"):
    
    neuronName = self.neuronName[neuronID]
    swcFile = self.networkConfig[neuronName]["morphology"]
    
    if(swcFile in self.neuronCache):
      morphology = self.neuronCache[swcFile]
    else:
      # !!! Should I use NeuronMorphology, or maybe the loadSWC in ConvertNetwork
      #                      -- maybe move loadSWC out to a separate file
      
      morphology = NeuronMorphology(name=neuronID,
                                    swc_filename=swcFile,
                                    axonStumpIDFlag=self.axonStumpIDFlag)
      self.neuronCache[swcFile] = morphology

    # morphology.dend -- 0-2: x,y,z 3: r, 4: dist to soma
    d = morphology.dend[:,4]
    Px = eval(Pdist)

    if(type(Px) in (int, float)):
      # If Px is a constant, we need to set it for all points
      Px *= np.ones(d.shape)

    Pcomp = (Px[morphology.dendLinks[:,0]] + Px[morphology.dendLinks[:,1]])/2
    compLen = morphology.compartmentLength(compType="dend")
    
      
    # Multiply by length, so longer compartments are proportionally more likely to be
    # connected to
    Pcomp = np.multiply(Pcomp,compLen)

    # Randomize locations. Here we then sort the locations, this is to make
    # it quicker to locate where they are (one pass in the Pxcumsum array)
    Pxcumsum = np.cumsum(Pcomp)
    x = np.random.uniform(low=0,high=Pxcumsum[-1],size=nLocations)
    xsort = np.sort(x)

    nComps = len(d)
    lastValidIdx = nComps - 1 
    compIdx = 0
    
    secID = np.zeros((nLocations,),dtype=int)
    secX = np.zeros((nLocations,))
    compCoords = np.zeros((nLocations,3))
    
    
    for ix,xval in enumerate(xsort):
      while(xval > Pxcumsum[compIdx] and compIdx < lastValidIdx):
        compIdx += 1
        
      secID[ix] = morphology.dendSecID[compIdx]
      secX[ix] = np.random.rand()*(morphology.dendSecX[compIdx,1] \
                                   - morphology.dendSecX[compIdx,0]) \
                 + morphology.dendSecX[compIdx,0]

      compCoords[ix,:3] = (morphology.dend[morphology.dendLinks[compIdx,0],:3] \
                        + morphology.dend[morphology.dendLinks[compIdx,1],:3])/2
      
        
    # The probability of picking a compartment is dependent on a distance
    # dependent probability function, and the length of the compartment

    # Plot synapses and segments, as a verification
    if(False):
      import matplotlib.pyplot as plt
      plt.hist(d[compID],label="synapses")
      plt.hist(d,label="segments")
      plt.legend(loc='upper right')
      plt.show()
    
    return(compCoords,secID,secX)

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
    self.dView.execute("nl.setSeed(workerSeed[0])",block=True)

    self.writeLog("New worker seeds: " + str(workerSeeds))
  
  ############################################################################

  def setupParallell(self):
    
    import os
    SlurmJobID = os.getenv("SLURM_JOBID")

    if(SlurmJobID is None):
      self.SlurmID = 0
    else:
      self.SlurmID = int(SlurmJobID)

    self.writeLog("IPYTHON_PROFILE = " + str(os.getenv('IPYTHON_PROFILE')))
    
    if(os.getenv('IPYTHON_PROFILE') is not None):
      from ipyparallel import Client
      self.rc = Client(profile=os.getenv('IPYTHON_PROFILE'))

      # http://davidmasad.com/blog/simulation-with-ipyparallel/
      # http://people.duke.edu/~ccc14/sta-663-2016/19C_IPyParallel.html
      self.writeLog("Client IDs: " + str(self.rc.ids))
      self.dView = self.rc[:] # Direct view into clients
      self.lbView = self.rc.load_balanced_view()

      if(self.logFile is not None):
        logFileName = self.logFile.name
        engineLogFile = [logFileName + "-" \
                         + str(x) for x in range(0,len(self.dView))]
      else:
        engineLogFile = [None for x in range(0,len(self.dView))]
    else:
      self.writeLog("No IPYTHON_PROFILE enviroment variable set, running in serial")
      self.dView = None
      self.lbView = None
      return

    with self.dView.sync_imports():
      from snudda.input import SnuddaInput

    self.dView.push({ "inputConfigFile" : self.inputConfigFile,
                      "networkConfigFile" : self.networkConfigFile,
                      "positionFile" : self.positionFile,
                      "spikeDataFileName" : self.spikeDataFileName,
                      "isMaster" : False,
                      "time" : self.time})

    self.writeLog("Scattering engineLogFile = " + str(engineLogFile))
      
    self.dView.scatter('logFileName',engineLogFile,block=True)
   
    
    self.writeLog("nl = SnuddaInput(inputConfigFile='" + self.inputConfigFile \
                  + "',networkConfigFile='" + self.networkConfigFile \
                  + "',positionFile='" + self.positionFile\
                  + "',spikeDataFileName='" + self.spikeDataFileName \
                  + "',isMaster=False " \
                  + ",time=" +str(self.time) + ",logFile=logFileName[0])")
    
    cmdStr = 'global nl; nl = SnuddaInput(inputConfigFile=inputConfigFile,networkConfigFile=networkConfigFile,positionFile=positionFile,spikeDataFileName=spikeDataFileName,isMaster=isMaster,time=time,logFile=logFileName[0])'
    
    self.dView.execute(cmdStr,block=True)

    self.newWorkerSeeds(self.dView)

    self.writeLog("Workers set up")
    
  ############################################################################
    
  # Function for debugging

  def dumpToRandomFile(self,filePrefix,dataToDump):
    import uuid
    tmp = open("save/" + filePrefix + "-file-" + str(uuid.uuid4()),'w')
    tmp.write(str(dataToDump))
    tmp.close()

  ############################################################################

  def checkSorted(self):

    # Just a double check that the spikes are not jumbled
    
    for neuronID in self.neuronInput:
      for inputType in self.neuronInput[neuronID]:
        if(inputType == "VirtualNeuron"):
          s = self.neuronInput[neuronID][inputType]["spikes"]
          assert (np.diff(s) >= 0).all(), \
            str(neuronID) + " " + inputType + ": Spikes must be in order"
        else:
          for spikes in self.neuronInput[neuronID][inputType]["spikes"]:
            assert len(spikes) == 0 or spikes[0] >= 0
            assert (np.diff(spikes) >= 0).all(), \
              str(neuronID) + " " + inputType + ": Spikes must be in order"

  ############################################################################
          
  def plotSpikes(self,neuronID=None):

    self.writeLog("Plotting spikes for neuronID: " + str(neuronID))
    
    if(neuronID is None):
      neuronID = self.neuronInput

    spikeTimes = []
    
    for nID in neuronID:
      for inputType in self.neuronInput[nID]:
        for spikes in self.neuronInput[nID][inputType]["spikes"]:
          spikeTimes.append(spikes)

    self.rasterPlot(spikeTimes)
          
  ############################################################################

  def readHDF5info(self,hdf5File):
    self.writeLog("Loading HDF5-file: " + hdf5File)

    try:
      with h5py.File(hdf5File,'r') as f:
        self.networkConfigFile = f["meta"]["configFile"][()]
        self.positionFile = f["meta"]["positionFile"][()]
        self.networkSlurmID = int(f["meta/SlurmID"][()])
        
        self.axonStumpIDFlag = f["meta/axonStumpIDFlag"][()]
    except Exception as e:
      self.writeLog("Error in readHDF5info: " + str(e))
      self.writeLog("Opening: " + hdf5File)

      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
       
      import pdb
      pdb.set_trace()
 
          
  ############################################################################


  def findLatestFile(self):

    files = glob('save/network-connect-voxel-pruned-synapse-file-*.hdf5')

    modTime = [os.path.getmtime(f) for f in files]
    idx = np.argsort(modTime)

    self.writeLog("Using the newest file: " + files[idx[-1]])
    
    return files[idx[-1]]
  
  
  ############################################################################

  def makeInputHelperParallel(self,args):

    try:
    
      neuronID,inputType,freq,start,end,synapseDensity,nSpikeTrains,Pkeep,populationUnitSpikes,jitterDt,populationUnitID,conductance,correlation,modFile,parameterFile,parameterList = args

      return self.makeInputHelperSerial(neuronID = neuronID,
                                        inputType = inputType,
                                        freq = freq,
                                        start = start,
                                        end = end,
                                        synapseDensity = synapseDensity,
                                        nSpikeTrains = nSpikeTrains,
                                        Pkeep = Pkeep,
                                        populationUnitSpikes = populationUnitSpikes,
                                        jitterDt = jitterDt,
                                        populationUnitID = populationUnitID,
                                        conductance = conductance,
                                        correlation = correlation,
                                        modFile = modFile,
                                        parameterFile = parameterFile,
                                        parameterList = parameterList)

    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()
      

  ############################################################################

  # Normally specify synapseDensity which then sets number of inputs
  # ie leave nSpikeTrains as None. If nSpikeTrains is set, that will then
  # scale synapseDensity to get the requested number of inputs (approximately)

  # For virtual neurons nSpikeTrains must be set, as it defines their activity
  
  def makeInputHelperSerial(self,
                            neuronID,
                            inputType,
                            freq,
                            start,
                            end,
                            synapseDensity,
                            nSpikeTrains,
                            Pkeep,
                            populationUnitSpikes,
                            jitterDt,
                            populationUnitID,
                            conductance,
                            correlation,
                            modFile,
                            parameterFile,
                            parameterList):
                            
    # First, find out how many inputs and where, based on morphology and
    # synapse density

    try:

      timeRange = (start,end)
      
      if(inputType.lower() == "VirtualNeuron".lower()):
        # This specifies activity of a virtual neuron
        loc = None
        conductance = None

        assert nSpikeTrains is None or nSpikeTrains == 1, \
          "Virtual neuron " + self.neuronName[neuronID] \
          + " should have only one spike train, fix nSpikeTrains in config"
        
        spikes = self.makeCorrelatedSpikes(freq=freq,
                                           timeRange=timeRange,
                                           nSpikeTrains=1,
                                           Pkeep=Pkeep,
                                           populationUnitSpikes=populationUnitSpikes,
                                           jitterDt=jitterDt)
        nInputs = 1
      else:
        
        # x,y,z, secID, secX    
        inputLoc = self.dendriteInputLocations(neuronID=neuronID,
                                               synapseDensity=synapseDensity,
                                               nSpikeTrains=nSpikeTrains)

        nInputs = inputLoc[0].shape[0]
        print("Generating " + str(nInputs) + " inputs for " \
              + self.neuronName[neuronID])
        
        # OBS, nInputs might differ slightly from nSpikeTrains if that is given
        spikes = self.makeCorrelatedSpikes(freq=freq,
                                           timeRange=timeRange,
                                           nSpikeTrains=nInputs,
                                           Pkeep=Pkeep,
                                           populationUnitSpikes=populationUnitSpikes,
                                           jitterDt=jitterDt)        
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()

      
    # We need to pick which parameter set to use for the input also
    parameterID = np.random.randint(1e6,size=nInputs)
      
    # We need to keep track of the neuronID, since it will all be jumbled
    # when doing asynchronous prallellisation
    return (neuronID, inputType, spikes, inputLoc, synapseDensity,freq,
            jitterDt,populationUnitID,conductance,correlation,
            timeRange,
            modFile,parameterFile,parameterList,parameterID)


      
  ############################################################################
  
  def makeInputHelperSerialOLD(self,
                            neuronID,
                            inputType,
                            freq,
                            start,
                            end,
                            nSpikeTrains,
                            Pkeep,
                            populationUnitSpikes,
                            jitterDt,
                            location,
                            populationUnitID,
                            conductance,
                            correlation,
                            modFile,
                            parameterFile,
                            parameterList,
                            parameterID):
    try:

      assert False, "Depricated use makeInputHelperSerial"
    
      timeRange = (start,end)
      
      spikes = self.makeCorrelatedSpikes(freq=freq,
                                         timeRange=timeRange,
                                         nSpikeTrains=nSpikeTrains,
                                         Pkeep=Pkeep,
                                         populationUnitSpikes=populationUnitSpikes,
                                         jitterDt=jitterDt)

      if(inputType.lower() == "VirtualNeuron".lower()):
        loc = None
        conductance = None
      else:
        loc = self.randomCompartment(neuronID,nSpikeTrains,location)
      
    except Exception as e:
      import uuid
      import traceback
      tstr = traceback.format_exc()
      tmp = open("save/tmp-log-file-" + str(uuid.uuid4()),'w')
      tmp.write("Exception: " + str(e))
      tmp.write("Trace:" + tstr)
      tmp.close()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()
    
    # We need to keep track of the neuronID, since it will all be jumbled
    # when doing asynchronous prallellisation
    return (neuronID, inputType, spikes, loc, freq,
            jitterDt,populationUnitID,conductance,correlation,location,
            timeRange,
            modFile,parameterFile,parameterList,parameterID)

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

  
if __name__ == "__main__":

  print("Please do not call this file directly, use snudda.py")
  exit(-1)
