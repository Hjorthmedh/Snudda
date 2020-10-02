#
# This code reads the network created by Network_connect.py and set it
# up in memory
#
# mpiexec -n 4 python snudda_simulate.py
#
#
# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union's Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).

#
############################################################################

# Plot all sections
# [neuron.h.psection(x) for x in neuron.h.allsec()]



from mpi4py import MPI # This must be imported before neuron, to run parallel
from neuron import h #, gui
import neuron
import h5py
import json
import timeit

import bluepyopt.ephys as ephys
from snudda.Neuron_model_extended import NeuronModel
import pickle
# from Network_place_neurons import NetworkPlaceNeurons
import numpy as np
from snudda.NrnSimulatorParallel import NrnSimulatorParallel

from glob import glob
import re
import os

# !!! Need to gracefully handle the situation where there are more workers than
# number of neurons, currently we get problem when adding the voltage saving

# !!! Have added code for dopamine modulation of neuron intrinsic channels
#     need to add dopamine modulation for the synaptic channels also !!!



##############################################################################

# If simulationConfig is set, those values override other values

class SnuddaSimulate(object):

  def __init__(self, networkFile, inputFile=None,
               verbose=True, logFile=None, \
               disableGapJunctions=True,
               simulationConfig=None):

    self.verbose = verbose
    self.logFile = logFile

    self.networkFile = networkFile
    self.inputFile = inputFile

    if(simulationConfig):
      simInfo = json.load(simulationConfig)

      if("networkFile" in simInfo):
        self.networkFile = networkFile

      if("inputFile" in simInfo):
        self.inputFile = inputFile

      if(logFile in simInfo):
        self.logFile = open(logFile,"w")


    if(type(self.logFile) == str):
      self.logFile = open(self.logFile,"w")

    self.writeLog("Using networkFile: " + str(networkFile))
    self.writeLog("Using inputFile: " + str(inputFile))

    if(self.logFile is not None):
      self.writeLog("Using logFile: " + str(self.logFile.name))



    # !!! What value to use for synaptic weight and synapse delay?
    # !!! different for AMPA and GABA?
    self.synapseWeight = 10.0 # microsiemens
    self.synapseDelay = 1      # ms
    self.spikeThreshold = -20
    self.axonSpeed = 0.8 # Tepper and Lee 2007, Wilson 1986, Wilson 1990
                         # refs taken from Damodaran et al 2013

    self.disableGapJunctions = disableGapJunctions

    self.synapseTypeLookup = { 1 : "GABA", 2: "AMPA_NMDA", 3: "GapJunction" }

    self.neurons = {}
    self.sim = None
    self.neuronNodes = []

    self.virtualNeurons = {}

    self.netConList = [] # Avoid premature garbage collection
    self.synapseList = []
    self.iStim = []
    self.vClampList = []
    self.gapJunctionList = []
    self.externalStim = dict([])
    self.tSave = []
    self.vSave = []
    self.vKey = []
    self.iSave = []
    self.iKey = []

    self.inputData = None

    self.gapJunctionNextGid = 0 # Are these gids separate from cell gids?

    self.tSpikes = h.Vector()
    self.idSpikes = h.Vector()


    # Make sure the output dir exists, so we dont fail at end because we
    # cant write file
    self.createDir("save/traces")

    self.pc = h.ParallelContext()

    # self.writeLog("I am node " + str(int(self.pc.id())))


    # We need to initialise random streams, see Lytton el at 2016 (p2072)

    self.loadNetworkInfo(networkFile)

    self.checkMemoryStatus()
    self.distributeNeurons()
    self.setupNeurons()
    self.checkMemoryStatus()
    self.pc.barrier()

#    for i in range(0,self.nNeurons):
#      print("Node : " + str(int(self.pc.id())) + " cell " + str(i) + " status " + str(self.pc.gid_exists(i)))


    self.connectNetwork()
    self.checkMemoryStatus()
    self.pc.barrier()

    # Do we need blocking call here, to make sure all neurons are setup
    # before we try and connect them



    # READ ABOUT PARALLEL NEURON
# https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/network/parcon.html#paralleltransfer

  ############################################################################

  def loadNetworkInfo(self, networkFile, config_file=None):

    self.networkFile = networkFile

    if(False):
      # We need to check if a split network file exists
      splitFile = networkFile.replace('save/','save/TEMP/').replace('.hdf5', '-%d.hdf5') % int(self.pc.id())

      import os.path
      if(os.path.isfile(splitFile)):
        networkFile = splitFile
      else:
        self.writeLog("Unable to find " + splitFile + " using " \
                      + networkFile)

    self.writeLog("Worker " + str(int(self.pc.id())) \
                  + ": Loading network from " + networkFile)

    from snudda.load import SnuddaLoad
    self.snuddaLoader = SnuddaLoad(networkFile)
    self.network_info = self.snuddaLoader.data

    self.synapses = self.network_info["synapses"]
    self.gapJunctions = self.network_info["gapJunctions"]

    # We are only passed information about neurons on our node if
    # SplitConnectionFile was run, so need to use nNeurons to know
    # how many neurons in total
    self.nNeurons = self.network_info["nNeurons"]

    if(config_file is None):
      config_file = self.getPath(self.network_info["configFile"])

    self.config_file = config_file
    self.writeLog("Loading config file " + config_file)

    # Add checks to see that config file and networkFile matches

    import json
    with open(config_file,'r') as config_file:
      self.config = json.load(config_file)

    # I do not know if the gap junction GIDs are a separate entity from the
    # neuron cell GIDs, so to be on safe side, let's make sure they
    # do not overlap
    self.gapJunctionNextGid = self.nNeurons + 100000000

    # Make a bool array indicating if cells are virtual or not
    self.isVirtualNeuron = [n["virtualNeuron"] \
                            for n in self.network_info["neurons"]]




  ############################################################################

  def distributeNeurons(self):
    # This code is run on all workers, will generate different lists on each
    self.writeLog("Distributing neurons.")

    self.neuronID = range(int(self.pc.id()),self.nNeurons,int(self.pc.nhost()))

    self.neuronNodes = [x%int(self.pc.nhost()) for x in range(0,self.nNeurons)]

    if(False):
      self.writeLog("Node " + str(int(self.pc.id())) + " handling neurons: " \
                    + ' '.join(map(str, self.neuronID)))

  ############################################################################

  def destroy(self):
    for neuron in self.neurons:
      neuron.destroy(sim=self.sim)

  ############################################################################

  # This requires self.sim to be defined

  def loadSynapseParameters(self):

    # We need to load all the synapse parameters
    self.synapseParameters = dict([])

    for (preType,postType) in self.network_info["connectivityDistributions"]:

      synData = self.network_info["connectivityDistributions"][preType,postType]

      for synType in synData:

        synapseTypeID = synData[synType]["channelModelID"]
        infoDict = synData[synType]

        if(synapseTypeID == 3):
          # Gap junctions, skip parameters
          continue

        #import pdb
        #pdb.set_trace()

        if("channelParameters" in infoDict \
           and infoDict["channelParameters"] is not None):
          channelParamDict = infoDict["channelParameters"].copy()
          modFile = channelParamDict["modFile"]

          try:
            evalStr = "self.sim.neuron.h." + modFile
            channelModule = eval(evalStr)

          except:
            self.writeLog("Are your NEURON modfiles correctly compiled?")
            
            import traceback
            tstr = traceback.format_exc()
            print(tstr)
            import pdb
            pdb.set_trace()

          try:
            # These are not variables to set in the modFile
            if("modFile" in channelParamDict):
              del channelParamDict["modFile"]

            if("parameterFile" in channelParamDict):
              del channelParamDict["parameterFile"]
              
          except:
            print(str(channelParamDict))
            
            import traceback
            tstr = traceback.format_exc()
            print(tstr)
            import pdb
            pdb.set_trace()

            
        else:
          channelParamDict = dict()
          modFile = None

          assert False, "No channel module specified for " \
            + str(preType) + "->" + str(postType) + " synapses, type ID= " \
            + str(synapseTypeID)
          channelModule = None


        if("parameterFile" in infoDict \
           and infoDict["parameterFile"] is not None):
          parFile = self.getPath(infoDict["parameterFile"])
          parDataDict = json.load(open(parFile,'r'))

          # Save data as a list, we dont need the keys
          parData = []
          for pd in parDataDict:
            if("synapse" in parDataDict[pd]):

              # Add channel parameters specified in network file, however
              # any values in the synapse parameter file will overwrite them
              pDict = channelParamDict.copy()
              for x in parDataDict[pd]["synapse"]:
                pDict[x] = parDataDict[pd]["synapse"][x]

              parData.append(pDict)
            else:
              self.writeLog("WARNING: Old data format in parameter file " \
                            + str(parFile))

              pDict = channelParamDict.copy()
              for x in parDataDict[pd]:
                pDict[x] = parDataDict[pd][x]

              parData.append(pDict)
        elif(len(channelParamDict) > 0):
          parData = [channelParamDict]
        else:
          parData = None

        self.synapseParameters[synapseTypeID] = (channelModule,parData)


  ############################################################################

  def setupNeurons(self):

    self.writeLog("Setup neurons")

    # self.sim = ephys.simulators.NrnSimulator(cvode_active=False)
    #self.sim = NrnSimulatorParallel()
    self.sim = NrnSimulatorParallel(cvode_active=False)

    # We need to load all the synapse parameters
    self.loadSynapseParameters()

    # The neurons this node is responsible for is in self.neuronID
    for ID in self.neuronID:

      name = self.network_info["neurons"][ID]["name"]

      config = self.config["Neurons"][name]
      
      morph = self.getPath(config["morphology"])
      param = self.getPath(config["parameters"])
      mech = self.getPath(config["mechanisms"])

      if("modulation" in config):
        modulation = self.getPath(config["modulation"])
      else:
        modulation = None
      
      # Obs, neurons is a dictionary
      if(self.network_info["neurons"][ID]["virtualNeuron"]):

        if(self.inputData is None):
          self.writeLog("Using " + self.inputFile + " for virtual neurons")
          self.inputData = h5py.File(self.getPath(self.inputFile),'r')

          name = self.network_info["neurons"][ID]["name"]
          spikes = self.inputData["input"][ID]["activity"]["spikes"][:,0]

        # Creating NEURON VecStim and vector
        # https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3125
        vs = h.VecStim()
        v = h.Vector(spikes.size)
        v.from_python(spikes)
        vs.play(v)

        self.virtualNeurons[ID] = dict([])
        self.virtualNeurons[ID]["spikes"] = (v,vs,spikes)
        self.virtualNeurons[ID]["name"] = name

        self.pc.set_gid2node(ID, int(self.pc.id()))

        nc = h.NetCon(vs,None)
        self.pc.cell(ID,nc,1) # The 1 means broadcast spikes to other machines

      else:
        # A real neuron (not a virtual neuron that just provides input)
        parameterID = self.network_info["neurons"][ID]["parameterID"]
        modulationID = self.network_info["neurons"][ID]["modulationID"]

        self.neurons[ID] = NeuronModel(param_file=param,
                                       morph_file=morph,
                                       mech_file=mech,
                                       cell_name=name,
                                       modulation_file=modulation,
                                       parameterID=parameterID,
                                       modulationID=modulationID)

        # Register ID as belonging to this worker node
        self.pc.set_gid2node(ID, int(self.pc.id()))

        if(True or False):
          self.writeLog("Node " + str(int(self.pc.id())) + " - cell " \
                        + str(ID) + " " + name)

        # We need to instantiate the cell
        try:
          self.neurons[ID].instantiate(sim=self.sim)

          self.setRestingVoltage(ID)

        except:
          import traceback
          tstr = traceback.format_exc()
          print(tstr)
          import pdb
          pdb.set_trace()


        # !!! DIRTY FIX for
        # https://github.com/BlueBrain/BluePyOpt/blob/master/bluepyopt/ephys/morphologies.py
        # This is likely the offending line, that pushes a segment to the stack
        # --> sim.neuron.h.execute('create axon[2]', icell)

        self.writeLog("!!! Popping extra segment from neuron -- temp fix!")
        h.execute("pop_section()")

        # !!! END OF DIRTY FIX

        # !!! Connect a netcon and register it, taken from ballandstick's
        #     connect2target function
        nc = h.NetCon(self.neurons[ID].icell.axon[0](0.5)._ref_v,
                      None,
                      sec = self.neurons[ID].icell.axon[0])
        nc.threshold = 10

        self.pc.cell(ID,nc,1) # The 1 means broadcast spikes to other machines
        # self.pc.outputcell(ID) # -- not needed, cell was called with a 1
        # self.netConList.append(nc) -- Not needed according to Lytton et al 2016

        # Record all spikes
        self.pc.spike_record(ID,self.tSpikes,self.idSpikes)

  ############################################################################

  def connectNetwork(self):

    self.pc.barrier()

    # Add synapses
    self.connectNetworkSynapses()

    # Add gap junctions
    if(self.disableGapJunctions):
      self.writeLog("!!! Gap junctions disabled.")
    else:
      self.writeLog("Adding gap junctions.")
      #self.connectNetworkGapJunctions()

      self.connectNetworkGapJunctionsLOCAL()
      self.pc.setup_transfer()
    self.pc.barrier()

  ############################################################################

  def connectNetworkSynapses(self):

    self.writeLog("connectNetworkSynapses")

    # This loops through all the synapses, and connects the relevant ones
    nextRow = 0
    # nextRowSet = [ fromRow, toRow ) -- ie range(fromRow,toRow)
    nextRowSet = self.findNextSynapseGroup(nextRow)

    while(nextRowSet is not None):

      # Add the synapses to the neuron
      self.connectNeuronSynapses(startRow=nextRowSet[0],endRow=nextRowSet[1])

      # Find the next group of synapses
      nextRow = nextRowSet[1] # 2nd number was not included in range
      nextRowSet = self.findNextSynapseGroup(nextRow)

  ############################################################################

  # This function starts at nextRow, then returns all synapses onto
  # a neuron which is located on the worker

  # This works for synapses, but it will not work for gap junctions, because
  # we need to connect the gap junctions from both sides

  # --- perhaps rewrite this as an iterator

  def findNextSynapseGroup(self,nextRow=0,connectionType="synapses"):

    if(connectionType == "synapses"):
      synapses = self.synapses
    elif(connectionType == "gapjunctions"):
      synapses = self.gapJunctions
    else:
      self.writeLog("!!! findNextSynapseGroup: Unknown connectionType: " \
                    + connectionType)
      import pdb
      pdb.set_trace()

    try:
      nSynRows = synapses.shape[0]
    except:
      self.writeLog("findNextSynapseGroup: If synapses was not loaded into memory, your problem is probably that the HDF5 file that holds the synapses were closed. Sorry.")
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()


    if(nextRow >= nSynRows):
      # No more synapses to get
      return None

    # The synapse matrix is sorted on destID, ascending order
    # We also assume that self.neuronID is sorted in ascending order

    startRow = None
    notOurID = None

    while(startRow is None):

      # What is the next destination ID
      nextID = synapses[nextRow,1]

      # Is the next ID ours?
      if(nextID in self.neuronID):
        foundRow = True
        startRow = nextRow
        ourID = nextID
        continue
      else:
        notOurID = nextID

      while(nextRow < nSynRows and\
            synapses[nextRow,1] == notOurID):
        nextRow += 1

      if(nextRow >= nSynRows):
        # No more synapses to get
        return None


    # Next find the last of the rows with this ID
    endRow = startRow

    while(endRow < nSynRows \
          and synapses[endRow,1] == ourID):
      endRow += 1

    return (startRow,endRow)

  ############################################################################

  # Processing the range(startRow,endRow) (ie, to endRow-1)

  def connectNeuronSynapses(self,startRow,endRow):

    sourceIDs = self.synapses[startRow:endRow,0]
    destID = self.synapses[startRow,1]
    assert (self.synapses[startRow:endRow,1] == destID).all()

    # Double check mapping
    assert self.pc.gid2cell(destID) == self.neurons[destID].icell, \
      "GID mismatch: " + str(self.pc.gid2cell(destID)) \
      + " != " + str(self.neurons[destID].icell)

    synapseTypeID = self.synapses[startRow:endRow,6]
    axonDistance = self.synapses[startRow:endRow,7] # Obs in micrometers

    secID=self.synapses[startRow:endRow,9]
    dendSections = self.neurons[destID].mapIDtoCompartment(secID)
    secX = self.synapses[startRow:endRow,10]/1000.0 # Convert to number 0-1

    # conductances are stored in pS (because we use INTs),
    # Neuron wants it in microsiemens??!
    conductance = self.synapses[startRow:endRow,11]*1e-6
    parameterID = self.synapses[startRow:endRow,12]

    voxelCoords = self.synapses[startRow:endRow,2:5]
    self.verifySynapsePlacement(dendSections,secX,destID,voxelCoords)

    for (srcID,section,sectionX,sTypeID,axonDist,cond,pID) \
      in zip(sourceIDs,dendSections,secX,synapseTypeID,
             axonDistance,conductance,parameterID):

      try:
        # !!!
        self.addSynapse(cellIDsource=srcID,
                        dendCompartment=section,
                        sectionDist=sectionX,
                        synapseTypeID=sTypeID,
                        axonDist=axonDist,
                        conductance=cond,
                        parameterID=pID)
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)
        import pdb
        pdb.set_trace()


  ############################################################################

  # OBS!! The src and dest lists can be different length
  #
  # src are all the gap junctions where the source compartment are
  # on the local worker.
  # dest are the gap junctions where the dest compartment are on the
  # local worker
  # The same GJ might appear in both src and dest lists, but at different rows

  def findLocalGapJunctions(self):

    # If the gap junction matrix is too large to fit in memory then
    # this will need to be optimised

    self.writeLog("Finding node local gap junctions...")

    GJidxA = np.where([x in self.neuronID \
                       for x in self.gapJunctions[:,0]])[0]

    GJidxB = np.where([x in self.neuronID \
                       for x in self.gapJunctions[:,1]])[0]

    # GJIDoffset = self.network_info["GJIDoffset"]
    GJIDoffset = 100*self.nNeurons
    GJGIDsrcA = GJIDoffset + 2*GJidxA
    GJGIDsrcB = GJIDoffset + 2*GJidxB+1

    GJGIDdestA = GJIDoffset + 2*GJidxA + 1
    GJGIDdestB = GJIDoffset + 2*GJidxB + 0

    neuronIDA = self.gapJunctions[GJidxA,0]
    neuronIDB = self.gapJunctions[GJidxB,1]

    segIDA = self.gapJunctions[GJidxA,2]
    segIDB = self.gapJunctions[GJidxB,3]

    compartmentA = [self.neurons[x].mapIDtoCompartment([y])[0] \
                    for (x,y) in zip(neuronIDA,segIDA)]
    compartmentB = [self.neurons[x].mapIDtoCompartment([y])[0] \
                    for (x,y) in zip(neuronIDB,segIDB)]

    segXA = self.gapJunctions[GJidxA,4] / 1000.0
    segXB = self.gapJunctions[GJidxB,5] / 1000.0

    # Since we had ints we stored pS, but Neuron wants microsiemens
    condA = self.gapJunctions[GJidxA,10] * 1e-6
    condB = self.gapJunctions[GJidxB,10] * 1e-6

    # Merge the two lists together

    GJidx = np.concatenate([GJidxA,GJidxB])
    GJGIDsrc = np.concatenate([GJGIDsrcA,GJGIDsrcB])
    GJGIDdest = np.concatenate([GJGIDdestA,GJGIDdestB])
    neuronID = np.concatenate([neuronIDA,neuronIDB])
    segID = np.concatenate([segIDA,segIDB])
    compartment = np.concatenate([compartmentA,compartmentB])
    segX = np.concatenate([segXA,segXB])
    cond = np.concatenate([condA,condB])

    return (neuronID,compartment,segX,GJGIDsrc,GJGIDdest,cond)

  ############################################################################

  # We can only do half the setup of the gap junction if it is split between
  # two workers.

  def connectNetworkGapJunctionsLOCAL(self):

    self.writeLog("connectNetworkGapJunctionsLOCAL")

    (neuronID,compartment,segX,GJGIDsrc,GJGIDdest,cond) \
     = self.findLocalGapJunctions()

    #import pdb
    #pdb.set_trace()

    try:
      # WHY??!
      # ValueError: too many values to unpack (expected 6)

      for nID,comp,sX,GIDsrc,GIDdest,g \
          in zip(neuronID,compartment,segX,GJGIDsrc,GJGIDdest,cond):

        self.addGapJunction(section=comp,
                            sectionDist=sX,
                            GIDsourceGJ=GIDsrc,
                            GIDdestGJ=GIDdest,
                            gGapJunction=g)

    except:
      import traceback
      tstr = traceback.format_exc()
      print(tstr)
      import pdb
      pdb.set_trace()


  ############################################################################

  def connectNetworkGapJunctions(self):

    self.writeLog("connectNetworkGapJunctions")

    self.writeLog("!!! Please verify connectNeuronGapJunctions function, that currents go bidirectionally")

    # This loops through all the synapses, and connects the relevant ones
    nextRow = 0
    # nextRowSet = [ fromRow, toRow ) -- ie range(fromRow,toRow)
    nextRowSet = self.findNextSynapseGroup(nextRow,
                                           connectionType="gapjunctions")

    while(nextRowSet is not None):

      # Add the synapses to the neuron
      self.connectNeuronGapJunctions(startRow=nextRowSet[0],
                                     endRow=nextRowSet[1])

      # Find the next group of synapses
      nextRow = nextRowSet[1] # 2nd number was not included in range
      nextRowSet = self.findNextSynapseGroup(nextRow,
                                             connectionType="gapjunctions")


  ############################################################################

  # Verify that the gap junctions conducts currents in both directions

  def connectNeuronGapJunctions(self,startRow,endRow):

    sourceID = self.gapJunctions[startRow:endRow,0]
    destID = self.gapJunctions[startRow,1]

    assert (self.gapJunctions[startRow:endRow,1] == destID).all()

    # Double check mapping
    assert self.pc.gid2cell(destID) == self.neurons[destID].icell, \
      "GID mismatch: " + str(self.pc.gid2cell(destID)) \
      + " != " + str(self.neurons[destID].icell)

    sourceSecID = self.gapJunctions[startRow:endRow,2]
    destSecID = self.gapJunctions[startRow:endRow,3]

    # !!! Double check we get number between 0.0 and 1.0
    sourceSecX = self.gapJunctions[startRow:endRow,4]*1e-4
    destSecX = self.gapJunctions[startRow:endRow,5]*1e-4

    # conductances are stored in pS, Neuron wants it in microsiements??!
    # (reason for not storing SI units is that we use INTs)
    conductance = self.gapJunctions[startRow:endRow,10]*1e-6

    destLoc = self.neurons[destID].mapIDtoCompartment(destSecID)

    for (srcID,srcSecID,srcSecX,dstLoc,dstSecX,rowIdx) \
        in zip(sourceID,sourceSecID,sourceSecX,
               destLoc,destSecX,range(startRow,endRow)):

      srcLoc = self.neurons[srcID].mapIDtoCompartment([srcSecID])[0]

      # Change the src and dest ID to be based on the row idx
      # to avoid overlaps
      GJsrcID = rowIdx * 2 + 10*self.nNeurons
      GJdestID = GJsrcID + 1

      print("rowIdx = " + str(rowIdx))
      print("GJsrcID = " + str(GJsrcID))
      print("GJdestID = " + str(GJdestID))

      try:
        self.addGapJunction(section=srcLoc,
                            sectionDist=srcSecX,
                            GIDsourceGJ=GJsrcID,
                            GIDdestGJ=GJdestID)

        self.addGapJunction(section=dstLoc,
                            sectionDist=dstSecX,
                            GIDsourceGJ=GJdestID,
                            GIDdestGJ=GJsrcID)
      except:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)
        import pdb
        pdb.set_trace()



############################################################################

  def findGapJunctionCompartments(self):

    allLoc = dict([])

    origGJcoords = self.network_info["origGJCoords"]

    for ID in self.neuronID:

      idxGJ1 = np.where(np.logical_and(self.network_info["synapses"][:,0]==ID,\
                              self.network_info["synapses"][:,5] == 3))

      idxGJ2 = np.where(np.logical_and(self.network_info["synapses"][:,2]==ID,\
                              self.network_info["synapses"][:,5] == 3))

      # Coordinates on either side of the gap junction
      GJcoords1 = np.array([origGJcoords[x][0:3] for x in idxGJ1[0]])
      GJcoords2 = np.array([origGJcoords[x][3:6] for x in idxGJ2[0]])

      GJid1 = np.array([origGJcoords[x][6:8] for x in idxGJ1[0]])
      GJid2 = np.array([origGJcoords[x][6:8] for x in idxGJ2[0]])

      if(GJcoords1.shape[0] == 0):
        GJcoords = GJcoords2
      elif(GJcoords2.shape[0] == 0):
        GJcoords = GJcoords1
      else:
        GJcoords = np.concatenate([GJcoords1,GJcoords2],axis=0)

      GJlocType = 4*np.ones(shape=(GJcoords.shape[0],1))

      lenGJ1 = len(GJcoords1)

      #import pdb
      #pdb.set_trace()

      if(GJcoords.shape[0] > 0):

        self.writeLog("Looking for " + str(GJcoords.shape[0]) +" gap junctions")

        # Get the compartment location of each coordinate
        GJdendLoc = self.neurons[ID].findDendCompartment(GJcoords,
                                                         GJlocType,
                                                         self.sim)
        GJdendLoc1 = GJdendLoc[:lenGJ1]
        GJdendLoc2 = GJdendLoc[lenGJ1:]

        assert(GJcoords1.shape[0] == len(GJdendLoc1))
        assert(GJcoords2.shape[0] == len(GJdendLoc2))

        for (idx,loc,id1) in zip(idxGJ1[0],GJdendLoc1,GJid1):
          allLoc[(idx,1)] = (loc,id1[0],id1[1])

        for (idx,loc,id2) in zip(idxGJ2[0],GJdendLoc2,GJid2):
          allLoc[(idx,2)] = (loc,id2[0],id2[1])

    return allLoc


  ############################################################################

  def addSynapse(self, cellIDsource, dendCompartment, sectionDist, conductance,
                 parameterID,synapseTypeID,axonDist=None):

    # You can not locate a point process at
    # position 0 or 1 if it needs an ion
    if(sectionDist == 0.0):
      sectionDist = 0.01
    if(sectionDist == 1.0):
      sectionDist = 0.99

    (channelModule,parData) = self.synapseParameters[synapseTypeID]

    syn = channelModule(dendCompartment(sectionDist))

    if(parData is not None):
      # Picking one of the parameter sets stored in parData
      parID = parameterID % len(parData)

      parSet = parData[parID]
      for par in parSet:
        if(par == "expdata" or par == "cond"):
          # expdata is not a parameter, and cond we take from synapse matrix
          continue

        try:
          # Can be value, or a tuple/list, if so second value is scale factor
          # for SI -> natural units conversion
          val = parSet[par]
          
          # Do we need to convert from SI to natural units?
          if(type(val) == tuple or type(val) == list):
             val = val[0] * val[1]
             
          setattr(syn,par,val)

          #evalStr = "syn." + par + "=" + str(parSet[par])
          # self.writeLog("Updating synapse: " + evalStr)
          # !!! Can we avoid an eval here, it is soooo SLOW
          #exec(evalStr)
        except:
          import traceback
          tstr = traceback.format_exc()
          print(tstr)
          import pdb
          pdb.set_trace()


    # Just create a default expsyn for test, will need to create proper GABA
    # synapses later
    #if(synapseType == 'ExpSyn'):
    #  syn = self.sim.neuron.h.ExpSyn(dendCompartment(sectionDist))
    #elif(synapseType == 'GABA'):
    #  syn = self.sim.neuron.h.tmGabaA(dendCompartment(sectionDist))
    #elif(synapseType == "AMPA_NMDA"):
    #  syn = self.sim.neuron.h.tmGlut(dendCompartment(sectionDist))
    #else:
    #  self.writeLog("Synapse type not implemented: ", synapseType)
    #  import pdb
    #  pdb.set_trace()

    if(axonDist is not None):
      # axon dist is in micrometer, want delay in ms
      synapseDelay = (1e3*1e-6*axonDist)/self.axonSpeed + self.synapseDelay
    else:
      synapseDelay = self.synapseDelay

    if(False):
      self.writeLog("Synapse delay: " + str(synapseDelay) + " ms")

    # What do we do if the GID does not exist?
    # print("GID exists:" + str(self.pc.gid_exists(cellIDsource)))

    if(self.isVirtualNeuron[cellIDsource]):
      # Source is a virtual neuron, need to read and connect input
      srcName = self.network_info["neurons"][cellIDsource]["name"]

      # self.writeLog("Connecting " + srcName + " to " + str(dendCompartment))

      # !!! OLD CODE, WRONG? DEL IT
      # (v,vs,spikes) = self.virtualNeurons[cellIDsource]["spikes"]
      # nc = h.NetCon(vs,syn)

      nc = self.pc.gid_connect(cellIDsource, syn)
      nc.weight[0] = conductance
      nc.delay = synapseDelay
      nc.threshold = self.spikeThreshold

      # Prevent garbage collection in python
      self.netConList.append(nc)
      self.synapseList.append(syn)


    else:
      # print("GID connect " + str(cellIDsource) + " syn: " + str(syn))
      # print("GID exists:" + str(self.pc.gid_exists(cellIDsource)))

      nc = self.pc.gid_connect(cellIDsource, syn)
      nc.weight[0] = conductance
      nc.delay = synapseDelay
      nc.threshold = self.spikeThreshold

      self.netConList.append(nc)
      self.synapseList.append(syn)

    return syn

  ############################################################################

  # Add one gap junction to specific location

  def addGapJunction(self, \
                     section, sectionDist, \
                     GIDsourceGJ, GIDdestGJ, \
                     gGapJunction=0.5e-9, \
                     GID=None): # GID unused??

    # There was a bug in neuron that affected gap junctions in parallel
    # simulations -- fixed in neuron 7.7    

    #self.writeLog("Adding src = " + str(GIDsourceGJ) + ", dest = " + str(GIDdestGJ))

    # If neuron complains, make sure you have par_ggap.mod
    GJ = h.gGapPar(section(sectionDist))
    self.gapJunctionList.append(GJ)

    self.pc.target_var(GJ._ref_vgap, GIDdestGJ)

    # !!! The line below sometimes gives this error:
    # /cfs/klemming/nobackup/h/hjorth/ChINopt/model/x86_64/special: source var gid already in use: 17124416
    # --- ok can replicate error if create 200 FS in small volume...
    # --- need to fix. HMMM how, src and dest gid has to be unique for each GJ?
    self.pc.source_var(section(sectionDist)._ref_v, GIDsourceGJ,sec=section)

    GJ.g = gGapJunction
    #print("Setting conductance: " + str(GJ.g))


  ############################################################################

  def findAllLoc(self):

    # These are the location of all post-synaptic synapses targeted
    # by neurons on this node
    destLoc = []

    # These are the post-synaptic synapses which are located on this node
    # which needs to be sent to other nodes


    # nNeurons = len(self.network_info["neuron_info"])
    # nodeLoc =  [[] for i in range(nNeurons)]
    nodeLoc =  [[] for i in range(int(self.pc.nhost()))]
    totSynapseCtr = 0

    # For each neuron on the node, find where all the synapses are located
    for cellID in self.neuronID:
      # Count how many synapses added on this neuron
      synapseCtr = 0

      # All synapses targeting this neuron
      idx = (self.network_info["synapses"][:,2] == cellID)
      synapseInfo = self.network_info["synapses"][idx,:]

      # Locations of all synapses (in the original SWC coordinate frame)
      synapseInputLoc = self.network_info["origSynapseCoords"][idx,:]

      sourceCellID = synapseInfo[:,0]
      destCellID = synapseInfo[:,2]
      locType = synapseInfo[:,4]
      synapseType = synapseInfo[:,5]

      dendLoc = self.neurons[cellID].findDendCompartment(synapseInputLoc,
                                                         locType,
                                                         self.sim)

      # We need to store the synapses on the target cells so we can find them
      # later when we connect the network together

      for srcID,destID,dLoc,sType \
          in zip(sourceCellID, destCellID, dendLoc, synapseType):

        # We need to add the synapse
        synType=self.synapseTypeLookup[sType]
        dendCompartment=dLoc[0]
        sectionDist=dLoc[1]

        if(synType == 'ExpSyn'):
          syn = self.sim.neuron.h.ExpSyn(dendCompartment(sectionDist))
        elif(synType == 'GABA'):
          syn = self.sim.neuron.h.tmGabaA(dendCompartment(sectionDist))
        elif(synType == "AMPA_NMDA"):
          syn = self.sim.neuron.h.tmGlut(dendCompartment(sectionDist))
        else:
          self.writeLog("Synapse type not implemented: ", synapseType)
          import pdb
          pdb.set_trace()

        self.synapseList.append(syn)
        self.neurons[cellID].icell.synlist.append(syn)

        synID = synapseCtr
        synapseCtr += 1
        totSynapseCtr += 1

        # Sanity check
        try:
          assert(len(self.synapseList) == totSynapseCtr)
          assert(len(self.neurons[cellID].icell.synlist) == synapseCtr)
        except Exception as e:
          self.writeLog("Sanity check: " + str(e))
          import pdb
          pdb.set_trace()

        nodeID = srcID % int(self.pc.nhost())
        assert(cellID == destID)

        nodeLoc[nodeID].append([srcID,destID,[str(dLoc[0]),dLoc[1]],int(synID)])

    self.writeLog("nhosts = " + str(int(self.pc.nhost())))
    self.writeLog("len(nodeLoc) = " + str(len(nodeLoc)))

    # import pdb
    # pdb.set_trace()

    self.writeLog("About to transfer data")
    # self.writeLog("nodeLoc = " + str(nodeLoc))

    # Transfer all the data between nodes

    self.pc.barrier()
    data = self.pc.py_alltoall(nodeLoc)
    self.pc.barrier()
    self.writeLog("All data transferred")

    # This is a list of lists, one element per node in the simulation
    # each element is a list with synapses, srcID, destID, dendLoc,synapseID
    return data

  ############################################################################

  # Wilson 2007 - GABAergic inhibition in the neostriatum
  # 80% of synapses in Striatum are glutamatergic
  # Estimated 10000 glutamate and 2000 GABA synapses per MS,
  # 325 dopamine synapses per MS
  # Wilson 1996 - 10000 spines per MS = 10000 glutamatergic inputs

  # Ingham et al 1998, approx 1 glutamatergic synapse per 0.92 mum3
  # --> ~11000 glutamate synapses per MS
  # Kemp 1971 -- The synaptic organization of the caudate nucleus (85% glu)


  def addExternalInput(self,inputFile=None):

    if(inputFile is None):
      inputFile = self.inputFile


    self.writeLog("Adding external (cortical, thalamic) input from " \
                  + inputFile)

    self.inputData = h5py.File(inputFile,'r')

    for neuronID, neuron in self.neurons.items():

      #if(neuronID != 0):
      #  self.writeLog("Skipping input temporarilly")
      #  continue

      # !!! WE ALSO NEED TO HANDLE modFile and parameterFile parameters that are
      # in inputData

      self.externalStim[neuronID] = []
      name = neuron.name

      if(str(neuronID) not in self.inputData["input"]):
        self.writeLog("Warning - No input specified for " + name)
        continue

      for inputType in self.inputData["input"][str(neuronID)]:
        neuronInput = self.inputData["input"][str(neuronID)][inputType]

        locType = 1*np.ones((neuronInput["sectionID"].shape[0],)) # Axon-Dend

        sections = self.neurons[neuronID].mapIDtoCompartment(neuronInput["sectionID"])

        # Setting individual parameters for synapses
        modFile = neuronInput["modFile"][()]
        paramList = json.loads(neuronInput["parameterList"][()])

        evalStr = "self.sim.neuron.h." + modFile
        channelModule = eval(evalStr)

        for inputID,(section,sectionX,paramID,nSpikes) \
            in enumerate(zip(sections,
                             neuronInput["sectionX"],
                             neuronInput["parameterID"],
                             neuronInput["nSpikes"])):
          # We need to find cellID (int) from neuronID (string, eg. MSD1_3)

          idx = inputID
          spikes = neuronInput["spikes"][inputID,:nSpikes] * 1e3 # Neuron uses ms
          assert (spikes >= 0).all(), \
            "Negative spike times for neuron " + str(neuronID) + " " + inputType

          if(False):
            print("Neuron " + str(neuron.name) + " receive " + str(len(spikes)) + " spikes from " + str(inputID))
            if(len(spikes) > 0):
              print("First spike at " + str(spikes[0]) + " ms")

          # Creating NEURON VecStim and vector
          # https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3125
          #import pdb
          #pdb.set_trace()
          try:
            vs = h.VecStim()
            v = h.Vector(spikes.size)
            v.from_python(spikes)
            vs.play(v)
          except:
            print("!!! If you see this, make sure that vecevent.mod is included in nrnivmodl compilation")

            import traceback
            tstr = traceback.format_exc()
            print(tstr)
            import pdb
            pdb.set_trace()


          # NEURON: You can not locate a point process at position 0 or 1
          # if it needs an ion
          if(sectionX == 0.0):
            sectionX = 0.01
          elif(sectionX == 1.0):
            sectionX = 0.99

          # !!! Parameters for the tmGlut should be possible to set in the
          # input specification !!!
          #syn = self.sim.neuron.h.tmGlut(section(sectionX))
          syn = channelModule(section(sectionX))
          nc = h.NetCon(vs,syn)

          nc.delay = 0.0
          # Should weight be between 0 and 1, or in microsiemens?
          nc.weight[0] = neuronInput["conductance"][()] * 1e6 # !! what is unit? microsiemens?
          nc.threshold = 0.1

          if(False):
            print("Weight: " + str(nc.weight[0]))

          # Get the modifications of synapse parameters, specific to
          # this synapse
          if(paramList is not None and len(paramList) > 0):
            try:
              synParams = paramList[paramID % len(paramList)]["synapse"]
            except:
              import traceback
              tstr = traceback.format_exc()
              print(tstr)
              import pdb
              pdb.set_trace()

            for par in synParams:
              if(par == "expdata"):
                # Not a parameter
                continue

              if(par == "cond"):
                # Ignoring cond value specified for synapse, using the
                # one specified in the input information instead
                continue

              try:

                evalStr = "syn." + par + "=" + str(synParams[par])
                # self.writeLog("Updating synapse: " + evalStr)
                # !!! Can we avoid an eval here, it is soooo SLOW
                exec(evalStr)

              except:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)
                import pdb
                pdb.set_trace()




          # !!! Set parameters in synParams

          # Need to save references, otherwise they will be freed
          # So sorry, but that is how neuron is
          self.externalStim[neuronID].append((v,vs,nc,syn,spikes))

          # ps = h.PatternStim()

          # HOW DO WE USE PATTERNSTIM?


  ############################################################################

  def setRestingVoltage(self,neuronID,restVolt=None):

    if(restVolt is None):
      # If no resting voltage is given, extract it from parameters
      restVolt = [x for x in self.neurons[neuronID].parameters \
                  if x["param_name"] == "v_init"][0]["value"]
      self.writeLog("Neuron " + self.neurons[neuronID].name \
                    + " resting voltage = " + str(restVolt))

    soma = [x for x in self.neurons[neuronID].icell.soma]
    axon = [x for x in self.neurons[neuronID].icell.axon]
    dend = [x for x in self.neurons[neuronID].icell.dend]

    cell = soma+axon+dend

    for sec in cell:
      for seg in sec.allseg():
        seg.v = restVolt

  ############################################################################

  def addVirtualNeuronInput(self):

    self.writeLog("Adding inputs from virtual neurons")

    assert False, "addVirtualNeuronInput not implemented"    


  ############################################################################

  # This adds external input (presumably cortical and thalamic)
  def addExternalInputOLD(self,nInputs=50,freq=5.0):
    self.writeLog("Adding external (cortical, thalamic) input (OLD VERSION)")
    assert False, "This uses neuron generated spikes, please use new function"

    for neuronID in self.neurons:
      neuron = self.neurons[neuronID]

      self.externalStim[neuronID] = []

      for i in range(0,nInputs):
        #if(neuronID != 1): # !!! TEST, only input to one neuron to see synapses work
        #  continue
        try:
          randComp = np.random.choice(neuron.icell.dend)
          randDist = np.random.random(1)

          netStim = self.sim.neuron.h.NetStim()
          netStim.start = 0
          netStim.interval = 1000.0/freq # units ms :(

          # self.writeLog("Interval " + str(netStim.interval) + " ms")

          netStim.noise = 1.0
          netStim.number = 10000

          self.externalStim[neuronID].append(netStim)

          # self.writeLog("RandComp: " + str(randComp) \
          #               + "dist: " + str(randDist[0]))

          syn = self.sim.neuron.h.tmGlut(randComp(randDist[0]))
          nc = self.sim.neuron.h.NetCon(netStim,syn)
          nc.delay = 1
          nc.weight[0] = 0.1
          nc.threshold = 0.1

          self.netConList.append(nc)
          self.synapseList.append(syn)
        except Exception as e:
          self.writeLog("Error! " + str(e))
          import pdb
          pdb.set_trace()

  ############################################################################

  # This code uses PatternStim to send input to neurons from file

  # inputFile is a csv file
  # neuronID, synapseLocation, spikeFile


#  def addExternalInput(self, inputFile=None):
#
#    if(inputFile is None):
#      inputFile = self.inputFile
#
#    for neuronID, neuron in self.neurons.items():
#      # We must store tvect and idvect, to avoid them being garbage collected
#      # same with PatternStim
#
#      # How should we specify where the input is stored, and where
#      # on the neuron it should be placed?
#      for inputFiles in neuron.input:
#
#        ps = h.PatternStim()
#
#        #!!! CODE NOT DONE, FIGURE OUT
#        assert(False)

  ############################################################################

  def centreNeurons(self,sideLen=None,neuronID=None):
    if(neuronID is None):
      neuronID = self.neuronID

    if(sideLen is None):
      return neuronID

    cID = []

    positions = self.network_info["neuronPositions"]

    centrePos = np.min(positions,axis=0)

    for nid in neuronID:
      # pos = self.network_info["neurons"][nid]["position"]
      pos = positions[nid,:]

      if(abs(pos[0]-centrePos[0]) <= sideLen
         and abs(pos[1]-centrePos[1]) <= sideLen
         and abs(pos[2]-centrePos[2]) <= sideLen):
        cID.append(nid)

    print("Centering: Keeping " + str(len(cID)) + "/" + str(len(neuronID)))

    return cID


  ############################################################################

  def addRecordingOfType(self,neuronType,nNeurons=None):

    cellID = self.snuddaLoader.getCellIDofType(neuronType=neuronType,
                                               nNeurons=nNeurons)

    self.addRecording(cellID)

  ############################################################################

  def addVoltageClamp(self,cellID,voltage,duration,res=1e-9,saveIflag=False):

    if(type(cellID) not in [list,np.ndarray]):
      cellID = [cellID]

    if(type(voltage) not in [list,np.ndarray]):
      voltage = [voltage for x in cellID]

    if(type(duration) not in [list,np.ndarray]):
      duration = [duration for x in cellID]

    if(type(res) not in [list,np.ndarray]):
      res = [res for x in cellID]

    if(saveIflag and (len(self.tSave) == 0 or self.tSave is None)):
      self.tSave = self.sim.neuron.h.Vector()
      self.tSave.record(self.sim.neuron.h._ref_t)

    for cID,v,rs,dur in zip(cellID,voltage,res,duration):

      try:
        if(not(cID in self.neuronID)):
          # Not in the list of neuronID on the worker, skip it
          continue
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)
        import pdb
        pdb.set_trace()


      self.writeLog("Adding voltage clamp to " + str(cID))
      s = self.neurons[cID].icell.soma[0]
      vc = neuron.h.SEClamp(s(0.5))
      vc.rs = rs
      vc.amp1 = v*1e3
      vc.dur1 = dur*1e3

      self.writeLog("Resistance: " + str(rs) \
                    + ", voltage: " + str(vc.amp1) + "mV")

      self.vClampList.append(vc)

      if(saveIflag):
        cur = self.sim.neuron.h.Vector()
        cur.record(vc._ref_i)
        self.iSave.append(cur)
        self.iKey.append(cID)

  ############################################################################

  def addRecording(self,cellID=None,sideLen=None):
    self.writeLog("Adding somatic recordings")

    if(cellID is None):
      cellID = self.neuronID

    # Does nothing if sideLen is not specified (otherwise, give neurons in
    # the centre)
    cellID = self.centreNeurons(sideLen=sideLen,neuronID=cellID)

    # Only include neuron IDs on this worker, ie those in self.neuronID
    # (filtering in the if statement)
    cells = dict((k,self.neurons[k]) \
                 for k in cellID if (not self.isVirtualNeuron[k] \
                                     and k in self.neuronID))

    if(len(self.tSave) == 0 or self.tSave is None):
      self.tSave = self.sim.neuron.h.Vector()
      self.tSave.record(self.sim.neuron.h._ref_t)

    for cellKey in cells:
      cell = cells[cellKey]
      try:
        v = self.sim.neuron.h.Vector()
        #import pdb
        #pdb.set_trace()
        v.record(getattr(cell.icell.soma[0](0.5),'_ref_v'))
        self.vSave.append(v)
        self.vKey.append(cellKey)
      except Exception as e:
        self.writeLog("Error: " + str(e))
        import pdb
        pdb.set_trace()

  ############################################################################


  def run(self,t=1000.0,holdV=None):

    self.setupPrintSimTime(t)

    startTime = timeit.default_timer()

    # If we want to use a non-default initialisation voltage, we need to
    # explicitly set: h.v_init
    #self.sim.neuron.h.v_init = -78
    #self.sim.neuron.h.finitialize(-78)
    if(holdV is None):
      self.sim.neuron.h.finitialize()
    else:
      self.writeLog("User override for holding voltage: " \
                    + str(holdV*1e3) + " mV")
      self.sim.neuron.h.finitialize(holdV*1e3)

    # Asked on neuron, check answer:
    # https://www.neuron.yale.edu/phpBB/viewtopic.php?f=2&t=4161&p=18021

    # Make sure all processes are synchronised
    self.pc.barrier()
    self.writeLog("Running simulation for " + str(t/1000) + " s")
    # self.sim.psolve(t)
    self.sim.run(t,dt = 0.025)
    self.pc.barrier()
    self.writeLog("Simulation done.")

    endTime = timeit.default_timer()
    self.writeLog("Simulation run time: " \
                  + str(endTime - startTime) + " s")

  ############################################################################

  def plot(self):
    import matplotlib.pyplot as pyplot
    pyplot.figure()
    for v in self.vSave:
      pyplot.plot(self.tSave,v)

    pyplot.xlabel('Time (ms)')
    pyplot.ylabel('Voltage (mV)')
    pyplot.show()
    from os.path import basename
    name = basename(self.network_info["configFile"])
    pyplot.savefig('figures/Network-voltage-trace' + name + '.pdf')

  ############################################################################

  def getSpikes(self):
    spiketrain.netconvecs_to_listoflists(self.tSpikes,self.idSpikes)

  ############################################################################

  def writeSpikes(self,outputFile=None):

    if(outputFile is None):
      outputFile = self.getSpikeFileName()

    self.writeLog("Writing spike times to " + outputFile)

    for i in range(int(self.pc.nhost())):
      self.pc.barrier() # sync all processes
      if(i == int(self.pc.id())):
        if(i == 0):
          mode = 'w'
        else:
          mode = 'a'
        with open(outputFile,mode) as spikeFile:
          for (t,id) in zip(self.tSpikes,self.idSpikes):
            spikeFile.write('%.3f\t%d\n' %(t,id))
      self.pc.barrier()

  ############################################################################

  # secList is a list of sections
  # secXList is a list of X values 0 to 1.0
  # destID is the ID of the neuron receiving synapse (one value!)
  # voxel coords are the voxel that the synapse is in

  # We want to check that voxel coords transformed to local coordinate system
  # of neuron matches with where neuron places the synapse

  def verifySynapsePlacement(self,secList,secXList,destID,voxelCoords):

    # print("Running verify synapse placement")

    simulationOrigo = self.network_info["simulationOrigo"]
    voxelSize = self.network_info["voxelSize"]
    neuronPosition = self.network_info["neurons"][destID]["position"]
    neuronRotation = self.network_info["neurons"][destID]["rotation"]

    # Transform voxel coordinates to local neuron coordinates to match neuron
    synapsePos = (voxelSize*voxelCoords+simulationOrigo-neuronPosition)*1e6

    try:
      synPosNrn = np.zeros((len(secList),3))

      for i,(sec,secX) in enumerate(zip(secList,secXList)):
        nPoints = h.n3d(sec=sec)
        arcLen = h.arc3d(nPoints-1,sec=sec)
        idx = int(np.round(secX*(nPoints-1)))
        arcLenX = h.arc3d(idx,sec=sec)

        #print("X : " + str(secX) + " = " + str(arcLenX/arcLen) + " ???")

        synPosNrn[i,0] = h.x3d(idx,sec=sec)
        synPosNrn[i,1] = h.y3d(idx,sec=sec)
        synPosNrn[i,2] = h.z3d(idx,sec=sec)
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()


    # We need to rotate the neuron to match the big simulation
    # !!! OBS, this assumes that some is in 0,0,0 local coordinates
    synPosNrnRot = np.transpose(np.matmul(neuronRotation, \
                                          np.transpose(synPosNrn)))

    synMismatch = np.sqrt(np.sum((synPosNrnRot - synapsePos)**2,axis=1))

    badThreshold = 50
    nBad = np.sum(synMismatch > badThreshold)


    if(nBad > 0):
      # If this happens, check that Neuron does not warn for removing sections
      # due to having only one point
      self.writeLog("!!! Found " + str(nBad) + " synapses on " \
                    + self.network_info["neurons"][destID]["name"] \
                    + "( " + str(destID) + ") " \
                    " that are further than " + str(badThreshold) + "mum away."\
                    + " morphology: " \
                    + self.network_info["neurons"][destID]["morphology"])

      ### DEBUG PLOT!!!

      if(True):
        import matplotlib.pyplot as plt
        plt.figure()

        somaDist = np.sqrt(np.sum(synapsePos**2,axis=1))
        plt.scatter(somaDist*1e6,synMismatch)
        plt.ion()
        plt.show()
        plt.title(self.network_info["neurons"][destID]["name"])

        from mpl_toolkits.mplot3d import Axes3D
        fig=plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(synapsePos[:,0],
                   synapsePos[:,1],
                   synapsePos[:,2],color="red")
        ax.scatter(synPosNrnRot[:,0],
                   synPosNrnRot[:,1],
                   synPosNrnRot[:,2],color="black",s=50)

        if(False):
          # Draw neuron
          allSec = [x for x in neuron.h.allsec()  if "axon" not in str(x)]
          for x in np.linspace(0,1,10):
            secPos =  np.array([[h.x3d(x,sec=sec),
                                 h.y3d(x,sec=sec),
                                 h.z3d(x,sec=sec)] \
                                for sec in allSec])

            ax.scatter(secPos[:,0],secPos[:,1],secPos[:,2],color="blue")

        import pdb
        pdb.set_trace()

    #voxelCoords *

    # Get local neuron position
    #self.neurons["position"]

    #for sec,secX
    #h.x3d(

  ############################################################################

# File format for csv voltage file:
# -1,t0,t1,t2,t3 ... (time)
# cellID,v0,v1,v2,v3, ... (voltage for cell #ID)
# repeat

  def writeVoltage(self,outputFile="save/traces/network-voltage",
                   downSampling=20):
    for i in range(int(self.pc.nhost())):
      self.pc.barrier()

      if(i == int(self.pc.id())):
        if(i == 0):
          mode = 'w'
        else:
          mode = 'a'

        with open(outputFile,mode) as voltageFile:
          if(mode == 'w'):
            voltageFile.write('-1') # Indiciate that first column is time

            for tIdx in range(0,len(self.tSave),downSampling):
              voltageFile.write(',%.4f' % self.tSave[tIdx])

          for vID, voltage in zip(self.vKey,self.vSave):
            voltageFile.write('\n%d' % vID)

            for vIdx in range(0,len(voltage),downSampling):
              voltageFile.write(',%.4f' % voltage[vIdx])

      self.pc.barrier()

  ############################################################################

  # File format for csv current file:
  # -1,t0,t1,t2,t3 ... (time)
  # cellID,i0,i1,i2,i3, ... (current for cell #ID)
  # repeat

  def writeCurrent(self,outputFile="save/traces/network-current",
                   downSampling=20):
    for i in range(int(self.pc.nhost())):
      self.pc.barrier()

      if(i == int(self.pc.id())):
        if(i == 0):
          mode = 'w'
        else:
          mode = 'a'

        with open(outputFile,mode) as currentFile:
          if(mode == 'w'):
            currentFile.write('-1') # Indiciate that first column is time

            for tIdx in range(0,len(self.tSave),downSampling):
              currentFile.write(',%.4f' % self.tSave[tIdx])

          for iID, cur in zip(self.iKey,self.iSave):
            currentFile.write('\n%d' % iID)

            for iIdx in range(0,len(cur),downSampling):
              currentFile.write(',%.4f' % cur[iIdx])

      self.pc.barrier()

##############################################################################

  def writeLog(self,text,flush=True):
    if(self.logFile is not None):
      self.logFile.write(text + "\n")
      print(text)
      if(flush):
        self.logFile.flush()
    else:
      if(self.verbose):
        print(text)

############################################################################

  def createDir(self,dirName):
    if(not os.path.isdir(dirName)):
      print("Creating " + str(dirName))
      try:
        os.makedirs(dirName)
      except:
        print("Failed to create dir. Already exists?")

############################################################################

  def addCurrentInjection(self,neuronID,startTime,endTime,amplitude):

    if neuronID not in self.neuronID:
      # The neuron ID does not exist on this worker
      return

    assert endTime > startTime, \
      "addCurrentInection: End time must be after start time"

    curStim = self.sim.neuron.h.IClamp(0.5,
                                       sec=self.neurons[neuronID].icell.soma[0])
    curStim.delay = startTime*1e3
    curStim.dur = (endTime-startTime)*1e3
    curStim.amp = amplitude*1e9 # What is units of amp?? nA??

    self.iStim.append(curStim)

  ############################################################################

  def getSpikeFileName(self):

    spikeFile = os.path.basename(self.networkFile) + "/simulation/spike-data.txt"
    return spikeFile


  ############################################################################

  def getVoltFileName(self):

    voltFile = os.path.basename(self.networkFile) + "/simulation/simulation-volt.txt"

    return voltFile

  ############################################################################

  # Use event handling

  def setupPrintSimTime(self,tMax):

    # Only have the first node print time estimates
    if(self.pc.id() == 0):
      self.tMax = tMax
      self.simStartTime = timeit.default_timer()
      self.fihTime = h.FInitializeHandler((self._setupPrintSimTimeHelper,tMax))

  ############################################################################

  def _setupPrintSimTimeHelper(self,tMax):
    updatePoints = np.arange(tMax/100., tMax, tMax/100. )
    for t in updatePoints:
      h.cvode.event(t, self.printSimTime)

  ############################################################################

  def printSimTime(self):
    curTime = timeit.default_timer()
    elapsedTime = curTime - self.simStartTime
    fractionDone = h.t/self.tMax
    timeLeft = elapsedTime * ((self.tMax - h.t)/ h.t)

    self.writeLog("%.0f%% done. Elapsed: %.1f s, estimated time left: %.1f s" \
                  % (fractionDone*100, elapsedTime, timeLeft))

  ############################################################################

  def checkMemoryStatus(self,threshold=0.1):
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

    self.writeLog(str(self.pc.id()) + ": Memory status: " \
                  + str(int(memoryRatio * 100)) + "% free")

    return memoryRatio < threshold

  ############################################################################

  def setDopamineModulation(self,sec,transientVector=[]):
    
    channelList = { 'spn':  ['naf_ms', 'kas_ms', 'kaf_ms', 'kir_ms',
                             'cal12_ms', 'cal13_ms', 'can_ms', 'car_ms'],
                    'fs':   ['kir_fs', 'kas_fs', 'kaf_fs', 'naf_fs'],
                    'chin': ['na_ch','na2_ch','kv4_ch','kir2_ch',
                             'hcn12_ch','cap_ch'],
                    'lts':  ['na3_lts','hd_lts'] }

    for cellType in channelList:
      for seg in sec:
        for mech in seg:
          if(mech.name in channelList[cellType]):
            if(len(transientVector) == 0):
              mech.damod = 1
            else:
              transientVector.play(mech._ref_damod,
                                   self.sim.neuron.h.dt)

  ############################################################################

  def applyDopamine(self,cellID=None,transientVector=[]):
    
    if(cellID is None):
      cellID = self.neuronID

    cells = dict((k,self.neurons[k]) \
                 for k in cellID if not self.isVirtualNeuron[k])

    for c in cells.values():
      for comp in [c.icell.dend, c.icell.axon, c.icell.soma]:
        for sec in comp:
          self.setDopamineModulation(sec,transientVector)

  ############################################################################

  def getPath(self,pathStr):

    return pathStr.replace("$DATA", os.path.dirname(__file__) + "/data")
          
  ############################################################################

def findLatestFile(fileMask):

  files = glob(fileMask)

  modTime = [os.path.getmtime(f) for f in files]
  idx = np.argsort(modTime)

  return files[idx[-1]]

############################################################################

#
# Test code to run a simulation

if __name__ == "__main__":

  # Please see the wrapper script snudda.py which shows how to generate
  # a network, and how to then simulate it using this script
  import sys
  if('-python' in sys.argv):
    print("Network_simulate.py called through nrniv, fixing arguments")
    pythonidx = sys.argv.index('-python')
    if(len(sys.argv) > pythonidx):
      sys.argv = sys.argv[pythonidx+1:]

  import argparse
  parser = argparse.ArgumentParser(description="Simulate network generated by Snudda")
  parser.add_argument("networkFile",help="Network model (HDF5)")
  parser.add_argument("inputFile",help="Network input (HDF5)")
  parser.add_argument("--spikesOut","--spikesout",
                      default=None,
                      help="Name of spike output file (csv)")
  parser.add_argument("--voltOut","--voltout",
                      default=None,
                      help="Name of voltage output file (csv)")
  parser.add_argument("--disableGJ",action="store_true",
                      help="Disable gap junctions")
  parser.add_argument("--time",type=float,default=1.5,
                      help="Duration of simulation in seconds")
  parser.add_argument("--verbose",action="store_true")

  # If called through "nrniv -python Network_simulate.py ..." then argparse
  # gets confused by -python flag, and we need to ignore it
  # parser.add_argument("-python",help=argparse.SUPPRESS,
  #                    action="store_true")

  args = parser.parse_args()
  networkDataFile = args.networkFile
  inputFile = args.inputFile
  logFile = os.path.dirname(args.networkFile) + "/network-simulation-log.txt"

  saveDir = os.path.dirname(args.networkFile) + "/simulation/"

  if(not os.path.exists(saveDir)):
    print("Creating directory " + saveDir)
    os.makedirs(saveDir, exist_ok=True)

  # Get the SlurmID, used in default file names
  SlurmID = os.getenv('SLURM_JOBID')

  if(SlurmID is None):
    digits = re.findall(r'\d+', inputFile)
    # Second to last digit is slurmID of old run, reuse that
    try:
      SlurmID = digits[-2]
    except:
      print("Failed to auto detect SlurmID, defaulting to 666")
      SlurmID = str(666)

  if(args.voltOut is None):
    # Do not save neuron soma voltage
    voltFile = None
  else:
    # Save neuron voltage
    if(args.voltOut == "default"):
      voltFile = saveDir + 'network-voltage-' + SlurmID + '.csv'
    else:
      voltFile = args.voltOut

  if(args.spikesOut is None or args.spikesOut == "default"):
    spikesFile = saveDir + 'network-output-spikes-' + SlurmID + '.txt'
  else:
    spikesFile = args.spikesOut

  start = timeit.default_timer()

  disableGJ = args.disableGJ
  # assert disableGJ, "Please use --disableGJ for now, need to test code"
  if(disableGJ):
    print("!!! WE HAVE DISABLED GAP JUNCTIONS !!!")

  pc = h.ParallelContext()

  sim = SnuddaSimulate(networkFile=networkDataFile,
                       inputFile=inputFile,
                       disableGapJunctions=disableGJ,
                       logFile=logFile,
                       verbose=args.verbose)

  sim.addExternalInput()
  sim.checkMemoryStatus()

  if(voltFile is not None):
    sim.addRecording(sideLen=None) # Side len let you record from a subset
    #sim.addRecordingOfType("dSPN",5) # Side len let you record from a subset
    #sim.addRecordingOfType("dSPN",2)
    #sim.addRecordingOfType("iSPN",2)
    #sim.addRecordingOfType("FSN",2)
    #sim.addRecordingOfType("LTS",2)
    #sim.addRecordingOfType("ChIN",2)

  tSim = args.time*1000 # Convert from s to ms for Neuron simulator

  sim.checkMemoryStatus()
  print("Running simulation for " + str(tSim) + " ms.")
  sim.run(tSim) # In milliseconds

  print("Simulation done, saving output")
  if(spikesFile is not None):
    sim.writeSpikes(spikesFile)

  if(voltFile is not None):
    sim.writeVoltage(voltFile)


  stop = timeit.default_timer()
  if(sim.pc.id() == 0):
    print("Program run time: " + str(stop - start ))

  # sim.plot()
  exit(0)

# Check this code example
# Why are spikes not propagated from one neuron to another
# https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=188544&file=%2FLyttonEtAl2016%2FREADME.html#tabs-2
