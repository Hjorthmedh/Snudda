# snudda_place.py
#
# Johannes Hjorth, Royal Institute of Technology (KTH)
# Human Brain Project 2019
#
# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union’s Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907, No 945539
# (Human Brain Project SGA1, SGA2, SGA3).

#

import numpy as np
import os
from collections import OrderedDict
import h5py
import json

from .Neuron_morphology import NeuronMorphology
from .RegionMesh import RegionMesh

class SnuddaPlace(object):

  ''' This code places all neurons in space, but does not setup their
      connectivity. That is done in another script. '''

  ############################################################################

  # dMin=15e-6, # MS usually less than 15 micrometers diameter

  def __init__(self,
               config_file=None,
               verbose=False,
               logFile = None,
               dView=None,
               lbView=None,
               h5libver="latest"):

    self.verbose = verbose
    self.logFile = logFile

    self.dView = dView
    self.lbView = lbView

    self.h5libver = h5libver
    self.writeLog("Using hdf5 version: " + str(self.h5libver))

    # List of all neurons
    self.neurons = []
    self.neuronPrototypes = {}

    # This defines the neuron units/channels. The dictionary lists all the
    # members of each unit, the neuronChannel gives the individual neurons
    # channel membership
    self.nPopulationUnits = 1
    self.populationUnitPlacementMethod = "random"
    self.populationUnits = dict([])
    self.populationUnit = None

    # These are the dimensions of our space, dMin also needs a "padding"
    # region outside the space where fake neurons are placed. This is to
    # avoid boundary effects, without padding we would get too high density
    # at the edges
    self.volume = dict([])

    self.config_file = config_file
    self.readConfig()


  ############################################################################

  def writeLog(self,text):
    if(self.logFile is not None):
      self.logFile.write(text + "\n")
      print(text)
    else:
      if(self.verbose):
        print(text)

  ############################################################################

  def addNeurons(self,
                 swc_filename,
                 nNeurons,
                 param_data=None,
                 mech_filename=None,
                 modulation=None,                 
                 name="Unnamed",
                 hoc=None,
                 volumeID=None,
                 rotationMode="random",
                 virtualNeuron=False,
                 axonDensity=None):

    assert volumeID is not None, "You must specify a volume for neuron " + name

    nm = NeuronMorphology(swc_filename=swc_filename,
                          param_data=param_data,
                          mech_filename=mech_filename,
                          name=name,
                          hoc=hoc,
                          virtualNeuron=virtualNeuron)

    neuronType = name.split("_")[0]
    neuronCoords = self.volume[volumeID]["mesh"].placeNeurons(nNeurons,
                                                              neuronType)

    firstAdded = True

    for coords in neuronCoords:
      # We set loadMorphology = False, to preserve memory
      # Only morphology loaded for nm then, to get axon and dend
      # radius needed for connectivity

      # Pick a random parameterset
      # parameter.json can be a list of lists, this allows you to select the
      # parameterset randomly
      # modulation.json is similarly formatted, pick a parameter set here
      parameterID = np.random.randint(1000000)
      modulationID = np.random.randint(1000000)      
      
      if(rotationMode=="random"):
        rotation = nm.randRotationMatrix()
      elif(rotationMode is None or rotationMode == ""):
        self.writeLog("Rotation mode: None (disabled) for " + name)
        rotation = np.eye(3)
      else:
        self.writeLog("Unknown rotation mode: " + str(rotationMode) \
                      + ", valid modes '' or 'random'." )
        assert False, "Unknown rotation mode: " + str(rotationMode)

      n = nm.clone(position=coords,
                   rotation=rotation,
                   loadMorphology=False,
                   parameterID=parameterID,
                   modulationID=modulationID)

      # self.writeLog("Place " + str(self.cellPos[i,:]))

      n.neuronID = len(self.neurons)
      n.volumeID = volumeID

      assert axonDensity is None or len(n.axon) == 0, \
        "!!! ERROR: Neuron: " + str(n.name) + " has both axon and axon density."

      n.axonDensity = axonDensity
      self.neurons.append(n)

      # This info is used by workers to speed things up
      if(firstAdded):
        firstAdded = False
        self.neuronPrototypes[n.name] = n



  ############################################################################

  def readConfig(self, config_file=None):

    if(config_file is None):
      config_file = self.config_file
    
    if(config_file is None):
      self.writeLog("No configuration file specified")
      exit(-1)

    if(not os.path.exists(config_file)):
      self.writeLog("Config file does not exist: " + str(config_file))
      self.writeLog("Run snudda init <your directory> first")
      exit(-1)

    self.writeLog("Parsing configuration file " + config_file)

    cfg_file = open(config_file,'r')

    try:
      config = json.load(cfg_file,object_pairs_hook=OrderedDict)
    finally:
      cfg_file.close()

    if(self.logFile is None):
      meshLogFileName = "mesh-log.txt"
    else:
      meshLogFileName = self.logFile.name + "-mesh"
    meshLogFile = open(meshLogFileName,'wt')

    # First handle volume definitions
    volumeDef = config["Volume"]

    for volumeID, volDef in volumeDef.items():

      self.volume[volumeID] = volDef

      if("meshFile" in volDef):

        assert "dMin" in volDef, "You must specify dMin if using a mesh" \
          + " for volume " + str(volumeID)

        if("meshBinWidth" not in volDef):
          self.writeLog("No meshBinWidth specified, using 1e-4")
          meshBinWidth = 1e-4
        else:
          meshBinWidth = volDef["meshBinWidth"]

        self.writeLog("Using meshBinWidth " + str(meshBinWidth))

        if("-cube-mesh-" in volDef["meshFile"]):
          self.writeLog("Cube mesh, switching to serial processing.")
          dView = None
          lbView = None
        else:
          dView = self.dView
          lbView = self.lbView

        self.volume[volumeID]["mesh"] \
          =  RegionMesh(volDef["meshFile"],
                        dView=dView,
                        lbView=lbView,
                        raytraceBorders=False,
                        dMin=volDef["dMin"],
                        binWidth=meshBinWidth,
                        logFile=meshLogFile)

      self.writeLog("Using dimensions from config file")

    if("PopulationUnits" in config):
      self.populationUnitPlacementMethod = config["PopulationUnits"]["method"]
      self.nPopulationUnits = config["PopulationUnits"]["nPopulationUnits"]

      if(self.populationUnitPlacementMethod == "populationUnitSpheres"):
        self.populationUnitRadius = config["PopulationUnits"]["radius"]
        self.populationUnitCentres = config["PopulationUnits"]["centres"]
      
    assert "Neurons" in config, \
      "No neurons defined. Is this config file old format?"

    # Read in the neurons
    for name, definition in config["Neurons"].items():

      try:
        neuronName = name
        morph = definition["morphology"]
        param = definition["parameters"]
        mech = definition["mechanisms"]

        if("modulation" in definition):
          modulation = definition["modulation"]
        else:
          # Modulation optional
          modulation = None                  
        
        num = definition["num"]
        volumeID = definition["volumeID"]

        if("neuronType" in definition):
          # type is "neuron" or "virtual" (provides input only)
          modelType = definition["neuronType"]
        else:
          modelType = "neuron"


        if 'hoc' in definition:
          hoc = definition["hoc"]
        else:
          hoc = None

        if(modelType == "virtual"):
          # Virtual neurons gets spikes from a file
          mech = ""
          hoc = ""
          virtualNeuron=True
        else:
          virtualNeuron = False

        rotationMode = definition["rotationMode"]

        if("axonDensity" in definition):
          axonDensity = definition["axonDensity"]
        else:
          axonDensity = None

        self.writeLog("Adding: " + str(num) + " " + str(neuronName))
        self.addNeurons(name=neuronName,
                        swc_filename=morph,
                        param_data=param,
                        mech_filename=mech,
                        modulation=modulation,                        
                        nNeurons=num,
                        hoc=hoc,
                        volumeID=volumeID,
                        virtualNeuron=virtualNeuron,
                        rotationMode=rotationMode,
                        axonDensity=axonDensity)

      except:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)
        print("SnuddaPlace: problem reading config file")
        import pdb
        pdb.set_trace()

    self.config_file = config_file

    # We reorder neurons, sorting their IDs after position
    self.sortNeurons()

    if(self.populationUnitPlacementMethod is not None):
      self.definePopulationUnits(method=self.populationUnitPlacementMethod)

  ############################################################################

  def allNeuronPositions(self):
    nNeurons = len(self.neurons)
    pos = np.zeros((nNeurons,3))

    for i in range(0,nNeurons):
      pos[i,:] = self.neurons[i].position

    return pos

  ############################################################################

  def allNeuronRotations(self):
    nNeurons = len(self.neurons)
    rot = np.zeros((nNeurons,3,3))

    for i in range(0,nNeurons):
      rot[i,:,:] = self.neurons[i].rotation

    return rot

  ############################################################################

  def allNeuronNames(self):
    return map(lambda x: x.name, self.neurons)

  ############################################################################

  def writeDataHDF5(self,file_name):

    self.writeLog("Writing data to HDF5 file: " + file_name)

    posFile = h5py.File(file_name,"w",libver=self.h5libver)

    with open(self.config_file,'r') as cfg_file:
      config = json.load(cfg_file,object_pairs_hook=OrderedDict)


    # Meta data
    saveMetaData = [(self.config_file,"configFile"),
                    (json.dumps(config),"config")]

    metaGroup = posFile.create_group("meta")

    for data,dataName in saveMetaData:
      metaGroup.create_dataset(dataName,data=data)

    networkGroup = posFile.create_group("network")

    # Neuron information
    neuronGroup = networkGroup.create_group("neurons")


    # If the name list is longer than 20 chars, increase S20
    nameList = [n.name.encode("ascii","ignore") for n in self.neurons]
    strType = 'S'+str(max(1,max([len(x) for x in nameList])))
    neuronGroup.create_dataset("name", (len(nameList),), strType, nameList,
                               compression="gzip")

    neuronIDlist = np.arange(len(self.neurons))
    neuronGroup.create_dataset("neuronID",(len(neuronIDlist),), \
                               'int',neuronIDlist)

    volumeIDlist = [n.volumeID.encode("ascii","ignore") \
                    for n in self.neurons]
    strTypeVID = 'S' + str(max(1,max([len(x) for x in volumeIDlist])))

    neuronGroup.create_dataset("volumeID", \
                               (len(volumeIDlist),),strTypeVID,volumeIDlist,
                               compression="gzip")

    hocList = [n.hoc.encode("ascii","ignore") for n in self.neurons]
    maxHocLen = max([len(x) for x in hocList])
    maxHocLen = max(maxHocLen,10) # In case there are none
    neuronGroup.create_dataset("hoc",(len(hocList),),'S'+str(maxHocLen), hocList,
                               compression="gzip")

    swcList = [n.swc_filename.encode("ascii","ignore") for n in self.neurons]
    maxSwcLen = max([len(x) for x in swcList])
    neuronGroup.create_dataset("morphology",(len(swcList),),'S'+str(maxSwcLen),swcList,
                               compression="gzip")

    virtualNeuronList = np.array([n.virtualNeuron for n in self.neurons],
                                 dtype=bool)
    virtualNeuron = neuronGroup.create_dataset("virtualNeuron",
                                               data=virtualNeuronList)

    # Create dataset, filled further down
    neuronPosition = neuronGroup.create_dataset("position",\
                                                (len(self.neurons),3),\
                                                "float",
                                                compression="gzip")
    neuronRotation = neuronGroup.create_dataset("rotation",\
                                                (len(self.neurons),9),\
                                                "float",\
                                                compression="gzip")

    neuronDendRadius = neuronGroup.create_dataset("maxDendRadius", \
                                                  (len(self.neurons),),\
                                                  "float",
                                                  compression="gzip")

    neuronAxonRadius = neuronGroup.create_dataset("maxAxonRadius", \
                                                  (len(self.neurons),),\
                                                  "float",
                                                  compression="gzip")

    neuronParamID = neuronGroup.create_dataset("parameterID",
                                               (len(self.neurons),),\
                                               "int",
                                               compression="gzip")
    neuronModulationID = neuronGroup.create_dataset("modulationID",
                                                    (len(self.neurons),),\
                                                    "int",
                                                    compression="gzip")

    
    for (i,n) in enumerate(self.neurons):
      neuronPosition[i]     = n.position
      neuronRotation[i]     = n.rotation.reshape(1,9)
      neuronDendRadius[i]   = n.maxDendRadius
      neuronAxonRadius[i]   = n.maxAxonRadius
      neuronParamID[i]      = n.parameterID
      neuronModulationID[i] = n.modulationID

    # Store input information
    neuronGroup.create_dataset("populationUnitID", data=self.populationUnit,dtype=int)
    neuronGroup.create_dataset("nPopulationUnits", data=self.nPopulationUnits,dtype=int)

    if(self.populationUnitPlacementMethod is not None):
      neuronGroup.create_dataset("populationUnitPlacementMethod", data=self.populationUnitPlacementMethod)
    else:
      neuronGroup.create_dataset("populationUnitPlacementMethod", data="")

    # Variable for axon density "r", "xyz" or "" (No axon density)
    axonDensityType = [n.axonDensity[0].encode("ascii","ignore") \
                       if n.axonDensity is not None \
                       else b"" \
                       for n in self.neurons]

    adStrType2 = "S"+str(max(1,max([len(x) if x is not None else 1 \
                                    for x in axonDensityType])))
    neuronGroup.create_dataset("axonDensityType", (len(axonDensityType),),
                               adStrType2,data=axonDensityType,
                               compression="gzip")


    axonDensity = [n.axonDensity[1].encode("ascii","ignore") \
                   if n.axonDensity is not None \
                   else b"" \
                   for n in self.neurons]
    adStrType = "S"+str(max(1,max([len(x) if x is not None else 1 \
                                   for x in axonDensity])))

    try:
      neuronGroup.create_dataset("axonDensity", (len(axonDensity),),
                                 adStrType,data=axonDensity,compression="gzip")
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()


    axonDensityRadius = [n.axonDensity[2] \
                         if n.axonDensity is not None \
                         and n.axonDensity[0] == "r" \
                         else np.nan for n in self.neurons]

    try:
      neuronGroup.create_dataset("axonDensityRadius",data=axonDensityRadius)
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()

    # We also need to save axonDensityBoundsXYZ, and nAxon points for the
    # non-spherical axon density option

    axonDensityBoundsXYZ = np.nan*np.zeros((len(self.neurons),6))

    for ni, n in enumerate(self.neurons):

      if(n.axonDensity is None):
        # No axon density specified, skip
        continue

      if(n.axonDensity[0] == "xyz"):

        try:
          axonDensityBoundsXYZ[ni,:] = np.array(n.axonDensity[2])
        except:
          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)

          self.writeLog("Incorrect density string: " + str(n.axonDensity))

    neuronGroup.create_dataset("axonDensityBoundsXYZ",data=axonDensityBoundsXYZ)

    posFile.close()


  ############################################################################

  def definePopulationUnits(self,method="random",nPopulationUnits=None):

    if(nPopulationUnits is None):
      nPopulationUnits = self.nPopulationUnits

    if(method == "random"):
      self.randomLabeling()
    elif(method == "populationUnitSpheres"):
      self.populationUnitSpheresLabeling(self.populationUnitCentres,self.populationUnitRadius)
    else:
      self.populationUnit = np.zeros((len(self.neurons),),dtype=int)
      self.populationUnits = dict([])
      
  ############################################################################

  def randomLabeling(self,nPopulationUnits=None):

    if(nPopulationUnits is None):
      nPopulationUnits = self.nPopulationUnits

    self.populationUnit = np.random.randint(nPopulationUnits,size=len(self.neurons))

    self.populationUnits = dict([])

    for i in range(0,nPopulationUnits):
      self.populationUnits[i] = np.where(self.populationUnit == i)[0]

  ############################################################################
  
  def populationUnitSpheresLabeling(self,populationUnitCentres,populationUnitRadius):
  
    xyz = self.allNeuronPositions()

   
    centres = np.array(populationUnitCentres)
    self.populationUnit = np.zeros((xyz.shape[0],),dtype=int)
    

    for (ctr, pos) in enumerate(xyz):
      d = [np.linalg.norm(pos-c) for c in centres]
      idx = np.argsort(d)

      if(d[idx[0]] <= populationUnitRadius):
        self.populationUnit[ctr] = idx[0] + 1 # We reserve 0 for no channel


    nPopulationUnits = np.max(self.populationUnit)+1

    for i in range(0,nPopulationUnits):
      # Channel 0 is unassigned, no channel, poor homeless neurons!
      self.populationUnits[i] = np.where(self.populationUnit == i)[0]

  ############################################################################

  def sortNeurons(self):

    # This changes the neuron IDs so the neurons are sorted along x,y or z
    xyz = self.allNeuronPositions()

    sortIdx = np.lexsort(xyz[:,[2,1,0]].transpose()) # x, y, z sort order

    self.writeLog("Re-sorting the neuron IDs after location")

    for newIdx,oldIdx in enumerate(sortIdx):
      self.neurons[oldIdx].neuronID = newIdx

    self.neurons = [self.neurons[x] for x in sortIdx]

    for idx,n in enumerate(self.neurons):
      assert idx == self.neurons[idx].neuronID, \
        "Something went wrong with sorting"

  ############################################################################

if __name__ == "__main__":

  assert False, "Please use snudda.py place networks/yournetwork"

  if(os.getenv('IPYTHON_PROFILE') is not None):
    from ipyparallel import Client
    rc = Client(profile=os.getenv('IPYTHON_PROFILE'),
                # sshserver='127.0.0.1',
                debug=False)
    print('Client IDs: ' + str(rc.ids))

    # http://davidmasad.com/blog/simulation-with-ipyparallel/
    # http://people.duke.edu/~ccc14/sta-663-2016/19C_IPyParallel.html
    print("Client IDs: " + str(rc.ids))
    dView = rc.direct_view(targets='all') # rc[:] # Direct view into clients
    lbView = rc.load_balanced_view(targets='all')
  else:
    print("No IPYTHON_PROFILE enviroment variable set, running in serial")
    dView = None
    lbView = None


  npn = SnuddaPlace(config_file="config/Network-striatum-mesh-v9-population-Units-10000-10.json",verbose=True,dView=dView,lbView=lbView)


  # Should we renumber neuron order, so that the neurons in the same hyper
  # voxel have numbers close to each other. This could speed up the merge
  # of files (but could slow down other things?)

  npn.writeDataHDF5('save/npn-test-save.hdf5')
