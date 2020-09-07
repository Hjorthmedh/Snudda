# Rewriting the create network config file to make it more general

# !!! Currently writing full path for the files in the snudda data directory
#     this makes moving between computers difficult. Fix...

#
# Add a function so that $SNUDDADATA refers to the base datapath for snudda
#

import numpy as np
import os.path
import glob
import collections
from .CreateCubeMesh import CreateCubeMesh
from .CreateSliceMesh import CreateSliceMesh

import json

class SnuddaInit(object):

  def __init__(self,structDef,configName,nPopulationUnits=1,PopulationUnitCentres="[[]]",PopulationUnitRadius=None):

    print("CreateConfig")

    self.networkData = collections.OrderedDict([])
    self.networkData["Volume"] = collections.OrderedDict([])
    self.nTotal = 0
    self.configName = configName

    if(configName is not None):
      self.basePath = os.path.dirname(configName)
    else:
      self.basePath = ""

    self.dataPath = os.path.dirname(__file__) + "/data"

    # Population Units here refer to processing units, where the neurons within a Population Unit
    # might have different connectivity than neurons belonging to different population Units
    self.networkData["PopulationUnits"] = collections.OrderedDict([])
    self.networkData["PopulationUnits"]["nPopulationUnits"] = nPopulationUnits

    useRandomPopulationUnits = False
    
    if(useRandomPopulationUnits):
      self.networkData["PopulationUnits"]["method"] = "random"
      
    else:
      self.networkData["PopulationUnits"]["method"] = "populationUnitSpheres"

      #Centre of striatum mesh is [3540e-6,4645e-6,5081e-6] - the population units will be shifted according to this coordinate

      self.networkData["PopulationUnits"]["centres"] = eval(PopulationUnitCentres)
      try:
        if(len(self.networkData["PopulationUnits"]["centres"]) == nPopulationUnits):
            pass
        else:
            raise ValueError
      except ValueError:
         print("The number of Population Units does not equal the number of centres.")
         import pdb
         pdb.set_trace()

      self.networkData["PopulationUnits"]["radius"] = PopulationUnitRadius*1e-6
      
      print("Overriding the number of population units")

      self.networkData["PopulationUnits"]["nPopulationUnits"] \
        = len(self.networkData["PopulationUnits"]["centres"])

      
      
    self.networkData["Connectivity"] = dict([])
    self.networkData["Neurons"] = dict([])

    # self.neuronTargets = collections.OrderedDict([])


    print("Using " + str(nPopulationUnits) + " Population Units")

    if(nPopulationUnits > 1):
        from scipy import spatial

        distanceBetweenPopulationUnits = spatial.distance.cdist(self.networkData["PopulationUnits"]["centres"],
                                                                self.networkData["PopulationUnits"]["centres"],
                                                                metric="euclidean")[np.triu_indices(len(self.networkData["PopulationUnits"]["centres"]),k = 1)]

        print("Using radius of Population Unit Sphere " + str(np.ceil(self.networkData["PopulationUnits"]["radius"]*1e6)) + " microns")

        print("Using distance between Population Unit Centres"  + str(distanceBetweenPopulationUnits*1e6) + " microns")

    
    self.nPopulationUnits = nPopulationUnits #5

    structFunc = { "Striatum" : self.defineStriatum,
                   "GPe" : self.defineGPe,
                   "GPi" : self.defineGPi,
                   "STN" : self.defineSTN,
                   "SNr" : self.defineSNr,
                   "Cortex" : self.defineCortex,
                   "Thalamus" : self.defineThalamus }

    if(structDef):

      for sn in structDef:
        print("Adding " + sn + " with " + str(structDef[sn]) + " neurons")
        structFunc[sn](nNeurons=structDef[sn])

      # Only write JSON file if the structDef was not empty
      self.writeJSON(self.configName)
    else:
      print("No structDef defined, not writing JSON file in init")

  ############################################################################

  # meshBinWidth is used when voxelising the mesh to determine which part of
  # space is inside the mesh. For smaller structures we might need to use
  # a smaller meshBinWidth than the default 1e-4

  def defineStructure(self,
                      structName,
                      structMesh,
                      dMin=15e-6,
                      structCentre=None,
                      sideLen=None,
                      sliceDepth=None,
                      meshBinWidth=None):

    if(structMesh == "cube"):
      assert sliceDepth is None, \
        "defineStructure: sliceDepth is not used for cubes, please set to None"
      assert sideLen is not None, \
        "defineStructure: cube needs sideLen specified"
      assert structCentre is not None, \
        "defineStructuer: cube needs a structCentre"

      structMesh = self.basePath + "/mesh/" + structName \
        + "-cube-mesh-" + str(sideLen) + ".obj"

      if(meshBinWidth is None):
        meshBinWidth = sideLen/3.0
        print("Setting meshBinWidth to " + str(meshBinWidth))

      CreateCubeMesh(fileName=structMesh,
                                    centrePoint=structCentre,
                                    sideLen=sideLen,
                                    description=structName + " cube mesh" \
                                    + ", centre = " + str(structCentre) \
                                    + ", side = " + str(sideLen))

    elif(structMesh == "slice"):

      structMesh = self.basePath + "/mesh/" + structName \
        + "-slice-mesh-150mum-depth.obj"

      # 2019-11-26 : Anya said that her sagital striatal slices
      # were 2.36 x 2.36 mm. So that can be an upper limit

      if(sideLen is None):
        sideLen = 200e-6

      if(sliceDepth is None):
        sliceDepth = 150e-6

      print("Using slice depth: " + str(sliceDepth))

      if(meshBinWidth is None):
        meshBinWidth = np.minimum(sideLen,sliceDepth)/3.0
        print("Setting meshBinWidth to " + str(meshBinWidth))

      CreateSliceMesh(fileName=structMesh,
                                      centrePoint=np.array([0,0,0]),
                                      xLen=sideLen,
                                      yLen=sideLen,
                                      zLen=sliceDepth,
                                      description=structName + " slice mesh")

    assert structName not in self.networkData["Volume"], \
      "defineStruct: Volume " + structName + " is already defined."

    self.networkData["Volume"][structName] = \
      self.defineVolume(dMin=dMin,
                        meshFile=structMesh,
                        meshBinWidth=meshBinWidth)

  ############################################################################

  def defineVolume(self,meshFile=None,dMin=15e-6,meshBinWidth=1e-4):

    vol = dict([])
    vol["type"] = "mesh"
    vol["dMin"] = dMin
    vol["meshFile"] = meshFile
    vol["meshBinWidth"] = meshBinWidth

    return vol

  ############################################################################

  # conductance and conductanceStd -- allow variation of conductances
  #
  # channelParamDictionary = dictionary specifying other parameters, such as
  #                   fascilitation and depression of AMPA/NMDA channels etc

  def addNeuronTarget(self,neuronName, targetName, connectionType,
                      distPruning,
                      f1,
                      softMax,
                      mu2,
                      a3,
                      distPruning_other=None,
                      f1_other=None,
                      softMax_other=None,
                      mu2_other=None,
                      a3_other=None,
                      conductance=[1.0e-9, 0],
                      modFile=None,
                      parameterFile=None,
                      channelParamDictionary=None):

    #if(connectionType == "GapJunction"):
    #  assert f1 is None and softMax is None and mu2 is None and a3 is None\
    #    and f1_other is None and softMax_other is None and mu2_other is None \
    #    and a3_other is None, \
    #    "addNeuronTarget: " + str(neuronName) \
    #    + ", pruning not currently available for gap junctions"

    if(parameterFile is not None):
      if(channelParamDictionary is None):
        channelParamDictionary = dict([])

      channelParamDictionary["parameterFile"] = parameterFile

    if(modFile is not None):
      if(channelParamDictionary is None):
        channelParamDictionary = dict([])

      channelParamDictionary["modFile"] = modFile


    if(type(conductance) == list):
      cond = conductance[0]
      condStd = conductance[1]
    else:
      cond = conductance
      condStd = 0

    conInfo = dict([])
    conInfo["conductance"] = [cond,condStd] # Mean, Std
    conInfo["channelParameters"] = channelParamDictionary
    pruningInfo = dict([])
    pruningInfo["f1"] = f1
    pruningInfo["softMax"] = softMax
    pruningInfo["mu2"] = mu2
    pruningInfo["a3"] = a3
    pruningInfo["distPruning"] = distPruning
    conInfo["pruning"] = pruningInfo

    # pruneInfo = (distPruning,f1,softMax,mu2,a3)

    if(distPruning_other is not None
       or f1_other is not None
       or softMax_other is not None
       or mu2_other is not None
       or a3_other is not None):

      # If any of the other varibles is set,
      # then all other "other" variables that are
      # not set assume default values

      if(distPruning_other is None):
        distPruning_other = distPruning

      if(f1_other is None):
        f1_other = f1

      if(softMax_other is None):
        softMax_other = softMax

      if(mu2_other is None):
        mu2_other = mu2

      if(a3_other is None):
        a3_other = a3

      #pruneInfo_other = (distPruning_other,
      #                   f1_other,
      #                   softMax_other,
      #                   mu2_other,
      #                   a3_other)

      pruningInfoOther = dict([])
      pruningInfoOther["f1"] = f1_other
      pruningInfoOther["softMax"] = softMax
      pruningInfoOther["mu2"] = mu2_other
      pruningInfoOther["a3"] = a3_other
      pruningInfoOther["distPruning"] = distPruning_other

      # Different pruning rules for within and between neuron units
      conInfo["pruningOther"] = pruningInfoOther

      #targetInfo = [targetName,
      #              [connectionType,cond,condStd,channelParamDictionary],
      #              pruneInfo, pruneInfo_other]
    #else:
    #  # All targets of same type are treated equally, no channels
    #  targetInfo = [targetName,
    #                [connectionType,cond,condStd,channelParamDictionary],
    #                pruneInfo]
    #
    # if(neuronName not in self.neuronTargets):
    #   self.neuronTargets[neuronName] = []

    #import pdb
    #pdb.set_trace()

    # Just make sure we are not specifying the same output twice
    # Can have GJ and synapse to same target, so moved this check to
    # Network_connect_voxel.py
    #assert targetName not in [x[0] for x in self.neuronTargets[neuronName]], \
    #  "Error " + neuronName + " already has output for " + targetName \
    #   + " specified (DUPLICATE!)"

    # self.neuronTargets[neuronName].append(targetInfo)

    # New format for connection info, now stored as dictionary
    # Json did not like typles in keys, so we separate by comma
    ntKey = neuronName + "," + targetName
    if( ntKey not in self.networkData["Connectivity"] ):
      self.networkData["Connectivity"][ntKey] = dict([])

    self.networkData["Connectivity"][ntKey][connectionType] = conInfo

  ############################################################################


  # modelType is "neuron" or "virtual" (= just provides input to network)
  # For axonDensity when it is "xyz" we assume that soma is at 0,0,0

  # neuronDir contains all the neurons in separate directories
  # Each of those directories have a config and morphology subdirectory

  def addNeurons(self,name,
                 neuronDir,
                 numNeurons, \
                 axonDensity=None,
                 modelType="neuron",
                 volumeID=None,
                 rotationMode="random"):

    if(numNeurons <= 0):
      return

    if(axonDensity is not None):
      if(axonDensity[0] == "r"):
        # Verify axon density function
        r = np.linspace(0,axonDensity[2],10)
        try:
          eval(axonDensity[1])
        except:
          print("!!! Axon density failed test: " + str(axonDensity))
          print("Inparameter: r = 1-D array of radius in meter")
      elif(axonDensity[0] == "xyz"):
        x = np.linspace(axonDensity[2][0], axonDensity[2][1],10)
        y = np.linspace(axonDensity[2][2], axonDensity[2][3],10)
        z = np.linspace(axonDensity[2][4], axonDensity[2][5],10)
        try:
          eval(axonDensity[1])
        except:
          print("!!! Axon density failed test: " + str(axonDensity))
          print("Inparameters: x,y,z three 1-D arrays (units in meter)")
          import traceback
          tstr = traceback.format_exc()

          import pdb
          pdb.set_trace()


        print("Checking boundaries, to make sure P is not too high")
        x = np.zeros((8,1))
        y = np.zeros((8,1))
        z = np.zeros((8,1))
        ctr = 0
        for xx in axonDensity[2][0:2]:
          for yy in axonDensity[2][2:4]:
            for zz in axonDensity[2][4:6]:
              x[ctr] = xx
              y[ctr] = yy
              z[ctr] = zz
              ctr += 1

        Pcorner = eval(axonDensity[1])*(3e-6**3)

        for P,xx,yy,zz in zip(Pcorner,x,y,z):
          print(name + " axon density P(" + str(xx) + "," + str(yy) \
                    + "," + str(zz) + ") = " + str(P))

        if((Pcorner > 0.01).any()):
          print("Axon density too high at boundary!!")
          print("Please increase bounding box")
          exit(-1)

        #print(str(axonDensity[3]) + " " + str(name) \
        #      + " axon points to place")

    print("Adding neurons: " + str(name) + " from dir " + str(neuronDir))

    # Find which neurons are available in neuronDir
    dirList = glob.glob(neuronDir + "/*")
    neuronFileList = []

    assert len(dirList) > 0, "Neuron dir " + str(neuronDir) + " is empty!"

    for d in dirList:

      if(os.path.isdir(d)):
        parFile = d + "/parameters.json"
        mechFile = d + "/mechanisms.json"
        modulationFile = d + "/modulation.json"
        if(not os.path.exists(modulationFile)):
          modulationFile = None
        
        swcFile = glob.glob(d + "/*swc")
        hocFile = glob.glob(d + "/*hoc")

        assert len(swcFile) == 1, "Morph dir " + d \
          + " should contain one swc file"

        assert len(hocFile) <= 1, "Morph dir " + d \
          + " contains more than one hoc file"

        if(len(hocFile) == 0):
          hocFile = [None]

        neuronFileList.append((d,
                               swcFile[0],
                               parFile,
                               mechFile,
                               modulationFile,
                               hocFile[0]))

    # First check how many unique cells we hava available, then we
    # calculate how many of each to use in simulation
    nInd = len(neuronFileList)
    nOfEachInd = np.zeros((nInd,))
    nOfEachInd[:] = int(numNeurons/nInd)
    stillToAdd = int(numNeurons - np.sum(nOfEachInd))
    addIdx = np.random.permutation(nInd)[0:stillToAdd]
    nOfEachInd[addIdx] += 1

    # Add the neurons to config

    for ctr, ((nrnDir,swcFile,parFile,mechFile,modulationFile,hocFile),num) \
        in enumerate(zip(neuronFileList,nOfEachInd)):

      if(int(num) == 0):
        continue

      uniqueName = name + "_" + str(ctr)
      cellData = dict([])

      if(not os.path.isfile(parFile) and modelType is not "virtual"):
        print("Parameter file not found: " + str(parFile))

      if(not os.path.isfile(mechFile) and modelType is not "virtual"):
        print("Mechanism file not found: " + str(mechFile))

      if(hocFile is not None and not os.path.isfile(hocFile)):
        print("Hoc file not found: " + str(hocFile))

      cellData["morphology"] = swcFile
      cellData["parameters"] = parFile
      cellData["mechanisms"] = mechFile

      if(modulationFile is not None):
        # Modulation is optional
        cellData["modulation"] = modulationFile      
      
      cellData["num"] = int(num)
      cellData["hoc"] = hocFile

      cellData["neuronType"] = modelType
      cellData["rotationMode"] = rotationMode
      cellData["volumeID"] = volumeID

      if(axonDensity is not None):
        cellData["axonDensity"] = axonDensity

      self.networkData["Neurons"][uniqueName] = cellData


  ############################################################################

  def writeJSON(self,filename):

    # !!! Dont need to do this anymore
    ## We need to copy over the target data to each neuron
    #for n in self.networkData:
    #  if(n in ["Volume","Units"]):
    #    # Non-neuron keywords, skip
    #    continue
    #
    #  nType = n.split("_")[0]
    #  if(nType in self.neuronTargets):
    #    self.networkData[n]["targets"] = self.neuronTargets[nType]
    #  else:
    #    print("No targets defined for " + str(nType))

    # import pdb
    # pdb.set_trace()

    print("Writing " + filename)

    import json
    with open(filename,'w') as f:
      json.dump(self.networkData,f,indent=4)

  ############################################################################

  # Normally: nNeurons = number of neurons set, then the fractions specified
  # fMSD1, fMSD2, fFS, fChIN, fLTS are used to  calculate the number of neurons
  # of each type.
  #
  # If nNeurons is set to None, then we use nMSD1, nMSD2, nFS, nChIN, nLTS
  # to set the number of neurons. This is useful if you want to create a small
  # test network.

  # Divide by fTot since we are not including all neurons and we want the
  # proportions to sum to 1.0 (f means fraction)

  def defineStriatum(self,nNeurons=None,
                     fMSD1=0.475,
                     fMSD2=0.475,
                     fFS=0.013,
                     fChIN=0.011,
                     fLTS=0.007,
                     nMSD1=None,
                     nMSD2=None,
                     nFS=None,
                     nChIN=None,
                     nLTS=None,
                     volumeType=None,
                     sideLen=None,
                     sliceDepth=None,
                     cellSpecDir=None,
                     neuronDensity=80500):

    getVal = lambda x : 0 if x is None else x
    if(nNeurons is None):
      self.nMSD1 = getVal(nMSD1)
      self.nMSD2 = getVal(nMSD2)
      self.nFS   = getVal(nFS)
      self.nChIN = getVal(nChIN)
      self.nLTS  = getVal(nLTS)

      self.nTotal += self.nFS + self.nMSD1 + self.nMSD2 + self.nChIN + self.nLTS
      nNeurons = self.nTotal

      if(self.nTotal <= 0):
        # No neurons specified, skipping structure
        return
    else:
      if(nNeurons <= 0):
        # No neurons specified, skipping structure
        return

      fTot = fMSD1 + fMSD2 + fFS + fChIN + fLTS

      self.nFS = np.round(fFS*nNeurons/fTot)
      self.nMSD1 = np.round(fMSD1*nNeurons/fTot)
      self.nMSD2 = np.round(fMSD2*nNeurons/fTot)
      self.nChIN = np.round(fChIN*nNeurons/fTot)
      self.nLTS = np.round(fLTS*nNeurons/fTot)

      self.nTotal += self.nFS + self.nMSD1 + self.nMSD2 + self.nChIN + self.nLTS

      if(abs(nNeurons - self.nTotal) > 5):
        print("Striatum should have " + str(nNeurons) + " but " + str(self.nTotal) \
              + " are being requested, check fractions set for defineStriatum.")

    if(volumeType == "mouseStriatum"):
      self.defineStructure(structName="Striatum",
                           structMesh=self.dataPath + "/mesh/Striatum-mesh.obj",
                           meshBinWidth=1e-4)

    elif(volumeType == "slice"):
      self.defineStructure(structName="Striatum",
                           structMesh="slice",
                           sideLen=sideLen)

    elif(nNeurons <= 1e6): #1e6
      print("Using cube for striatum")
      # 1.73 million neurons, volume of allen striatal mesh is 21.5mm3
      striatumVolume = 1e-9*(nNeurons)/neuronDensity # 80.5e3
      striatumSideLen = striatumVolume ** (1./3)
      striatumCentre = np.array([3540e-6,4645e-6,5081e-6])

      if(nNeurons < 500):
        meshBinWidth = striatumSideLen
      elif(nNeurons < 5000):
        meshBinWidth = striatumSideLen / 5
      else:
        meshBinWidth = striatumSideLen / 10

      # Reduced striatum, due to few neurons
      self.defineStructure(structName="Striatum",
                           structMesh="cube",
                           structCentre=striatumCentre,
                           sideLen=striatumSideLen,
                           meshBinWidth=meshBinWidth)

                           # meshBinWidth=5e-5)

      #import pdb
      #pdb.set_trace()
    else:
      # Default, full size striatum
      self.defineStructure(structName="Striatum",
                           structMesh=self.dataPath + "/mesh/Striatum-mesh.obj",
                           meshBinWidth=1e-4)
      
    if(cellSpecDir is None):
      csDir = self.dataPath + "/cellspecs-v2"
      #csDir = self.dataPath + "/cellspecs-var"      
    else:
      csDir = cellSpecDir

    FSdir   = csDir + "/fs"
    MSD1dir = csDir + "/dspn"
    MSD2dir = csDir + "/ispn"
    ChINdir = csDir + "/chin"
    LTSdir  = csDir + "/lts"


    self.regSize = 5

    if(self.nPopulationUnits == 1):
      self.populationUnitMSNmodifier = 1
    else:
      print("!!! OBS, modifying probaiblities within and between channe")
      self.populationUnitMSNmodifier = 0.2 # 0.2 = 20% within, 2 = 2x higher within
      print("populationUnitMSNmodifier: " + str(self.populationUnitMSNmodifier))



    # Add the neurons

    self.addNeurons(name="FSN",neuronDir=FSdir,
                    numNeurons=self.nFS, \
                    volumeID="Striatum")

    self.addNeurons(name="dSPN",neuronDir=MSD1dir,
                    numNeurons=self.nMSD1, \
                    volumeID="Striatum")

    self.addNeurons(name="iSPN",neuronDir=MSD2dir,
                    numNeurons=self.nMSD2, \
                    volumeID="Striatum")


    # ChIN axon density,
    # We start with the axon length per unit volume, then we scale it
    # to synapses per unit volume
    # This will then be used to precompute a lookup table
    # Guestimated density from Suzuki 2001, J Neurosci, figure 1bb
    # see directory morphology/ChINdensityEstimate


    # "I Janickova et al. 2017 så har de 2018 varicosities i en area på 655 um²,
    # deras slices är 70 um tjocka och om man antar att det inte är några
    # varicositites som täcker varandra så är volym-densiteten/mm³: 4.4*10⁷/mm3"
    # 1.7e6/24*0.01 = 708 ChIN per mm3
    # 4.4e7 / 708 = 62000 varicosities per ChIN
    #
    # 325 ChIN synapser per MS
    # 2-5 ChIN per MS
    # --> 65-160 synapser between a ChIN-MS pair
    # --> Each ChIN connect to 400 - 950 MS
    #
    # Number of MS within 350 micrometer radius 4*pi*(350e-6)^3/3*1.76e6/24e-9
    # --> 13100 MS reachable by ChIN at most (or rather number of MS somas
    # within radius of axonal arbour)
    # -->  3-7% connectivity probability??

    # ChINaxonDensity = ("6*5000*1e12/3*np.exp(-d/60e-6)",350e-6)

    # func type, density function, max axon radius
    ChINaxonDensity = ("r", "5000*1e12/3*np.exp(-r/120e-6)",350e-6)
    # !!! TEST
    #ChINaxonDensity = ("xyz", "2*5000*1e12/3*np.exp(-np.sqrt(x**2+y**2+z**2)/120e-6)",[-350e-6,350e-6,-350e-6,350e-6,-350e-6,350e-6])


    self.addNeurons(name="ChIN",neuronDir=ChINdir,
                    numNeurons=self.nChIN, \
                    axonDensity=ChINaxonDensity,
                    volumeID="Striatum")

    ############################################################################

    # Add LTS neuron

    # OBS, the SWC coordinates assume that the soma is centred at 0,0,0
    # Func type, Density function, [[xmin,xmax,ymin,ymax,zmin,zmax]], nAxonPoints
    # LTSaxonDensity = ("xyz", "np.exp(-(((x-500e-6)/300e-6)**2 + ((y-0)/100e-6)**2 + ((z-0)/50e-6)**2))",[150e-6,1000e-6,-200e-6,200e-6,-100e-6,100e-6],1000)
    # LTSaxonDensity = ("xyz", "np.exp(-(((x-750e-6)/500e-6)**2 + ((y-0)/100e-6)**2 + ((z-0)/50e-6)**2))",[-100e-6,2000e-6,-200e-6,200e-6,-100e-6,100e-6],1000)
    # LTSaxonDensity = ("xyz", "2000*1e12*np.exp(-(((x-750e-6)/500e-6)**2 + ((y-0)/100e-6)**2 + ((z-0)/50e-6)**2))",[-100e-6,2000e-6,-200e-6,200e-6,-100e-6,100e-6])
    # LTSaxonDensity = ("xyz", "2000*1e12*np.exp(-(((x-300e-6)/300e-6)**2 + ((y-0)/10e-6)**2 + ((z-0)/10e-6)**2))",[-300e-6,900e-6,-30e-6,30e-6,-30e-6,30e-6])


    # See plotLTSdensity.py
    # LTSaxonDensity = ("xyz", "10*3000*1e12*(0.25*np.exp(-(((x-200e-6)/100e-6)**2 + ((y-0)/50e-6)**2 + ((z-0)/20e-6)**2)) + 1*np.exp(-(((x-300e-6)/300e-6)**2 + ((y-0)/10e-6)**2 + ((z-0)/10e-6)**2)) + 1*np.exp(-(((x-700e-6)/100e-6)**2 + ((y-0)/15e-6)**2 + ((z-0)/15e-6)**2)))",[-200e-6,900e-6,-100e-6,100e-6,-30e-6,30e-6])

    LTSaxonDensity = ("xyz", "12*3000*1e12*( 0.25*np.exp(-(((x-200e-6)/100e-6)**2 + ((y-0)/50e-6)**2 + ((z-0)/30e-6)**2)) + 1*np.exp(-(((x-300e-6)/300e-6)**2 + ((y-0)/15e-6)**2 + ((z-0)/10e-6)**2)) + 1*np.exp(-(((x-700e-6)/100e-6)**2 + ((y-0)/15e-6)**2 + ((z-0)/15e-6)**2)) )",[-200e-6,900e-6,-100e-6,100e-6,-30e-6,30e-6])

    # !!! Remember to update bounding box

    self.addNeurons(name="LTS",neuronDir=LTSdir,
                    numNeurons=self.nLTS, \
                    axonDensity=LTSaxonDensity,
                    volumeID="Striatum")


    ## Define FS targets

    # Szydlowski SN, Pollak Dorocic I, Planert H, Carlen M, Meletis K,
    # Silberberg G (2013) Target selectivity of feedforward inhibition
    # by striatal fast-spiking interneurons. J Neurosci
    # --> FS does not target ChIN

    #FSDistDepPruning = "np.exp(-(0.3*d/60e-6)**2)"
    FSDistDepPruning = "np.exp(-(0.5*d/60e-6)**2)" # updated 2019-10-31
    # Temp disable dist dep pruning
    # FSDistDepPruning = None
    FSgGABA = [1.1e-9, 1.5e-9] # cond (1nS Gittis et al 2010), condStd
    FStoLTSgGABA = [1.1e-10, 1.5e-10] # cond (1nS Gittis et al 2010), condStd
    FSgGapJunction = [0.5e-9, 0.1e-9]
    # (gap junctions: 0.5nS, P=0.3 -- Galarreta Hestrin 2002, Koos Tepper 1999)
    # total 8.4nS ?? Gittis et al 2010??

    # File with FS->FS parameters (dont have one yet)
    pfFSFS = None # Gittis 2010?
    pfFSLTS = None

    #pfFSdSPN = "synapses/v1/trace_table.txt-FD-model-parameters.json"
    #pfFSiSPN = "synapses/v1/trace_table.txt-FI-model-parameters.json"
    pfFSdSPN = self.dataPath + "/synapses/v2/PlanertFitting-FD-tmgaba-fit.json"
    pfFSiSPN = self.dataPath + "/synapses/v2/PlanertFitting-FI-tmgaba-fit.json"


    # Increased from a3=0.1 to a3=0.7 to match FS-FS connectivity from Gittis
    self.addNeuronTarget(neuronName="FSN",
                         targetName="FSN",
                         connectionType="GABA",
                         distPruning=None,
                         f1=0.15, softMax=5, mu2=2, a3=1,
                         conductance=FSgGABA,
                         parameterFile=pfFSFS,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : (1.33e-3,1e3),
                                                 "tau2" : (5.7e-3,1e3) })
    # !!! Double check that channelParamDictionary works, and SI units gets
    # converted to natural units

    self.addNeuronTarget(neuronName="FSN",
                         targetName="dSPN",
                         connectionType="GABA",
                         distPruning=FSDistDepPruning,
                         f1=0.5, softMax=5, mu2=2, a3=1.0,
                         conductance=FSgGABA,
                         parameterFile=pfFSdSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : (1.2e-3,1e3),
                                                 "tau2" : (8e-3,1e3) })

    self.addNeuronTarget(neuronName="FSN",
                         targetName="iSPN",
                         connectionType="GABA",
                         distPruning=FSDistDepPruning,
                         f1=0.5, softMax=5, mu2=2, a3=0.9,
                         conductance=FSgGABA,
                         parameterFile=pfFSiSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : (1.2e-3,1e3),
                                                 "tau2" : (8e-3,1e3) })

    self.addNeuronTarget(neuronName="FSN",
                         targetName="LTS",
                         connectionType="GABA",
                         distPruning=None,
                         f1=0.15, softMax=3, mu2=2,a3=1.0,
                         conductance=FStoLTSgGABA,
                         parameterFile=pfFSLTS,
                         modFile="tmGabaA",
                         channelParamDictionary=None)



    # FS-FS gap junction, currently without pruning
    if(True):
      self.addNeuronTarget(neuronName="FSN",
                           targetName="FSN",
                           connectionType="GapJunction",
                           distPruning=None,
                           f1=0.7, softMax=8, mu2=2, a3=1.0,
                           conductance=FSgGapJunction,
                           channelParamDictionary=None)


    ## Define MSD1 targets

    # 3e-6 voxel method
    MSP11 = 1.0 #0.55
    MSP12 = 1.0 #0.20


    # Taverna 2008, fig 3E&F:
    # D1D1 22.6+/-3pS per synapse, 37+/-15 synapses (approx)
    # D2D1 24.6+/-6pS per synapse, 75+/-30 synapses (approx)
    # D2D2 24+/-1.5pS per synapse, 78+/-11 synapses (approx)

    # !!! But Taverna 2008 analyse aggregates all synapses into a conductance
    # measure?? if so, we need to divide the values by 3 or 4.
    #

    # !!! UPDATE: Assume 24pS per channel, and 10 channels per synapse

    MSD1gGABA = [0.24e-9, 0.1e-9]
    # Koos, Tepper 1999 says max 0.75nS?
    MSD1GABAfailRate = 0.7 # Taverna 2008, figure 2

    # OLD: Previously: 23pA * 50 receptors = 1.15e-9 -- Taverna 2008, fig3
    # OLD: std ~ +/- 8 receptors, we used before:  [1.15e-9, 0.18e-9]


    P11withinUnit = MSP11 * self.populationUnitMSNmodifier
    P11betweenUnit = MSP11 *(1 +(1-self.populationUnitMSNmodifier) / self.nPopulationUnits)
    P12withinUnit = MSP12 * self.populationUnitMSNmodifier
    P12betweenUnit = MSP12 *(1 +(1-self.populationUnitMSNmodifier) / self.nPopulationUnits)

    #pfdSPNdSPN = "synapses/v1/trace_table.txt-DD-model-parameters.json"
    #pfdSPNiSPN = "synapses/v1/trace_table.txt-DI-model-parameters.json"
    pfdSPNdSPN = self.dataPath +"/synapses/v2/PlanertFitting-DD-tmgaba-fit.json"
    pfdSPNiSPN = self.dataPath +"/synapses/v2/PlanertFitting-DI-tmgaba-fit.json"
    pfdSPNChIN = None


    # Argument for distance dependent SPN-SPN synapses:
    # Koos, Tepper, Wilson 2004 -- SPN-SPN more distally

    # From this paper, https://www.frontiersin.org/articles/10.3389/fnana.2010.00150/full,
    #
    # This is in contrast to the axon collateral synapses between SPNs
    # (Tunstall et al., 2002), which typically evoke significantly
    # smaller IPSPs/IPSCs than FSI-evoked synaptic responses when
    # recorded somatically (Koós et al., 2004; Tepper et al., 2004,
    # 2008; Gustafson et al., 2006) due to a combination of
    # predominantly distal synaptic locations (88%; Wilson and Groves,
    # 1980) and relatively few synaptic (2–3) connections made by each
    # SPN on each postsynaptic SPN (Koós et al., 2004)
    #
    # Also, In Kai's Thesis on the first page, He used this reference,
    # https://www.sciencedirect.com/science/article/pii/S0166223612001191?via%3Dihub,
    #


    SPN2SPNdistDepPruning = "1-np.exp(-(0.4*d/60e-6)**2)" # With Taverna conductances, we see that the response is much stronger than Planert 2010. We try to introduce distance dependent pruning to see if removing strong proximal synapses will give a better match to experimental data.

    SPN2ChINDistDepPruning = "1-np.exp(-(0.4*d/60e-6)**2)" # Chuhma about 20pA response from 10% SPN, we need to reduce activity, try dist dep pruning (already so few synapses and connectivity)

    # old f1 = 0.15
    self.addNeuronTarget(neuronName="dSPN",
                         targetName="dSPN",
                         connectionType="GABA",
                         distPruning=SPN2SPNdistDepPruning,
                         f1=0.38, softMax=3, mu2=2.4,
                         a3=P11withinUnit,
                         a3_other=P11betweenUnit,
                         conductance=MSD1gGABA,
                         parameterFile=pfdSPNdSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : (1.3e-3,1e3),
                                                 "tau2" : (12.4e-3,1e3),
                                                 "failRate" : MSD1GABAfailRate})

    # old f1 = 0.15
    self.addNeuronTarget(neuronName="dSPN",
                         targetName="iSPN",
                         connectionType="GABA",
                         distPruning=SPN2SPNdistDepPruning,
                         f1=0.20, softMax=3, mu2=2.4,
                         a3=P12withinUnit,
                         a3_other=P12betweenUnit,
                         conductance=MSD1gGABA,
                         parameterFile=pfdSPNiSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : (1.3e-3,1e3),
                                                 "tau2" : (12.4e-3,1e3),
                                                 "failRate" : MSD1GABAfailRate})

    # Doig, Magill, Apicella, Bolam, Sharott 2014:
    # 5166 +/- 285 GABA synapses on ChIN (antag att 95% av dem är från MS?)
    # 2859 +/- Assymetrical (Glut) synapses on ChIN


    # Set a3 pruning to 0.1, to remove 90% of connected pairs
    # removed softMax = 3 (want to get 5000 MSD1+D2 synapses on ChIN)

    self.addNeuronTarget(neuronName="dSPN",
                         targetName="ChIN",
                         connectionType="GABA",
                         distPruning=SPN2ChINDistDepPruning,
                         f1=0.1, softMax=3, mu2=2.4,a3=0.1,
                         conductance=MSD1gGABA,
                         parameterFile=pfdSPNChIN,
                         modFile="tmGabaA",
                         channelParamDictionary={"failRate" : MSD1GABAfailRate})



    ## Define MSD2 targets


    # 3e-6 voxel method
    MSP21 = 1.0 #0.50
    MSP22 = 1.0 # 0.95

    # OLD: 24pA * 51 receptors = 1.15e-9 -- Taverna 2008, fig3
    # OLD: std ~ +/- 10 receptors [1.24e-9, 0.24e-9]

    # Taverna 2008, fig 3E&F:
    # D1D1 22.6+/-3pS per synapse, 37+/-15 synapses (approx)
    # D2D1 24.6+/-6pS per synapse, 75+/-30 synapses (approx)
    # D2D2 24+/-1.5pS per synapse, 78+/-11 synapses (approx)

    # !!! But Taverna 2008 analyse aggregates all synapses into a conductance
    # measure?? if so, we need to divide the values by 3 or 4.
    #

    # !!! UPDATE: Assume 24pS per channel, and 10 channels per synapse
    # Because in Taverna 2008 iSPN has more receptors in total, we increase
    # softMax from 3 to 4

    MSD2gGABA = [0.24e-9, 0.1e-9]
    MSD2GABAfailRate = 0.4 # Taverna 2008, 2mM


    # Voxel method
    P21withinUnit = MSP21 * self.populationUnitMSNmodifier
    P21betweenUnit = MSP21 *(1 +(1-self.populationUnitMSNmodifier) / self.nPopulationUnits)
    P22withinUnit = MSP22 * self.populationUnitMSNmodifier
    P22betweenUnit = MSP22 *(1 +(1-self.populationUnitMSNmodifier) / self.nPopulationUnits)

    pfiSPNdSPN = self.dataPath +"/synapses/v2/PlanertFitting-ID-tmgaba-fit.json"
    pfiSPNiSPN = self.dataPath +"/synapses/v2/PlanertFitting-II-tmgaba-fit.json"
    pfiSPNChIN = None

    # GABA decay från Taverna 2008

    # old f1 = 0.15
    self.addNeuronTarget(neuronName="iSPN",
                         targetName="dSPN",
                         connectionType="GABA",
                         distPruning=SPN2SPNdistDepPruning,
                         f1=0.3, softMax=4, mu2=2.4,
                         a3=P21withinUnit,
                         a3_other=P21betweenUnit,
                         conductance=MSD2gGABA,
                         parameterFile=pfiSPNdSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : (1.3e-3,1e3),
                                                 "tau2" : (12.4e-3,1e3),
                                                 "failRate" : MSD2GABAfailRate})

    # old f1 = 0.15
    self.addNeuronTarget(neuronName="iSPN",
                         targetName="iSPN",
                         connectionType="GABA",
                         distPruning=SPN2SPNdistDepPruning,
                         f1=0.55, softMax=4, mu2=2.4,
                         a3=P22withinUnit,
                         a3_other=P22betweenUnit,
                         conductance=MSD2gGABA,
                         parameterFile=pfiSPNiSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : (1.3e-3,1e3),
                                                 "tau2" : (12.4e-3,1e3),
                                                 "failRate" : MSD2GABAfailRate})


    # See comment for dSPN to ChIN
    self.addNeuronTarget(neuronName="iSPN",
                         targetName="ChIN",
                         connectionType="GABA",
                         distPruning=SPN2ChINDistDepPruning,
                         f1=0.1, softMax=3, mu2=2.4,a3=0.1,
                         conductance=MSD2gGABA,
                         parameterFile=pfiSPNChIN,
                         modFile="tmGabaA",
                         channelParamDictionary={"failRate" : MSD2GABAfailRate})



    ## Define ChIN targets

    # Nelson AB, Hammack N, Yang CF, Shah NM, Seal RP, Kreitzer AC
    # (2014) Striatal choliner- gic interneurons Drive GABA release
    # from dopamine terminals. Neuron
    # Mamaligas, Ford 2016 -- connectivity, 2-5ChIN per MS (in slice)

    ChINgGABA = 1e-9 # If just one value given, then gSTD = 0
    ChINgACh = 1e-9 # FIXME

    # Run 1142 -- No mu2
    # Run 1150 -- Mu2 2.4
    # Run 1153 -- Mu2 D1: 5, D2: 10 (för att testa fler värden)

    # Guzman et al 2003 "Dopaminergic Modulation of Axon Collaterals Interconnecting Spiny Neurons of the Rat Striatum"
    # 325 ChIN inputs per MS (2500 * 0.13)

    # Do ChIN co-release GABA?!! otherwise should be ACh

    pfChINdSPN = None
    pfChINiSPN = None
    pfChINLTS = None

    # !!! SET RELEASE TO GABA FOR NOW
    
    # ================================================================
    # commenting gabaergic ChIN -> SPN connections Feb. 25th 2020 (RL)
    
    if(False):
        self.addNeuronTarget(neuronName="ChIN",
                             targetName="dSPN",
                             connectionType="GABA",
                             distPruning=None,
                             f1=0.5, softMax=10, mu2=15,a3=0.1, # SM 15
                             conductance=ChINgGABA,
                             parameterFile=pfChINdSPN,
                             modFile="tmGabaA",
                             channelParamDictionary=None)

        # TEST SETTING THIS TO ACh (SHOULD BE GABA), will this change?
        # !!!
        
        self.addNeuronTarget(neuronName="ChIN",
                             targetName="iSPN",
                             connectionType="GABA",
                             distPruning=None,
                             f1=0.5, softMax=10, mu2=10,a3=0.1, # SM 12
                             conductance=ChINgGABA,
                             parameterFile=pfChINiSPN,
                             modFile="tmGabaA",
                             channelParamDictionary=None)
    # ================================================================
    
    # We got an increasing connection distribution with distance, looks fishy
    # !!! Should be ACh, lets try set it to GABA and see if that changes things
    # --- trying same pruning as for ChIN to MSD2
    if(False):
      self.addNeuronTarget(neuronName="ChIN",
                           targetName="LTS",
                           connectionType="ACh",
                           distPruning=None,
                           f1=0.5, softMax=None, mu2=10,a3=None, # SM 12
                           conductance=ChINgACh,
                           parameterFile=pfChINLTS,
                           modFile="ACh", # !!! DOES NOT YET EXIST --- FIXME
                           channelParamDictionary=None)


    # !!! USE SAME PARAMS FOR FS AS FOR MS??


    # ??? ChIN does not connect to FS and MS directly ???

    # Add targets for LTS neurons

    LTSgGABA = 1e-9 # !!! FIXME
    #LTSgNO = 1e-9

    LTSDistDepPruning = "1-np.exp(-(0.4*d/60e-6)**2)" # updated 2019-10-31

    # !!! Straub, Sabatini 2016
    # No LTS synapses within 70 micrometers of proximal MS dendrite
    # !!! ADD DISTANCE DEPENDENT PRUNING

    pfLTSdSPN = None
    pfLTSiSPN = None
    pfLTSChIN = None

    self.addNeuronTarget(neuronName="LTS",
                         targetName="dSPN",
                         connectionType="GABA",
                         distPruning=LTSDistDepPruning,
                         f1=1.0, softMax=15, mu2=3, a3=0.3,
                         conductance=LTSgGABA,
                         parameterFile=pfLTSdSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : (3e-3,1e3),
                                                 "tau2" : (38e-3,1e3) })
    # LTS -> SPN, rise time 3+/-0.1 ms, decay time 38+/-3.1 ms, Straub 2016

    self.addNeuronTarget(neuronName="LTS",
                         targetName="iSPN",
                         connectionType="GABA",
                         distPruning=LTSDistDepPruning,
                         f1=1.0, softMax=15, mu2=3, a3=0.3,
                         conductance=LTSgGABA,
                         parameterFile=pfLTSiSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : (3e-3,1e3),
                                                 "tau2" : (38e-3,1e3) })

    self.addNeuronTarget(neuronName="LTS",
                         targetName="ChIN",
                         connectionType="GABA", # also NO, nitric oxide
                         distPruning=None,
                         f1=0.5, softMax=10, mu2=3, a3=0.4,
                         conductance=LTSgGABA,
                         parameterFile=pfLTSChIN,
                         modFile="tmGabaA",
                         channelParamDictionary=None)


  ############################################################################

  def defineGPe(self,nNeurons):

    if(nNeurons <= 0):
      # No neurons specified, skipping structure
      return

    self.nGPeneurons = nNeurons
    self.nTotal += nNeurons

    self.defineStructure(structName="GPe",
                         structMesh="mesh/GPe-mesh.obj")

    # !!! Need to add targets for neurons in GPe

  ############################################################################

  def defineGPi(self,nNeurons):

    if(nNeurons <= 0):
      # No neurons specified, skipping structure
      return

    self.nGPineurons = nNeurons
    self.nTotal += nNeurons

    self.defineStructure(structName="GPi",
                         structMesh="mesh/GPi-mesh.obj")

    # !!! Need to add targets for neurons in GPi

  ############################################################################

  def defineSTN(self,nNeurons):

    if(nNeurons <= 0):
      # No neurons specified, skipping structure
      return

    self.nSTNneurons = nNeurons
    self.nTotal += nNeurons

    self.defineStructure(structName="STN",
                         structMesh="mesh/STN-mesh.obj")

    # !!! Need to add targets for neurons in STN

  ############################################################################

  def defineSNr(self,nNeurons):

    if(nNeurons <= 0):
      # No neurons, skipping
      return

    self.nSNrneurons = nNeurons
    self.nTotal += nNeurons

    self.defineStructure(structName="SNr",
                         structMesh="mesh/SNr-mesh.obj",
                         meshBinWidth=1e-4)

    # !!! Need to add targets for neurons in SNr

  ############################################################################

  def defineCortex(self,nNeurons):

    if(nNeurons <= 0):
      # No neurons specified, skipping structure
      return

    # Neurons with corticostriatal axons
    self.nCortex = nNeurons

    self.nTotal += nNeurons


    # Using start location of neuron  DOI: 10.25378/janelia.5521780 for centre
    # !!! If we use a larger mesh for cortex, we will need to reduce
    #     meshBinWidth to 1e-4 (or risk getting memory error)
    self.defineStructure(structName="Cortex",
                         structMesh="cube",
                         structCentre=np.array([7067e-6,3007e-6,2570e-6]),
                         sideLen=200e-6,
                         meshBinWidth=5e-5)


    CortexDir = "morphology/InputAxons/Cortex/Reg10/"

    # Add cortex axon

    self.addNeurons("CortexAxon",CortexDir, self.nCortex, \
                    modelType="virtual",
                    rotationMode="",
                    volumeID="Cortex")

    # Define targets

    CortexGlutCond = [1e-9,0.1e-9]

    # We should have both ipsi and contra, M1 and S1 input, for now
    # picking one
    cortexSynParMS = self.dataPath \
      + "/synapses/v2/M1RH_Analysis_190925.h5-parameters-MS.json"
    cortexSynParFS = self.dataPath \
      + "synapses/v2/M1RH_Analysis_190925.h5-parameters-FS.json"

    self.addNeuronTarget(neuronName="CortexAxon",
                         targetName="dSPN",
                         connectionType="AMPA_NMDA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,a3=None,
                         parameterFile=cortexSynParMS,
                         modFile="tmGlut",
                         conductance=CortexGlutCond,
                         channelParamDictionary=None)

    self.addNeuronTarget(neuronName="CortexAxon",
                         targetName="iSPN",
                         connectionType="AMPA_NMDA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,a3=None,
                         parameterFile=cortexSynParMS,
                         modFile="tmGlut",
                         conductance=CortexGlutCond,
                         channelParamDictionary=None)

    self.addNeuronTarget(neuronName="CortexAxon",
                         targetName="FSN",
                         connectionType="AMPA_NMDA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,a3=None,
                         parameterFile=cortexSynParFS,
                         modFile="tmGlut",
                         conductance=CortexGlutCond,
                         channelParamDictionary=None)


    # !!! No input for LTS and ChIN right now...

  ############################################################################

  def defineThalamus(self,nNeurons):

    if(nNeurons <= 0):
      # No neurons specified, skipping structure
      return

    # Neurons with thalamustriatal axons
    self.nThalamus = nNeurons

    self.nTotal += nNeurons

    # Using start location of neuron DOI: 10.25378/janelia.5521765 for centre
    self.defineStructure(structName="Thalamus",
                         structMesh="cube",
                         structCentre=np.array([4997e-6,4260e-6,7019e-6]),
                         sideLen=200e-6,
                         meshBinWidth=5e-5)

    # Define neurons

    ThalamusDir = "morphology/InputAxons/Thalamus/Reg10/"

    self.addNeurons("ThalamusAxon",ThalamusDir, self.nThalamus, \
                    modelType="virtual",
                    rotationMode="",
                    volumeID="Thalamus")


    # Define targets

    thalamusSynParMS = self.dataPath \
      + "synapses/v2/TH_Analysis_191001.h5-parameters-MS.json"
    thalamusSynParFS = self.dataPath \
      + "synapses/v2/TH_Analysis_191001.h5-parameters-FS.json"


    ThalamusGlutCond = [1e-9,0.1e-9]

    self.addNeuronTarget(neuronName="ThalamusAxon",
                         targetName="dSPN",
                         connectionType="AMPA_NMDA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,a3=None,
                         conductance=ThalamusGlutCond,
                         parameterFile=thalamusSynParMs,
                         modFile="tmGlut",
                         channelParamDictionary=None)

    self.addNeuronTarget(neuronName="ThalamusAxon",
                         targetName="iSPN",
                         connectionType="AMPA_NMDA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,a3=None,
                         conductance=ThalamusGlutCond,
                         parameterFile=thalamusSynParMS,
                         modFile="tmGlut",
                         channelParamDictionary=None)

    # Picked D1 parameters, lack
    self.addNeuronTarget(neuronName="ThalamusAxon",
                         targetName="FSN",
                         connectionType="AMPA_NMDA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,a3=None,
                         conductance=ThalamusGlutCond,
                         parameterFile=thalamusSynParFS,
                         modFile="tmGlut",
                         channelParamDictionary=None)



  ############################################################################

if __name__ == "__main__":

  fullStriatum = True #False #True

  # Striatum has about 1.73 million neurons in mouse

  # Rat data (Oorschot 1996 J Comp Neurol 366)
  # --> mouse estimated from scaling down to mouse from rat
  # Striatum: 2.79M --> 1.73M
  # GPe: 46,000 --> 28500
  # GPi: 3,200 --> 2000
  # STN: 13,600 --> 8400
  # SNRpc : 7,200 --> 4500
  # SNRpr : 26,300 --> 16300

  # --> SNr = 20800


 # Nd1=Nd2=9493
 # Nfsi=400
 # Nstn=97
 # Nta=82
 # Nti=247
 # Nsnr=189


  if(fullStriatum):
    structDef = { "Striatum" : 1730000,
                  "GPe" : 28500,
                  "GPi" : 2000,
                  "SNr" : 20800,
                  "STN" : 8400,
                  "Cortex" : 1,
                  "Thalamus" : 1}

    # !!! TEMP, only do stratium for now
    structDef = { "Striatum" : 1730000,
                  "GPe" : 0,
                  "GPi" : 0,
                  "SNr" : 0,
                  "STN" : 0,
                  "Cortex" : 0,
                  "Thalamus" : 0}


  else:
    structDef = { "Striatum" : 100000,
                  "GPe" : 0,
                  "GPi" : 0,
                  "SNr" : 0,
                  "STN" : 0,
                  "Cortex" : 0,
                  "Thalamus" : 0}


  nTotals = 0
  for x in structDef:
    nTotals += structDef[x]

  fName = "config/basal-ganglia-config-" + str(nTotals) + ".json"

  SnuddaInit(structDef=structDef,configName=fName,nPopulationUnits=1)
