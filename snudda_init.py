# Rewriting the create network config file to make it more general

import numpy as np
import os.path
import glob
import collections
import CreateCubeMesh

import json

class SnuddaInit(object):  
    
  def __init__(self,structDef,configName,nChannels=1):

    print("CreateConfig")
    
    self.networkData = collections.OrderedDict([])
    self.networkData["Volume"] = collections.OrderedDict([])
    self.nTotal = 0
    self.configName = configName

    # Channels here refer to processing units, where the neurons within a channel
    # might have different connectivity than neurons belonging to different channels
    self.networkData["Channels"] = collections.OrderedDict([])
    self.networkData["Channels"]["nChannels"] = nChannels
    self.networkData["Channels"]["method"] = "random"

    self.neuronTargets = collections.OrderedDict([])

    print("Using " + str(nChannels) + " functional channels")
    self.nChannels = nChannels #5

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
                      meshBinWidth=1e-4):

    if(structMesh == "cube"):
      assert sideLen is not None, \
        "defineStructure: cube needs sideLen specified"
      assert structCentre is not None, \
        "defineStructuer: cube needs a structCentre"
      
      structMesh = "mesh/" + structName + "-cube-mesh-" + str(sideLen) + ".obj"

      if(meshBinWidth is None):
        meshBinWidth = sideLen/20.0
        print("Setting meshBinWidth to " + str(meshBinWidth))
      
      CreateCubeMesh.CreateCubeMesh(fileName=structMesh,
                                    centrePoint=structCentre,
                                    sideLen=sideLen,
                                    description=structName + " cube mesh" \
                                    + ", centre = " + str(structCentre) \
                                    + ", side = " + str(sideLen))

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

    if(connectionType == "GapJunction"):
      assert f1 is None and softMax is None and mu2 is None and a3 is None\
        and f1_other is None and softMax_other is None and mu2_other is None \
        and a3_other is None, \
        "addNeuronTarget: " + str(neuronName) \
        + ", pruning not currently available for gap junctions"

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
    
    pruneInfo = (distPruning,f1,softMax,mu2,a3)
    
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
      
      pruneInfo_other = (distPruning_other,
                         f1_other,
                         softMax_other,
                         mu2_other,
                         a3_other)
      
      # Different pruning rules for within and between neuron channels
      targetInfo = [targetName,
                    [connectionType,cond,condStd,channelParamDictionary],
                    pruneInfo, pruneInfo_other]
    else:
      # All targets of same type are treated equally, no channels 
      targetInfo = [targetName,
                    [connectionType,cond,condStd,channelParamDictionary],
                    pruneInfo]

    if(neuronName not in self.neuronTargets):
      self.neuronTargets[neuronName] = []

    #import pdb
    #pdb.set_trace()
      
    # Just make sure we are not specifying the same output twice
    # Can have GJ and synapse to same target, so moved this check to
    # Network_connect_voxel.py
    #assert targetName not in [x[0] for x in self.neuronTargets[neuronName]], \
    #  "Error " + neuronName + " already has output for " + targetName \
    #   + " specified (DUPLICATE!)"
      
    self.neuronTargets[neuronName].append(targetInfo)

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

        swcFile = glob.glob(d + "/*swc")
        hocFile = glob.glob(d + "/*hoc")

        assert len(swcFile) == 1, "Morph dir " + d \
          + " should contain one swc file"

        assert len(hocFile) <= 1, "Morph dir " + d \
          + " contains more than one hoc file"
        
        if(len(hocFile) == 0):
          hocFile = [None]
          
        neuronFileList.append((d,swcFile[0],parFile,mechFile,hocFile[0]))
        
    # First check how many unique cells we hava available, then we
    # calculate how many of each to use in simulation
    nInd = len(neuronFileList)
    nOfEachInd = np.zeros((nInd,))
    nOfEachInd[:] = int(numNeurons/nInd)
    stillToAdd = int(numNeurons - np.sum(nOfEachInd))
    addIdx = np.random.permutation(nInd)[0:stillToAdd]
    nOfEachInd[addIdx] += 1

    # Add the neurons to config
    
    for ctr, ((nrnDir,swcFile,parFile,mechFile,hocFile),num) \
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
      cellData["num"] = int(num)
      cellData["hoc"] = hocFile
      
      cellData["neuronType"] = modelType
      cellData["rotationMode"] = rotationMode
      cellData["volumeID"] = volumeID

      if(axonDensity is not None):
        cellData["axonDensity"] = axonDensity
      
      self.networkData[uniqueName] = cellData
       

  ############################################################################

  def writeJSON(self,filename):

    # We need to copy over the target data to each neuron
    for n in self.networkData:
      if(n in ["Volume","Channels"]):
        # Non-neuron keywords, skip
        continue

      nType = n.split("_")[0]      
      if(nType in self.neuronTargets):
        self.networkData[n]["targets"] = self.neuronTargets[nType]
      else:
        print("No targets defined for " + str(nType))
        
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
                     nLTS=None):

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

    if(nNeurons <= 1e6): #1e6
      print("Using cube for striatum")
      # 1.73 million neurons, volume of allen striatal mesh is 21.5mm3
      striatumVolume = 1e-9*(nNeurons)/80.5e3 
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
                           structMesh="mesh/Striatum-mesh.obj",
                           meshBinWidth=1e-4)

    FSdir   = "cellspecs/fs"
    MSD1dir = "cellspecs/dspn"
    MSD2dir = "cellspecs/ispn"
    ChINdir = "cellspecs/chin"
    LTSdir  = "cellspecs/lts"

    #FSdir   = "cellspecs.sfn/fs"
    #MSD1dir = "cellspecs.sfn/dspn"
    #MSD2dir = "cellspecs.sfn/ispn"
    #ChINdir = "cellspecs.sfn/chin"
    #LTSdir  = "cellspecs.sfn/lts"
    
    self.regSize = 5

    if(self.nChannels == 1):
      self.channelMSNmodifier = 1
    else:
      print("!!! OBS, modifying probaiblities within and between channels")
      self.channelMSNmodifier = 0.2 # 0.2 = 20% within, 2 = 2x higher within
      print("channelMSNmodifier: " + str(self.channelMSNmodifier))

    

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
    ChINaxonDensity = ("r", "2*5000*1e12/3*np.exp(-r/120e-6)",350e-6)
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
    LTSaxonDensity = ("xyz", "3000*1e12*(0.25*np.exp(-(((x-200e-6)/100e-6)**2 + ((y-0)/50e-6)**2 + ((z-0)/20e-6)**2)) + 1*np.exp(-(((x-300e-6)/300e-6)**2 + ((y-0)/10e-6)**2 + ((z-0)/10e-6)**2)) + 1*np.exp(-(((x-700e-6)/100e-6)**2 + ((y-0)/15e-6)**2 + ((z-0)/15e-6)**2)))",[-200e-6,900e-6,-100e-6,100e-6,-30e-6,30e-6])    

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
    
    FSDistDepPruning = "np.exp(-(0.3*d/60e-6)**2)"
    # Temp disable dist dep pruning
    # FSDistDepPruning = None
    FSgGABA = [1.1e-9, 1.5e-9] # cond (1nS Gittis et al 2010), condStd
    FSgGapJunction = [0.5e-9, 0.1e-9]
    # (gap junctions: 0.5nS, P=0.3 -- Galarreta Hestrin 2002, Koos Tepper 1999)
    # total 8.4nS ?? Gittis et al 2010??
    
    # File with FS->FS parameters (dont have one yet)
    pfFSFS = None # Gittis 2010?
    pfFSLTS = None

    pfFSdSPN = "synapses/trace_table.txt-FD-model-parameters.json"
    pfFSiSPN = "synapses/trace_table.txt-FI-model-parameters.json"

    
    # Increased from a3=0.1 to a3=0.7 to match FS-FS connectivity from Gittis
    self.addNeuronTarget(neuronName="FSN",
                         targetName="FSN",
                         connectionType="GABA",
                         distPruning=None,
                         f1=1, softMax=8, mu2=2, a3=0.7,
                         conductance=FSgGABA,
                         parameterFile=pfFSFS,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : 1.33e-3,
                                                 "tau2" : 5.7e-3 })
    # !!! Double check that channelParamDictionary works, and SI units gets
    # converted to natural units

    self.addNeuronTarget(neuronName="FSN",
                         targetName="dSPN",
                         connectionType="GABA",
                         distPruning=FSDistDepPruning,
                         f1=1, softMax=8, mu2=2, a3=None, # mu2 was 2
                         conductance=FSgGABA,
                         parameterFile=pfFSdSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : 1.2e-3,
                                                 "tau2" : 8e-3 })

    self.addNeuronTarget(neuronName="FSN",
                         targetName="iSPN",
                         connectionType="GABA",
                         distPruning=FSDistDepPruning,
                         f1=1, softMax=8, mu2=2, a3=None, # mu2 was 2
                         conductance=FSgGABA,
                         parameterFile=pfFSiSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : 1.2e-3,
                                                 "tau2" : 8e-3 })

    self.addNeuronTarget(neuronName="FSN",
                         targetName="LTS",
                         connectionType="GABA",
                         distPruning=None,
                         f1=1.0, softMax=8, mu2=2,a3=0.63,
                         conductance=FSgGABA,
                         parameterFile=pfFSLTS,
                         modFile="tmGabaA",
                         channelParamDictionary=None)

    

    # FS-FS gap junction, currently without pruning
    if(True):
      self.addNeuronTarget(neuronName="FSN",
                           targetName="FSN",
                           connectionType="GapJunction",
                           distPruning=None,
                           f1=None, softMax=None, mu2=None, a3=None,
                           conductance=FSgGapJunction,
                           channelParamDictionary=None)

    
    ## Define MSD1 targets

    # 2e-6 voxel method
    # MSP11 = 0.25
    # MSP12 = 0.17

    # 3e-6 voxel method
    MSP11 = 0.17
    MSP12 = 0.085 #0.1 # 0.14 then 0.16 old

    # 23pA * 50 receptors = 1.15e-9 -- Taverna 2008, fig3
    # std ~ +/- 8 receptors
    MSD1gGABA = [1.15e-9, 0.18e-9]
    # Koos, Tepper 1999 says max 0.75nS?
    
    
    P11withinChannel = MSP11 * self.channelMSNmodifier
    P11betweenChannel = MSP11 *(1 +(1-self.channelMSNmodifier) / self.nChannels)
    P12withinChannel = MSP12 * self.channelMSNmodifier
    P12betweenChannel = MSP12 *(1 +(1-self.channelMSNmodifier) / self.nChannels)

    pfdSPNdSPN = "synapses/trace_table.txt-DD-model-parameters.json"
    pfdSPNiSPN = "synapses/trace_table.txt-DI-model-parameters.json"
    pfdSPNChIN = None
    
    self.addNeuronTarget(neuronName="dSPN",
                         targetName="dSPN",
                         connectionType="GABA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,
                         a3=P11withinChannel,
                         a3_other=P11betweenChannel,
                         conductance=MSD1gGABA,
                         parameterFile=pfdSPNdSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : 1.3e-3,
                                                 "tau2" : 12.4e-3 })

    self.addNeuronTarget(neuronName="dSPN",
                         targetName="iSPN",
                         connectionType="GABA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,
                         a3=P12withinChannel,
                         a3_other=P12betweenChannel,
                         conductance=MSD1gGABA,
                         parameterFile=pfdSPNiSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : 1.3e-3,
                                                 "tau2" : 12.4e-3 })

    # Doig, Magill, Apicella, Bolam, Sharott 2014:
    # 5166 +/- 285 GABA synapses on ChIN (antag att 95% av dem är från MS?)
    # 2859 +/- Assymetrical (Glut) synapses on ChIN


    # Set a3 pruning to 0.1, to remove 90% of connected pairs
    # removed softMax = 3 (want to get 5000 MSD1+D2 synapses on ChIN)
    
    self.addNeuronTarget(neuronName="dSPN",
                         targetName="ChIN",
                         connectionType="GABA",
                         distPruning=None,
                         f1=1.0, softMax=None, mu2=2.4,a3=0.1,
                         conductance=MSD1gGABA,
                         parameterFile=pfdSPNChIN,
                         modFile="tmGabaA",
                         channelParamDictionary=None)

                         

    ## Define MSD2 targets

    # 2e-6 voxels
    # MSP21 = 0.4
    # MSP22 = 0.8

    # 3e-6 voxel method
    MSP21 = 0.23
    MSP22 = 0.4

    # 24pA * 51 receptors = 1.15e-9 -- Taverna 2008, fig3
    # std ~ +/- 10 receptors
    MSD2gGABA = [1.24e-9, 0.24e-9]
    
    # Voxel method 
    P21withinChannel = MSP21 * self.channelMSNmodifier
    P21betweenChannel = MSP21 *(1 +(1-self.channelMSNmodifier) / self.nChannels)
    P22withinChannel = MSP22 * self.channelMSNmodifier
    P22betweenChannel = MSP22 *(1 +(1-self.channelMSNmodifier) / self.nChannels)

    pfiSPNdSPN = "synapses/trace_table.txt-ID-model-parameters.json"
    pfiSPNiSPN = "synapses/trace_table.txt-II-model-parameters.json"    
    pfiSPNChIN = None

    
    self.addNeuronTarget(neuronName="iSPN",
                         targetName="dSPN",
                         connectionType="GABA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,
                         a3=P21withinChannel,
                         a3_other=P21betweenChannel,
                         conductance=MSD2gGABA,
                         parameterFile=pfiSPNdSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : 1.3e-3,
                                                 "tau2" : 12.4e-3 })

    self.addNeuronTarget(neuronName="iSPN",
                         targetName="iSPN",
                         connectionType="GABA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,
                         a3=P22withinChannel,
                         a3_other=P22betweenChannel,
                         conductance=MSD2gGABA,
                         parameterFile=pfiSPNiSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : 1.3e-3,
                                                 "tau2" : 12.4e-3 })

    self.addNeuronTarget(neuronName="iSPN",
                         targetName="ChIN",
                         connectionType="GABA",
                         distPruning=None,
                         f1=1.0, softMax=None, mu2=2.4,a3=0.1,
                         conductance=MSD2gGABA,
                         parameterFile=pfiSPNChIN,
                         modFile="tmGabaA",
                         channelParamDictionary=None)

    
    
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
    
    self.addNeuronTarget(neuronName="ChIN",
                         targetName="dSPN",
                         connectionType="GABA",
                         distPruning=None,
                         f1=1.0, softMax=20, mu2=15,a3=None,
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
                         f1=1.0, softMax=20, mu2=10,a3=None,
                         conductance=ChINgGABA,
                         parameterFile=pfChINiSPN,
                         modFile="tmGabaA",
                         channelParamDictionary=None)

    # We got an increasing connection distribution with distance, looks fishy
    # !!! Should be ACh, lets try set it to GABA and see if that changes things
    # --- trying same pruning as for ChIN to MSD2
    self.addNeuronTarget(neuronName="ChIN",
                         targetName="LTS",
                         connectionType="ACh",
                         distPruning=None,
                         f1=1.0, softMax=20, mu2=10,a3=None,
                         conductance=ChINgACh,
                         parameterFile=pfChINLTS,
                         modFile="tmGabaA",
                         channelParamDictionary=None)
    

    # !!! USE SAME PARAMS FOR FS AS FOR MS??
    

    # ??? ChIN does not connect to FS and MS directly ???

    # Add targets for LTS neurons

    LTSgGABA = 1e-9 # !!! FIXME
    LTSgNO = 1e-9

    # !!! Straub, Sabatini 2016
    # No LTS synapses within 70 micrometers of proximal MS dendrite
    # !!! ADD DISTANCE DEPENDENT PRUNING

    pfLTSdSPN = None
    pfLTSiSPN = None
    pfLTSChIN = None
    
    self.addNeuronTarget(neuronName="LTS",
                         targetName="dSPN",
                         connectionType="GABA",
                         distPruning=None,
                         f1=1.0, softMax=10, mu2=2, a3=None,
                         conductance=LTSgGABA,
                         parameterFile=pfLTSdSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : 3e-3,
                                                 "tau2" : 38e-3 })
    # LTS -> SPN, rise time 3+/-0.1 ms, decay time 38+/-3.1 ms, Straub 2016
    
    self.addNeuronTarget(neuronName="LTS",
                         targetName="iSPN",
                         connectionType="GABA",
                         distPruning=None,
                         f1=1.0, softMax=10, mu2=2, a3=None,
                         conductance=LTSgGABA,
                         parameterFile=pfLTSiSPN,
                         modFile="tmGabaA",
                         channelParamDictionary={"tau1" : 3e-3,
                                                 "tau2" : 38e-3 })

    self.addNeuronTarget(neuronName="LTS",
                         targetName="ChIN",
                         connectionType="GABA", # NO, nitric oxide
                         distPruning=None,
                         f1=1.0, softMax=10, mu2=2, a3=None,
                         conductance=LTSgNO,
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
                         structMesh="mesh/SNr-mesh.obj")
    
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
    cortextSynPar = "synapses/M1_Analysis_RH_extr_UP.pxp-traceList-MSND1-require-H20-model-parameters.json"
    
    self.addNeuronTarget(neuronName="CortexAxon",
                         targetName="dSPN",
                         connectionType="AMPA_NMDA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,a3=None,
                         parameterFile=cortexSynPar,
                         modFile="tmGlut",
                         conductance=CortexGlutCond,
                         channelParamDictionary=None)

    self.addNeuronTarget(neuronName="CortexAxon",
                         targetName="iSPN",
                         connectionType="AMPA_NMDA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,a3=None,
                         parameterFile=cortexSynPar,
                         modFile="tmGlut",
                         conductance=CortexGlutCond,
                         channelParamDictionary=None)

    self.addNeuronTarget(neuronName="CortexAxon",
                         targetName="FSN",
                         connectionType="AMPA_NMDA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,a3=None,
                         parameterFile=cortexSynPar,
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

    thalamusSynParD1 = "synapses/TH_Analysis_extr_UP.pxp-traceList-MSND1-require-H20-model-parameters.json"
    thalamusSynParD2 = "synapses/TH_Analysis_extr_UP.pxp-traceList-MSND2-require-H20-model-parameters.json"
    thalamusSynParFS = "synapses/TH_Analysis_extr_UP.pxp-traceList-FS-require-H20-model-parameters.json"

    
    ThalamusGlutCond = [1e-9,0.1e-9]
    
    self.addNeuronTarget(neuronName="ThalamusAxon",
                         targetName="dSPN",
                         connectionType="AMPA_NMDA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,a3=None,
                         conductance=ThalamusGlutCond,
                         parameterFile=thalamusSynParD1,
                         modFile="tmGlut",
                         channelParamDictionary=None)

    self.addNeuronTarget(neuronName="ThalamusAxon",
                         targetName="iSPN",
                         connectionType="AMPA_NMDA",
                         distPruning=None,
                         f1=1.0, softMax=3, mu2=2.4,a3=None,
                         conductance=ThalamusGlutCond,
                         parameterFile=thalamusSynParD2,
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
    
  SnuddaInit(structDef=structDef,configName=fName,nChannels=1)
