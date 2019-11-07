# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union’s Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).
#
# etc...
#

import numpy as np

class NeuronMorphology(object):

  # axonStumpIDFlag should be True if running Network_simulate.py
  # it should be False if we are running Neurodamus simulation.
  
  def __init__(self,
               name = None,
               position = np.zeros((1,3)),
               rotation = None, #np.eye(3),
               swc_filename = None,
               param_filename = None,
               param_data = None,
               mech_filename = None,
               verbose=False,
               loadMorphology=True,
               hoc=None,
               colour=None,
               useCache=True,
               pickleVersion=-1,
               logFile = None,
               virtualNeuron = False,
               axonStumpIDFlag = False ):

    self.cacheVersion = 0.9
    
    self.position = np.copy(np.array(position))
    
    if(rotation is not None):
      self.rotation = np.copy(np.array(rotation))
    else:
      self.rotation = None

    self.soma = []
    self.axon = []
    self.dend = []  # 0,1,2: x,y,z  3: radie, 4: dist to soma (all in meters)


    self.axonDensityType = None
    self.dendDensity = None
    self.axonDensity = None
    self.axonDensityBoundsXYZ = None
    
    self.voxelSize = 5e6
    self.densityBinSize = 10e-6

    self.loadMorphology = loadMorphology

    # This tells if axon is indexed fully or if only as a stump
    # this affects sectionID for both axon and dendrites
    self.axonStumpIDFlag = axonStumpIDFlag
    
    # Meta data
    self.name = name      
    self.swc_filename = swc_filename
    self.param_filename = param_filename
    self.param_data = param_data
    self.mech_filename = mech_filename
    self.verbose = verbose
    self.useCache = useCache
    self.pickleVersion = pickleVersion    
    self.logFile = logFile
    self.virtualNeuron=virtualNeuron

    self.rotatedFlag = False
    
    self.cache_filename = swc_filename.replace('.swc','-cache.pickle')
    assert(self.cache_filename != swc_filename)

    # This is used for Neurodamus, which instantiates through hoc files
    if(hoc is None):
      hoc = ""
      
    self.hoc = hoc
    
    # This is useful when determining connectivity, to exclude pairs
    self.maxAxonRadius = 0
    self.maxDendRadius = 0
    
    # Telling how the different points link together into lines
    self.axonLinks = [] # These should never be changed after CLONE
    self.dendLinks = []

    self.dendSecID = []
    self.dendSecX = []

    if(colour is None):
      self.colour = np.random.random((3,))
    else:
      self.colour = colour
    
    if(loadMorphology):
      # This loads, rotates and places neuron
      self.loadNeuronMorphology()


      
  ############################################################################
    
  def loadNeuronMorphology(self):

    if(self.useCache):
      if(self.cacheExist()):
        # Load existing cache
        try:
          self.loadCache() 
        except Exception as e:
          
          import traceback
          tstr = traceback.format_exc()
          print(tstr)
          
          self.writeLog("!!! Failed to read cache file, loading: " \
                        + self.swc_filename)
          self.loadSWC(self.swc_filename)
          self.saveCache()

      else:
        self.writeLog("No cache found, create it.")
        # Load SWC and save cache file
        self.loadSWC(self.swc_filename)
        self.saveCache()
    else:
      # Load SWC file
      self.writeLog("Ignoring old cache, rewriting new cache file")
      self.loadSWC(self.swc_filename)
      self.saveCache()

    self.place() # Updates position and rotation

    # Remove axonStumpIDFlag completely later...
    assert not self.axonStumpIDFlag, \
      "axonStumpFlag is depricated, should be off"
    
  ############################################################################

  def clone(self,
            loadMorphology=None, # True or False, None = same as parent
            position=np.zeros((1,3)),
            rotation=None):

    if(loadMorphology is None):
      loadMorphology = self.loadMorphology

    # If these are explicitly set to None, reuse to original coordinates
    # and rotation
    if(position is None):
      position=self.position

    if(rotation is None):
      rotation=self.rotation
      
    newNeuron = NeuronMorphology(name=self.name,
                                 position=position,
                                 rotation=rotation,
                                 swc_filename=self.swc_filename,
                                 param_filename=self.param_filename,
                                 param_data=self.param_data,
                                 mech_filename=self.mech_filename,
                                 verbose=self.verbose,
                                 loadMorphology=False,
                                 hoc=self.hoc,
                                 virtualNeuron = self.virtualNeuron)

    if(loadMorphology):
      # Set the flag
      newNeuron.loadMorphology = loadMorphology

      # Copy the data
      newNeuron.axon = np.copy(self.axon)
      newNeuron.dend = np.copy(self.dend)
      newNeuron.soma = np.copy(self.soma)

      # Warn the user if the neuron is already rotated
      newNeuron.rotatedFlag = self.rotatedFlag
      
      newNeuron.place()

      # These dont change either, so skip np.copy
      newNeuron.axonLinks = self.axonLinks
      newNeuron.dendLinks = self.dendLinks
      newNeuron.dendSecX = self.dendSecX
      newNeuron.dendSecID = self.dendSecID

      newNeuron.axonStumpIDFlag = self.axonStumpIDFlag
      
    newNeuron.maxAxonRadius = self.maxAxonRadius
    newNeuron.maxDendRadius = self.maxDendRadius

    if(self.dendDensity is not None):
      newNeuron.dendDensity = self.dendDensity

    if(self.axonDensity is not None):
      newNeuron.axonDensity = self.axonDensity

    newNeuron.voxelSize = self.voxelSize

    if(self.axonDensityType is not None):
      newNeuron.axonDensityType = self.axonDensityType

    if(self.axonDensityBoundsXYZ is not None):
       newNeuron.axonDensityBoundsXYZ = self.axonDensityBoundsXYZ
    
    return newNeuron

  ############################################################################
  
  def writeLog(self,text):
    if(self.logFile is not None):
      self.logFile.write(text + "\n")
      print(text)
    else:
      if(self.verbose):
        print(text)
        
  ############################################################################
  
  # http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
  def randRotationMatrix(self, deflection=1.0, randnums=None):
    """
    Creates a random rotation matrix.
    
    deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
    rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
    """

    # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
  
    if randnums is None:
      randnums = np.random.uniform(size=(3,))
    
    theta, phi, z = randnums
    
    theta = theta * 2.0*deflection*np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0*np.pi  # For direction of pole deflection.
    z = z * 2.0*deflection  # For magnitude of pole deflection.
    
    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
  
    r = np.sqrt(z)
    Vx, Vy, Vz = V = (
      np.sin(phi) * r,
      np.cos(phi) * r,
      np.sqrt(2.0 - z)
    )
    
    st = np.sin(theta)
    ct = np.cos(theta)
    
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
    
    # Construct the rotation matrix  ( V Transpose(V) - I ) R.  
    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M

  ############################################################################
  
  # We can specify a position and rotation
  def place(self, rotation=None, position=None):

    
    if(self.rotatedFlag):
      self.writeLog("!!! WARNING, rotating a rotated neuron...")
    
    if(rotation is None):
      rotation = self.rotation
    elif(type(rotation) is not np.ndarray):
      rotation = np.array(rotation)

    if(position is None):
      position = self.position
    elif(type(position) is not np.ndarray):
      position = np.array(position)

    # print("Place called! pos = " + str(position) + ", rot = " + str(rotation))

      
    # rotation = self.randRotationMatrix()
      
    # We subtract soma before rotating to centre neuron
    if(rotation is not None):

      self.rotatedFlag = True
      
      if(len(self.axon) > 0):
        self.axon[:,0:3] = \
          np.transpose(np.matmul(rotation, \
                                 np.transpose(self.axon[:,0:3] \
                                              - self.soma[0,0:3])))        
      if(len(self.dend) > 0):
        self.dend[:,0:3] = \
          np.transpose(np.matmul(rotation, \
                                 np.transpose(self.dend[:,0:3] \
                                              - self.soma[0,0:3])))
      if(len(self.soma) > 0):
        self.soma[:,0:3] = \
          np.transpose(np.matmul(rotation, \
                                 np.transpose(self.soma[:,0:3] \
                                              - self.soma[0,0:3])))

        

    # Place neuron in correct position
    if(len(self.axon) > 0):
      self.axon[:,0:3] = self.axon[:,0:3] - self.soma[0,0:3] + position

    if(len(self.dend) > 0):
      self.dend[:,0:3] = self.dend[:,0:3] - self.soma[0,0:3] + position

    if(len(self.soma) > 0):
      self.soma[:,0:3] = self.soma[:,0:3] - self.soma[0,0:3] + position
    
    # Track rotation and location
    self.rotation = rotation
    self.position = position
    
    # Plot neuron post rotation
    if(False):
      self.plotNeuron()

    return self

  ############################################################################

  def saveCache(self, cacheFile=None):

    if(cacheFile is None):
      cacheFile = self.cache_filename

    assert not self.rotatedFlag, \
      "saveCache: The neuron should not be rotated when saving cache"
    
    morph = dict([])

    morph["swc_filename"] = self.swc_filename
    morph["soma"] = self.soma
    morph["axon"] = self.axon
    morph["dend"] = self.dend
    morph["axonLinks"] = self.axonLinks
    morph["dendLinks"] = self.dendLinks
    morph["dendSecX"] = self.dendSecX
    morph["dendSecID"] = self.dendSecID
    morph["axonStumpIDFlag"] = self.axonStumpIDFlag
    morph["maxAxonRadius"] = self.maxAxonRadius
    morph["maxDendRadius"] = self.maxDendRadius
    morph["dendDensity"] = self.dendDensity
    morph["axonDensity"] = self.axonDensity
    morph["version"] = self.cacheVersion
    
    assert(cacheFile != self.swc_filename)
    print("Saving cache file: " + cacheFile)
    
    import pickle
    with open(cacheFile,'wb') as cache_file:
      pickle.dump(morph, cache_file, self.pickleVersion)

  ############################################################################

  def cacheExist(self, cacheFile=None):

    if(cacheFile is None):
      cacheFile = self.cache_filename
    
    cacheFlag = False
    
    import os

    if(os.path.isfile(cacheFile)):
      
      swcTime = os.path.getmtime(self.swc_filename) 
      cacheTime = os.path.getmtime(cacheFile) 

      if(cacheTime > swcTime):
        print("Found cache file: " + cacheFile)
        cacheFlag = True
      else:
        print("Found old cache file: " + cacheFile)

    else:
      print("No cache file found.")
        
    return cacheFlag
  
  ############################################################################
    
  def loadCache(self, cacheFile=None):

    if(cacheFile is None):
      cacheFile = self.cache_filename
    
    import pickle
    with open(cacheFile,'rb') as cache_file:
      morph = pickle.load(cache_file)
    
    assert(self.swc_filename == morph["swc_filename"])
    assert self.axonStumpIDFlag == morph["axonStumpIDFlag"], \
    "axonStumpIDFlag must match cached version"
    
    # axonStumpIDFlag affects the section ID for the dendrites (and axon)
    # True when running Network_simulate.py and False if running Neurodamus.
    
    # self.axonStumpIDFlag = morph["axonStumpIDFlag"] # True or False
    
    self.soma = np.copy(morph["soma"])
    self.axon = np.copy(morph["axon"])
    self.dend = np.copy(morph["dend"])

    self.axonLinks = morph["axonLinks"]
    self.dendLinks = morph["dendLinks"]
    self.dendSecX = morph["dendSecX"]
    self.dendSecID = morph["dendSecID"]

    assert morph["version"] == self.cacheVersion, \
      "Cache version mismatch, regenerating cache"

    
    self.maxAxonRadius = morph["maxAxonRadius"]
    self.maxDendRadius = morph["maxDendRadius"]

    if(morph["dendDensity"] is not None):
      self.dendDensity = morph["dendDensity"]
    else:
      self.dendDensity = None
      
    if(morph["axonDensity"] is not None):
      self.axonDensity = morph["axonDensity"]
    else:
      self.axonDensity = None
      
    # Place neuron -- Do not place neuron, loadNeuronMorphology does that
    # self.place()

      
  ############################################################################

  # self.actionStumpIDFlag only affects the section ID.
  # If it is set to False, all sectionID are computed normally
  # if it is set to True, each axon will have the same sectionID throughout
  # if there are multiple axons they will have separate sectionIDs
  
  def loadSWC(self,swcFile):
    
    with open(swcFile,'r') as f:
      lines = f.readlines()

    comp_type = { 1: "soma", 2: "axon", 3: "dend", 4: "apic" }

    swcVals = np.zeros(shape=(len(lines),7))
   
    nComps = 0
    for ss in lines:
      if(ss[0] != '#'):
        swcVals[nComps,:] = [float(s) for s in ss.split()]
        nComps = nComps + 1
        
    # swcVals -- 0: compID, 1: type, 2,3,4: xyz coords, 5: radius, 6: parentID
    assert (1 <= swcVals[:nComps,1]).all() \
      and (swcVals[:nComps,1] <= 4).all(), \
      "loadMorphology: Only types 1,2,3,4 are supported: " + str(swcFile)

    # Subtract 1 from ID and parentID, so we get easier indexing
    swcVals[:,0] -= 1
    swcVals[:,6] -= 1

    swcVals[:,2:6] *= 1e-6 # Convert to meter x,y,z, radie

    # Columns:
    # 0: ID, 1,2,3: x,y,z 4: radie, 5: type, 6: parent, 7: somaDist,
    # 8: nodeParent, 9: childCount, 10: sectionID, 11: sectionLen,
    # 12: segmentLen
    
    # -- careful with sectionID and sectionLen at branch points,
    #    they belong to the parent section
    # -- also dont confuse sectionLen and segmentLen (the latter is
    #    for the segment, which is a part of the larger section)
    points = np.zeros((nComps,13))
    points[:nComps,0] = swcVals[:nComps,0]     # ID
    points[:nComps,1:5] = swcVals[:nComps,2:6] # x,y,z,r    
    points[:nComps,5] = swcVals[:nComps,1]     # type
    points[:nComps,6] = swcVals[:nComps,6]     # parent

    
    assert points[0,5] == 1, \
      "First compartment must be a soma: " + str(swcFile)
    
    # Create list of the links,
    # exclude soma -> first comp link (should be within soma radius)
    # Columns: 0: ID1, 1: ID2, 2: sectionID, 3: sectionX0, 4: sectionX1
    # 5: nodeParent, 6:type

    links = np.zeros((nComps,7))
    
    linkIdx = 0
    for idx in range(0,nComps):
      ID0 = int(points[idx,6]) # parent
      ID1 = int(points[idx,0]) # point

      if(ID0 <= 0):
        # No parent or soma is parent, skip link
        continue
      
      links[linkIdx,0:2] = [ID0,ID1]
      links[linkIdx,5] = points[idx,5]

      linkIdx += 1

    # Trim link list
    links = links[:linkIdx,:]

    # Count children each node has
    for idx in range(1,nComps):
      try:
        # Increment parents child counter
        points[int(points[idx,6]),9] += 1
      except:
        print("Are there gaps in the numbering of the compartments in the SWC file: " + str(swcFile))
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)
        import pdb
        pdb.set_trace()
      
    # Make sure soma has a child count > 1 --- no we dont want some as node 
    # if(points[0,9] == 0):
    #   points[0,9] = 100 

    # Also make sure all points with soma as parent get child count > 1
    # (Child Count > 1 ==> start or end of segment)
    somaChildIdx = np.where(points[:,6] == 0)[0]
    points[somaChildIdx,9] += 50
    
    # Mark node parent, and assign sectionID
    # -- this is used to set sectionID for links, the link
    # will use the end points sectionID
    # !!! Make sure sectionID is correct, and match what Neuron uses internally

    # Nodes are branch points (> 1 child), or end points (0 children)
    # and not soma
    nodeIdx = np.where((points[:,9] != 1) & (points[:,5] != 1))[0]    

    # soma is section 0, but we dont include connection soma to first node
    # so let the first dend node be 0, since the section ID is taken from
    # the child ID
    sectionID=1

    # Sonata specifies first axon, then basal, then apical sections
    axonIdx = nodeIdx[np.where(points[nodeIdx,5] == 2)[0]]
    basalIdx = nodeIdx[np.where(points[nodeIdx,5] == 3)[0]]
    apicalIdx = nodeIdx[np.where(points[nodeIdx,5] == 4)[0]]    

    # If simulation will use an axon stump, where each axon branch is shortened
    # to a stump with the same section ID, then we need to make sure the
    # numbering is correct for the dendrites.

    # Update, set axonID to -1
    for nIdx in axonIdx:
      points[nIdx,10] = -1

    # Set soma ID to 0
    points[0,10] = 0

    # Calculate sectionID for dendrites
    sectionID = 1

    # Axon dealt with, only loop over dendrites next
    nodeLoopList = [basalIdx,apicalIdx]
        
    for idxList in nodeLoopList:
      for nIdx in idxList:
        if(points[nIdx,6] > 0):
          # Set section ID, exclude soma, and compartments bordering to soma
          points[nIdx,10] = sectionID
          sectionID += 1
          
    # Assign node parents
    for nIdx in nodeIdx:

      # Find node parent
      parentIdx = int(points[nIdx,6])
      # While one child (= no node), keep following parent
      # But stop if parent is soma, or if grandparent is soma
      # !!! Here last link node to soma is not included in neurite morphology
      #     since we assume it is inside the soma
      while(points[parentIdx,9] == 1 \
            and parentIdx > 0 \
            and points[parentIdx,6] > 0):
        parentIdx = int(points[parentIdx,6])

      nodeParentIdx = parentIdx
      points[nIdx,8] = nodeParentIdx

      sectionID = points[nIdx,10]
      parentIdx = int(points[nIdx,6])
      while(points[parentIdx,9] == 1 and parentIdx > 0):
        points[parentIdx,8] = nodeParentIdx
        assert points[parentIdx,10] == 0, "SectionID should be unset prior"
        points[parentIdx,10] = sectionID
        parentIdx = int(points[parentIdx,6])
        
    for idx in range(1,nComps):      
      parentIdx = int(points[idx,6])
      
      # Calculate soma dist (and also save segLen)
      segLen = np.sqrt(np.sum((points[idx,1:4] - points[parentIdx,1:4])**2))
      points[idx,7] = points[parentIdx,7] + segLen
      points[idx,12] = segLen

      
    # Calculate section length (length between nodes)
    for idx in nodeIdx:
      nodeParentIdx = int(points[idx,8])
        
      # Difference in soma distance is section length
      sectionLen = points[idx,7] - points[nodeParentIdx,7]
      points[idx,11] = sectionLen

      if(sectionLen == 0):
        self.writeLog("Section length is zero --- !!! ")
        import pdb
        pdb.set_trace()
      
      prevIdx = int(points[idx,6])
      while(prevIdx > nodeParentIdx): 
        points[prevIdx,11] = sectionLen
        prevIdx = int(points[prevIdx,6])
        
    # Calculate sectionX
    for idx in range(0,links.shape[0]):
      ID0 = int(links[idx,0])
      ID1 = int(links[idx,1])
      links[idx,2] = points[ID1,10] # Section ID from point (not parent)

      nodeParent = int(points[ID1,8])
      nodeParentSomaDist = points[nodeParent,7]
      sectionLen = points[ID1,11]

      # segX0 and segX1
      links[idx,3] = (points[ID0,7] - nodeParentSomaDist)/sectionLen
      links[idx,4] = (points[ID1,7] - nodeParentSomaDist)/sectionLen

      links[idx,5] = nodeParent
      links[idx,6] = points[ID0,5] # type (use parent,
                                   # to avoid soma to dend link)
      
    # Store the soma, axon, dend and links in the object

    self.soma = np.zeros((1,4))
    self.soma[0,:] = swcVals[0,2:6] # save x,y,z,r

    dendIdx = np.where((points[:,5] == 3) | (points[:,5] == 4))[0]
    axonIdx = np.where(points[:,5] == 2)[0]

    dendLinkIdx = np.where((links[:,6] == 3) | (links[:,6] == 4))[0]
    axonLinkIdx = np.where(links[:,6] == 2)[0]

    # 0,1,2: x,y,z  3: radie, 4: dist to soma
    self.dend = np.zeros((len(dendIdx),5))
    self.axon = np.zeros((len(axonIdx),5))

    self.dendLinks = np.zeros((len(dendLinkIdx),2),dtype=int) # ID0,ID1
    self.axonLinks = np.zeros((len(axonLinkIdx),2),dtype=int) # ID0,ID1

    self.dendSecID = np.zeros((len(dendLinkIdx),),dtype=int) # SectionID
    self.dendSecX  = np.zeros((len(dendLinkIdx),2)) # SecX0, SecX1
    
    dendLookup = dict([])
    axonLookup = dict([])

    for idx in range(0,len(dendIdx)):
      dendLookup[dendIdx[idx]] = idx

    for idx in range(0,len(axonIdx)):
      axonLookup[axonIdx[idx]] = idx

    for idx,dIdx in enumerate(dendIdx):
      self.dend[idx,0:4] = points[dIdx,1:5] # x,y,z,r
      self.dend[idx,4] = points[dIdx,7] # dist to soma

    for idx,aIdx in enumerate(axonIdx):
      self.axon[idx,0:4] = points[aIdx,1:5] # x,y,z,r
      self.axon[idx,4] = points[aIdx,7] # dist to soma

    for idx,dIdx in enumerate(dendLinkIdx):
      self.dendLinks[idx,0] = dendLookup[int(links[dIdx,0])] # ID0 - parent
      self.dendLinks[idx,1] = dendLookup[int(links[dIdx,1])] # ID1

      self.dendSecID[idx] = links[dIdx,2]
      self.dendSecX[idx,:] = links[dIdx,3:5] 
      
    for idx,aIdx in enumerate(axonLinkIdx):
      self.axonLinks[idx,0] = axonLookup[links[aIdx,0]]
      self.axonLinks[idx,1] = axonLookup[links[aIdx,1]]
      # We alsoe have sectionID, secX0 and secX1 saved in links[:,2:5]
      # if needed in the future
      
    if(False):
      print("Inspect self.dend and axon")
      import pdb
      pdb.set_trace()

    if(self.virtualNeuron):
      # For virtual neurons, skip the dendrites (save space)
      self.dend = np.zeros((0,self.dend.shape[1]))
      self.dendLinks = np.zeros((0,2))
      self.dendSecID = np.zeros((0,))
      self.dendSecX = np.zeros((0,2))
      
    # self.dendriteDensity() # -- depricated
    self.findRadius()
    self.place()

    if(False):
      print("Debug plot")
      self.debugPlot()
      import pdb
      pdb.set_trace()

  ############################################################################

  def findRadius(self):

    if(len(self.axon) > 0):
      self.maxAxonRadius = \
        np.max(np.linalg.norm(self.axon[:,0:3]-self.soma[0,0:3], axis=1))

    if(len(self.dend) > 0):
      self.maxDendRadius = \
        np.max(np.linalg.norm(self.dend[:,0:3]-self.soma[0,0:3], axis=1))
    else:
      self.maxDendRadius = 0

    if(self.verbose):
      print("Max axon radius = " + str(self.maxAxonRadius))
      print("Max dend radius = " + str(self.maxDendRadius))
      
  ############################################################################

  def plotNeuron(self,axis=None,plotAxon=True,plotDendrite=True,lineStyle='-',alpha=1.0,plotOrigo=np.array([0,0,0]),plotScale=1.0):
    
    if(self.verbose):
      print("Plotting neuron " + self.swc_filename)
      
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    if(axis is None):
      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')
    else:
      ax = axis

    if(len(self.axon) > 0 and plotAxon):
      axLinks = []

      for row in self.axonLinks[:,:2].astype(int):
        if(len(axLinks) == 0):
          axLinks = list(row)
        elif(row[0] == axLinks[-1]):
          axLinks.append(row[1])
        elif(row[1] == axLinks[-1]):
          axLinks.append(row[0])
        else:
          ax.plot((self.axon[axLinks,0]-plotOrigo[0])*plotScale,
                  (self.axon[axLinks,1]-plotOrigo[1])*plotScale,
                  (self.axon[axLinks,2]-plotOrigo[2])*plotScale,
                  linestyle=lineStyle,
                  marker=',',
                  alpha=alpha,
                  c=self.colour)

          axLinks = list(row)

      if(len(axLinks) > 0):
          ax.plot((self.axon[axLinks,0]-plotOrigo[0])*plotScale,
                  (self.axon[axLinks,1]-plotOrigo[1])*plotScale,
                  (self.axon[axLinks,2]-plotOrigo[2])*plotScale,
                  linestyle=lineStyle,
                  marker=',',
                  alpha=alpha,
                  c=self.colour)

    if(plotDendrite):
      dendLinks = []
      for row in self.dendLinks[:,:2].astype(int):
        if(len(dendLinks) == 0):
          dendLinks = list(row)
        elif(row[0] == dendLinks[-1]):
          dendLinks.append(row[1])
        elif(row[1] == dendLinks[-1]):
          dendLinks.append(row[0])
        else:
          ax.plot((self.dend[dendLinks,0]-plotOrigo[0])*plotScale,
                  (self.dend[dendLinks,1]-plotOrigo[1])*plotScale,
                  (self.dend[dendLinks,2]-plotOrigo[2])*plotScale,
                  linestyle=lineStyle,
                  marker=',',
                  alpha=alpha,
                  c=self.colour)

          dendLinks = list(row)

      if(len(dendLinks) > 0):
        ax.plot((self.dend[dendLinks,0]-plotOrigo[0])*plotScale,
                (self.dend[dendLinks,1]-plotOrigo[1])*plotScale,
                (self.dend[dendLinks,2]-plotOrigo[2])*plotScale,
                linestyle=lineStyle,
                marker=',',
                alpha=alpha,
                c=self.colour)
        
          
    if(len(self.soma) > 0):
      ax.scatter((self.soma[:,0]-plotOrigo[0])*plotScale,
                 (self.soma[:,1]-plotOrigo[1])*plotScale,
                 (self.soma[:,2]-plotOrigo[2])*plotScale,
                 c=self.colour,alpha=alpha)
      
    plt.axis('equal')
    plt.ion()
    plt.show()
    plt.draw()
    plt.pause(0.001)

    return ax

  ############################################################################

  # !!! Is this function depricated?
  
  def dendriteDensity(self):

    assert False, "Depricated function dendriteDensity?"
    
    if(len(self.dend) == 0):
      assert self.virtualNeuron, \
        "No dendrites in " + str(self.name) \
        + ". Only virtual neurons are allowed to have no dendrites!"
 
      # No dendrites, neuron is virtual
      return None
    
    if(True or self.dendDensity is None):
      
      # Calculate all the segment lengths
      dendSegmentLength = \
        np.sum(((self.dend[self.dendLinks[:,0].astype(int),:][:,0:3]
                 -self.dend[self.dendLinks[:,1].astype(int),:][:,0:3])**2),
               axis=-1) ** 0.5     
      
      # Calculate all segment centres distances to soma
      dendSegmentDist = \
        np.sum((((self.dend[self.dendLinks[:,0].astype(int),:][:,0:3]
                  +self.dend[self.dendLinks[:,1].astype(int),:][:,0:3])/2
                 - self.soma[0,0:3])**2),
               axis=-1) ** 0.5

      binSize = self.densityBinSize
      maxDist = np.max(dendSegmentDist)
      nBins = int(np.ceil(maxDist/binSize))+1
      
      self.dendDensity = np.zeros((nBins,1))
      
      for sd,sl in zip(dendSegmentDist,dendSegmentLength):
        idx = int(np.floor(sd/binSize))
        self.dendDensity[idx] += sl

      # Divide by volume of shell to get density
      for i in range(0,nBins):
        self.dendDensity[i] /= 4*np.pi/3*(((i+1)*binSize)**3 - (i*binSize)**3)
                      
    return (self.dendDensity,self.densityBinSize)

  ############################################################################

  def setAxonVoxelRadialDensity(self,density,maxAxonRadius):

    print("Only saving equation now")

    self.axonDensityType = "r"
    self.axonDensity = density
    self.maxAxonRadius = maxAxonRadius

  ############################################################################

  def setAxonVoxelXYZDensity(self,
                             density,
                             axonDensityBoundsXYZ):

    print("Only saving equation now")

    self.axonDensityType = "xyz"
    self.axonDensity = density
    self.axonDensityBoundsXYZ = axonDensityBoundsXYZ
      
  ############################################################################

  def compartmentLength(self,compType="dend"):
    if(compType == "dend"):
      links = self.dendLinks
      coords = self.dend
    elif(compType == "axon"):
      links = self.axonLinks
      coords = self.axon
    else:
      assert False, "Unknown compartment type: " + str(compType) \
        + ", valid types are 'axon' and 'dend'"
    
    compLen = np.linalg.norm(coords[links[:,0],:][:,:3] \
                             - coords[links[:,1],:][:,:3], \
                             axis=1)

    return compLen
    
  ############################################################################

  # !!! Add plot functon for density

  # Now density is specified as an equation of d, this func needs to be rewritten
  
  def plotDensityOLD(self):

    d = [self.densityBinSize*x*1e6 for x in range(0,len(self.dendDensity))]
    
    import matplotlib.pyplot as plt
    plt.figure()
    plt.step(d,self.dendDensity*1e-12)
    plt.xlabel('Distance from soma (mum)')
    plt.ylabel('Density (mum/mum3)')

    if(self.axonDensity is not None):
      da = [self.densityBinSize*x*1e6 for x in range(0,len(self.axonDensity))]
      plt.figure()
      plt.step(da,self.axonDensity*1e-12)
      plt.xlabel('Distance from soma (mum)')
      plt.ylabel('Axon density (mum/mum3)')
     
    
  ############################################################################

  def debugPlot(self,waitFlag=True,plotStep=1,plotAxonFlag=False):

    ax = self.plotNeuron(plotAxon=plotAxonFlag)

    if(plotAxonFlag):
      for a in self.axonLinks:
        x0 = self.axon[int(a[0]),0:3]
        x1 = self.axon[int(a[1]),0:3]
        x = (x0 + x1) / 2

        #ax.text(x=x0[0],y=x0[1],z=x0[2],s=str(np.around(a[3],2)),color='blue')
        #ax.text(x=x1[0],y=x1[1],z=x1[2],s=str(np.around(a[4],2)),color='red')      
        ax.text(x=x[0],y=x[1],z=x[2],s=str(a[2]),color='black')

        print("ID: " + str(a[2]))      
        input(" ")

    ctr = 0
    for (d,dID,dX) in zip(self.dendLinks,self.dendSecID,self.dendSecX):
      x0 = self.dend[int(d[0]),0:3]
      x1 = self.dend[int(d[1]),0:3]
      x = (x0 + x1) / 2

      #ax.text(x=x0[0],y=x0[1],z=x0[2],s=str(np.around(dX[0],2)),color='blue')
      #ax.text(x=x1[0],y=x1[1],z=x1[2],s=str(np.around(dX[1],2)),color='red')

      if(ctr % plotStep == 0):
        ax.text(x=x[0],y=x[1],z=x[2],s=str(dID),color='black')
      ctr += 1
        
      print("ID: " + str(dID) + " X = " + str(np.around(dX[0],2)) + " - " \
            + str(np.around(dX[1],2)))

      if(waitFlag):
        input(" ")

    return ax
      
  ############################################################################

if __name__ == "__main__":
  # The lines below are just for testing purposes
  
  # nm = NeuronMorphology(swc_filename='morphology/network/FSN-BE79B-3ak-compact-5.swc',verbose=True,clusterFlag=True,nClustersDend=-1,nClustersAxon=-1)

  # fName = "cellspecs/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20190508/WT-0728MSN01-cor-rep-ax.swc"
  fName = "cellspecs/lts/LTS_Experiment-9862_20181211/Experiment-9862corrected-cor-rep.swc"
  
  nm = NeuronMorphology(swc_filename=fName,verbose=True,useCache=False)
  
  nm.place(rotation=nm.randRotationMatrix(),position=np.array([0,0,0]))

  nm.debugPlot()
  

  # nm.setAxonDensity("3e9*np.exp(-d/100e-6)",300e-6)
  
  # nm.plotDensity()
  
  ax1 = nm.plotNeuron()

  print("In main function")
  import pdb
  pdb.set_trace()
  
  nm2 = nm.clone(rotation=nm.randRotationMatrix(),position=np.array([0.001,0.001,0.001]))
  nm2.plotNeuron(ax1)

  nm2.plotDensity()

  
  # raw_input("Test")

  import pdb
  pdb.set_trace()

    
  
