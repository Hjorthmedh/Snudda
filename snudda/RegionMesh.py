import numpy as np
import scipy
from scipy import ndimage
import re
import os
import pickle
import timeit


class RegionMesh(object):

  ############################################################################
  
  def __init__(self,fileName,dView = None, lbView = None, role="master",
               useCache=True,pickleVersion=-1,raytraceBorders=True,
               dMin=15e-6,binWidth=1e-4,logFileName=None,logFile=None):

    self.dView = dView
    self.lbView = lbView
    
    self.role = role
    self.workersInitialised = False

    self.verbose = True
    
    if(logFile is not None):
      self.logFile = logFile
      self.logFileName = logFile.name
    elif(logFileName is not None and len(logFileName) > 0):
      self.logFile = open(logFileName,'wt')
      self.logFileName = logFileName
    else:
      self.logFile = None
      self.logFileName = None

    #self.binWidth = 5e-4
    self.binWidth = binWidth # 1e-4 # 5e-5 # 1e-4    
    self.padding = max(self.binWidth,dMin)
    
    # This determines if we ray trace the border voxels, for finer detail
    # or not (activating this is SLOW)
    self.raytraceBorders = raytraceBorders
    
    if(raytraceBorders):
      self.writeLog("Ray tracing points in border voxels, this is slow.")
      rtStr = "-RTB"
    else:
      rtStr = ""
    
    # binWidth 5e-4 (94s) --> 10.8 % border voxels
    # binWidth 2.5e-4 (577) --> 6.2 % border voxels
    # binWidth 1e-4 (8090s) --> 2.7 % border voxels
    # binWidth 0.5e-4 (??? s) --> ?? % border voxels

    self.fileName = fileName
    self.cacheFile = fileName + "-" + str(int(1e6*self.binWidth)) + rtStr + "-cache.pickle"
    self.pickleVersion = pickleVersion
    
    self.debugFlag = False

    cacheLoaded = False
    
    if(useCache and self.cacheExist()):
      try:
        self.loadCache()
        cacheLoaded = True
      except:
        self.writeLog("Failed to load cache.")
        cacheLoaded = False

    if(not cacheLoaded):
      self.loadMesh(fileName=fileName)

      tA = timeit.default_timer()
      
      self.preCompute()
      self.setupVoxelFilter()
      
      tB = timeit.default_timer()

      self.writeLog("Calculation time: " + str(tB-tA) + " s")

    self.setupVoxelList()

    self.setupPlaceNeurons(dMin=dMin)

    self.writeLog("Inner voxel bin volume: " \
                  + str(np.round(self.innerVoxelVolume()*1e9,1)) + " mm³")
    # !!! RUN THIS
    # self.verifyInside(100)


  ############################################################################

  def markBorders(self):

    self.writeLog("Marking borders")
    
    for row in self.meshVec:
      ic = np.floor((row - self.minCoord)/self.binWidth)
      self.voxelMaskBorder[int(ic[0]),int(ic[1]),int(ic[2])] = 1
          
    maxN = 0
    
    for row in self.meshFaces:
      coord1 = self.meshVec[row[0],:]
      coord2 = self.meshVec[row[1],:]
      coord3 = self.meshVec[row[2],:]

      vx = coord2-coord1
      vy = coord3-coord1
      
      dx = np.linalg.norm(coord2-coord1)
      dy = np.linalg.norm(coord3-coord1)

      nx = int(2*np.ceil(dx/self.binWidth))+1
      ny = int(2*np.ceil(dy/self.binWidth))+1

      maxN = max(maxN,max(nx,ny))
      
      for xStep in np.linspace(0,1,nx):
        for yStep in np.linspace(0,1-xStep,ny):
          xPoint = coord1 + vx*xStep + vy*yStep
          ic = np.floor((xPoint - self.minCoord)/self.binWidth)
          self.voxelMaskBorder[int(ic[0]),int(ic[1]),int(ic[2])] = 1
          
    self.writeLog("maxN = " + str(maxN))
    

  ############################################################################

  def setupParallel(self,dView=None):

    if(dView is None):
      
      if(self.dView is None):
        self.writeLog("Running in serial")
        return
      
      dView = self.dView

    self.writeLog("Setting up parallel")
    
    assert self.role == "master", \
      "setupParallel should only be called by master"

    if(self.workersInitialised):
      self.writeLog("setupParallel: workers already initialised")
      return

    with dView.sync_imports():
      from RegionMesh import RegionMesh
    
    try:
      self.writeLog("Setting up RegionMesh on workers")
      dView.push({"fileName":self.fileName}, block=True)

      if(self.logFileName is not None):
        engineLogFile = [self.logFileName + "-" \
                         + str(x) for x in range(0,len(dView))]
      else:
        engineLogFile = [[] for x in range(0,len(dView))]
        
      dView.scatter('logFileName',engineLogFile,block=True)

      cmdStr = "sm = RegionMesh(fileName=fileName,role='worker',logFileName=logFileName[0])"
      
      dView.execute(cmdStr,block=True)
      self.writeLog("Worker RegionMesh setup done.")

    except Exception as e:
      import uuid
      import traceback
      tstr = traceback.format_exc()
      tmp = open("save/tmp-striatum-mesh-log-file-" + str(uuid.uuid4()),'w')
      tmp.write("Exception: " + str(e))
      tmp.write("Trace:" + tstr)
      tmp.close()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()
      

    self.workersInitialised = True
    
  ############################################################################

  def cacheExist(self):

    cacheFlag = False

    if(os.path.isfile(self.cacheFile)):

      objTime = os.path.getmtime(self.fileName)
      cacheTime = os.path.getmtime(self.cacheFile)

      if(cacheTime > objTime):
        self.writeLog("Found cache file " + self.cacheFile)
        cacheFlag = True
      else:
        self.writeLog("Found old cache file (" + str(self.cacheFile) + "), ignoring.")
    else:
      self.writeLog("No mesh cache file found (" + str(self.cacheFile) + ")")

    return cacheFlag

  ############################################################################

  def loadCache(self):

    self.writeLog("Loading cache file: " + self.cacheFile)
    
    with open(self.cacheFile,'rb') as f:
      data = pickle.load(f)

    self.meshVec   = data["meshVec"]
    self.meshFaces = data["meshFaces"]
    self.meshNorm  = data["meshNorm"]
    self.minCoord  = data["minCoord"]
    self.maxCoord  = data["maxCoord"]
    self.voxelMaskInner  = data["voxelMaskInner"]
    self.voxelMaskBorder = data["voxelMaskBorder"]

    self.pointOut = self.maxCoord + 1e-2
    
    assert self.binWidth == data["binWidth"], \
      "Mismatch binWidth: " + str(self.binWidth) +" vs " + str(data["binWidth"])
    assert self.padding == data["padding"], \
      "Missmatch padding: " + str(self.padding) + " vs " + str(data["padding"])

    assert self.raytraceBorders == data["raytraceBorders"], \
      "Missmatch raytraceBorders: " + str(self.raytraceBorders) \
      + " vs " + str(data["raytraceBorders"])
    
    self.nBins = self.voxelMaskInner.shape
    self.preCompute()
    
  ############################################################################

  def saveCache(self):

    if(self.role != "master"):
      return
    
    data = dict([])
    data["meshVec"]   = self.meshVec
    data["meshFaces"] = self.meshFaces
    data["meshNorm"]  = self.meshNorm
    data["minCoord"]  = self.minCoord
    data["maxCoord"]  = self.maxCoord
    data["voxelMaskInner"]  = self.voxelMaskInner
    data["voxelMaskBorder"] = self.voxelMaskBorder      
    data["padding"] = self.padding
    data["binWidth"] = self.binWidth
    data["raytraceBorders"] = self.raytraceBorders
    
    self.writeLog("Saving mesh cache file " + self.cacheFile)
    with open(self.cacheFile,'wb') as f:
      pickle.dump(data,f,self.pickleVersion)
    
  ############################################################################
  
  def loadMesh(self,fileName):

    self.fileName = fileName
    
    allVec = []
    allFaces = []
    allNorm = []

    # https://stackoverflow.com/questions/4703390/how-to-extract-a-floating-number-from-a-string
    numeric_const_pattern = r"""
    [-+]? # optional sign
    (?:
    (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
    |
    (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
    )
    # followed by optional exponent part if desired
    (?: [Ee] [+-]? \d+ ) ?
    """
    rx = re.compile(numeric_const_pattern, re.VERBOSE)
    
    with open(fileName,'rt') as f:
      for row in f:
        if(row[0:2] == 'v '):
          digits = rx.findall(row)
          # Convert to SI units
          allVec.append([float(d)*1e-6 for d in digits])
          
        if(row[0:2] == 'f '):
          # Only take first value of each triplet 1/?/? 2/?/? 3/?/?
          digits = re.findall(r'f\s+(\d+)/\d*/\d*\s+(\d+)//\d*\s+(\d+)//\d*',
                              row)
          # Subtract one, to get python indexing
          try:
            allFaces.append([int(d)-1 for d in digits[0]])
          except:
            self.writeLog("Problem with reading digits")
            self.writeLog(row + "\nread: " + str(digits))
            import pdb
            pdb.set_trace()

        if(row[0:2] == 'vn'):
          digits = rx.findall(row)
          allNorm.append([float(d) for d in digits])
          
      self.meshVec = np.zeros((len(allVec),3))
      self.meshFaces = np.zeros((len(allFaces),3),dtype=int)
      self.meshNorm = np.zeros((len(allNorm),3))      

      for ir,row in enumerate(allVec):
        self.meshVec[ir,:] = row
      
      for ir,row in enumerate(allFaces):
        self.meshFaces[ir,] = row
        
      for ir,row in enumerate(allNorm):
        self.meshNorm[ir,:] = row

      try:
        self.minCoord = np.min(self.meshVec,axis=0)
        self.maxCoord = np.max(self.meshVec,axis=0)
      except:
        self.writeLog("Problem calculating min and max coords")
        import pdb
        pdb.set_trace()

      # Used by ray casting when checking if another point is interior
      self.pointOut = self.maxCoord + 1e-2
      
  ############################################################################

  def preCompute(self):

    i0 = self.meshFaces[:,0]
    i1 = self.meshFaces[:,1]    
    i2 = self.meshFaces[:,2]
    
    self.u = self.meshVec[i1,:] - self.meshVec[i0,:]
    self.v = self.meshVec[i2,:] - self.meshVec[i0,:]

    self.v0 = self.meshVec[i0,:]

    self.uv = np.sum(np.multiply(self.u,self.v),axis=1)
    self.vv = np.sum(np.multiply(self.v,self.v),axis=1)
    self.uu = np.sum(np.multiply(self.u,self.u),axis=1)
    
    self.denom = np.multiply(self.uv,self.uv) - np.multiply(self.uu,self.vv)

    # Normal of triangle
    self.nrm = np.cross(self.u,self.v)

    # We need to normalise it
    nl = np.repeat(np.reshape(self.uv,[self.uv.shape[0],1]),3,axis=1)
    self.nrm = np.divide(self.nrm,nl)
    
    
  ############################################################################

  def setupVoxelFilter(self):

    
    if(self.role == "master"):
      self.setupParallel()
    
    self.minCoord = np.floor((np.min(self.meshVec,axis=0)
                              - self.padding)/self.binWidth )*self.binWidth
    self.maxCoord = np.ceil((np.max(self.meshVec,axis=0)
                             + self.padding)/self.binWidth)*self.binWidth

    self.nBins = np.array(np.ceil((self.maxCoord-self.minCoord)/self.binWidth + 1),
                          dtype=int)


    self.writeLog("Voxel mask: " + str(self.nBins[0]) + "x" + str(self.nBins[1]) \
          + "x" + str(self.nBins[2]))

    self.voxelMaskInner = np.zeros(self.nBins,dtype=bool)
    self.voxelMaskBorder = np.zeros(self.nBins,dtype=bool)    

    # All voxels with a mesh point in them are "border voxels"
    # For all remaining points, check if inner or outer
    # If raytraceBorders is false, we dont raytace for border points
    # at run time, this gives a bit jagged edges, but is MUCH faster
    # when placing cells (no ray tracing then)
    
    if(self.raytraceBorders):
      self.markBorders()
      
    nBinsTotal = self.nBins[0]*self.nBins[1]*self.nBins[2]
    iterCtr = 0

    # This second part is only run by the master, it calls the workers
    # to perform part of the computation

    if(self.role == "master"):
      # This should only be done by master

      if(self.dView is None):

        # No workers, do all work ourselves
        # The worker function adds a dimension (so gather works in parallel
        # case), here we just need to reshape results.
        vmInner = self._voxelMaskHelper(range(0,self.nBins[0]))
        self.voxelMaskInner = np.reshape(vmInner,self.nBins)

        # !!! Old code, now using the _voxelMaskHelper function instead
        #
        #for ix in range(0,self.nBins[0]):
        #  for iy in range(0,self.nBins[1]):
        #    print(str(iterCtr) + "/" + str(nBinsTotal))
        #    for iz in range(0,self.nBins[2]):
        #      iterCtr += 1
        #  
        #      if(not self.voxelMaskBorder[ix,iy,iz]):
        #        # Inner or outer point, check centre
        #        xyz = np.array([self.minCoord[0] + (ix+0.5)*self.binWidth,
        #                        self.minCoord[1] + (iy+0.5)*self.binWidth,
        #                        self.minCoord[2] + (iz+0.5)*self.binWidth])
        #        
        #        self.voxelMaskInner[ix,iy,iz] = self.rayCasting(xyz)

      else:
        # Distribute the work to the workers
        # Randomize order, to spread work load a bit better
        allX = np.random.permutation(np.arange(0,self.nBins[0]))

        try:
          self.dView.scatter("xRange",allX,block=True)
          self.writeLog("Starting parallel job")
          self.dView.execute("innerMask = sm._voxelMaskHelper(xRange)",block=True)
          self.writeLog("Gathering results")
          innerMask = self.dView.gather("innerMask",block=True)
        except Exception as e:
          self.writeLog("Oh no, something failed with parallel meshbuidling")
          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)
          self.writeToRandomFile(tstr)
          import pdb
          pdb.set_trace()
          
        for m in innerMask:
          self.voxelMaskInner = np.logical_or(self.voxelMaskInner,m)

    self.writeLog("Fraction of border voxels: " \
          + str(np.sum(self.voxelMaskBorder)
                /np.prod(self.voxelMaskBorder.shape)))

    self.saveCache()


    if(np.sum(self.voxelMaskInner) == 0):
      self.writeLog("Warning no inner voxels in mesh, check your meshBinWidth")
      self.writeLog("mesh file: " + self.fileName)

  #  if(self.role == "master"):
  #    import pdb
  #    pdb.set_trace()

  ############################################################################

  def writeToRandomFile(self,text):
    
    import uuid
    tmp = open("save/regmesh-tmp-log-file-" + str(uuid.uuid4()),'w')
    tmp.write(text)
    tmp.close()
    print(text)

  
  ############################################################################

  def checkInside(self,coords):
    idx = np.array(np.floor((coords - self.minCoord)/self.binWidth),dtype=int)
    
    if(self.voxelMaskInner[idx[0],idx[1],idx[2]]):
      # We know it is an inner voxel
      return True
    elif(self.voxelMaskBorder[idx[0],idx[1],idx[2]]):
      # We are in a border voxel, need to ray cast this
      return self.rayCasting(coords)
    else:
      # We are outside structure
      return False

    
  ############################################################################

  def _voxelMaskHelper(self, xRange):

    try:

      # Need the extra dimension at the top to make gather work
      vmInner = np.zeros((1,self.nBins[0],self.nBins[1],self.nBins[2]),
                         dtype=bool)
    
      for ix in xRange:
        self.writeLog("Processing x = " + str(ix))
        
        for iy in range(0,self.nBins[1]):
          for iz in range(0,self.nBins[2]):
            if(not self.voxelMaskBorder[ix,iy,iz]):
              # Inner or outer point, check centre
              xyz = np.array([self.minCoord[0] + (ix+0.5)*self.binWidth,
                              self.minCoord[1] + (iy+0.5)*self.binWidth,
                              self.minCoord[2] + (iz+0.5)*self.binWidth])

              vmInner[0,ix,iy,iz] = self.rayCasting(xyz)

    except:

      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      self.writeToRandomFile(tstr)
      import pdb
      pdb.set_trace()
  
    return vmInner

  ############################################################################

  # Cast a ray, see how many times it intersects the triangles of the structure
  # if an odd number of intersections, then it is an inner point
  #
  # Based on: https://www.erikrotteveel.com/python/three-dimensional-ray-tracing-in-python/
  
  def rayCasting(self,point):

    nTri = self.meshFaces.shape[0]

    P = self.pointOut-point
    # rn = nominator, rd = denominator
    rn = np.sum(np.multiply(self.nrm,self.v0 - point),axis=1)
    rd = np.dot(self.nrm,P)
    
    #r = np.divide(rn,rd)
    
    intersectCount = 0

    for i in range(0,nTri):
      
      if(rd[i] == 0):
        if(rn[i] == 0):
          # Parallel and lies in the plane
          rI = 0.0
        else:
          # Parallel to plane, but outside. Mark by -1 to avoid counting
          rI = -1
      else:
        rI = rn[i] / rd[i]

      if(rI >= 0 and rI <= 1):
        # Crosses the plane, but is it within triangle?

        w = point + rI * P - self.v0[i,:]
        
        si = (self.uv[i] * np.inner(w,self.v[i,:]) \
              - self.vv[i] * np.inner(w,self.u[i,:])) / self.denom[i]
        
        if(si < 0 or si > 1):
          # outside of triangle
          continue

        ti = (self.uv[i] * np.inner(w,self.u[i,:]) \
              - self.uu[i] * np.inner(w,self.v[i,:]))/self.denom[i]

        if(ti < 0 or (si + ti) > 1):
          # outside of trinagle
          continue

        # print("intersects face i = " + str(i))
        # print("si = " +str(si) + ", ti = " + str(ti))
        
        intersectCount += 1

        if(self.debugFlag):
          import pdb
          pdb.set_trace()
        
        
    # print("intersection count = " + str(intersectCount))
    
    return np.mod(intersectCount,2) == 1

  ############################################################################

  def verifyInside(self,nPoints=1000):

    if(self.role != "master"):
      return
    
    xTest = np.random.uniform(self.minCoord[0],self.maxCoord[0],nPoints)
    yTest = np.random.uniform(self.minCoord[1],self.maxCoord[1],nPoints)
    zTest = np.random.uniform(self.minCoord[2],self.maxCoord[2],nPoints)

    self.writeLog("Verifying our method")
    
    for i in range(0,nPoints):

      self.writeLog(str(i) + "/" + str(nPoints))
      coords = np.array([xTest[i],yTest[i],zTest[i]])
      
      try:
        assert(self.checkInside(coords) == self.rayCasting(coords))
      except:
        self.writeLog("Mismatch for coords: " + str(coords))
        self.writeLog("Cached: " + str(self.checkInside(coords)))
        self.writeLog("RC: " + str(self.rayCasting(coords)))
        import pdb
        pdb.set_trace()


  ############################################################################

  def setupPlaceNeurons(self,dMin=None):

    self.writeLog("Setup place neurons")
    
    self.dMin = dMin
    
    self.maxRand = 1000000  
    self.randCtr = self.maxRand + 1

    self.randomPool = np.zeros((self.maxRand,3))

    self.maxNeurons = 3000000

    self.neuronCoords = np.zeros((self.maxNeurons,3))
    self.neuronCtr = 0

    self.paddingCoords = np.zeros((self.maxNeurons,3))
    self.paddingCtr = 0

    self.neuronTypes = dict([])
    self.neuronType = np.zeros((self.maxNeurons,))
    self.nextNeuronType = 1

    self.maxReject = 100e6
    self.rejectCtr = 0

    self.nBins = self.voxelMaskInner.shape

    self.updatePaddingMask()
    self.updateRandomPool()

    self.writeLog("Setup done")

  ############################################################################

  def setupVoxelList(self):

    self.writeLog("Setup voxel list")
    
    self.maxNeuronsVoxel = 10000

    self.voxelNextNeuron = np.zeros(self.nBins,dtype=int)
    self.voxelNeurons = np.zeros((self.nBins[0],self.nBins[1],
                                  self.nBins[2],self.maxNeuronsVoxel,3))

  ############################################################################

  def updatePaddingMask(self):

    self.writeLog("Update padding mask")
    
    # Only save padding for border voxels, and voxels nearby
    nDist = int(np.ceil(2*self.dMin/self.binWidth))
    s = ndimage.generate_binary_structure(3,3)
    
    if(np.sum(self.voxelMaskBorder)):
      # If we do raytracing then the border exists
      self.voxelMaskPadding = np.copy(self.voxelMaskBorder)

      self.voxelMaskPadding = \
            scipy.ndimage.binary_dilation(self.voxelMaskPadding,
                                          structure=s,
                                          iterations=nDist)
    else:
      # No ray tracing, we need to create a padding region
      
      dilatedMask = scipy.ndimage.binary_dilation(self.voxelMaskInner,
                                                  structure=s,
                                                  iterations=nDist)

      self.voxelMaskPadding = \
          np.logical_xor(dilatedMask,self.voxelMaskInner)

      
  ############################################################################

  def checkPaddingZone(self,coords):

    idx = np.array(np.floor((coords - self.minCoord)/self.binWidth),dtype=int)

    return self.voxelMaskPadding[idx[0],idx[1],idx[2]]
    
  ############################################################################

  def updateRandomPool(self):

    if(self.randCtr >= self.maxRand):

      self.writeLog("Regenerating new random pool")
      for i in range(0,3):
        self.randomPool[:,i] = np.random.uniform(low=self.minCoord[i],
                                                 high=self.maxCoord[i],
                                                 size=self.maxRand)

      self.randCtr = 0
    
  ############################################################################

  # neuronType is included since in future versions we might want varying
  # density depending on the type of neuron
  
  def placeNeurons(self, nCells, neuronType=None, dMin=None):

    if(dMin is None):
      dMin = self.dMin

    try:
      dMin2 = dMin ** 2
    except:
      self.writeLog("Make sure dMin is set. dMin = " + str(dMin))

    # If this is not fullfullid, then we need to update the range values below
    assert 2*dMin < self.binWidth, \
      "2*dMin (2 * " + str(dMin) + ") must be smaller than binWidth (" + str(self.binWidth) + ")"

    if(neuronType in self.neuronTypes):
      neuronID = self.neuronTypes[neuronType]
    else:
      neuronID = self.nextNeuronType
      self.nextNeuronType += 1
      
      self.neuronTypes[neuronType] = neuronID

    
    startCtr = self.neuronCtr
    endCtr = self.neuronCtr + nCells

    tA = timeit.default_timer()

    while(self.neuronCtr < endCtr \
          and self.rejectCtr < self.maxReject):

      putativeLoc = self.randomPool[self.randCtr,:]
      self.randCtr += 1

      if(self.randCtr % 100000 == 0):
        self.writeLog("Neurons: " + str(self.neuronCtr) \
              + " Rejected: " + str(self.rejectCtr) \
              + " Padding: " + str(self.paddingCtr))

        if(self.neuronCtr == 0):
          self.writeLog("No neurons placed, check why!")
          import pdb
          pdb.set_trace()
      
      if(self.randCtr >= self.maxRand):
        self.updateRandomPool()

      minDist2 = 1e6
      
      insideFlag = self.checkInside(coords=putativeLoc)
      
      if(not insideFlag and not self.checkPaddingZone(putativeLoc)):
        # We are outside, go to next neuron
        self.rejectCtr += 1
        continue

      
      # Check that we are not too close to existing points
      # Only check the neighbouring voxels, to speed things up
      voxelIdx = np.array(np.floor((putativeLoc - self.minCoord)
                                   / self.binWidth),dtype=int)

      voxelIdxList = []
      voxelIdxList.append(voxelIdx)

      borderVoxel = np.zeros((3,),dtype=int)
      
      for idx in range(0,3):
        if((putativeLoc[idx] - self.minCoord[idx]) % self.binWidth < dMin):
          borderVoxel[idx] = -1
          newIdx = np.copy(voxelIdx)
          newIdx[idx] -= 1
          voxelIdxList.append(newIdx)

        elif((putativeLoc[idx] - self.minCoord[idx]) % self.binWidth > self.binWidth - dMin):
          borderVoxel[idx] = 1
          newIdx = np.copy(voxelIdx)
          newIdx[idx] += 1
          voxelIdxList.append(newIdx)

      nBorder = np.sum(np.abs(borderVoxel))
      
      if(nBorder == 2):
        # Along one of the lines, need to check diagonal voxel also
        voxelIdxList.append(voxelIdx + borderVoxel)
      elif(nBorder == 3):
        # Close to corner, need to check 8 voxels in total (ouch!)
        voxelIdxList.append(voxelIdx + [borderVoxel[0],borderVoxel[1],0])
        voxelIdxList.append(voxelIdx + [borderVoxel[0],0,borderVoxel[2]])
        voxelIdxList.append(voxelIdx + [0,borderVoxel[1],borderVoxel[2]])
        voxelIdxList.append(voxelIdx + borderVoxel)
        
      minDist2 = 1e6
      for voxIdx in voxelIdxList:
        if((voxIdx<0).any() or (voxIdx > self.nBins).any()):
          # Voxel outside bounds, ignore
          continue

        if(minDist2 < dMin2):
          # No need to calculate further, we are too close
          break

        if(self.voxelNextNeuron[voxIdx[0],voxIdx[1],voxIdx[2]] > 0):
          tmp = self.voxelNeurons[voxIdx[0],voxIdx[1],voxIdx[2],
                                  0:self.voxelNextNeuron[voxIdx[0],
                                                         voxIdx[1],
                                                         voxIdx[2]],
                                  :] \
                                  - putativeLoc
        
            
          minDist2 = min(minDist2, np.min(np.sum(np.square(tmp), axis=1)))
          
      if(dMin2 < minDist2):
        # Ok neuron is not too close to any neighbours

        saveInVoxelList = False

        if(insideFlag):
          # We are inside, add to inside points
          self.neuronCoords[self.neuronCtr,:] = putativeLoc
          self.neuronType[self.neuronCtr] = neuronID
          self.neuronCtr += 1
          #self.writeLog("Placed neuron " + str(self.neuronCtr))
        else:
          self.paddingCtr += 1

        # Also save the point in the specific voxel, this way we can ignore
        # lots of distance comparisons
        self.voxelNeurons[voxelIdx[0],voxelIdx[1],voxelIdx[2],
                          self.voxelNextNeuron[voxelIdx[0],
                                               voxelIdx[1],
                                               voxelIdx[2]],
                          :] = putativeLoc
          
        self.voxelNextNeuron[voxelIdx[0],voxelIdx[1],voxelIdx[2]] += 1
          
      else:
        self.rejectCtr += 1
        
    tB = timeit.default_timer()
    self.writeLog("Placed " + str(nCells) + " in " + str(tB-tA) + " s")
  
        
    return self.neuronCoords[startCtr:endCtr,:]
        

      # Store point in global list, but also reference it in specific voxel
      
  
  ############################################################################

  # shape = "cube" or "sphere"
  # radius = radius of sphere, or half the length of side of the cube
  
  def getSubset(self,centre,radius=None,nNeurons=None,shape="cube",
                returnIdxFlag=False):

    assert ((radius is None) ^ (nNeurons is None)), \
      "Specify one of radius or nNeurons."

    coords = self.neuronCoords[:self.neuronCtr,:]
    nrnType = self.neuronType[:self.neuronCtr,:]
    
    if(shape == "cube"):
      dist = np.amax(abs(coords - centre),axis=1)
    elif(shape == "sphere"):
      dist = np.sqrt(np.sum(np.square(coords - centre),axis=1))
    else:
      assert False, "Unknown shape " + str(shape) + " use cube or sphere"


    sortedIdx = np.argsort(dist)

    if(radius is not None):
      self.writeLog("Using radius " + str(radius))
      idx = sortedIdx[np.where(dist[sortedIdx] <= radius)]
    else:
      self.writeLog("Selecting " + str(nNeurons) + " closest neuron")
      idx = sortedIdx[:nNeurons]

    # Next we need to return them in the order they were originally sorted
    keepMask = np.zeros((self.neuronCtr,),dtype=bool)
    keepMask[idx] = True

    if(returnIdxFlag):
      return np.where(keepMask)
    else:
      return (coords[keepMask,:],nrnType[keepMask,:])
  
  ############################################################################
  
  def plotStruct(self,pdfName=None):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(self.meshVec[:,0],
               self.meshVec[:,1],
               self.meshVec[:,2],
               'black')

    plt.ion()
    plt.show()
    plt.pause(0.001)

    ax.view_init(150,40); plt.draw()
    plt.axis("off")
    
    if(pdfName is not None):
      plt.savefig(pdfName)
    
  ############################################################################

  def plotNeurons(self,plotIdx=None,pdfName=None):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    if(plotIdx is None):
      plotIdx = range(0,self.neuronCtr)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(self.neuronCoords[plotIdx,0],
               self.neuronCoords[plotIdx,1],
               self.neuronCoords[plotIdx,2],
               'black')

    plt.ion()
    plt.show()
    plt.pause(0.001)

    ax.view_init(150,40); plt.draw()
    plt.axis("off")
    
    #import pdb
    #pdb.set_trace()
    
    if(pdfName is not None):
      plt.savefig(pdfName)

  ############################################################################

  def testPlot(self):

    # !!! NEXT ADD dMIN TO THIS
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


    nPoints = 300

    xTest = np.random.uniform(self.minCoord[0],self.maxCoord[0],nPoints)
    yTest = np.random.uniform(self.minCoord[1],self.maxCoord[1],nPoints)
    zTest = np.random.uniform(self.minCoord[2],self.maxCoord[2],nPoints)
    
    for i in range(0,nPoints):
      
      self.writeLog("Checking " + str(i+1) + "/" + str(nPoints))
      
      if(self.rayCasting(np.array([xTest[i],yTest[i],zTest[i]]))):
        self.writeLog("Inside!")
        color = 'red'
        ax.scatter(xTest[i],yTest[i],zTest[i],color=color)      
      else:
        color = 'black'

      # ax.scatter(xTest[i],yTest[i],zTest[i],color=color)
    
    plt.show()
    plt.pause(0.001)       

  ############################################################################

  def testPlotCached(self):

    # !!! NEXT ADD dMIN TO THIS
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


    nPoints = 1000

    xTest = np.random.uniform(self.minCoord[0],self.maxCoord[0],nPoints)
    yTest = np.random.uniform(self.minCoord[1],self.maxCoord[1],nPoints)
    zTest = np.random.uniform(self.minCoord[2],self.maxCoord[2],nPoints)
    
    for i in range(0,nPoints):
      

      coord = np.array([xTest[i],yTest[i],zTest[i]])
      rayVal = self.rayCasting(coord)
      cachedVal = self.checkInside(coord)

      if(rayVal == cachedVal):
        self.writeLog("Checking " + str(i+1) + "/" + str(nPoints))
        
      elif(rayVal):
        # Inside, but cached was wrong
        color = 'red'
        ax.scatter(xTest[i],yTest[i],zTest[i],color=color)      
      else:
        # Outside, but cached wrong
        color = 'blue'
        ax.scatter(xTest[i],yTest[i],zTest[i],color=color)
    
    plt.show()
    plt.pause(0.001)       

  ############################################################################

  def verifyDmin(self):

    self.writeLog("Verifying that dMin constraint is met")

    minDist = np.zeros((self.neuronCtr,))

    if(self.neuronCtr < 100000):
      neuronRange = range(0,self.neuronCtr)
    else:
      self.writeLog("Too many to check all, picking random neurons to check")
      neuronRange = np.random.randint(0,self.neuronCtr,size=(100000,1))
      self.writeLog(str(neuronRange))
      
    ctr = 0

    for iNeuron in range(0,self.neuronCtr):
      
      ctr = ctr + 1
      if(ctr % 10000 == 0):
        self.writeLog(str(ctr) + "/" + str(len(neuronRange)))
      
      d = np.sqrt(np.sum(np.square(self.neuronCoords[:self.neuronCtr,:]
                                   - self.neuronCoords[iNeuron,:]),
                         axis=1))
      d[iNeuron] = 1e6 # Dont count self distance
      minDist[iNeuron] = np.min(d)

    self.writeLog("Closest neighbour, min = " + str(np.min(minDist)))

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    if(True):
      fig = plt.figure()
      plt.hist(minDist)
      plt.title("Closest neighbour " + str(np.min(minDist)))
      plt.ylabel("Count")
      plt.xlabel("Distance")
      plt.ion()
      plt.show()
      plt.pause(0.001)

    # Also plot where the neurons are relative to the voxel boundaries
    nr = [x for x in neuronRange]
    badIdx = [nr[i] for i in np.where(minDist < self.dMin)[0]]

    fig2 = plt.figure()
    ax = fig2.add_subplot(111, projection='3d')
    
    xc = (self.neuronCoords[badIdx,:] - self.minCoord) % self.binWidth
    ax.scatter(xc[:,0],xc[:,1],xc[:,2],'black')
    plt.title("Bad neuron locations: " + str(xc.shape[0]) )
    plt.axis('tight')
    plt.xlabel('X')
    plt.ylabel('Y')
    #plt.zlabel('Z')
    plt.ion()
    plt.show()
    plt.pause(0.001)
    
    try:
      assert self.dMin <= np.min(minDist), "dMin criteria not fullfilled: " \
        + str(np.min(minDist)) + " < dMin = " + str(self.dMin)
    except:
      self.writeLog("dMin not fullfilled")
      import pdb
      pdb.set_trace()
      
  ############################################################################

  def simpleTestCase(self):

    self.writeLog("This redefines the object, please restart after")

    self.fileName = "cube.obj"
    self.cacheFile = "cube.obj-cached.pickle"
    self.loadMesh(self.fileName)
    self.preCompute()
    self.setupVoxelFilter()

    self.pointOut = np.array([0.5,0.5,-10.5])*1e-6

    for i in range(0,1000):
      testPoint = np.random.uniform(0,1e-6,3) \
                  + np.array([0,0,0])*1e-6
    
    #testPoint = np.array([0.5,0.61,0.5])*1e-6
    
      #self.debugFlag = True
      insideFlag = self.rayCasting(testPoint)

      if(not insideFlag):
        self.writeLog("Test point = " + str(testPoint))
        self.writeLog("wrong!")
        import pdb
        pdb.set_trace()
      

    self.writeLog("All correct")
    self.writeLog("Debug mode")
    import pdb
    pdb.set_trace()
  
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

  def innerVoxelVolume(self):

    return np.sum(self.voxelMaskInner)*(self.binWidth**3)
  
  ############################################################################
  
if __name__ == "__main__":
  
  #sm = RegionMesh("cube.obj",useCache=False)
  #sm.simpleTestCase()
  
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
    
  meshFile = 'mesh/striatum-mesh.obj'
  #meshFile = "mesh/cortex-mesh-200.obj"
  sm = RegionMesh(meshFile,dView=dView,lbView=lbView,
                    raytraceBorders=False)

  # import cProfile
  # cProfile.run("neuronPos = sm.placeNeurons(1000)")

  # sm.plotStruct()

  
  nNeurons = 1730000
  neuronPos = sm.placeNeurons(nNeurons)
  # sm.verifyDmin()
  sm.plotNeurons(pdfName="figures/striatum-fig-somas.png")

  if(False):
    cent = np.array([0.0045,0.0050,0.007])
    # idx = sm.getSubset(cent,radius=1e-3,shape="cube",returnIdxFlag=True)
    # idx = sm.getSubset(cent,radius=1e-3,shape="sphere",returnIdxFlag=True)
    idx = sm.getSubset(cent,nNeurons=1000,shape="cube",returnIdxFlag=True)
    sm.plotNeurons(plotIdx=idx)

  sm.plotStruct(pdfName="figures/striatum-fig-struct.png")
  
  # sm.testPlot()
  # sm.testPlotCached()
  
  #tp = (sm.minCoord + sm.maxCoord)/2
  #sm.rayCasting(np.array(tp))
  
  import pdb
  pdb.set_trace()

  if(dView is not None):
    rc.shutdown(hub=True)

  
