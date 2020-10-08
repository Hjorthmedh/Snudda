import numpy as np
import h5py
import timeit
import json
import sys
import os
from glob import glob

class SnuddaLoad(object):

  ############################################################################

  def __init__(self,network_file,loadSynapses=True):

    if(network_file == "last"):
      network_file = self.findLatestFile()

    # This variable will only be set if the synapses are not kept in
    # memory so we can access them later, otherwise the hdf5 file is
    # automatically closed
    self.hdf5File = None
      
    self.config = None
    self.data = self.loadHDF5(network_file,loadSynapses)
    self.network_file = network_file


  ############################################################################

  def __del__(self):

    if(self.hdf5File is not None):
      try:
        self.hdf5File.close()
      except:
        print("Unable to close HDF5, alread closed?")

  ############################################################################

  def loadHDF5(self, network_file, loadSynapses=True, loadMorph=True):
    print("Loading " + network_file)

    startTime = timeit.default_timer()
    data = dict([])

    f = h5py.File(network_file,'r')

    # with h5py.File(network_file,'r') as f:
    if(True): # Need f open when loadSynapses = False, "with" doesnt work then

      if("config" in f):
        print("Loading config data from HDF5")
        data["config"] = f["config"][()]
        self.config = json.loads(f["config"][()])

      # Added so this code can also load the position file, which
      # does not have the network group yet
      if("network/synapses" in f):
        data["nNeurons"] = f["network/neurons/neuronID"].shape[0]
        data["nSynapses"] = f["network/synapses"].shape[0]
        data["nGapJunctions"] = f["network/gapJunctions"].shape[0]

        if(data["nSynapses"] > 100e6):
          print(str(data["nSynapses"]) + \
                " synapses, which is a lot, not loading them into memory!")
          loadSynapses = False

        # Depricated ??
        if("network/GJIDoffset" in f):
          data["GJIDoffset"] = f["network/GJIDoffset"][()]

        if("network/hyperVoxelIDs" in f):
          data["hyperVoxelIDs"] = f["network/hyperVoxelIDs"][()]

        if(loadSynapses):
          data["synapses"] = f["network/synapses"][:]
          data["gapJunctions"] = f["network/gapJunctions"][:]

          # !!! Convert from voxel idx to coordinates
          data["synapseCoords"] = f["network/synapses"][:,2:5] \
                                * f["meta/voxelSize"][()] \
                                + f["meta/simulationOrigo"][()]
        else:
          # Point the data structure to the synapses and gap junctions on file
          # This will be slower, and only work while the file is open
          data["synapses"] = f["network/synapses"]
          data["gapJunctions"] = f["network/gapJunctions"]

          # We need to keep f alive, since we did not load synapses into
          # the memory
          self.hdf5File = f

          # data["origSynapseCoords"] = f["network/origSynapseCoords"][:]
          # gatheredSynapses = f["network/origGJCoords"][:]
          # data["origGJCoords"] = self.extractSynapseCoords(gatheredSynapses)
      else:
        data["nNeurons"] = f["network/neurons/neuronID"].shape[0]
        assert data["nNeurons"] == f["network/neurons/neuronID"][-1] + 1, \
          "Internal error, something fishy with number of neurons found"



      configFile = f["meta/configFile"][()]
      if(type(configFile) == bytes):
        configFile = configFile.decode()
      data["configFile"] = configFile

      if("meta/positionFile" in f):
        positionFile = f["meta/positionFile"][()]

        if(type(positionFile) == bytes):
          positionFile = positionFile.decode()

        data["positionFile"] = positionFile

      if("meta/SlurmID" in f):
        if(type(f["meta/SlurmID"][()]) == bytes):
          data["SlurmID"] = int(f["meta/SlurmID"][()].decode())
        else:
          data["SlurmID"] = int(f["meta/SlurmID"][()])

      else:
        print("No SlurmID set, using -1")
        data["SlurmID"] = -1

      if("meta/simulationOrigo" in f):
        data["simulationOrigo"] = f["meta/simulationOrigo"][()]

      if("meta/voxelSize" in f):
        data["voxelSize"] = f["meta/voxelSize"][()]

      if("meta/axonStumpIDFlag" in f):
        data["axonStumpIDFlag"] = f["meta/axonStumpIDFlag"][()]

      data["neurons"] = self.extractNeurons(f)

      # This is for old format, update for new format
      if("parameters" in f):
        # print("Parameters found, loading")
        data["synapseRange"] = f["parameters/synapseRange"][()]
        data["gapJunctionRange"] = f["parameters/gapJunctionRange"][()]
        data["minSynapseSpacing"] = f["parameters/minSynapseSpacing"][()]


      data["neuronPositions"] = f["network/neurons/position"][()]

      if("nPopulationUnits" in f["network/neurons"]):
        data["nPopulationUnits"] = f["network/neurons/nPopulationUnits"][()]
        data["populationUnit"] = f["network/neurons/populationUnitID"][()]
        data["populationUnitPlacementMethod"] = f["network/neurons/populationUnitPlacementMethod"][()]
      else:
        print("No Population Units detected.")
        data["nPopulationUnits"] = 0
        data["populationUnit"] = np.zeros(data["nNeurons"])
        data["populationUnitPlacementMethod"] = "none"

      if(loadMorph and "morphologies" in f):
        data["morph"] = dict([])

        for name in f["morphologies"].keys():

          data["morph"][name] = { "swc" :
                                  f["morphologies"][name]["swc"][()],
                                  "location" :
                                  f["morphologies"][name]["location"][()] }


      data["connectivityDistributions"] = dict([])
      #data["connectivityDistributionsGJ"] = dict([])

      if("connectivityDistributions" in f["meta"]):
        origConnectivityDistributions = \
          json.loads(f["meta/connectivityDistributions"][()])

        for keys in origConnectivityDistributions:
          (preType,postType) = keys.split("$$")
          data["connectivityDistributions"][preType,postType] \
            = origConnectivityDistributions[keys]


#      if("connectivityDistributionsGJ" in f["meta"]):
#        origConnectivityDistributionsGJ = \
#          json.loads(f["meta/connectivityDistributionsGJ"][()])
#
#        for keys in origConnectivityDistributionsGJ:
#          (preType,postType) = keys.split("$$")
#          data["connectivityDistributionsGJ"][preType,postType] \
#            = origConnectivityDistributionsGJ[keys]

      if("synapses" in data):
        if("gapJunctions" in data):
          print(str(len(data["neurons"])) + " neurons with " \
                + str(data["synapses"].shape[0]) + " synapses" \
                + " and " + str(data["gapJunctions"].shape[0]) \
                + " gap junctions")
        else:
          print(str(len(data["neurons"])) + " neurons with " \
                + str(data["synapses"].shape[0]) + " synapses")

      print("Load done. " + str(timeit.default_timer() - startTime))

    if(loadSynapses):
      f.close()
    else:
      self.hdf5File = f

    return data

  ############################################################################

  def extractSynapseCoords(self,gatheredSynapses):

    synCoords = dict([])

    for row in gatheredSynapses:
      synCoords[row[8]] = row[0:8]

    return synCoords

  ############################################################################

  def extractNeurons(self,HDF5file):

    if("parameterID" not in HDF5file["network/neurons"]):
      return self.extractNeuronsOLD(HDF5file)
      
    neurons = []


    for name,neuronID,hoc,pos,rot,dendR,axonR,virtual,vID, \
        axonDensityType, axonDensity,axonDensityRadius, \
        axonDensityBoundsXYZ, \
        morph,parameterID,modulationID \
        in zip(HDF5file["network/neurons/name"][:],
               HDF5file["network/neurons/neuronID"][:],
               HDF5file["network/neurons/hoc"][:],
               HDF5file["network/neurons/position"][()],
               HDF5file["network/neurons/rotation"][()],
               HDF5file["network/neurons/maxDendRadius"][:],
               HDF5file["network/neurons/maxAxonRadius"][:],
               HDF5file["network/neurons/virtualNeuron"][:],
               HDF5file["network/neurons/volumeID"][:],
               HDF5file["network/neurons/axonDensityType"][:],
               HDF5file["network/neurons/axonDensity"][:],
               HDF5file["network/neurons/axonDensityRadius"][:],
               HDF5file["network/neurons/axonDensityBoundsXYZ"][:],
               HDF5file["network/neurons/morphology"][:],
               HDF5file["network/neurons/parameterID"][:],
               HDF5file["network/neurons/modulationID"][:]):

      n = dict([])

      if(type(name) == np.ndarray):
        # Old version of savefiles give different output
        name = name[0]
        neuronID = neuronID[0]
        hoc = hoc[0]
        dendR = dendR[0]
        axonR = axonR[0]

      if(type(name) in [bytes, np.bytes_] ):
        n["name"] = name.decode()
      else:
        n["name"] = name

      if(morph is not None):
        if(type(morph) in [bytes,np.bytes_] ):
          n["morphology"] = morph.decode()
        else:
          n["morphology"] = morph

      # Naming convention is TYPE_X, where XX is a number starting from 0
      n["type"] = n["name"].split("_")[0]

      n["neuronID"] = neuronID

      if(type(vID) in [bytes,np.bytes_] ):
        n["volumeID"] = vID.decode()
      else:
        n["volumeID"] = vID

      if(type(hoc) in [bytes, np.bytes_]):
        n["hoc"] = hoc.decode()
      else:
        n["hoc"] = hoc

      n["position"] = pos
      n["rotation"] = rot.reshape(3,3)
      n["maxDendRadius"] = dendR
      n["maxAxonRadius"] = axonR
      n["virtualNeuron"] = virtual

      if(len(axonDensityType) == 0):
        n["axonDensityType"] = None
      elif(type(axonDensityType) in [bytes, np.bytes_]):
        n["axonDensityType"] = axonDensityType.decode()
      else:
        n["axonDensityType"] = axonDensityType

      if(len(axonDensity) > 0):
        if(type(axonDensity) in [bytes,np.bytes_] ):
          n["axonDensity"] = axonDensity.decode()
        else:
          n["axonDensity"] = axonDensity
      else:
        n["axonDensity"] = None


      if(n["axonDensityType"] == "xyz"):
        n["axonDensityBoundsXYZ"] = axonDensityBoundsXYZ
      else:
        n["axonDensityBoundsXYZ"] = None

      n["axonDensityRadius"] = axonDensityRadius

      n["parameterID"] = parameterID
      n["modulationID"] = modulationID
      
      neurons.append(n)

    return neurons


  
  ############################################################################


  # OLD version does not include parameterID and modulationID
  
  def extractNeuronsOLD(self,HDF5file):

    neurons = []

    for name,neuronID,hoc,pos,rot,dendR,axonR,virtual,vID, \
        axonDensityType, axonDensity,axonDensityRadius, \
        axonDensityBoundsXYZ, \
        morph \
        in zip(HDF5file["network/neurons/name"][:],
               HDF5file["network/neurons/neuronID"][:],
               HDF5file["network/neurons/hoc"][:],
               HDF5file["network/neurons/position"][()],
               HDF5file["network/neurons/rotation"][()],
               HDF5file["network/neurons/maxDendRadius"][:],
               HDF5file["network/neurons/maxAxonRadius"][:],
               HDF5file["network/neurons/virtualNeuron"][:],
               HDF5file["network/neurons/volumeID"][:],
               HDF5file["network/neurons/axonDensityType"][:],
               HDF5file["network/neurons/axonDensity"][:],
               HDF5file["network/neurons/axonDensityRadius"][:],
               HDF5file["network/neurons/axonDensityBoundsXYZ"][:],
               HDF5file["network/neurons/morphology"][:]):

      n = dict([])

      if(type(name) == np.ndarray):
        # Old version of savefiles give different output
        name = name[0]
        neuronID = neuronID[0]
        hoc = hoc[0]
        dendR = dendR[0]
        axonR = axonR[0]

      if(type(name) in [bytes, np.bytes_] ):
        n["name"] = name.decode()
      else:
        n["name"] = name

      if(morph is not None):
        if(type(morph) in [bytes,np.bytes_] ):
          n["morphology"] = morph.decode()
        else:
          n["morphology"] = morph

      # Naming convention is TYPE_X, where XX is a number starting from 0
      n["type"] = n["name"].split("_")[0]

      n["neuronID"] = neuronID

      if(type(vID) in [bytes,np.bytes_] ):
        n["volumeID"] = vID.decode()
      else:
        n["volumeID"] = vID

      if(type(hoc) in [bytes, np.bytes_]):
        n["hoc"] = hoc.decode()
      else:
        n["hoc"] = hoc

      n["position"] = pos
      n["rotation"] = rot.reshape(3,3)
      n["maxDendRadius"] = dendR
      n["maxAxonRadius"] = axonR
      n["virtualNeuron"] = virtual

      if(len(axonDensityType) == 0):
        n["axonDensityType"] = None
      elif(type(axonDensityType) in [bytes, np.bytes_]):
        n["axonDensityType"] = axonDensityType.decode()
      else:
        n["axonDensityType"] = axonDensityType

      if(len(axonDensity) > 0):
        if(type(axonDensity) in [bytes,np.bytes_] ):
          n["axonDensity"] = axonDensity.decode()
        else:
          n["axonDensity"] = axonDensity
      else:
        n["axonDensity"] = None


      if(n["axonDensityType"] == "xyz"):
        n["axonDensityBoundsXYZ"] = axonDensityBoundsXYZ
      else:
        n["axonDensityBoundsXYZ"] = None

      n["axonDensityRadius"] = axonDensityRadius
      
      neurons.append(n)

    return neurons

  ############################################################################

  
  def loadConfigFile(self):

    if(self.config is None):
      configFile = self.data["configFile"]
      self.config = json.load(open(configFile,'r'))

  ############################################################################

  def loadNeuron(self,neuronID):

    neuronInfo = self.data["neurons"][neuronID]
    self.loadConfigFile()

    prototypeInfo = self.config[neuronInfo["name"]]

    from Neuron_morphology import NeuronMorphology
    neuron = NeuronMorphology(name=neuronInfo["name"],
                              position=neuronInfo["position"],
                              rotation=neuronInfo["rotation"],
                              swc_filename=prototypeInfo["morphology"],
                              mech_filename=prototypeInfo["mechanisms"],
                              loadMorphology=True)

    return neuron

  ############################################################################

  def synapseIterator(self,chunkSize=1000000,dataType="synapses"):

    typeStrDict = { "synapses" : "network/synapses",
                    "gapJunctions" : "network/gapJunctions" }
    dataStr = typeStrDict[dataType]

    with h5py.File(self.network_file,'r') as f:
      nRows = f[dataStr].shape[0]
      if(nRows == 0):
        # No synapses
        return

      chunkSize = min(nRows,chunkSize)
      nSteps = int(np.ceil(nRows/chunkSize))

      rowStart = 0

      for rowEnd in np.linspace(chunkSize,nRows,nSteps,dtype=int):

        synapses = f[dataStr][rowStart:rowEnd,:]
        rowStart = rowEnd

        yield synapses

  ############################################################################

  def gapJunctionIterator(self,chunkSize=1000000):

    with h5py.File(self.network_file,'r') as f:
      nRows = f["network/gapJunctions"].shape[0]
      chunkSize = min(nRows,chunkSize)

      nSteps = int(np.ceil(nRows/chunkSize))
      rowStart = 0

      for rowEnd in np.linspace(chunkSize,nRows,nSteps,dtype=int):

        gapJunctions = f["network/gapJunctions"][rowStart:rowEnd,:]
        rowStart = rowEnd

        yield gapJunctions


  ############################################################################

  def _rowEvalPostPre(self,row,nNeurons):
    return row[1]*nNeurons + row[0]

  def _rowEvalPost(self,row,nNeurons):
    return row[1]

  ############################################################################

  def findSynapsesSLOW(self,preID,nMax=1000000):

    print("Finding synapses originating from " + str(preID) + ", this is slow")

    synapses = np.zeros((nMax,13),dtype=np.int32)
    synCtr = 0

    if(np.issubdtype(type(preID),np.integer)):
      for synList in self.synapseIterator():
        for syn in synList:
          if(syn[0] == preID):
            synapses[synCtr,:] = syn
            synCtr += 1

    else:
      for synList in self.synapseIterator():
        for syn in synList:
          if(syn[0] in preID):
            synapses[synCtr,:] = syn
            synCtr += 1

    with h5py.File(self.network_file,'r') as f:
      synapseCoords = synapses[:,2:5][:synCtr,:] \
                      * f["meta/voxelSize"][()] \
                      + f["meta/simulationOrigo"][()]


    return (synapses[:synCtr,:],synapseCoords)

  ############################################################################

  # Either give preID and postID, or just postID

  def findSynapses(self,preID=None,postID=None,silent=True):

    if(postID is None):
      return self.findSynapsesSLOW(preID=preID)


    with h5py.File(self.network_file,'r') as f:

      assert postID is not None, "Must specify at least postID"

      nRows = f["network/synapses"].shape[0]
      nNeurons = f["network/neurons/neuronID"].shape[0]

      if(preID is None):
        rowEval = self._rowEvalPost
        valTarget = postID
      else:
        rowEval = self._rowEvalPostPre
        valTarget = postID * nNeurons + preID

      idxA1 = 0
      idxA2 = nRows-1

      idxFound = None

      # We use idxA1 and idxA2 as upper and lower range within which we
      # hope to find one of the synapses. Once we found a synapse row
      # we go up and down in matrix to find the range of the synapses
      # matching the requested condition. This works because the matrix is
      # sorted on postID, and then preID if postID matches

      if(rowEval(f["network/synapses"][idxA1,:],nNeurons) == valTarget):
        idxFound = idxA1

      if(rowEval(f["network/synapses"][idxA2,:],nNeurons) == valTarget):
        idxFound = idxA2

      while(idxA1 < idxA2 and idxFound is None):

        idxNext = int(np.round((idxA1 + idxA2)/2))
        valNext = rowEval(f["network/synapses"][idxNext,:],nNeurons)

        #print("synRow = " + str( f["network/synapses"][idxNext,:]))
        #print("valTarget = " + str(valTarget) + " , valNext = " + str(valNext))
        #print("idxNext = " + str(idxNext), " - " + str(idxA1) + " " + str(idxA2))

        if(valNext < valTarget):
          idxA1 = idxNext
        elif(valNext > valTarget):
          idxA2 = idxNext
        else:
          # We found a hit
          idxFound = idxNext
          break

      if(idxFound is None):
        # No synapses found
        print("No synapses found")
        return None, None

      # Find start of synapse range
      idxB1 = idxFound
      valB1 = rowEval(f["network/synapses"][idxB1-1,:],nNeurons)

      while(valB1 == valTarget and idxB1 > 0):
        idxB1 -= 1
        valB1 = rowEval(f["network/synapses"][idxB1-1,:],nNeurons)

      # Find end of synapse range
      idxB2 = idxFound

      if(idxB2 + 1 < f["network/synapses"].shape[0]):
        valB2 = rowEval(f["network/synapses"][idxB2+1,:],nNeurons)

        while(valB2 == valTarget and idxB2+1 < f["network/synapses"].shape[0]):
          idxB2 += 1
          valB2 = rowEval(f["network/synapses"][idxB2+1,:],nNeurons)

      synapses = f["network/synapses"][idxB1:idxB2+1,:].copy()

      if(not silent):
        print("Synapse range, first " + str(idxB1) + ", last " + str(idxB2))
        print(str(synapses))

      # Calculate coordinates
      synapseCoords = synapses[:,2:5] \
                      * f["meta/voxelSize"][()] \
                      + f["meta/simulationOrigo"][()]

      return synapses,synapseCoords


  ############################################################################

  # !!! DEPRICATED

  def synapseOrigCoordIterator(self,chunkSize=1000000):

    assert False, "This code is now depricated."

    with h5py.File(self.network_file,'r') as f:
      nRows = f["network/synapses"].shape[0]
      chunkSize = min(nRows,chunkSize)

      nSteps = int(np.ceil(nRows/chunkSize))
      rowStart = 0

      for rowEnd in np.linspace(chunkSize,nRows,nSteps,dtype=int):

        synapses = f["network/synapses"][rowStart:rowEnd,:]
        #origCoords = f["network/origSynapseCoords"][rowStart:rowEnd,:]
        rowStart = rowEnd

        # OBS, origCoords is the coordinates in the original swc file
        # coordinate frame (not in the slice coordinates)
        yield synapses #,origCoords

  ############################################################################

  def findLatestFile(self):

    files = glob('save/network-connect-voxel-pruned-synapse-file-*.hdf5')

    modTime = [os.path.getmtime(f) for f in files]
    idx = np.argsort(modTime)

    print("Using the newest file: " + files[idx[-1]])

    return files[idx[-1]]

  ############################################################################

  # Returns cellID of all neurons of neuronType

  def getCellIDofType(self,neuronType,nNeurons=None,randomPermute=False):


    cellID = [x["neuronID"] for x in self.data["neurons"] \
              if x["type"] == neuronType]


    if(nNeurons is not None):
      if(randomPermute):
        # Do not use this if you have a simulation with multiple
        # workers... they might randomize differently, and you might
        # get more or less neurons in total than you wanted
        keepIdx = np.random.permutation(len(cellID))[:nNeurons]
        cellID = np.array([cellID[x] for x in keepIdx])
      else:
        cellID = np.array([cellID[x] for x in range(nNeurons)])

      if(len(cellID) < nNeurons):
        print("getCellIDofType: wanted " + str(nNeurons) \
              + " only got " + str(len(cellID)) \
              + " neurons of type " + str(neuronType))

    # Double check that all of the same type
    assert np.array([self.data["neurons"][x]["type"] == neuronType \
                     for x in cellID]).all()

    return cellID

  ############################################################################


if __name__ == "__main__":

  from argparse import ArgumentParser
  parser = ArgumentParser(description="Load snudda network file (hdf5)")
  parser.add_argument("networkFile", help="Network file (hdf5)",type=str)
  parser.add_argument("--listN", help="Lists neurons in network",
                      action="store_true")
  parser.add_argument("--listT", type=str,
                      help="List neurons of type, --listT ? list the types.",
                      default=None)
  parser.add_argument("--listPre", help="List pre synaptic neurons",
                      type=int)
  parser.add_argument("--listPost", help="List post synaptic neurons (slow)",
                      type=int)
  parser.add_argument("--keepOpen", help="This prevents loading of synapses to memory, and keeps HDF5 file open", action="store_true")

  args = parser.parse_args()

  if(args.keepOpen):
    loadSynapses=False
  else:
    loadSynapses=True

  nl = SnuddaLoad(args.networkFile,loadSynapses=loadSynapses)

  if(args.listN):
    print("Neurons in network: ")

    for nid,name,pos in [(x["neuronID"],x["name"],x["position"])\
                         for x in nl.data["neurons"]]:
      print("%d : %s  (x: %f, y: %f, z: %f)" % (nid,name,pos[0],pos[1],pos[2]))

  if(args.listT is not None):
    if(args.listT == "?"):
      print("List neuron types in network:")

      nTypes = np.unique([x["type"] for x in nl.data["neurons"]])
      for nt in nTypes:
        num = len([x["type"] for x in nl.data["neurons"] if x["type"] == nt])
        print(nt + " (" + str(num) + " total)")

    else:
      print("Neurons of type " + args.listT + ":")
      nOfType = [(x["neuronID"], x["name"]) for x in nl.data["neurons"] \
                 if x["type"] == args.listT]
      for nid,name in nOfType:
        print("%d : %s" % (nid,name))

  if(args.listPre):
    print("List neurons pre-synaptic to neuronID = " + str(args.listPre) \
          + " (" + str(nl.data["neurons"][args.listPre]["name"]) + ")" )
    synapses = nl.findSynapses(postID=args.listPre)
    preID = np.unique(synapses[0][:,0])

    for nid,name in [(x["neuronID"], x["name"]) for x in nl.data["neurons"] \
                     if x["neuronID"] in preID]:
      nSyn = np.sum(synapses[0][:,0] == nid)
      print("%d : %s (%d synapses)" % (nid,name,nSyn))

  if(args.listPost):
    print("List neurons post-synaptic to neuronID = " + str(args.listPost) \
          + " (" + str(nl.data["neurons"][args.listPost]["name"]) + ")" )
    synapses = nl.findSynapses(preID=args.listPost)
    postID = np.unique(synapses[0][:,1])

    for nid,name in [(x["neuronID"], x["name"]) for x in nl.data["neurons"] \
                     if x["neuronID"] in postID]:
      nSyn = np.sum(synapses[0][:,1] == nid)
      print("%d : %s (%d synapses)" % (nid,name,nSyn))

    # List neurons of network


  #syn = nl.findSynapses(22,5)

  #syn2 = nl.findSynapses(postID=5)



  # cellID = nl.getCellIDofType(neuronType="FSN")
