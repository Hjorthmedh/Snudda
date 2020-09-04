#!/usr/bin/env python3

# A wrapper script for the touch detection algorithm
#
# Usage:
#
# snudda init <networkPath> --size XXX
# -- Creates an a json config file
#
# snudda place <networkPath>
# -- Cell placement within volumes specified
#
# snudda detect <networkPath> [--hvsize hyperVoxelSize]
# -- Touch detection of putative synapses
#
# snudda prune <networkPath> [--mergeonly]
# -- Prune the synapses
#
# snudda input <networkPath> [--input yourInputConfig]
#
# snudda export <networkPath>
# -- Export to SONATA format (optional)
#
# snudda simulate <networkPath>
#
# snudda analyse <networkPath>
#
#
# snudda help me

# Johannes Hjorth, Royal Institute of Technology (KTH)
# Human Brain Project 2019

#
# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union’s Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).
#


import os
import sys
import timeit
import numpy as np
import zmq
import pkg_resources

def get_data_file(*dirs):
  path = os.path.join("data", *dirs)
  if not pkg_resources.resource_exists(__package__, path):
    raise FileNotFoundError("Data file '{}' not found".format(path))
  return pkg_resources.resource_filename(__package__, path)

class Snudda(object):

  ############################################################################

  def __init__(self,networkPath):

    if(networkPath[-1] == "/"):
      self.networkPath = networkPath[:-1]
    else:
      self.networkPath = networkPath

    # Add current dir to python path
    sys.path.append(os.getcwd())

    self.start = timeit.default_timer()

  ############################################################################

  def helpInfo(self,args):
    from .snudda_help import snudda_help_text
    print(snudda_help_text())


  ############################################################################

  def initConfig(self,args):
    # self.networkPath = args.path
    print("Creating config file")
    print("Network path: " + str(self.networkPath))

    assert args.size is not None, \
      "You need to speicfy --size when initialising config for network2"

    from .init import SnuddaInit
    structDef = { "Striatum" : args.size,
                  "GPe" : 0,
                  "GPi" : 0,
                  "SNr" : 0,
                  "STN" : 0,
                  "Cortex" : 0,
                  "Thalamus" : 0}
    # Cortex and thalamus axons disabled right now, set to 1 to include one

    if not args.overwrite:
      assert not os.path.exists(self.networkPath), \
        "Network path " + str(self.networkPath) + " already exists" \
        + " (aborting to prevent accidental overwriting)"

    self.makeDirIfNeeded(self.networkPath)

    nPopulationUnits = args.nPopulationUnits
    PopulationUnitCentres = args.PopulationUnitCentres
    PopulationUnitRadius = args.PopulationUnitRadius

    configFile = self.networkPath + "/network-config.json"
    SnuddaInit(structDef=structDef,
               configName=configFile,
               nPopulationUnits=nPopulationUnits,PopulationUnitCentres=PopulationUnitCentres,PopulationUnitRadius=PopulationUnitRadius)

    if(args.size > 1e5):
      print("Make sure there is enough disk space in " + str(self.networkPath))
      print("Large networks take up ALOT of space")


  ############################################################################

  def placeNeurons(self,args):
    # self.networkPath = args.path
    print("Placing neurons")
    print("Network path: " + str(self.networkPath))


    configFile = self.networkPath + "/network-config.json"
    positionFile = self.networkPath + "/network-neuron-positions.hdf5"
    logFileName = self.networkPath + "/log/logFile-place-neurons.txt"

    self.setupLogFile(logFileName) # sets self.logFile
    self.setupParallel() # sets self.dView and self.lbView

    from .place import SnuddaPlace

    if(args.h5legacy):
      h5libver = "earliest"
    else:
      h5libver = "latest" # default

    npn = SnuddaPlace(config_file=configFile,
                      logFile=self.logFile,
                      verbose=True,
                      dView=self.dView,
                      h5libver=h5libver)

    npn.writeDataHDF5(positionFile)

    self.stopParallel()
    self.closeLogFile()


  ############################################################################

  def touchDetection(self,args):
    # self.networkPath = args.path
    print("Touch detection")
    print("Network path: " + str(self.networkPath))

    if(args.hvsize is not None):
      hyperVoxelSize = int(args.hvsize)
    else:
      hyperVoxelSize = 100

    if(args.volumeID is not None):
      volumeID = args.volumeID
    else:
      volumeID = None

    logDir = self.networkPath + "/log"

    configFile = self.networkPath + "/network-config.json"
    positionFile = self.networkPath + "/network-neuron-positions.hdf5"
    logFileName = self.networkPath + "/log/logFile-touch-detection.txt"
    saveFile = self.networkPath + "/voxels/network-putative-synapses.hdf5"


    voxelDir = self.networkPath + "/voxels"
    self.makeDirIfNeeded(voxelDir)

    self.setupLogFile(logFileName) # sets self.logFile
    self.setupParallel() # sets self.dView and self.lbView

    if(args.h5legacy):
      h5libver = "earliest"
    else:
      h5libver = "latest" # default

    from .detect import SnuddaDetect

    if(args.cont):
      # Continue previous run
      print("Continuing previous touch detection")

      ncv = SnuddaDetect(configFile=configFile,
                         positionFile=positionFile,
                         logFile=self.logFile,
                         saveFile=saveFile,
                         SlurmID=self.SlurmID,
                         volumeID=volumeID,
                         rc=self.rc,
                         hyperVoxelSize=hyperVoxelSize,
                         h5libver=h5libver,
                         restartDetectionFlag=False)


    else:
      ncv = SnuddaDetect(configFile=configFile,
                         positionFile=positionFile,
                         logFile=self.logFile,
                         saveFile=saveFile,
                         SlurmID=self.SlurmID,
                         volumeID=volumeID,
                         rc=self.rc,
                         h5libver=h5libver,
                         hyperVoxelSize=hyperVoxelSize)

    self.stopParallel()
    self.closeLogFile()

  ############################################################################

  def pruneSynapses(self,args):
    # self.networkPath = args.path
    print("Prune synapses")
    print("Network path: " + str(self.networkPath))

    from .prune import SnuddaPrune

    logFileName = self.networkPath + "/log/logFile-synapse-pruning.txt"

    workLog = self.networkPath + "/log/network-detect-worklog.hdf5"

    self.setupLogFile(logFileName) # sets self.logFile
    self.setupParallel() # sets self.dView and self.lbView

    # Optionally set this
    scratchPath = None

    if(args.mergeonly):
      preMergeOnly = True
    else:
      preMergeOnly = False

    print("preMergeOnly : " + str(preMergeOnly))

    if(args.h5legacy):
      h5libver = "earliest"
    else:
      h5libver = "latest" # default

    ncvp = SnuddaPrune(workHistoryFile=workLog,
                       logFile=self.logFile,
                       logFileName=logFileName,
                       dView=self.dView, lbView=self.lbView,
                       scratchPath=scratchPath,
                       h5libver=h5libver,
                       preMergeOnly=preMergeOnly)

    self.stopParallel()
    self.closeLogFile()


  ############################################################################

  def setupInput(self,args):

    from .input import SnuddaInput

    print("Setting up inputs, assuming input.json exists")
    logFileName = self.networkPath + "/log/logFile-setup-input.log"
    self.setupLogFile(logFileName) # sets self.logFile
    self.setupParallel() # sets self.dView and self.lbView

    if "input" in args:
      inputConfig = args.input
    else:
      inputConfig = self.networkPath + "/input.json"

    if(not os.path.isfile(inputConfig)):
      print("Missing input config file: " + str(inputConfig))
      return

    if(args.networkFile):
      networkFile = args.networkFile
    else:
      networkFile = self.networkPath \
        + "/network-pruned-synapses.hdf5"

    if(args.inputFile):
      spikeFile = args.inputFile
    else:
      spikeFile = self.networkPath + "/input-spikes.hdf5"

    if(args.time):
      inputTime = args.time

    print("Writing input spikes to " + spikeFile)

    ni = SnuddaInput(inputConfigFile=inputConfig,
                     HDF5networkFile=networkFile,
                     spikeDataFileName=spikeFile,
                     time=inputTime,
                     logFile=self.logFile)

    self.stopParallel()
    self.closeLogFile()

  ############################################################################

  def exportToSONATA(self,args):

    from ConvertNetwork import ConvertNetwork

    print("Exporting to SONATA format")
    print("Network path: " + str(self.networkPath))

    if(args.networkFile):
      networkFile = args.networkFile
    else:
      networkFile = self.networkPath \
        + "/network-pruned-synapses.hdf5"

    if(args.inputFile):
      inputFile = args.inputFile
    else:
      inputFile = self.networkPath + "/input-spikes.hdf5"


    outDir = self.networkPath + "/SONATA/"

    cn = ConvertNetwork(networkFile=networkFile,
                        inputFile=inputFile,
                        outDir=outDir)

  ############################################################################

  def simulate(self,args):

    start = timeit.default_timer()
    
    from .simulate import SnuddaSimulate
    
    if(args.networkFile):
      networkFile = args.networkFile
    else:
      networkFile = self.networkPath \
        + "/network-pruned-synapses.hdf5"

    if(args.inputFile):
      inputFile = args.inputFile
    else:
      inputFile = self.networkPath + "/input-spikes.hdf5"

    self.makeDirIfNeeded(self.networkPath + "/simulation")

    print("Using input file " + inputFile)

    #nWorkers = args.ncores
    #print("Using " + str(nWorkers) + " workers for neuron")

    # Problems with nested symbolic links when the second one is a relative
    # path going beyond the original base path
    if(args.mechDir is None):
      mechDir = os.path.dirname(networkFile) + "/mechanisms"

      # !!! problem with paths, testing to create mechanism dir in current dir
      mechDir = "mechanisms"
      
      if(not os.path.exists(mechDir)):
        mDir = os.path.dirname(__file__) + "/data/cellspecs/mechanisms"
        os.symlink(mDir,mechDir)
    else:
      mechDir = args.mechDir

    # !!! These are saved in current directory x86_64
    # --- problem since nrnivmodl seems to want a relative path...

    makeModsStr = "nrnivmodl " +  mechDir
    if(not os.path.exists('x86_64')):
      print("Please first run: " + makeModsStr)
      exit(-1)
      # I was having problems when running nrnivmodl in the script, but
      # running it manually in bash works... WHY?!!
      
    # os.system(makeModsStr)

    saveDir = os.path.dirname(networkFile) + "/simulation/"

    if(not os.path.exists(saveDir)):
      print("Creating directory " + saveDir)
      os.makedirs(saveDir, exist_ok=True)

    # Get the SlurmID, used in default file names
    SlurmID = os.getenv('SLURM_JOBID')
  
    if(SlurmID is None):
      SlurmID = str(666)


    print("args: " + str(args))
    
    if(args.voltOut is not None):
      # Save neuron voltage
      if(args.voltOut == "default"):
        voltFile = saveDir + 'network-voltage-' + SlurmID + '.csv'
      else:
        voltFile = args.voltOut
    else:
      voltFile = None
    
    if(args.spikesOut is None or args.spikesOut == "default"):
      spikesFile = saveDir + 'network-output-spikes-' + SlurmID + '.txt'
    else:
      spikesFile = args.spikesOut
      
    disableGJ = args.disableGJ
    if(disableGJ):
      print("!!! WE HAVE DISABLED GAP JUNCTIONS !!!")
      
    logFile = os.path.dirname(networkFile) \
      + "/log/network-simulation-log.txt"

    logDir = os.path.dirname(networkFile) + "/log"
    if(not os.path.exists(logDir)):
      print("Creating directory " + logDir)
      os.makedirs(logDir, exist_ok=True)
    
    
    from mpi4py import MPI # This must be imported before neuron, to run parallel
    from neuron import h #, gui
      
    pc = h.ParallelContext()
    
    sim = SnuddaSimulate(networkFile=networkFile,
                         inputFile=inputFile,
                         disableGapJunctions=disableGJ,
                         logFile=logFile,
                         verbose=args.verbose)

    sim.addExternalInput()
    
    sim.checkMemoryStatus()

    if(voltFile is not None):
      sim.addRecording(sideLen=None) # Side len let you record from a subset
      #sim.addRecordingOfType("dSPN",5) # Side len let you record from a subset

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

    
    #cmdStr = "nrnivmodl " + mechDir + " && mpiexec -n " + str(nWorkers) + " -map-by socket:OVERSUBSCRIBE python3 " + os.path.dirname(__file__) + " simulate.py " + networkFile + " " + inputFile + " --time " + str(args.time)

    #if(args.voltOut is not None):
    #  cmdStr += " --voltOut " + args.voltOut

    #os.system(cmdStr)

  ############################################################################

  def analyse(self,args):

    print("Add analysis code here, see Network_analyse.py")

  ############################################################################

  def setupParallel(self):
    self.SlurmID = os.getenv('SLURM_JOBID')

    if(self.SlurmID is None):
      self.SlurmID = self.nextRunID()
    else:
      self.SlurmID = int(self.SlurmID)

    self.logFile.write("Using SlurmID: " + str(self.SlurmID))

    if(os.getenv('IPYTHON_PROFILE') is not None):

      self.logFile.write('Creating ipyparallel client\n')

      from ipyparallel import Client
      #self.rc = Client(profile=os.getenv('IPYTHON_PROFILE'),
      #            # sshserver='127.0.0.1',
      #            debug=False)

      ufile = os.getenv('IPYTHONDIR') + "/profile_" \
              + os.getenv('IPYTHON_PROFILE') \
              + "/security/ipcontroller-client.json"
      self.rc = Client(url_file=ufile, timeout=120, debug=False)

      self.logFile.write('Client IDs: ' + str(self.rc.ids))

      # http://davidmasad.com/blog/simulation-with-ipyparallel/
      # http://people.duke.edu/~ccc14/sta-663-2016/19C_IPyParallel.html
      self.dView = self.rc.direct_view(targets='all') # rc[:] # Direct view into clients
      self.lbView = self.rc.load_balanced_view(targets='all')

      # Define nc globally
      # self.dView.execute("nc = None",block=True)
    else:
      self.logFile.write("No IPYTHON_PROFILE enviroment variable set, running in serial")
      self.dView = None
      self.lbView = None
      self.rc = None

  ############################################################################

  def stopParallel(self):

    # Disable this function, keep the pool running for now
    return

    if(self.rc is not None):
      print("Stopping ipyparallel")
      self.rc.shutdown(hub=True)

  ############################################################################

  def setupLogFile(self, logFileName):
    dataDir = os.path.dirname(logFileName)

    self.makeDirIfNeeded(dataDir)

    try:
      self.logFile = open(logFileName,'w')
      self.logFile.write('Starting log file\n')
    except:
      print("Unable to set up log file " + str(logFileName))

  ############################################################################

  def closeLogFile(self):

    stop = timeit.default_timer()

    print("\nProgram run time: " + str(stop - self.start ))

    self.logFile.write("Program run time: " + str(stop - self.start))
    self.logFile.write("End of log. Closing file.")
    self.logFile.close()

  ##############################################################################


  def nextRunID(self):

    import pickle

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

  def makeDirIfNeeded(self,dirPath):

    if(not os.path.exists(dirPath)):
      print("Creating missing directory " + dirPath)
      try:
        os.makedirs(dirPath)
        print("Created directory " + dirPath)
      except:
        print("Failed to create dir " + dirPath)
        
##############################################################################


if __name__ == "__main__":

  # This is fix to handle if user calles python from within neuron
  import sys
  if('-python' in sys.argv):
    print("Network_simulate.py called through nrniv, fixing arguments")
    pythonidx = sys.argv.index('-python')
    if(len(sys.argv) > pythonidx):
      sys.argv = sys.argv[pythonidx+1:]

  
  from .cli import snudda_cli
  snudda_cli()
