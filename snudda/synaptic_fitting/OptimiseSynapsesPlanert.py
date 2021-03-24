import os
import numpy as np
import time


from OptimiseSynapses import OptimiseSynapses


class OptimiseSynapsesPlanert(OptimiseSynapses):

  def __init__(self,traceList,dView=None,logFileName="logs/planert.txt",
               role="master"):

    self.traceList = traceList
    
    # Create object
    super().__init__(fileName=None,synapseType="gaba",loadCache=False,
                     dView=dView,logFileName=logFileName,role=role)

    self.cacheFileName = traceList + "-parameters.json"
    self.data = {}
    self.nameList = []
    
    # Next load the data
    self.loadVonly(traceList)


    
    # import pdb
    # pdb.set_trace()
    
  ############################################################################
    
  def loadVonly(self,traceList):


    print("Reading traceList: " + str(traceList))
    self.fileName = traceList
    self.traceList = None

    with open(traceList) as f:
      tt = f.readlines()

    dataPath = os.path.dirname(traceList)
      
    for row in tt:
      rowdata = row.split()

      tracefile = rowdata[0]
      dt = float(rowdata[1]) * 1e-6
      timeToFirst = float(rowdata[2]) * 1e-3
      freq = float(rowdata[3])
      nSpikes = int(rowdata[4])
      delayToLast = float(rowdata[5]) * 1e-3

      name = tracefile.replace(".txt","")
      self.nameList.append(name)
      
      tFile = dataPath + "/" + tracefile

      if(not os.path.isfile(tFile)):
        tFileAlt = tFile + "_amp.dat"
        if(os.path.isfile(tFileAlt)):
          # We do not have the file, but we have the extracted amplitudes
          dd = np.genfromtxt(tFileAlt)

          self.data[name] = dict([])
          self.data[name]["tPeak"] = dd[:,0]
          self.data[name]["peakHeight"] = dd[:,1]
          print("Loaded only peaks from " + tFileAlt)
          continue
        else:
          print("File not found: " + tFile)
          continue
      else:
        print("Loading " + tFile)
      
      volt = np.genfromtxt(tFile)
      
      self.data[name] = dict([])
      
      self.data[name]["nPoints"] = nSpikes #  +1 for control spike?
      self.data[name]["time"] = dt*np.arange(0,len(volt))
      self.data[name]["voltage"] = volt

      print("Loaded " + tracefile)
      
      # Find peaks...
      peakTrigger = [timeToFirst]
      for idx in range(1,nSpikes):
        peakTrigger.append(peakTrigger[-1] + 1.0/freq)

      peakTrigger.append(peakTrigger[-1] + delayToLast)

      pTimes = np.array(peakTrigger)

      pWindow = 1.0/(2*freq)*np.ones(pTimes.shape)
      pWindow[-1] *= 5

      self.setPeakInWindow(name,pTimes,pWindow)
      
    return freq      

  ############################################################################
  
  # Look for peaks within a short window (pWindow) after specific time points
  
  def setPeakInWindow(self,name,pTimes,pWindow):

    if(not "time" in self.data[name] or not "voltage" in self.data[name]):
      print("Time or voltage missing from " + str(name))
      return
    
    t = self.data[name]["time"]
    v = self.data[name]["voltage"]
    
    peakIdx = []
    peakTime = []
    peakValue = []
    
    for pt,pw in zip(pTimes,pWindow):
      tStart = pt
      tEnd = pt + pw

      tIdx = np.where(np.logical_and(tStart <= t,t <= tEnd))[0]
      pIdx = tIdx[np.argmax(v[tIdx])]

      peakIdx.append(pIdx)
      peakTime.append(t[pIdx])
      peakValue.append(v[pIdx])

    self.data[name]["peakIdx"] = np.array(peakIdx)
    self.data[name]["tPeak"] = np.array(peakTime)
    self.data[name]["vPeak"] = np.array(peakValue)

  ############################################################################

  def getData(self,dataType,cellID):

    # Ignore dataType, used in data format from Yvonne
    name = self.nameList[cellID]

    return self.data[name]["voltage"],self.data[name]["time"]

  ############################################################################

  def getModelBounds(self,cellID):

    # U, tauR, tauF, tauRatio, cond (obs, tau = tauRatio * tauR)
    modelBounds = ([1e-3,1e-4,1e-4,0, 1e-5],
                   [1.0,2,2,0.9999999,1e-1])

    return modelBounds

  ############################################################################

  # We do not know these properties for the neuron, so assume some values
  # that are not too far off for MS and FS
  
  def getCellProperties(self,dataType,cellID):

    inputResSteadyState = 60e6
    tauDelta = 5e-3

    v,t = self.getData(dataType,cellID)
    baselineDepol = np.mean(v[np.where(t < 0.05)][0])

    return (inputResSteadyState,tauDelta,baselineDepol)
  
  ############################################################################

  def getStimTime(self,dataType,cellID,firstSpike=None,delayToLast=None):

    # Ignoring all parameters but cell ID
    
    stimTime = np.array([0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,1.0])
    
    return stimTime

  ############################################################################

  def getCellType(self,cellID):

    return self.nameList[cellID].split("_")[0]

  ############################################################################

  def setupParallel(self,dView=None):

    assert self.role == "master", "Only master should call setupParallel"

    if(dView is None):
      self.writeLog("No dView, no parallel")
      return
    else:
      self.dView = dView

    if(self.parallelSetupFlag):
      # Already setup servants
      return
    
    with self.dView.sync_imports():
      from run_little_synapse_run import RunLittleSynapseRun
      from OptimiseSynapses import NumpyEncoder
      from OptimiseSynapsesPlanert import OptimiseSynapsesPlanert

    self.writeLog("Setting up workers: " \
          + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    # Create unique log file names for the workers
    if(self.logFileName is not None):
      engineLogFile = [self.logFileName + "-" \
                       + str(x) for x in range(0,len(self.dView))]
    else:
      engineLogFile = [[] for x in range(0,len(self.dView))]

    nWorkers = len(self.dView)
    self.dView.scatter("engineLogFile",engineLogFile) 
    
    self.dView.push({"traceList": self.traceList,
                     "role" : "servant"})

    cmdStr = "ly = OptimiseSynapsesPlanert(traceList=traceList,role=role,logFileName=engineLogFile[0])"
    self.dView.execute(cmdStr,block=True)
    self.parallelSetupFlag = True

  
  ############################################################################
    
if __name__ == "__main__":

  assert False, "This code uses the old incomplete Planert data, for the article we instead resorted to surrogate data. See Planert2010.py and Planert2010part2.py"
  
  print("IPYTHON_PROFILE = " + str(os.getenv('IPYTHON_PROFILE')))
  
  if(os.getenv('IPYTHON_PROFILE') is not None \
     or os.getenv('SLURMID') is not None):
    from ipyparallel import Client
    rc = Client(profile=os.getenv('IPYTHON_PROFILE'),
                debug=False)
    print('Client IDs: ' + str(rc.ids))

    # http://davidmasad.com/blog/simulation-with-ipyparallel/
    # http://people.duke.edu/~ccc14/sta-663-2016/19C_IPyParallel.html
    dView = rc.direct_view(targets='all') # rc[:] # Direct view into clients
    lbView = rc.load_balanced_view(targets='all')
  else:
    dView = None  
  
  traceList = "DATA/Planert2010/d1d2conns/traces/trace_table.txt"
  logFile = "logs/Planert-log.txt"

  osp = OptimiseSynapsesPlanert(traceList,dView=dView,
                                logFileName=logFile)
  
  osp.parallel_optimise_cells("planert",
                              np.arange(0,len(osp.nameList)))

  ops.save_parameter_cache()
  
  # osp.optimiseCell("planert",0)

  #import pdb
  #pdb.set_trace()
