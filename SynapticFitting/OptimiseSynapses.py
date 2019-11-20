import os
import re
import numpy as np
import h5py
import scipy
import scipy.optimize
import matplotlib.pyplot as plt
import neuron
import json
import copy
import time

# !!! Add check that if the voltage is 0 or 5, then the trace is skipped entirely

from RunLittleSynapseRun import RunLittleSynapseRun

############################################################################

# JSON can not handle numpy arrays or numpy ints/floats, fix with a
# custom serialiser

class NumpyEncoder(json.JSONEncoder):
  def default(self, obj):
    print("NumpyEncoder: " + str(type(obj)))
      
    if isinstance(obj, np.integer):
      return int(obj)
    elif isinstance(obj, np.floating):
      return float(obj)
    elif isinstance(obj, np.ndarray):
      return obj.tolist()
    else:
      #return super(NumpyEncoder, self).default(obj)
      return json.JSONEncoder.default(self, obj)

############################################################################
          
class OptimiseSynapses(object):

  ############################################################################
  
  def __init__(self, fileName, synapseType="glut",loadCache=True,
               role="master",dView=None,verbose=True,logFileName=None,
               optMethod="swarm",prettyPlot=True):

    # Parallel execution role, "master" or "servant"
    self.role = role
    self.parallelSetupFlag = False # Set to True when servants are done setup
    self.dView = dView
    self.verbose = verbose
    self.logFileName = logFileName
    self.optMethod = optMethod
    self.nSmoothing = 200 # How many smoothing points do we use?
    self.simTime = 1.5
    
    self.debugParsFlag = False
    self.debugPars = []

    self.prettyPlot = prettyPlot
    
    print("Init optMethod = " + str(optMethod))
    
    if(self.logFileName is not None and len(self.logFileName) > 0):
      print("Log file: " + self.logFileName)
      self.logFile = open(self.logFileName,'w')
    else:
      self.logFile = None
      
    self.figResolution = 300

    self.fileName = fileName
    self.cacheFileName = str(self.fileName) + "-parameters.json"
    self.loadCache = loadCache
    self.synapseType = synapseType

    self.rsrSynapseModel = None
    self.rsrDeltaModel = None

    self.modelInfo = None

    if(loadCache):
      self.loadParameterCache()
    else:
      self.parameterCache = dict([])

    if(fileName is not None):
      self.hFile = h5py.File(fileName,"r")
      # self.hFile.swmr_mode = True # Allow multiple reads from h5py file
    
      self.writeLog("Loading " + str(fileName))

    if(self.role == "master"):
      self.setupParallel(dView=dView)


  ############################################################################
    
  def __delete__(self):

    if(self.hFile is not None):
      self.hFile.close()

    # Save the parameter cache before closing
    self.saveParameterCache()

    if(self.logFile is not None):
      self.logFile.close()
      self.logFile = None
    
  ############################################################################

  def saveParameterCache(self):

    if(self.role != "master"):
      self.writeLog("No servants are allowed to write output to json, ignoring call.")
      return
    
    self.writeLog("Saving parameters to cache file: " + str(self.cacheFileName))

    try:
      with open(self.cacheFileName,"w") as f:
        json.dump(self.parameterCache,f,indent=2,cls=NumpyEncoder)
        f.close()
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      
      self.writeLog("Failed to save cache file ... " + str(self.cacheFileName))
      
  ############################################################################

  def loadParameterCache(self):

    if(os.path.exists(self.cacheFileName)):
      try:
        self.writeLog("Loading cache file " + str(self.cacheFileName))
        with open(self.cacheFileName,"r") as f:
          tmpDict = json.load(f)
          
          self.parameterCache = dict([])
          for k in tmpDict:
            self.parameterCache[int(k)] = tmpDict[k]

          f.close()
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)
        
        self.writeLog("Unable to open " + str(self.cacheFileName))
        self.parameterCache = dict([])
    else:
      # No cache file to load, create empty dictionary
      self.parameterCache = dict([])

  ############################################################################

  def getSubDictionary(self,cellID):

    try:
      if(type(cellID) == list or type(cellID) == np.ndarray):
        retDict = dict([])

        for c in cellID:
          retDict[c] = copy.deepcopy(self.parameterCache[c])

        return retDict
    
      else:
        return self.parameterCache[cellID]
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()
      
  ############################################################################
  
  def addParameterCache(self,cellID,name,value):

    if(cellID not in self.parameterCache):
      self.parameterCache[int(cellID)] = dict([])
    
    self.parameterCache[int(cellID)][name] = value

  ############################################################################

  def getParameterCache(self,cellID,name):

    if(cellID in self.parameterCache and name in self.parameterCache[cellID]):
      return self.parameterCache[cellID][name]
    else:
      return None
    
  ############################################################################

  def getData(self,dataType,cellID=None):

    if(cellID is None):
      data = self.hFile[dataType].value.copy()
    else:
      data = self.hFile[dataType][:,cellID].copy()
      
    # data: 0 = no recording, 5 = recording but no response
      
    if("IGORWaveScaling" in self.hFile[dataType].attrs):
      tStep = self.hFile[dataType].attrs["IGORWaveScaling"][1,0]

      # 0 = no recording
      # 5 = has recording, no respons
      
      assert 0 < tStep and tStep < 1e-2, " Double check tStep = " + str(tStep)
      # Which variable contains the start time? Do not know, so assume it is 0
    
      nPoints = data.shape[0]
      t = 0 + tStep * np.arange(0,nPoints)

      return (data,t)

    else:
      
      return data

  ############################################################################

  def getCellType(self,cellID):

    # If the text is byte encoded, decode it first
    allCellTypes = [ x if type(x) != bytes else x.decode() \
                     for x in self.getData("w_cell").flatten()]

    return allCellTypes[cellID]
  
  ############################################################################
  
  def getCellID(self,dataType,cellType):

    # If the text is byte encoded, decode it first
    allCellTypes = [ x if type(x) != bytes else x.decode() \
                     for x in self.getData("w_cell").flatten()]
    
    
    # Just a check that our assumption on the ID is correct
    # tmp = self.getData("w_num").astype(int)
    # assert (tmp == np.arange(0,len(tmp)+1)).all(), "Numering mis-match!"
     
    cellIDs = [x for x,y in enumerate(allCellTypes) \
               if cellType.lower() in y.lower()]

    validIDs = self.getValidCellID(dataType)

    validCellIDs = np.array(list(set(cellIDs).intersection(validIDs)))
    
    return validCellIDs

  ############################################################################

  def getAllCellID(self):

    return np.arange(0,len(self.getData("w_cell").flatten()),dtype=int)

  ############################################################################

  def getValidCellID(self,dataType):

    allIdx = self.getAllCellID()

    helper = lambda x : self.checkValidData(dataType, x)
    validFlag = [x for x in map(helper,allIdx)]

    self.writeLog("Valid data: " + str(sum(validFlag)))
    
    return allIdx[validFlag]

  ############################################################################

  def getUserID(self,IDstring,dataType=None):

    IDlist = [int(x) for x in IDstring.split(",")]

    if(dataType is not None):
      validID = getValidCellID(dataType)
      validIDlist = [x for x in IDlist if x in validID]
      return validIDlist

    else:
      return IDlist
  
  ############################################################################

  # parDict is the parameters that are associated with cellID


  
  def plotData(self,dataType,cellID=None,params={},show=True,skipTime=0.05,
               prettyPlot=None):

    if(prettyPlot is None):
      prettyPlot = self.prettyPlot
    
    (data,t) = self.getData(dataType,cellID)

    if(data is None):
      self.writeLog(dataType + " " + str(cellID) + ": Nothing to plot")
      return
    
    somaParams = self.getParameterCache(cellID,"soma")

    vPlot = None
    if(somaParams is not None):

      for p in somaParams:
        params[p] = somaParams[p]

      synapseParams = self.getParameterCache(cellID,"synapse")
      if(synapseParams is not None):
        
        for p in synapseParams:
          params[p] = synapseParams[p]
      
        # Get experimental data
        (inputResSteadyState,tauDelta,baselineDepol) \
          = self.getCellProperties(dataType,cellID)

        tStim = self.getStimTime(dataType,cellID)
    
        plotModel = RunLittleSynapseRun(stimTimes=tStim,
                                        holdingVoltage=baselineDepol,
                                        synapseType=self.synapseType,
                                        params=params,
                                        time=self.simTime)
    
        (tPlot,vPlot,iPlot) = plotModel.run2(pars=params)

    tIdx = np.where(skipTime <= t)[0]
        
    plt.figure()

    plt.plot(t[tIdx]*1e3,data[tIdx]*1e3,'r-')
    if(vPlot is not None):
      t2Idx = np.where(skipTime <= tPlot)[0]
      plt.plot(tPlot[t2Idx]*1e3,vPlot[t2Idx]*1e3,'k-')
    # plt.title(dataType + " " + str(cellID))
    cellType = self.getCellType(cellID)
    if(not prettyPlot):
      titleStr = cellType + " " + str(cellID)

      if("nmda_ratio" in params):
        titleStr += "\nU=%.3g, tauR=%.3g, tauF=%.3g, tau=%.3g,\ncond=%.3g, nmda_ratio=%.3g" \
          % (params["U"],
             params["tauR"],
             params["tauF"],
             params["tau"],
             params["cond"],
             params["nmda_ratio"])
      else:
        titleStr += "\nU=%.3g, tauR=%.3g, tauF=%.3g, tau=%.3g,\ncond=%.3g" \
          % (params["U"],
             params["tauR"],
             params["tauF"],
             params["tau"],
             params["cond"])
        

      plt.title(titleStr)

    if(prettyPlot):
      # Draw scalebars
      vScaleX = 1300
      #vMax = np.max(vPlot[np.where(tPlot > 0.050)[0]])
      vBase = vPlot[-1]
      vScaleY1 = vBase*1e3+2
      vScaleY2 = vBase*1e3+1
      tScaleY  = vBase*1e3+1
      tScaleX1 = 1300
      tScaleX2 = 1400
      
      plt.plot([vScaleX,vScaleX],[vScaleY1,vScaleY2],color="black")
      plt.plot([tScaleX1,tScaleX2],[tScaleY,tScaleY],color="black")

      # Mark optogenetical stimulation
      yHeight = float(np.diff(plt.ylim()))/13
      
      tStim = self.getStimTime(dataType,cellID)*1e3
      yStimMarker1 = vPlot[-1]*1e3-1.5*yHeight
      yStimMarker2 = vPlot[-1]*1e3-2.5*yHeight
      for ts in tStim:
        plt.plot([ts,ts],[yStimMarker1,yStimMarker2],color="cyan")      

      plt.axis("off")

      #import pdb
      #pdb.set_trace()
      
      
    plt.xlabel("Time (ms)")
    plt.ylabel("Volt (mV)")
    
    if(not os.path.exists("figures/")):
      os.makedirs("figures/")
    
    baseName = os.path.splitext(os.path.basename(self.fileName))[0]
    figName = "figures/" + baseName + "-"  + dataType \
      + "-" + cellType + "-" + str(cellID) + ".pdf"
    plt.savefig(figName,dpi=self.figResolution)

    if(show):
      plt.ion()
      plt.show()
    else:
      plt.ioff()
      plt.close()
      
    
  ############################################################################

  def getCellProperties(self,dataType,cellID):

    data = self.getData("w_ecode",cellID)
    inputResSteadyState = data[32]*1e6
    tauDelta = data[40]*1e-3
    # baselineDepol = data[58]*1e-3
    baselineDepol = self.getExpBaseline(dataType,cellID)
    
    return (inputResSteadyState,tauDelta,baselineDepol)

  ############################################################################

  def extractInputResTau(self, t,v,curAmp,curStart,curEnd,baseStart,baseEnd):

    # We assume SI units
    tIdxBase = np.where(np.logical_and(baseStart < t, t < baseEnd))[0]
    vBase = np.mean([v[x] for x in tIdxBase])

    tIdxPeak = np.where(np.logical_and(curStart < t, t < curEnd))[0]
    vPeak = np.min(v[tIdxPeak])
    #vPeak = np.max([v[x] for x in tIdxPeak])

    assert np.abs(vPeak-vBase) > np.abs(np.max(v[tIdxPeak]) - vBase), \
      "The code assumes a hyperpolarising pulse, not a peak maximum"
    
    RM = (vPeak-vBase)/curAmp

    idxPostPulse = np.where(curEnd < t)[0]
    idxMaxPostPulse = idxPostPulse[np.argmax(v[idxPostPulse])]
    tMaxPostPulse = t[idxMaxPostPulse]
    
    tIdxDecay = np.where(np.logical_and(curEnd < t,t<tMaxPostPulse))[0]

    decayFunc = lambda x,a,b,c : a*np.exp(-x/b) + c

    tABfit = t[tIdxDecay] - t[tIdxDecay[0]]
    vABfit = v[tIdxDecay]

    p0d = [-0.005,0.01,-0.06]

    if(np.isnan(vABfit).any() or np.isinf(vABfit).any()):
      self.writeLog("We have inifinite or nan values in the voltage")
      import pdb
      pdb.set_trace()
    
    try:
      fitParams,pcov = scipy.optimize.curve_fit(decayFunc,tABfit,vABfit,
                                                p0=p0d)               
      tau = fitParams[1]

    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)

      plt.figure()
      plt.plot(tABfit,vABfit)
      plt.ion()
      plt.show()
      
      import pdb
      pdb.set_trace()        
      
      
    if(False):
      plt.figure()
      plt.plot(tABfit,vABfit,'-r')
      plt.plot(tABfit,decayFunc(tABfit,
                                   fitParams[0],
                                   fitParams[1],
                                   fitParams[2]))
      plt.xlabel("t")
      plt.ylabel("v")
      plt.title("Tau extraction")
      plt.ion()
      plt.show()                  

    # self.writeLog("RM = " + str(RM) + " tau = " + str(tau))

    # Return membrane resistance and tau
    return (RM,tau) 
    

  ############################################################################

  def setupModelDeltaPulse(self,dataType,cellID):

    # Experimental data from file (only baselineDepol needed at this point)
    (inputResSteadyState,tauDelta,baselineDepol) \
      = self.getCellProperties(dataType,cellID)

    # Initial guesses for parameters of soma
    somaDiam = 16e-6     # m
    somaGleak = 0.0003e4 # S/m2

    #(1.0 / tauDelta) / (np.pi*2*16*16) # !!! Check me!!

    params = { "somaDiameter" : somaDiam, "somaGleak" : somaGleak }
    self.rsrDeltaModel = RunLittleSynapseRun(stimTimes=np.array([]),
                                             holdingVoltage=baselineDepol,
                                             synapseType=self.synapseType,
                                             params=params,
                                             time=self.simTime)
    # tau = R * C

    # This is the short hyperpolarising current pulse that we use to
    # estimate tau and input resistance
    # SI units
    self.curInj = -0.1e-9
    self.curStart = 0.5
    self.curDuration = 2e-3
    self.curEnd = self.curStart + self.curDuration
    self.baseStart = 0.2
    self.baseEnd = 0.3

    # Neuron uses natural units !
    self.rsrDeltaModel.IClamp2 = neuron.h.IClamp(self.rsrDeltaModel.soma(0.5))
    self.rsrDeltaModel.IClamp2.amp = self.curInj*1e9
    self.rsrDeltaModel.IClamp2.delay = self.curStart*1e3
    self.rsrDeltaModel.IClamp2.dur = self.curDuration*1e3
    neuron.h.finitialize(self.rsrDeltaModel.holdingVoltage*1e3) # OLD : -70

  ############################################################################

  # x is somaDiameter, Gleak ...
  
  def optHelperFunction(self,x,inputResSteadyState,tauDelta,plotResults=False):

    try:
      somaDiameter,somaGleak = x
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)

      import pdb
      pdb.set_trace()

    print("somaDiameter = " +str(somaDiameter))
    print("somaGleak = " + str(somaGleak))
    
    # Set soma diameter and Gleak (convert to natural units!)
    self.rsrDeltaModel.soma.L = somaDiameter*1e6
    self.rsrDeltaModel.soma.diam = somaDiameter*1e6
    self.rsrDeltaModel.soma.gl_hh = somaGleak *1e-4
    
    # We probably dont need to disable the short hyperpolarising pulse...
    
    # Update the holding current (use same holding voltage as before by default)
    self.rsrDeltaModel.updateHoldingCurrent()

    # Run the simulation    
    neuron.h.tstop = self.rsrDeltaModel.time *1e3# Must set tstop
    neuron.h.run()

    (RM,tau) =self.extractInputResTau(t=np.array(self.rsrDeltaModel.tSave)*1e-3,
                                      v=np.array(self.rsrDeltaModel.vSave)*1e-3,
                                      curAmp=self.curInj,
                                      curStart=self.curStart,
                                      curEnd=self.curEnd,
                                      baseStart=self.baseStart,
                                      baseEnd=self.baseEnd)

    # Return input resistance and tau    
    optError = np.abs((RM-inputResSteadyState)*(tauDelta-tau))

    if(plotResults):
      plt.figure()
      plt.plot(self.rsrDeltaModel.tSave,
               self.rsrDeltaModel.vSave)
      plt.xlabel("Time (ms)")
      plt.ylabel("Volt (mV)")
      plt.title("RM = " + str(RM) + ", tau = " + str(tau))
      plt.ion()
      plt.show()
    
    return optError
    
  ############################################################################
  
  def fitCellProperties(self,dataType,cellID,plotResults=True):

    params = self.getParameterCache(cellID,"soma")
    
    if(params is not None):
      return params
    
    import neuron

    # Experimental data from file
    (inputResSteadyState,tauDelta,baselineDepol) \
      = self.getCellProperties(dataType,cellID)

    # Setup neuron model
    self.setupModelDeltaPulse(dataType,cellID)

    p0 = [16e-6,0.0003e4]
    res = scipy.optimize.minimize(self.optHelperFunction,p0,
                                  args=(inputResSteadyState,tauDelta),
                                  bounds=[(5e-6,50e-6),(0.0003e3,0.0003e5)])

    somaDiameter = res.x[0]
    somaGleak = res.x[1]
    
    self.writeLog("Soma Diameter: " + str(somaDiameter))
    self.writeLog("Soma Gleak: " + str(somaGleak))

    somaDict = { "somaDiameter" : somaDiameter,
                 "somaGleak" : somaGleak }
    self.addParameterCache(cellID,"soma", somaDict)

    self.saveParameterCache()

    if(plotResults and len(res) == 2):
      self.optHelperFunction(res,
                             inputResSteadyState=inputResSteadyState,
                             tauDelta=tauDelta,
                             plotResults=True)
    
    # Remove old model
    self.rsrDeltaModel = None
 
    return somaDict
    
  ############################################################################

  
  def getPeakIdx(self, dataType,cellID,
                 firstSpike=0.1,delayToLast=None):

    pTime = self.getStimTime(dataType=dataType,
                              cellID=cellID,
                              firstSpike=firstSpike,
                              delayToLast=delayToLast)

    freq = 1.0/(pTime[1]-pTime[0])
    pWindow = 1.0/(2*freq)*np.ones(pTime.shape)
    pWindow[-1] *= 5
    
    peakInfo = self.findPeaksHelper(cellID=cellID,
                                    dataType=dataType,
                                    pTime=pTime,
                                    pWindow=pWindow)

    return peakInfo["peakIdx"]

  ############################################################################

  def getPeakIdx2(self,stimTime,time,volt):

    freq = 1.0/(stimTime[1]-stimTime[0])
    pWindow = 1.0/(2*freq)*np.ones(stimTime.shape)
    pWindow[-1] *= 5
    
    peakInfo = self.findPeaksHelper(pTime=stimTime,
                                    pWindow=pWindow,
                                    time=time,
                                    volt=volt)

    return peakInfo["peakIdx"]
   
   
  
  ############################################################################
    
  def getStimTime(self,dataType,cellID,
                  firstSpike=0.1,delayToLast=None):
    
    try:
      freq = float(re.findall(r'H\d+',dataType)[0][1:])
    except:
      self.writeLog("Unable to extract frequency from " + str(dataType))
      import pdb
      pdb.set_trace()

    if(delayToLast is None):
      if(freq == 20.0):
        delayToLast = 0.55
      else:
        delayToLast = 0.5
      
    pTime = 0.1 + np.arange(0,8)*1.0/freq
    pTime = np.append(pTime,pTime[-1] + delayToLast)

    return pTime
    
  ############################################################################

  # Find peaks within pStart[i] and pStart[i]+pWindow[i]
  # The value is not the amplitude of the peak, just the voltage at the peak
  
  def findPeaksHelper(self,pTime,pWindow,
                      cellID=None,dataType=None,
                      time=None,volt=None):

    if(cellID is not None):
      assert dataType is not None, "If you give cellID you need dataType"
      peakData = self.getParameterCache(cellID,"peaks")
      if(peakData is not None):
        return peakData
    
      (volt,time) = self.getData(dataType,cellID) 
    else:
      assert volt is not None and time is not None, \
        "Either use cellID to get time and volt of experimental data, or send time and volt explicitly"
      
    peakIdx = []
    peakTime = []
    peakVolt = []
    
    for pt,pw in zip(pTime,pWindow):
      tStart = pt
      tEnd = pt + pw

      tIdx = np.where(np.logical_and(tStart <= time,time <= tEnd))[0]

      if(self.synapseType == "glut"):
        pIdx = tIdx[np.argmax(volt[tIdx])]
      elif(self.synapseType == "gaba"):
        # We assume that neuron is more depolarised than -65, ie gaba is
        # also depolarising
        pIdx = tIdx[np.argmax(volt[tIdx])]
      else:
        self.writeLog("Unknown synapse type : " +str(self.synapseType))
        import pdb
        pdb.set_trace()
        
      peakIdx.append(int(pIdx))
      peakTime.append(time[pIdx])
      peakVolt.append(volt[pIdx])
      
    # Save to cache -- obs peakVolt is NOT amplitude of peak, just volt

    peakDict = { "peakIdx" : np.array(peakIdx),
                 "peakTime" : np.array(peakTime),
                 "peakVolt" : np.array(peakVolt)} # NOT AMPLITUDE

    if(cellID is not None):
      self.addParameterCache(cellID,"peaks",peakDict)
                           
    return peakDict
  # (peakIdx,peakTime,peakVolt)

  ############################################################################

  def findTraceHeights(self,time,volt,peakIdx):

    decayFunc = lambda x,a,b,c : a*np.exp(-x/b) + c
    
    vBase = np.mean(volt[int(0.3*peakIdx[0]):int(0.8*peakIdx[0])])

    peakHeight = np.zeros((len(peakIdx,)))
    peakHeight[0] = volt[peakIdx[0]] - vBase

    decayFits = []

    
    for idxB in range(1,len(peakIdx)):

      if(peakHeight[0] > 0):
        if(idxB < len(peakIdx) -1):        
          p0d = [0.06,0.05,-0.074]
        else:
          p0d = [1e-5,100,-0.074]

          if(self.synapseType == "gaba"):
            p0d = [1e-8,10000,-0.0798]  
      else:
        # In some cases for GABA we had really fast decay back
        if(idxB < len(peakIdx) -1):        
          p0d = [-0.06,0.05,-0.0798]
        else:
          p0d = [-1e-5,1e5,-0.0798]

      idxA = idxB -1

      peakIdxA = peakIdx[idxB-1] # Prior peak
      peakIdxB = peakIdx[idxB] # Next peak

      if(idxB < len(peakIdx) -1):
        # Not the last spike
      
        idxStart = int(peakIdxA*0.9 + peakIdxB*0.1)
        idxEnd = int(peakIdxA*0.1 + peakIdxB*0.9)
      else:
        # Last spike, use only last half of decay trace
        idxStart = int(peakIdxA*0.5 + peakIdxB*0.5)
        idxEnd = int(peakIdxA*0.05 + peakIdxB*0.95)

      try:
        assert idxStart < idxEnd
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)

        plt.figure()
        plt.plot(time,volt)
        plt.xlabel("Time (error plot)")
        plt.ylabel("Volt (error plot)")
        plt.ion()
        plt.show()
        plt.title("ERROR!!!")
        import pdb
        pdb.set_trace()
        
      tAB = time[idxStart:idxEnd]
      vAB = volt[idxStart:idxEnd]        

      tABfit = tAB - tAB[0]
      vABfit = vAB
     
      try:

        try:
          fitParams,pcov = scipy.optimize.curve_fit(decayFunc,tABfit,vABfit,
                                                    p0=p0d)
        except:
          import traceback
          tstr = traceback.format_exc()
          self.writeLog(tstr)
          
          self.writeLog("!!! Failed to converge, trying with smaller decay constant")
          p0d[1] *= 0.01
          fitParams,pcov = scipy.optimize.curve_fit(decayFunc,tABfit,vABfit,
                                                    p0=p0d)
          
        tB = time[peakIdxB] - tAB[0]
        vBaseB = decayFunc(tB,fitParams[0],fitParams[1],fitParams[2])
        
        peakHeight[idxB] = volt[peakIdxB] - vBaseB

        vFit = decayFunc(tAB-tAB[0],fitParams[0],fitParams[1],fitParams[2])
        decayFits.append((tAB,vFit))
        
      except:
        
        self.writeLog("Check that the threshold in the peak detection before is OK")
        #self.plot(name)
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)

        if(True):
          plt.figure()
          plt.plot(tAB,vAB,'r')
          plt.title("Error in findTraceHeights")
          plt.xlabel("time")
          plt.ylabel("volt")
          #plt.plot(tAB,vFit,'k-')
          plt.ion()
          plt.show()

        import pdb
        pdb.set_trace()

    #import pdb
    #pdb.set_trace()
    
        
    return peakHeight.copy(), decayFits,vBase
        
  ############################################################################

  def setupModelSynapseFitting(self,dataType,cellID,params={}):

    # If the delta model existed, clear it
    #self.rsrDeltaModel = None
    
    tStim = self.getStimTime(dataType,cellID)
    
    somaParams = self.getParameterCache(cellID,"soma")
    assert somaParams is not None, "You must fit soma parameters first"

    # !!! Should we also try and load synapse parameters from cache here??
    for p in somaParams:
      params[p] = somaParams[p]

    # Get experimental data
    (inputResSteadyState,tauDelta,baselineDepol) \
      = self.getCellProperties(dataType,cellID)

    # !!! We need to get the baseline depolarisation in another way

    self.rsrSynapseModel = RunLittleSynapseRun(stimTimes=tStim,
                                               holdingVoltage=baselineDepol,
                                               synapseType=self.synapseType,
                                               params=params,
                                               time=self.simTime)

  ############################################################################

  def neuronSynapseSwarmHelper(self,pars,tSpikes,peakHeight,
                               smoothExpTrace):

    if(self.debugParsFlag):
      self.debugPars.append(pars)
    
    res = np.zeros((pars.shape[0]))
    
    for idx,p in enumerate(pars):
      peakH,tSim,vSim = self._neuronSynapseSwarmHelper(p,tSpikes)

      # Calculating error in peak height
      hDiff = np.abs(peakH - peakHeight)
      hDiff[0] *= 3
      hDiff[-1] *= 3
      hError = np.sum(hDiff)/len(hDiff)

      # Calculate error in decay fit
      simTrace,simTime = self.smoothingTrace(vSim,self.nSmoothing,
                                             time=tSim,
                                             startTime=self.decayStartFit,
                                             endTime=self.decayEndFit)

      # We only want to use the bit of the trace after max
      idxMax = np.argmax(smoothExpTrace)

      # We divide by number of points in vector, to get the average deviation
      # then we multiply by 10000 to get an error comparable to the others
      decayError = np.sum((smoothExpTrace[idxMax:] \
                           - simTrace[idxMax:])**2) \
                           /(self.nSmoothing-idxMax+1) * 2000

      if(False):
        plt.figure()
        plt.plot(smoothExpTrace[idxMax:])
        plt.plot(simTrace[idxMax:])
        plt.ion()
        plt.show()
        import pdb
        pdb.set_trace()
      
      res[idx] = hError + decayError
      
    return res

  ############################################################################

  def smoothingTrace(self,originalTrace,nParts,time=None,
                     startTime=None,endTime=None):

    if(time is not None):
      tFlag = np.ones((len(originalTrace),),dtype=bool)
      
      if(endTime is not None):
        tFlag[np.where(time > endTime)[0]] = False

      if(startTime is not None):
         tFlag[np.where(time < startTime)[0]] = False
         
      #tIdx = np.where(tFlag)[0]
      trace = originalTrace[tFlag]
      t = time[tFlag]
    else:
      trace = originalTrace
      t = time
    
    N = int(np.round(len(trace)/nParts))
    
    smoothTrace = np.convolve(trace, np.ones((N,))/N, mode='valid')

    idx = np.linspace(0,len(smoothTrace)-1,num=nParts,dtype=int)
    
    return smoothTrace[idx], t[idx]
  
  ############################################################################
  
  def _neuronSynapseSwarmHelper(self,
                               pars,
                               tSpikes):

    if(self.synapseType == "glut"):
      U,tauR,tauF,tauRatio,cond,nmdaRatio = pars
      tau = tauR*tauRatio
      params = { "nmda_ratio" : nmdaRatio }
    elif(self.synapseType == "gaba"):
      U,tauR,tauF,tauRatio,cond = pars
      tau = tauR*tauRatio
      params = {}
    else:
      self.writeLog("Unknown synapse type: " + str(self.synapseType))
      import pdb
      pdb.set_trace()
    
    peakHeights,tSim,vSim = self.synapseModelNeuron(tSpikes,U,
                                                    tauR,tauF,cond,tau,
                                                    params=params,
                                                    returnTrace=True)
    
    return peakHeights,tSim,vSim
  
  ############################################################################
  
  def neuronSynapseHelper(self,
                          tSpike,U,tauR,tauF,
                          tauRatio=None,
                          cond = 1e-7, tau=None):

    assert tauRatio is None or tau is None, \
      "Only one of tau and tauRatio should be defined"
    
    if(tauRatio is not None):
      tau = tauR*tauRatio
    elif(tau is None):
      assert False, "tau or tauRatio must be specified"

    peakHeights = self.synapseModelNeuron(tSpike,U,tauR,tauF,cond,tau)
    
    return peakHeights
    
   ############################################################################
  
  def synapseModelNeuron(self,tSpike,U,tauR,tauF,cond,tau,
                         params = {},
                         returnTrace=False):

    # self.writeLog("Running neuron model")
    
    assert self.rsrSynapseModel is not None, \
      "!!! Need to call setupModelSynapseFitting first"

    # Should we make a copy of params, to not destroy it? ;)
    params["U"]    = U
    params["tauR"] = tauR
    params["tauF"] = tauF
    params["cond"] = cond
    params["tau"]  = tau

    #self.writeLog("params=" + str(params))
    
    (tSim,vSim,iSim) = \
      self.rsrSynapseModel.run2(pars=params)

    if(tSim.shape != vSim.shape):
      self.writeLog("Shape are different, why?!")
      import pdb
      pdb.set_trace()            
    
    peakIdx = self.getPeakIdx2(time=tSim,volt=vSim,stimTime=tSpike)
    peakHeight,decayFits,vBase = self.findTraceHeights(tSim,vSim,peakIdx)

    if(returnTrace):
      return peakHeight,tSim,vSim
    else:
      return peakHeight
   
  ############################################################################

  def getExpBaseline(self,dataType,cellID,tBefore=0.08):

     (volt,time) = self.getData(dataType,cellID)

     idx = np.where(time < tBefore)[0]

     return np.mean(volt[idx])
  
  
  ############################################################################

  def getModelBounds(self,cellID):

    cellType = self.getCellType(cellID)

    if("FS" in cellType.upper()):
      # U, tauR, tauF, tauRatio, cond (obs, tau = tauRatio * tauR), nmda_ratio
      modelBounds = ([1e-3,1e-4,1e-4,0, 1e-5,0.1],
                     [1.0,2,2,0.9999999,1e-1,1])

    elif("MS" in cellType.upper()):
      # U, tauR, tauF, tauRatio, cond (obs, tau = tauRatio * tauR), nmda_ratio
      modelBounds = ([1e-3,1e-4,1e-4,0, 1e-5,0.3],
                     [1.0,2,2,0.9999999,1e-1,6])

    elif("CHAT" in cellType.upper()):
      # U, tauR, tauF, tauRatio, cond (obs, tau = tauRatio * tauR), nmda_ratio
      modelBounds = ([1e-3,1e-4,1e-4,0, 1e-5,2],
                     [1.0,2,2,0.9999999,1e-1,8])

    elif("LTS" in cellType.upper()):
      # U, tauR, tauF, tauRatio, cond (obs, tau = tauRatio * tauR), nmda_ratio
      modelBounds = ([1e-3,1e-4,1e-4,0, 1e-5,2],
                     [1.0,2,2,0.9999999,1e-1,8])
    else:
      self.writeLog("Unknown celltype in " + str(cellType))
      import pdb
      pdb.set_trace()

    return modelBounds
    
  
  ############################################################################
  
  def fitTrace(self,dataType,cellID,optMethod=None):
    
    if(optMethod is None):
      optMethod = self.optMethod
    
    (volt,time) = self.getData(dataType,cellID)

    # This should NOT be empty, it should have AMPA/NMDA specific parameters
    # that we do not optimise for...
    params = dict([])
    
    synapseModel = self.setupModelSynapseFitting(dataType=dataType,
                                                 cellID=cellID,
                                                 params=params)




    peakIdx = self.getPeakIdx(dataType,cellID)
    tSpikes = time[peakIdx]
    
    sigma = np.ones(len(peakIdx))
    sigma[-1] = 1./3

    stimTime =  self.getStimTime(dataType=dataType,
                                 cellID=cellID)
    
    peakHeight,decayFits,vBase = self.findTraceHeights(time,volt,peakIdx)

    # U, tauR, tauF, tauRatio, cond (obs, tau = tauRatio * tauR)
    modelBounds = ([1e-3,1e-4,1e-4,0, 1e-5],[1.0,2,2,0.9999999,1e-1])
    # !!! bestRandom should also see these bounds...

    try:
      
      if(optMethod == "stupid"):

        assert False, "Use the swarm method!"
        # U, tauR, tauF, tauRatio, cond (obs, tau = tauRatio * tauR)
        modelBounds = ([1e-3,1e-4,1e-4,0, 1e-5],[1.0,2,2,0.9999999,1e-1])

        # This one gives nmda_ratio as parameter also
        #modelBounds = self.getModelBounds(cellID)
        
        # Call best random, to find a good set of starting points
        p0 = self.bestRandom(synapseModel=synapseModel,
                             cellID=cellID,
                             tPeak = stimTime,
                             hPeak = peakHeight,
                             modelBounds=modelBounds)


        fitParams,pcov = scipy.optimize.curve_fit(self.neuronSynapseHelper,
                                                  stimTime,peakHeight,
                                                  sigma=sigma,
                                                  absolute_sigma=False,
                                                  p0=p0,
                                                  bounds=modelBounds)

        # tau < tauR, so we use tauRatio for optimisation
        fitParams[3] *= fitParams[1] # tau = tauR * tauRatio

        modelHeights = self.neuronSynapseHelper(stimTime,
                                                U=fitParams[0],
                                                tauR=fitParams[1],
                                                tauF=fitParams[2],
                                                tau=fitParams[3],
                                                cond=fitParams[4])


        self.writeLog("Parameters: U = %.3g, tauR = %.3g, tauF = %.3g, tau = %.3g, cond = %3.g" % tuple(fitParams))
        
      elif(optMethod=="swarm"):

        if(self.debugParsFlag):
          self.debugPars = []

        # Increase upper cond limit to 1? -- M1LH 11 might need it
          
        # U, tauR, tauF, tauRatio, cond (obs, tau = tauRatio * tauR), nmda_ratio
        #modelBounds = ([1e-3,1e-4,1e-4,0, 1e-5,0.3],
        #               [1.0,2,2,0.9999999,1e-1,4])

        modelBounds = self.getModelBounds(cellID)

        self.decayStartFit = 1.0
        self.decayEndFit = 1.2
         
        smoothVolt,smoothTime = self.smoothingTrace(volt,self.nSmoothing,
                                                 time=time,
                                                 startTime = self.decayStartFit,
                                                 endTime = self.decayEndFit)

        import pyswarms as ps
        from pyswarms.utils.functions import single_obj as fx

        options = {"c1":0.5, "c2":0.3, "w":0.9} # default
        
        # Increased from 200 to 300, to 400
        if(self.synapseType == "glut"):
          nPart = 400
        else:
          nPart = 50
          
        optimizer = ps.single.GlobalBestPSO(n_particles=nPart,
                                            dimensions=len(modelBounds[0]),
                                            options=options,
                                            bounds=modelBounds)

        # Set iter to 50, lowered to 10
        cost, fitParams = optimizer.optimize(self.neuronSynapseSwarmHelper,
                                             iters=10,
                                             tSpikes = stimTime,
                                             peakHeight = peakHeight,
                                             smoothExpTrace=smoothVolt)

        modelHeights,tSim,vSim = \
          self._neuronSynapseSwarmHelper(fitParams,stimTime)
        
        # tau < tauR, so we use tauRatio for optimisation
        fitParams[3] *= fitParams[1] # tau = tauR * tauRatio

        if(len(fitParams) == 6):
          self.writeLog("Parameters: U = %.3g, tauR = %.3g, tauF = %.3g, tau = %.3g, cond = %3.g, nmda_ratio = %.3g" % tuple(fitParams))
        elif(len(fitParams) == 5):
          self.writeLog("Parameters: U = %.3g, tauR = %.3g, tauF = %.3g, tau = %.3g, cond = %3.g" % tuple(fitParams))
        else:
          self.writeLog("Unknown parameter format: " + str(fitParams))
          import pdb
          pdb.set_trace()
      else:
        self.writeLog("NO METHOD SELECTED!!!")
        exit(-1)
        
      self.writeLog("peakHeight = " + str(peakHeight))
      self.writeLog("modelHeights = " + str(modelHeights))
            
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()


    if(len(fitParams) == 6):
      self.addParameterCache(cellID,"synapse", \
                             { "U" : fitParams[0],
                               "tauR" : fitParams[1],
                               "tauF" : fitParams[2],
                               "tau"  : fitParams[3],
                               "cond" : fitParams[4],
                               "nmda_ratio" : fitParams[5] })
    elif(len(fitParams) == 5):
      self.addParameterCache(cellID,"synapse", \
                             { "U" : fitParams[0],
                               "tauR" : fitParams[1],
                               "tauF" : fitParams[2],
                               "tau"  : fitParams[3],
                               "cond" : fitParams[4] })
    else:
      print("Unknown fitParams format, unable to save paramters = " \
            +str(fitParams))
      import pdb
      pdb.set_trace()
      

    self.addParameterCache(cellID,"meta", \
                           { "type" : self.getCellType(cellID),
                             "experiment" : dataType,
                             "datafile" : self.fileName })

    if(self.debugParsFlag):
      self.plotDebugPars()
    
    # Deallocate Neuron model
    self.rsrSynapseModel = None
    
  ############################################################################

  def bestRandom(self,synapseModel,cellID,
                 tPeak,hPeak,
                 modelBounds,
                 nTrials=5000,loadParamsFlag=False):

    assert nTrials >= 1, "nTrials should be a positive integer"

    minError = np.inf
    minPar = None

    for idx in range(0,nTrials):
      if(idx % 100 == 0):
        self.writeLog("Pre-trial : " + str(idx) + "/" + str(nTrials))

      if(idx == 0 and loadParamsFlag):
        # If we should load params then do so first
        pars = self.getParameterCache(cellID,"synapse")
      else:
        pars = None
        
      if(pars is not None):
        U    = pars["U"]
        tauR = pars["tauR"]
        tauF = pars["tauF"]
        tau  = pars ["tau"]
        cond = pars["cond"]
      else:
        U = np.random.uniform(modelBounds[0][0],modelBounds[1][0])
        tauF = np.random.uniform(modelBounds[0][1],modelBounds[1][1])
        tauR = np.random.uniform(modelBounds[0][2],modelBounds[1][2])
        tau = tauR * np.random.uniform(modelBounds[0][3],modelBounds[1][3])
        cond = np.random.uniform(modelBounds[0][4],modelBounds[1][4])

      try:
        peakHeights = self.synapseModelNeuron(tPeak,U,tauR,tauF,tau,cond)
      
        error = np.abs(peakHeights - hPeak)
        error[0] *= 3
        error[1] *= 2
        error[-1] *= 3
        error = np.sum(error)
      
        if(error < minError):
          minError = error
          minPar = np.array([U,tauR,tauF,tau/tauR,cond])

      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)
        
        # For big runs we do no want to give up. Let's try again...
        continue

        
    return minPar

  ############################################################################

  def optimiseCell(self,dataType,cellID):

    try:
      cellID = int(cellID)

      self.writeLog("Optimising " + str(dataType) + " ID = " + str(cellID))
    
      self.fitCellProperties(dataType,cellID)
      self.fitTrace(dataType,cellID)
      self.plotData(dataType,cellID,show=False)

      # We need to extract the dictionary associated with cellID and
      # return the result so that the main function can plug it into the
      # central parameter cache -- OR DO WE?

      return self.getSubDictionary(cellID)
    
    except:
      
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      self.addParameterCache(cellID,"error",tstr)
      

  ############################################################################

  def parallelOptimiseCells(self,dataType,cellIDlist=None):

    if(self.role == "master"):

      if(cellIDlist is None):
        cellIDlist = self.getValidCellID(dataType)

      self.printCellProperties(dataType,cellIDlist)
        
      if(self.dView is None):
        self.writeLog("We are in serial mode, use serial code..")
        for c in cellIDlist:
          self.optimiseCell(dataType,c)

        # Save parameters to file
        self.saveParameterCache()
        return

      self.writeLog("Optimising in parallel")
      
      self.setupParallel(self.dView)
      self.dView.scatter("cellIDlist",cellIDlist,block=True)
      self.dView.push({"dataType" : dataType},block=True)
      
      try:

        cmdStr = "res = ly.parallelOptimiseCells(dataType=dataType,cellIDlist=cellIDlist)"
        self.writeLog("Starting execution")
        self.dView.execute(cmdStr,block=True)
        import pdb
        pdb.set_trace()
        self.writeLog("Execution finished, gathering results")
        #res = self.dView.gather("res",block=True)
        res = self.dView["res"]
      
        # res now contains a list of dictionaries, each dictionary from one worker, we need
        # to merge these
        self.writeLog("Results gathered, merging dictionary")
        for m in res:
          #self.parameterCache.update(m)
          for k in m:
            self.parameterCache[int(k)] = m[k]
          
        self.saveParameterCache()

        self.writeLog("Done.")
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)
        import pdb
        pdb.set_trace()
        
    else:
      try:
        # Servant runs here
        for cellID in cellIDlist:
          self.optimiseCell(dataType,cellID)

        return self.getSubDictionary(cellIDlist)
      
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)
        import pdb
        pdb.set_trace()
     
  ############################################################################

  def printCellProperties(self,dataType,cellID):

    inputResList = []
    tauDeltaList = []
    baselineDepolList = []
    
    for cid in cellID:

      inputResSteadyState,tauDelta,baselineDepol \
        = self.getCellProperties(dataType,cid)

      inputResList.append(inputResSteadyState)
      tauDeltaList.append(tauDelta)
      baselineDepolList.append(baselineDepol)

    inputRes = np.array(inputResList)
    tauDelta = np.array(tauDeltaList)
    baselineDepol = np.array(baselineDepolList)
      
    plt.figure()
    plt.scatter(inputRes*1e-6,tauDelta*1e3)
    plt.xlabel("Input resistance (Mohm)")
    plt.ylabel("Tau (ms)")
    plt.ion()
    plt.show()

    plt.figure()
    plt.scatter(inputRes*1e-6,baselineDepol*1e3)
    plt.xlabel("Input resistance (Mohm)")
    plt.ylabel("Baseline depolarisation (mV)")
    plt.ion()
    plt.show()

    
  
  ############################################################################
  
  def checkValidData(self,dataType,cellID):

    data = self.hFile[dataType][:,cellID]
    goodDataFlag = True
    
    if((data == 0).all()):
      self.writeLog("No recording for " + str(dataType) + " ID = " + str(cellID))
      goodDataFlag = False
      
    if((data == 5).all()):
      self.writeLog("No response for " + str(dataType) + " ID = " + str(cellID))
      goodDataFlag = False

    (inputResSteadyState,tauDelta,baselineDepol) = \
      self.getCellProperties(dataType,cellID)

    cellType = self.getCellType(cellID)
    
    if("LTS" in cellType.upper()):
      inputResMaxLimit = 1000e6
    else:
      inputResMaxLimit = 500e6
    
    if(np.isnan(inputResSteadyState) or np.isnan(tauDelta) \
       or inputResSteadyState > inputResMaxLimit
       or inputResSteadyState < 5e6
       or tauDelta > 5):
      self.writeLog("Bad values for input resistance and/or tauDelta for " \
                    + str(dataType) + " ID = " +str(cellID) \
                    + "(" + str(self.getCellType(cellID)) + ")")
      goodDataFlag = False
      
    return goodDataFlag
      
  ############################################################################
  
  def listDataTypes(self):

    return [x for x in self.hFile]

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
      from RunLittleSynapseRun import RunLittleSynapseRun
      from OptimiseSynapses import NumpyEncoder
      from OptimiseSynapses import OptimiseSynapses

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
    
    self.dView.push({"fileName": self.fileName,
                     "synapseType" : self.synapseType,
                     "loadCache" : self.loadCache,
                     "role" : "servant"})

    cmdStr = "ly = OptimiseSynapses(fileName=fileName, synapseType=synapseType,loadCache=loadCache,role=role,logFileName=engineLogFile[0])"
    self.dView.execute(cmdStr,block=True)
    self.parallelSetupFlag = True

  ############################################################################

  def writeLog(self,text,flush=True): # Change flush to False in future, debug
    if(self.logFile is not None):
      self.logFile.write(text + "\n")
      
      if(self.verbose):
        print(text)
        
      if(flush):
        self.logFile.flush()
    else:
      if(self.verbose):
        print(text)  

  ############################################################################

  def plotDebugPars(self):

    try:
      nIter = len(self.debugPars)
      nPoints = self.debugPars[0].shape[0]
      nPars = self.debugPars[0].shape[1]

      assert nPars == 6, "Should be six parameters"

      Uall      = np.zeros((nPoints,nIter))
      tauR      = np.zeros((nPoints,nIter))
      tauF      = np.zeros((nPoints,nIter))
      tauRatio  = np.zeros((nPoints,nIter))
      cond      = np.zeros((nPoints,nIter))
      nmdaRatio = np.zeros((nPoints,nIter))

      for ctr,par in enumerate(self.debugPars):
        Uall[:,ctr]      = par[:,0]
        tauR[:,ctr]      = par[:,1]
        tauF[:,ctr]      = par[:,2]
        tauRatio[:,ctr]  = par[:,3]
        cond[:,ctr]      = par[:,4]
        nmdaRatio[:,ctr] = par[:,5]
      
      
      plt.figure()
      plt.plot(Uall,cond,'-')
      plt.xlabel("U")
      plt.ylabel("cond")

      plt.figure()
      plt.plot(tauR,tauF,'-')
      plt.xlabel("tauR")
      plt.ylabel("tauF")

      plt.figure()
      plt.plot(tauRatio,nmdaRatio)
      plt.xlabel("tauRatio")
      plt.ylabel("nmdaRatio")

      plt.figure()
      plt.plot(Uall)
      plt.xlabel("Uall")

      plt.figure()
      plt.plot(nmdaRatio)
      plt.xlabel("nmdaRatio")
      
      plt.ion()
      plt.show()
      
      import pdb
      pdb.set_trace()

    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)

      import pdb
      pdb.set_trace()
      
  ############################################################################
      
if __name__ == "__main__":

  import argparse
  parser = argparse.ArgumentParser(description="Extract synaptic parameters from electrophysiological recordings")
  parser.add_argument("file", help="HDF5 file to read in")
  parser.add_argument("--st", help="Synapse type (glut or gaba)",
                      choices=["glut","gaba"])
  parser.add_argument("--type", help="Cell types",
                      choices=["MSN","FS","LTS","CHAT"])
  parser.add_argument("--id",
                      help="Cell ID, comma separated no spaces, eg. 1,2,3,7")
  #parser.add_argument("--optMethod",
  #                    help="Optimisation method",
  #                    choices=["stupid","swarm"],
  #                    default="swarm")
  parser.add_argument("--plot",action="store_true",
                      help="plotting previous optimised model")
  parser.add_argument("--prettyplot",action="store_true",
                      help="plotting traces for article")
  
  args = parser.parse_args()

  optMethod = "swarm" # args.optMethod
  
  print("Reading file : " + args.file)
  print("Synapse type : " + args.st)
  print("Optimisation method : " + optMethod)

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

  
  logFileName = "logs/" + os.path.basename(args.file) + "-log.txt"
  if(not os.path.exists("logs/")):
    os.makedirs("logs/")
     
  
  
  # "DATA/Yvonne2019/M1RH_Analysis_190925.h5"
  ly = OptimiseSynapses(args.file,synapseType=args.st,dView=dView,role="master",
                        logFileName=logFileName,optMethod=optMethod)
  
  if(args.plot or args.prettyplot):

    if(args.prettyplot):
      prettyPlotFlag = True
    else:
      prettyPlotFlag = False
    
    assert args.id is not None, "You need to specify which trace(s) to plot"
    for id in ly.getUserID(args.id):
      print("Only plotting data for ID " + str(int(id)))

      ly.plotData("GBZ_CC_H20",int(id),show=True,prettyPlot=prettyPlotFlag)
    exit(0)
    
  if(args.id is not None):
    # ID has priority over type
    userID = ly.getUserID(args.id)
    print("Optimising only " + str(userID))
    ly.parallelOptimiseCells("GBZ_CC_H20", userID)

  elif(args.type):
    print("Optimising only " + args.type + " in " + "GBZ_CC_H20")
    ly.parallelOptimiseCells("GBZ_CC_H20",
                             ly.getCellID("GBZ_CC_H20",args.type))
  else:
    ly.parallelOptimiseCells("GBZ_CC_H20")

    
  # (d2,t2) = ly.getData("GBZ_CC_H20",21)
  # xx = ly.getCellID("MSN")



  if(False):
    optHelper = lambda x : ly.optimiseCell( "GBZ_CC_H20", x)

    for idx in ly.getAllCellID():
      try:
        if(not ly.checkValidData("GBZ_CC_H20", idx)):
          # Skip invalid trace
          continue
        ly.optimiseCell("GBZ_CC_H20", idx)
      except:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)
        print("!!! Something went wrong. Skipping neuron " + str(idx))
  
    ly.saveParameterCache()
  
  print("Did we get all the parameters in the file?")
  #import pdb
  #pdb.set_trace()
