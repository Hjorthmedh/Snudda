import os
import re
import numpy as np
import h5py
import scipy
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib
import neuron
import json
import copy
import time

# TODO 2020-10-19
#
# We need to make sure params contains the nmda_ratio and other parameters
# that are important for mod file (but not optimised for at this stage)
# Please fix.
#
#
# python3 OptimiseSynapsesFull.py DATA/Yvonne2020/M1RH-ipsi_MSN_D1_20Hz_depressing.json --synapseParameters ../data/synapses/v3/M1RH-ipsi_D1-MSN.json --st glut
#
#
# TODO 2020-10-09
#
# We are no longer reading data from Yvonne's HDF5 directly, instead we are
# reading data from a JSON file, that Ilaria created from Igor exports.
#
# Consequences and TODO:
# - We no longer need CellID to identify which cell we are optmising,
#   each JSON file only contains one dataset (and it is an averaged dataset)
# - Holding voltage is no longer extractable from data, we need to set it
# - Create a JSON optmisation parameter file which contains the holding voltage
#   to be used, as well as the modelbounds (currently in getModelBounds)
# - The JSON data is now loaded into self.data, go through all functions
#   remove references to CellID, and extract the data directly from the
#   self.data variable.
#
#

# TODO 2020-07-15
#
# Make it so that the code can run in serial mode also, to simplify debugging
#

#
# TODO 2020-07-02
# -- We just wrote parallelOptimiseSingleCell -- need to make sure we call it
#    the function will optimise one cellID, using all workers available
#    need to debug the code, to make sure it works... have fun! :)
#

#
# TODO 2020-06-16
# -- We need to make sure one neuron can be optimised by all workers
#    collectively, right now one worker does one cell alone
#
# -- Determine synapse locations, pass it to all workers
# -- Best random is easy to parallelise, just do the work, then gather it at
#    the master node.
# -- Difficult: How to parallise scipy.optimize.minimize
#    Possible idea: let func to optimize handle vectors, and then each
#    position in vector is sent to one worker.

#
#
# TODO (depricated):
# 1. Remove the old soma optimisation code, not needed anymore
# 2. Place the inital synapses, then find out the sectionID,X,
# 3. Pass sectionID, sectionX to all workers, so they use same synapselocations
# 4. Optimise.
#

# !!! Add check that if the voltage is 0 or 5, then the trace is skipped entirely

from RunSynapseRun import RunSynapseRun

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
          
class OptimiseSynapsesFull(object):

  ############################################################################

  # optMethod not used anymore, potentially allow user to set sobol or refine

  # datafile is JSON file from Ilaria's Igor extraction
  
  def __init__(self, datafile, synapseType="glut",loadCache=True,
               role="master",dView=None,verbose=True,logFileName=None,
               optMethod="sobol",prettyPlot=False,
               model_bounds="model_bounds.json",
               neuronSetFile="neuronSet.json",
               synapseParameterFile=None):

    # Parallel execution role, "master" or "servant"
    self.role = role
    
    self.parallelSetupFlag = False # Set to True when servants are done setup
    self.dView = dView
    self.verbose = verbose
    self.logFileName = logFileName
    self.optMethod = optMethod
    self.nSmoothing = 200 # How many smoothing points do we use?
    self.simTime = 1.8
    self.neuronSetFile = neuronSetFile
    
    self.debugParsFlag = False
    self.debugPars = []
    self.cellProperties = None

    self.prettyPlot = prettyPlot
    
    print("Init optMethod = " + str(optMethod))
    
    if(self.logFileName is not None and len(self.logFileName) > 0):
      print("Log file: " + self.logFileName)
      self.logFile = open(self.logFileName,'w')
    else:
      self.logFile = None
      
    self.figResolution = 300

    self.datafile = datafile

    self.writeLog("Loading " + str(datafile))
    with open(datafile,"r") as f:
      self.data = json.load(f)
      
      self.volt = np.array(self.data["data"]["mean_norm_trace"])
      self.sampleFreq = self.data["metadata"]["sample_freq"]
      
      dt = 1/self.sampleFreq
      self.time = 0 + dt * np.arange(0,len(self.volt))

      self.stimTime = np.array(self.data["metadata"]["stim_time"]) * 1e-3 # ms
      
      self.cellType = self.data["metadata"]["cell_type"]

    self.synapseParameterFile = synapseParameterFile
    
    if synapseParameterFile:
      with open(synapseParameterFile,'r') as f:
        print(f"Reading synapse parameters from {synapseParameterFile}")
        tmp = json.load(f)
        self.synapseParameters = tmp["data"]
    else:
      self.synapseParameters = {}
      
    self.cacheFileName = str(self.datafile) + "-parameters-full.json"
    self.loadCache = loadCache
    self.synapseType = synapseType

    self.rsrSynapseModel = None
    self.rsrDeltaModel = None

    self.modelInfo = None

    with open(model_bounds,'r') as f:
      print(f"Loading model bounds from {model_bounds}")
      self.model_bounds = json.load(f)
    
    if(loadCache):
      self.loadParameterCache()
    else:
      self.parameterCache = dict([])


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
  
  def addParameterCache(self,name,value):

    self.parameterCache[name] = value

  ############################################################################

  def getParameterCache(self,name):

    if name in self.parameterCache:
      return self.parameterCache[name]
    else:
      return None
    

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


  
  def plotData(self,
               params=None,
               show=True,
               skipTime=0.3,
               prettyPlot=None):
      

    if params is None:
      params = self.synapseParameters
    
    if(prettyPlot is None):
      prettyPlot = self.prettyPlot

    if(prettyPlot):
      matplotlib.rcParams.update({'font.size': 24})
    else:
      matplotlib.rcParams.update({'font.size': 5})
      
    if(self.volt is None):
      self.writeLog(dataType + " " + str(cellID) + ": Nothing to plot")
      return


    
    bestParams = self.getParameterCache("param")
    synapsePositionOverride = (self.getParameterCache("sectionID"),
                               self.getParameterCache("sectionX"))
    minError = self.getParameterCache("error")

    vPlot = None

    if(bestParams is not None):
      U, tauR, tauF, tauRatio, cond = bestParams

      params = { "U" : U,
                 "tauR" : tauR,
                 "tauF" : tauF,
                 "cond" : cond,
                 "tau" : tauR * tauRatio}
    
    

      plotModel = self.setupModel(params=params,
                                  synapsePositionOverride \
                                    = synapsePositionOverride)
    
      (tPlot,vPlot,iPlot) = plotModel.run2(pars=params)

    tIdx = np.where(skipTime <= t)[0]
        
    plt.figure()

    plt.plot(self.time[tIdx]*1e3,self.data[tIdx]*1e3,'r-')
    if(vPlot is not None):
      t2Idx = np.where(skipTime <= tPlot)[0]
      plt.plot(tPlot[t2Idx]*1e3,vPlot[t2Idx]*1e3,'k-')
    # plt.title(dataType + " " + str(cellID))

    if(not prettyPlot):
      titleStr = self.cellType
      
      titleStr += "\nU=%.3g, tauR=%.3g, tauF=%.3g, tau=%.3g,\ncond=%.3g" \
          % (params["U"],
             params["tauR"],
             params["tauF"],
             params["tau"],
             params["cond"])
        

      plt.title(titleStr)

    if(prettyPlot):
      # Draw scalebars
      vScaleX = 1200
      #vMax = np.max(vPlot[np.where(tPlot > 0.050)[0]])
      vBase = vPlot[-1]
      yScaleBar = vBase*1e3 + float(np.diff(plt.ylim()))/4
      vScaleY1 = yScaleBar+1
      vScaleY2 = yScaleBar
      tScaleY  = yScaleBar
      tScaleX1 = vScaleX
      tScaleX2 = vScaleX + 100
      
      plt.plot([vScaleX,vScaleX],[vScaleY1,vScaleY2],color="black")
      plt.plot([tScaleX1,tScaleX2],[tScaleY,tScaleY],color="black")

      plt.text(vScaleX-100, vScaleY2 + 0.20*float(np.diff(plt.ylim())),\
               ("%.0f" % (vScaleY1-vScaleY2)) + " mV",
                     rotation=90)
      plt.text(vScaleX,vScaleY2 - float(np.diff(plt.ylim()))/10,
               ("%.0f" % (tScaleX2 - tScaleX1) + " ms"))
      
      # Mark optogenetical stimulation
      yHeight = float(np.diff(plt.ylim()))/13
      
      tStim = self.getStimTime()*1e3
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
    
    baseName = os.path.splitext(os.path.basename(self.datafile))[0]
    figName = "figures/" + baseName + ".pdf"
    plt.savefig(figName,dpi=self.figResolution)

    if(show):
      plt.ion()
      plt.show()
    else:
      plt.ioff()
      plt.close()
      
    
  ############################################################################

  def getCellProperties(self):

    if(self.cellProperties is None):
      with open(self.neuronSetFile,'r') as f:
        self.cellProperties = json.load(f)

    cellType = self.data["metadata"]["cell_type"]

    return self.cellProperties[cellType].copy()

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

  
  def getPeakIdx(self):

    pTime = np.array(self.data["metadata"]["stim_time"])*1e-3
    freq = self.data["metadata"]["freq"]

    assert np.abs(1.0-freq/(pTime[1]-pTime[0])) < 0.01, "frequency mismatch"
    
    pWindow = 1.0/(2*freq)*np.ones(pTime.shape)
    pWindow[-1] *= 5
    
    peakInfo = self.findPeaksHelper(pTime=pTime,
                                    pWindow=pWindow,
                                    time=self.time,
                                    volt=self.volt)

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

  # Find peaks within pStart[i] and pStart[i]+pWindow[i]
  # The value is not the amplitude of the peak, just the voltage at the peak
  
  def findPeaksHelper(self,pTime,pWindow,
                      time=None,volt=None):
      
    peakIdx = []
    peakTime = []
    peakVolt = []
    
    for pt,pw in zip(pTime,pWindow):
      tStart = pt
      tEnd = pt + pw

      tIdx = np.where(np.logical_and(tStart <= time,time <= tEnd))[0]
      assert len(tIdx) > 0, f"No time points within {tStart} and {tEnd}" 
      
      if(self.synapseType == "glut"):
        pIdx = tIdx[np.argmax(volt[tIdx])]
      elif(self.synapseType == "gaba"):
        # We assume that neuron is more depolarised than -65, ie gaba is
        # also depolarising
        pIdx = tIdx[np.argmax(volt[tIdx])]
      else:
        self.writeLog("Unknown synapse type : " + str(self.synapseType))
        import pdb
        pdb.set_trace()
        
      peakIdx.append(int(pIdx))
      peakTime.append(time[pIdx])
      peakVolt.append(volt[pIdx])
      
    # Save to cache -- obs peakVolt is NOT amplitude of peak, just volt

    peakDict = { "peakIdx" : np.array(peakIdx),
                 "peakTime" : np.array(peakTime),
                 "peakVolt" : np.array(peakVolt)} # NOT AMPLITUDE

#    if(cellID is not None):
#      self.addParameterCache("peaks",peakDict)
                           
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
        idxEnd = int(peakIdxA*0.05 + peakIdxB*0.85) # might need 0.85 as last

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

  def setupModel(self,params={},
                 synapseDensityOverride=None,
                 nSynapsesOverride=None,
                 synapsePositionOverride=None):
    
    tStim = self.stimTime

    # Read the info needed to setup the neuron hosting the synapses
    cProp = self.getCellProperties()

    if(synapsePositionOverride is not None):
      synapseSectionID,synapseSectionX = synapsePositionOverride
    else:
      synapseSectionID,synapseSectionX = None,None
      
    if(synapseDensityOverride is not None):
      synapseDensity = synapseDensityOverride
    else:
      synapseDensity = cProp["synapseDensity"]
      
    if(nSynapsesOverride is not None):
      nSynapses = nSynapsesOverride
    else:
      nSynapses = cProp["nSynapses"]
      
    # !!! We need to get the baseline depolarisation in another way

    self.rsrSynapseModel = \
      RunSynapseRun(neuronMorphology=cProp["neuronMorphology"],
                    neuronMechanisms=cProp["neuronMechanisms"],
                    neuronParameters=cProp["neuronParameters"],
                    neuronModulation=cProp["neuronModulation"],
                    stimTimes=tStim,
                    nSynapses=nSynapses,
                    synapseDensity=synapseDensity,
                    holdingVoltage=cProp["baselineVoltage"],
                    synapseType=self.synapseType,
                    params=params,
                    time=self.simTime,
                    logFile=self.logFile,
                    synapseSectionID=synapseSectionID,
                    synapseSectionX=synapseSectionX)


    return self.rsrSynapseModel
    
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

    U,tauR,tauF,tauRatio,cond = pars
    tau = tauR*tauRatio
    
    peakHeights,tSim,vSim = self.runModel(tSpikes,U,
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

    peakHeights = self.runModel(tSpike,U,tauR,tauF,cond,tau)
    
    return peakHeights

   ############################################################################

   # !!! OBS tauRatio is inparameter
   
  def neuronSynapseHelperGlut(self,tSpike, 
                              U, tauR, tauF, tauRatio, cond,
                              smoothExpTrace8, smoothExpTrace9, expPeakHeight,
                              returnType="peaks"):

    if(self.debugParsFlag):
      self.debugPars.append([U, tauR, tauF, tauRatio, cond])

    params = self.synaseParameters
    tau = tauR * tauRatio
    
    peakH,tSim,vSim = self.runModel(tSpike,U,
                                    tauR,tauF,cond,tau,
                                    params=params,
                                    returnTrace=True)
    
    
    # Calculate error in decay fit
    simTrace8,simTime8 = self.smoothingTrace(vSim,self.nSmoothing,
                                             time=tSim,
                                             startTime=self.decayStartFit8,
                                             endTime=self.decayEndFit8)

    simTrace9,simTime9 = self.smoothingTrace(vSim,self.nSmoothing,
                                             time=tSim,
                                             startTime=self.decayStartFit9,
                                             endTime=self.decayEndFit9)

    # We only want to use the bit of the trace after max
    idxMax8 = np.argmax(smoothExpTrace8)
    idxMax9 = np.argmax(smoothExpTrace9)    
    
    # Calculating error in peak height
    hDiff = np.abs(peakH - expPeakHeight)
    hDiff[0] *= 3
    hDiff[-2] *= 2    
    hDiff[-1] *= 3

    # This is to prevent the model spiking
    spikePenalty = np.sum(peakH > 0.03)*1
    
    hError = np.sum(hDiff)/len(hDiff)

    decayError8 = np.sum((smoothExpTrace8[idxMax8:] \
                          - simTrace8[idxMax8:])**2) \
                         /(self.nSmoothing-idxMax8+1) * 10000

    decayError9 = np.sum((smoothExpTrace9[idxMax9:] \
                         - simTrace9[idxMax9:])**2) \
                         /(self.nSmoothing-idxMax9+1) * 10000

    fitError = hError + decayError8 + decayError9 + spikePenalty
      

    if(spikePenalty > 0):
      self.writeLog("Action potential detected in trace. Penalising!")
      
    
    if(False):
      peakBase = vSim[-1]
      plt.figure()
      plt.plot(tSim,vSim,'k-')
      plt.plot(simTime8,simTrace8,'y--')
      plt.plot(simTime8,smoothExpTrace8,'r--')
      plt.plot(simTime9,simTrace9,'y--')
      plt.plot(simTime9,smoothExpTrace9,'r--')
      
      for tp,expH,modH in zip(tSpike,expPeakHeight,peakH):
        plt.plot([tp,tp],[peakBase,expH+peakBase],'r-',linewidth=3)
        plt.plot([tp,tp],[peakBase,modH+peakBase],'b-')
      plt.title("hE = %g, dE8 = %g, dE9 = %g" \
                % (hError,decayError8,decayError9))
        
      plt.ion()
      plt.show()
      


    if(returnType == "peaks"):
      return peakH
    elif(returnType == "error"):
      return fitError
    elif(returnType == "full"):
      return fitError,peakH,tSim,vSim
    else:
      assert False, "Unknown return type: " + str(returnType)
    
   ############################################################################
  
  def runModel(self,tSpike,U,tauR,tauF,cond,tau,
               params = {},
               returnTrace=False):

    # self.writeLog("Running neuron model")
    
    assert self.rsrSynapseModel is not None, \
      "!!! Need to call setupModel first"

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

  # This should read from a JSON file instead
  
  def getModelBounds(self):
    
    cellType = self.data["metadata"]["cell_type"]

    mb = self.model_bounds[cellType]

    paramList = ["U", "tauR", "tauF", "tauRatio", "cond"]
    lowerBound = [mb[x][0] for x in paramList]
    upperBound = [mb[x][1] for x in paramList]

    return (lowerBound,upperBound)
    
  
  ############################################################################

  def sobolScan(self,synapseModel,
                tStim,hPeak,
                modelBounds,
                smoothExpTrace8, smoothExpTrace9,
                nTrials=6,loadParamsFlag=False,
                parameterSets = None,
                returnMinError=False):

    assert self.synapseType == "glut", \
      "GABA synapse not supported yet in new version"

    if(parameterSets is None):
      parameterSets = self.setupParameterSet(modelBounds,nTrials)
    

    # zip(*xxx) unzips xxx -- cool.
    USobol,tauRSobol,tauFSobol,tauRatioSobol,condSobol \
      = zip(*parameterSets)
    
    # tauSobol = np.multiply(tauRatioSobol,tauRSobol)

    minPars = None
    minError = np.inf
      
    if(loadParamsFlag):
      # If we should load params then do so first
      minPars = self.getParameterCache("synapse")

      # --- now parameters are read from cache, but we can in future
      # have them read from a work-log with parameters to do etc
      
      # What was error of the cached parameterset
      if(minPars is not None):
        minError = self.neuronSynapseHelperGlut(tStim,
                                                U=minPars[0],
                                                tauR=minPars[1],
                                                tauF=minPars[2],
                                                tauRatio=minPars[3]/minPars[1],
                                                cond=minPars[4],
                                                smoothExpTrace8=smoothExpTrace8,
                                                smoothExpTrace9=smoothExpTrace9,
                                                expPeakHeight=hPeak,
                                                returnType="error")

    idx = 0
        
    for U,tauR,tauF,tauRatio,cond \
        in zip(USobol, tauRSobol, tauFSobol, tauRatioSobol, \
               condSobol):

      idx += 1
      if(idx % 50 == 0):
        self.writeLog("%d / %d : minError = %g" % (idx, len(USobol),minError))
        self.writeLog(str(minPar))
   
      error = self.neuronSynapseHelperGlut(tStim,U,tauR,tauF,tauRatio,
                                           cond,
                                           smoothExpTrace8=smoothExpTrace8,
                                           smoothExpTrace9=smoothExpTrace9,
                                           expPeakHeight=hPeak,
                                           returnType="error")
      try:      
        if(error < minError):
          minError = error
          minPar = np.array([U,tauR,tauF,tauRatio,cond])

          # TODO, write intermediate results to file, in case of a crash...
          
      except:
        import traceback
        tstr = traceback.format_exc()
        self.writeLog(tstr)

        import pdb
        pdb.set_trace()
        
        # For big runs we do no want to give up. Let's try again...
        continue
    

    if(returnMinError):
      return minPar,minError
    else:
      return minPar
      
  ############################################################################
  
  def bestRandom(self,synapseModel,
                 tStim,hPeak,
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
        pars = self.getParameterCache("synapse")
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
        tauR = np.random.uniform(modelBounds[0][1],modelBounds[1][1])
        tauF = np.random.uniform(modelBounds[0][2],modelBounds[1][2])
        tau = tauR * np.random.uniform(modelBounds[0][3],modelBounds[1][3])
        cond = np.random.uniform(modelBounds[0][4],modelBounds[1][4])

      try:
        peakHeights = self.runModel(tStim,U,tauR,tauF,tau,cond)
      
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

  def optimiseCell(self):

    assert False, "obsolete"
    
    try:
      # self.fitCellProperties(dataType,cellID)
      self.fitTrace()
      self.plotData(show=False)

      # We need to extract the dictionary associated with cellID and
      # return the result so that the main function can plug it into the
      # central parameter cache -- OR DO WE?

      return self.getSubDictionary(cellID)
    
    except:
      
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      self.addParameterCache("error",tstr)
      

  ############################################################################

  def parallelOptimiseSingleCell(self,nTrials=10000, \
                                 postOpt=False):

    # !!! Future improvement. Allow continuation of old optimisation by
    # reading synapse location and old parameter set, so that is not thrown away
    
    if(self.role == "master"):

      # 1. Setup workers
      params = self.synapseParameters

      # 2. Setup one cell to optimise, randomise synapse positions
      synapseModel = self.setupModel(params=params)

      #(volt,time) = self.getData(dataType,cellID)
      peakIdx = self.getPeakIdx2(stimTime=self.stimTime,
                                 time=self.time,
                                 volt=self.volt)
      tSpikes = self.time[peakIdx]
    
      sigma = np.ones(len(peakIdx))
      sigma[-1] = 1./3

      peakHeight,decayFits,vBase = self.findTraceHeights(self.time,
                                                         self.volt,
                                                         peakIdx)
      
      # 2b. Create list of all parameter points to investigate
      modelBounds = self.getModelBounds()
      parameterPoints = self.setupParameterSet(modelBounds,nTrials)
      
      # 3. Send synapse positions to all workers, and split parameter points
      #    between workers

      if(self.dView is not None):
        self.setupParallel(self.dView)

        self.dView.scatter("parameterPoints",parameterPoints,block=True)
      
        self.dView.push( { "params"           : params,
                           "synapseSectionID" : synapseModel.synapseSectionID,
                           "synapseSectionX"  : synapseModel.synapseSectionX,
                           "modelBounds"      : modelBounds,
                           "stimTime"         : self.stimTime,
                           "peakHeight"       : peakHeight },
                         block=True)

        cmdStrSetup = \
          "ly.sobolWorkerSetup(params=params," \
          + "synapsePositionOverride=(synapseSectionID,synapseSectionX))"

        self.dView.execute(cmdStrSetup,block=True)
      
        cmdStr = "res = ly.sobolScan(synapseModel=ly.synapseModel, \
                                     tStim = stimTime, \
                                     hPeak = peakHeight, \
                                     parameterSets=parameterPoints, \
                                     modelBounds=modelBounds, \
                                     smoothExpTrace8=ly.smoothExpVolt8, \
                                     smoothExpTrace9=ly.smoothExpVolt9, \
                                     returnMinError=True)"

        self.writeLog("Executing workers, bang bang")
        self.dView.execute(cmdStr,block=True)
      
        # 5. Gather worker data
        self.writeLog("Gathering results from workers")
        res = self.dView["res"]

        # * unpacks res variable
        parSets,parError = zip(*res)

        minErrorIdx = np.argsort(parError)

        import pdb
        pdb.set_trace()
        
        # We save parameter set, synapse locations, error value
        bestPar = (parSets[minErrorIdx[0]],
                   (synapseModel.synapseSectionID,synapseModel.synapseSectionX),
                   parError[minErrorIdx[0]])

        
      else:

        # No dView, run in serial mode...
        # !!!
        self.sobolWorkerSetup(params=params,
                              synapsePositionOverride \
                                = (synapseModel.synapseSectionID,
                                   synapseModel.synapseSectionX))

        parSet,parError = self.sobolScan(synapseModel = synapseModel,
                                         tStim = self.stimTime, \
                                         hPeak = peakHeight, \
                                         modelBounds=modelBounds, \
                                         smoothExpTrace8=ly.smoothExpVolt8, \
                                         smoothExpTrace9=ly.smoothExpVolt9, \
                                         returnMinError=True)

        
        
        
        bestPar = (parSet,
                   (synapseModel.synapseSectionID,synapseModel.synapseSectionX),
                   parError)
        
 
      self.addParameterCache("param", bestPar[0])
      self.addParameterCache("sectionID", synapseModel.synapseSectionID)
      self.addParameterCache("sectionX", synapseModel.synapseSectionX)
      self.addParameterCache("error", bestPar[2])
      

      self.saveParameterCache()
      self.writeLog(f"Sobol search done. Best parameter {bestPar}")

      if(postOpt):
        # This updates parameters and saves new parameter cache
        self.getRefinedParameters()
        
      
  ############################################################################

  # This sets up the model also, so can be run in a self-contained loop
  # We might later want to let the workers do this, but then they cant
  # write to cache --- THAT WILL LEAD TO DATA CORRUPTION!!
  
  def getRefinedParameters(self):

    assert self.role == "master", \
      "You do not want to run this on workers in parallel, " \
      + " since it writes directly to parameter cache. " \
      + " That could lead to corrupted data."
    
    # Load parameters from disk
    self.loadParameterCache()
    modelBounds = self.getModelBounds()

    startPar = self.getParameterCache("param")
    sectionX = self.getParameterCache("sectionX")
    sectionID = self.getParameterCache("sectionID")    
    startParErrorVal = self.getParameterCache("error")

    synapsePositionOverride = (sectionID,sectionX)

    # Make sure we have correct taus etc for synapse
    params = self.synapseParameters
    
    peakIdx = self.getPeakIdx2(stimTime=self.stimTime,
                               time=self.time,
                               volt=self.volt)
    tSpikes = time[peakIdx]

    
    peakHeight,decayFits,vBase = self.findTraceHeights(self.time,
                                                       self.volt,
                                                       peakIdx)
    
    self.sobolWorkerSetup(params, \
                          synapsePositionOverride = synapsePositionOverride)
    
    func = lambda x : \
           self.neuronSynapseHelperGlut(tSpike=self.stimTime,
                                        U=x[0],
                                        tauR=x[1],
                                        tauF=x[2],
                                        tauRatio=x[3],
                                        cond=x[4],
                                        smoothExpTrace8=self.smoothExpVolt8,
                                        smoothExpTrace9=self.smoothExpVolt9,
                                        expPeakHeight=peakHeight,
                                        returnType="error")

    mBounds = [x for x in zip(modelBounds[0],modelBounds[1])]
    startPar = self.getParameterCache("param")
    
    res = scipy.optimize.minimize(func,
                                  x0=startPar,
                                  bounds=mBounds)

    fitParams = res.x
    minError = res.fun

    if(minError >= startParErrorVal):
      print("Refinement failed. Sobol parameters are better match than new fitting")
      # Dont overwrite the old parameters
      
    else:
      self.addParameterCache("param", fitParams)
      self.addParameterCache("error", minError)

      self.saveParameterCache()

      print(f"Old error: {startParErrorVal}, New error: {minError}")
                    
  ############################################################################

  def sobolWorkerSetup(self,params,
                       synapsePositionOverride=None):

    self.synapseModel = self.setupModel(params=params,
                                        synapsePositionOverride \
                                            = synapsePositionOverride)

    self.decayStartFit8 = 0.45
    self.decayEndFit8   = 0.8

    self.decayStartFit9 = 1.0
    self.decayEndFit9   = 1.3

    self.smoothExpVolt8,self.smoothExpTime8 \
      = self.smoothingTrace(self.volt,self.nSmoothing,
                            time=self.time,
                            startTime = self.decayStartFit8,
                            endTime = self.decayEndFit8)
        
    self.smoothExpVolt9,self.smoothExpTime9 \
      = self.smoothingTrace(self.volt,self.nSmoothing,
                            time=self.time,
                            startTime = self.decayStartFit9,
                            endTime = self.decayEndFit9)
    
    
      
  ############################################################################

  def setupParameterSet(self,modelBounds,nSets):
    
    import chaospy
    distribution = chaospy.J(chaospy.Uniform(modelBounds[0][0],
                                             modelBounds[1][0]),
                             chaospy.Uniform(modelBounds[0][1],
                                             modelBounds[1][1]),
                             chaospy.Uniform(modelBounds[0][2],
                                             modelBounds[1][2]),
                             chaospy.Uniform(modelBounds[0][3],
                                             modelBounds[1][3]),
                             chaospy.Uniform(modelBounds[0][4],
                                             modelBounds[1][4]))
    
    USobol,tauRSobol,tauFSobol,tauRatioSobol,condSobol \
      = distribution.sample(nSets, rule="sobol")

    parameterSets = [x for x in zip(USobol, \
                                    tauRSobol,tauFSobol,tauRatioSobol, \
                                    condSobol)]
    
    return parameterSets
      
  ############################################################################

  def parallelOptimiseCells(self,dataType,cellIDlist=None):

    assert False, "Only doing one cell at a time"
    
    if(self.role == "master"):

      if(cellIDlist is None):
        cellIDlist = self.getValidCellID(dataType)

      # self.printCellProperties(dataType,cellIDlist)
        
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
      from RunSynapseRun import RunSynapseRun
      from OptimiseSynapsesFull import NumpyEncoder
      from OptimiseSynapsesFull import OptimiseSynapsesFull

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
    
    self.dView.push({"datafile": self.datafile,
                     "synapseType" : self.synapseType,
                     "synapseparameters" : self.synapseParameterFile,
                     "loadCache" : self.loadCache,
                     "role" : "servant"})

    cmdStr = "ly = OptimiseSynapsesFull(datafile=datafile, synapseParameterFile=synapseparameters, synapseType=synapseType,loadCache=loadCache,role=role,logFileName=engineLogFile[0])"
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

      for ctr,par in enumerate(self.debugPars):
        Uall[:,ctr]      = par[:,0]
        tauR[:,ctr]      = par[:,1]
        tauF[:,ctr]      = par[:,2]
        tauRatio[:,ctr]  = par[:,3]
        cond[:,ctr]      = par[:,4]
      
      
      plt.figure()
      plt.plot(Uall,cond,'-')
      plt.xlabel("U")
      plt.ylabel("cond")

      plt.figure()
      plt.plot(tauR,tauF,'-')
      plt.xlabel("tauR")
      plt.ylabel("tauF")

      plt.figure()
      plt.plot(Uall)
      plt.xlabel("Uall")
      
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
  parser.add_argument("datafile", help="JSON DATA file")
  parser.add_argument("--synapseParameters", help="Static synapse parameters (JSON)",
                      default=None)
  parser.add_argument("--st", help="Synapse type (glut or gaba)",
                      choices=["glut","gaba"])
  parser.add_argument("--optMethod",
                      help="Optimisation method",
                      choices=["sobol","stupid","swarm"],
                      default="sobol")
  parser.add_argument("--plot",action="store_true",
                      help="plotting previous optimised model")
  parser.add_argument("--prettyplot",action="store_true",
                      help="plotting traces for article")
  
  args = parser.parse_args()

  optMethod = args.optMethod
  
  print("Reading file : " + args.datafile)
  print("Synapse type : " + args.st)
  print("Synapse params :" + args.synapseParameters)
  print("Optimisation method : " + optMethod)

  print("IPYTHON_PROFILE = " + str(os.getenv('IPYTHON_PROFILE')))
  
  if(os.getenv('IPYTHON_PROFILE') is not None \
     or os.getenv('SLURMID') is not None):
    from ipyparallel import Client
    rc = Client(profile=os.getenv('IPYTHON_PROFILE'),
                debug=False)

    # http://davidmasad.com/blog/simulation-with-ipyparallel/
    # http://people.duke.edu/~ccc14/sta-663-2016/19C_IPyParallel.html
    dView = rc.direct_view(targets='all') # rc[:] # Direct view into clients
    lbView = rc.load_balanced_view(targets='all')
  else:
    dView = None

  
  logFileName = "logs/" + os.path.basename(args.datafile) + "-log.txt"
  if(not os.path.exists("logs/")):
    os.makedirs("logs/")
     
  
  
  # "DATA/Yvonne2019/M1RH_Analysis_190925.h5"
  ly = OptimiseSynapsesFull(datafile=args.datafile,
                            synapseParameterFile=args.synapseParameters,
                            synapseType=args.st,dView=dView,
                            role="master",
                            logFileName=logFileName,optMethod=optMethod)
  
  if(args.plot or args.prettyplot):

    if(args.prettyplot):
      prettyPlotFlag = True
    else:
      prettyPlotFlag = False

    ly.plotData(show=True,prettyPlot=prettyPlotFlag)

    exit(0)
    
    
  ly.parallelOptimiseSingleCell(nTrials=12)

