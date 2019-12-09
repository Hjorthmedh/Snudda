import numpy as np
import json
import os
import sys
import scipy
import scipy.optimize

import pyswarms as ps
from pyswarms.utils.functions import single_obj as fx
from RunLittleSynapseRun import RunLittleSynapseRun
import matplotlib.pyplot as plt
import matplotlib

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

##############################################################################
    
class Planert2010part2(object):
  
  def __init__(self,dataType,cellID=None,prettyPlot=True,fitData=True):

    self.dataLegend = { "II" : "iSPN to iSPN",
                        "ID" : "iSPN to dSPN",
                        "DD" : "dSPN to dSPN",
                        "DI" : "dSPN to iSPN",
                        "FD" : "FSN to dSPN",
                        "FI" : "FSN to iSPN" }
    
    # Since we reuse code that was parallel, need to say that this is master
    self.role = "master"
    self.figResolution = 300
    self.prettyPlot = prettyPlot

    fileName = "DATA/Planert2010/PlanertFitting-" + dataType + "-cache.json"
    print("Loading data " + str(dataType) + ": " + str(fileName))

    with open(fileName,"r") as f:
      self.data = json.load(f)

    if(cellID is None):
      cellID = np.arange(0,len(self.data["simAmp"]))
      
    # Also load old parameters if they exist?
    self.cacheFileName = "DATA/Planert2010/PlanertFitting-" \
      + dataType + "-tmgaba-fit.json"

    self.loadParameterCache()

    if(fitData):
      # Save parameters after fitting
    
      if(type(cellID) in [list,np.ndarray]):
        for cID in cellID:
          self.fitPeaks(dataType,cID)
      else:
        self.fitPeaks(dataType,int(cellID))

      self.saveParameterCache()
    else:
      self.setupModelSynapseFitting(dataType)
      
      if(type(cellID) in [list,np.ndarray]):
        for cID in cellID:
          self.plotData(dataType,cID,runSim=True,show=False)
      else:
        self.plotData(dataType,int(cellID),runSim=True,show=False)
      

  ############################################################################

  def __delete__(self):
    # Save cache file before exiting
    self.saveParameterCache()
    

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

  def loadParameterCache(self):

    if(os.path.exists(self.cacheFileName)):
      try:
        print("Loading cache file " + str(self.cacheFileName))
        with open(self.cacheFileName,"r") as f:
          tmpDict = json.load(f)
          
          self.parameterCache = dict([])
          for k in tmpDict:
            self.parameterCache[int(k)] = tmpDict[k]

          f.close()
      except:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)
        
        print("Unable to open " + str(self.cacheFileName))
        self.parameterCache = dict([])
    else:
      # No cache file to load, create empty dictionary
      self.parameterCache = dict([])  

  ############################################################################

  def saveParameterCache(self):

    if(self.role != "master"):
      print("No servants are allowed to write output to json, ignoring call.")
      return
    
    print("Saving parameters to cache file: " + str(self.cacheFileName))

    try:
      with open(self.cacheFileName,"w") as f:
        json.dump(self.parameterCache,f,indent=2,cls=NumpyEncoder)
        f.close()
    except:
      import traceback
      tstr = traceback.format_exc()
      print(tstr)
      
      print("Failed to save cache file ... " + str(self.cacheFileName))
      
      
  ############################################################################
  
  def fitPeaks(self,dataType,cellID):

    print("Optimising for ID = " + str(cellID))
    
    # Get peaks
    peakHeight = np.array(self.data["simAmp"][cellID])
    tStim = np.array(self.data["tStim"])

    assert len(tStim) == len(peakHeight), \
      "Inconsistent data lengths for fitPeaks"
    sigma = np.ones(len(peakHeight))
    
    # Setup neuron model
    self.setupModelSynapseFitting(dataType)

    # U, tauR, tauF, tauRatio, cond (obs, tau = tauRatio * tauR)
    modelBounds = ([1e-3,1e-4,1e-4,0, 1e-5],
                   [1.0,2,2,0.9999999,1e-1])

    
    # Pyswarm options
    options = {"c1":0.5, "c2":0.3, "w":0.9} # default

    nParticles = 200
    nIterations = 10

    optimizer = ps.single.GlobalBestPSO(n_particles=nParticles,
                                        dimensions=len(modelBounds[0]),
                                        options=options,
                                        bounds=modelBounds)

    cost, fitParams = optimizer.optimize(self.neuronSynapseSwarmHelperPlanert,
                                         iters=nIterations,
                                         tStim = tStim,
                                         peakHeight = peakHeight)

    modelHeights,tSim,vSim = \
      self._neuronSynapseSwarmHelper(fitParams,tStim)
    
    # tau < tauR, so we use tauRatio for optimisation
    fitParams[3] *= fitParams[1] # tau = tauR * tauRatio

    print("Parameters: U = %.3g, tauR = %.3g, tauF = %.3g, tau = %.3g, cond = %3.g" % tuple(fitParams))
    
    print("peakHeight = " + str(peakHeight))
    print("modelHeights = " + str(modelHeights))
    
    self.addParameterCache(cellID,"synapse", \
                           { "U" : fitParams[0],
                             "tauR" : fitParams[1],
                             "tauF" : fitParams[2],
                             "tau"  : fitParams[3],
                             "cond" : fitParams[4] })

    self.addParameterCache(cellID,"surrogatePeaks",peakHeight)
    self.addParameterCache(cellID,"fittedPeaks",modelHeights)

    self.plotData(dataType=dataType,cellID=cellID,
                  tStim=tStim,
                  surrogatePeaks=peakHeight,
                  modelPeaks=modelHeights,
                  tSim=tSim,
                  vSim=vSim,
                  params=fitParams)
                  

  ############################################################################

  def plotData(self,dataType,cellID,
               tStim=None,
               surrogatePeaks=None,
               modelPeaks=None,
               tSim=None,vSim=None,
               modelParams=None,
               show=True,
               tSkip=0.01,
               prettyPlot=None,
               runSim=False):

    matplotlib.rcParams.update({'font.size': 24})
    
    if(surrogatePeaks is None):
      surrogatePeaks = np.array(self.data["simAmp"][cellID])

    if(tStim is None):
      tStim = np.array(self.data["tStim"])

    if(modelParams is None):
      pDict = self.getParameterCache(cellID,"synapse")
      modelParams = [pDict[x] for x in ["U","tauR","tauF","tau","cond"]]
      
    if(runSim):

      assert modelParams is not None, \
        "plotData: modelParams must be given if runSim = True"
      if(tSim is not None or vSim is not None):
        print("Ignoring tSim and vSim when runSim = True")
      
      print("Run simulation for " + str(dataType) + " " + str(cellID))
      U,tauR,tauF,tau,cond = modelParams

      modelPeaks,tSim,vSim = self.synapseModelNeuron(tStim,U,
                                                     tauR,tauF,cond,tau,
                                                     params={},
                                                     returnTrace=True)
      
    if(prettyPlot is None):
      prettyPlot = self.prettyPlot
    
    plt.figure()
    if(tSim is not None):

      tIdx = np.where(tSim > tSkip)[0][0]
      
      vBase = np.mean(vSim[-50:-1]*1e3)

      for t,sp,mp in zip(tStim*1e3,surrogatePeaks*1e3,modelPeaks*1e3):

        idx = np.argmin(np.abs(tSim-t))
        plt.plot([t,t],[vBase,vBase+sp],color=(1,0.31,0),linewidth=3)
        #plt.plot([t,t],[vBase,vBase+sp],color="red",linewidth=3)        
        plt.plot([t,t],[vBase,vBase+mp],color="blue")

      plt.plot(1e3*tSim[tIdx:],1e3*vSim[tIdx:],color="black")
        
    else:
      for t,sp,mp in zip(tStim*1e3,surrogatePeaks*1e3,modelPeaks*1e3):

        #plt.plot([t,t],[0,sp],color="red",linewidth=3)
        plt.plot([t,t],[0,sp],color=(1,0.31,0),linewidth=3)        
        plt.plot([t,t],[0,mp],color="blue")

      vBase = 0 # No sim data, base it at 0
      
    if(prettyPlot):
      # Draw scalebars
      vScaleX = 1300
      #vMax = np.max(vPlot[np.where(tPlot > 0.050)[0]])
      vScaleY1 = vBase+0.15
      vScaleY2 = vBase+0.05
      tScaleY  = vBase+0.05
      tScaleX1 = 1300
      tScaleX2 = 1400
      
      plt.plot([vScaleX,vScaleX],[vScaleY1,vScaleY2],color="black")
      plt.plot([tScaleX1,tScaleX2],[tScaleY,tScaleY],color="black")
     
        
    if(not prettyPlot):
      if(modelParams is not None):
        titleStr = "\nU=%.3g, tauR=%.3g, tauF=%.3g, tau=%.3g,\ncond=%.3g" \
          % (modelParams[0],
             modelParams[1],
             modelParams[2],
             modelParams[3],
             modelParams[4])
        plt.title(titleStr)

      plt.xlabel("Time (ms)")
      plt.ylabel("Volt (mV)")

    else:
      if(dataType in self.dataLegend):
        plt.title(self.dataLegend[dataType])

      plt.axis("off")

    if(not os.path.exists("figures/")):
      os.makedirs("figures/")

    if(prettyPlot):
      figName = "figures/PlanertFitting-" + dataType \
        + "-" + str(cellID) + "-noaxis.pdf"      
    else:
      figName = "figures/PlanertFitting-" + dataType \
        + "-" + str(cellID) + ".pdf"
    plt.savefig(figName,dpi=self.figResolution)

    if(show):
      plt.ion()
      plt.show()
    else:
      plt.ioff()
      plt.close()
     
    
  ############################################################################

  def neuronSynapseSwarmHelperPlanert(self,pars,tStim,peakHeight):

    res = np.zeros((pars.shape[0]))

    for idx,p in enumerate(pars):
      peakH,tSim,vSim = self._neuronSynapseSwarmHelper(p,tStim)

      # Calculating error in peak height
      hDiff = np.abs(peakH - peakHeight)
      hDiff[0] *= 3
      hDiff[-1] *= 3
      hError = np.sum(hDiff)/len(hDiff)
   
      res[idx] = hError
      
    return res
      
  ############################################################################
    
  def _neuronSynapseSwarmHelper(self,
                               pars,
                               tSpikes):

    U,tauR,tauF,tauRatio,cond = pars
    tau = tauR*tauRatio
    params = {}
    
    peakHeights,tSim,vSim = self.synapseModelNeuron(tSpikes,U,
                                                    tauR,tauF,cond,tau,
                                                    params=params,
                                                    returnTrace=True)
    
    return peakHeights,tSim,vSim
   
  ############################################################################
  
  def setupModelSynapseFitting(self,dataType,params={}):

    # If the delta model existed, clear it
    #self.rsrDeltaModel = None

    somaDiameter = 20e-6
    somaGleak = 3

    params = { "somaDiameter" : somaDiameter, "somaGleak" : somaGleak }

    
    tStim = np.array(self.data["tStim"])
    maxTime = np.max(tStim) + 0.5
    
    baselineDepol = -80e-3

    self.rsrSynapseModel = RunLittleSynapseRun(stimTimes=tStim,
                                               holdingVoltage=baselineDepol,
                                               synapseType="GABA",
                                               params=params,
                                               time=maxTime)

   ############################################################################
  
  def synapseModelNeuron(self,tSpike,U,tauR,tauF,cond,tau,
                         params = {},
                         returnTrace=False):

    # print("Running neuron model")
    
    assert self.rsrSynapseModel is not None, \
      "!!! Need to call setupModelSynapseFitting first"

    # Should we make a copy of params, to not destroy it? ;)
    params["U"]    = U
    params["tauR"] = tauR
    params["tauF"] = tauF
    params["cond"] = cond
    params["tau"]  = tau

    #print("params=" + str(params))
    
    (tSim,vSim,iSim) = \
      self.rsrSynapseModel.run2(pars=params)

    if(tSim.shape != vSim.shape):
      print("Shape are different, why?!")
      import pdb
      pdb.set_trace()            
    
    peakIdx = self.getPeakIdx2(time=tSim,volt=vSim,stimTime=tSpike)
    peakHeight,decayFits,vBase = self.findTraceHeights(tSim,vSim,peakIdx)
    
    if(returnTrace):
      return peakHeight,tSim,vSim
    else:
      return peakHeight

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

      # We assume that neuron is more depolarised than -65, ie gaba is
      # also depolarising
      pIdx = tIdx[np.argmax(volt[tIdx])]
        
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

  # We are using a simplified function that skips decay fits
  
  def findTraceHeights(self,time,volt,peakIdx):

    decayFunc = lambda x,a,b,c : a*np.exp(-x/b) + c
    
    vBase = np.mean(volt[int(0.3*peakIdx[0]):int(0.8*peakIdx[0])])

    peakHeight = np.zeros((len(peakIdx,)))
    peakHeight[0] = volt[peakIdx[0]] - vBase

    decayFits = []

    
    for idxB in range(1,len(peakIdx)):

      if(peakHeight[0] > 0):
        if(idxB < len(peakIdx) -1):        
          p0d = [0.06,-0.05,-0.074]
        else:
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
      
        idxStart = int(peakIdxA*0.2 + peakIdxB*0.8)
        idxEnd = int(peakIdxA*0.08 + peakIdxB*0.92)
      else:
        # Last spike, use only last half of decay trace
        idxStart = int(peakIdxA*0.6 + peakIdxB*0.4)
        idxEnd = int(peakIdxA*0.08 + peakIdxB*0.92)

      try:
        assert idxStart < idxEnd
      except:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)

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
          print(tstr)
          
          print("!!! Failed to converge, trying with smaller decay constant")
          p0d[1] *= 0.01
          fitParams,pcov = scipy.optimize.curve_fit(decayFunc,tABfit,vABfit,
                                                    p0=p0d)
          
        tB = time[peakIdxB] - tAB[0]
        vBaseB = decayFunc(tB,fitParams[0],fitParams[1],fitParams[2])
        
        peakHeight[idxB] = volt[peakIdxB] - vBaseB

        vFit = decayFunc(tAB-tAB[0],fitParams[0],fitParams[1],fitParams[2])
        decayFits.append((tAB,vFit))

        ################################################################

        if(False):
          plt.figure()
          plt.plot(tAB,vAB,'r')
          plt.title("Error in findTraceHeights")
          plt.xlabel("time")
          plt.ylabel("volt")
          plt.plot(tAB,vFit,'k-')
          plt.ion()
          plt.show()

          import pdb
          pdb.set_trace()

        
        ########################################

      except:
        
        print("Check that the threshold in the peak detection before is OK")
        #self.plot(name)
        import traceback
        tstr = traceback.format_exc()
        print(tstr)

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
  
    
  ############################################################################

if __name__ == "__main__":

  import argparse
  parser = argparse.ArgumentParser(description="Fit GABA model to peaks")
  parser.add_argument("dataType",choices=["DD","ID","DI","II","FI","FD"])
  parser.add_argument("--plotOnly",help="Only plot, no new fitting",
                      action="store_true")
  args = parser.parse_args()

  if(args.plotOnly):
    fitData = False
  else:
    fitData = True
    
  pp = Planert2010part2(dataType=args.dataType,cellID=None,
                        prettyPlot=True,fitData=fitData)
  #pp = Planert2010part2(dataType=args.dataType,cellID=0)
          
