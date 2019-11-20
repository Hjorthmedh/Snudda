# Plan:
#
# Randomise dtau, ftau, U -- check that it is within the ranges specified
# for pair pulse ratio, and recovery test ratio
#
# Fit the neuron synapse model to the modelled traces
#
#

import json
import numpy as np
import scipy.stats as stats

import matplotlib.pyplot as plt

################################################################################

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

################################################################################    

class Planert2010(object):


  def __init__(self,parType,num=10):
    

    self.defExpData()
    pars= self.pickRandomParam(parType=parType,num=num)
    amps = None

    if(True):
      # Just checking that original starting distributions are ok
      self.plotModelParams(parType,pars)
    
    amps = self.checkDataMatch(parType=parType,pars=pars,amps=amps)
    
    for ctr in range(0,3):
      pars,amps = self.pruneData(parType=parType,pars=pars,amps=amps)
      self.checkDataMatch(parType=parType,pars=pars,amps=amps)

    self.saveData(parType,pars,amps)

    self.plotModelParams(parType,pars)
      

  ############################################################################

  def pruneData(self,parType,pars,amps=None):

    if(amps.shape[0] == 0):
      print("No data left?")
      import pdb
      pdb.set_trace()
    
    pprMean,pprStd = self.conInfo[parType]["ppr"]
    rtrMean,rtrStd = self.conInfo[parType]["rtr"]    

    if(amps is None):
      amps = self.synapseModelWrapper(self.tStim,pars)

    ppr = np.divide(amps[:,1],amps[:,0])
    rtr = np.divide(amps[:,-1],amps[:,0])

    nBins = 20
    pprMax = np.max(ppr)
    rtrMax = np.max(rtr)
        
    # Which bin does each ppr belong to
    pprIdx = (nBins*ppr / (pprMax+1e-6)).astype(int)
    rtrIdx = (nBins*rtr / (rtrMax+1e-6)).astype(int)

    pprCount = np.zeros((nBins,))
    rtrCount = np.zeros((nBins,))

    for p,r in zip(pprIdx,rtrIdx):
      pprCount[p] += 1
      rtrCount[r] += 1
    
    pprBinWidth = pprMax/(nBins-1)
    rtrBinWidth = rtrMax/(nBins-1)
    
    pprCentre = pprBinWidth/2 + np.arange(0,nBins)*pprBinWidth
    rtrCentre = rtrBinWidth/2 + np.arange(0,nBins)*rtrBinWidth

    pprDensity = stats.norm.pdf(pprCentre,pprMean,pprStd)
    rtrDensity = stats.norm.pdf(rtrCentre,rtrMean,rtrStd)

    # We pick random param sets, and see if that PPR and RTR are over
    # represented, if so we remove it.

    nLeft = len(pars["u"])
    keepFlag = np.ones((nLeft,),dtype=bool)
    done = False

    idxRand = np.random.permutation(np.where(keepFlag)[0])

    for idx in idxRand:
      pprBin = pprIdx[idx]
      rtrBin = rtrIdx[idx]
      pprP = pprCount[pprBin] / nLeft / pprBinWidth
      rtrP = rtrCount[rtrBin] / nLeft / rtrBinWidth
      
      #if(pprP > pprDensity[pprBin] or rtrP > rtrDensity[rtrBin] \
      #   or pprBin == 0 or rtrBin == 0):
      if(pprP + rtrP > pprDensity[pprBin] + rtrDensity[rtrBin] \
         or pprBin == 0 or rtrBin == 0 \
         or pprDensity[pprBin] < 1e-4 \
         or rtrDensity[rtrBin] < 1e-4):
        # This point is over represented, lets remove it
        keepFlag[idx] = False
        nLeft -= 1
        pprCount[pprIdx[idx]] -= 1
        rtrCount[rtrIdx[idx]] -= 1
      
    pars["u"]    = pars["u"][keepFlag]
    pars["ftau"] = pars["ftau"][keepFlag]
    pars["dtau"] = pars["dtau"][keepFlag]    
    pars["amp"]  = pars["amp"][keepFlag]    

    print("Reduced " + str(len(keepFlag)) + " down to " + str(sum(keepFlag)))
    
    return pars,amps[keepFlag,:]

  ############################################################################

  # This prunes the parameters to get back their distribution to match
  
  def pruneParameters(self,parType,pars):

    (dtauMean,dtauStd) = self.conInfo[dataType]["dtau"]
    (ftauMean,ftauStd) = self.conInfo[dataType]["ftau"]    
    (uMean,uStd) = self.conInfo[dataType]["U"]

    uPar = pars["u"]
    dtauPar = pars["dtau"]
    ftauPar = pars["ftau"]

     #!!! TODO
    
    
  ############################################################################

  def checkDataMatch(self,parType,pars,amps=None):
  
    pprMean,pprStd = self.conInfo[parType]["ppr"]
    rtrMean,rtrStd = self.conInfo[parType]["rtr"]    

    if(amps is None):
      amps = self.synapseModelWrapper(self.tStim,pars)

    ppr = np.divide(amps[:,1],amps[:,0])
    rtr = np.divide(amps[:,-1],amps[:,0])


    # self.plotSpikes(amps,self.tSpike,pars,0)
    
    
    if(True):
      self.plotDataMatch(values=ppr,
                         mean=pprMean,
                         std=pprStd,
                         title=parType,
                         xlabel="PPR")

      self.plotDataMatch(values=rtr,
                         mean=rtrMean,
                         std=rtrStd,
                         title=parType,
                         xlabel="RTR")

    return amps
      
  ############################################################################

  def plotDataMatch(self,values,mean,std,title,xlabel,nBins=20):
    
    print(str(xlabel) + " Mean : " + str(mean) + ", std : " + str(std))
    print("Data : " + str(np.mean(values)) + ", std : " + str(np.std(values)))
    nPoints = 100
    x = np.linspace(mean - 3*std, mean + 3*std, nPoints)

    plt.figure()
    plt.plot(x, stats.norm.pdf(x, mean, std))
    plt.hist(values,density=True,bins=nBins)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ion()
    plt.show()
    plt.tight_layout()    
    
    figName = "DATA/Planert2010/figs/Planert2010-surrogate-data-pre-fit-" \
      + str(title) + "-" + str(xlabel) + ".pdf"
    plt.savefig(figName)

  ############################################################################

  def plotSpikes(self,amps,tStim,pars=None,idx=0):

    
    plt.figure()
    for a,t in zip(amps[idx,:],tStim):
      plt.plot([t,t],[0,a],color="black")
    plt.xlabel("Time")
    plt.ylabel("Amplitude")  
    
    if(pars is not None):
      titleStr = "U = %.3f, ftau = %.3f, dtau = %.3f, amp = %.3f" \
        % (pars["u"][idx], pars["ftau"][idx], \
           pars["dtau"][idx], pars["amp"][idx])
      plt.title(titleStr)
      
    plt.ion()
    plt.show()
    
  ############################################################################

  # parType = "DD", "DI", "ID", "II", "FD" or "FI"
  
  def pickRandomParam(self,parType,num):

    print("Picking " + str(num) + " " + parType + " parameter sets...")
    
    # Need to randomise amp, dtau, ftau and U

    ampMean,ampStd = self.conInfo[parType]["amp"]
    dtauMean,dtauStd = self.conInfo[parType]["dtau"]
    ftauMean,ftauStd = self.conInfo[parType]["ftau"]
    uMean,uStd = self.conInfo[parType]["U"]        

    amp  = (np.random.normal(size=num,scale=ampStd)  + ampMean)*1e-3
    dtau = (np.random.normal(size=num,scale=dtauStd) + dtauMean)*1e-3
    ftau = (np.random.normal(size=num,scale=ftauStd) + ftauMean)*1e-3
    u    = np.random.normal(size=num,scale=uStd)    + uMean

    # If any of the conditions are true, we are outside OK range,
    # so only OK if all conditions are false, and the result sums to False
    okIdx = False == (np.array(amp < 0) \
                       + np.array(dtau < 0) \
                       + np.array(ftau < 0) \
                       + np.array(u < 0) \
                       + np.array(u > 1))

    
    # Then we need to verify that the ppr and rtr are within range

    print("Npoints = " + str(np.sum(okIdx)))

    pars = {"amp"  : amp[okIdx],
            "dtau" : dtau[okIdx],
            "ftau" : ftau[okIdx],
            "u"    : u[okIdx]}
    
    return pars
    
  ############################################################################

  def synapseModelWrapper(self,tStim,pars):

    print("Running " + str(len(pars["u"])) + " simulations")

    amps = np.zeros((len(pars["u"]),len(tStim)))

    for i,(u,dtau,ftau,Asc) in enumerate(zip(pars["u"],
                                             pars["dtau"],
                                             pars["ftau"],
                                             pars["amp"])):
      if(i > 0 and i % 100000 == 0):
        # Print progress...
        print(str(i) + "/" + str(amps.shape[0]))
        
      amps[i,:] = self.synapseModel(tStim,U=u,dtau=dtau,ftau=ftau,Asc=Asc)


    return amps
      
  ############################################################################
  
  def synapseModel(self,tStim,U,dtau,ftau,Asc):

    nSpikes = len(tStim)
    u = np.zeros((nSpikes,))
    r = np.zeros((nSpikes,))
    ISI = np.diff(tStim)

    # Init
    u[0] = U
    r[0] = 1
      
    for idx,ii in enumerate(ISI):
      
      u[idx+1] = u[idx]*np.exp(-ISI[idx]/ftau) \
                  + U*(1-u[idx]*np.exp(-ISI[idx]/ftau))
      
      r[idx+1] = r[idx]*(1-u[idx])*np.exp(-ISI[idx]/dtau) \
                 + 1 - np.exp(-ISI[idx]/dtau)

    amp = np.multiply(u,r)*Asc

    return amp
  
    
  ############################################################################
  
  def defExpData(self):
    self.conInfo = dict()

    self.conInfo["DD"] = { "P"    : (3,43),
                           "amp"  : (0.24, 0.15),
                           "ppr"  : (0.91, 0.63),
                           "rtr"  : (1.23, 0.50),
                           "dtau" : ( 192, 114 ),
                           "ftau" : (1266, 1427),
                           "U"    : (0.39, 0.22)}

    self.conInfo["DI"] = { "P" : (3,66),
                           "amp"  : (0.33,0.15),
                           "ppr"  : (0.84,0.3),
                           "rtr"  : (1.16,0.29),
                           "dtau" : (96,9),
                           "ftau" : (313,363),
                           "U"    : (0.46,0.24)}

    self.conInfo["ID"] = { "P" : (10,80),
                           "amp"  : (0.27,0.09),
                           "ppr"  : (1.1,0.6),
                           "rtr"  : (1.51,0.64),
                           "dtau" : (365,471),
                           "ftau" : (570,783),
                           "U"    : (0.36,0.18)}

    self.conInfo["II"] = { "P" : (7,31),
                           "amp"  : (0.45,0.44),
                           "ppr"  : (0.95,0.48),
                           "rtr"  : (1.39,0.69),
                           "dtau" : (149,90),
                           "ftau" : (1462,1800),
                           "U"    : (0.34,0.19)}

    self.conInfo["FD"] = { "P" : (8,9),
                           "amp"  : (4.8,4.9),
                           "ppr"  : (0.62,0.12),
                           "rtr"  : (0.72,0.08),
                           "dtau" : (740,350),
                           "ftau" : (3.1,2.4),
                           "U"    : (0.24,0.07)}

    self.conInfo["FI"] = { "P" : (6,9),
                          "amp"  : (3.1,4.1),
                           "ppr"  : (0.66,0.14),
                           "rtr"  : (0.63,0.19),
                           "dtau" : (850,500),
                           "ftau" : (4.5,2.7),
                           "U"    : (0.23,0.07)}

    self.tStim = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.9]

  ############################################################################

  def saveData(self,parType,pars,amp):
    fileName="DATA/Planert2010/PlanertFitting-"+parType+"-cache.json"

    data = dict()
    
    for p in pars:
      data[p] = pars[p]
      
    data["simAmp"] = amp
    data["tStim"] = self.tStim
    
    with open(fileName,"w") as f:
      json.dump(data,f,indent=2,cls=NumpyEncoder)

  ############################################################################

  def loadData(self,parType):
    fileName="Data/Planert2010/PlanertFitting-"+parType+"-cache.json"

    if(not os.path.exists(fileName)):
      return None,None
    
    with open(fileName,"r") as f:
      data = json.load(f)

    pars = dict([])
    for p in data:
      if(p in ["u","dtau","ftau","amp"]):
        pars[p] = data[p]
      elif(p != "simAmp"):
        print("Unknown parameter : " + str(p))

    amps = data["simAmp"]

    return pars, amps
        
  ############################################################################

  def plotModelParams(self,dataType,pars):

    (dtauMean,dtauStd) = self.conInfo[dataType]["dtau"]
    (ftauMean,ftauStd) = self.conInfo[dataType]["ftau"]    
    (uMean,uStd) = self.conInfo[dataType]["U"]

    self._plotModelParamsHelper(dataType,dtauMean*1e-3,dtauStd*1e-3,\
                                pars["dtau"],"dtau")
    self._plotModelParamsHelper(dataType,ftauMean*1e-3,ftauStd*1e-3,
                                pars["ftau"],"ftau")
    self._plotModelParamsHelper(dataType,uMean,uStd,pars["u"],"U")        
        
  ############################################################################
  
  def _plotModelParamsHelper(self,dataType,dataMean,dataStd,dataPoints,xlabel):

    xMin = np.min(dataPoints)
    xMax = np.max(dataPoints)

    x = np.linspace(xMin,xMax,100)
    xDensity = stats.norm.pdf(x,dataMean,dataStd) # / (x[1]-x[0])

    plt.figure()
    plt.hist(dataPoints,bins=20,density=True)
    plt.plot(x,xDensity)
    plt.xlabel(xlabel)
    plt.ylabel("Density")

    plt.tight_layout()

    
    figName = "DATA/Planert2010/figs/Surrogate-variables-distribution-" \
      + dataType + "-" + xlabel + ".pdf"
    plt.savefig(figName)
    
  ############################################################################

if __name__ == "__main__":

  nRuns = 1000000

  pp = Planert2010(parType="FI",num=nRuns)

  import pdb
  pdb.set_trace()
  
  pp = Planert2010(parType="FD",num=nRuns)
  
  pp = Planert2010(parType="DD",num=nRuns)
  pp = Planert2010(parType="DI",num=nRuns)
  pp = Planert2010(parType="ID",num=nRuns)
  pp = Planert2010(parType="II",num=nRuns)  

  
  import pdb
  pdb.set_trace()
