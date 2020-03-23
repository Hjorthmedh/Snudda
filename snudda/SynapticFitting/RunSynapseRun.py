import neuron
import numpy as np

from snudda.NrnSimulatorParallel import NrnSimulatorParallel

import bluepyopt.ephys as ephys
from snudda.Neuron_model_extended import NeuronModel
from snudda.Neuron_morphology import NeuronMorphology

# Plot all sections
# [neuron.h.psection(x) for x in neuron.h.allsec()]

##############################################################################

class RunSynapseRun(object):
  
  ############################################################################
  
  def __init__(self,
               neuronMorphology,
               neuronMechanisms,
               neuronParameters,
               neuronModulation,
               stimTimes,
               synapseDensity,
               nSynapses,
               neuronParameterID=0, # Which param set in parameter file to use
               neuronModulationID=0,
               holdingVoltage=-70e-3,
               synapseType='glut',
               params={},
               time=2.0):
    

    #plotting neuronMorphology
    '''
    import neurom as nm
    nrn=nm.load_neuron(neuronMorphology)
    from neurom import viewer
    fig, ax = viewer.draw(nrn) #2D plot
    fig.show() 
    fig, ax = viewer.draw(nrn, mode='3d')
    fig.show()
    '''
    #import pdb
    #pdb.set_trace()

    #plotting synapseDensity
    #import pdb
    #pdb.set_trace()

    print("Holding voltage: " + str(holdingVoltage) + " V")
    print("Stim times: " + str(stimTimes) + " s")
    print("Synapse type: " + str(synapseType))

    self.time = time
    self.synapses = []
    self.IClamp = None
    self.ncSyn = None

    # Done in NrnSimulatorParallel
    # neuron.h.load_file('stdrun.hoc')

    self.sim = NrnSimulatorParallel(cvode_active=False)
    
    # Should we use weak reference for garbage collection? (weakref package)

    # We load the neuron morphology object also, used to place synapses
    self.morphology = NeuronMorphology(swc_filename = neuronMorphology)

    # We need to setup the Neuron model
    self.neuron = NeuronModel(param_file=neuronParameters,
                              morph_file=neuronMorphology,
                              mech_file=neuronMechanisms,
                              cell_name="OptimisationNeuron",
                              modulation_file=neuronModulation,
                              parameterID=neuronParameterID,
                              modulationID=neuronModulationID)

    self.neuron.instantiate(sim=self.sim)
    self.setRestingVoltage(holdingVoltage)
    
    neuron.h.celsius = 35
      
    # gnabar_hh: The maximum specific sodium channel conductance [Default value = 0.120 S/cm2]
    # gkbar_hh: The maximum specific potassium channel conductance [Default value = 0.036 S/cm2]
    # gl_hh: The maximum specific leakage conductance [Default value = 0.0003 S/cm2]
    # ena: The reversal potential for the sodium channel [Default value = 50 mV]
    # ek: The reversal potential for the potassium channel [Default value = -77 mV]
    # el_hh: The reversal potential for the leakage channel [Default value = -54.3 mV]

    # We need to set the params also
    self.params = params
    self.defaultCond = 5e-7

    self.addSynapseDensity(synapseType,synapseDensity,nSynapses)
          
    self.stimTimes = stimTimes*1e3

    # Assumes input in seconds (SI units)
    self.connectInputToSynapses(stimTimes)
    
    self.somaRecord()
    self.synapseCurrentRecord()

    self.updateHoldingCurrent(holdingVoltage)

    #import pdb
    #pdb.set_trace()  
    
  ############################################################################

  def __del__(self):

    # This should not be needed but...
    self.neuron = None
    self.morphology = None
    self.IClamp = None
    self.VClamp = None
    self.vSave  = None
    self.tSave  = None
    self.iSave  = None
    self.ncSyn  = None
    
    self.vecStim    = None
    self.stimVector = None
    self.littleSynapse = None
    
  ############################################################################

  def updateHoldingCurrent(self,holdingVoltage=None):

    print("Updating holding current, might take a bit of time")
    
    if(holdingVoltage is None):
      holdingVoltage = self.holdingVoltage
    else:
      self.holdingVoltage = holdingVoltage

    # Disable old iClamp temporarilly
    if(self.IClamp is not None):
      self.IClamp.amp = 0
    
    # Setup a temporary VClamp
    self.VClamp = neuron.h.SEClamp(self.neuron.icell.soma[0](0.5))
    self.VClamp.rs = 1e-9
    self.VClamp.amp1 = holdingVoltage*1e3
    self.VClamp.dur1 = self.time*2*1e3
    # print("VClamp duration: " + str(self.VClamp.dur1))

    neuron.h.finitialize(self.holdingVoltage*1e3)
    # !!! There is a WEIRD neuron bug, that if this tstop here is
    # different from duration of simulation, then the *SECOND* time
    # a model is initialised we get the length of tSave set by this
    # value, and not by the tStop of that simulation --- go figure!
    self.setRestingVoltage(self.holdingVoltage*1e3)
    
    neuron.h.tstop = self.time*1e3 # Must set tstop
    neuron.h.run()

    if(False):
      import matplotlib.pyplot as plt
      plt.figure()
      plt.plot(self.tSave,self.vSave)
      plt.title("Holding voltage should be " \
                + str(self.holdingVoltage*1e3) + "mV")
      plt.xlabel("time (ms)")
      plt.ylabel("volt (mV)")
      plt.ion()
      plt.show()
          
      import pdb
      pdb.set_trace()    
    
    cur = float(self.VClamp.i)

    # Remove VClamp
    self.VClamp = None

    if(self.IClamp is None):
      self.IClamp = neuron.h.IClamp(self.neuron.icell.soma[0](0.5))
    
    # Update current on iClamp
    self.IClamp.amp = cur
    self.IClamp.dur = 2*self.time*1e3

    self.setRestingVoltage(self.holdingVoltage*1e3)

    print("Holding voltage " + str(self.holdingVoltage*1e3) + "mV, " \
          + "IClamp amp = " + str(cur) + "nA" )

    
  ############################################################################

  def setStimTimes(self,stimTimes):

    if(len(stimTimes) != len(self.stimTimes) \
       or (stimTimes != self.stimTimes).any()):

      print("Setting stim times to " + str(stimTimes) + " s")
      self.stimVector = neuron.h.Vector(stimTimes*1e3)
      self.stimTimes = stimTimes*1e3


  ############################################################################

  def addSynapseDensity(self,synapseType,
                        synapseDensity,
                        nSynapses=None):

    inputCoords,sectionID,sectionX = self.morphology.dendriteInputLocations( \
                                               synapseDensity=synapseDensity,
                                               nLocations=nSynapses)
    self.synapseLocations=inputCoords

    sections = self.neuron.mapIDtoCompartment(sectionID)
    
    for s,sx in zip(sections,sectionX):
      self.addSynapse(synapseType,s,sx,self.params)

  ############################################################################
    
  def addSynapse(self,synapseType,section,sectionX,params):

    sectionX = np.maximum(sectionX,1e-6) # Cant be 0 or 1
    
    try:
      if(synapseType.lower() == 'glut'):
        syn = neuron.h.tmGlut(section(sectionX))
      elif(synapseType.lower() == "gaba"):
        syn = neuron.h.tmGabaA(section(sectionX))
      else:
        assert False, "Unknown synapse type: " + str(synapseType)

      self.synapses.append(syn)
        
    except:
      import traceback
      tstr = traceback.format_exc()
      print(tstr)

      print("Did you remember to run nrnivmodl first, " \
            + "to generate channels mod files?")
      exit(-1)


    for p in params:

      # Conductance is set separately, it is a property of netcon
      if(p == "cond"):
        self.defaultCond = params["cond"]
        continue

      val = self.SItoNaturalUnits(p,params[p])
        
      setattr(syn,p,val)
      print("Setting parameters: " + str(p) + " = " + str(val) \
            + " (neuron natural units)")
      
      
  ############################################################################

  def connectInputToSynapses(self,stimTimes):

    print("Stimulation times (s): " + str(stimTimes))
    
    self.ncSyn = []
    
    self.stimVector = neuron.h.Vector(stimTimes*1e3)
    self.vecStim  = neuron.h.VecStim()
    self.vecStim.play(self.stimVector)

    for syn in self.synapses:
      ncs = neuron.h.NetCon(self.vecStim,syn)
      ncs.delay = 0
      ncs.threshold = 0
      self.ncSyn.append(ncs)

  ############################################################################

  def somaRecord(self):

    self.tSave = neuron.h.Vector()
    self.tSave.record(neuron.h._ref_t)

    self.vSave = neuron.h.Vector()
    self.vSave.record(self.neuron.icell.soma[0](0.5)._ref_v)

  ############################################################################

  def synapseCurrentRecord(self):

    self.iSave = []
    
    for syn in self.synapses:
      iRec = neuron.h.Vector()
      iRec.record(syn._ref_i)
      self.iSave.append(iRec)

  ############################################################################
    
  def run(self,tau,tauR,tauF,U,cond=None,time=None):

    if(time is None):
      time = self.time
    else:
      self.time = time

    if(cond is None):
      cond = self.defaultCond
      
    # print(vars())
    
    #print("Running: tau: %.3g, tauR: %.3g, tauF: %.3g, U: %.3g, cond: %.3g\n" \
    #      % (tau,tauR,tauF,U,cond))
    
    
    # Convert from SI units to natural units that Neuron uses
    for ncs in self.ncSyn:
      ncs.weight[0] = 1*cond*1e6

    for syn in self.synapses:
      syn.tau = tau*1e3
      syn.tauR = tauR*1e3
      syn.tauF = tauF*1e3
      syn.U = U
    
    # print(self.littleSynapse.tau)


    # print("Initialise voltage to " + str(self.holdingVoltage*1e3) + " mV")
    neuron.h.finitialize(self.holdingVoltage*1e3) # OLD : -70
    neuron.h.tstop = time*1e3

    
    neuron.h.run()

    #self.tSave.resize()
    #self.vSave.resize()
    #self.iSave.resize()
    
    if(np.array(self.tSave).shape != np.array(self.vSave).shape):
      print("THIS IS WRONG, should be same shape!!")
      print("size t = " + str(np.array(self.tSave).shape))
      print("size v = " + str(np.array(self.vSave).shape))      
      import pdb
      pdb.set_trace()
    
    # print("Last V = " + str(self.vSave[-1]*1e-3))
    
    # Convert back to SI units
    return (np.array(self.tSave)*1e-3,
            np.array(self.vSave)*1e-3,
            np.array(self.iSave)*1e-9)

  ############################################################################

  def setRestingVoltage(self,restVolt):

    soma = [x for x in self.neuron.icell.soma]
    axon = [x for x in self.neuron.icell.axon]
    dend = [x for x in self.neuron.icell.dend]

    cell = soma+axon+dend

    for sec in cell:
      for seg in sec.allseg():
        seg.v = restVolt
  
  ############################################################################

  # I wish Neuron would use natural units...
  
  def SItoNaturalUnits(self,varName,value):

    convFactor = { "U" : 1.0,
                   "tauR" : 1e3,
                   "tauF" : 1e3,
                   "cond" : 1e6,
                   "tau" : 1e3,
                   "nmda_ratio" : 1.0 }

    if varName not in convFactor:
      print("Missing conversion fractor for " + str(varName) \
            + ". Please update SItoNaturalUnits function.")
      print("convFactor = " + str(convFactor))
      import pdb
      pdb.set_trace()

    try:
      return value*convFactor[varName]
    except:
      import traceback
      tstr = traceback.format_exc()
      print(tstr)
      import pdb
      pdb.set_trace()
    
    
    
  ############################################################################
  def plot(self):
    
    ax=self.morphology.plotNeuron()

    ax.scatter(self.synapseLocations[:,0],
               self.synapseLocations[:,1],
               self.synapseLocations[:,2],
               color='black', 
               marker='o',
               s=20)
    plt.ion()    


  ############################################################################
  
  # pars = { "tau1" : 0.25e-3 }
  # The code converts to natural units internally, if your parameter is missing
  # then update SItoNaturalUnits

  # OBS, soma parameters are ignored by run2 (they should be set during setup)
  
  def run2(self,pars,time=None,cond=1e-8):

    if(time is None):
      time = self.time
    else:
      self.time = time

    for p in pars:
      if(p in ["somaDiameter","somaCM","somaGleak"]):
        # Ignore soma parameters at this stage, should be set at setup
        continue
      
      if(p == "cond"):
        cond = self.SItoNaturalUnits(p,pars[p])
      else:
        v = self.SItoNaturalUnits(p,pars[p])
        setattr(self.littleSynapse,p,v)
        
    neuron.h.finitialize(self.holdingVoltage*1e3)  
    self.ncSyn.weight[0] = cond
    #print("Synapse conductance: " + str(cond) + " uS")
    #print("Verification of conductance: " + str(self.ncSyn.weight[0]))

    neuron.h.tstop = time*1e3
    neuron.h.run()

    # Convert results back to SI units
    return (np.array(self.tSave)*1e-3,
            np.array(self.vSave)*1e-3,
            np.array(self.iSave)*1e-9)
  
##############################################################################
  
if __name__== "__main__":
  
  import matplotlib.pyplot as plt
  import time
  start = time.time() 
  
  stimTimes = np.array([0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,1.0])
  tau  = 10e-3
  tauR = 200e-3
  tauF = 900e-3
  U    = 0.5

  tau = 0.1
  tauR = 0.1
  tauF = 0.1
  U = 0.1
  cond = 1e-9
  
  
  rlsr = RunSynapseRun(stimTimes=stimTimes,
                       synapseType="glut",
                       neuronMorphology="../data/cellspecs/dspn/str-dspn-e150917_c9_d1-mWT-1215MSN03-v20190521/WT-1215MSN03-cor-rep-ax2.swc",
                       neuronParameters="../data/cellspecs/dspn/str-dspn-e150917_c9_d1-mWT-1215MSN03-v20190521/parameters.json",
                       neuronMechanisms="../data/cellspecs/dspn/str-dspn-e150917_c9_d1-mWT-1215MSN03-v20190521/mechanisms.json",
                       neuronModulation="../data/cellspecs/dspn/str-dspn-e150917_c9_d1-mWT-1215MSN03-v20190521/modulation.json",
                       synapseDensity="0.05/(1+np.exp(-(d-30e-6)/5e-6))",
                       nSynapses=20,
                       neuronParameterID=0,
                       neuronModulationID=0)

  # I would like to plot the morphology
  ax = rlsr.plot()

  for i in range(3):
    tS,vS,iS = rlsr.run(tau*i,tauR,tauF,U,cond)
    plt.figure()
    plt.plot(tS,vS)
    plt.xlabel("Time (s)")
    plt.ylabel("Voltage (V)")
    plt.figure()
    
    for iiS in iS:
      plt.plot(tS,iiS)
      
    plt.xlabel("Time (s)")
    plt.ylabel("Current (A)")
  
  end = time.time()
  print(end - start)

  plt.ion()
  plt.show()

  import pdb
  pdb.set_trace()
