import neuron
import numpy as np

import bluepyopt.ephys as ephys
from snudda.simulate.nrn_simulator_parallel import NrnSimulatorParallel

from snudda.neurons.neuron_model_extended import NeuronModel
from snudda.neurons.neuron_morphology import NeuronMorphology

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
               synapseSectionID=None, # if given, nSynapses is ignored
               synapseSectionX=None,  # same # of elements as synapseSectionID
               neuronParameterID=0, # Which param set in parameter file to use
               neuronModulationID=0,
               holdingVoltage=-70e-3,
               synapseType='glut',
               params={},
               time=2,
               logFile=None,
               verbose=True,
               rng=None):

    self.logFile = logFile # File pointer
    self.verbose = verbose

    self.writeLog("Holding voltage: " + str(holdingVoltage) + " V")
    self.writeLog("Stim times: " + str(stimTimes) + " s")
    self.writeLog("Synapse type: " + str(synapseType))

    self.time = time
    self.synapses = []
    self.IClamp = None
    self.ncSyn = None

    # Done in NrnSimulatorParallel
    # neuron.h.load_file('stdrun.hoc')

    self.sim = NrnSimulatorParallel(cvode_active=False)
    
    # Should we use weak reference for garbage collection? (weakref package)

    # We load the neuron morphology object also, used to place synapses
    self.writeLog("Using morphology: " + str(neuronMorphology))    
    self.morphology = NeuronMorphology(swc_filename = neuronMorphology)

    # We need to setup the Neuron model
    self.neuron = NeuronModel(param_file=neuronParameters,
                              morph_path=neuronMorphology,
                              mech_file=neuronMechanisms,
                              cell_name="OptimisationNeuron",
                              modulation_file=neuronModulation,
                              parameter_id=neuronParameterID,
                              modulation_id=neuronModulationID)


    self.neuron.instantiate(sim=self.sim)
    self.setRestingVoltage(holdingVoltage*1e3)
    
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
    self.rng = rng

    assert rng, " Code has been updated, rng must be a random number generator, eg. rng = np.random.default_rng(random_seed)"

    self.addSynapseDensity(synapseType,synapseDensity,nSynapses,
                           synapseSectionID,synapseSectionX,
                           rng=self.rng)
          
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

    self.writeLog("Updating holding current, might take a bit of time")
    
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
    # self.writeLog("VClamp duration: " + str(self.VClamp.dur1))

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
      plt.plot(self.t_save, self.v_save)
      plt.title("Holding voltage should be " \
                + str(self.holding_voltage * 1e3) + "mV")
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
      self.IClamp = neuron.h.i_clamp(self.neuron.icell.soma[0](0.5))
    
    # Update current on iClamp
    self.IClamp.amp = cur
    self.IClamp.dur = 2*self.time*1e3

    self.setRestingVoltage(self.holdingVoltage*1e3)

    self.writeLog("Holding voltage " + str(self.holdingVoltage*1e3) + "mV, " \
          + "IClamp amp = " + str(cur) + "nA" )

    
  ############################################################################

  def setStimTimes(self,stimTimes):

    if(len(stimTimes) != len(self.stimTimes) \
       or (stimTimes != self.stimTimes).any()):


      print("Setting stim times to " + str(stimTimes) + " s")
      self.writeLog("Setting stim times to " + str(stimTimes) + " s")
      self.stimVector = neuron.h.Vector(stimTimes*1e3)
      self.stimTimes = stimTimes*1e3


  ############################################################################

  def addSynapseDensity(self,synapseType,
                        synapseDensity,
                        rng,
                        nSynapses=None,
                        plotFlag=False,
                        sectionID=None,
                        sectionX=None):


    if(plotFlag):

      assert sectionID is None and sectionX is None, \
        "Can not plot if sectionID and sectionX are given"
      
    
      inputCoords,sectionID,sectionX,densityFunction,distSynSoma =\
                                    self.morphology.dendrite_input_locations( \
                                                synapse_density=synapseDensity,
                                                num_locations=nSynapses,
                                                return_density=True,
                                                rng=rng)

      self.synapseLocations=inputCoords
      distFromSoma = self.morphology.dend[:,4]

      #plot density function
      plt.figure()
      plt.plot(distFromSoma*1e6,densityFunction,'o')
      plt.xlabel('distance from soma $(\mu m)$')
      plt.title('density function')
      plt.ion()
      plt.show()

      #plot histogram - distance synapses from soma
      plt.figure()
      plt.hist(distSynSoma*1e6,edgecolor='gray', bins=distSynSoma.size)
      plt.xlabel('distance from soma $(\mu m)$')
      plt.ylabel('frequency')
      plt.title('synaptic distance from soma')
      plt.ion()
      plt.show()

    elif(sectionID is None or sectionX is None):

      inputCoords,sectionID,sectionX =\
        self.morphology.dendrite_input_locations(synapse_density=synapseDensity,
                                                 num_locations=nSynapses,
                                                 rng=rng)

      
      self.synapseLocations=inputCoords

    else:
      self.synapseLocations = "Unknown, you specified sectionX and sectionID, "\
                            + "but not synapse xyz coordinates."
      

    # This is so we can find out where the synapses were placed
    self.synapseSectionID = sectionID
    self.synapseSectionX = sectionX
      
    sections = self.neuron.map_id_to_compartment(sectionID)
    
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
      self.writeLog(tstr)

      self.writeLog("Did you remember to run nrnivmodl first, " \
            + "to generate channels mod files?")
      exit(-1)


    for p in params:

      # Conductance is set separately, it is a property of netcon
      if(p == "cond"):
        self.defaultCond = params["cond"]
        continue

      val = self.SItoNaturalUnits(p,params[p])
        
      setattr(syn,p,val)
      self.writeLog("Setting parameters: " + str(p) + " = " + str(val) \
            + " (neuron natural units)")
      
      
  ############################################################################

  def connectInputToSynapses(self,stimTimes):

    self.writeLog("Stimulation times (s): " + str(stimTimes))
    
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

    assert False, "This is the old run method"
    
    if(time is None):
      time = self.time
    else:
      self.time = time

    if(cond is None):
      cond = self.default_cond
      
    # print(vars())
    
    #print("Running: tau: %.3g, tauR: %.3g, tauF: %.3g, U: %.3g, cond: %.3g\n" \
    #      % (tau,tauR,tauF,U,cond))
    
    
    # Convert from SI units to natural units that Neuron uses
    for ncs in self.nc_syn:
      ncs.weight[0] = 1*cond*1e6

    for syn in self.synapses:
      syn.tau = tau*1e3
      syn.tau_r = tauR * 1e3
      syn.tau_f = tauF * 1e3
      syn.u = U
    
    # print(self.littleSynapse.tau)


    # print("Initialise voltage to " + str(self.holdingVoltage*1e3) + " mV")
    neuron.h.finitialize(self.holding_voltage * 1e3) # OLD : -70
    neuron.h.tstop = time*1e3

    
    neuron.h.run()

    #self.tSave.resize()
    #self.vSave.resize()
    #self.iSave.resize()
    
    if(np.array(self.t_save).shape != np.array(self.v_save).shape):
      self.write_log("THIS IS WRONG, should be same shape!!")
      self.write_log("size t = " + str(np.array(self.t_save).shape))
      self.write_log("size v = " + str(np.array(self.v_save).shape))
      import pdb
      pdb.set_trace()
    
    # print("Last V = " + str(self.vSave[-1]*1e-3))
    
    # Convert back to SI units
    return (np.array(self.t_save) * 1e-3,
            np.array(self.v_save) * 1e-3,
            np.array(self.i_save) * 1e-9)

  ############################################################################

  def setRestingVoltage(self,restVolt):

    self.writeLog("Setting resting voltage to %.3f mV" % restVolt)
    
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
      self.writeLog("Missing conversion fractor for " + str(varName) \
            + ". Please update SItoNaturalUnits function.")
      self.writeLog("convFactor = " + str(convFactor))
      import pdb
      pdb.set_trace()

    try:
      return value*convFactor[varName]
    except:
      import traceback
      tstr = traceback.format_exc()
      self.writeLog(tstr)
      import pdb
      pdb.set_trace()
    
    
    
  ############################################################################
  def plot(self):

    ax=self.morphology.plot_neuron(axon_colour='red', dend_colour='blue', soma_colour='green')

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

    self.writeLog("Running with pars: " + str(pars))
    
    if(time is None):
      time = self.time
    else:
      self.time = time

    for p in pars:
      
      if(p == "cond"):
        cond = self.SItoNaturalUnits(p,pars[p])
          
      else:
        v = self.SItoNaturalUnits(p,pars[p])
        for s in self.synapses:
          setattr(s,p,v)
        
    neuron.h.finitialize(self.holdingVoltage*1e3)  
    for ncs in self.ncSyn:
      ncs.weight[0] = cond

    # self.writeLog("Synapse conductance: " + str(cond) + " uS")

    self.setRestingVoltage(self.holdingVoltage*1e3)

#    print("Check voltage")
#    import pdb
#    pdb.set_trace()

    neuron.h.v_hold = self.holdingVoltage * 1e3
    neuron.h.tstop = time*1e3
    neuron.h.run()

    # Convert results back to SI units
    return (np.array(self.tSave)*1e-3,
            np.array(self.vSave)*1e-3,
            np.array(self.iSave)*1e-9)

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

  
##############################################################################
  
if __name__== "__main__":
  
  import matplotlib.pyplot as plt
  import time
  start = time.time() 
  
  stimTimes = np.array([0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,1.0]) + 0.3

  #synapseDensity="0.05/(1+np.exp(-(d-30e-6)/5e-6))"  
  # synapseDensity="0.05/(1+np.exp(-(d-150e-6)/5e-6))"

  synapseDensity="np.exp(-((d-50e-6)/15e-6)**2)"  
  holdingVoltage = -79.5e-3
  
  # !!! OBS changed number of synapses, before was 20 -- cond 1nS
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
                       synapseDensity=synapseDensity,
                       nSynapses=10,
                       neuronParameterID=0,
                       neuronModulationID=0,
                       holdingVoltage=holdingVoltage)

  # I would like to plot the morphology
  # ax = rlsr.plot()

  #import pdb
  #pdb.set_trace()

  params = { "U" : 0.5,
             "tauR" : 0.6,
             "tauF" : 2.0, #1.1
             "tau" : 0.08,
             "cond" : 1e-9,
             "nmda_ratio" : 0.5 }
  

  tS,vS,iS = rlsr.run2(pars=params)
  plt.figure()
  plt.plot(tS,vS)
  plt.xlabel("Time (s)")
  plt.ylabel("Voltage (V)")

  if(True):
    expData = 'DATA/YvonneJohansson2019/M1LH_Analysis_191001.h5'
    expCellID = 131

    import LoadExpData

    led = LoadExpData.LoadExpData(expData)
    cellType = led.get_cell_type(expCellID)
    plt.title("ExpData is cellType : " + str(cellType))
    expV,expT = led.getData("GBZ_CC_H20",expCellID)
    plt.plot(expT,expV,'r-')
  
  if(False):
    plt.figure()
  
    for iiS in iS:
      plt.plot(tS,iiS)
      
    plt.xlabel("Time (s)")
    plt.ylabel("Current (A)")

  plt.ion()
  plt.show()

    
  end = time.time()
  print(end - start)


  import pdb
  pdb.set_trace()
