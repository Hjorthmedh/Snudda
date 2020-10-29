import neuron
import numpy as np

# Plot all sections
# [neuron.h.psection(x) for x in neuron.h.allsec()]

##############################################################################

class RunLittleSynapseRun(object):
  
  ############################################################################
  
  def __init__(self,
               stimTimes,
               holdingVoltage=-70e-3,
               synapseType='glut',
               params={},
               time=2.0):

    # This just clears the sections, just to be on the safe side
    # !!! Didnt work... still problem
    #for s in neuron.h.allsec():
    #  import pdb
    #  pdb.set_trace()
    #  neuron.h.delete_section(sec=s)
    
    print("Holding voltage: " + str(holdingVoltage) + " V")
    print("Stim times: " + str(stimTimes) + " s")
    print("Synapse type: " + str(synapseType))

    self.time = time
    
    neuron.h.load_file('stdrun.hoc')

    # Should we use weak reference for garbage collection? (weakref package)

    self.soma = neuron.h.Section(name='soma')
    self.soma.insert('hh')
    neuron.h.celsius = 35
      
    # gnabar_hh: The maximum specific sodium channel conductance [Default value = 0.120 S/cm2]
    # gkbar_hh: The maximum specific potassium channel conductance [Default value = 0.036 S/cm2]
    # gl_hh: The maximum specific leakage conductance [Default value = 0.0003 S/cm2]
    # ena: The reversal potential for the sodium channel [Default value = 50 mV]
    # ek: The reversal potential for the potassium channel [Default value = -77 mV]
    # el_hh: The reversal potential for the leakage channel [Default value = -54.3 mV]

    try:
      if(synapseType.lower() == 'glut'):
        self.littleSynapse=neuron.h.tmGlut(self.soma(0.5))
      elif(synapseType.lower() == "gaba"):
        self.littleSynapse=neuron.h.tmGabaA(self.soma(0.5))
      else:
        assert False, "Unknown synapse type: " + str(synapseType)

    except:
      import traceback
      tstr = traceback.format_exc()
      print(tstr)

      print("Did you remember to run nrnivmodl first, to generate channels mod files?")
      exit(-1)

      
    # We need to set the params also
    self.params = params

    self.defaultCond = 5e-7
    
    for p in self.params:

      # Conductance is set separately, it is a property of netcon
      if(p == "cond"):
        self.defaultCond = self.params["cond"]
        continue

      if(p == "somaDiameter"):
        print("Setting soma diameter to " \
              + str(self.params["somaDiameter"]*1e6) + " mum")
        self.soma.L = self.params["somaDiameter"]*1e6
        self.soma.diam = self.params["somaDiameter"]*1e6
        continue
        
      if(p == "somaCM"):
        print("Setting soma CM to " \
              + str(self.params["somaCM"]*1e2) + " uf/cm2")
        self.soma.cm = self.params["somaCM"]*1e2
        continue

      if(p == "somaGleak"):
        print("Setting GL_HH in soma to " \
              + str(self.params["somaGleak"]*1e-4) + " S/cm2")
        self.soma.gl_hh = self.params["somaGleak"]*1e4
        continue
        
      val = self.SItoNaturalUnits(p,self.params[p])
        
      setattr(self.littleSynapse,p,val)
      print("Setting parameters: " + str(p) + " = " + str(val) \
            + " (neuron natural units)")
      
    self.vecStim  = neuron.h.VecStim()
    
    self.stimTimes = stimTimes*1e3
      
    self.stimVector = neuron.h.Vector(stimTimes*1e3)
    self.vecStim.play(self.stimVector)
    
    self.ncSyn = neuron.h.NetCon(self.vecStim,self.littleSynapse)
    self.ncSyn.delay = 0
    self.ncSyn.threshold = 0
    
    self.tSave = neuron.h.Vector()
    self.tSave.record(neuron.h._ref_t)

    self.vSave = neuron.h.Vector()
    self.vSave.record(self.soma(0.5)._ref_v)

    
    # Channel current save
    self.iSave = neuron.h.Vector()
    self.iSave.record(self.littleSynapse._ref_i)

    self.holdingVoltage = holdingVoltage
    
    self.VClamp = neuron.h.SEClamp(self.soma(0.5))

    self.vClampRS = 1e-9 # 1e-6
    self.VClamp.rs = self.vClampRS
    self.VClamp.amp1 = holdingVoltage*1e3
    self.VClamp.dur1 = self.time*2*1e3
    # print("VClamp duration: " + str(self.VClamp.dur1))

    neuron.h.finitialize(self.holdingVoltage*1e3)
    # !!! There is a WEIRD neuron bug, that if this tstop here is
    # different from duration of simulation, then the *SECOND* time
    # a model is initialised we get the length of tSave set by this
    # value, and not by the tStop of that simulation --- go figure!
    neuron.h.tstop = self.time*1e3 # Must set tstop
    neuron.h.run()
    
    #neuron.h.tstop = 2001
    #neuron.h.finitialize(self.holdingVoltage*1e3)
    #neuron.h.run()
    
    if(False):
      import matplotlib.pyplot as plt
      plt.figure()
      plt.plot(self.t_save, self.v_save)
      plt.title("Holding voltage should be " \
                + str(self.holdingVoltage*1e3) + "mV")
      plt.xlabel("time (ms)")
      plt.ylabel("volt (mV)")
      plt.ion()
      plt.show()
      import pdb
      pdb.set_trace()
    
      
      import pdb
      pdb.set_trace()
    
    cur = float(self.VClamp.i)

    # Remove VClamp
    self.VClamp = None
    
    self.IClamp = neuron.h.IClamp(self.soma(0.5))
    self.IClamp.amp = cur #nA
    self.IClamp.dur = 2*self.time*1e3

    print("Holding voltage " + str(holdingVoltage*1e3) + "mV, " \
          + "IClamp amp = " + str(cur) + "nA" )

  ############################################################################

  def __del__(self):

    # This should not be needed but...
    self.soma   = None
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

    if(holdingVoltage is None):
      holdingVoltage = self.holdingVoltage
    else:
      self.holdingVoltage = holdingVoltage

    # Disable old iClamp temporarilly
    self.IClamp.amp = 0

    # Setup a temporary VClamp
    self.VClamp = neuron.h.SEClamp(self.soma(0.5))
    self.VClamp.rs = self.vClampRS
    self.VClamp.amp1 = holdingVoltage*1e3
    self.VClamp.dur1 = self.time*2*1e3
    # print("VClamp duration: " + str(self.VClamp.dur1))

    neuron.h.finitialize(self.holdingVoltage*1e3)
    # !!! There is a WEIRD neuron bug, that if this tstop here is
    # different from duration of simulation, then the *SECOND* time
    # a model is initialised we get the length of tSave set by this
    # value, and not by the tStop of that simulation --- go figure!
    neuron.h.tstop = self.time*1e3 # Must set tstop
    neuron.h.run()
   
    cur = float(self.VClamp.i)

    # Remove VClamp
    self.VClamp = None

    # Update current on iClamp
    self.IClamp.amp = cur
    self.IClamp.dur = 2*self.time*1e3

    
  ############################################################################

  def setStimTimes(self,stimTimes):

    if(len(stimTimes) != len(self.stimTimes) \
       or (stimTimes != self.stimTimes).any()):

      print("Setting stim times to " + str(stimTimes) + " s")
      self.stimVector = neuron.h.Vector(stimTimes*1e3)
      self.stimTimes = stimTimes*1e3

    
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
    self.ncSyn.weight[0] = 1*cond*1e6    
    self.littleSynapse.tau = tau*1e3
    self.littleSynapse.tauR = tauR*1e3
    self.littleSynapse.tauF = tauF*1e3
    self.littleSynapse.U = U
    
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
  cond = 1e-7
  
  
  rlsr = RunLittleSynapseRun(stimTimes,synapseType="glut")
  
  for i in range(10):
    tS,vS,iS = rlsr.run(tau*i,tauR,tauF,U,cond)
    plt.figure()
    plt.plot(tS,vS)
    plt.xlabel("Time (s)")
    plt.ylabel("Voltage (V)")
    plt.figure()
    plt.plot(tS,iS)
    plt.xlabel("Time (s)")
    plt.ylabel("Current (A)")

  #tS,vS = rlsr.run(tau*1,tauR,tauF,U,cond)
  #plt.plot(tS,vS)
    
    
  plt.xlabel("Time")
  plt.ylabel("Voltage")
  plt.title("How tau affects dynamics")
  
  end = time.time()
  print(end - start)

  plt.ion()
  plt.show()

  import pdb
  pdb.set_trace()
