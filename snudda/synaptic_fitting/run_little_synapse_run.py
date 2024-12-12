import neuron
import numpy as np
import sys


# Plot all sections
# [neuron.h.psection(x) for x in neuron.h.allsec()]

##############################################################################

class RunLittleSynapseRun(object):

    ############################################################################

    def __init__(self,
                 stim_times,
                 holding_voltage=-70e-3,
                 synapse_type='glut',
                 params={},
                 time=2.0):

        # This just clears the sections, just to be on the safe side
        # !!! Didnt work... still problem
        # for s in neuron.h.allsec():
        #  import pdb
        #  pdb.set_trace()
        #  neuron.h.delete_section(sec=s)

        print(f"Holding voltage: {holding_voltage} V")
        print(f"Stim times: {stim_times} s")
        print(f"Synapse type: {synapse_type}")

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
            if synapse_type.lower() == 'glut':
                self.little_synapse = neuron.h.tmGlut(self.soma(0.5))
            elif synapse_type.lower() == "gaba":
                self.little_synapse = neuron.h.tmGabaA(self.soma(0.5))
            else:
                assert False, f"Unknown synapse type: {synapse_type}"

        except:
            import traceback
            tstr = traceback.format_exc()
            print(tstr)

            print("Did you remember to run nrnivmodl first, to generate channels mod files?")
            sys.exit(-1)

        # We need to set the params also
        self.params = params

        self.default_cond = 5e-7

        for p in self.params:

            # Conductance is set separately, it is a property of netcon
            if p == "cond":
                self.default_cond = self.params["cond"]
                continue

            if p == "somaDiameter":
                print(f"Setting soma diameter to {self.params['somaDiameter'] * 1e6} mum")
                self.soma.L = self.params["somaDiameter"] * 1e6
                self.soma.diam = self.params["somaDiameter"] * 1e6
                continue

            if p == "somaCM":
                print(f"Setting soma CM to {self.params['somaCM'] * 1e2} uf/cm2")
                self.soma.cm = self.params["somaCM"] * 1e2
                continue

            if p == "somaGleak":
                print("Setting GL_HH in soma to {self.params['somaGleak'] * 1e-4} S/cm2")
                self.soma.gl_hh = self.params["somaGleak"] * 1e4
                continue

            val = self.si_to_natural_units(p, self.params[p])

            setattr(self.little_synapse, p, val)
            print(f"Setting parameters: {p} = {val} (neuron natural units)")

        self.vec_stim = neuron.h.VecStim()

        self.stim_times = stim_times * 1e3

        self.stim_vector = neuron.h.Vector(stim_times * 1e3)
        self.vec_stim.play(self.stim_vector)

        self.nc_syn = neuron.h.NetCon(self.vec_stim, self.little_synapse)
        self.nc_syn.delay = 0
        self.nc_syn.threshold = 0

        self.t_save = neuron.h.Vector()
        self.t_save.record(neuron.h._ref_t)

        self.v_save = neuron.h.Vector()
        self.v_save.record(self.soma(0.5)._ref_v)

        # Channel current save
        self.i_save = neuron.h.Vector()
        self.i_save.record(self.little_synapse._ref_i)

        self.holding_voltage = holding_voltage

        self.v_clamp = neuron.h.SEClamp(self.soma(0.5))

        self.v_clamp_rs = 1e-9  # 1e-6
        self.v_clamp.rs = self.v_clamp_rs
        self.v_clamp.amp1 = holding_voltage * 1e3
        self.v_clamp.dur1 = self.time * 2 * 1e3
        # print("VClamp duration: " + str(self.VClamp.dur1))

        neuron.h.finitialize(self.holding_voltage * 1e3)
        # !!! There is a WEIRD neuron bug, that if this tstop here is
        # different from duration of simulation, then the *SECOND* time
        # a model is initialised we get the length of tSave set by this
        # value, and not by the tStop of that simulation --- go figure!
        neuron.h.tstop = self.time * 1e3  # Must set tstop
        neuron.h.run()

        # neuron.h.tstop = 2001
        # neuron.h.finitialize(self.holdingVoltage*1e3)
        # neuron.h.run()

        if False:
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

            import pdb
            pdb.set_trace()

        cur = float(self.v_clamp.i)

        # Remove VClamp
        self.v_clamp = None

        self.i_clamp = neuron.h.IClamp(self.soma(0.5))
        self.i_clamp.amp = cur  # nA
        self.i_clamp.dur = 2 * self.time * 1e3

        print(f"Holding voltage {holding_voltage * 1e3} mV, IClamp amp = {cur} nA")

    ############################################################################

    def __del__(self):

        # This should not be needed but...
        self.soma = None
        self.i_clamp = None
        self.v_clamp = None
        self.v_save = None
        self.t_save = None
        self.i_save = None
        self.nc_syn = None

        self.vec_stim = None
        self.stim_vector = None
        self.little_synapse = None

    ############################################################################

    def update_holding_current(self, holding_voltage=None):

        if holding_voltage is None:
            holding_voltage = self.holding_voltage
        else:
            self.holding_voltage = holding_voltage

        # Disable old iClamp temporarily
        self.i_clamp.amp = 0

        # Setup a temporary VClamp
        self.v_clamp = neuron.h.SEClamp(self.soma(0.5))
        self.v_clamp.rs = self.v_clamp_rs
        self.v_clamp.amp1 = holding_voltage * 1e3
        self.v_clamp.dur1 = self.time * 2 * 1e3
        # print("VClamp duration: " + str(self.VClamp.dur1))

        neuron.h.finitialize(self.holding_voltage * 1e3)
        # !!! There is a WEIRD neuron bug, that if this tstop here is
        # different from duration of simulation, then the *SECOND* time
        # a model is initialised we get the length of tSave set by this
        # value, and not by the tStop of that simulation --- go figure!
        neuron.h.tstop = self.time * 1e3  # Must set tstop
        neuron.h.run()

        cur = float(self.v_clamp.i)

        # Remove VClamp
        self.v_clamp = None

        # Update current on iClamp
        self.i_clamp.amp = cur
        self.i_clamp.dur = 2 * self.time * 1e3

    ############################################################################

    def set_stim_times(self, stim_times):

        if len(stim_times) != len(self.stim_times) or (stim_times != self.stim_times).any():
            print(f"Setting stim times to {stim_times} s")
            self.stim_vector = neuron.h.Vector(stim_times * 1e3)
            self.stim_times = stim_times * 1e3

    ############################################################################

    def run(self, tau, tau_r, tau_f, u, cond=None, time=None):

        if time is None:
            time = self.time
        else:
            self.time = time

        if cond is None:
            cond = self.default_cond

        # print(vars())

        # print("Running: tau: %.3g, tauR: %.3g, tauF: %.3g, U: %.3g, cond: %.3g\n" \
        #      % (tau,tauR,tauF,U,cond))

        # Convert from SI units to natural units that Neuron uses
        self.nc_syn.weight[0] = 1 * cond * 1e6
        self.little_synapse.tau = tau * 1e3
        self.little_synapse.tauR = tau_r * 1e3
        self.little_synapse.tauF = tau_f * 1e3
        self.little_synapse.U = u

        # print(self.littleSynapse.tau)

        # print("Initialise voltage to " + str(self.holdingVoltage*1e3) + " mV")
        neuron.h.finitialize(self.holding_voltage * 1e3)  # OLD : -70
        neuron.h.tstop = time * 1e3

        neuron.h.run()

        # self.tSave.resize()
        # self.vSave.resize()
        # self.iSave.resize()

        if np.array(self.t_save).shape != np.array(self.v_save).shape:
            print("THIS IS WRONG, should be same shape!!")
            print(f"size t = {np.array(self.t_save).shape}")
            print(f"size v = {np.array(self.v_save).shape}")
            import pdb
            pdb.set_trace()

        # print("Last V = " + str(self.vSave[-1]*1e-3))

        # Convert back to SI units
        return (np.array(self.t_save) * 1e-3,
                np.array(self.v_save) * 1e-3,
                np.array(self.i_save) * 1e-9)

    ############################################################################

    # I wish Neuron would use natural units...

    def si_to_natural_units(self, var_name, value):

        conv_factor = {"U": 1.0,
                       "tauR": 1e3,
                       "tauF": 1e3,
                       "cond": 1e6,
                       "tau": 1e3,
                       "nmda_ratio": 1.0}

        if var_name not in conv_factor:
            print(f"Missing conversion fractor for {var_name}. Please update SItoNaturalUnits function.")
            print("convFactor = " + str(conv_factor))
            import pdb
            pdb.set_trace()

        try:
            return value * conv_factor[var_name]
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

    def run2(self, pars, time=None, cond=1e-8):

        if time is None:
            time = self.time
        else:
            self.time = time

        for p in pars:
            if p in ["somaDiameter", "somaCM", "somaGleak"]:
                # Ignore soma parameters at this stage, should be set at setup
                continue

            if p == "cond":
                cond = self.si_to_natural_units(p, pars[p])
            else:
                v = self.si_to_natural_units(p, pars[p])
                setattr(self.little_synapse, p, v)

        neuron.h.finitialize(self.holding_voltage * 1e3)
        self.nc_syn.weight[0] = cond
        # print("Synapse conductance: " + str(cond) + " uS")
        # print("Verification of conductance: " + str(self.ncSyn.weight[0]))

        neuron.h.tstop = time * 1e3
        neuron.h.run()

        # Convert results back to SI units
        return (np.array(self.t_save) * 1e-3,
                np.array(self.v_save) * 1e-3,
                np.array(self.i_save) * 1e-9)


##############################################################################

if __name__ == "__main__":

    import matplotlib.pyplot as plt
    import time

    start = time.time()

    stim_times = np.array([0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 1.0])
    tau = 10e-3
    tau_r = 200e-3
    tau_f = 900e-3
    u = 0.5

    tau = 0.1
    tau_r = 0.1
    tau_f = 0.1
    u = 0.1
    cond = 1e-7

    rlsr = RunLittleSynapseRun(stim_times, synapse_type="glut")

    for i in range(10):
        t_s, v_s, i_s = rlsr.run(tau * i, tau_r, tau_f, u, cond)
        plt.figure()
        plt.plot(t_s, v_s)
        plt.xlabel("Time (s)")
        plt.ylabel("Voltage (V)")
        plt.figure()
        plt.plot(t_s, i_s)
        plt.xlabel("Time (s)")
        plt.ylabel("Current (A)")

    # tS,vS = rlsr.run(tau*1,tauR,tauF,U,cond)
    # plt.plot(tS,vS)

    plt.xlabel("Time")
    plt.ylabel("Voltage")
    plt.title("How tau affects dynamics")

    end = time.time()
    print(end - start)

    plt.ion()
    plt.show()

    import pdb

    pdb.set_trace()
