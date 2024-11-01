import sys
import os.path

import neuron
import numpy as np

from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.simulate.nrn_simulator_parallel import NrnSimulatorParallel

import bluepyopt.ephys as ephys
from snudda.neurons.neuron_model_extended import NeuronModel
from snudda.neurons.neuron_morphology import NeuronMorphology


# Plot all sections
# [neuron.h.psection(x) for x in neuron.h.allsec()]

##############################################################################


class RunSynapseRun(object):

    def __init__(self,
                 stim_times,
                 synapse_density,
                 num_synapses,
                 neuron_parameter_key,
                 neuron_morphology_key,
                 neuron_path=None,
                 neuron_morphology=None,
                 neuron_mechanisms=None,
                 neuron_parameters=None,
                 # neuron_modulation=None,
                 synapse_section_id=None,  # if given, nSynapses is ignored
                 synapse_section_x=None,  # same # of elements as synapseSectionID
                 # neuron_parameter_id=0,  # Which param set in parameter file to use
                 # neuron_morphology_id=0,
                 # neuron_modulation_id=0,
                 holding_voltage=-70e-3,
                 holding_current=None,
                 synapse_type='glut',
                 params={},
                 time=2,
                 random_seed=None,
                 log_file=None,
                 verbose=True):

        self.log_file = log_file  # File pointer
        self.verbose = verbose
        self.rng = np.random.default_rng(random_seed)

        self.write_log(f"Holding voltage: {holding_voltage} V")
        self.write_log(f"Stim times: {stim_times} s")
        self.write_log(f"Synapse type: {synapse_type}")

        self.time = time
        self.synapses = []
        self.i_clamp = None
        self.nc_syn = None

        self.i_save = None
        self.t_save = None
        self.v_save = None

        self.stim_vector = None
        self.vec_stim = None
        self.synapse_section_id = None
        self.synapse_section_x = None

        if neuron_path is None and neuron_parameters is not None:
            neuron_path = os.path.dirname(neuron_parameters)

        self.neuron_path = neuron_path

        if neuron_parameters is None and self.neuron_path is not None:
            neuron_parameters = os.path.join(self.neuron_path, "parameters.json")

        if neuron_morphology is None and self.neuron_path is not None:
            neuron_morphology = os.path.join(self.neuron_path, "morphology")

        if neuron_mechanisms is None and self.neuron_path is not None:
            neuron_mechanisms = os.path.join(self.neuron_path, "mechanisms.json")

        # if neuron_modulation is None and self.neuron_path is not None:
        #     neuron_modulation = os.path.join(self.neuron_path, "modulation.json")

        # Done in NrnSimulatorParallel
        # neuron.h.load_file('stdrun.hoc')

        self.sim = NrnSimulatorParallel(cvode_active=False)

        # Should we use weak reference for garbage collection? (weakref package)

        # We load the neuron morphology object also, used to place synapses
        self.write_log(f"Using morphology: {neuron_morphology}")
        neuron_prototype = NeuronPrototype(neuron_path=neuron_path,
                                           neuron_name="OptimisationNeuron")
        self.morphology = neuron_prototype.clone(parameter_key=neuron_parameter_key,
                                                 morphology_key=neuron_morphology_key)
        # self.morphology = NeuronMorphology(swc_filename=neuron_morphology)

        # We need to setup the Neuron model
        self.neuron = NeuronModel(param_file=neuron_parameters,
                                  morph_path=neuron_morphology,
                                  mech_file=neuron_mechanisms,
                                  cell_name="OptimisationNeuron",
                                  # modulation_file=neuron_modulation,
                                  parameter_key=neuron_parameter_key,
                                  morphology_key=neuron_morphology_key)

        self.neuron.instantiate(sim=self.sim)
        self.set_resting_voltage(holding_voltage * 1e3)

        neuron.h.celsius = 35

        # gnabar_hh: The maximum specific sodium channel conductance [Default value = 0.120 S/cm2]
        # gkbar_hh: The maximum specific potassium channel conductance [Default value = 0.036 S/cm2]
        # gl_hh: The maximum specific leakage conductance [Default value = 0.0003 S/cm2]
        # ena: The reversal potential for the sodium channel [Default value = 50 mV]
        # ek: The reversal potential for the potassium channel [Default value = -77 mV]
        # el_hh: The reversal potential for the leakage channel [Default value = -54.3 mV]

        # We need to set the params also
        self.params = params
        self.default_cond = 5e-7

        # This returns (section,sectionX) so we can reuse it if needed
        self.synapse_positions = self.add_synapse_density(synapse_type=synapse_type,
                                                          synapse_density=synapse_density,
                                                          num_synapses=num_synapses,
                                                          section_id=synapse_section_id,
                                                          section_x=synapse_section_x)

        self.stim_times = np.array(stim_times)

        # Assumes input in seconds (SI units)
        self.connect_input_to_synapses(self.stim_times)

        self.soma_record()
        self.synapse_current_record()

        self.holding_current = self.update_holding_current(holding_voltage=holding_voltage,
                                                           holding_current=holding_current)

        self.write_log("RunSynapseRun: Init done.")

        # import pdb
        # pdb.set_trace()

    ############################################################################

    def __del__(self):

        # This should not be needed but...
        self.neuron = None
        self.morphology = None
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

    def update_holding_current(self, holding_voltage=None, holding_current=None):

        if holding_voltage is None:
            holding_voltage = self.holding_voltage
        else:
            self.holding_voltage = holding_voltage

        if holding_current is not None:
            if self.i_clamp is None:
                self.i_clamp = neuron.h.IClamp(self.neuron.icell.soma[0](0.5))

            # Update current on iClamp
            self.i_clamp.amp = holding_current * 1e9  # Convert SI -> nA for NEURON
            self.i_clamp.dur = 2 * self.time * 1e3

            self.set_resting_voltage(self.holding_voltage * 1e3)
            self.write_log(f"Set holding current {holding_current}A and holding voltage {holding_voltage}V,"
                           f" until {self.i_clamp.dur} ms")
            return holding_current

        self.write_log("Updating holding current, might take a bit of time")

        # Disable old iClamp temporarily
        if self.i_clamp is not None:
            self.i_clamp.amp = 0

        # Setup a temporary VClamp
        self.v_clamp = neuron.h.SEClamp(self.neuron.icell.soma[0](0.5))
        self.v_clamp.rs = 1e-9
        self.v_clamp.amp1 = holding_voltage * 1e3
        self.v_clamp.dur1 = self.time * 2 * 1e3
        # self.writeLog("VClamp duration: " + str(self.VClamp.dur1))

        neuron.h.finitialize(self.holding_voltage * 1e3)
        # !!! There is a WEIRD neuron bug, that if this tstop here is
        # different from duration of simulation, then the *SECOND* time
        # a model is initialised we get the length of tSave set by this
        # value, and not by the tStop of that simulation --- go figure!
        self.set_resting_voltage(self.holding_voltage * 1e3)

        neuron.h.tstop = self.time * 1e3  # Must set tstop
        neuron.h.run()

        if False:
            import matplotlib.pyplot as plt
            plt.figure()
            plt.plot(self.t_save, self.v_save)
            plt.title("Holding voltage should be " + str(self.holding_voltage * 1e3) + "mV")
            plt.xlabel("time (ms)")
            plt.ylabel("volt (mV)")
            plt.ion()
            plt.show()

            import pdb
            pdb.set_trace()

        cur = float(self.v_clamp.i)

        # Remove VClamp
        self.v_clamp = None

        if self.i_clamp is None:
            self.i_clamp = neuron.h.IClamp(self.neuron.icell.soma[0](0.5))

        # Update current on iClamp
        self.i_clamp.amp = cur
        self.i_clamp.dur = 2 * self.time * 1e3

        self.set_resting_voltage(self.holding_voltage * 1e3)

        self.write_log(f"Holding voltage {self.holding_voltage * 1e3} mV, IClamp amp = {cur} nA")

        return cur * 1e-9  # Convert to SI units

    ############################################################################

    def set_stim_times(self, stim_times):

        if len(stim_times) != len(self.stim_times) or (stim_times != self.stim_times).any():
            print(f"Setting stim times to {stim_times} s")
            self.write_log(f"Setting stim times to {stim_times} s")
            self.stim_vector = neuron.h.Vector(stim_times * 1e3)
            self.stim_times = stim_times * 1e3

    ############################################################################

    def add_synapse_density(self, synapse_type,
                            synapse_density,
                            num_synapses=None,
                            plot_flag=False,
                            section_id=None,
                            section_x=None):

        if plot_flag:

            assert section_id is None and section_x is None, \
                "Can not plot if sectionID and sectionX are given"

            input_coords, section_id, section_x, density_function, dist_syn_soma = \
                self.morphology.dendrite_input_locations(synapse_density_str=synapse_density,
                                                         num_locations=num_synapses,
                                                         return_density=True)

            self.synapse_locations = input_coords
            dist_from_soma = self.morphology.dend[:, 4]

            # plot density function
            plt.figure()
            plt.plot(dist_from_soma * 1e6, density_function, 'o')
            plt.xlabel('distance from soma $(\mu m)$')
            plt.title('density function')
            plt.ion()
            plt.show()

            # plot histogram - distance synapses from soma
            plt.figure()
            plt.hist(dist_syn_soma * 1e6, edgecolor='gray', bins=dist_syn_soma.size)
            plt.xlabel('distance from soma $(\mu m)$')
            plt.ylabel('frequency')
            plt.title('synaptic distance from soma')
            plt.ion()
            plt.show()

        elif section_id is None or section_x is None:

            input_coords, section_id, section_x, dist_syn_soma = \
                self.morphology.dendrite_input_locations(synapse_density_str=synapse_density,
                                                         num_locations=num_synapses,
                                                         rng=self.rng)

            self.synapse_locations = input_coords

        else:
            self.synapse_locations = "Unknown, you specified sectionX and sectionID, " \
                                     + "but not synapse xyz coordinates."

        # This is so we can find out where the synapses were placed
        self.synapse_section_id = section_id
        self.synapse_section_x = section_x

        sections = self.neuron.map_id_to_compartment(section_id)

        for s, sx in zip(sections, section_x):
            self.add_synapse(synapse_type, s, sx, self.params)

        # Return synapse position if we want to reuse them elsewhere
        return sections, section_x

    ############################################################################

    def add_synapse(self, synapse_type, section, section_x, params):

        section_x = np.maximum(section_x, 1e-6)  # Cant be 0 or 1

        try:
            if synapse_type.lower() == 'glut':
                syn = neuron.h.tmGlut_double(section(section_x))
            elif synapse_type.lower() == "gaba":
                syn = neuron.h.tmGabaA(section(section_x))
            else:
                assert False, f"Unknown synapse type: {synapse_type}"

            self.synapses.append(syn)

        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)

            self.write_log("Did you remember to run nrnivmodl first, to generate channels mod files?")
            sys.exit(-1)

        for p in params:

            # Conductance is set separately, it is a property of netcon
            if p == "cond":
                self.default_cond = params["cond"]
                continue

            val = self.si_to_natural_units(p, params[p])

            setattr(syn, p, val)
            self.write_log(f"Setting parameters: {p} = {val} (neuron natural units)")

    ############################################################################

    def connect_input_to_synapses(self, stim_times):

        self.write_log(f"Stimulation times (s): {stim_times}")
        self.nc_syn = []

        self.stim_vector = neuron.h.Vector(stim_times * 1e3)
        self.vec_stim = neuron.h.VecStim()
        self.vec_stim.play(self.stim_vector)

        for syn in self.synapses:
            ncs = neuron.h.NetCon(self.vec_stim, syn)
            ncs.delay = 0
            ncs.threshold = 0
            self.nc_syn.append(ncs)

    ############################################################################

    def soma_record(self):

        self.t_save = neuron.h.Vector()
        self.t_save.record(neuron.h._ref_t)

        self.v_save = neuron.h.Vector()
        self.v_save.record(self.neuron.icell.soma[0](0.5)._ref_v)

    ############################################################################

    def synapse_current_record(self):

        self.i_save = []

        for syn in self.synapses:
            i_rec = neuron.h.Vector()
            i_rec.record(syn._ref_i)
            self.i_save.append(i_rec)

    ############################################################################

    def run(self, tau, tau_r, tau_f, u, cond=None, time=None):

        assert False, "This is the old run method"

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
        for ncs in self.nc_syn:
            ncs.weight[0] = 1 * cond * 1e6

        for syn in self.synapses:
            syn.tau = tau * 1e3
            syn.tau_r = tau_r * 1e3
            syn.tau_f = tau_f * 1e3
            syn.u = u

        # print(self.littleSynapse.tau)

        # print("Initialise voltage to " + str(self.holdingVoltage*1e3) + " mV")
        neuron.h.finitialize(self.holding_voltage * 1e3)  # OLD : -70
        neuron.h.tstop = time * 1e3

        neuron.h.run()

        # self.tSave.resize()
        # self.vSave.resize()
        # self.iSave.resize()

        if np.array(self.t_save).shape != np.array(self.v_save).shape:
            self.write_log("THIS IS WRONG, should be same shape!!")
            self.write_log(f"size t = {np.array(self.t_save).shape}")
            self.write_log(f"size v = {np.array(self.v_save).shape}")
            import pdb
            pdb.set_trace()

        # print("Last V = " + str(self.vSave[-1]*1e-3))

        # Convert back to SI units
        return (np.array(self.t_save) * 1e-3,
                np.array(self.v_save) * 1e-3,
                np.array(self.i_save) * 1e-9)

    ############################################################################

    def set_resting_voltage(self, rest_volt):

        self.write_log("Setting resting voltage to %.3f mV" % rest_volt)

        soma = [x for x in self.neuron.icell.soma]
        axon = [x for x in self.neuron.icell.axon]
        dend = [x for x in self.neuron.icell.dend]

        cell = soma + axon + dend

        for sec in cell:
            for seg in sec.allseg():
                seg.v = rest_volt

    ############################################################################

    # I wish Neuron would use natural units...

    def si_to_natural_units(self, var_name, value):

        conv_factor = {"U": 1.0,
                       "tauR": 1e3,
                       "tauF": 1e3,
                       "cond": 1e6,
                       "tau": 1e3,
                       "nmda_ratio": 1.0,

                       "tau1_ampa": 1.0,  # Ilaria's file has ms already
                       "tau2_ampa": 1.0,  # Ilaria's file has ms already
                       "tau3_ampa": 1.0,  # Ilaria's file has ms already
                       "tau1_nmda": 1.0,  # Ilaria's file has ms already
                       "tau2_nmda": 1.0,  # Ilaria's file has ms already
                       "tau3_nmda": 1.0,  # Ilaria's file has ms already
                       "tpeak_ampa": 1.0,  # Ilaria's file has ms already
                       "tpeak_nmda": 1.0,  # Ilaria's file has ms already :/
                       "ratio_nmda": 1.0,

                       "I2_ampa": 1.0,  # Ilaria's file has ms already
                       "I3_ampa": 1.0,  # Ilaria's file has ms already
                       "I2_nmda": 1.0,  # Ilaria's file has ms already
                       "I3_nmda": 1.0,  # Ilaria's file has ms already
                       "factor_ampa": 1.0,  # Ilaria's file has ms already
                       "factor_nmda": 1.0  # Ilaria's file has ms already
                       }

        if var_name not in conv_factor:
            self.write_log("Missing conversion fractor for " + str(var_name) \
                           + ". Please update SItoNaturalUnits function.")
            self.write_log("convFactor = " + str(conv_factor))
            import pdb
            pdb.set_trace()

        try:
            return value * conv_factor[var_name]
        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)
            import pdb
            pdb.set_trace()

    ############################################################################
    def plot(self):

        ax = self.morphology.plot_neuron(axon_colour='red', dend_colour='blue', soma_colour='green')

        ax.scatter(self.synapse_locations[:, 0],
                   self.synapse_locations[:, 1],
                   self.synapse_locations[:, 2],
                   color='black',
                   marker='o',
                   s=20)
        plt.ion()

        ############################################################################

    # pars = { "tau1" : 0.25e-3 }
    # The code converts to natural units internally, if your parameter is missing
    # then update SItoNaturalUnits

    # OBS, soma parameters are ignored by run2 (they should be set during setup)

    def run2(self, pars, time=None, cond=1e-8):

        self.write_log(f"Running with pars: {pars}")

        if time is None:
            time = self.time
        else:
            self.time = time

        for p in pars:

            if p == "cond":
                cond = self.si_to_natural_units(p, pars[p])

            else:
                v = self.si_to_natural_units(p, pars[p])
                for s in self.synapses:
                    setattr(s, p, v)

        neuron.h.finitialize(self.holding_voltage * 1e3)
        for ncs in self.nc_syn:
            ncs.weight[0] = cond

        self.set_resting_voltage(self.holding_voltage * 1e3)

        neuron.h.v_init = self.holding_voltage * 1e3
        neuron.h.tstop = time * 1e3
        self.write_log("About to start NEURON... stay safe")
        neuron.h.run()
        self.write_log("NEURON actually completed?!")

        # Convert results back to SI units
        return (np.array(self.t_save) * 1e-3,
                np.array(self.v_save) * 1e-3,
                np.array(self.i_save) * 1e-9)

    ############################################################################

    def write_log(self, text, flush=True):  # Change flush to False in future, debug
        if self.log_file is not None:
            self.log_file.write(text + "\n")

            if self.verbose:
                print(text)

            if flush:
                self.log_file.flush()
        else:
            if self.verbose:
                print(text)

            ##############################################################################


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    import time

    start = time.time()

    stim_times = np.array([0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 1.0]) + 0.3

    # synapseDensity="0.05/(1+np.exp(-(d-30e-6)/5e-6))"
    # synapseDensity="0.05/(1+np.exp(-(d-150e-6)/5e-6))"

    synapse_density = "np.exp(-((d-50e-6)/15e-6)**2)"
    holding_voltage = -79.5e-3

    # !!! OBS changed number of synapses, before was 20 -- cond 1nS
    tau = 10e-3
    tau_r = 200e-3
    tauF = 900e-3
    U = 0.5

    tau = 0.1
    tau_r = 0.1
    tauF = 0.1
    U = 0.1
    cond = 1e-9

    rlsr = RunSynapseRun(stim_times=stim_times,
                         synapse_type="glut",
                         neuron_morphology="../data/cellspecs/dspn/str-dspn-e150917_c9_d1-mWT-1215MSN03-v20190521/WT-1215MSN03-cor-rep-ax2.swc",
                         neuron_parameters="../data/cellspecs/dspn/str-dspn-e150917_c9_d1-mWT-1215MSN03-v20190521/parameters.json",
                         neuron_mechanisms="../data/cellspecs/dspn/str-dspn-e150917_c9_d1-mWT-1215MSN03-v20190521/mechanisms.json",
                         synapse_density=synapse_density,
                         num_synapses=10,
                         neuron_parameter_key=0,
                         # neuron_modulation_key=0,
                         holding_voltage=holding_voltage)

    # I would like to plot the morphology
    # ax = rlsr.plot()

    # import pdb
    # pdb.set_trace()

    params = {"U": 0.5,
              "tauR": 0.6,
              "tauF": 2.0,  # 1.1
              "tau": 0.08,
              "cond": 1e-9,
              "nmda_ratio": 0.5}

    tS, vS, iS = rlsr.run2(pars=params)
    plt.figure()
    plt.plot(tS, vS)
    plt.xlabel("Time (s)")
    plt.ylabel("Voltage (V)")

    if True:
        expData = 'DATA/YvonneJohansson2019/M1LH_Analysis_191001.h5'
        expCellID = 131

        import LoadExpData

        led = LoadExpData.LoadExpData(expData)
        cellType = led.get_cell_type(expCellID)
        plt.title("ExpData is cellType : " + str(cellType))
        expV, expT = led.getData("GBZ_CC_H20", expCellID)
        plt.plot(expT, expV, 'r-')

    if False:
        plt.figure()

        for iiS in iS:
            plt.plot(tS, iiS)

        plt.xlabel("Time (s)")
        plt.ylabel("Current (A)")

    plt.ion()
    plt.show()

    end = time.time()
    print(end - start)

    import pdb

    pdb.set_trace()
