import neuron
import numpy as np

from snudda.simulate.nrn_simulator_parallel import NrnSimulatorParallel
from snudda.neurons.neuron_prototype import NeuronPrototype
import bluepyopt.ephys as ephys
from snudda.neurons.neuron_model_extended import NeuronModel
from snudda.neurons.neuron_morphology import NeuronMorphology
import matplotlib.pyplot as plt


# Plot all sections
# [neuron.h.psection(x) for x in neuron.h.allsec()]

##############################################################################


class SynapseClass:

    def __init__(self,
                 neuron_morphology,
                 neuron_mechanisms,
                 neuron_parameters,
                 neuron_modulation,
                 stim_times=None,
                 synapse_density=None,
                 num_synapses=None,
                 synapse_section_id=None,  # if given, nSynapses is ignored
                 synapse_section_x=None,  # same # of elements as synapseSectionID
                 neuron_parameter_id=0,  # Which param set in parameter file to use
                 neuron_modulation_id=0,
                 holding_voltage=-70e-3,
                 holding_current=None,
                 synapse_type=None,
                 params={},
                 time=2,
                 random_seed=None,
                 log_file=None,
                 verbose=True):

        self.synapse_type = synapse_type
        self.stim_times = stim_times
        self.synapse_density = synapse_density
        self.num_synapses = num_synapses

        self.neuron_parameter_id = neuron_parameter_id
        self.neuron_modulation_id = neuron_modulation_id
        self.neuron_morphology = neuron_morphology
        self.neuron_mechanisms = neuron_mechanisms
        self.neuron_parameters = neuron_parameters
        self.neuron_modulation = neuron_modulation
        self.holding_voltage = holding_voltage
        self.holding_current = holding_current
        self.log_file = log_file  # File pointer
        self.verbose = verbose
        self.rng = np.random.default_rng(random_seed)
        self.params = params
        self.default_cond = 5e-7
        self.size_of_noise = 0.1

        self.write_log("Holding voltage: " + str(holding_voltage) + " V")
        self.write_log("Stim times: " + str(stim_times) + " s")
        self.write_log("Synapse type: " + str(synapse_type))

        self.time = time
        self.internal_synapses = []
        self.external_synapses = []
        self.vec_list = list()
        self.i_clamp = None

        self.i_save = None
        self.t_save = None
        self.v_save = None
        self.nc_syn =list()

        self.stim_vector = None
        self.vec_stim = None
        self.synapse_locations = list()

        # Done in NrnSimulatorParallel
        # neuron.h.load_file('stdrun.hoc')

        self.sim = NrnSimulatorParallel(cvode_active=False)
        self.sim.neuron.h.celsius = 35

        # Should we use weak reference for garbage collection? (weakref package)

        # We load the neuron morphology object also, used to place synapses
        #self.write_log(f"Using morphology: {neuron_morphology}")
        #self.morphology = NeuronMorphology(swc_filename=neuron_morphology)
        #

        self.nm = NeuronPrototype(neuron_name="JJJ",
                                         neuron_path=None,
                                        morphology_path=self.neuron_morphology,
                                    parameter_path=self.neuron_parameters,
                             mechanism_path=self.neuron_mechanisms,
                             virtual_neuron=False,
                             axon_stump_id_flag=False)
        self.nm.instantiate()

        swc_file = self.nm.get_morphology(self.neuron_parameter_id,0)
        self.morphology = NeuronMorphology(swc_filename=swc_file, verbose=True, use_cache=False)


        # We need to setup the Neuron model
        self.neuron = NeuronModel(param_file=self.neuron_parameters,
                                  morph_path=self.neuron_morphology,
                                  mech_file=self.neuron_mechanisms,
                                  cell_name="OptimisationNeuron",
                                  modulation_file=self.neuron_modulation,
                                  morphology_id=0,
                                  parameter_id=self.neuron_parameter_id,
                                  modulation_id=self.neuron_modulation_id)
    def return_external_synapses(self):

        return self.external_synapses

    def return_internal_synapses(self):

        return self.internal_synapses

    def old_setup(self, synapse_section_id, synapse_section_x):

        self.instantiate()

        self.set_rest_voltage(self.holding_voltage)

        self.add_synapse_distributions(self.synapse_type,self.synapse_density,self.num_synapses,synapse_section_id,
                                       synapse_section_x)

        self.add_stim_times_to_synapse_distributions(self.stim_times)

        self.record_voltage()
        self.record_currents()

        self.set_holding_current(self.holding_voltage,self.holding_current)

        self.write_log("RunSynapseRun: Init done.")

    def instantiate(self):

        self.neuron.instantiate(sim=self.sim)

    def set_rest_voltage(self,holding_voltage):


        self.set_resting_voltage(holding_voltage * 1e3)


    def add_synapse_distributions(self,synapse_type,synapse_density=None,num_synapses=None,synapse_section_id=None,synapse_section_x=None):

        _, _,synapses, synapse_positions = self.add_synapse_density(synapse_type=synapse_type,
                                                          synapse_density=synapse_density,
                                                          num_synapses=num_synapses,
                                                          section_id=synapse_section_id,
                                                          section_x=synapse_section_x)
        self.internal_synapses = synapses

        return synapse_positions

    def add_stim_times_to_synapse_distributions(self, stim_times,synapses,weight):

        stim_times = np.array(stim_times)

        # Assumes input in seconds (SI units)
        with_noise = self.connect_input_to_synapses(stim_times,synapses,weight)

        return with_noise

    def record_voltage(self):

        self.soma_record()

    def record_currents(self):

        self.synapse_current_record()

    def set_holding_current(self, holding_voltage, holding_current):

        self.holding_current = self.update_holding_current(holding_voltage=holding_voltage,
                                                           holding_current=holding_current)

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
            self.write_log(f"Set holding current {holding_current}A and holding voltage {holding_voltage}V")
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
                            synapse_density=None,
                            num_synapses=None,
                            section_id=None,
                            section_x=None):

        synapse_locations = list()

        if section_id is None or section_x is None:

            input_coords, section_id, section_x, dist_syn_soma = \
                self.morphology.dendrite_input_locations(synapse_density=synapse_density,
                                                         num_locations=num_synapses,
                                                         rng=self.rng)

            synapse_locations = list(input_coords)

        '''
        Do we need to save it???
        # This is so we can find out where the synapses were placed
        self.section_id.append(section_id)
        self.synapse_section_x = section_x
        '''

        sections = self.neuron.map_id_to_compartment(section_id)

        synapses = list()

        for s, sx in zip(sections, section_x):

            number_along = int(np.floor(s.n3d() * s(sx).x) - 1)

            number_along = 0 if number_along < 0 else number_along

            input_coords = np.array([s.x3d(number_along),s.y3d(number_along),s.z3d(number_along)])

            syn = self.add_synapse(synapse_type, s, sx)
            synapse_locations.append(np.array(input_coords)*1e-6)

            synapses.append(syn)

        # Return synapse position if we want to reuse them elsewhere

        return sections, section_x, synapses, synapse_locations

    ############################################################################

    def add_synapse(self, synapse_type, section, section_x):

        section_x = np.maximum(section_x, 1e-6)  # Cant be 0 or 1

        try:
            if synapse_type.lower() == 'glut':
                syn = neuron.h.tmGlut_double(section(section_x))
            elif synapse_type.lower() == 'single_glut':
                syn = neuron.h.tmGlut(section(section_x))
            elif synapse_type.lower() == "gaba":
                syn = neuron.h.tmGabaA(section(section_x))
            else:
                assert False, f"Unknown synapse type: {synapse_type}"

        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)

            self.write_log("Did you remember to run nrnivmodl first, to generate channels mod files?")
            sys.exit(-1)

        return syn


    ############################################################################

    def set_size_of_noise(self,size):

        self.size_of_noise = size
    def connect_input_to_synapses(self, stim_times,synapses,weight):

        self.write_log(f"Stimulation times Original (s): {stim_times}")

        print(f"Size of noise is {self.size_of_noise}")
        stim_times_with_noise = list()
        for syn in synapses:
            noise = np.random.normal(1, self.size_of_noise, stim_times.shape) * 1e2

            stim_times_n = stim_times * 1e3 + noise.astype(int)
            stim_times_n = stim_times_n.astype(int)
            stim_times_with_noise.append(stim_times_n)

            if not(all(i >= 0 for i in stim_times_n)):
                print(stim_times_n)
                raise ValueError('Change size of noise or move start times - '
                                 'some values are negative')
            stim_vector = neuron.h.Vector(list(stim_times_n))
            vec_stim = neuron.h.VecStim()
            vec_stim.play(stim_vector)
            ncs = neuron.h.NetCon(vec_stim, syn)
            ncs.delay = 0
            ncs.threshold = 0
            ncs.weight[0] = weight
            self.nc_syn.append(ncs)
            self.vec_list.append([stim_vector,vec_stim])

        return stim_times_with_noise

    ############################################################################

    def soma_record(self):

        self.t_save = neuron.h.Vector()
        self.t_save.record(neuron.h._ref_t)

        self.v_save = neuron.h.Vector()
        self.v_save.record(self.neuron.icell.soma[0](0.5)._ref_v)

    ############################################################################

    def synapse_current_record(self,synapses):

        self.i_save = []

        for syn in synapses:
            i_rec = neuron.h.Vector()
            i_rec.record(syn._ref_i)
            self.i_save.append(i_rec)


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
                       "tau1" : 1e3,
                       "tau2" : 1e3,
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
    def plot_synapse_locations(self, synapse_locations,colour=None):

        synapse_locations = np.array(synapse_locations)
        ax = self.morphology.plot_neuron(axon_colour='grey', dend_colour='blue', soma_colour='green')

        ax.scatter(synapse_locations[:, 0],
                   synapse_locations[:, 1],
                   synapse_locations[:, 2],
                   color=colour,
                   marker='o',
                   s=20)
        plt.show()

    def plot_interal_and_external_synapses(self,internal,external,colors=None):

        internal = np.array(internal)
        ax = self.morphology.plot_neuron(axon_colour='grey', dend_colour='blue', soma_colour='green')

        ax.scatter(internal[:, 0],
                   internal[:, 1],
                   internal[:, 2],
                   color=colors[0],
                   marker='o',
                   s=20)

        external = np.array(external)
        ax.scatter(external[:, 0],
                   external[:, 1],
                   external[:, 2],
                   color=colors[1],
                   marker='o',
                   s=20)

        plt.show()

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
                         neuron_modulation="../data/cellspecs/dspn/str-dspn-e150917_c9_d1-mWT-1215MSN03-v20190521/modulation.json",
                         synapse_density=synapse_density,
                         num_synapses=10,
                         neuron_parameter_id=0,
                         neuron_modulation_id=0,
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
