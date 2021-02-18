# This code writes the input spikes for the neuron simulation --
#
#
# If nInputs is given then synapseDensity is scaled to give approximately
# that total number of synapses, otherwise it is used without scaling.
# see config/input-tinytest-v2.json for example config.
#

#
# !!!! Change how data is stored, many small datasets is inefficient
#

# Smith, Galvan, ..., Bolam 2014 -- Bra info om thalamic inputs, CM/PF
#


import numpy as np
import h5py
import json
import os
from glob import glob
import itertools

import matplotlib.pyplot as plt

from .neuron_morphology import NeuronMorphology
from .load import SnuddaLoad

nl = None

# When specifying vectors for start and end time, they should normally not overlap
# if we want to allow that, set time_interval_overlap_warning = False


class SnuddaInput(object):

    def __init__(self, spike_data_filename, input_config_file,
                 network_config_file=None,
                 position_file=None,
                 hdf5_network_file=None,
                 time=10.0,
                 is_master=True,
                 h5libver="latest",
                 rc=None,
                 random_seed=None,
                 time_interval_overlap_warning=True,
                 logfile=None,
                 verbose=True):

        if type(logfile) == str:
            self.logfile = open(logfile, "w")
        else:
            self.logfile = logfile

        self.verbose = verbose
        self.rc = rc

        self.time_interval_overlap_warning = time_interval_overlap_warning
        self.input_info = None
        self.population_unit_spikes = None
        self.network_info = None
        self.neuron_info = None
        self.num_population_units = None
        self.population_unit_id = None
        self.neuron_name = None
        self.neuron_id = None
        self.neuron_type = None
        self.d_view = None
        self.lb_view = None
        self.rc = None
        self.network_slurm_id = None
        self.network_config = None
        self.neuron_input = None
        self.slurm_id = None

        self.write_log("Time = " + str(time))

        self.random_seed = random_seed

        self.h5libver = h5libver
        self.write_log("Using hdf5 version " + str(h5libver))

        if hdf5_network_file is not None:
            assert network_config_file is None and position_file is None, \
                "If HDF5networkFile specified then positionFile " + \
                "and networkConfigFile should be left empty."

            if hdf5_network_file == "last":
                hdf5_network_file = self.find_latest_file()

            self.hdf5_network_file = hdf5_network_file
            self.read_hdf5_info(hdf5_network_file)
        else:
            self.network_config_file = network_config_file
            self.position_file = position_file

            # self.writeLog("Assuming axonStumpIDFlag is True (Running Network_simulate.py)")
            self.axon_stump_id_flag = False

        self.input_config_file = input_config_file

        if spike_data_filename is None:
            spike_data_filename = "save/input-spikes-" + str(self.network_slurm_id) + ".hdf5"
        self.spike_data_filename = spike_data_filename

        self.time = time  # How long time to generate inputs for

        self.neuron_cache = dict([])

        self.is_master = is_master

    def generate(self):

        # Read in the input configuration information from JSON file
        self.read_input_config_file()

        # Read the network position file
        self.read_neuron_positions()

        # Read the network config file -- This also reads random seed
        self.read_network_config_file()

        # Only the master node should start the work
        if self.is_master:
            # Initialises lbView and dView (load balance, and direct view)
            self.setup_parallel()

            # Make the "master input" for each channel
            # TODO: Having common spikes within a population unit caused problems, a lot of neurons fired at same
            # time. To avoid this the population_spike_list is cleared in make_neuron_input_parallel
            # Should rewrite this later.
            rng = self.get_master_node_rng()
            self.make_population_unit_spike_trains(rng=rng)

            # Generate the actual input spikes, and the locations
            # stored in self.neuronInput dictionary

            self.make_neuron_input_parallel()

            # Write spikes to disk, HDF5 format
            self.write_hdf5()

            # Verify correlation --- THIS IS VERY VERY SLOW
            # self.verifyCorrelation()

            self.check_sorted()

        # !!! TODO

        # 1. Define what the within correlation, and between correlation should be
        #    for each neuron type. Also what input frequency should we have for each
        #    neuron. --- Does it depend on size of dendritic tree?
        #    Store the info in an internal dict.

        # 2. Read the position file, so we know what neurons are in the network

        # 3. Create the "master input" for each population unit.

        # 4. Mix the master input with random input, for each neuron, to create
        #    the appropriate correlations

        # 5. Randomize which compartments each synaptic input should be on

        # 6. Verify correlation of input

        # 7. Write to disk

        # If more than one worker node, then we need to split the data
        # into multiple files
        # self.nWorkers=nWorkers

    ############################################################################

    def write_hdf5(self):

        self.write_log("Writing spikes to " + self.spike_data_filename)

        out_file = h5py.File(self.spike_data_filename, 'w', libver=self.h5libver)
        config_data = out_file.create_dataset("config", data=json.dumps(self.input_info, indent=4))
        input_group = out_file.create_group("input")

        for neuron_id in self.neuron_input:

            nid_group = input_group.create_group(str(neuron_id))

            neuron_type = self.neuron_type[neuron_id]
            # nName = self.neuronName[neuronID]

            for input_type in self.neuron_input[neuron_id]:

                if input_type[0] == '!':
                    self.write_log(f"Disabling input {input_type} for neuron {neuron_id} "
                                   f" (input_type was commented with ! before name)")
                    continue

                if input_type.lower() != "VirtualNeuron".lower():
                    it_group = nid_group.create_group(input_type)

                    neuron_in = self.neuron_input[neuron_id][input_type]
                    spike_mat, num_spikes = self.create_spike_matrix(neuron_in["spikes"])

                    it_group.create_dataset("spikes", data=spike_mat)
                    it_group.create_dataset("nSpikes", data=num_spikes)

                    it_group.create_dataset("sectionID", data=neuron_in["location"][1])
                    it_group.create_dataset("sectionX", data=neuron_in["location"][2])
                    it_group.create_dataset("distanceToSoma", data=neuron_in["location"][3])

                    it_group.create_dataset("freq", data=neuron_in["freq"])
                    it_group.create_dataset("correlation", data=neuron_in["correlation"])
                    it_group.create_dataset("jitter", data=neuron_in["jitter"])
                    it_group.create_dataset("synapseDensity",
                                            data=neuron_in["synapseDensity"])
                    it_group.create_dataset("start", data=neuron_in["start"])
                    it_group.create_dataset("end", data=neuron_in["end"])
                    it_group.create_dataset("conductance", data=neuron_in["conductance"])

                    population_unit_id = neuron_in["populationUnitID"]
                    it_group.create_dataset("populationUnitID", data=population_unit_id)

                    # TODO: What to do with population_unit_spikes, should we have mandatory jittering for them?
                    if neuron_type in self.population_unit_spikes:
                        chan_spikes = self.population_unit_spikes[neuron_type][input_type][population_unit_id]
                    else:
                        chan_spikes = np.array([])
                    it_group.create_dataset("populationUnitSpikes", data=chan_spikes)

                    it_group.create_dataset("generator", data=neuron_in["generator"])

                    it_group.create_dataset("modFile", data=neuron_in["modFile"])
                    if neuron_in["parameterFile"]:
                        it_group.create_dataset("parameterFile", data=neuron_in["parameterFile"])
                    # We need to convert this to string to be able to save it
                    it_group.create_dataset("parameterList", data=json.dumps(neuron_in["parameterList"]))
                    it_group.create_dataset("parameterID", data=neuron_in["parameterID"])

                else:

                    # Input is activity of a virtual neuron
                    a_group = nid_group.create_group("activity")
                    spikes = self.neuron_input[neuron_id][input_type]["spikes"]

                    a_group.create_dataset("spikes", data=spikes)
                    generator = self.neuron_input[neuron_id][input_type]["generator"]
                    a_group.create_dataset("generator", data=generator)

        out_file.close()

    ############################################################################

    @staticmethod
    def create_spike_matrix(spikes):

        if len(spikes) == 0:
            return np.zeros((0, 0)), 0

        num_input_trains = len(spikes)
        num_spikes = np.array([len(x) for x in spikes])
        max_len = max(num_spikes)

        spike_mat = -1 * np.ones((num_input_trains, max_len))
        for idx, st in enumerate(spikes):
            n = st.shape[0]
            spike_mat[idx, :n] = st

        return spike_mat, num_spikes

    ############################################################################

    # Reads from self.inputConfigFile

    def read_input_config_file(self):

        self.write_log("Loading input configuration from " + str(self.input_config_file))

        with open(self.input_config_file, 'rt') as f:
            self.input_info = json.load(f)

        for neuron_type in self.input_info:
            for input_type in self.input_info[neuron_type]:
                if "parameterFile" in self.input_info[neuron_type][input_type]:
                    par_file = self.input_info[neuron_type][input_type]["parameterFile"]

                    # Allow user to use $DATA to refer to snudda data directory
                    par_file = par_file.replace("$DATA", os.path.join(os.path.dirname(__file__), "data"))

                    with open(par_file, 'r') as f:
                        par_data_dict = json.load(f)

                    # Read in parameters into a list
                    par_data = []
                    for pd in par_data_dict:
                        par_data.append(par_data_dict[pd])
                else:
                    par_data = None

                self.input_info[neuron_type][input_type]["parameterList"] = par_data

    ############################################################################

    # Each synaptic input will contain a fraction of population unit spikes, which are
    # taken from a stream of spikes unique to that particular population unit
    # This function generates these correlated spikes

    def make_population_unit_spike_trains(self, rng, num_population_units=None):

        self.write_log("Running makePopulationUnitSpikeTrains")

        if num_population_units is None:
            num_population_units = self.num_population_units

        self.population_unit_spikes = dict([])

        for cell_type in self.input_info:

            self.population_unit_spikes[cell_type] = dict([])

            for input_type in self.input_info[cell_type]:

                if self.input_info[cell_type][input_type]["generator"] == "poisson":

                    freq = self.input_info[cell_type][input_type]["frequency"]
                    self.population_unit_spikes[cell_type][input_type] = dict([])

                    if "start" in self.input_info[cell_type][input_type]:
                        start_time = self.input_info[cell_type][input_type]["start"]
                    else:
                        start_time = 0

                    if "end" in self.input_info[cell_type][input_type]:
                        end_time = self.input_info[cell_type][input_type]["end"]
                    else:
                        end_time = self.time

                    if "populationUnitID" in self.input_info[cell_type][input_type]:
                        pop_unit_list = \
                            self.input_info[cell_type][input_type]["populationUnitID"]

                        if type(pop_unit_list) != list:
                            pop_unit_list = [pop_unit_list]
                    else:
                        pop_unit_list = range(0, num_population_units)

                    for idxPopUnit in pop_unit_list:
                        self.population_unit_spikes[cell_type][input_type][idxPopUnit] = \
                            self.generate_spikes(freq=freq, time_range=(start_time, end_time), rng=rng)

        return self.population_unit_spikes

    ############################################################################

    def make_neuron_input_parallel(self):

        self.write_log("Running makeNeuronInputParallel")

        self.neuron_input = dict([])

        neuron_id_list = []
        input_type_list = []
        freq_list = []
        start_list = []
        end_list = []
        synapse_density_list = []
        num_inputs_list = []
        p_keep_list = []
        population_unit_spikes_list = []
        jitter_dt_list = []
        population_unit_id_list = []
        conductance_list = []
        correlation_list = []

        mod_file_list = []
        parameter_file_list = []
        parameter_list_list = []

        for (neuron_id, neuron_name, neuron_type, populationUnitID) \
                in zip(self.neuron_id, self.neuron_name, self.neuron_type, self.population_unit_id):

            self.neuron_input[neuron_id] = dict([])

            # The input can be specified using neuron_id, neuron_name or neuron_type
            if str(neuron_id) in self.input_info:
                input_info = self.input_info[str(neuron_id)]
            elif neuron_name in self.input_info:
                input_info = self.input_info[neuron_name]
            elif neuron_type in self.input_info:
                input_info = self.input_info[neuron_type]
            else:
                self.write_log(f"!!! Warning, no synaptic input for neuron ID {neuron_id}, "
                               f"name {neuron_name} or type {neuron_type}")
                continue

            # if a number --> use a specific neuron with that given ID
            # if dSPN --> use neuron_type dSPN
            # if dSPN_3 --> use specific neuron morphology corresponding to dSPN_3

            for input_type in input_info:

                if input_type[0] == '!':
                    self.write_log(f"Disabling input {input_type} for neuron {neuron_name} "
                                   f" (input_type was commented with ! before name)")
                    continue

                input_inf = input_info[input_type]

                if "populationUnitID" in input_inf:
                    pop_unit_id = input_inf["populationUnitID"]

                    if type(pop_unit_id) == list and populationUnitID not in pop_unit_id:
                        # We have a list of functional channels, but this neuron
                        # does not belong to a functional channel in that list
                        continue
                    elif populationUnitID != pop_unit_id:
                        # We have a single functional channel, but this neuron is not
                        # in that functional channel
                        continue

                self.neuron_input[neuron_id][input_type] = dict([])

                if input_inf["generator"] == "poisson":
                    neuron_id_list.append(neuron_id)
                    input_type_list.append(input_type)
                    freq_list.append(input_inf["frequency"])

                    if "populationUnitCorrelation" in input_inf:
                        p_keep_list.append(np.sqrt(input_inf["populationUnitCorrelation"]))
                    else:
                        p_keep_list.append(0)

                    if "jitter" in input_inf:
                        jitter_dt_list.append(input_inf["jitter"])
                    else:
                        jitter_dt_list.append(None)

                    if "start" in input_inf:
                        start_list.append(input_inf["start"])
                    else:
                        start_list.append(0.0)  # Default start at beginning

                    if "end" in input_inf:
                        end_list.append(input_inf["end"])
                    else:
                        end_list.append(self.time)

                    if input_type.lower() == "VirtualNeuron".lower():
                        # Virtual neurons spikes specify their activity, location and conductance not used
                        cond = None
                        n_inp = 1

                        mod_file = None
                        parameter_file = None
                        parameter_list = None
                    else:
                        assert "location" not in input_inf, \
                            "Location in input config has been replaced with synapseDensity"
                        cond = input_inf["conductance"]

                        if "nInputs" in input_inf:

                            # If a dictionary, then extract the info for the relevant neuron
                            if type(input_inf["nInputs"]) == dict:
                                if neuron_name in input_inf["nInputs"]:
                                    n_inp = input_inf["nInputs"][neuron_name]
                                elif neuron_type in input_inf["nInputs"]:
                                    n_inp = input_inf["nInputs"][neuron_type]
                                else:
                                    n_inp = None
                            else:
                                n_inp = input_inf["nInputs"]

                        else:
                            n_inp = None

                        mod_file = input_inf["modFile"]
                        if type(mod_file) in [bytes, np.bytes_]:
                            mod_file = mod_file.decode()

                        if "parameterFile" in input_inf:
                            parameter_file = input_inf["parameterFile"]
                        else:
                            parameter_file = None

                        if "parameterList" in input_inf:
                            parameter_list = input_inf["parameterList"]
                        else:
                            parameter_list = None

                    if "synapseDensity" in input_inf:
                        synapse_density = input_inf["synapseDensity"]
                    else:
                        synapse_density = "1"

                    synapse_density_list.append(synapse_density)
                    num_inputs_list.append(n_inp)

                    population_unit_id_list.append(populationUnitID)
                    conductance_list.append(cond)

                    if "populationUnitCorrelation" in input_inf:
                        correlation_list.append(input_inf["populationUnitCorrelation"])
                    else:
                        correlation_list.append(0)

                    if neuron_type in self.population_unit_spikes:
                        # TODO: These are cleared at the end anyway, so not currently used. Remove altogether?
                        c_spikes = self.population_unit_spikes[neuron_type][input_type][populationUnitID]
                        population_unit_spikes_list.append(c_spikes)
                    else:
                        print(f"No population spikes specified for neuron type {neuron_type}")
                        population_unit_spikes_list.append(np.array([]))

                    mod_file_list.append(mod_file)
                    parameter_file_list.append(parameter_file)
                    parameter_list_list.append(parameter_list)

                elif input_inf["generator"] == "csv":
                    csv_file = input_inf["csvFile"] % neuron_id

                    self.neuron_input[neuron_id][input_type]["spikes"] \
                        = np.genfromtxt(csv_file, delimiter=',')
                    self.neuron_input[neuron_id][input_type]["generator"] = "csv"

                else:
                    self.write_log("Unknown input generator: " + input_inf["generator"]
                                   + " for " + str(neuron_id))

        seed_list = self.generate_seeds(num_states=len(neuron_id_list))

        # The old code had so that all neurons within a population unit shared the same
        # mother process, which caused them all to activate at the same time
        # with high probability. By setting channelSpikeList to None we disable it
        self.write_log(
            "Clearing populationUnitSpikesList, thus all neurons will have their own mother process for each input")
        population_unit_spikes_list = [None for x in population_unit_spikes_list]
        amr = None

        # Lets try and swap self.lbView for self.dView
        if self.d_view is not None:

            # self.writeLog("Sending jobs to workers, using lbView")
            self.write_log("Sending jobs to workers, using dView")

            # Changed the logic, the old input helper needed a global
            # variable to be visible, but it was not always so in its scope

            input_list = list(zip(neuron_id_list,
                                  input_type_list,
                                  freq_list,
                                  start_list,
                                  end_list,
                                  synapse_density_list,
                                  num_inputs_list,
                                  p_keep_list,
                                  population_unit_spikes_list,
                                  jitter_dt_list,
                                  population_unit_id_list,
                                  conductance_list,
                                  correlation_list,
                                  mod_file_list,
                                  parameter_file_list,
                                  parameter_list_list,
                                  seed_list))

            self.d_view.scatter("input_list", input_list, block=True)
            cmd_str = "inpt = list(map(nl.make_input_helper_parallel,input_list))"
            self.d_view.execute(cmd_str, block=True)

            inpt = self.d_view["inpt"]
            amr = list(itertools.chain.from_iterable(inpt))

        else:
            # If no lbView then we run it in serial
            self.write_log("Running input generation in serial")
            amr = map(self.make_input_helper_serial,
                      neuron_id_list,
                      input_type_list,
                      freq_list,
                      start_list,
                      end_list,
                      synapse_density_list,
                      num_inputs_list,
                      p_keep_list,
                      population_unit_spikes_list,
                      jitter_dt_list,
                      population_unit_id_list,
                      conductance_list,
                      correlation_list,
                      mod_file_list,
                      parameter_file_list,
                      parameter_list_list,
                      seed_list)

        # Gather the spikes that were generated in parallel
        for neuron_id, input_type, spikes, loc, synapse_density, frq, \
            jdt, p_uid, cond, corr, timeRange, \
                mod_file, paramFile, paramList, paramID in amr:
            self.write_log("Gathering " + str(neuron_id) + " - " + str(input_type))
            self.neuron_input[neuron_id][input_type]["spikes"] = spikes

            if input_type.lower() != "VirtualNeuron".lower():
                # Virtual neurons have no location of their input, as the "input"
                # specifies the spike times of the virtual neuron itself
                self.neuron_input[neuron_id][input_type]["location"] = loc
                self.neuron_input[neuron_id][input_type]["synapseDensity"] = synapse_density
                self.neuron_input[neuron_id][input_type]["conductance"] = cond

            self.neuron_input[neuron_id][input_type]["freq"] = frq
            self.neuron_input[neuron_id][input_type]["correlation"] = corr
            self.neuron_input[neuron_id][input_type]["jitter"] = jdt
            self.neuron_input[neuron_id][input_type]["start"] = timeRange[0]
            self.neuron_input[neuron_id][input_type]["end"] = timeRange[1]
            self.neuron_input[neuron_id][input_type]["populationUnitID"] = p_uid

            assert p_uid == self.population_unit_id[neuron_id], \
                "Internal error: Neuron should belong to the functional channel " \
                + "that input is generated for"

            self.neuron_input[neuron_id][input_type]["generator"] = "poisson"
            self.neuron_input[neuron_id][input_type]["modFile"] = mod_file
            self.neuron_input[neuron_id][input_type]["parameterFile"] = paramFile
            self.neuron_input[neuron_id][input_type]["parameterList"] = paramList
            self.neuron_input[neuron_id][input_type]["parameterID"] = paramID

        return self.neuron_input

    ############################################################################

    def generate_spikes_helper(self, frequencies, time_ranges, rng):
        t_spikes = []

        for f, t_start, t_end in zip(frequencies, time_ranges[0], time_ranges[1]):
            t_spikes.append(self.generate_spikes(f, (t_start, t_end), rng))

        # Double check correct dimension
        return np.sort(np.concatenate(t_spikes))

    def generate_spikes(self, freq, time_range, rng):
        # This generates poisson spikes with frequency freq, for a given time range

        if type(time_range[0]) == list:

            if type(freq) != list:
                freq_list = [freq for t in time_range[0]]
                freq = freq_list

            assert len(time_range[0]) == len(time_range[1]) == len(freq), \
                (f"Frequency, start and end time vectors need to be of same length." 
                    f"\nfreq: {freq}\nstart: {time_range[0]}\nend:{time_range[1]}")

            if self.time_interval_overlap_warning:
                assert (np.array(time_range[0][1:]) - np.array(time_range[1][0:-1]) >= 0).all(), \
                    f"Time range should not overlap: start: {time_range[0]}, end: {time_range[1]}"

            return self.generate_spikes_helper(frequencies=freq, time_ranges=time_range, rng=rng)

        # https://stackoverflow.com/questions/5148635/how-to-simulate-poisson-arrival
        start_time = time_range[0]
        end_time = time_range[1]
        duration = end_time - start_time

        assert duration > 0, f"Start time {start_time} and end time {end_time} incorrect (duration > 0 required)"

        if freq > 0:
            t_diff = -np.log(1.0 - rng.random(int(np.ceil(max(1, freq * duration))))) / freq

            t_spikes = [start_time + np.cumsum(t_diff)]

            # Is last spike after end of duration
            while t_spikes[-1][-1] <= end_time:
                t_diff = -np.log(1.0 - rng.random(int(np.ceil(freq * duration * 0.1)))) / freq
                t_spikes.append(t_spikes[-1][-1] + np.cumsum(t_diff))

            # Prune away any spikes after end
            if len(t_spikes[-1]) > 0:
                t_spikes[-1] = t_spikes[-1][t_spikes[-1] <= end_time]

            # Return spike times
            return np.concatenate(t_spikes)
        else:
            # Frequency was 0 or negative(!)
            assert not freq < 0, "Negative frequency specified."
            return np.array([])

    ############################################################################

    # This takes a list of spike trains and returns a single spike train
    # including all spikes

    @staticmethod
    def mix_spikes(spikes):

        return np.sort(np.concatenate(spikes))

    ############################################################################

    @staticmethod
    def cull_spikes(spikes, p_keep, rng):

        return spikes[rng.random(spikes.shape) < p_keep]

    ############################################################################

    # timeRange --- (start,end time) of spike train
    # freq -- frequency of spike train
    # nSpikeTrains -- number of spike trains to generate
    # Pkeep -- fraction of channel spikes to include in spike train
    # retpopUnitSpikes -- if true, returns tuple with second item population unit spikes
    #                  if false, just return spikes
    # populationUnitSpikes --- if None, new population unit spikes will be generated
    #                  (population unit Spikes are the spikes shared between correlated
    #                   spike trains)

    def make_correlated_spikes(self, freq, time_range, num_spike_trains, p_keep, rng,
                               population_unit_spikes=None, 
                               ret_pop_unit_spikes=False, jitter_dt=None):

        assert (0 <= p_keep <= 1)

        if population_unit_spikes is None:
            population_unit_spikes = self.generate_spikes(freq, time_range, rng=rng)

        if type(freq) == list:
            unique_freq = [f * (1 - p_keep) for f in freq]
        else:
            unique_freq = freq * (1 - p_keep)
        spike_trains = []

        for i in range(0, num_spike_trains):
            t_unique = self.generate_spikes(unique_freq, time_range, rng)
            t_population_unit = self.cull_spikes(population_unit_spikes, p_keep, rng)

            spike_trains.append(self.mix_spikes([t_unique, t_population_unit]))

        # if(False):
        # self.verifyCorrelation(spikeTrains=spikeTrains) # THIS STEP IS VERY VERY SLOW

        if jitter_dt is not None:
            spike_trains = self.jitter_spikes(spike_trains, jitter_dt, time_range=time_range, rng=rng)

        if ret_pop_unit_spikes:
            return spike_trains, population_unit_spikes
        else:
            return spike_trains

    ############################################################################

    def make_uncorrelated_spikes(self, freq, t_start, t_end, n_spike_trains, rng):

        spike_trains = []

        for i in range(0, n_spike_trains):
            spike_trains.append(self.generate_spikes(freq, (t_start, t_end), rng))

        return spike_trains

    ############################################################################

    # If a timeRange (start,endtime) is given then all spike times will
    # be modulo duration, so if we jitter and they go to before start time,
    # they wrap around and appear at end of the timeline

    @staticmethod
    def jitter_spikes(spike_trains, dt, rng, time_range=None):

        jittered_spikes = []

        for i in range(0, len(spike_trains)):
            spikes = spike_trains[i] + rng.normal(0, dt, spike_trains[i].shape)

            # No modulo time jittering if list of times specified
            if time_range is not None and type(time_range[0]) != list:
                start = time_range[0]
                end = time_range[1]
                spikes = np.mod(spikes - start, end - start) + start

            s = np.sort(spikes)
            # Remove any spikes that happened to go negative
            s = s[np.where(s >= 0)]
            jittered_spikes.append(s)

        return jittered_spikes

    ############################################################################

    # Plot spikes as a raster plot, for debugging and visualisation purposes

    @staticmethod
    def raster_plot(spike_times,
                    mark_spikes=None, mark_idx=None,
                    title=None, fig_file=None, fig=None):

        if fig is None:
            fig = plt.figure()
        # ax = plt.gca()

        for i, spikes in enumerate(spike_times):
            plt.vlines(spikes, i + 1.5, i + 0.5, color="black")

        plt.ylim(0.5, len(spike_times) + 0.5)

        if mark_spikes is not None and mark_idx is not None:
            for i, spikes in zip(mark_idx, mark_spikes):
                plt.vlines(spikes, i + 1.5, i + 0.5, color="red")

            plt.ylim(min(0.5, min(mark_idx) - 0.5),
                     max(max(mark_idx) + 0.5, len(spike_times)) + 0.5)

        plt.xlabel("Time")
        plt.ylabel("Inputs")

        plt.ion()
        plt.show()

        if title is not None:
            plt.title(title)

        fig.show()

        if fig_file is not None:
            plt.savefig(fig_file)

        return fig

    ############################################################################

    def read_neuron_positions(self):

        self.write_log("Reading neuron postions")

        pos_info = SnuddaLoad(self.position_file).data
        self.network_info = pos_info
        self.neuron_info = pos_info["neurons"]

        #    import pdb
        #    pdb.set_trace()

        # Make sure the position file matches the network config file
        assert (pos_info["configFile"] == self.network_config_file)

        self.num_population_units = pos_info["nPopulationUnits"]
        self.population_unit_id = pos_info["populationUnit"]

        self.neuron_id = [n["neuronID"] for n in self.neuron_info]
        self.neuron_name = [n["name"] for n in self.neuron_info]
        self.neuron_type = [n["type"] for n in self.neuron_info]
        # self.nInputs =  [n["nInputs"] for n in self.neuronInfo]

    ############################################################################

    def read_network_config_file(self):

        self.write_log("Reading config file " + str(self.network_config_file))

        with open(self.network_config_file, 'r') as f:
            self.network_config = json.load(f)

        # This also loads random seed from config file while we have it open
        if self.random_seed is None:
            if "RandomSeed" in self.network_config and "input" in self.network_config["RandomSeed"]:
                self.random_seed = self.network_config["RandomSeed"]["input"]
                self.write_log(f"Reading random see from config file: {self.random_seed}")
            else:
                # No random seed given, invent one
                self.random_seed = 1004
                self.write_log(f"No random seed provided, using: {self.random_seed}")
        else:
            self.write_log(f"Using random seed provided by command line: {self.random_seed}")

    def generate_seeds(self, num_states):

        ss = np.random.SeedSequence(self.random_seed)
        all_seeds = ss.generate_state(num_states+1)

        return all_seeds[1:]  # First seed in sequence is reserved for master

    def get_master_node_rng(self):

        ss = np.random.SeedSequence(self.random_seed)
        master_node_seed = ss.generate_state(1)
        return np.random.default_rng(master_node_seed)

    ############################################################################

    def verify_correlation(self, spike_trains, expected_corr=None, dt=0):

        # THIS FUNCTION IS VERY VERY SLOW

        corr_vec = []

        for si, s in enumerate(spike_trains):
            for s2i, s2 in enumerate(spike_trains):
                if si == s2i:
                    # No self comparison
                    continue

                corr_vec.append(self.estimate_correlation(s, s2, dt=dt))

        # print("corr = " + str(corrVec))
        self.write_log("meanCorr = " + str(np.mean(corr_vec)))

    ############################################################################

    @staticmethod
    def estimate_correlation(spikes_a, spikes_b, dt=0):

        n_spikes_a = len(spikes_a)
        corr_spikes = 0

        for t in spikes_a:
            if np.min(abs(spikes_b - t)) <= dt:
                corr_spikes += 1

        return corr_spikes / float(n_spikes_a)

    ############################################################################

    # inputDensity = f(d) where d is micrometers from soma,
    #                unit of f is synapses/micrometer

    # !!! Returns input locations only on dendrites, not on soma

    def dendrite_input_locations(self,
                                 neuron_id,
                                 rng,
                                 synapse_density="1",
                                 num_spike_trains=None):

        neuron_name = self.neuron_name[neuron_id]
        swc_file = self.network_config["Neurons"][neuron_name]["morphology"]

        if swc_file in self.neuron_cache:
            morphology = self.neuron_cache[swc_file]
        else:
            morphology = NeuronMorphology(name=neuron_id,
                                          swc_filename=swc_file,
                                          axon_stump_id_flag=self.axon_stump_id_flag)
            self.neuron_cache[swc_file] = morphology

        return morphology.dendrite_input_locations(synapse_density=synapse_density,
                                                   num_locations=num_spike_trains,
                                                   rng=rng)

    ############################################################################

    # Returns random dendrite compartments
    # Pdist = function of d (micrometers), e.g. "1", "2*d", ... etc

    def random_compartment(self, neuron_id, num_locations, rng, p_dist="1"):

        assert False, "obsolete, remove! Use dendrite_input_locations in NeuronMorphology"

        neuron_name = self.neuron_name[neuron_id]
        swc_file = self.network_config[neuron_name]["morphology"]

        if swc_file in self.neuron_cache:
            morphology = self.neuron_cache[swc_file]
        else:
            # !!! Should I use NeuronMorphology, or maybe the loadSWC in ConvertNetwork
            #                      -- maybe move loadSWC out to a separate file

            morphology = NeuronMorphology(name=neuron_id,
                                          swc_filename=swc_file,
                                          axon_stump_id_flag=self.axon_stump_id_flag)
            self.neuron_cache[swc_file] = morphology

        # morphology.dend -- 0-2: x,y,z 3: r, 4: dist to soma
        d = morphology.dend[:, 4]
        p_x = eval(p_dist)

        if type(p_x) in (int, float):
            # If Px is a constant, we need to set it for all points
            p_x *= np.ones(d.shape)

        p_comp = (p_x[morphology.dend_links[:, 0]] + p_x[morphology.dend_links[:, 1]]) / 2
        comp_len = morphology.compartment_length(comp_type="dend")

        # Multiply by length, so longer compartments are proportionally more likely to be
        # connected to
        p_comp = np.multiply(p_comp, comp_len)

        # Randomize locations. Here we then sort the locations, this is to make
        # it quicker to locate where they are (one pass in the Pxcumsum array)
        p_x_cumsum = np.cumsum(p_comp)
        x = rng.uniform(low=0, high=p_x_cumsum[-1], size=num_locations)
        x_sort = np.sort(x)

        n_comps = len(d)
        last_valid_idx = n_comps - 1
        comp_idx = 0

        sec_id = np.zeros((num_locations,), dtype=int)
        sec_x = np.zeros((num_locations,))
        comp_coords = np.zeros((num_locations, 3))

        for ix, xval in enumerate(x_sort):
            while xval > p_x_cumsum[comp_idx] and comp_idx < last_valid_idx:
                comp_idx += 1

            sec_id[ix] = morphology.dend_sec_id[comp_idx]
            sec_x[ix] = (rng.random() * (morphology.dend_sec_x[comp_idx, 1] - morphology.dend_sec_x[comp_idx, 0])
                         + morphology.dend_sec_x[comp_idx, 0])

            comp_coords[ix, :3] = (morphology.dend[morphology.dend_links[comp_idx, 0], :3]
                                   + morphology.dend[morphology.dend_links[comp_idx, 1], :3]) / 2

        # The probability of picking a compartment is dependent on a distance
        # dependent probability function, and the length of the compartment

        # Plot synapses and segments, as a verification
        if False:
            plt.hist(d[compID], label="synapses")
            plt.hist(d, label="segments")
            plt.legend(loc='upper right')
            plt.show()

        return comp_coords, sec_id, sec_x

    ############################################################################

    def setup_parallel(self):

        slurm_job_id = os.getenv("SLURM_JOBID")

        if slurm_job_id is None:
            self.slurm_id = 0
        else:
            self.slurm_id = int(slurm_job_id)

        if self.rc is not None:
            # http://davidmasad.com/blog/simulation-with-ipyparallel/
            # http://people.duke.edu/~ccc14/sta-663-2016/19C_IPyParallel.html
            self.write_log("Client IDs: " + str(self.rc.ids))
            self.d_view = self.rc[:]  # Direct view into clients
            self.lb_view = self.rc.load_balanced_view()

            if self.logfile is not None:
                log_filename = self.logfile.name
                engine_logfile = [log_filename + "-" + str(x) for x in range(0, len(self.d_view))]
            else:
                engine_logfile = [None for x in range(0, len(self.d_view))]
        else:
            self.write_log("Running in serial")
            self.d_view = None
            self.lb_view = None
            return

        with self.d_view.sync_imports():
            from snudda.input import SnuddaInput

        self.d_view.push({"input_config_file": self.input_config_file,
                          "network_config_file": self.network_config_file,
                          "position_file": self.position_file,
                          "spike_data_filename": self.spike_data_filename,
                          "is_master": False,
                          "time": self.time,
                          "random_seed": self.random_seed})

        self.write_log("Scattering engineLogFile = " + str(engine_logfile))

        self.d_view.scatter('log_filename', engine_logfile, block=True)

        self.write_log(f"nl = SnuddaInput(input_config_file='{self.input_config_file}'"
                       f",network_config_file='{self.network_config_file}'"
                       f",position_file='{self.position_file}'"
                       f",spike_data_filename='{self.spike_data_filename}'"
                       f",is_master=False "
                       f", random_seed={self.random_seed}"
                       f",time={self.time}, logfile=log_filename[0])")

        cmd_str = ("global nl; nl = SnuddaInput(input_config_file=input_config_file," 
                   "network_config_file=network_config_file,position_file=position_file," 
                   " spike_data_filename=spike_data_filename,is_master=is_master,time=time," 
                   " random_seed=random_seed, logfile=log_filename[0])")

        self.d_view.execute(cmd_str, block=True)

        self.write_log("Workers set up")

    ############################################################################

    # Function for debugging

    @staticmethod
    def dump_to_random_file(file_prefix, data_to_dump):
        import uuid
        tmp = open("save/" + file_prefix + "-file-" + str(uuid.uuid4()), 'w')
        tmp.write(str(data_to_dump))
        tmp.close()

    ############################################################################

    def check_sorted(self):

        # Just a double check that the spikes are not jumbled

        for neuron_id in self.neuron_input:
            for input_type in self.neuron_input[neuron_id]:
                if input_type == "VirtualNeuron":
                    s = self.neuron_input[neuron_id][input_type]["spikes"]
                    assert (np.diff(s) >= 0).all(), \
                        str(neuron_id) + " " + input_type + ": Spikes must be in order"
                else:
                    for spikes in self.neuron_input[neuron_id][input_type]["spikes"]:
                        assert len(spikes) == 0 or spikes[0] >= 0
                        assert (np.diff(spikes) >= 0).all(), \
                            str(neuron_id) + " " + input_type + ": Spikes must be in order"

    ############################################################################

    def plot_spikes(self, neuron_id=None):

        self.write_log("Plotting spikes for neuronID: " + str(neuron_id))

        if neuron_id is None:
            neuron_id = self.neuron_input

        spike_times = []

        for nID in neuron_id:
            for inputType in self.neuron_input[nID]:
                for spikes in self.neuron_input[nID][inputType]["spikes"]:
                    spike_times.append(spikes)

        self.raster_plot(spike_times)

    ############################################################################

    def read_hdf5_info(self, hdf5_file):
        self.write_log("Loading HDF5-file: " + hdf5_file)

        try:
            with h5py.File(hdf5_file, 'r') as f:
                self.network_config_file = SnuddaLoad.to_str(f["meta"]["configFile"][()])
                self.position_file = SnuddaLoad.to_str(f["meta"]["positionFile"][()])
                self.network_slurm_id = int(f["meta/SlurmID"][()])

                self.axon_stump_id_flag = f["meta/axonStumpIDFlag"][()]
        except Exception as e:
            self.write_log("Error in readHDF5info: " + str(e))
            self.write_log("Opening: " + hdf5_file)

            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)
            os.sys.exit(-1)

    ############################################################################

    def find_latest_file(self):

        files = glob('save/network-connect-voxel-pruned-synapse-file-*.hdf5')

        mod_time = [os.path.getmtime(f) for f in files]
        idx = np.argsort(mod_time)

        self.write_log("Using the newest file: " + files[idx[-1]])

        return files[idx[-1]]

    ############################################################################

    def make_input_helper_parallel(self, args):

        try:

            neuron_id, input_type, freq, start, end, synapse_density, num_spike_trains, p_keep, \
                population_unit_spikes, jitter_dt, population_unit_id, conductance, correlation, mod_file, \
                parameter_file, parameter_list, random_seed = args

            return self.make_input_helper_serial(neuron_id=neuron_id,
                                                 input_type=input_type,
                                                 freq=freq,
                                                 t_start=start,
                                                 t_end=end,
                                                 synapse_density=synapse_density,
                                                 num_spike_trains=num_spike_trains,
                                                 p_keep=p_keep,
                                                 population_unit_spikes=population_unit_spikes,
                                                 jitter_dt=jitter_dt,
                                                 population_unit_id=population_unit_id,
                                                 conductance=conductance,
                                                 correlation=correlation,
                                                 mod_file=mod_file,
                                                 parameter_file=parameter_file,
                                                 parameter_list=parameter_list,
                                                 random_seed=random_seed)

        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)
            import pdb
            pdb.set_trace()

    ############################################################################

    # Normally specify synapseDensity which then sets number of inputs
    # ie leave nSpikeTrains as None. If nSpikeTrains is set, that will then
    # scale synapseDensity to get the requested number of inputs (approximately)

    # For virtual neurons nSpikeTrains must be set, as it defines their activity

    def make_input_helper_serial(self,
                                 neuron_id,
                                 input_type,
                                 freq,
                                 t_start,
                                 t_end,
                                 synapse_density,
                                 num_spike_trains,
                                 p_keep,
                                 population_unit_spikes,
                                 jitter_dt,
                                 population_unit_id,
                                 conductance,
                                 correlation,
                                 mod_file,
                                 parameter_file,
                                 parameter_list,
                                 random_seed):

        # First, find out how many inputs and where, based on morphology and
        # synapse density

        time_range = (t_start, t_end)

        rng = np.random.default_rng(random_seed)

        if input_type.lower() == "VirtualNeuron".lower():
            # This specifies activity of a virtual neuron
            conductance = None

            assert num_spike_trains is None or num_spike_trains == 1, \
                "Virtual neuron " + self.neuron_name[neuron_id] \
                + " should have only one spike train, fix nSpikeTrains in config"

            # Virtual neurons input handled through touch detection
            input_loc = None

            spikes = self.make_correlated_spikes(freq=freq,
                                                 time_range=time_range,
                                                 num_spike_trains=1,
                                                 p_keep=p_keep,
                                                 population_unit_spikes=population_unit_spikes,
                                                 jitter_dt=jitter_dt,
                                                 rng=rng)
            num_inputs = 1
        else:

            # (x,y,z), secID, secX, dist_to_soma
            input_loc = self.dendrite_input_locations(neuron_id=neuron_id,
                                                      synapse_density=synapse_density,
                                                      num_spike_trains=num_spike_trains,
                                                      rng=rng)

            num_inputs = input_loc[0].shape[0]
            print("Generating " + str(num_inputs) + " inputs for " + self.neuron_name[neuron_id])

            # OBS, nInputs might differ slightly from nSpikeTrains if that is given
            spikes = self.make_correlated_spikes(freq=freq,
                                                 time_range=time_range,
                                                 num_spike_trains=num_inputs,
                                                 p_keep=p_keep,
                                                 population_unit_spikes=population_unit_spikes,
                                                 jitter_dt=jitter_dt,
                                                 rng=rng)

        # We need to pick which parameter set to use for the input also
        parameter_id = rng.integers(1e6, size=num_inputs)

        # We need to keep track of the neuronID, since it will all be jumbled
        # when doing asynchronous parallelisation
        return (neuron_id, input_type, spikes, input_loc, synapse_density, freq,
                jitter_dt, population_unit_id, conductance, correlation,
                time_range,
                mod_file, parameter_file, parameter_list, parameter_id)

    ############################################################################

    def write_log(self, text, flush=True):  # Change flush to False in future, debug
        if self.logfile is not None:
            self.logfile.write(text + "\n")
            print(text)
            if flush:
                self.logfile.flush()
        else:
            if self.verbose:
                print(text)

    ############################################################################


if __name__ == "__main__":
    print("Please do not call this file directly, use snudda command line")
    os.sys.exit(-1)
