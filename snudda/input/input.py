# This code writes the input spikes for the NEURON simulation --
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

# TODO: Randomise conductance for the inputs and store it, use it later when adding external synapses in simulate.py

import json
import os
import sys
from collections import OrderedDict

import h5py
import numexpr
import numpy as np

from snudda.utils.snudda_path import get_snudda_data
from snudda.input.time_varying_input import TimeVaryingInput
from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.utils.load import SnuddaLoad
from snudda.utils.snudda_path import snudda_parse_path

nl = None


# When specifying vectors for start and end time, they should normally not overlap
# if we want to allow that, set time_interval_overlap_warning = False


class SnuddaInput(object):
    """ Generates input for the simulation. """

    def __init__(self,
                 network_path=None,
                 snudda_data=None,
                 input_config_file=None,
                 spike_data_filename=None,
                 hdf5_network_file=None,
                 time=10.0,
                 is_master=True,
                 h5libver="latest",
                 rc=None,
                 random_seed=None,
                 time_interval_overlap_warning=True,
                 logfile=None,
                 verbose=False,
                 use_meta_input=True):

        """
        Constructor.

        Args:
            network_path (str): Path to network directory
            snudda_data (str): Path to Snudda Data
            input_config_file (str): Path to input config file, default input.json in network_path
            spike_data_filename (str): Path to output file, default input-spikes.hdf5
            hdf5_network_file (str): Path to network file, default network-synapses.hdf5
            time (float): Duration of simulation to generate input for, default 10 seconds
            is_master (bool): "master" or "worker"
            h5libver (str): Version of HDF5 library to use, default "latest"
            rc: ipyparallel remote client
            random_seed (int): Random seed for input generation
            time_interval_overlap_warning (bool): Warn if input intervals specified overlap
            logfile (str): Log file
            verbose (bool): Print logging
        """

        if type(logfile) == str:
            self.logfile = open(logfile, "w")
        else:
            self.logfile = logfile

        self.verbose = verbose
        self.rc = rc

        if network_path:
            self.network_path = network_path
        elif hdf5_network_file:
            self.network_path = os.path.dirname(hdf5_network_file)
        elif input_config_file:
            self.network_path = os.path.dirname(input_config_file)
        else:
            self.network_path = None

        self.snudda_data = get_snudda_data(snudda_data=snudda_data,
                                           network_path=self.network_path)

        if input_config_file:
            self.input_config_file = input_config_file
        elif self.network_path:
            self.input_config_file = os.path.join(self.network_path, "input.json")
        else:
            self.input_config_file = None

        if spike_data_filename:
            self.spike_data_filename = spike_data_filename
        elif self.network_path:
            self.spike_data_filename = os.path.join(self.network_path, "input-spikes.hdf5")
        else:
            self.spike_data_filename = None

        if hdf5_network_file:
            self.hdf5_network_file = hdf5_network_file
        elif self.network_path:
            self.hdf5_network_file = os.path.join(self.network_path, "network-synapses.hdf5")
        else:
            self.hdf5_network_file = None

        self.time_interval_overlap_warning = time_interval_overlap_warning
        self.input_info = None
        self.population_unit_spikes = None
        self.all_population_units = None  # List of all population units in simulation

        self.num_population_units = None
        self.population_unit_id = None
        self.neuron_name = None
        self.neuron_id = None
        self.neuron_type = None
        self.d_view = None
        self.network_config = None
        self.neuron_input = None
        self.slurm_id = None

        self.snudda_load = None
        self.network_data = None
        self.neuron_info = None

        self.network_config_file = None
        self.position_file = None

        self.axon_stump_id_flag = None
        self.network_slurm_id = None
        self.population_unit_id = []

        self.use_meta_input = use_meta_input

        self.neuron_id = []
        self.neuron_name = []
        self.neuron_type = []

        if self.hdf5_network_file:
            self.load_network(self.hdf5_network_file)
        else:
            print("No network file specified, use load_network to load network info")

        if time:
            self.time = time  # How long time to generate inputs for
        else:
            self.time = 10

        self.write_log(f"Time = {time}")

        self.random_seed = random_seed

        self.h5libver = h5libver
        self.write_log(f"Using hdf5 version {h5libver}")

        self.neuron_cache = dict([])

        self.is_master = is_master

    def load_network(self, hdf5_network_file=None):

        if hdf5_network_file is None:
            hdf5_network_file = self.hdf5_network_file

        self.snudda_load = SnuddaLoad(hdf5_network_file)
        self.network_data = self.snudda_load.data
        self.neuron_info = self.network_data["neurons"]

        self.network_config_file = self.network_data["configFile"]
        self.position_file = self.network_data["positionFile"]

        self.axon_stump_id_flag = self.network_data["axonStumpIDFlag"]
        self.network_slurm_id = self.network_data["SlurmID"]
        self.population_unit_id = self.network_data["populationUnit"]

        self.neuron_id = [n["neuronID"] for n in self.network_data["neurons"]]
        self.neuron_name = [n["name"] for n in self.network_data["neurons"]]
        self.neuron_type = [n["type"] for n in self.network_data["neurons"]]

    def generate(self):

        """ Generates input for network. """

        # Read in the input configuration information from JSON file
        self.read_input_config_file()

        # Read the network config file -- This also reads random seed
        self.read_network_config_file()

        # Only the master node should start the work
        if self.is_master:
            # Initialises lbView and dView (load balance, and direct view)
            self.setup_parallel()

            # Make the "master input" for each channel
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

        """ Writes input spikes to HDF5 file. """

        self.write_log(f"Writing spikes to {self.spike_data_filename}", force_print=True)

        out_file = h5py.File(self.spike_data_filename, 'w', libver=self.h5libver)
        out_file.create_dataset("config", data=json.dumps(self.input_info, indent=4))
        input_group = out_file.create_group("input")

        for neuron_id in self.neuron_input:

            nid_group = input_group.create_group(str(neuron_id))

            neuron_type = self.neuron_type[neuron_id]

            for input_type in self.neuron_input[neuron_id]:

                if input_type[0] == '!':
                    self.write_log(f"Disabling input {input_type} for neuron {neuron_id} "
                                   f" (input_type was commented with ! before name)")
                    continue

                if input_type.lower() != "VirtualNeuron".lower():
                    it_group = nid_group.create_group(input_type)

                    neuron_in = self.neuron_input[neuron_id][input_type]
                    spike_mat, num_spikes = self.create_spike_matrix(neuron_in["spikes"])

                    it_group.create_dataset("spikes", data=spike_mat, compression="gzip", dtype=np.float32)
                    it_group.create_dataset("nSpikes", data=num_spikes, dtype=np.int32)

                    it_group.create_dataset("sectionID", data=neuron_in["location"][1].astype(int),
                                            compression="gzip", dtype=np.int16)
                    it_group.create_dataset("sectionX", data=neuron_in["location"][2],
                                            compression="gzip", dtype=np.float16)
                    it_group.create_dataset("distanceToSoma", data=neuron_in["location"][3],
                                            compression="gzip", dtype=np.float16)

                    it_group.create_dataset("freq", data=neuron_in["freq"])
                    it_group.create_dataset("correlation", data=neuron_in["correlation"])

                    if "jitter" in neuron_in and neuron_in["jitter"]:
                        it_group.create_dataset("jitter", data=neuron_in["jitter"])

                    if "synapseDensity" in neuron_in and neuron_in["synapseDensity"]:
                        it_group.create_dataset("synapseDensity", data=neuron_in["synapseDensity"])

                    it_group.create_dataset("start", data=neuron_in["start"])
                    it_group.create_dataset("end", data=neuron_in["end"])
                    it_group.create_dataset("conductance", data=neuron_in["conductance"])

                    population_unit_id = int(neuron_in["populationUnitID"])
                    it_group.create_dataset("populationUnitID", data=population_unit_id)

                    # TODO: What to do with population_unit_spikes, should we have mandatory jittering for them?

                    # population_unit_id = 0 means not population unit membership, so no population spikes available
                    if neuron_type in self.population_unit_spikes and population_unit_id > 0 \
                            and input_type in self.population_unit_spikes[neuron_type]:
                        chan_spikes = self.population_unit_spikes[neuron_type][input_type][population_unit_id]
                    else:
                        chan_spikes = np.array([])

                    it_group.create_dataset("populationUnitSpikes", data=chan_spikes, compression="gzip",
                                            dtype=np.float32)

                    it_group.create_dataset("generator", data=neuron_in["generator"])

                    it_group.create_dataset("modFile", data=neuron_in["modFile"])
                    if neuron_in["parameterFile"]:
                        it_group.create_dataset("parameterFile", data=neuron_in["parameterFile"])
                    # We need to convert this to string to be able to save it
                    it_group.create_dataset("parameterList", data=json.dumps(neuron_in["parameterList"]))
                    it_group.create_dataset("parameterID", data=neuron_in["parameterID"], dtype=np.int32)

                else:

                    # Input is activity of a virtual neuron
                    a_group = nid_group.create_group("activity")
                    spikes = self.neuron_input[neuron_id][input_type]["spikes"]

                    a_group.create_dataset("spikes", data=spikes, compression="gzip")
                    generator = self.neuron_input[neuron_id][input_type]["generator"]
                    a_group.create_dataset("generator", data=generator)

        out_file.close()

    ############################################################################

    @staticmethod
    def create_spike_matrix(spikes):

        """ Creates a spike matrix from a list of spikes. """

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

        """ Read input configuration from JSON file. """

        self.write_log(f"Loading input configuration from {self.input_config_file}")

        with open(snudda_parse_path(self.input_config_file, self.snudda_data), 'rt') as f:
            self.input_info = json.load(f, object_pairs_hook=OrderedDict)

        for neuron_type in self.input_info:
            for input_type in self.input_info[neuron_type]:
                if "parameterFile" in self.input_info[neuron_type][input_type]:
                    # Allow user to use $DATA to refer to snudda data directory
                    par_file = snudda_parse_path(self.input_info[neuron_type][input_type]["parameterFile"],
                                                 self.snudda_data)

                    with open(par_file, 'r') as f:
                        par_data_dict = json.load(f, object_pairs_hook=OrderedDict)

                    # Read in parameters into a list
                    par_data = []
                    for pd in par_data_dict:

                        if "parameterList" in self.input_info[neuron_type][input_type]:
                            for par_key, par_d in self.input_info[neuron_type][input_type]["parameterList"].items():
                                print(f"Overriding {par_key} with value {par_d} for {neuron_type}:{input_type}")
                                par_data_dict[pd]["synapse"][par_key] = par_d

                        par_data.append(par_data_dict[pd])
                elif "parameterList" in self.input_info[neuron_type][input_type]:

                    par_data = [{"synapse": self.input_info[neuron_type][input_type]["parameterList"]}]
                else:
                    par_data = None

                self.input_info[neuron_type][input_type]["parameterList"] = par_data

    ############################################################################

    # Each synaptic input will contain a fraction of population unit spikes, which are
    # taken from a stream of spikes unique to that particular population unit
    # This function generates these correlated spikes

    def make_population_unit_spike_trains(self, rng):

        """
        Generate population unit spike trains.
        Each synaptic input will contain a fraction of population unit spikes, which are
        taken from a stream of spikes unique to that particular population unit
        This function generates these correlated spikes
        """

        self.write_log("Running makePopulationUnitSpikeTrains")

        self.population_unit_spikes = dict([])

        for cell_type in self.input_info:

            self.population_unit_spikes[cell_type] = dict([])

            for input_type in self.input_info[cell_type]:

                if "start" in self.input_info[cell_type][input_type]:
                    start_time = np.array(self.input_info[cell_type][input_type]["start"])
                else:
                    start_time = 0

                if "end" in self.input_info[cell_type][input_type]:
                    end_time = np.array(self.input_info[cell_type][input_type]["end"])
                else:
                    end_time = self.time

                if "populationUnitID" in self.input_info[cell_type][input_type]:
                    pop_unit_list = self.input_info[cell_type][input_type]["populationUnitID"]

                    if type(pop_unit_list) != list:
                        pop_unit_list = [pop_unit_list]
                else:
                    pop_unit_list = self.all_population_units

                # Handle Poisson input
                if self.input_info[cell_type][input_type]["generator"] == "poisson":

                    freq = self.input_info[cell_type][input_type]["frequency"]
                    self.population_unit_spikes[cell_type][input_type] = dict([])

                    for idx_pop_unit in pop_unit_list:
                        self.population_unit_spikes[cell_type][input_type][idx_pop_unit] = \
                            self.generate_poisson_spikes(freq=freq, time_range=(start_time, end_time), rng=rng)

                # Handle frequency function
                elif self.input_info[cell_type][input_type]["generator"] == "frequency_function":

                    frequency_function = self.input_info[cell_type][input_type]["frequency"]
                    self.population_unit_spikes[cell_type][input_type] = dict([])

                    for idx_pop_unit in pop_unit_list:
                        try:
                            self.population_unit_spikes[cell_type][input_type][idx_pop_unit] = \
                                self.generate_spikes_function(frequency_function=frequency_function,
                                                              time_range=(start_time, end_time),
                                                              rng=rng)
                        except:
                            import traceback
                            print(traceback.format_exc())
                            import pdb
                            pdb.set_trace()
                else:
                    assert False, f"Unknown input generator {self.input_info[cell_type][input_type]['generator']} " \
                                  f"for cell_type {cell_type}, input_type {input_type}"

                if "setMotherSpikes" in self.input_info[cell_type][input_type]:
                    self.write_log(f"Warning, overwriting mother spikes for {cell_type} {input_type} with user defined spikes")

                    for idx_pop_unit in pop_unit_list:
                        # User defined mother spikes
                        self.population_unit_spikes[cell_type][input_type][idx_pop_unit] = \
                            np.array(self.input_info[cell_type][input_type]["setMotherSpikes"])

                if "addMotherSpikes" in self.input_info[cell_type][input_type]:
                    self.write_log(f"Adding user defined extra spikes to mother process for {cell_type} {input_type} -- but not for population unit 0")

                    for idx_pop_unit in pop_unit_list:
                        self.population_unit_spikes[cell_type][input_type][idx_pop_unit] = \
                            np.sort(np.concatenate((self.population_unit_spikes[cell_type][input_type][idx_pop_unit],
                                                   np.array(self.input_info[cell_type][input_type]["addMotherSpikes"]))))

        return self.population_unit_spikes

    ############################################################################

    def make_neuron_input_parallel(self):

        """ Generate input, able to run in parallel if rc (Remote Client) has been provided at initialisation."""

        self.write_log("Running make_neuron_input_parallel")

        self.neuron_input = dict([])

        neuron_id_list = []
        input_type_list = []
        freq_list = []
        start_list = []
        end_list = []
        synapse_density_list = []
        num_inputs_list = []
        population_unit_spikes_list = []
        jitter_dt_list = []
        population_unit_id_list = []
        conductance_list = []
        correlation_list = []

        mod_file_list = []
        parameter_file_list = []
        parameter_list_list = []
        cluster_size_list = []
        cluster_spread_list = []

        generator_list = []
        population_unit_fraction_list = []

        dendrite_location_override_list = []
        if self.use_meta_input:
            self.write_log("Input from meta.json will be used")
        else:
            self.write_log("Input from meta.json will NOT be used")

        for (neuron_id, neuron_name, neuron_type, population_unit_id) \
                in zip(self.neuron_id, self.neuron_name, self.neuron_type, self.population_unit_id):

            self.neuron_input[neuron_id] = dict([])

            # The input can be specified using neuron_id, neuron_name or neuron_type
            if str(neuron_id) in self.input_info:
                input_info = self.input_info[str(neuron_id)].copy()
            elif neuron_name in self.input_info:
                input_info = self.input_info[neuron_name].copy()
            elif neuron_type in self.input_info:
                input_info = self.input_info[neuron_type].copy()
            else:
                input_info = dict()

            # if a number --> use a specific neuron with that given ID
            # if dSPN --> use neuron_type dSPN
            # if dSPN_3 --> use specific neuron morphology corresponding to dSPN_3

            # Also see if we have additional input specified in the meta.json file for the neuron?

            # Add baseline activity:
            #  1. From neuron_id derive the parameter_id and morphology_id
            #  2. Using parameter_id, morphology_id check if the meta.json has any additional input specified
            #  3. Add the input to input_info

            parameter_key = self.network_data["neurons"][neuron_id]["parameterKey"]
            morphology_key = self.network_data["neurons"][neuron_id]["morphologyKey"]
            neuron_path = snudda_parse_path(self.network_data["neurons"][neuron_id]["neuronPath"], self.snudda_data)
            meta_path = os.path.join(neuron_path, "meta.json")

            if self.use_meta_input and os.path.exists(meta_path):
                with open(meta_path, "r") as f:
                    meta_data = json.load(f)

                if parameter_key in meta_data and morphology_key in meta_data[parameter_key] \
                        and "input" in meta_data[parameter_key][morphology_key]:

                    for inp_name, inp_data in meta_data[parameter_key][morphology_key]["input"].items():
                        if inp_name in input_info:

                            self.write_log(f"!!! Warning, combining definition of {inp_name} input for neuron "
                                           f"{self.network_data['neurons'][neuron_id]['name']} {neuron_id} "
                                           f"(meta modified by input_config)",
                                           force_print=True)

                            old_info = input_info[inp_name]
                            new_info = inp_data.copy()

                            for key, data in old_info.items():
                                new_info[key] = data

                            input_info[inp_name] = new_info.copy()

                        else:
                            input_info[inp_name] = inp_data.copy()

            if len(input_info) == 0:
                self.write_log(f"!!! Warning, no synaptic input for neuron ID {neuron_id}, "
                               f"name {neuron_name} or type {neuron_type}")

            for input_type in input_info:

                if input_type[0] == '!':
                    self.write_log(f"Disabling input {input_type} for neuron {neuron_name} "
                                   f" (input_type was commented with ! before name)")
                    continue

                input_inf = input_info[input_type].copy()

                if "populationUnitID" in input_inf:
                    pop_unit_id = input_inf["populationUnitID"]

                    if type(pop_unit_id) in [list, np.ndarray] and population_unit_id not in pop_unit_id:
                        # We have a list of functional channels, but this neuron
                        # does not belong to a functional channel in that list
                        continue
                    elif population_unit_id != int(pop_unit_id):
                        # We have a single functional channel, but this neuron is not
                        # in that functional channel
                        continue
                    else:
                        pop_unit_id = int(pop_unit_id)
                else:
                    pop_unit_id = None

                self.neuron_input[neuron_id][input_type] = dict([])

                if input_inf["generator"] == "csv":
                    csv_file = snudda_parse_path(input_inf["csvFile"] % neuron_id, self.snudda_data)

                    self.neuron_input[neuron_id][input_type]["spikes"] \
                        = np.genfromtxt(csv_file, delimiter=',')
                    self.neuron_input[neuron_id][input_type]["generator"] = "csv"

                    # Done for CSV input
                    continue

                # These parameters are shared between "poisson" and "frequency_function"

                neuron_id_list.append(neuron_id)
                input_type_list.append(input_type)

                if "jitter" in input_inf:
                    jitter_dt_list.append(input_inf["jitter"])
                else:
                    jitter_dt_list.append(None)

                if "start" in input_inf:
                    start_list.append(np.array(input_inf["start"]))
                else:
                    start_list.append(0.0)  # Default start at beginning

                if "end" in input_inf:
                    end_list.append(np.array(input_inf["end"]))
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

                    # If your code fails here, it might be that you are trying to override the background input
                    # frequency, but have the incorrect name of that input (check the meta.json file)
                    cond = input_inf["conductance"]

                    if "nInputs" in input_inf:

                        # TODO: We need to read this from meta.json
                        dir_name = snudda_parse_path(os.path.basename(neuron_path), self.snudda_data)

                        # If a dictionary, then extract the info for the relevant neuron
                        # Priority order is:
                        # 1. Morphology key, 2: neuron directory name,
                        # 3: Neuron name (note this can change if additional neurons are added to neuron type dir)
                        # 4: Neuron type

                        if type(input_inf["nInputs"]) == OrderedDict:
                            if morphology_key in input_inf["nInputs"]:
                                n_inp = input_inf["nInputs"][morphology_key]
                            elif dir_name in input_inf["nInputs"]:
                                n_inp = input_inf["nInputs"][dir_name]
                            elif neuron_name in input_inf["nInputs"]:
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

                population_unit_id_list.append(population_unit_id)
                conductance_list.append(cond)

                if "populationUnitCorrelation" in input_inf:
                    correlation_list.append(input_inf["populationUnitCorrelation"])
                else:
                    correlation_list.append(0)

                if "populationUnitCorrelationFraction" in input_inf:
                    population_unit_fraction_list.append(np.array(input_inf["populationUnitCorrelationFraction"]))
                else:
                    population_unit_fraction_list.append(1)

                if (neuron_type in self.population_unit_spikes
                        and input_type in self.population_unit_spikes[neuron_type]
                        and population_unit_id in self.population_unit_spikes[neuron_type][input_type]):

                    # TODO: Currently only correlated within a neuron type for a given population unit
                    #       should the spikes be shared between all neuron types in that population unit?
                    c_spikes = self.population_unit_spikes[neuron_type][input_type][population_unit_id]
                    population_unit_spikes_list.append(c_spikes)
                else:
                    # self.write_log(f"No population spikes specified for neuron type {neuron_type}")
                    population_unit_spikes_list.append(None)

                mod_file_list.append(mod_file)
                parameter_file_list.append(parameter_file)
                parameter_list_list.append(parameter_list)

                if "clusterSize" in input_inf:
                    cluster_size = input_inf["clusterSize"]
                else:
                    cluster_size = None

                if "clusterSpread" in input_inf:
                    cluster_spread = input_inf["clusterSpread"]
                else:
                    cluster_spread = 20e-6

                cluster_size_list.append(cluster_size)
                cluster_spread_list.append(cluster_spread)

                if "dendriteLocation" in input_info:
                    assert "morphologyKey" in input_info, \
                        f"If you specify dendriteLocation you must also specify morphologyKey"

                    assert morphology_key == self.network_data["neurons"][neuron_id]["morphologyKey"], \
                        f"Neuron {neuron_id} has morphology_key " \
                        f"{self.network_data['neurons'][neuron_id]['morphologyKey']}" \
                        f"which does not match what is specified in input JSON file: {morphology_key}"

                    dend_location = input_info["dendriteLocation"]
                else:
                    dend_location = None

                dendrite_location_override_list.append(dend_location)

                if input_inf["generator"] == "poisson":
                    freq_list.append(input_inf["frequency"])
                    generator_list.append("poisson")

                elif input_inf["generator"] == "frequency_function":
                    freq_list.append(input_inf["frequency"])
                    generator_list.append("frequency_function")

                else:
                    self.write_log(f"Unknown input generator: {input_inf['generator']} for {neuron_id}", is_error=True)
                    assert False, f"Unknown input generator {input_inf['generator']}"

        seed_list = self.generate_seeds(num_states=len(neuron_id_list))

        amr = None

        assert len(neuron_id_list) == len(input_type_list) == len(freq_list)\
            == len(start_list) == len(end_list) == len(synapse_density_list) == len(num_inputs_list)\
            == len(num_inputs_list) == len(population_unit_spikes_list) == len(jitter_dt_list)\
            == len(population_unit_id_list) == len(conductance_list) == len(correlation_list)\
            == len(mod_file_list) == len(parameter_file_list) == len(parameter_list_list)\
            == len(seed_list) == len(cluster_size_list) == len(cluster_spread_list)\
            == len(dendrite_location_override_list) == len(generator_list) == len(population_unit_fraction_list),\
            "Internal error, input lists length missmatch"

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
                                  population_unit_spikes_list,
                                  jitter_dt_list,
                                  population_unit_id_list,
                                  conductance_list,
                                  correlation_list,
                                  mod_file_list,
                                  parameter_file_list,
                                  parameter_list_list,
                                  seed_list,
                                  cluster_size_list,
                                  cluster_spread_list,
                                  dendrite_location_override_list,
                                  generator_list,
                                  population_unit_fraction_list))

            self.d_view.scatter("input_list", input_list, block=True)
            cmd_str = "inpt = list(map(nl.make_input_helper_parallel,input_list))"

            self.write_log("Calling workers to generate input in parallel")
            self.d_view.execute(cmd_str, block=True)
            self.d_view.execute("nl.write_log('Execution done on workers')")

            self.write_log("Execution done")

            # On this line it stalls... WHY?
            # inpt = self.d_view["inpt"]
            amr = self.d_view.gather("inpt", block=True)
            self.write_log("Results received")

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
                      population_unit_spikes_list,
                      jitter_dt_list,
                      population_unit_id_list,
                      conductance_list,
                      correlation_list,
                      mod_file_list,
                      parameter_file_list,
                      parameter_list_list,
                      seed_list,
                      cluster_size_list,
                      cluster_spread_list,
                      dendrite_location_override_list,
                      generator_list,
                      population_unit_fraction_list)

        # Gather the spikes that were generated in parallel
        for neuron_id, input_type, spikes, loc, synapse_density, frq, \
            jdt, p_uid, cond, corr, timeRange, mod_file, param_file, param_list, param_id in amr:

            self.write_log(f"Gathering {neuron_id} - {input_type}")
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
            self.neuron_input[neuron_id][input_type]["parameterFile"] = param_file
            self.neuron_input[neuron_id][input_type]["parameterList"] = param_list
            self.neuron_input[neuron_id][input_type]["parameterID"] = param_id

        return self.neuron_input

    ############################################################################

    def generate_spikes_helper(self, frequency, time_range, rng, input_generator=None):

        if input_generator == "poisson":
            spikes = self.generate_poisson_spikes(freq=frequency, time_range=time_range, rng=rng)
        elif input_generator == "frequency_function":
            spikes = self.generate_spikes_function(frequency_function=frequency,
                                                   time_range=time_range, rng=rng)
        else:
            assert False, f"Unknown input_generator {input_generator}"

        return spikes

    def generate_poisson_spikes_helper(self, frequencies, time_ranges, rng):

        """
        Generates spike trains with given frequencies within time_ranges, using rng stream.

        Args:
             frequencies (list): List of frequencies
             time_ranges (list): List of tuples with start and end time for each frequency range
             rng: Numpy random stream
        """

        t_spikes = []

        for f, t_start, t_end in zip(frequencies, time_ranges[0], time_ranges[1]):
            t_spikes.append(self.generate_poisson_spikes(f, (t_start, t_end), rng=rng))

        # Double check correct dimension
        return np.sort(np.concatenate(t_spikes))

    def generate_poisson_spikes(self, freq, time_range, rng):
        # This generates poisson spikes with frequency freq, for a given time range

        assert np.size(freq) == np.size(time_range[0]) or np.size(freq) == 1

        if np.size(time_range[0]) > 1:

            if np.size(freq) == 1:
                freq = np.full(np.size(time_range[0]), freq)

            assert len(time_range[0]) == len(time_range[1]) == len(freq), \
                (f"Frequency, start and end time vectors need to be of same length."
                 f"\nfreq: {freq}\nstart: {time_range[0]}\nend:{time_range[1]}")

            if self.time_interval_overlap_warning:
                assert (np.array(time_range[0][1:]) - np.array(time_range[1][0:-1]) >= 0).all(), \
                    f"Time range should not overlap: start: {time_range[0]}, end: {time_range[1]}"

            return self.generate_poisson_spikes_helper(frequencies=freq, time_ranges=time_range, rng=rng)

        # https://stackoverflow.com/questions/5148635/how-to-simulate-poisson-arrival
        start_time = time_range[0]
        end_time = time_range[1]
        duration = end_time - start_time

        assert duration > 0, f"Start time = {start_time} and end time = {end_time} incorrect (duration > 0 required)"

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

    def generate_spikes_function_helper(self, frequencies, time_ranges, rng, dt, p_keep=1):

        """
        Generates spike trains with given frequencies within time_ranges, using rng stream.

        Args:
             frequencies (list): List of frequencies
             time_ranges (list): List of tuples with start and end time for each frequency range
             rng: Numpy random stream
             dt: timestep
        """

        if np.size(frequencies) == np.size(time_ranges[0]):
            frequency_list = frequencies
        else:
            frequency_list = np.full(np.size(time_ranges[0]), frequencies)

        if np.size(p_keep) == np.size(time_ranges[0]):
            p_keep_list = p_keep
        else:
            p_keep_list = np.full(np.size(time_ranges[0]), p_keep)

        t_spikes = []

        for freq, t_start, t_end, p_k in zip(frequency_list, time_ranges[0], time_ranges[1], p_keep_list):
            t_spikes.append(self.generate_spikes_function(freq, (t_start, t_end), rng=rng, dt=dt, p_keep=p_k))

        try:
            spikes = np.sort(np.concatenate(t_spikes))
        except:
            import traceback
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()

        # Double check correct dimension
        return spikes

    def generate_spikes_function(self, frequency_function, time_range, rng, dt=1e-4, p_keep=1):

        # TODO: Replace this with the code in time_varying_input.py

        """
        Generates frequency based on frequency_function.

        Args
            frequency_function: vector based python function taking t as argument, returning momentary frequency
                                if it is not a python then numexpr.evaluate is run on it (with t as argument)
                                OBS: t passed to the function is 0 at the stimultion start time, e.g. for a stimulus
                                that starts at time 4s seconds, f(t=0) is calculated, and at the end 5s f(t=1) is calculated.
            time_range: Interval of time to generate spikes for
            rng: Numpy rng object
            dt: timestep
        """

        if np.size(time_range[0]) > 1:
            return self.generate_spikes_function_helper(frequencies=frequency_function,
                                                        time_ranges=time_range,
                                                        rng=rng, dt=dt, p_keep=p_keep)

        assert 0 <= p_keep <= 1, \
            f"Error: p_keep = {p_keep}, valid range 0-1. If p_keep is a list, " \
            f"then time_ranges must be two lists, ie. (start_times, end_times)"

        if callable(frequency_function):
            func = lambda t: frequency_function(t) * p_keep
        else:
            try:
                func_str = f"{frequency_function}*{p_keep}"
                func = lambda t: numexpr.evaluate(func_str)
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

        # TODO: Utilise the n_spike trains better
        spikes = TimeVaryingInput.generate_spikes(frequency_function=func,
                                                  start_time=time_range[0], end_time=time_range[1],
                                                  n_spike_trains=1, rng=rng)[0].T

        return spikes

    ############################################################################

    # This takes a list of spike trains and returns a single spike train
    # including all spikes

    @staticmethod
    def mix_spikes(spikes):

        """ Mixes spikes in list of spike trains into one sorted spike train. """

        return np.sort(np.concatenate(spikes))

    @staticmethod
    def mix_fraction_of_spikes_OLD(spikes_a, spikes_b, fraction_a, fraction_b, rng):

        """ Picks fraction_a of spikes_a and fraction_b of spikes_b and returns sorted spike train

        Args:
            spikes_a (np.array) : Spike train A
            spikes_b (np.array) : Spike train B
            fraction_a (float) : Fraction of spikes in train A picked, e.g 0.4 means 40% of spikes are picked
            fraction_b (float) : Fraction of spikes in train B picked
            rng : Numpy rng object
        """

        len_a = np.size(spikes_a) * fraction_a
        len_b = np.size(spikes_b) * fraction_b

        len_a_rand = int(np.floor(len_a) + (len_a % 1 > rng.uniform()))
        len_b_rand = int(np.floor(len_b) + (len_b % 1 > rng.uniform()))

        idx_a = rng.choice(np.size(spikes_a), size=len_a_rand, replace=False)
        idx_b = rng.choice(np.size(spikes_b), size=len_b_rand, replace=False)

        return np.sort(np.concatenate([spikes_a[idx_a], spikes_b[idx_b]]))


    @staticmethod
    def mix_fraction_of_spikes(spikes_a, spikes_b, fraction_a, fraction_b, rng, time_range=None):

        """ Picks fraction_a of spikes_a and fraction_b of spikes_b and returns sorted spike train

        Args:
            spikes_a (np.array) : Spike train A
            spikes_b (np.array) : Spike train B
            fraction_a (float) : Fraction of spikes in train A picked, e.g 0.4 means 40% of spikes are picked
            fraction_b (float) : Fraction of spikes in train B picked
            rng : Numpy rng object
            time_range : (start_times, end_times) for the different fractions
        """

        p_keep_a = np.zeros((np.size(spikes_a),))
        p_keep_b = np.zeros((np.size(spikes_b),))

        if time_range is None:
            assert np.size(fraction_a) == np.size(fraction_b) == 1

            assert 0 <= fraction_a <= 1 and 0 <= fraction_b <= 1

            p_keep_a[:] = fraction_a
            p_keep_b[:] = fraction_b
        else:
            assert len(time_range) == 2
            assert np.ndim(time_range[0]) == np.ndim(time_range[1])

            if np.ndim(time_range[0]) == 0:
                time_range = (np.array([time_range[0]]), np.array([time_range[1]]))

            if np.ndim(fraction_a) == 0:
                fraction_a = np.full(time_range[0].shape, fraction_a)
            else:
                fraction_a = np.array(fraction_a)

            if np.ndim(fraction_b) == 0:
                fraction_b = np.full(time_range[0].shape, fraction_b)
            else:
                fraction_b = np.array(fraction_b)

            assert np.size(fraction_a) == np.size(fraction_b) == np.size(time_range[0]) == np.size(time_range[1]), \
                f"Lengths must match for time_range start {time_range[0]}, end {time_range[1]}, " \
                f"fraction_a {fraction_a} and fraction_b {fraction_b}"
            assert np.logical_and(0 <= fraction_a, fraction_a <= 1).all() \
                and np.logical_and(0 <= fraction_b, fraction_b <= 1).all(), \
                f"Fractions must be between 0 and 1: {fraction_a}, {fraction_b}"

            try:
                for start, end, f_a, f_b in zip(*time_range, fraction_a, fraction_b):
                    idx_a = np.where(np.logical_and(start <= spikes_a, spikes_a <= end))[0]
                    idx_b = np.where(np.logical_and(start <= spikes_b, spikes_b <= end))[0]
                    p_keep_a[idx_a] = f_a
                    p_keep_b[idx_b] = f_b

            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

        keep_idx_a = np.where(p_keep_a >= rng.uniform(size=p_keep_a.shape))[0]
        keep_idx_b = np.where(p_keep_b >= rng.uniform(size=p_keep_b.shape))[0]

        return np.sort(np.concatenate([spikes_a[keep_idx_a], spikes_b[keep_idx_b]]))

    ############################################################################

    @staticmethod
    def cull_spikes(spikes, p_keep, rng, time_range=None):

        """
        Keeps a fraction of all spikes.

        Args:
            spikes: Spike train
            p_keep: Probability to keep each spike
            rng: Numpy random number stream
            time_range: If p_keep is vector, this specifies which part of those ranges each p_keep is for
        """

        if time_range is None:
            assert np.size(p_keep) == 1, f"If not time_range is given then p_keep must be a scalar. p_keep = {p_keep}"
            return spikes[rng.random(spikes.shape) < p_keep]
        else:
            if np.size(time_range[0]) == 1:
                old_time_range = time_range
                time_range = (np.array([time_range[0]]), np.array([time_range[1]]))

            if np.size(p_keep) == 1:
                p_keep = np.full(np.size(time_range[0]), p_keep)
            
            p_keep_spikes = np.zeros(spikes.shape)
            
            try:
                for p_k, start, end in zip(p_keep, time_range[0], time_range[1]):
                    idx = np.where(np.logical_and(start <= spikes, spikes <= end))[0]
                    p_keep_spikes[idx] = p_k
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()
              
        return spikes[rng.random(spikes.shape) < p_keep_spikes]

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

    def make_correlated_spikes(self,
                               freq, time_range, num_spike_trains, p_keep, rng,
                               population_unit_spikes=None,
                               ret_pop_unit_spikes=False, jitter_dt=None,
                               input_generator=None):

        """
        Make correlated spikes.

        Args:
            freq (float or str): frequency of spike train
            time_range (tuple): start time, end time of spike train
            num_spike_trains (int): number of spike trains to generate
            p_keep (float or list of floats): fraction of shared channel spikes to include in spike train, p_keep=1 (100% correlated)
            rng: Numpy random number stream
            population_unit_spikes
            ret_pop_unit_spikes (bool): if false, returns only spikes,
                                        if true returns (spikes, population unit spikes)
            jitter_dt (float): amount to jitter all spikes
            input_generator (str) : "poisson" (default) or "frequency_functon"
        """

        assert np.all(np.logical_and(0 <= p_keep, p_keep <= 1)), f"p_keep = {p_keep} should be between 0 and 1"

        if population_unit_spikes is None:
            population_unit_spikes = self.generate_spikes_helper(freq, time_range, rng=rng,
                                                                 input_generator=input_generator)
        spike_trains = []

        if input_generator == "poisson":
            pop_freq = np.multiply(freq, 1 - p_keep)
        elif input_generator == "frequency_function":
            pop_freq = freq
        else:
            assert False, f"Unknown input_generator {input_generator}"


        # TODO: For frequency_functions we can skip parts of the loop and directly request n_spike_trains
        #       That would speed up the frequency_function generation call

        for i in range(0, num_spike_trains):
            t_unique = self.generate_spikes_helper(frequency=pop_freq, time_range=time_range, rng=rng,
                                                   input_generator=input_generator)
            t_population_unit = self.cull_spikes(spikes=population_unit_spikes,
                                                 p_keep=p_keep, rng=rng,
                                                 time_range=time_range)

            spike_trains.append(SnuddaInput.mix_spikes([t_unique, t_population_unit]))

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

        """
        Generate uncorrelated spikes.

        Args:
            freq: frequency
            t_start: start time
            t_end: end time
            n_spike_trains: number of spike trains to generate
            rng: numpy random number stream
        """

        spike_trains = []

        for i in range(0, n_spike_trains):
            spike_trains.append(self.generate_poisson_spikes(freq, (t_start, t_end), rng))

        return spike_trains

    ############################################################################

    # If a timeRange (start,endtime) is given then all spike times will
    # be modulo duration, so if we jitter and they go to before start time,
    # they wrap around and appear at end of the timeline

    @staticmethod
    def jitter_spikes(spike_trains, dt, rng, time_range=None):

        """
        Jitter spikes in a spike train.

        If a time_range (start,end_time) is given then all spike times will
        be modulo duration, so if we jitter and they go to before start time,
        they wrap around and appear at end of the timeline


        Args:
            spike_trains: spike times
            dt: amount of jitter
            rng: Numpy random stream
            time_range (tuple): (start, end) see comment above about wrapping around edges.

        """

        jittered_spikes = []

        for i in range(0, len(spike_trains)):
            spikes = spike_trains[i] + rng.normal(0, dt, spike_trains[i].shape)

            # No modulo time jittering if list of times specified
            if time_range is not None and np.size(time_range[0]) == 1:
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

        """
        Raster plot of spike trains.

        Args:
            spike_times
            mark_spikes: list of spikes to mark
            mark_idx: index of neuron with spikes to mark
            title: title of plot
            fig_file: path to figure
            fig: matplotlib figure object
        """

        import matplotlib.pyplot as plt

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

    def read_network_config_file(self):

        """ Read network configuration JSON file."""

        self.write_log(f"Reading config file {self.network_config_file}")

        with open(self.network_config_file, 'r') as f:
            self.network_config = json.load(f, object_pairs_hook=OrderedDict)

        # This also loads random seed from config file while we have it open
        if self.random_seed is None:
            if "RandomSeed" in self.network_config and "input" in self.network_config["RandomSeed"]:
                self.random_seed = self.network_config["RandomSeed"]["input"]
                self.write_log(f"Reading random seed from config file: {self.random_seed}")
            else:
                # No random seed given, invent one
                self.random_seed = 1004
                self.write_log(f"No random seed provided, using: {self.random_seed}")
        else:
            self.write_log(f"Using random seed provided by command line: {self.random_seed}")

        if "PopulationUnits" in self.network_config:

            all_id = []
            for volume in self.network_config["PopulationUnits"]:
                if "unitID" in self.network_config["PopulationUnits"][volume]:
                    all_id += self.network_config["PopulationUnits"][volume]["unitID"]

            all_id = sorted(all_id)

            if "AllUnitID" in self.network_config["PopulationUnits"]:
                self.all_population_units = sorted(self.network_config["PopulationUnits"]["AllUnitID"])
                assert all_id == self.all_population_units, \
                    (f"Inconsistency: AllUnitID = {self.all_population_units}, "
                     f"but all units in unitID blocks = {all_id}")
            else:
                self.write_log("Missing AllUnitID tag, deriving it from unitID tag for volumes")
                self.all_population_units = all_id

        else:
            self.all_population_units = [0]

    def generate_seeds(self, num_states):

        """ From the master seed, generate a seed sequence for inputs. """

        ss = np.random.SeedSequence(self.random_seed)
        all_seeds = ss.generate_state(num_states + 1)

        return all_seeds[1:]  # First seed in sequence is reserved for master

    def get_master_node_rng(self):

        """ Get random number for master node, from master seed. """

        ss = np.random.SeedSequence(self.random_seed)
        master_node_seed = ss.generate_state(1)
        return np.random.default_rng(master_node_seed)

    ############################################################################

    def verify_correlation(self, spike_trains, dt=0):

        """
        Verify correlation. This function is slow.

        Args:
            spike_trains
            dt
        """

        # THIS FUNCTION IS VERY VERY SLOW

        corr_vec = []

        for si, s in enumerate(spike_trains):
            for s2i, s2 in enumerate(spike_trains):
                if si == s2i:
                    # No self comparison
                    continue

                corr_vec.append(self.estimate_correlation(s, s2, dt=dt))

        # print("corr = " + str(corrVec))
        self.write_log(f"meanCorr = {np.mean(corr_vec)}")

    ############################################################################

    @staticmethod
    def estimate_correlation(spikes_a, spikes_b, dt=0):

        """
        Estimate correlation between spikes_a and spikes_b, assuming correlation window of dt.

        Args:
            spikes_a
            spikes_b
            dt
        """

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
                                 synapse_density=None,
                                 num_spike_trains=None,
                                 cluster_size=None,
                                 cluster_spread=30e-6):

        """
        Return dendrite input location.

        Args:
            neuron_id: Neuron ID
            rng: Numpy random number stream
            synapse_density (str): Distance function f(d)
            num_spike_trains (int): Number of spike trains
            cluster_size (int): Size of each synaptic cluster (None = No clustering)
            cluster_spread (float): Spread of cluster along dendrite (in meters)
        """

        if synapse_density is None:
            synapse_density = "1"

        neuron_name = self.neuron_name[neuron_id]

        # self.write_log(f"self.network_config = {self.network_config}")
        # self.write_log(f"self.network_config['Neurons'] = {self.network_config['Neurons']}")
        # self.write_log(f"self.network_config['Neurons'][neuron_name] = {self.network_config['Neurons'][neuron_name]}")

        morphology_path = self.network_config["Neurons"][neuron_name]["morphology"]
        parameters_path = self.network_config["Neurons"][neuron_name]["parameters"]

        if "modulation" in self.network_config["Neurons"][neuron_name]:
            modulation_path = self.network_config["Neurons"][neuron_name]["modulation"]
        else:
            modulation_path = None

        mechanisms_path = self.network_config["Neurons"][neuron_name]["mechanisms"]

        parameter_id = self.neuron_info[neuron_id]["parameterID"]
        morphology_id = self.neuron_info[neuron_id]["morphologyID"]
        modulation_id = self.neuron_info[neuron_id]["modulationID"]

        parameter_key = self.neuron_info[neuron_id]["parameterKey"]
        morphology_key = self.neuron_info[neuron_id]["morphologyKey"]
        modulation_key = self.neuron_info[neuron_id]["modulationKey"]

        if neuron_name in self.neuron_cache:
            self.write_log(f"About to clone cache of {neuron_name}.")
            # Since we do not care about location of neuron in space, we can use get_cache_original
            morphology = self.neuron_cache[neuron_name].clone(parameter_id=parameter_id,
                                                              morphology_id=morphology_id,
                                                              parameter_key=parameter_key,
                                                              morphology_key=morphology_key,
                                                              modulation_id=None, position=None, rotation=None,
                                                              get_cache_original=True)
        else:
            self.write_log(f"Creating prototype {neuron_name}")
            morphology_prototype = NeuronPrototype(neuron_name=neuron_name,
                                                   snudda_data=self.snudda_data,
                                                   morphology_path=morphology_path,
                                                   parameter_path=parameters_path,
                                                   modulation_path=modulation_path,
                                                   mechanism_path=mechanisms_path,
                                                   neuron_path=None,
                                                   axon_stump_id_flag=self.axon_stump_id_flag)
            self.neuron_cache[neuron_name] = morphology_prototype
            morphology = morphology_prototype.clone(parameter_id=parameter_id,
                                                    morphology_id=morphology_id,
                                                    parameter_key=parameter_key,
                                                    morphology_key=morphology_key,
                                                    modulation_id=None, position=None, rotation=None,
                                                    get_cache_original=True)

        self.write_log(f"morphology = {morphology}")

        # input_info = self.neuron_cache[neuron_name].get_input_parameters(parameter_id=parameter_id,
        #                                                                  morphology_id=morphology_id,
        #                                                                  parameter_key=parameter_key,
        #                                                                  morphology_key=morphology_key)

        return morphology.dendrite_input_locations(synapse_density=synapse_density,
                                                   num_locations=num_spike_trains,
                                                   rng=rng,
                                                   cluster_size=cluster_size,
                                                   cluster_spread=cluster_spread)

    ############################################################################

    def setup_parallel(self):

        """ Setup worker nodes for parallel execution. """

        slurm_job_id = os.getenv("SLURM_JOBID")

        if slurm_job_id is None:
            self.slurm_id = 0
        else:
            self.slurm_id = int(slurm_job_id)

        if self.rc is not None:
            # http://davidmasad.com/blog/simulation-with-ipyparallel/
            # http://people.duke.edu/~ccc14/sta-663-2016/19C_IPyParallel.html
            self.write_log(f"Client IDs: {self.rc.ids}")
            self.d_view = self.rc.direct_view(targets='all')

            if self.logfile is not None:
                log_filename = self.logfile.name
                engine_logfile = [log_filename + "-" + str(x) for x in range(0, len(self.d_view))]
            else:
                engine_logfile = [None for x in range(0, len(self.d_view))]
        else:
            self.write_log("Running in serial")
            self.d_view = None
            return

        with self.d_view.sync_imports():
            from snudda.input.input import SnuddaInput

        self.d_view.push({"network_path": self.network_path,
                          "input_config_file": self.input_config_file,
                          "spike_data_filename": self.spike_data_filename,
                          "hdf5_network_file": self.hdf5_network_file,
                          "snudda_data": self.snudda_data,
                          "is_master": False,
                          "time": self.time,
                          "h5libver": self.h5libver,
                          "random_seed": self.random_seed,
                          "use_meta_input": self.use_meta_input,
                          "verbose": self.verbose,
                          "time_interval_overlap_warning": self.time_interval_overlap_warning})

        self.write_log(f"Scattering engineLogFile = {engine_logfile}")

        self.d_view.scatter('log_filename', engine_logfile, block=True)

        self.write_log(f"nl = SnuddaInput(network_path={self.network_path}"
                       f", snudda_data='{self.snudda_data}'"
                       f", input_config_file='{self.input_config_file}'"
                       f", spike_data_filename='{self.spike_data_filename}'"
                       f", hdf5_nework_file={self.hdf5_network_file}"
                       f", h5libver={self.h5libver}"
                       f", is_master=False "
                       f", random_seed={self.random_seed}"
                       f", use_meta_input={self.use_meta_input}"
                       f", verbose={self.verbose}"
                       f", time_interval_overlap_warning={self.time_interval_overlap_warning}"
                       f", time={self.time}, logfile='{log_filename[0]}')")

        cmd_str = ("global nl; nl = SnuddaInput(network_path=network_path, "
                   "snudda_data=snudda_data, "
                   "input_config_file=input_config_file, "
                   "spike_data_filename=spike_data_filename, "
                   "hdf5_network_file=hdf5_network_file, "
                   "use_meta_input=use_meta_input, "
                   "is_master=is_master, time=time, "
                   "h5libver=h5libver, "
                   "verbose=verbose, "
                   "time_interval_overlap_warning=time_interval_overlap_warning, "
                   "random_seed=random_seed, logfile=log_filename[0])")

        self.d_view.execute(cmd_str, block=True)

        self.write_log("Read network config on workers")
        cmd_str3 = "nl.read_network_config_file()"
        self.d_view.execute(cmd_str3, block=True)

        self.write_log("Workers set up")

    ############################################################################

    def check_sorted(self):

        """ Checks that spikes are in chronological order. """

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

        """ Plot spikes for neuron_id """

        self.write_log(f"Plotting spikes for neuronID: {neuron_id}")

        if neuron_id is None:
            neuron_id = self.neuron_input

        spike_times = []

        for nID in neuron_id:
            for inputType in self.neuron_input[nID]:
                for spikes in self.neuron_input[nID][inputType]["spikes"]:
                    spike_times.append(spikes)

        self.raster_plot(spike_times)

    ############################################################################

    def make_input_helper_parallel(self, args):

        """ Helper function for parallel input generation."""

        try:

            neuron_id, input_type, freq, start, end, synapse_density, num_spike_trains, \
            population_unit_spikes, jitter_dt, population_unit_id, conductance, correlation, mod_file, \
            parameter_file, parameter_list, random_seed, cluster_size, cluster_spread, \
            dendrite_location_override, input_generator, population_unit_fraction = args

            return self.make_input_helper_serial(neuron_id=neuron_id,
                                                 input_type=input_type,
                                                 freq=freq,
                                                 t_start=start,
                                                 t_end=end,
                                                 synapse_density=synapse_density,
                                                 num_spike_trains=num_spike_trains,
                                                 population_unit_spikes=population_unit_spikes,
                                                 jitter_dt=jitter_dt,
                                                 population_unit_id=population_unit_id,
                                                 conductance=conductance,
                                                 correlation=correlation,
                                                 mod_file=mod_file,
                                                 parameter_file=parameter_file,
                                                 parameter_list=parameter_list,
                                                 random_seed=random_seed,
                                                 cluster_size=cluster_size,
                                                 cluster_spread=cluster_spread,
                                                 dendrite_location=dendrite_location_override,
                                                 input_generator=input_generator,
                                                 population_unit_fraction=population_unit_fraction)

        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr, is_error=True)
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
                                 population_unit_spikes,
                                 jitter_dt,
                                 population_unit_id,
                                 conductance,
                                 correlation,
                                 mod_file,
                                 parameter_file,
                                 parameter_list,
                                 random_seed,
                                 cluster_size=None,
                                 cluster_spread=None,
                                 dendrite_location=None,
                                 input_generator=None,
                                 population_unit_fraction=1):

        """
        Generate poisson input.

        Args:
            neuron_id (int): Neuron ID to generate input for
            input_type: Input type
            freq: Frequency of input
            t_start: Start time of input
            t_end: End time of input
            synapse_density: Density function f(d), d=distance to soma along dendrite
            num_spike_trains: Number of spike trains
            jitter_dt: Amount of time to jitter all spikes
            population_unit_spikes: Population unit spikes
            population_unit_id: Population unit ID
            conductance: Conductance
            correlation: correlation
            mod_file: Mod file
            parameter_file: Parameter file for input synapses
            parameter_list: Parameter list (to inline parameters, instead of reading from file)
            random_seed: Random seed.
            cluster_size: Input synapse cluster size
            cluster_spread: Spread of cluster along dendrite (in meters)
            dendrite_location: Override location of dendrites, list of (sec_id, sec_x) tuples.
            input_generator: "poisson" or "frequency_function"
            population_unit_fraction: Fraction of population unit spikes used, 1.0=all correlation within population unit, 0.0 = only correlation within the particular neuron
        """

    # First, find out how many inputs and where, based on morphology and
        # synapse density

        time_range = (t_start, t_end)

        rng = np.random.default_rng(random_seed)

        if input_type.lower() == "VirtualNeuron".lower():
            # This specifies activity of a virtual neuron
            conductance = None

            assert num_spike_trains is None or num_spike_trains == 1, \
                (f"Virtual neuron {self.neuron_name[neuron_id]}"
                 f" should have only one spike train, fix nSpikeTrains in config")

            # Virtual neurons input handled through touch detection
            input_loc = None

            num_inputs = 1
            p_keep = np.divide(1, (num_inputs - np.sqrt(correlation) * (num_inputs - 1)))

            # !!! Pass the input_generator
            spikes = self.make_correlated_spikes(freq=freq,
                                                 time_range=time_range,
                                                 num_spike_trains=1,
                                                 p_keep=p_keep,
                                                 population_unit_spikes=population_unit_spikes,
                                                 jitter_dt=jitter_dt,
                                                 rng=rng,
                                                 input_generator=input_generator)
        else:

            if dendrite_location:
                self.write_log(f"Overriding input location for {input_type} on neuron_id={neuron_id}")
                sec_id, sec_x = zip(*dendrite_location)

                # TODO: Calculate the correct x,y,z and distance to soma
                x = y = z = dist_to_soma = np.zeros((len(sec_id),))
                input_loc = [(x, y, z), np.array(sec_id), np.array(sec_x), dist_to_soma]

            else:
                # (x,y,z), secID, secX, dist_to_soma
                input_loc = self.dendrite_input_locations(neuron_id=neuron_id,
                                                          synapse_density=synapse_density,
                                                          num_spike_trains=num_spike_trains,
                                                          rng=rng,
                                                          cluster_size=cluster_size,
                                                          cluster_spread=cluster_spread)

            num_inputs = input_loc[0].shape[0]

            if num_inputs > 0:
                p_keep = np.divide(1, (num_inputs - np.sqrt(correlation) * (num_inputs - 1)))
            else:
                p_keep = 0

            if population_unit_spikes is not None:
                neuron_correlated_spikes = self.generate_spikes_helper(frequency=freq, time_range=time_range, rng=rng,
                                                                       input_generator=input_generator)

                mother_spikes = SnuddaInput.mix_fraction_of_spikes(population_unit_spikes, neuron_correlated_spikes,
                                                                   population_unit_fraction, 1-population_unit_fraction,
                                                                   rng=rng, time_range=time_range)
            else:
                mother_spikes = population_unit_spikes

            self.write_log(f"Generating {num_inputs} inputs (correlation={correlation}, p_keep={p_keep}, "
                           f"population_unit_fraction={population_unit_fraction}) "
                           f"for {self.neuron_name[neuron_id]} ({neuron_id})")

            # OBS, n_inputs might differ slightly from n_spike_trains if that is given
            spikes = self.make_correlated_spikes(freq=freq,
                                                 time_range=time_range,
                                                 num_spike_trains=num_inputs,
                                                 p_keep=p_keep,
                                                 population_unit_spikes=mother_spikes,
                                                 jitter_dt=jitter_dt,
                                                 rng=rng,
                                                 input_generator=input_generator)

        # We need to pick which parameter set to use for the input also
        parameter_id = rng.integers(1e6, size=num_inputs)

        # We need to keep track of the neuronID, since it will all be jumbled
        # when doing asynchronous parallelisation
        return (neuron_id, input_type, spikes, input_loc, synapse_density, freq,
                jitter_dt, population_unit_id, conductance, correlation,
                time_range,
                mod_file, parameter_file, parameter_list, parameter_id)

    ############################################################################

    def write_log(self, text, flush=True, is_error=False, force_print=False):  # Change flush to False in future, debug

        """
        Writes to log file. Use setup_log first. Text is only written to screen if self.verbose=True,
        or is_error = True, or force_print = True.

        test (str) : Text to write
        flush (bool) : Should all writes be flushed to disk directly?
        is_error (bool) : Is this an error, always written.
        force_print (bool) : Force printing, even if self.verbose=False.
        """

        if self.logfile is not None:
            self.logfile.write(text + "\n")
            if flush:
                self.logfile.flush()

        if self.verbose or is_error or force_print:
            print(text, flush=True)

    ############################################################################


if __name__ == "__main__":
    print("Please do not call this file directly, use snudda command line")
    sys.exit(-1)
