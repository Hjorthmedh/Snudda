#
# This code reads the network created by Network_connect.py and set it
# up in memory
#
# mpiexec -n 4 python snudda_simulate.py
#
#
# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union's Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).

#
############################################################################

import copy
import gc
import json
import os
import re
import time
import timeit
# Plot all sections
# [neuron.h.psection(x) for x in neuron.h.allsec()]
from collections import OrderedDict

import h5py
import neuron
import numpy as np

from mpi4py import MPI  # This must be imported before neuron, to run parallel
from neuron import h  # , gui

import snudda.utils.memory
from snudda.neurons.neuromodulation_extracellular import ExtracellularNeuromodulation
from snudda.neurons.neuron_model_extended import NeuronModel
from snudda.simulate.nrn_simulator_parallel import NrnSimulatorParallel
from snudda.simulate.save_network_recording import SnuddaSaveNetworkRecordings
# If simulationConfig is set, those values override other values
from snudda.utils.load import SnuddaLoad
from snudda.utils.snudda_path import snudda_parse_path, get_snudda_data


# !!! Need to gracefully handle the situation where there are more workers than
# number of neurons, currently we get problem when adding the voltage saving
##############################################################################


class SnuddaSimulate(object):
    """ Simulate network """

    def __init__(self,
                 network_path=None,
                 network_file=None,
                 snudda_data=None,
                 input_file=None,
                 output_file=None,
                 verbose=False,
                 log_file=None,
                 disable_synapses=None,
                 disable_gap_junctions=None,
                 sample_dt=None,
                 use_rxd_neuromodulation=None,
                 simulation_config=None):

        """
        Constructor

        Args:
            network_path (str): Path to network directory
            network_file (str, optional): Path to network file, default network-synapses.hdf5 in network_path
            snudda_data (str, optional): Path to Snudda Data
            input_file (str, optional): Path to input file, default input-spikes.hdf5 in network_path
            output_file (str, optional): Path to output file, default simulation/output.hdf5
            verbose (bool): Extra printouts
            log_file (str): Path to logfile
            disable_gap_junctions (bool): Disable gap junctions, default False
            disable_synapses (bool): Disable synapses, default False
            simulation_config (str, optional): Path to config file with simulation info (including network_path)

        """

        self.verbose = verbose
        self.log_file = log_file
        self.sim_info = dict()

        if network_path:
            self.network_path = network_path
        elif network_file:
            self.network_path = os.path.dirname(network_file)
        else:
            self.write_log("No network network_path or network_file specified")
            self.network_path = None

        self.snudda_data = get_snudda_data(snudda_data=snudda_data,
                                           network_path=self.network_path)

        if not network_file and self.network_path is not None:
            alt_network_file = os.path.join(self.network_path, "network-synapses.hdf5")
            if os.path.isfile(alt_network_file):
                self.network_file = alt_network_file
            else:
                self.network_file = None
        else:
            self.network_file = network_file

        if self.network_file is None:
            self.write_log(f"Warning: no network_file defined.", is_error=True)

        if not input_file:
            default_input_file = os.path.join(self.network_path, "input-spikes.hdf5")

            if os.path.exists(default_input_file):
                self.input_file = default_input_file
            else:
                self.write_log("Warning: No external synaptic input file given!")
                self.input_file = None
        else:
            self.input_file = input_file

        if output_file is not None:
            self.output_file = output_file
        else:
            self.output_file = os.path.join(self.network_path, "simulation", "output.hdf5")

        # Init
        self.snudda_loader = None
        self.network_info = None
        self.synapses = None
        self.gap_junctions = None
        self.num_neurons = None
        self.config_file = None
        self.config = None
        self.is_virtual_neuron = None
        self.neuron_id = None
        self.neuron_id_on_node = None
        self.synapse_parameters = None
        self.use_rxd_neuromodulation = use_rxd_neuromodulation
        self.bath_application = dict()
        self.extracellular_regions = dict()  # TODO: This needs to be set when extracellular space is defined!

        self.sim_start_time = 0
        self.fih_time = None
        self.last_sim_report_time = 0
        self.sample_dt = sample_dt  # None means all values, 0.0005 means 0.5ms time resolution in saved files
        self.sim_dt = 0.025
        self.use_cvode = False

        self.pc = h.ParallelContext()

        self.print_error_once = dict()

        self.disable_synapses = False
        self.disable_gap_junctions = False

        self.current_injection_info = dict()
        self.current_clamps = dict()

        self.node_id = int(self.pc.id())
        self.total_nodes = int(self.pc.nhost())

        self.post_init_mods = dict()

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        print(f"MPI Rank: {rank}, Size: {size} -- NEURON: This is node {self.node_id} out of {self.total_nodes}")

        if simulation_config:

            print(f"Reading config: {simulation_config}")

            if type(simulation_config) == dict:
                self.sim_info = simulation_config

            elif os.path.isfile(simulation_config):
                print(f"Loading simulation_config from {simulation_config}")
                with open(simulation_config, "r") as f:
                    self.sim_info = json.load(f, object_pairs_hook=OrderedDict)
            else:
                print(f"Unable to find simulation_config file: {simulation_config}")
                import sys
                sys.exit(-1)

            if "log_file" in self.sim_info:

                if self.total_nodes > 1:
                    log_file_name = f"{self.sim_info['log_file']}-{self.node_id}"
                else:
                    log_file_name = self.sim_info["log_file"]

                self.log_file = open(log_file_name, "w")
                self.write_log(f"Using log file {self.log_file}")

            elif isinstance(self.log_file, str):
                self.log_file = open(self.log_file, "w")

            if "network_path" in self.sim_info:
                self.network_path = self.sim_info["network_path"]
                self.write_log(f"Network path: {self.network_file}")

            if "network_file" in self.sim_info:
                self.network_file = self.sim_info["network_file"]
                self.write_log(f"Network file: {self.network_file}")

            if "input_file" in self.sim_info:
                self.input_file = self.sim_info["input_file"]
                self.write_log(f"Input file: {self.input_file}")

            if "output_file" in self.sim_info:
                self.output_file = self.sim_info["output_file"]
                self.write_log(f"Output file: {self.output_file}")

            if "use_cvode" in self.sim_info:
                self.use_cvode = self.sim_info["use_cvode"]
                if self.use_cvode:
                    self.write_log(f"Using CVODE for NEURON simulation")

            if "sim_dt" in self.sim_info:
                self.sim_dt = self.sim_info["sim_dt"] * 1e3  # OBS, converted to ms for NEURON
                self.write_log(f"Setting simulation dt={self.sim_dt} ms")

            if "sample_dt" in self.sim_info:
                self.sample_dt = self.sim_info["sample_dt"]
                self.write_log(f"Sample dt: {self.sample_dt}")

            if "disable_synapses" in self.sim_info:
                self.disable_synapses = self.sim_info["disable_synapses"]
                self.write_log(f"Disable synapses: {self.disable_synapses}")

            if "disable_gap_junctions" in self.sim_info:
                self.disable_gap_junctions = self.sim_info["disable_gap_junctions"]
                self.write_log(f"Disable gap junctions: {self.disable_gap_junctions}")

            if "verbose" in self.sim_info:
                self.verbose = self.sim_info["verbose"]

            if "snudda_data" in self.sim_info:
                # Do not change this unless you know what you are doing
                self.snudda_data = self.sim_info["snudda_data"]
                self.write_log(f"Reading snudda_data {self.snudda_data} from config file", force_print=True)

            if "current_injection_file" in self.sim_info:
                current_file = self.sim_info["current_injection_file"]
                if not os.path.isfile(current_file):
                    raise ValueError(f"No such current injection file {current_file}")

                with open(current_file, "rt") as f:
                    self.current_injection_info = json.load(f)

                self.write_log(f"Current injection file: {current_file}")

            if "current_injection_info" in self.sim_info:
                # This is merged with current injection info read from file (above)
                self.current_injection_info |= self.sim_info["current_injection_info"]
                self.write_log(f"Updating current_injection_info from config file")

            if "post_init_modifications" in self.sim_info:
                self.post_init_mods = self.sim_info["post_init_modifications"]

        if self.log_file is None:
            self.log_file = os.path.join(self.network_path, "log", "simulation-log.txt")

        if isinstance(self.log_file, str):
            log_dir_name = os.path.dirname(self.log_file)
            self.create_dir(log_dir_name)

            self.log_file += f'-{int(self.pc.id())}'
            self.log_file = open(self.log_file, "w")

        # Earliest point we can write to log file, it needs to be opened.

        self.write_log(f"Using network_file: {self.network_file}")
        self.write_log(f"Using input_file: {self.input_file}")
        self.write_log(f"Using output_file: {self.output_file}")

        if self.log_file is not None:
            self.write_log(f"Using logFile: {self.log_file.name}")

        # !!! What values to use for synaptic weight and synapse delay?
        # !!! different for AMPA and GABA?
        self.synapse_weight = 10.0  # microsiemens
        self.synapse_delay = 1  # ms
        self.spike_threshold = -20  # TODO: Let each neuron type have individual spike threshold, based on what config file says.
        self.axon_speed = 0.8  # Tepper and Lee 2007, Wilson 1986, Wilson 1990
        # refs taken from Damodaran et al 2013

        if disable_synapses is not None:
            self.disable_synapses = disable_synapses

        if disable_gap_junctions is not None:
            self.disable_gap_junctions = disable_gap_junctions

        self.synapse_type_lookup = {1: "GABA", 2: "AMPA_NMDA", 3: "gap_junction"}

        self.neurons = {}
        self.sim = None
        self.neuron_nodes = []  # Is this used?

        self.virtual_neurons = {}

        self.synapse_dict = dict()  # Avoid premature garbage collection
        self.i_stim = []
        self.v_clamp_list = []
        self.gap_junction_dict = dict()
        self.external_stim = dict()
        self.t_save = []
        self.i_save = []
        self.i_key = []

        self.input_data = None

        self.gap_junction_next_gid = 0  # Are these gids separate from cell gids?

        self.check_id_recordings = []  # Prevent segmentation fault due to garbage collection of spike id

        # Make sure the output dir exists, so we don't fail at end because we cant write file
        self.create_dir(os.path.join("save", "traces"))

        self.conv_factor = {"tauR": 1e3,
                            "tauF": 1e3,
                            "tau": 1e3}

        # We need to initialise random streams, see Lytton el at 2016 (p2072)

        if self.network_file is not None:
            self.load_network_info(self.network_file)
        else:
            self.write_log("No network path or file specified, not loading network.")

        self.record = SnuddaSaveNetworkRecordings(output_file=self.output_file, network_data=self.network_info,
                                                  sample_dt=self.sample_dt, node_id=self.node_id)
        self.record.add_unit(data_type="voltage", target_unit="V", conversion_factor=1e-3)
        self.record.add_unit(data_type="synaptic_current", target_unit="A", conversion_factor=1e-9)
        self.record.add_unit(data_type="spikes", target_unit="s", conversion_factor=1e-3)
        self.record.add_unit(data_type="time", target_unit="s", conversion_factor=1e-3)
        # TODO: Add more units as needed https://www.neuron.yale.edu/neuron/static/docs/units/unitchart.html

    # def __del__(self):
    #     print("Destructor called -- explicitly deleting neurons from NEURON")
    #     for n in self.neurons.values():
    #         for k, v in n.__dict__.items():
    #            del v

    def setup(self):

        """ Setup simulation """
        self.check_memory_status()
        self.distribute_neurons()
        self.pc.barrier()

        # This activates RxD functionality, if needed
        self.pre_setup_parse_sim_info()

        # TODO: Setup extracellular space
        self.setup_extracellular_region()

        self.setup_neurons()
        self.check_memory_status()
        self.pc.barrier()

        self.connect_network()
        self.check_memory_status()
        self.pc.barrier()

        self.setup_parse_sim_info()

    def pre_setup_parse_sim_info(self):

        if self.use_rxd_neuromodulation is None:
            if "use_rxd_neuromodulation" in self.sim_info:
                self.use_rxd_neuromodulation = self.sim_info["use_rxd_neuromodulation"]
            else:
                # Setting default to False
                self.use_rxd_neuromodulation = False

        if self.use_rxd_neuromodulation:
            self.write_log(f"RxD for neuromodulation: {'ENABLED' if self.use_rxd_neuromodulation else 'DiSABLED'}.")

        if "rxd_enable_extracellular" in self.sim_info:
            import neuron.rxd as rxd
            rxd.options.enable.extracellular = self.sim_info["rxd_enable_extracellular"]

    def setup_parse_sim_info(self):

        if self.sim_info:

            if "record_all_soma" in self.sim_info and self.sim_info["record_all_soma"]:
                self.add_volt_recording_soma()

            if "record_soma" in self.sim_info:
                record_soma_cell_id = np.array(self.sim_info["record_soma"], dtype=int)
                self.add_volt_recording_soma(cell_id=record_soma_cell_id)

            if "record_all_compartments" in self.sim_info:
                raise DeprecationWarning("record_all_compartments deprecated, please use record_voltate_all_compartments")
                record_comp_cell_id = np.array(self.sim_info["record_all_compartments"], dtype=int)
                self.add_volt_recording_all(cell_id=record_comp_cell_id)

            if "record_voltage_all_compartments" in self.sim_info:
                record_comp_cell_id = np.array(self.sim_info["record_voltage_all_compartments"], dtype=int)
                self.add_volt_recording_all(cell_id=record_comp_cell_id)

            if "record_all_synapses" in self.sim_info:
                raise DeprecationWarning("record_all_synapses deprecated, please use record_current_all_synapses")
                record_syn_cell_id = np.array(self.sim_info["record_all_synapses"], dtype=int)
                self.add_synapse_current_recording_all(record_syn_cell_id)

            if "record_current_all_synapses" in self.sim_info:
                record_syn_cell_id = np.array(self.sim_info["record_current_all_synapses"], dtype=int)
                self.add_synapse_current_recording_all(record_syn_cell_id)

            if "record_rxd_species_concentration_all_compartments" in self.sim_info and self.use_rxd_neuromodulation:

                rec_info = self.sim_info["record_rxd_species_concentration_all_compartments"]

                if type(rec_info[0]) == str:
                    # Only one species
                    rec_info = [rec_info]

                for record_rx_species, record_rxd_neuron_id in rec_info:
                    self.write_log(f"Recording {record_rx_species} from neuron_id = {record_rxd_neuron_id}")
                    for rxd_neuron_id in record_rxd_neuron_id:
                        self.add_rxd_internal_concentration_recording_all(record_rx_species, rxd_neuron_id)

            if "record_rxd_species_all" in self.sim_info and self.use_rxd_neuromodulation:
                rxd_record_neuron_id = self.sim_info["record_rxd_species_all"]

                self.add_rxd_internal_concentration_recording_all_species(neuron_id=rxd_record_neuron_id)

            if "record_rxd_species_soma" in self.sim_info and self.use_rxd_neuromodulation:
                rxd_record_neuron_id = self.sim_info["record_rxd_species_soma"]

                self.add_rxd_internal_concentration_recording_all_species(neuron_id=rxd_record_neuron_id,
                                                                          include_dendrites=False)

            if "record_rxd_extracellular" in self.sim_info and self.use_rxd_neuromodulation:
                for region_name, ecs_data_parent in self.sim_info["record_rxd_extracellular"].items():
                    for ecs_data in ecs_data_parent:
                        for species in ecs_data["species"]:
                            for xyz in ecs_data["location"]:
                                self.add_rxd_extracellular_concentration_recording(species=species,
                                                                                   volume_id=region_name,
                                                                                   xyz=xyz)

            if "record_density_mechanism" in self.sim_info:
                record_info = self.sim_info["record_density_mechanism"]

                for key_name, density_data in record_info.items():
                    density_mechanism_name, variable_name = key_name.split(".")
                    neuron_id = density_data["neuron_id"]
                    sec_id = density_data["section_id"]
                    sec_x = density_data["section_x"]

                    if not (len(neuron_id) == len(sec_id) == len(sec_x)):
                        raise ValueError(f"neuron_id, section_id and section_x must be lists of same length in experiment config.\n"
                                         f"{record_info = }")

                    for nid, sid, sex in zip(neuron_id, sec_id, sec_x):

                        self.add_density_mechanism_recording(neuron_id=nid,
                                                             sec_id=sid, sec_x=sex,
                                                             density_mechanism=density_mechanism_name,
                                                             variable=variable_name)
            if "bath_application" in self.sim_info:
                for species_name, bath_info in self.sim_info["bath_application"].items():

                    bath_time = np.array(bath_info["time"])
                    bath_conc = np.array(bath_info["concentration"])
                    neuron_id = bath_info.get("neuron_id", None)
                    interpolate_bath = bath_info.get("interpolate", False)

                    self.add_bath_application(species_name=species_name,
                                              concentration=bath_conc,
                                              time=bath_time,
                                              neuron_id=neuron_id,
                                              interpolate=interpolate_bath)

            # Add any current injections that are specified
            self.parse_current_injection_info()

        # Do we need blocking call here, to make sure all neurons are setup
        # before we try and connect them

        # READ ABOUT PARALLEL NEURON

        # https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/network/parcon.html#paralleltransfer

    ############################################################################

    def load_mechanisms(self):

        """ Load the mechanisms. """

        try:
            if os.path.exists("nrnmech.dll"):
                self.write_log(f"Loading nrnmech.dll")
                h.nrn_load_dll("nrnmech.dll")
            elif os.path.exists("x86_64"):
                self.write_log(f"Loading x86_64/.libs/libnrnmech.so")
                h.nrn_load_dll("x86_64/.libs/libnrnmech.so")
            elif os.path.exists("aarch64"):
                self.write_log(f"Loading aarch64/.libs/libnrnmech.so")
                h.nrn_load_dll("aarch64/.libs/libnrnmech.so")
            elif os.path.exists("arm64"):
                self.write_log(f"Loading arm64/.libs/libnrnmech.so")
                h.nrn_load_dll("arm64/.libs/libnrnmech.so")
            else:
                self.write_log("No compiled mechanisms found. If you use custom mechanisms you need to run nrnivmodl")
        except:
            import traceback
            self.write_log(f"Error while loading mechanisms:\n{traceback.format_exc()}")

    ############################################################################

    def load_network_info(self, network_file=None, config_file=None):

        """
        Load network info from network data file (default network-synapses.hdf5)

        Args:
            network_file (str, optional): Path to network data, default network-synapses.hdf5
            config_file (str, optional): Path to network config file, default network-config.json

        """

        if network_file:
            self.network_file = network_file
        else:
            network_file = self.network_file

        self.write_log(f"Worker {int(self.pc.id())} : Loading network from {network_file}")

        from snudda.utils.load import SnuddaLoad
        self.snudda_loader = SnuddaLoad(network_file)
        self.network_info = self.snudda_loader.data

        self.synapses = self.network_info["synapses"]
        self.gap_junctions = self.network_info["gap_junctions"]

        # We are only passed information about neurons on our node if
        # SplitConnectionFile was run, so need to use nNeurons to know
        # how many neurons in total
        self.num_neurons = self.network_info["num_neurons"]

        if config_file is None:
            config_file = snudda_parse_path(self.network_info["config_file"], self.snudda_data)

        self.config_file = config_file
        self.write_log(f"Loading config file {config_file}")

        # Add checks to see that config file and network_file matches

        import json
        with open(config_file, 'r') as cf:
            self.config = json.load(cf, object_pairs_hook=OrderedDict)

        # I do not know if the gap junction GIDs are a separate entity from the
        # neuron cell GIDs, so to be on safe side, let's make sure they
        # do not overlap
        self.gap_junction_next_gid = self.num_neurons + 100000000

        # Make a bool array indicating if cells are virtual or not
        self.is_virtual_neuron = np.array([n["virtual_neuron"] for n in self.network_info["neurons"]], dtype=bool)

    ############################################################################

    def distribute_neurons_OLD(self):

        """
        Distribute neurons between workers.
        This code is run on all workers, will generate different lists on each
        """

        # This code is run on all workers, will generate different lists on each
        self.write_log("Distributing neurons.")

        assert self.num_neurons >= int(self.pc.nhost()), \
            f"Do not allocate more workers ({int(self.pc.nhost())}) than there are neurons ({self.num_neurons})."

        self.neuron_id = range(int(self.pc.id()), self.num_neurons, int(self.pc.nhost()))

        self.neuron_id_on_node = np.zeros((self.num_neurons,), dtype=bool)
        self.neuron_id_on_node[self.neuron_id] = True

        # TODO: Change to these ranges: range_borders = np.linspace(0, num_neurons, n_workers + 1).astype(int)
        #       will be faster, because of new numbering of neurons.

        range_borders = np.linspace(0, self.num_neurons, self.pc.nhost() + 1).astype(int)
        start_pos = range_borders[0]
        neuron_nodes = []
        for idx, end_pos in enumerate(range_borders[1:]):
            neuron_nodes += [idx] * (end_pos - start_pos)
            start_pos = end_pos

        self.neuron_nodes = neuron_nodes

    def distribute_neurons(self):
        # First distribute non-virtual neurons, then distribute the virtual neurons
        self.write_log("Distributing neurons (#2).")

        real_neuron_id = np.array([x["neuron_id"] for x in self.network_info["neurons"] if not x["virtual_neuron"]])
        virtual_neuron_id = np.array([x["neuron_id"] for x in self.network_info["neurons"] if x["virtual_neuron"]])

        r_idx = range(int(self.pc.id()), len(real_neuron_id), int(self.pc.nhost()))
        v_idx = range(int(self.pc.id()), len(virtual_neuron_id), int(self.pc.nhost()))

        neuron_id = list(set(real_neuron_id[r_idx]).union(set(virtual_neuron_id[v_idx])))
        self.neuron_id = np.array(neuron_id)

        self.neuron_id_on_node = np.zeros((self.num_neurons,), dtype=bool)
        self.neuron_id_on_node[self.neuron_id] = True

        self.write_log(f"Worker {int(self.pc.id())} has real neurons {real_neuron_id[r_idx] = }\n"
                       f"{virtual_neuron_id[v_idx] =}")

        assert len(real_neuron_id) >= int(self.pc.nhost()), \
            f"Do not allocate more workers ({int(self.pc.nhost())}) than there are real neurons ({len(real_neuron_id)})."

    ############################################################################

    # This requires self.sim to be defined

    def load_synapse_parameters(self):

        """
        Load synapse parameters. This requires self.sim to be defined.
        """

        # We need to load all the synapse parameters
        self.synapse_parameters = dict()

        for (pre_type, post_type) in self.network_info["connectivity_distributions"]:

            syn_data = self.network_info["connectivity_distributions"][pre_type, post_type]

            for syn_type in syn_data:

                synapse_type_id = syn_data[syn_type]["channel_model_id"]
                info_dict = syn_data[syn_type]

                if synapse_type_id == 3:
                    # Gap junctions, skip parameters
                    continue

                if "channel_parameters" in info_dict and info_dict["channel_parameters"] is not None:
                    channel_param_dict = copy.deepcopy(info_dict["channel_parameters"])
                    mod_file = channel_param_dict["mod_file"]

                    # TODO: Sanity check on the mod_file string
                    if mod_file:
                        eval_str = f"self.sim.neuron.h.{mod_file}"
                        channel_module = eval(eval_str)  # If this fails, check that NEURON modules are compiled
                    else:
                        self.write_log(f"Empty mod_file field for {pre_type} -> {post_type} synapses. This channel is IGNORED.", force_print=True)
                        channel_module = None

                    # These are not variables to set in the mod_file
                    if "mod_file" in channel_param_dict:
                        del channel_param_dict["mod_file"]

                    if "parameter_file" in channel_param_dict:
                        del channel_param_dict["parameter_file"]

                else:
                    assert False, (f"No channel module specified for {pre_type}->{post_type} synapses, "
                                   f"type ID={synapse_type_id}")

                if "parameter_file" in info_dict["channel_parameters"] \
                        and info_dict["channel_parameters"]["parameter_file"] is not None:
                    par_file = snudda_parse_path(info_dict["channel_parameters"]["parameter_file"], self.snudda_data)

                    with open(par_file, "r") as f:
                        par_data_dict = json.load(f, object_pairs_hook=OrderedDict)

                    # Save data as a list, we don't need the keys
                    par_data = []
                    for pd in par_data_dict:

                        if "synapse" in par_data_dict[pd]:
                            # Add channel parameters specified in network file, however
                            # any values in the synapse parameter file will overwrite them
                            p_dict = copy.deepcopy(channel_param_dict)
                            for x in par_data_dict[pd]["synapse"]:
                                p_dict[x] = par_data_dict[pd]["synapse"][x]

                            par_data.append(p_dict)

                        else:
                            self.write_log(f"WARNING: Old data format in parameter file {par_file}")

                            p_dict = copy.deepcopy(channel_param_dict)
                            for x in par_data_dict[pd]:
                                p_dict[x] = par_data_dict[pd][x]

                            par_data.append(p_dict)
                elif len(channel_param_dict) > 0:
                    par_data = [channel_param_dict]
                else:
                    par_data = None

                self.synapse_parameters[synapse_type_id] = (channel_module, par_data)

    ############################################################################

    def setup_extracellular_region(self):

        # What to do if self.sim_info is None?

        if self.sim_info is None or "rxd_enable_extracellular" not in self.sim_info \
                or not self.sim_info["rxd_enable_extracellular"]:
            # RxD extracellular not enabled

            self.write_log(f"RxD extracellular not enabled.")

            return

        for region_name, region_data in self.config["regions"].items():
            if "extracellular_space" in region_data:
                extracellular_info = region_data["extracellular_space"]

                extracellular_dx = extracellular_info.get("dx", None)
                if not isinstance(extracellular_dx, (list, tuple, np.ndarray)):
                    extracellular_dx = np.full(shape=(3, ), fill_value=extracellular_dx)

                extracellular_padding = extracellular_info.get("padding", None)
                extracellular_config = extracellular_info.get("extracellular_config", None)
                region_mesh = region_data.get("volume", {}).get("mesh_file", None)

                self.write_log(f"Setting up extracellular space for region {region_name}")

                self.extracellular_regions[region_name] = \
                    ExtracellularNeuromodulation(sim=self, volume_id=region_name,
                                                 padding=extracellular_padding, dx=extracellular_dx)

                self.extracellular_regions[region_name].load_json(config_path=extracellular_config)

        # TODO:
        # 1. Iterate through network_config file, to find out which regions have
        #    extracellular space defined. (maybe allow it to be other types of regions also)
        # 2. Setup each region, need to check how much padding is needed
        #    self.extracellular_region[region_name] = XXXXX
        # 3. Instantiate regions
        # 4. Check that neurons are able to couple to the regions

    ############################################################################

    def setup_neurons(self):

        """
        Instantiates neurons, saves them in self.neurons[neuronID]
        """

        self.write_log("Setup neurons")

        self.sim = NrnSimulatorParallel(cvode_active=self.use_cvode)

        # We need to load all the synapse parameters
        self.load_synapse_parameters()

        # The neurons this node is responsible for is in self.neuronID
        for ID in self.neuron_id:

            name = self.network_info["neurons"][ID]["name"]
            neuron_type = self.network_info["neurons"][ID]["type"]

            region = self.network_info["neurons"][ID]["volume_id"]
            position = self.network_info["neurons"][ID]["position"]
            rotation = self.network_info["neurons"][ID]["rotation"]

            # We need to get morphology from network_info, since it can now be redefined for bent morphologies
            morph = snudda_parse_path(self.network_info["neurons"][ID]["morphology"], self.snudda_data)

            neuron_path = snudda_parse_path(self.network_info["neurons"][ID]["neuron_path"], self.snudda_data)

            param = os.path.join(neuron_path, "parameters.json")
            mech = os.path.join(neuron_path, "mechanisms.json")

            if not self.use_rxd_neuromodulation:
                modulation = None

            elif "modulation" in self.network_info["neurons"][ID]:
                # TODO: This is old network_config location for "neurons", deprecate?
                raise DeprecationWarning(f"This is old location for 'neurons' in network_info, remove!")

                modulation = self.network_info["neurons"][ID]["modulation"]

                if not os.path.isfile(modulation):
                    raise ValueError(f"Missing modulation file {modulation} "
                                     f"for neuron {self.network_info['neurons'][ID]['name']}")

            elif neuron_type in self.network_info["config"]["regions"][region]["neurons"] and \
                    "modulation" in self.network_info["config"]["regions"][region]["neurons"][neuron_type]:
                modulation = self.network_info["config"]["regions"][region]["neurons"][neuron_type]["modulation"]

                if not os.path.isfile(modulation):
                    modulation = os.path.join(neuron_path, modulation)

            else:
                modulation = os.path.join(neuron_path, "modulation.json")

                if not os.path.isfile(modulation):
                    modulation = None

            if not self.use_rxd_neuromodulation:
                reaction_diffusion_file = None

            elif "reaction_diffusion_file" in self.network_info["neurons"][ID]:
                reaction_diffusion_file = SnuddaLoad.to_str(self.network_info["neurons"][ID]["reaction_diffusion_file"])

                if reaction_diffusion_file is not None and not os.path.isfile(reaction_diffusion_file):
                    raise ValueError(f"Missing RxD reaction diffusion file {reaction_diffusion_file} "
                                     f"for neuron {self.network_info['neurons'][ID]['name']}")

            else:
                reaction_diffusion_file = os.path.join(neuron_path, "reaction_diffusion.json")

                if not os.path.isfile(reaction_diffusion_file):
                    reaction_diffusion_file = None

            meta_file = snudda_parse_path(os.path.join(neuron_path, "meta.json"), self.snudda_data)
            axon_length = 60e-6
            axon_nseg_frequency = 40e-6
            replace_axon_diameter = None
            replace_axon_myelin_length = None
            replace_axon_myelin_diameter = None

            if os.path.isfile(meta_file):
                with open(meta_file, "r") as mf:
                    meta_data = json.load(mf)

                meta_parameter_key = self.network_info["neurons"][ID]["parameter_key"]
                meta_morphology_key = self.network_info["neurons"][ID]["morphology_key"]

                if meta_parameter_key in meta_data:
                    if meta_morphology_key in meta_data[meta_parameter_key]:
                        if "axon_stump" in meta_data[meta_parameter_key][meta_morphology_key]:
                            replace_info = meta_data[meta_parameter_key][meta_morphology_key]["axon_stump"]
                            # axon_length is normally a scalar, but if axon_diameter given
                            # then axon_length must also be a list
                            axon_length = replace_info.get("axon_length", 60e-6)
                            axon_nseg_frequency = replace_info.get("axon_nseg_frequency", 40e-6)

                            # Special treatment for Robert
                            replace_axon_diameter = replace_info.get("axon_diameter", None)

                            # Special treatment for Wilhelm
                            replace_axon_myelin_length = replace_info.get("axon_myelin_length", None)
                            replace_axon_myelin_diameter = replace_info.get("axon_myelin_diameter", None)

            # Obs, neurons is a dictionary
            if self.network_info["neurons"][ID]["virtual_neuron"]:

                if self.input_data is None:
                    self.write_log(f"Using {self.input_file} for virtual neurons")
                    self.input_data = h5py.File(snudda_parse_path(self.input_file, self.snudda_data), 'r')

                name = self.network_info["neurons"][ID]["name"]
                spikes = self.input_data["input"][str(ID)]["activity"]["spikes"][()] * 1e3  # s -> ms for NEURON

                # Creating NEURON VecStim and vector
                # https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3125
                vs = h.VecStim()
                v = h.Vector(spikes.size)
                v.from_python(spikes)
                vs.play(v)

                self.virtual_neurons[ID] = dict([])
                self.virtual_neurons[ID]["spikes"] = (v, vs, spikes)
                self.virtual_neurons[ID]["name"] = name

                self.pc.set_gid2node(ID, int(self.pc.id()))

                nc = h.NetCon(vs, None)
                self.pc.cell(ID, nc, 1)  # The 1 means broadcast spikes to other machines

            else:
                # A real neuron (not a virtual neuron that just provides input)
                parameter_key = self.network_info["neurons"][ID]["parameter_key"]
                morphology_key = self.network_info["neurons"][ID]["morphology_key"]

                # TODO: Modulation key currently has no USE -- deprecated? Remove?
                modulation_key = self.network_info["neurons"][ID]["modulation_key"]
                volume_id = self.network_info["neurons"][ID]["volume_id"]

                self.neurons[ID] = NeuronModel(param_file=param,
                                               morph_path=morph,
                                               mech_file=mech,
                                               cell_name=name,
                                               modulation_file=modulation,
                                               reaction_diffusion_file=reaction_diffusion_file,
                                               parameter_key=parameter_key,
                                               morphology_key=morphology_key,
                                               modulation_key=modulation_key,
                                               use_rxd_neuromodulation=self.use_rxd_neuromodulation,
                                               replace_axon_length=axon_length,
                                               replace_axon_nseg_frequency=axon_nseg_frequency,
                                               replace_axon_diameter=replace_axon_diameter,
                                               replace_axon_myelin_length=replace_axon_myelin_length,
                                               replace_axon_myelin_diameter=replace_axon_myelin_diameter,
                                               position=position,
                                               rotation=rotation,
                                               volume_id=volume_id)

                # Register ID as belonging to this worker node
                try:
                    # print(f"Debug:: Neuron {ID} on node {int(self.pc.id())}")
                    self.pc.set_gid2node(ID, int(self.pc.id()))
                except:
                    self.write_log(f"pc.set_gid2node failed ID = {ID}, {int(self.pc.id())}")
                    import traceback
                    self.write_log(traceback.format_exc())
                    import pdb
                    pdb.set_trace()

                self.write_log(f"Node {int(self.pc.id())} - cell {ID} {name}")

                # We need to instantiate the cell
                self.neurons[ID].instantiate(sim=self.sim, extracellular_regions=self.extracellular_regions)
                self.set_resting_voltage(ID)

                # !!! DIRTY FIX for
                # https://github.com/BlueBrain/BluePyOpt/blob/master/bluepyopt/ephys/morphologies.py
                # This is likely the offending line, that pushes a segment to the stack
                # --> sim.neuron.h.execute('create axon[2]', icell)

                self.write_log("!!! Popping extra segment from neuron -- temp fix!")
                h.execute("pop_section()")

                # !!! END OF DIRTY FIX

                # !!! Connect a netcon and register it, taken from ballandstick's connect2target function
                nc = h.NetCon(self.neurons[ID].icell.axon[0](0.5)._ref_v,
                              None,
                              sec=self.neurons[ID].icell.axon[0])
                nc.threshold = 10

                self.pc.cell(ID, nc, 1)  # The 1 means broadcast spikes to other machines

                # print(f"Spike threshold for neuron {ID} is {self.pc.threshold(ID)}")
                # self.pc.threshold(ID, new_threshold)  # if we want to update threshold, obs in mV

                # Record all spikes
                t_spikes = h.Vector()
                id_spikes = h.Vector()

                self.pc.spike_record(ID, t_spikes, id_spikes)
                self.check_id_recordings.append((ID, id_spikes))
                self.record.register_spike_data(neuron_id=ID, data=t_spikes, sec_id=-1, sec_x=0.5)

        # RxD is slow when doing species.nodes call, so we cache it...
        for neuron_id in self.neuron_id:
            if not self.network_info["neurons"][neuron_id]["virtual_neuron"] \
                    and self.neurons[neuron_id].modulation is not None:
                self.neurons[neuron_id].modulation.build_node_cache()

        # This allows us to modify ion channel conductance on the fly before runnigg simulation
        # Can be useful to e.g. increase KIR channel conductance
        self.post_initialisation_modifications()

    ############################################################################

    def connect_network(self):

        """ Connect neurons through synapses and gap junctions in network."""

        self.pc.barrier()

        # Add gap junctions
        if self.disable_gap_junctions:
            self.write_log("!!! Gap junctions disabled.", force_print=True)
        else:
            self.write_log("Adding gap junctions.")

            gap_junction_count = self.connect_network_gap_junctions_local()
            self.pc.setup_transfer()

            # Divide by two, since each gap junction is added twice (source and dest side)
            total_gap_junction_count = np.sum(self.pc.py_allgather(gap_junction_count)) / 2

            if total_gap_junction_count != self.gap_junctions.shape[0]:
                self.write_log(f"ERROR: Added only {total_gap_junction_count} out of {self.gap_junctions.shape[0]} gap junctions"
                               f"({np.sum(self.is_virtual_neuron)} virtual neurons)",
                               is_error=True)
            elif self.pc.id() == 0:
                self.write_log(f"Added {total_gap_junction_count} gap junctions to simulation ({self.gap_junctions.shape[0]} total)",
                               force_print=True)

        # Add synapses
        if self.disable_synapses:
            self.write_log("!!! Synapses disabled.", force_print=True)
        else:
            self.write_log("Adding synapses.")
            synapse_count = self.connect_network_synapses()
            self.write_log(f"Added {synapse_count} on worker {self.pc.id()}")
            total_synapse_count = np.sum(self.pc.py_allgather(synapse_count))

            if total_synapse_count != self.synapses.shape[0]:
                self.write_log(f"ERROR: Added only {total_synapse_count} out of {self.synapses.shape[0]} synapses! "
                               f"({np.sum(self.is_virtual_neuron)} virtual neurons)",
                               is_error=True)
            elif self.pc.id() == 0:
                self.write_log(f"Added {total_synapse_count} synapses to simulation ({self.synapses.shape[0]} total)",
                               force_print=True)

        self.pc.barrier()

    def set_new_output_file(self, new_output_file):
        self.write_log(f"Setting output file to {new_output_file}")
        self.output_file = new_output_file
        self.record.set_new_output_file(new_output_file)

    def reenable_synapses(self, new_output_file=None):

        assert self.disable_synapses, f"You can only reenable synapses if they were previously disabled"
        self.write_log("Re-enabling previously disabled synapses", force_print=True)

        self.disable_synapses = False
        self.connect_network_synapses()

        if new_output_file:
            self.set_new_output_file(new_output_file)

    def reenable_gap_junctions(self, new_output_file=None):
        assert self.disable_gap_junctions, f"You can only reenable gap junctions if they were previously disabled"
        self.write_log("Re-enabling previously disabled gap junctions", force_print=True)

        self.disable_gap_junctions = False
        self.connect_network_gap_junctions_local()
        self.pc.setup_transfer()

        if new_output_file:
            self.set_new_output_file(new_output_file)

    ############################################################################

    def connect_network_synapses(self):

        """ Connects neurons with synapses in network. """

        self.write_log("connect_network_synapses")

        # This loops through all the synapses, and connects the relevant ones
        next_row = 0
        next_row_set = self.find_next_synapse_group(next_row)

        synapse_count = 0

        while next_row_set is not None:
            # Add the synapses to the neuron
            synapse_count += self.connect_neuron_synapses(start_row=next_row_set[0], end_row=next_row_set[1])

            # Find the next group of synapses
            next_row = next_row_set[1]  # 2nd number was not included in range
            next_row_set = self.find_next_synapse_group(next_row)

        return synapse_count

    ############################################################################

    # This function starts at next_row, then returns all synapses onto a neuron which is located on the worker

    # This works for synapses, but it will not work for gap junctions, because
    # we need to connect the gap junctions from both sides

    # --- perhaps rewrite this as an iterator
    # --- also, see comment in loop below, to use boolean array to be faster

    def find_next_synapse_group(self, next_row=0):

        """
        Synapses are sorted by destination neuron (and then source neuron), this method starts from next_row
        and find the next range of synapses that have the same source and destination.

        Args:
            next_row (int): Row in the synapse matrix to start from
        """

        synapses = self.synapses

        try:
            num_syn_rows = synapses.shape[0]
        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr, is_error=True)

            assert False, ("find_next_synapse_group: If synapses was not loaded into memory, your problem is probably "
                           "that the HDF5 file that holds the synapses were closed. Sorry.")

        if next_row >= num_syn_rows:
            # No more synapses to get
            return None

        # The synapse matrix is sorted on dest_id, ascending order
        # We also assume that self.neuron_id is sorted in ascending order

        start_row = None
        our_id = None
        not_our_id = None  # used in the loop below, despite what pycharm thinks

        while start_row is None:

            # What is the next destination ID
            next_id = synapses[next_row, 1]

            # Is the next ID ours?
            # TODO: This can be speed up by instead having a bool array with 1 if neuron is in self.neuron_id and 0
            #       otherwise. That would prevent us from having to search in self.neuron_id
            # if next_id in self.neuron_id: -- OLD if statement, replaced with bool lookup below
            if self.neuron_id_on_node[next_id]:
                start_row = next_row
                our_id = next_id
                continue
            else:
                not_our_id = next_id

            # This loop just skips all synapses targeting not_our_id so we then can check next id
            while (next_row < num_syn_rows and
                   synapses[next_row, 1] == not_our_id):
                next_row += 1

            if next_row >= num_syn_rows:
                # No more synapses to get
                return None

        # Next find the last of the rows with this ID
        end_row = start_row

        while (end_row < num_syn_rows
               and synapses[end_row, 1] == our_id):
            end_row += 1

        return start_row, end_row

    ############################################################################

    # Processing the range(start_row,end_row) (ie, to end_row-1)

    def get_synapse_info(self, start_row, end_row):

        """
        Returns synapse information for synapses at start_row to end_row -1.

        Returns:
            (tuple):
                source_id_list (list of int): Presynaptic neuron ID
                dest_id (int): Destination neuron ID, obs single value, not a list
                dend_sections: Postsynaptic neuron sections
                sec_id (list of int): Section ID
                sec_x (list of float): Section X (between 0 and 1)
                synapse_type_id (list of int): Synapse type ID
                axon_distance (list of float):  Length of axon before synapse
                conductance (list of float): Conductance
                parameter_id (list of int): Synapse parameter ID

        """

        source_id_list = self.synapses[start_row:end_row, 0]
        source_id =self.synapses[start_row, 0]
        dest_id = self.synapses[start_row, 1]
        assert (self.synapses[start_row:end_row, 1] == dest_id).all()

        # Double check mapping (skip any synapses onto virtual neurons)
        assert self.is_virtual_neuron[dest_id] or self.pc.gid2cell(dest_id) == self.neurons[dest_id].icell, \
            f"GID mismatch: {self.pc.gid2cell(dest_id)} != {self.neurons[dest_id].icell}"

        synapse_type_id = self.synapses[start_row:end_row, 6]
        axon_distance = self.synapses[start_row:end_row, 7]  # Obs in micrometers

        sec_id = self.synapses[start_row:end_row, 9]
        sec_x = self.synapses[start_row:end_row, 10] / 1000.0  # Convert to number 0-1

        # Conductances are stored in pS (because we use INTs), NEURON wants it in microsiemens
        conductance = self.synapses[start_row:end_row, 11] * 1e-6
        parameter_id = self.synapses[start_row:end_row, 12]
        voxel_coords = self.synapses[start_row:end_row, 2:5]

        if not self.is_virtual_neuron[dest_id]:
            dend_sections = self.neurons[dest_id].map_id_to_compartment(sec_id)
            self.verify_synapse_placement(dend_sections, sec_x, dest_id, voxel_coords, source_id_list)
        else:
            dend_sections = None

        return source_id_list, dest_id, dend_sections, sec_id, sec_x, synapse_type_id, \
               axon_distance, conductance, parameter_id

    def connect_neuron_synapses(self, start_row, end_row):

        """ Connects the synapses present in the synapse matrix between start_row and end_row-1. """

        synapse_count = 0

        source_id_list, dest_id, dend_sections, sec_id, sec_x, synapse_type_id, \
        axon_distance, conductance, parameter_id = self.get_synapse_info(start_row=start_row, end_row=end_row)

        if dend_sections is None:
            # Target neuron was a virtual neuron, skip it
            return 0

        for (src_id, section, section_id, section_x, s_type_id, axon_dist, cond, p_id) \
                in zip(source_id_list, dend_sections, sec_id, sec_x, synapse_type_id,
                       axon_distance, conductance, parameter_id):

            try:
                self.add_synapse(cell_id_source=src_id,
                                 cell_id_dest=dest_id,
                                 dend_compartment=section,
                                 section_id=section_id,
                                 section_dist=section_x,
                                 synapse_type_id=s_type_id,
                                 axon_dist=axon_dist,
                                 conductance=cond,
                                 parameter_id=p_id)

                synapse_count += 1
            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr, is_error=True)
                import pdb
                pdb.set_trace()

        return synapse_count

    ############################################################################

    # OBS!! The src and dest lists can be different length
    #
    # src are all the gap junctions where the source compartment are
    # on the local worker.
    # dest are the gap junctions where the dest compartment are on the
    # local worker
    # The same GJ might appear in both src and dest lists, but at different rows

    def find_local_gap_junctions(self):

        """
        Locates the gap junctions present on the worker.

        Returns:
            (tuple):
                neuron_id: Neuron ID
                compartment: compartment
                seg_x: segment X
                gj_gid_src: gap junction source GID
                gj_gid_dest: gap junction destination GID
                cond: conductance
        """

        # If the gap junction matrix is too large to fit in memory then this will need to be optimised

        self.write_log("Finding node local gap junctions...")

        if self.gap_junctions.shape[0] == 0:
            return np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

        # gj_idx_a = np.where([x in self.neuron_id for x in self.gap_junctions[:, 0]])[0]
        # gj_idx_b = np.where([x in self.neuron_id for x in self.gap_junctions[:, 1]])[0]

        # We need to remove gap junctions where one or both of the neurons are virtual

        real_gj_idx = ~np.logical_or(self.is_virtual_neuron[self.gap_junctions[:, 0]],
                                     self.is_virtual_neuron[self.gap_junctions[:, 1]])

        if np.sum(real_gj_idx) == 0:
            return np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])

        gj_idx_a = np.where([self.neuron_id_on_node[x] for x in self.gap_junctions[real_gj_idx, 0]])[0]
        gj_idx_b = np.where([self.neuron_id_on_node[x] for x in self.gap_junctions[real_gj_idx, 1]])[0]

        gj_id_offset = 100 * self.num_neurons
        gj_gid_src_a = gj_id_offset + 2 * gj_idx_a
        gj_gid_src_b = gj_id_offset + 2 * gj_idx_b + 1

        gjg_id_dest_a = gj_id_offset + 2 * gj_idx_a + 1
        gj_gid_dest_b = gj_id_offset + 2 * gj_idx_b + 0

        neuron_id_a = self.gap_junctions[gj_idx_a, 0]
        neuron_id_b = self.gap_junctions[gj_idx_b, 1]

        seg_id_a = self.gap_junctions[gj_idx_a, 2]
        seg_id_b = self.gap_junctions[gj_idx_b, 3]

        try:
            compartment_a = [self.neurons[x].map_id_to_compartment([y])[0] for (x, y) in zip(neuron_id_a, seg_id_a)]
            compartment_b = [self.neurons[x].map_id_to_compartment([y])[0] for (x, y) in zip(neuron_id_b, seg_id_b)]
        except:
            import traceback
            self.write_log(traceback.format_exc(), is_error=True)
            import pdb
            pdb.set_trace()

        seg_xa = self.gap_junctions[gj_idx_a, 4] / 1000.0
        seg_xb = self.gap_junctions[gj_idx_b, 5] / 1000.0

        # Since we had ints we stored pS, but Neuron wants microsiemens
        cond_a = self.gap_junctions[gj_idx_a, 10] * 1e-6
        cond_b = self.gap_junctions[gj_idx_b, 10] * 1e-6

        # Merge the two lists together

        # gj_idx = np.concatenate([gj_idx_a, gj_idx_b])
        gj_gid_src = np.concatenate([gj_gid_src_a, gj_gid_src_b])
        gj_gid_dest = np.concatenate([gjg_id_dest_a, gj_gid_dest_b])
        neuron_id = np.concatenate([neuron_id_a, neuron_id_b])
        # seg_id = np.concatenate([seg_id_a, seg_id_b])
        compartment = np.concatenate([compartment_a, compartment_b])
        seg_x = np.concatenate([seg_xa, seg_xb])
        cond = np.concatenate([cond_a, cond_b])
        
        self.write_log(f"Found {len(neuron_id)} local gap junctions on node.")

        return neuron_id, compartment, seg_x, gj_gid_src, gj_gid_dest, cond

    ############################################################################

    # We can only do half the setup of the gap junction if it is split between two workers.

    def connect_network_gap_junctions_local(self):

        """ Setup gap junctions. Note that a gap junction is not complete until it has been setup by both workers
            that it is located on. """

        self.write_log("connect_network_gap_junctions_local")

        (neuron_id, compartment, seg_x, gj_gid_src, gj_gid_dest, cond) = self.find_local_gap_junctions()

        gap_junction_count = 0

        try:
            for nid, comp, s_x, gid_src, gid_dest, g \
                    in zip(neuron_id, compartment, seg_x, gj_gid_src, gj_gid_dest, cond):
                self.add_gap_junction(section=comp,
                                      section_dist=s_x,
                                      gid_source_gj=gid_src,
                                      gid_dest_gj=gid_dest,
                                      g_gap_junction=g,
                                      neuron_id=nid)

                gap_junction_count += 1

        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)
            import pdb
            pdb.set_trace()

        return gap_junction_count

    ############################################################################

    @staticmethod
    def get_synapse(channel_module, dend_compartment, section_dist):
        """ Helper method, returns channel_module(dend_compartment(section_dist)) """
        return channel_module(dend_compartment(section_dist))

    def add_synapse(self, cell_id_source, cell_id_dest, dend_compartment, section_id, section_dist, conductance,
                    parameter_id, synapse_type_id, axon_dist=None):

        """
        Add synapse

        Args:
            cell_id_source: Neuron ID of presynaptic neuron
            cell_id_dest
            dend_compartment: Dendrite compartment connecting to
            section_id
            section_dist: Section X
            conductance: Conductance of synapse
            parameter_id: Synapse parameter ID
            synapse_type_id: Synapse type ID
            axon_dist: Axon distance to presynaptic location

        """

        if self.is_virtual_neuron[cell_id_dest]:
            # The target neuron is a virtual neuron, do not add synapse
            return

        if conductance < 0:
            raise ValueError(f"Negative conductance found, this can be caused by specifying too large conductance values."
                             f"Remember that the synapse matrix is 32-bit int as pico siemens")

        # You can not locate a point process at endpoints (position 0.0 or 1.0) if it needs an ion
        if section_dist == 0.0:
            section_dist = 0.01
        if section_dist == 1.0:
            section_dist = 0.99

        (channel_module, par_data) = self.synapse_parameters[synapse_type_id]

        if channel_module is None:
            error_tag = ("channel_model_error", synapse_type_id)

            if error_tag not in self.print_error_once:
                error_message = (f"Warning: No channel module for {synapse_type_id} "
                                 f"between neuron {cell_id_source} and {cell_id_dest}, "
                                 f"did you miss specifying a mod file?")

                self.print_error_once[error_tag] = error_message
                self.write_log(error_message, is_error=True)

            # This channel was not implemented, skipping.
            return None

        syn = self.get_synapse(channel_module, dend_compartment, section_dist)

        weight_scale = 1

        if par_data is not None:
            # Picking one of the parameter sets stored in par_data
            par_id = parameter_id % len(par_data)

            par_set = par_data[par_id]

            for par in par_set:
                if par in ("expdata", "cond", "RxD"):
                    # expdata is not a parameter, and cond we take from synapse matrix
                    continue

                try:
                    # Can be value, or a tuple/list, if so second value is scale factor
                    # for SI -> natural units conversion
                    val = par_set[par]

                    # Do we need to convert from SI to natural units?
                    if type(val) == tuple or type(val) == list:
                        val_orig = val
                        val = val[0] * val[1]
                    else:
                        # If no list, we need to handle SI to natural units conversion automatically
                        val_orig = val
                        val = self.convert_to_natural_units(par, val)

                    if par in ["tau", "tauR"] and ((val < 0.01) or (10000 < val)):
                        self.write_log(f" !!! Warning: Converted from {val_orig} to {val} but expected "
                                       f"a value within [0.01, 10000) for neuron id {cell_id_source}. ", is_error=True)

                    setattr(syn, par, val)

                except:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)
                    import pdb
                    pdb.set_trace()

            if "RxD" in par_set and self.use_rxd_neuromodulation:
                species_name = par_set["RxD"]["species_name"]
                region = par_set["RxD"]["region"]

                if region in ("internal", "external"):
                    if section_id == -1:
                        region = f"soma_{region}"
                    else:
                        region = f"dend_{region}"

                weight_scale = par_set["RxD"].get("weight_scale", 1) * 1e-6  # (to compensate for 1e6 multiplication later)

                # If you have a RxD synapse it is good idea to set weight scale, especially
                # if your channel has valence 0, then cond variable is actually flux and needs to be in
                # number of molecules per second.
                flux_variable = par_set["RxD"]["flux_variable"]

                self.neurons[cell_id_dest].modulation.link_synapse(species_name=species_name,
                                                                   region=region,
                                                                   synapse=syn,
                                                                   flux_variable=flux_variable)

        if axon_dist is not None:
            # axon dist is in micrometer, want delay in ms
            synapse_delay = (1e3 * 1e-6 * axon_dist) / self.axon_speed + self.synapse_delay
        else:
            synapse_delay = self.synapse_delay

        # Prevent garbage collection in python
        if (cell_id_source, cell_id_dest) not in self.synapse_dict:
            self.synapse_dict[cell_id_source, cell_id_dest] = []

        nc = self.pc.gid_connect(cell_id_source, syn)
        nc.weight[0] = conductance * weight_scale
        nc.delay = synapse_delay
        nc.threshold = self.spike_threshold

        # This prevents garbage collection of syn and nc
        self.synapse_dict[cell_id_source, cell_id_dest].append((syn, nc, synapse_type_id, section_id))

        return syn

    ############################################################################

    # Add one gap junction to specific location

    def add_gap_junction(self,
                         section, section_dist,
                         gid_source_gj, gid_dest_gj,
                         g_gap_junction,
                         neuron_id):

        """
        Add gap junction.

        Args:
            section: Section to connect to
            section_dist: Section X
            gid_source_gj: GID of source gap junction
            gid_dest_gj: GID of destination gap junction
            g_gap_junction: Gap junction conductance
            neuron_id: ID of neuron, this is for book-keeping

        """

        # If neuron complains, make sure you have par_ggap.mod
        gj = h.gGapPar(section(section_dist))

        if neuron_id not in self.gap_junction_dict:
            self.gap_junction_dict[neuron_id] = [(gj, gid_source_gj, gid_dest_gj)]
        else:
            self.gap_junction_dict[neuron_id].append((gj, gid_source_gj, gid_dest_gj))

        # If you get a "NEURON: No source_var for target_var sid = 1301" error, then
        # make sure the neurons on both sides of the gap junction are included in the simulation

        # https://neuronsimulator.github.io/nrn/python/modelspec/programmatic/network/parcon.html?highlight=target_var
        # Update thanks to Lizhixin

        self.pc.target_var(gj, gj._ref_vgap, gid_dest_gj)
        self.pc.source_var(section(section_dist)._ref_v, gid_source_gj, sec=section)

        gj.g = g_gap_junction

    ############################################################################

    # Wilson 2007 - GABAergic inhibition in the neostriatum
    # 80% of synapses in Striatum are glutamatergic
    # Estimated 10000 glutamate and 2000 GABA synapses per MS,
    # 325 dopamine synapses per MS
    # Wilson 1996 - 10000 spines per MS = 10000 glutamatergic inputs

    # Ingham et al 1998, approx 1 glutamatergic synapse per 0.92 mum3
    # --> ~11000 glutamate synapses per MS
    # Kemp 1971 -- The synaptic organization of the caudate nucleus (85% glu)

    @staticmethod
    def get_external_input_synapse(channel_module, section, section_x):
        """ Helper method to return channel_module(section(section_x)) """
        return channel_module(section(section_x))

    def add_external_input(self, input_file=None):

        """ Adds external input from input_file to network. """

        if input_file is None:
            if self.input_file is None:
                self.write_log("No input file given, not adding external input!")
                return

            input_file = self.input_file

        self.write_log(f"Adding external (cortical, thalamic) input from {input_file}")

        self.input_data = h5py.File(input_file, 'r')

        for neuron_id, neuron in self.neurons.items():

            name = neuron.name

            if str(neuron_id) not in self.input_data["input"]:
                self.write_log(f"Warning - No input specified for {name}", is_error=True)
                continue

            for input_type in self.input_data["input"][str(neuron_id)]:

                self.external_stim[neuron_id, input_type] = []

                neuron_input = self.input_data["input"][str(neuron_id)][input_type]
                sections = self.neurons[neuron_id].map_id_to_compartment(neuron_input.attrs["section_id"])
                mod_file = SnuddaLoad.to_str(neuron_input.attrs["mod_file"])
                if "parameter_list" in neuron_input.attrs:
                    param_list = json.loads(neuron_input.attrs["parameter_list"], object_pairs_hook=OrderedDict)
                else:
                    param_list = None

                # TODO: Sanity check mod_file string
                eval_str = f"self.sim.neuron.h.{mod_file}"
                channel_module = eval(eval_str)

                rxd_species_name, rxd_flux_variable, rxd_region, rxd_weight_scale = self.get_rxd_external_input_parameters(neuron_input)

                for input_id, (section, section_id, section_x, param_id, n_spikes) \
                        in enumerate(zip(sections,
                                         neuron_input.attrs["section_id"],
                                         neuron_input.attrs["section_x"],
                                         neuron_input.attrs["parameter_id"],
                                         neuron_input["spikes"].attrs["num_spikes"])):

                    # We need to find cellID (int) from neuronID (string, eg. MSD1_3)

                    idx = input_id
                    spikes = neuron_input["spikes"][input_id, :n_spikes] * 1e3  # Neuron uses ms
                    assert (spikes >= 0).all(), f"Negative spike times for neuron {neuron_id} {input_type}"

                    # Creating NEURON VecStim and vector
                    # https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3125

                    try:
                        vs = h.VecStim()
                        v = h.Vector(spikes.size)
                        v.from_python(spikes)
                        vs.play(v)
                    except:
                        import traceback
                        tstr = traceback.format_exc()
                        self.write_log(tstr)

                        assert False, "!!! Make sure that vecevent.mod is included in nrnivmodl compilation"

                    # NEURON: You can not locate a point process at position 0 or 1 if it needs an ion
                    if section_x == 0.0:
                        section_x = 0.01
                    elif section_x == 1.0:
                        section_x = 0.99

                    syn = self.get_external_input_synapse(channel_module, section, section_x)
                    nc = h.NetCon(vs, syn)

                    nc.delay = 0.0
                    nc.weight[0] = neuron_input.attrs["conductance"][()] * rxd_weight_scale  # Neurons needs microsiemens
                    nc.threshold = 0.1

                    # Get the modifications of synapse parameters, specific to this synapse
                    if param_list is not None and len(param_list) > 0:
                        syn_params = param_list[param_id % len(param_list)]  # No longer need to take ["synapse"], only that part saved in hdf5

                        for par in syn_params:

                            if par == "RxD":
                                # RxD specific information
                                continue

                            if par == "expdata":
                                # Not a parameter
                                continue

                            if par == "cond":
                                # Ignoring cond value specified for synapse, using the
                                # one specified in the input information instead
                                continue

                            par_value = self.convert_to_natural_units(par, syn_params[par])

                            if par in ["tau", "tauR"]:
                                assert 0.01 <= par_value < 10000, \
                                    (f"Converting {self.neurons[neuron_id].name} {par}={syn_params[par]} "
                                     f"we get {par_value}, "
                                     f"but expected >= 0.01 and < 10000")
                            # print(f"Setting {par} to {par_value}.")
                            setattr(syn, par, par_value)

                    if rxd_species_name is not None:
                        if section_id == -1:
                            region_name = f"soma_{rxd_region}"
                        else:
                            region_name = f"dend_{rxd_region}"

                        self.neurons[neuron_id].modulation.link_synapse(species_name=rxd_species_name,
                                                                        region=region_name,
                                                                        synapse=syn,
                                                                        flux_variable=rxd_flux_variable)

                    # Need to save references, otherwise they will be freed
                    self.external_stim[neuron_id, input_type].append((v, vs, nc, syn, spikes, section_id, section_x))

    def get_rxd_external_input_parameters(self, neuron_input):

        if "RxD" in neuron_input.attrs.keys() and self.use_rxd_neuromodulation:

            rxd_dict = json.loads(neuron_input.attrs["RxD"])

            species_name = rxd_dict.get("species_name")
            flux_variable = rxd_dict.get("flux_variable")
            region = rxd_dict.get("region")
            weight_scale = rxd_dict.get("weight_scale", 1.0)

        else:
            species_name = None
            flux_variable = None
            region = None
            weight_scale = 1e6  # we need to do pS -> micro siemens

        return species_name, flux_variable, region, weight_scale

    ############################################################################

    def set_resting_voltage(self, neuron_id, rest_volt=None):

        """
        Sets resting voltage for neuron

        Args:
            neuron_id: Neuron ID (either int, or list of int)
            rest_volt: Resting voltage (either None = read from parameter files, float, or list of floats)
                       in SI units (volt), gets converted to mV before passing to NEURON

        """

        if type(neuron_id) != list:
            neuron_id_list = [neuron_id]
        else:
            neuron_id_list = neuron_id

        if rest_volt is None:
            rest_volt_list = [None] * len(neuron_id_list)
        elif type(rest_volt) != list:
            rest_volt_list = [rest_volt]
        else:
            rest_volt_list = rest_volt

        for neuron_id, rest_volt in zip(neuron_id_list, rest_volt_list):

            if neuron_id not in self.neuron_id:
                # This neuron is not on this worker, continue
                continue

            if rest_volt is None:
                # If no resting voltage is given, extract it from parameters
                # Note that the file has NEURON v_init from bluepyopt, in natural units, so convert to SI internally
                rest_volt = [x for x in self.neurons[neuron_id].parameters
                             if "param_name" in x and x["param_name"] == "v_init"][0]["value"] * 1e-3

            self.write_log(f"Neuron {self.neurons[neuron_id].name} ({neuron_id}) resting voltage = {rest_volt * 1e3}")

            soma = [x for x in self.neurons[neuron_id].icell.soma]
            axon = [x for x in self.neurons[neuron_id].icell.axon]
            # If there are no dendrites this will crash NEURON icell.dend has length == 1, but that is WRONG
            dend = [x for x in self.neurons[neuron_id].icell.dend]

            cell = soma + axon + dend

            for sec in cell:
                for seg in sec.allseg():
                    seg.v = rest_volt * 1e3

    ############################################################################

    def centre_neurons(self, side_len=None, neuron_id=None):

        """
        Returns neurons that are in the centre of the volume. Useful when simulating for example a cube,
        and we want to avoid neurons with border artifacts.
        """

        if neuron_id is None:
            neuron_id = self.neuron_id

        if side_len is None:
            return neuron_id

        c_id = []

        positions = self.network_info["neuron_positions"]

        centre_pos = np.min(positions, axis=0)

        for nid in neuron_id:

            pos = positions[nid, :]

            if (abs(pos[0] - centre_pos[0]) <= side_len
                    and abs(pos[1] - centre_pos[1]) <= side_len
                    and abs(pos[2] - centre_pos[2]) <= side_len):
                c_id.append(nid)

        self.write_log(f"Centering: Keeping {len(c_id)}/{len(neuron_id)}")

        return c_id

    ############################################################################

    def add_recording_of_type(self, neuron_type, num_neurons=None):

        """ Adds somatic voltage recording to num_neurons of neuron_type. """

        cell_id = self.snudda_loader.get_neuron_id_of_type(neuron_type=neuron_type, num_neurons=num_neurons)

        self.add_volt_recording_all(cell_id)

    ############################################################################

    def add_voltage_clamp(self, cell_id, voltage, duration, res=1e-9, save_i_flag=False):

        """
        Adds voltage clamp.

        Args:
            cell_id : Neuron ID
            voltage : Voltage
            duration : Duration
            res : Resistance
            save_i_flag : Should current be saved
        """

        if type(cell_id) not in [list, np.ndarray]:
            cell_id = [cell_id]

        if type(voltage) not in [list, np.ndarray]:
            voltage = [voltage for x in cell_id]

        if type(duration) not in [list, np.ndarray]:
            duration = [duration for x in cell_id]

        if type(res) not in [list, np.ndarray]:
            res = [res for x in cell_id]

        if save_i_flag and (len(self.t_save) == 0 or self.t_save is None):
            self.t_save = self.sim.neuron.h.Vector()
            self.t_save.record(self.sim.neuron.h._ref_t)

        for cID, v, rs, dur in zip(cell_id, voltage, res, duration):

            try:
                if not (cID in self.neuron_id):
                    # Not in the list of neuronID on the worker, skip it
                    continue
            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr, is_error=True)
                import pdb
                pdb.set_trace()

            self.write_log(f"Adding voltage clamp to {cID}")
            s = self.neurons[cID].icell.soma[0]
            vc = neuron.h.SEClamp(s(0.5))
            vc.rs = rs
            vc.amp1 = v * 1e3
            vc.dur1 = dur * 1e3

            self.write_log(f"Resistance: {rs}, voltage: {vc.amp1}mV")

            self.v_clamp_list.append(vc)

            if save_i_flag:
                cur = self.sim.neuron.h.Vector()
                cur.record(vc._ref_i)
                self.i_save.append(cur)
                self.i_key.append(cID)

    ############################################################################

    def add_volt_recording_all(self, cell_id=None, centre_only_flag=True, section_x=None):

        # centre_only_flag has priority over section_x parameter

        if cell_id is None:
            cell_id = self.neuron_id

        if isinstance(cell_id, (int, np.integer)):
            cell_id = [cell_id]

        local_cell_id = np.intersect1d(cell_id, list(self.neurons.keys()))

        for cid in local_cell_id:
            sec_id = [-1]  # Soma is sec_id = -1
            sec_x = [0.5]

            for sid, sec in enumerate(self.neurons[cid].icell.dend):

                assert int(sec.name().split('[')[-1].strip(']')) == sid, \
                    f"Internal error, assumed {sid} was section id of {sec.name()}"

                if centre_only_flag:
                    sec_id.append(sid)
                    sec_x.append(0.5)
                elif section_x:
                    sec_id.append(sid)
                    sec_x.append(section_x)
                else:
                    for seg in sec.allseg():
                        sec_id.append(sid)
                        sec_x.append(seg.x)

            self.add_volt_recording(cell_id=cid, sec_id=sec_id, sec_x=sec_x)

    def add_volt_recording_soma(self, cell_id=None):

        # Adds somatic voltage recording to neurons on worker (but not virtual neurons)

        if cell_id is None:
            cell_id = self.neuron_id

        for cid in cell_id:

            if cid in self.neurons.keys() and not self.is_virtual_neuron[cid]:
                self.add_volt_recording(cid, [-1], [0.5])  # Soma is sec_id = -1

    def add_volt_recording(self, cell_id: int, sec_id, sec_x):

        # cell_id cant be a list, but sec_id and sec_x must be iterable

        sections = self.neurons[cell_id].map_id_to_compartment(sec_id)

        self.pc.threshold(cell_id, self.spike_threshold)  # TODO: Set individual spike thresholds based on parameters
        for s, sx, sid in zip(sections, sec_x, sec_id):
            v = self.sim.neuron.h.Vector()
            v.record(getattr(s(sx), '_ref_v'))

            # From the Snudda synapse matrix. sec_id -1 is soma, sec_id >= 0 is dendrite, sec_id <= -2 is axon
            self.record.register_compartment_data(neuron_id=cell_id,
                                                  data_type="voltage", data=v,
                                                  sec_id=sid, sec_x=sx)

        if self.record.time is None:
            t_save = self.sim.neuron.h.Vector()
            t_save.record(self.sim.neuron.h._ref_t)
            self.record.register_time(time=t_save)

    def add_synapse_current_recording_all(self, dest_id=None, max_synapses=500):

        """
            Record all synaptic currents to neuron dest_id. If dest_id is None then all synaptic currents are recorded.
            The total number of currents recorded is capped at max_synapses.

            Args:
                dest_id (int or list of ints) : Postsynaptic neuron ID
                max_synapses (int) : Maximum number of synapses to record from (default: 500)
        """

        syn_ctr = 0

        if dest_id is None:
            self.write_log("Warning, recording ALL synapse currents.", force_print=True)

            for sid, did in self.synapse_dict.keys():
                if syn_ctr > max_synapses:
                    break

                syn_ctr += self.add_synapse_current_recording(sid, did)

        else:
            if isinstance(dest_id, (int, np.integer)):
                dest_id = [dest_id]

            for d_id in dest_id:
                for sid, did in self.synapse_dict.keys():
                    if did == d_id and syn_ctr < max_synapses:
                        syn_ctr += self.add_synapse_current_recording(sid, did)

        if syn_ctr > max_synapses:
            self.write_log(f"Warning: Not recording all synapse currents requested, capped at max_synapses={max_synapses}",
                           force_print=True)

    def add_density_mechanism_recording(self, density_mechanism: str, variable: str,
                                        neuron_id: int, sec_id: int, sec_x: float):

        """ Record density mechanism:

        Args:
            density_mechanism (str): Name of density mechanism, e.g. "pas"
            variable (str): Variable name to record, e.g. "i" (will then read _ref_i)
            neuron_id (int) : Id of neuron
            sec_id (int) : Section id of compartment
            sec_x (float) : Section x
        """

        if neuron_id not in self.neuron_id:
            # The neuron is not on this worker
            return

        segment = self.neurons[neuron_id].map_id_to_compartment(sec_id)(sec_x)
        mech = getattr(segment, density_mechanism)
        var = getattr(mech, f"_ref_{variable}")
        data = self.sim.neuron.h.Vector().record(var)

        self.record.register_compartment_data(data_type=f"{density_mechanism}.{variable}",
                                              neuron_id=neuron_id,
                                              data=data,
                                              sec_id=sec_id,
                                              sec_x=sec_x)

        if self.record.time is None:
            t_save = self.sim.neuron.h.Vector()
            t_save.record(self.sim.neuron.h._ref_t)
            self.record.register_time(time=t_save)

    def add_membrane_recording(self, variable, neuron_id, sec_id, sec_x):

        if neuron_id not in self.neuron_id:
            # The neuron is not on this worker
            return

        segment = self.neurons[neuron_id].map_id_to_compartment(sec_id)(sec_x)
        var = getattr(segment, f"_ref_{variable}")
        data = self.sim.neuron.h.Vector().record(var)

        self.record.register_compartment_data(data_type=f"membrane.{variable}",
                                              neuron_id=neuron_id,
                                              data=data,
                                              sec_id=sec_id,
                                              sec_x=sec_x)

    def get_internal_synapse_point_process(self, source_id, dest_id, synapse_type=None):

        synapse_list = self.synapse_dict.get((source_id, dest_id), [])

        channel_model_id = None

        if synapse_type is not None:
            pre_type = self.network_info["neurons"][source_id]["type"]
            post_type = self.network_info["neurons"][dest_id]["type"]

            if (pre_type, post_type) in self.network_info["connectivity_distributions"]:
                channel_model_id = self.network_info["connectivity_distributions"][(pre_type, post_type)][synapse_type]["channel_model_id"]

        s_list = [synapse_info for synapse_info in synapse_list
                  if channel_model_id is None or channel_model_id == synapse_info[2]]

        return s_list, pre_type, post_type

    def add_synapse_variable_recording(self, source_id, dest_id, variable, synapse_type=None):

        synapse_list, pre_type, post_type = self.get_internal_synapse_point_process(source_id=source_id,
                                                                                    dest_id=dest_id,
                                                                                    synapse_type=synapse_type)

        return self.add_point_process_variable_recording(point_process_list=synapse_list,
                                                         variable=variable,
                                                         post_synaptic_id=dest_id,
                                                         pre_synaptic_id=source_id,
                                                         name=f"{pre_type}_{post_type}_{synapse_type}")

    def add_gap_junction_current_recording(self, neuron_id, gj_idx=None):

        gj_list = self.gap_junction_dict.get(neuron_id, [])

        if gj_idx is not None:
            gj_list = [gj_list[gj_idx]]

        # syn, nc, synapse_type_id, sec_id
        point_process = [(x[0], x[0], 3, -100) for x in gj_list]

        return self.add_point_process_variable_recording(point_process_list=point_process,
                                                         variable="i",
                                                         post_synaptic_id=neuron_id,
                                                         pre_synaptic_id=-1,  # We need to add bookkeeping to track destination
                                                         name=f"gj_currents_{neuron_id}",
                                                         use_netcon_weight=False)

    def get_external_synapse_point_process(self, neuron_id, input_type):

        """ Returns point process of the external synapses on the neuron """

        external_input = self.external_stim[neuron_id, input_type]

        syn_list = [(x[3], x[2], x[5], x[6]) for x in external_input]   # sym and nc

        return syn_list

    def add_external_input_variable_recording(self, neuron_id, input_type, variable, name=""):

        # TODO: Verify this works...

        syn_list = self.get_external_synapse_point_process(neuron_id=neuron_id, input_type=input_type)
        syn_ctr = 0

        for syn, nc, sec_id, sec_x in syn_list:
            data = self.sim.neuron.h.Vector()
            data.record(getattr(syn, f"_ref_{variable}"))
            seg = syn.get_segment()

            # They are close, but not identical...
            # assert sec_x == seg.x, f"Internal error, {sec_x = } should be same as stored in {seg.x = }"

            self.record.register_synapse_data(neuron_id=neuron_id,
                                              data_type=f"{name}{'.' if len(name) > 0 else ''}{input_type}. {variable}", data=data,
                                              synapse_type=-1,  # Check what the real number is
                                              presynaptic_id=-1,  # External input
                                              sec_id=sec_id,
                                              sec_x=seg.x,
                                              cond=nc.weight[0])
            syn_ctr += 1

        return syn_ctr

    def add_point_process_variable_recording(self, point_process_list, variable,
                                             post_synaptic_id, pre_synaptic_id=-1,
                                             name="", use_netcon_weight=True):

        if not isinstance(point_process_list, list):
            point_process_list = list(point_process_list)

        syn_ctr = 0

        for syn, nc, synapse_type_id, sec_id in point_process_list:
            data = self.sim.neuron.h.Vector()
            data.record(getattr(syn, f"_ref_{variable}"))
            seg = syn.get_segment()

            if use_netcon_weight:
                cond = nc.weight[0]  # netcon object for synapses
            else:
                cond = nc.g  # For gap junctions nc is a gap junction object, and g conductance

            self.record.register_synapse_data(neuron_id=post_synaptic_id,
                                              data_type=f"{name}{'.' if len(name) > 0 else ''}{variable}", data=data,
                                              synapse_type=synapse_type_id,
                                              presynaptic_id=pre_synaptic_id,
                                              sec_id=sec_id,
                                              sec_x=seg.x,
                                              cond=cond)
            syn_ctr += 1

        return syn_ctr

    def add_synapse_current_recording(self, source_id, dest_id):

        assert (source_id, dest_id) in self.synapse_dict, f"No synapse between {source_id} and {dest_id}"

        synapse_info_list = self.synapse_dict[source_id, dest_id]
        syn_ctr = 0

        for syn, nc, synapse_type_id, sec_id in synapse_info_list:
            data = self.sim.neuron.h.Vector()
            data.record(syn._ref_i)
            seg = syn.get_segment()

            self.record.register_synapse_data(neuron_id=dest_id, data_type="synaptic_current", data=data,
                                              synapse_type=synapse_type_id,
                                              presynaptic_id=source_id,
                                              sec_id=sec_id,
                                              sec_x=seg.x,
                                              cond=nc.weight[0])
            syn_ctr += 1

        return syn_ctr

    def add_rxd_extracellular_concentration_recording(self, species: str, volume_id, xyz):

        if volume_id not in self.extracellular_regions:
            raise KeyError(f"{volume_id =} missing, available keys {list(self.extracellular_regions.keys())}")

        x, y, z = 1e6 * np.array(xyz)  # unit conversion

        if species not in self.extracellular_regions[volume_id].species:
            raise KeyError(f"{species =} not available in {volume_id}. Available spcies are: "
                           f"{list(self.extracellular_regions[volume_id].species.keys())}")

        spec = self.extracellular_regions[volume_id].species[species]["ecs"]
        ecs = self.extracellular_regions[volume_id].compartments["ecs"]
        conc_ref = spec[ecs].node_by_location(x, y, z)._ref_value
        vector = self.sim.neuron.h.Vector()
        vector.record(conc_ref)

        # Get the ijk index:
        e = self.extracellular_regions[volume_id].compartments["ecs"]

        if e._xlo <= x <= e._xhi and e._ylo <= y <= e._yhi and e._zlo <= z <= e._zhi:
            i = int((x - e._xlo) / e._dx[0])
            j = int((y - e._ylo) / e._dx[1])
            k = int((z - e._zlo) / e._dx[2])

        index_ijk = (i, j, k)

        # Convert back from RxD millimolar -> molar
        self.record.add_unit(data_type=species, target_unit="molar", conversion_factor=1e-3)

        self.record.register_extracellular_data(data=vector, data_type=species, index_ijk=index_ijk,
                                                ecs=self.extracellular_regions, region_name=volume_id)

        if self.record.time is None:
            t_save = self.sim.neuron.h.Vector()
            t_save.record(self.sim.neuron.h._ref_t)
            self.record.register_time(time=t_save)

        pass

    def add_rxd_concentration_recording(self, species: str, neuron_id: int, region, sec_id, sec_x):

        if not self.use_rxd_neuromodulation:
            self.write_log(f"add_rxd_concentration_recording:  not enabled, ignoring recording of {species} in neuron {neuron_id}")
            return

        if sec_id == -1:
            sec_type = "soma"
            neuron_sec_id = 0
        elif sec_id >= 0:
            sec_type = "dend"
            neuron_sec_id = sec_id
        else:
            sec_type = "axon"
            neuron_sec_id = 0
            self.write_log("Axon recordings currently not fully supported (using sec_id=0")

        if self.neurons[neuron_id].modulation is None:
            raise ValueError(f"No modulation specified for neuron {self.neurons[neuron_id].name} ({neuron_id})")

        try:
            segment = getattr(self.neurons[neuron_id].icell, sec_type)[neuron_sec_id](sec_x)
        except:
            import traceback
            self.write_log(traceback.format_exc())
            import pdb
            pdb.set_trace()

        conc_ref = self.neurons[neuron_id].modulation.species[species][region].nodes(segment)._ref_concentration

        vector = self.sim.neuron.h.Vector()
        vector.record(conc_ref)

        # Convert back from RxD millimolar -> molar
        self.record.add_unit(data_type=species, target_unit="molar", conversion_factor=1e-3)

        self.record.register_compartment_data(neuron_id=neuron_id,
                                              data_type=species,
                                              data=vector,
                                              sec_id=sec_id, sec_x=sec_x)

        if self.record.time is None:
            t_save = self.sim.neuron.h.Vector()
            t_save.record(self.sim.neuron.h._ref_t)
            self.record.register_time(time=t_save)

    def add_rxd_internal_concentration_recording_all(self, species, neuron_id):

        if neuron_id not in self.neuron_id:
            return

        if not self.use_rxd_neuromodulation:
            self.write_log(f"add_rxd_internal_concentration_recording_all:  not enabled, ignoring recording of {species} in neuron {neuron_id}")
            return

        # Add soma
        self.add_rxd_concentration_recording(species, neuron_id, "soma_internal", -1, 0.5)

        for sid, sec in enumerate(self.neurons[neuron_id].icell.dend):
            assert int(sec.name().split('[')[-1].strip(']')) == sid, \
                f"Internal error, assumed {sid} was section id of {sec.name()}"
            self.add_rxd_concentration_recording(species, neuron_id, "dend_internal", sid, 0.5)

    def add_rxd_internal_concentration_recording_all_species(self, neuron_id, include_dendrites=True):

        if self.verbose:
            self.write_log(f"Recording all RxD species from neurons: {neuron_id}")

        if isinstance(neuron_id, (list, np.ndarray)):
            for nid in neuron_id:
                self.add_rxd_internal_concentration_recording_all_species(neuron_id=nid)
            return

        if neuron_id not in self.neuron_id:
            return

        if not self.use_rxd_neuromodulation:
            self.write_log(f"add_rxd_internal_concentration_recording_all_species:  not enabled, ignoring recording of neuron {neuron_id}")
            return

        for species in self.neurons[neuron_id].modulation.species.keys():

            # Add soma
            self.add_rxd_concentration_recording(species, neuron_id, "soma_internal", -1, 0.5)

            if include_dendrites:
                for sid, sec in enumerate(self.neurons[neuron_id].icell.dend):
                    assert int(sec.name().split('[')[-1].strip(']')) == sid, \
                        f"Internal error, assumed {sid} was section id of {sec.name()}"
                    self.add_rxd_concentration_recording(species, neuron_id, "dend_internal", sid, 0.5)

    def add_bath_application(self, species_name, concentration, time, neuron_id=None, interpolate=True):

        if neuron_id is None:
            neuron_id = self.snudda_loader.get_neuron_id(include_virtual=False)

        conc_vect = self.sim.neuron.h.Vector(concentration * 1e3)  # SI to millimolar
        t_vect = self.sim.neuron.h.Vector(time * 1e3)  # s -> ms

        if self.verbose:
            self.write_log(f"Bath application t={time*1e3}ms, conc={concentration*1e3}")

        if species_name is self.bath_application:
            raise KeyError(f"Bath application already applied for {species_name}")

        self.bath_application[species_name] = (time, concentration, t_vect, conc_vect)

        for nid in neuron_id:
            if nid in self.neurons:
                n = self.neurons[nid]

                if n.modulation is not None:
                    n.modulation.concentration_from_vector(species_name=species_name,
                                                           concentration_vector=conc_vect,
                                                           time_vector=t_vect,
                                                           interpolate=interpolate)

    ############################################################################

    def sanity_check_play_vectors(self, sim_end_time):
        # TODO: Add additional checks that all bath application play vectors are long enough

        for species_name, bath_data in self.bath_application.items():
            bath_max_time = np.max(bath_data[0])

            if sim_end_time > bath_max_time:
                raise ValueError(f"Simulation duration {sim_end_time} is "
                                 f"longer than time vector for bath application of {species_name}")

    def post_initialisation_modifications(self):

        if len(self.post_init_mods) == 0:
            return

        self.write_log(f"Applying post_initialisation_modifications")
        # This allows us to modify ion channel conductance on the fly before runnigg simulation
        # Can be useful to e.g. increase KIR channel conductance

        for n in self.neurons.values():
            try:
                if n.type in self.post_init_mods:
                    for ion_channel, channel_mod_factor in self.post_init_mods[n.type].items():
                        for sec in n.icell.all:
                            for seg in sec:
                                channel = getattr(seg, ion_channel, None)
                                if channel is not None:
                                    setattr(channel, "gbar", getattr(channel, "gbar") * channel_mod_factor)
            except:
                import traceback
                self.write_log(traceback.format_exc(), is_error=True)
                import pdb
                pdb.set_trace()

    def run(self, t=None, hold_v=None):

        """ Run simulation. """

        start_time = timeit.default_timer()

        if self.is_virtual_neuron.all():
            self.write_log("ALL YOUR NEURONS IN THE SIMULATION ARE VIRTUAL")

        if t is None:
            if self.sim_info is not None and "time" in self.sim_info:
                t = self.sim_info["time"] * 1e3  # convert to ms for NEURON
            else:
                t = 1000.0

        self.sanity_check_play_vectors(sim_end_time=t*1e-3)

        if hold_v is None:
            if self.sim_info is not None and "hold_voltage" in self.sim_info:
                hold_v = self.sim_info["hold_voltage"]

        # If we want to use a non-default initialisation voltage, we need to
        # explicitly set: h.v_init
        # self.sim.neuron.h.v_init = -78
        # self.sim.neuron.h.finitialize(-78)

        # Important, when simulating gap junctions, neurons on both sides of the gap junction
        # must be simulated.

        if hold_v is None:
            self.sim.neuron.h.finitialize()
        else:
            self.write_log(f"User override for holding voltage: {hold_v * 1e3} mV")
            self.sim.neuron.h.finitialize(hold_v * 1e3)

        # Asked on neuron, check answer:
        # https://www.neuron.yale.edu/phpBB/viewtopic.php?f=2&t=4161&p=18021

        self.setup_print_sim_time(t)

        # Make sure all processes are synchronised
        self.pc.barrier()

        self.write_log(f"Running simulation for {t / 1000} s", force_print=True)

        self.sim.run(t, dt=self.sim_dt)
        self.pc.barrier()
        self.write_log("Simulation done.")

        end_time = timeit.default_timer()
        self.write_log(f"Simulation run time: {end_time - start_time:.1f} s", force_print=True)

    ############################################################################
    # export_to_core_neuron contributed by Zhixin, email: 2001210624@pku.edu.cn

    def export_to_core_neuron(self):

        """ Export network for faster optimised execution using core neuron."""

        core_data_path = os.path.join(self.network_path, "CoreNeuronModule")
        rank = self.pc.id()
        self.pc.barrier()

        self.write_log(f"Exporting core neuron module to {core_data_path}")
        if rank == 0 and core_data_path is not None and not os.path.isdir(core_data_path):
            os.mkdir(core_data_path)

        self.pc.barrier()
        h.secondorder = 1
        h.cvode.cache_efficient(1)
        self.pc.set_maxstep(10)
        h.stdinit()
        self.pc.nrnbbcore_write(core_data_path)
        self.pc.barrier()

    ############################################################################

    def verify_synapse_placement(self, sec_list, sec_x_list, dest_id, voxel_coords, source_id_list=None):

        """ Verifies synapse placement.

            We want to check that voxel coords transformed to local coordinate system
            of neuron matches with where neuron places the synapse

            Args:

                sec_list (list) : list of sections
                sec_x_list (list) : list of X values 0 to 1.0
                dest_id (int) : ID of the neuron receiving synapse (one value!)
                voxel_coords : voxel that the synapse is in
                source_id_list (np.array) : ID of sending neurons, not used, for debug

        """

        simulation_origo = self.network_info["simulation_origo"]
        voxel_size = self.network_info["voxel_size"]
        neuron_position = self.network_info["neurons"][dest_id]["position"]
        neuron_rotation = self.network_info["neurons"][dest_id]["rotation"]

        # Transform voxel coordinates to local neuron coordinates to match neuron
        # synapse_pos = (voxel_size * voxel_coords + simulation_origo - neuron_position) * 1e6

        # UPDATE 2025-01-24: Now the neuron NEURON coordinates should use the simulation coordinate system
        synapse_pos = (voxel_size * voxel_coords + simulation_origo) * 1e6

        syn_pos_nrn = np.zeros((len(sec_list), 3))
        old_sec = None
        norm_arc_dist = None

        for i, (sec, sec_x) in enumerate(zip(sec_list, sec_x_list)):

            # If statement is just so we don't recalculate the norm_arc_dist every time
            if old_sec is None or not sec.same(old_sec):
                num_points = int(h.n3d(sec=sec))
                arc_dist = np.array([sec.arc3d(x) for x in range(0, num_points)])
                norm_arc_dist = arc_dist / arc_dist[-1]
                old_sec = sec

            # Find the closest point
            closest_idx = np.argmin(np.abs(norm_arc_dist - sec_x))

            syn_pos_nrn[i, 0] = h.x3d(closest_idx, sec=sec)
            syn_pos_nrn[i, 1] = h.y3d(closest_idx, sec=sec)
            syn_pos_nrn[i, 2] = h.z3d(closest_idx, sec=sec)

        # We need to rotate the neuron to match the big simulation
        # !!! OBS, this assumes that soma is in 0,0,0 local coordinates
        # syn_pos_nrn_rot = np.transpose(np.matmul(neuron_rotation,
        #                                        np.transpose(syn_pos_nrn)))

        # UPDATED 2025-01-24: NEURON neuron coordinates are now updated to use the simulation coordinate frame
        syn_pos_nrn_rot = syn_pos_nrn

        syn_mismatch = np.sqrt(np.sum((syn_pos_nrn_rot - synapse_pos) ** 2, axis=1))

        bad_threshold = 22
        num_bad = np.sum(syn_mismatch > bad_threshold)

        if num_bad > 0:
            bad_idx = np.where(syn_mismatch > bad_threshold)[0]
            bad_sec_len = [sec_list[bi].L for bi in bad_idx]

            # If this happens, check that Neuron does not warn for removing sections
            # due to having only one point
            self.write_log(f"!!! Found {num_bad} synapses on "
                           f"{self.network_info['neurons'][dest_id]['name']} ({dest_id}) "
                           f" that are farther than {bad_threshold} mum away "
                           f" (out of {len(syn_mismatch)} synapses)"
                           f" Max found was {np.max(syn_mismatch):.0f} mum from expected location."
                           f" morphology: {self.network_info['neurons'][dest_id]['morphology']}\n"
                           f" Check that soma is centered at (0,0,0). Also check that the first dendritic"
                           f" compartment of each dendrite is not too far away from the soma, then NEURON "
                           f" adds an extra connecting compartment which messes up section IDs.",
                           is_error=True)
            self.write_log(f"Length of sections with bad synapses: {bad_sec_len}", is_error=True)

            for bi in bad_idx:
                closest_sec, closest_sec_x, min_dist = self.find_closest_point_on_neuron(neuron_id=dest_id,
                                                                                         synapse_xyz=synapse_pos[bi, :])

                self.write_log(f"Neuron id: {dest_id} Bad synapse {bi} on {sec_list[bi]} {sec_x_list[bi]}, "
                      f"closer match at {closest_sec} {closest_sec_x}, dist: {min_dist} (source: {source_id_list[bi]})")

            ### DEBUG PLOT!!!

            if True:
                import matplotlib.pyplot as plt
                # plt.figure()
                #
                # soma_dist = np.sqrt(np.sum(synapse_pos ** 2, axis=1))
                # plt.scatter(soma_dist * 1e6, syn_mismatch)
                # # plt.ion()
                # plt.show()
                # plt.title(self.network_info["neurons"][dest_id]["name"])

                fig = plt.figure()
                ax = fig.add_subplot(projection='3d')

                syn_size = np.full((synapse_pos.shape[0],), 10)
                syn_size[bad_idx] = 50

                col_list = ['orange' for x in range(synapse_pos.shape[0])]
                for b in bad_idx:
                    col_list[b] = 'red'

                ax.scatter(synapse_pos[:, 0],
                           synapse_pos[:, 1],
                           synapse_pos[:, 2], color=col_list, s=syn_size)
                ax.scatter(syn_pos_nrn_rot[:, 0],
                           syn_pos_nrn_rot[:, 1],
                           syn_pos_nrn_rot[:, 2], color="black", s=syn_size*2)

                if True:
                    # Draw neuron
                    # all_sec = [x for x in neuron.h.allsec() if "axon" not in str(x)]
                    all_sec = [s for s in self.neurons[dest_id].icell.dend]
                    for sec in all_sec:
                        x = np.array([h.x3d(i, sec=sec) for i in range(sec.n3d())])
                        y = np.array([h.y3d(i, sec=sec) for i in range(sec.n3d())])
                        z = np.array([h.z3d(i, sec=sec) for i in range(sec.n3d())])

                        xyz = np.matmul(neuron_rotation, np.array([x, y, z]))

                        ax.plot(xyz[0, :].T, xyz[1, :].T, xyz[2, :].T, 'b-')

                plt.pause(0.001)
                plt.savefig("DEBUG-plot-bad-synapse-placement.png", dpi=600)
                plt.ion()
                plt.show()
                plt.pause(0.001)

                import pdb
                pdb.set_trace()

    ############################################################################

    def find_closest_point_on_neuron(self, neuron_id, synapse_xyz):

        min_dist = np.inf
        closest_sec = None
        closest_sec_x = None

        neuron_rotation = self.network_info["neurons"][neuron_id]["rotation"]

        for sec in self.neurons[neuron_id].icell.dend:
            # Extract all coordinates
            n_points = int(h.n3d(sec=sec))
            xyz = np.zeros((n_points, 3))
            for i in range(0, n_points):
                xyz[i, 0] = h.x3d(i, sec=sec)
                xyz[i, 1] = h.y3d(i, sec=sec)
                xyz[i, 2] = h.z3d(i, sec=sec)

            xyz = np.transpose(np.matmul(neuron_rotation, np.transpose(xyz)))
            d = np.linalg.norm(xyz-synapse_xyz, axis=1)
            d_min_idx = np.argmin(d)

            if d[d_min_idx] < min_dist:
                min_dist = d[d_min_idx]
                closest_sec = sec
                closest_sec_x = sec.arc3d(d_min_idx) / sec.arc3d(n_points-1)

        return closest_sec, closest_sec_x, min_dist

    ############################################################################

    def write_output(self):

        self.record.write()

    ##############################################################################

    def write_log(self, text, flush=True, is_error=False, force_print=False):

        """
        Writes to log file. Use setup_log first. Text is only written to screen if self.verbose=True,
        or is_error = True, or force_print = True.

        test (str) : Text to write
        flush (bool) : Should all writes be flushed to disk directly?
        is_error (bool) : Is this an error, always written.
        force_print (bool) : Force printing, even if self.verbose=False.
        """

        if isinstance(self.log_file, str):
            self.log_file = open(self.log_file, "w")

        if self.log_file is not None:
            self.log_file.write(f"{text}\n")
            if flush:
                self.log_file.flush()

        if self.verbose or is_error or force_print:
            print(text)

    ############################################################################

    def create_dir(self, dir_name):

        """ Creates dir_name if needed. """
        if int(self.pc.id()) == 0:
            if not os.path.isdir(dir_name):
                self.write_log(f"Creating {dir_name} (on master node 0)")
                os.makedirs(dir_name)
        else:
            while not os.path.isdir(dir_name):
                self.write_log(f"Waiting 1 second for master node to create {dir_name}")
                time.sleep(1)

    ############################################################################

    def parse_current_injection_info(self):

        if self.current_injection_info and self.verbose:
            self.write_log(f"Parsing current_injection_info.")

        for neuron_id, cur_info in self.current_injection_info.items():

            if int(neuron_id) not in self.neurons:
                # Neuron not on this worker.
                continue

            time = np.array(cur_info["time"])
            cur_amp = np.array(cur_info["current"])
            neuron_id = int(neuron_id)

            # Default mode is to interpolate between given points (because otherwise we need a full vector)
            interpolate_flag = cur_info["interpolate"] if "interpolate" in cur_info else True

            if self.verbose:
                self.write_log(f"Adding current injection to neuron {neuron_id}: time = {time}, current = {cur_amp}, "
                               f"{'interpolate' if interpolate_flag else ''}")

            t_vec = neuron.h.Vector(time * 1e3)
            amp_vec = neuron.h.Vector(cur_amp * 1e9)

            if neuron_id not in self.current_clamps:
                self.current_clamps[neuron_id] = []

            i_clamp = self.sim.neuron.h.IClamp(0.5, sec=self.neurons[neuron_id].icell.soma[0])
            i_clamp.dur = 1e9
            amp_vec.play(i_clamp._ref_amp, t_vec, interpolate_flag)

            self.current_clamps[neuron_id].append((i_clamp, t_vec, amp_vec))

    def add_current_injection(self, neuron_id, start_time, end_time, amplitude):

        """
        Adds current injection to neuron_id starting at start_time, ending at end_time with amplitude.

        Args:
             neuron_id: Neuron ID
             start_time: Start time of current injection
             end_time: End time of current injection
             amplitude: Amplitude of current injection
        """

        if neuron_id not in self.neuron_id:
            # The neuron ID does not exist on this worker
            return

        assert end_time > start_time, "add_current_injection: End time must be after start time"

        cur_stim = self.sim.neuron.h.IClamp(0.5, sec=self.neurons[neuron_id].icell.soma[0])
        cur_stim.delay = start_time * 1e3
        cur_stim.dur = (end_time - start_time) * 1e3
        cur_stim.amp = amplitude * 1e9  # What is units of amp?? nA??

        self.i_stim.append(cur_stim)

    ############################################################################

    def add_current_pulses(self, neuron_id, start_times, end_times, amplitudes, amp_spread=None):

        """ Inject current pulses into a neuron

        Args:
            neuron_id (int) : ID of neuron receiving current injection
            start_times (list, ndarray) : List of start times
            end_times (list, ndarray) : List of end times
            amplitudes (list, ndarray) : List of amplitudes
            amp_spread (float) : Uncertainty in the amplitudes
        """

        assert type(neuron_id) == int, "add_current_pulses only works on one neuron_id at a time"

        if type(start_times) != np.ndarray:
            start_times = np.array(start_times)

        if type(end_times) != np.ndarray:
            end_times = np.array(end_times)

        if type(amplitudes) != np.ndarray:
            amplitudes = np.array(amplitudes)

        if amp_spread is not None:
            amplitudes += np.random.normal(loc=0, scale=amp_spread, size=amplitudes.shape)

        assert (end_times - start_times).all() > 0, \
            (f"All start times must be before corresponding end times: "
             f"\nStart times: {start_times}\nEnd times: {end_times}")

        if neuron_id not in self.neuron_id:
            return  # The neuron ID does not exist on this worker

        if len(amplitudes) == 1 and len(start_times) > 1:
            amplitudes = np.repeat(amplitudes[0], len(start_times))

        assert (end_times - start_times > 0).all(), \
            (f"End time must be after start time for each time pair"
             f"Start time {start_times}, End time {end_times}")

        all_times = np.concatenate([start_times, end_times])
        all_cur = np.concatenate([amplitudes, np.zeros(amplitudes.shape)])

        idx = np.argsort(all_times)

        all_times = all_times[idx]
        all_cur = all_cur[idx]

        t_vec = neuron.h.Vector(all_times * 1e3)
        amp_vec = neuron.h.Vector(all_cur * 1e9)

        i_clamp = self.sim.neuron.h.IClamp(0.5, sec=self.neurons[neuron_id].icell.soma[0])
        i_clamp.dur = 1e9

        amp_vec.play(i_clamp._ref_amp, t_vec)

        self.i_stim.append((i_clamp, t_vec, amp_vec))

    ############################################################################

    def add_noise(self, neuron_id, duration, start_time=0, noise_amp=0, noise_std=0.1e-9, dt=None):

        """
            Add noise starting at time 0, for duration.

            Args:
                neuron_id (int) : ID of neuro injected
                duration (float) : How long is the pulse (in seconds)
                start_time (float) : When does the pulse start (in seconds)
                noise_amp (float) : Mean mplitude of noise (in ampere)
                noise_std (float) : Standard deviation of noise (in ampere)
                dt (float) : How often does the noise change (in seconds, default = 0.001s)

        """
        if dt is None:
            dt = h.dt
        else:
            dt *= 1e3  # Convert to ms

        t_vec = h.Vector(np.linspace(start_time*1e3, duration*1e3, int(duration*1e3 / dt)))

        noise_current = np.random.normal(noise_amp*1e9, noise_std*1e9, len(t_vec))
        noise_current_vector = h.Vector()
        noise_current_vector.from_python(noise_current)

        i_clamp = self.sim.neuron.h.IClamp(0.5, sec=self.neurons[neuron_id].icell.soma[0])
        i_clamp.delay = 0.0
        i_clamp.dur = 1e9
        noise_current_vector.play(i_clamp._ref_amp, t_vec, True)

        self.i_stim.append((i_clamp, t_vec, noise_current_vector))

    ############################################################################

    def convert_to_natural_units(self, param_name, param_value):

        # TODO, move conversion list to separate file
        if param_name in self.conv_factor:
            val = param_value * self.conv_factor[param_name]
        else:
            val = param_value

        return val

    ############################################################################

    # Use event handling

    def setup_print_sim_time(self, t_max):
        """ Setup event handler to print progress, and estimate of simulation run time"""

        # Only have the first node print time estimates
        if self.pc.id() == 0:
            self.t_max = t_max
            self.sim_start_time = timeit.default_timer()
            self.fih_time = h.FInitializeHandler((self._setup_print_sim_time_helper, t_max))
            self.last_sim_report_time = 0

    ############################################################################

    def _setup_print_sim_time_helper(self, t_max):
        """ Helper method for printing simulation time during execution. """
        update_points = np.array([0.0, 0.01, 0.02, 0.03, 0.04, 0.05,
                                  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]) * t_max

        # update_points = np.arange(t_max / 100., t_max, t_max / 100.)
        for t in update_points:
            h.cvode.event(t, self.print_sim_time)

    ############################################################################

    def print_sim_time(self):
        """ Helper function, to print simulation time during execution. """

        if h.t == 0:
            # Starting simulation wall clock, for estimate of time left
            self.sim_start_time = timeit.default_timer()
            return

        cur_time = timeit.default_timer()
        elapsed_time = cur_time - self.sim_start_time
        fraction_done = h.t / self.t_max
        time_left = np.abs(elapsed_time * ((self.t_max - h.t) / h.t))

        # Do not print status update too often
        if cur_time - self.last_sim_report_time > 100 or fraction_done > 0.99:
            force_print = True
            self.last_sim_report_time = cur_time
        else:
            force_print = False

        self.write_log(f"{fraction_done * 100:>3.0f}% done. Elapsed: {elapsed_time:.1f} s, "
                       f"estimated time left: {time_left:.1f} s", force_print=force_print)

        # self.write_log("%.0f%% done. Elapsed: %.1f s, estimated time left: %.1f s"
        #                % (fraction_done * 100, elapsed_time, time_left), force_print=force_print)

    ############################################################################

    def check_memory_status(self, threshold=0.1):
        """ Checks memory status. """

        mem_available, mem_total = snudda.utils.memory.memory_status()

        memory_ratio = mem_available / mem_total

        if self.pc is not None:
            pc_id = self.pc.id()
        else:
            pc_id = "Unknown node"

        self.write_log(f"{pc_id} : Memory status: {int(memory_ratio * 100)}% free", force_print=True)

        return memory_ratio < threshold

    def __del__(self):
        if self is not None:
            self.clear_neuron()

        if self.log_file is not None:
            try:
                self.log_file.close()
            except:
                pass

    def clear_neuron(self):

        if self.pc is not None:
            self.pc.gid_clear()

        self.neurons = {}
        self.sim = None
        self.neuron_nodes = []  # Is this used?
        self.virtual_neurons = {}

        self.synapse_dict = dict()
        self.i_stim = []
        self.v_clamp_list = []
        self.gap_junction_dict = dict()
        self.external_stim = dict([])
        self.check_id_recordings = []
        self.pc = None

        self.t_save = []
        self.i_save = []
        self.i_key = []

        self.fih_time = None  # Event handler to print estimated time left

        self.record = None

        # We also need to clear everything setup by NEURON
        for sec in h.allsec():
            h.delete_section(sec=sec)

        gc.collect()

############################################################################

#
# Test code to run a simulation


if __name__ == "__main__":

    # Please see the wrapper script snudda.py which shows how to generate
    # a network, and how to then simulate it using this script
    import sys

    if '-python' in sys.argv:
        print("Network_simulate.py called through nrniv, fixing arguments")
        pythonidx = sys.argv.index('-python')
        if len(sys.argv) > pythonidx:
            sys.argv = sys.argv[pythonidx + 1:]

    import argparse

    parser = argparse.ArgumentParser(description="Simulate network generated by Snudda")
    parser.add_argument("networkFile", help="Network model (HDF5)")
    parser.add_argument("inputFile", help="Network input (HDF5)")
    parser.add_argument("--noVolt", "--novolt", action="store_false", dest="record_volt",
                        help="Exclude voltage when saving results, saves time and space.")
    parser.add_argument("--recordALLcompartments", dest="record_all_compartments", type=str, default=None)
    parser.add_argument("--recordALLsynapses", dest="record_all_synapses", type=str, default=None)
    parser.add_argument("--disableSyn", "--disableSynapses", action="store_true", dest="disable_synapses", default=None,
                        help="Disable synapses")
    parser.add_argument("--disableGJ", "--disableGapJunctions", action="store_true", dest="disable_gapjunctions", default=None,
                        help="Disable gap junctions")
    parser.add_argument("--time", type=float, default=None,
                        help="Duration of simulation in seconds (must be set, or specified in simulation_config)")
    parser.add_argument("--verbose", action="store_true")

    parser.add_argument("--outputFile", help="Output hdf5 file (from simulation)", dest="output_file", default=None)
    parser.add_argument("--simulation_config", type=str, default=None)


    # If called through "nrniv -python Network_simulate.py ..." then argparse
    # gets confused by -python flag, and we need to ignore it
    # parser.add_argument("-python",help=argparse.SUPPRESS,
    #                    action="store_true")

    args = parser.parse_args()
    network_data_file = args.networkFile
    input_file = args.inputFile
    log_file = os.path.join(os.path.dirname(args.networkFile), "log", "network-simulation-log.txt")
    save_dir = os.path.join(os.path.dirname(args.networkFile), "simulation")

    if not os.path.exists(save_dir):
        print(f"Creating directory {save_dir}")
        os.makedirs(save_dir, exist_ok=True)

    # Get the SlurmID, used in default file names
    slurm_id = os.getenv('SLURM_JOBID')

    if slurm_id is None:
        digits = re.findall(r'\d+', input_file)
        # Second to last digit is slurmID of old run, reuse that
        try:
            slurm_id = digits[-2]
        except:
            slurm_id = str(666)

    start = timeit.default_timer()

    if args.disable_gapjunctions:
        print("!!! WE HAVE DISABLED GAP JUNCTIONS !!!")

    if args.disable_synapses:
        print("!!! WE HAVE DISABLED SYNAPSES !!!")

    pc = h.ParallelContext()

    if args.output_file:
        output_file = args.output_file
    else:
        output_file = None

    sim = SnuddaSimulate(network_file=network_data_file,
                         input_file=input_file,
                         output_file=output_file,
                         disable_gap_junctions=args.disable_gapjunctions,
                         disable_synapses=args.disable_synapses,
                         log_file=log_file,
                         simulation_config=args.simulation_config,
                         verbose=args.verbose)
    sim.setup()
    sim.add_external_input()
    sim.check_memory_status()

    if args.record_volt:
        sim.write_log(f"Adding soma voltage recording to neurons.")
        sim.add_volt_recording_soma()
    else:
        sim.write_log(f"Somatic voltage not recorded.")

    if args.record_all_compartments:
        record_cell_id = np.array([int(x) for x in args.record_all_compartments.split(",")])
        sim.add_volt_recording_all(cell_id=record_cell_id, centre_only_flag=True)

    if args.record_all_synapses:
        record_cell_id = np.array([int(x) for x in args.record_all_synapses.split(",")])
        sim.add_synapse_current_recording_all(record_cell_id)

    if args.time:
        t_sim = args.time * 1000  # Convert from s to ms for Neuron simulator
    else:
        t_sim = None

    sim.check_memory_status()
    # print(f"Running simulation for {t_sim} ms.")
    sim.run(t_sim)  # In milliseconds

    sim.write_log("Simulation done, saving output")
    sim.write_output()

    stop = timeit.default_timer()
    if sim.pc.id() == 0:
        sim.write_log(f"Program run time: {(stop - start):.0f}")

    # sim.plot()
    sys.exit(0)

# Check this code example
# Why are spikes not propagated from one neuron to another
# https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=188544&file=%2FLyttonEtAl2016%2FREADME.html#tabs-2
