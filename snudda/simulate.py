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

# Plot all sections
# [neuron.h.psection(x) for x in neuron.h.allsec()]


from mpi4py import MPI  # This must be imported before neuron, to run parallel
from neuron import h  # , gui
import neuron
import h5py
import json
import timeit

import bluepyopt.ephys as ephys
from snudda.neuron_model_extended import NeuronModel
# from Network_place_neurons import NetworkPlaceNeurons
import numpy as np
from snudda.nrn_simulator_parallel import NrnSimulatorParallel

from glob import glob
import re
import os

import snudda.utils.memory

# !!! Need to gracefully handle the situation where there are more workers than
# number of neurons, currently we get problem when adding the voltage saving

# !!! Have added code for dopamine modulation of neuron intrinsic channels
#     need to add dopamine modulation for the synaptic channels also !!!

##############################################################################

# If simulationConfig is set, those values override other values
from snudda.load import SnuddaLoad


class SnuddaSimulate(object):

    def __init__(self, network_file, input_file=None,
                 verbose=True, log_file=None,
                 disable_gap_junctions=True,
                 simulation_config=None):

        self.verbose = verbose
        self.log_file = log_file

        self.network_file = network_file
        self.input_file = input_file

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
        self.synapse_parameters = None

        if simulation_config:
            sim_info = json.load(simulation_config)

            if "networkFile" in sim_info:
                self.network_file = network_file

            if "inputFile" in sim_info:
                self.input_file = input_file

            if log_file in sim_info:
                self.log_file = open(log_file, "w")

        if type(self.log_file) == str:
            self.log_file = open(self.log_file, "w")

        self.write_log("Using networkFile: " + str(network_file))
        self.write_log("Using inputFile: " + str(input_file))

        if self.log_file is not None:
            self.write_log("Using logFile: " + str(self.log_file.name))

        # !!! What value to use for synaptic weight and synapse delay?
        # !!! different for AMPA and GABA?
        self.synapse_weight = 10.0  # microsiemens
        self.synapse_delay = 1  # ms
        self.spike_threshold = -20
        self.axon_speed = 0.8  # Tepper and Lee 2007, Wilson 1986, Wilson 1990
        # refs taken from Damodaran et al 2013

        self.disable_gap_junctions = disable_gap_junctions

        self.synapse_type_lookup = {1: "GABA", 2: "AMPA_NMDA", 3: "GapJunction"}

        self.neurons = {}
        self.sim = None
        self.neuron_nodes = []

        self.virtual_neurons = {}

        self.net_con_list = []  # Avoid premature garbage collection
        self.synapse_list = []
        self.i_stim = []
        self.v_clamp_list = []
        self.gap_junction_list = []
        self.external_stim = dict([])
        self.t_save = []
        self.v_save = []
        self.v_key = []
        self.i_save = []
        self.i_key = []

        self.input_data = None

        self.gap_junction_next_gid = 0  # Are these gids separate from cell gids?

        self.t_spikes = h.Vector()
        self.id_spikes = h.Vector()

        # Make sure the output dir exists, so we don't fail at end because we
        # cant write file
        self.create_dir("save/traces")

        self.pc = h.ParallelContext()

        # self.writeLog("I am node " + str(int(self.pc.id())))

        # We need to initialise random streams, see Lytton el at 2016 (p2072)

        self.load_network_info(network_file)

        self.check_memory_status()
        self.distribute_neurons()
        self.setup_neurons()
        self.check_memory_status()
        self.pc.barrier()

        #    for i in range(0,self.nNeurons):
        #      print("Node : " + str(int(self.pc.id())) + " cell " + str(i) + " status " + str(self.pc.gid_exists(i)))

        self.connect_network()
        self.check_memory_status()
        self.pc.barrier()

        # Do we need blocking call here, to make sure all neurons are setup
        # before we try and connect them

        # READ ABOUT PARALLEL NEURON

    # https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/network/parcon.html#paralleltransfer

    ############################################################################

    def load_network_info(self, network_file, config_file=None):

        self.network_file = network_file

        if False:
            # We need to check if a split network file exists
            split_file = network_file.replace('save/', 'save/TEMP/').replace('.hdf5', '-%d.hdf5') % int(self.pc.id())

            import os.path
            if os.path.isfile(split_file):
                network_file = split_file
            else:
                self.write_log("Unable to find " + split_file + " using "
                               + network_file)

        self.write_log("Worker " + str(int(self.pc.id()))
                       + ": Loading network from " + network_file)

        from snudda.load import SnuddaLoad
        self.snudda_loader = SnuddaLoad(network_file)
        self.network_info = self.snudda_loader.data

        self.synapses = self.network_info["synapses"]
        self.gap_junctions = self.network_info["gapJunctions"]

        # We are only passed information about neurons on our node if
        # SplitConnectionFile was run, so need to use nNeurons to know
        # how many neurons in total
        self.num_neurons = self.network_info["nNeurons"]

        if config_file is None:
            config_file = self.get_path(self.network_info["configFile"])

        self.config_file = config_file
        self.write_log("Loading config file " + config_file)

        # Add checks to see that config file and networkFile matches

        import json
        with open(config_file, 'r') as config_file:
            self.config = json.load(config_file)

        # I do not know if the gap junction GIDs are a separate entity from the
        # neuron cell GIDs, so to be on safe side, let's make sure they
        # do not overlap
        self.gap_junction_next_gid = self.num_neurons + 100000000

        # Make a bool array indicating if cells are virtual or not
        self.is_virtual_neuron = [n["virtualNeuron"] for n in self.network_info["neurons"]]

    ############################################################################

    def distribute_neurons(self):
        # This code is run on all workers, will generate different lists on each
        self.write_log("Distributing neurons.")

        self.neuron_id = range(int(self.pc.id()), self.num_neurons, int(self.pc.nhost()))

        self.neuron_nodes = [x % int(self.pc.nhost()) for x in range(0, self.num_neurons)]

        # self.write_log("Node " + str(int(self.pc.id())) + " handling neurons: " + ' '.join(map(str, self.neuron_id)))

    ############################################################################

    def destroy(self):
        for n in self.neurons:
            n.destroy(sim=self.sim)

    ############################################################################

    # This requires self.sim to be defined

    def load_synapse_parameters(self):

        # We need to load all the synapse parameters
        self.synapse_parameters = dict([])

        for (preType, postType) in self.network_info["connectivityDistributions"]:

            syn_data = self.network_info["connectivityDistributions"][preType, postType]

            for synType in syn_data:

                synapse_type_id = syn_data[synType]["channelModelID"]
                info_dict = syn_data[synType]

                if synapse_type_id == 3:
                    # Gap junctions, skip parameters
                    continue

                if ("channelParameters" in info_dict
                        and info_dict["channelParameters"] is not None):
                    channel_param_dict = info_dict["channelParameters"].copy()
                    mod_file = channel_param_dict["modFile"]

                    eval_str = "self.sim.neuron.h." + mod_file
                    channel_module = eval(eval_str)  # If this fails, check that NEURON modules are compiled

                    # These are not variables to set in the modFile
                    if "modFile" in channel_param_dict:
                        del channel_param_dict["modFile"]

                    if "parameterFile" in channel_param_dict:
                        del channel_param_dict["parameterFile"]

                else:
                    assert False, "No channel module specified for " \
                                  + str(preType) + "->" + str(postType) + " synapses, type ID= " \
                                  + str(synapse_type_id)

                if "parameterFile" in info_dict and info_dict["parameterFile"] is not None:
                    par_file = self.get_path(info_dict["parameterFile"])
                    par_data_dict = json.load(open(par_file, 'r'))

                    # Save data as a list, we dont need the keys
                    par_data = []
                    for pd in par_data_dict:
                        if "synapse" in par_data_dict[pd]:

                            # Add channel parameters specified in network file, however
                            # any values in the synapse parameter file will overwrite them
                            p_dict = channel_param_dict.copy()
                            for x in par_data_dict[pd]["synapse"]:
                                p_dict[x] = par_data_dict[pd]["synapse"][x]

                            par_data.append(p_dict)
                        else:
                            self.write_log("WARNING: Old data format in parameter file " + str(par_file))

                            p_dict = channel_param_dict.copy()
                            for x in par_data_dict[pd]:
                                p_dict[x] = par_data_dict[pd][x]

                            par_data.append(p_dict)
                elif len(channel_param_dict) > 0:
                    par_data = [channel_param_dict]
                else:
                    par_data = None

                self.synapse_parameters[synapse_type_id] = (channel_module, par_data)

    ############################################################################

    def setup_neurons(self):

        self.write_log("Setup neurons")

        # self.sim = ephys.simulators.NrnSimulator(cvode_active=False)
        # self.sim = NrnSimulatorParallel()
        self.sim = NrnSimulatorParallel(cvode_active=False)

        # We need to load all the synapse parameters
        self.load_synapse_parameters()

        # The neurons this node is responsible for is in self.neuronID
        for ID in self.neuron_id:

            name = self.network_info["neurons"][ID]["name"]

            config = self.config["Neurons"][name]

            morph = self.get_path(config["morphology"])
            param = self.get_path(config["parameters"])
            mech = self.get_path(config["mechanisms"])

            if "modulation" in config:
                modulation = self.get_path(config["modulation"])
            else:
                modulation = None

            # Obs, neurons is a dictionary
            if self.network_info["neurons"][ID]["virtualNeuron"]:

                if self.input_data is None:
                    self.write_log("Using " + self.input_file + " for virtual neurons")
                    self.input_data = h5py.File(self.get_path(self.input_file), 'r')

                name = self.network_info["neurons"][ID]["name"]
                spikes = self.input_data["input"][ID]["activity"]["spikes"][:, 0]

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
                parameter_id = self.network_info["neurons"][ID]["parameterID"]
                modulation_id = self.network_info["neurons"][ID]["modulationID"]

                self.neurons[ID] = NeuronModel(param_file=param,
                                               morph_file=morph,
                                               mech_file=mech,
                                               cell_name=name,
                                               modulation_file=modulation,
                                               parameter_id=parameter_id,
                                               modulation_id=modulation_id)

                # Register ID as belonging to this worker node
                self.pc.set_gid2node(ID, int(self.pc.id()))

                if True or False:
                    self.write_log("Node " + str(int(self.pc.id())) + " - cell " + str(ID) + " " + name)

                # We need to instantiate the cell
                self.neurons[ID].instantiate(sim=self.sim)
                self.set_resting_voltage(ID)

                # !!! DIRTY FIX for
                # https://github.com/BlueBrain/BluePyOpt/blob/master/bluepyopt/ephys/morphologies.py
                # This is likely the offending line, that pushes a segment to the stack
                # --> sim.neuron.h.execute('create axon[2]', icell)

                self.write_log("!!! Popping extra segment from neuron -- temp fix!")
                h.execute("pop_section()")

                # !!! END OF DIRTY FIX

                # !!! Connect a netcon and register it, taken from ballandstick's
                #     connect2target function
                nc = h.NetCon(self.neurons[ID].icell.axon[0](0.5)._ref_v,
                              None,
                              sec=self.neurons[ID].icell.axon[0])
                nc.threshold = 10

                self.pc.cell(ID, nc, 1)  # The 1 means broadcast spikes to other machines
                # self.pc.outputcell(ID) # -- not needed, cell was called with a 1
                # self.netConList.append(nc) -- Not needed according to Lytton et al 2016

                # Record all spikes
                self.pc.spike_record(ID, self.t_spikes, self.id_spikes)

    ############################################################################

    def connect_network(self):

        self.pc.barrier()

        # Add synapses
        self.connect_network_synapses()

        # Add gap junctions
        if self.disable_gap_junctions:
            self.write_log("!!! Gap junctions disabled.")
        else:
            self.write_log("Adding gap junctions.")
            # self.connectNetworkGapJunctions()

            # TODO: Check difference with old non-local version
            self.connect_network_gap_junctions_local()
            self.pc.setup_transfer()
        self.pc.barrier()

    ############################################################################

    def connect_network_synapses(self):

        self.write_log("connectNetworkSynapses")

        # This loops through all the synapses, and connects the relevant ones
        next_row = 0
        # nextRowSet = [ fromRow, toRow ) -- ie range(fromRow,toRow)
        next_row_set = self.find_next_synapse_group(next_row)

        while next_row_set is not None:
            # Add the synapses to the neuron
            self.connect_neuron_synapses(start_row=next_row_set[0], end_row=next_row_set[1])

            # Find the next group of synapses
            next_row = next_row_set[1]  # 2nd number was not included in range
            next_row_set = self.find_next_synapse_group(next_row)

    ############################################################################

    # This function starts at nextRow, then returns all synapses onto
    # a neuron which is located on the worker

    # This works for synapses, but it will not work for gap junctions, because
    # we need to connect the gap junctions from both sides

    # --- perhaps rewrite this as an iterator

    def find_next_synapse_group(self, next_row=0, connection_type="synapses"):

        if connection_type == "synapses":
            synapses = self.synapses
        elif connection_type == "gapjunctions":
            synapses = self.gap_junctions
        else:
            assert False, "!!! find_next_synapse_group: Unknown connectionType: " + connection_type

        try:
            num_syn_rows = synapses.shape[0]
        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)

            assert False, "find_next_synapse_group: If synapses was not loaded into memory, your problem is probably " \
                          + "that the HDF5 file that holds the synapses were closed. Sorry."

        if next_row >= num_syn_rows:
            # No more synapses to get
            return None

        # The synapse matrix is sorted on destID, ascending order
        # We also assume that self.neuronID is sorted in ascending order

        start_row = None
        our_id = None
        not_our_id = None

        while start_row is None:

            # What is the next destination ID
            next_id = synapses[next_row, 1]

            # Is the next ID ours?
            if next_id in self.neuron_id:
                start_row = next_row
                our_id = next_id
                continue
            else:
                not_our_id = next_id

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

    # Processing the range(startRow,endRow) (ie, to endRow-1)

    def connect_neuron_synapses(self, start_row, end_row):

        source_id_list = self.synapses[start_row:end_row, 0]
        dest_id = self.synapses[start_row, 1]
        assert (self.synapses[start_row:end_row, 1] == dest_id).all()

        # Double check mapping
        assert self.pc.gid2cell(dest_id) == self.neurons[dest_id].icell, \
            "GID mismatch: " + str(self.pc.gid2cell(dest_id)) \
            + " != " + str(self.neurons[dest_id].icell)

        synapse_type_id = self.synapses[start_row:end_row, 6]
        axon_distance = self.synapses[start_row:end_row, 7]  # Obs in micrometers

        sec_id = self.synapses[start_row:end_row, 9]
        dend_sections = self.neurons[dest_id].map_id_to_compartment(sec_id)
        sec_x = self.synapses[start_row:end_row, 10] / 1000.0  # Convert to number 0-1

        # conductances are stored in pS (because we use INTs),
        # Neuron wants it in microsiemens??!
        conductance = self.synapses[start_row:end_row, 11] * 1e-6
        parameter_id = self.synapses[start_row:end_row, 12]

        voxel_coords = self.synapses[start_row:end_row, 2:5]
        self.verify_synapse_placement(dend_sections, sec_x, dest_id, voxel_coords)

        for (src_id, section, section_x, s_type_id, axon_dist, cond, p_id) \
                in zip(source_id_list, dend_sections, sec_x, synapse_type_id,
                       axon_distance, conductance, parameter_id):

            try:
                # !!!
                self.add_synapse(cell_id_source=src_id,
                                 dend_compartment=section,
                                 section_dist=section_x,
                                 synapse_type_id=s_type_id,
                                 axon_dist=axon_dist,
                                 conductance=cond,
                                 parameter_id=p_id)
            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr)
                import pdb
                pdb.set_trace()

    ############################################################################

    # OBS!! The src and dest lists can be different length
    #
    # src are all the gap junctions where the source compartment are
    # on the local worker.
    # dest are the gap junctions where the dest compartment are on the
    # local worker
    # The same GJ might appear in both src and dest lists, but at different rows

    def find_local_gap_junctions(self):

        # If the gap junction matrix is too large to fit in memory then
        # this will need to be optimised

        self.write_log("Finding node local gap junctions...")

        gj_idx_a = np.where([x in self.neuron_id for x in self.gap_junctions[:, 0]])[0]
        gj_idx_b = np.where([x in self.neuron_id for x in self.gap_junctions[:, 1]])[0]

        # GJIDoffset = self.network_info["GJIDoffset"]
        gj_id_offset = 100 * self.num_neurons
        gj_gid_src_a = gj_id_offset + 2 * gj_idx_a
        gj_gid_src_b = gj_id_offset + 2 * gj_idx_b + 1

        gjg_id_dest_a = gj_id_offset + 2 * gj_idx_a + 1
        gj_gid_dest_b = gj_id_offset + 2 * gj_idx_b + 0

        neuron_id_a = self.gap_junctions[gj_idx_a, 0]
        neuron_id_b = self.gap_junctions[gj_idx_b, 1]

        seg_id_a = self.gap_junctions[gj_idx_a, 2]
        seg_id_b = self.gap_junctions[gj_idx_b, 3]

        compartment_a = [self.neurons[x].map_id_to_compartment([y])[0] for (x, y) in zip(neuron_id_a, seg_id_a)]
        compartment_b = [self.neurons[x].map_id_to_compartment([y])[0] for (x, y) in zip(neuron_id_b, seg_id_b)]

        seg_xa = self.gap_junctions[gj_idx_a, 4] / 1000.0
        seg_xb = self.gap_junctions[gj_idx_b, 5] / 1000.0

        # Since we had ints we stored pS, but Neuron wants microsiemens
        cond_a = self.gap_junctions[gj_idx_a, 10] * 1e-6
        cond_b = self.gap_junctions[gj_idx_b, 10] * 1e-6

        # Merge the two lists together

        gj_idx = np.concatenate([gj_idx_a, gj_idx_b])
        gj_gid_src = np.concatenate([gj_gid_src_a, gj_gid_src_b])
        gj_gid_dest = np.concatenate([gjg_id_dest_a, gj_gid_dest_b])
        neuron_id = np.concatenate([neuron_id_a, neuron_id_b])
        seg_id = np.concatenate([seg_id_a, seg_id_b])
        compartment = np.concatenate([compartment_a, compartment_b])
        seg_x = np.concatenate([seg_xa, seg_xb])
        cond = np.concatenate([cond_a, cond_b])

        return neuron_id, compartment, seg_x, gj_gid_src, gj_gid_dest, cond

    ############################################################################

    # We can only do half the setup of the gap junction if it is split between
    # two workers.

    def connect_network_gap_junctions_local(self):

        self.write_log("connectNetworkGapJunctionsLOCAL")

        (neuron_id, compartment, seg_x, gj_gid_src, gj_gid_dest, cond) \
            = self.find_local_gap_junctions()

        # import pdb
        # pdb.set_trace()

        try:
            # WHY??!
            # ValueError: too many values to unpack (expected 6)

            for nid, comp, s_x, gid_src, gid_dest, g \
                    in zip(neuron_id, compartment, seg_x, gj_gid_src, gj_gid_dest, cond):
                self.add_gap_junction(section=comp,
                                      section_dist=s_x,
                                      gid_source_gj=gid_src,
                                      gid_dest_gj=gid_dest,
                                      g_gap_junction=g)

        except:
            import traceback
            tstr = traceback.format_exc()
            print(tstr)
            import pdb
            pdb.set_trace()

    ############################################################################

    def connect_network_gap_junctions(self):

        self.write_log("connectNetworkGapJunctions")

        self.write_log("!!! Please verify connectNeuronGapJunctions function, that currents go bidirectionally")

        # This loops through all the synapses, and connects the relevant ones
        next_row = 0
        # nextRowSet = [ fromRow, toRow ) -- ie range(fromRow,toRow)
        next_row_set = self.find_next_synapse_group(next_row,
                                                    connection_type="gapjunctions")

        while next_row_set is not None:
            # Add the synapses to the neuron
            self.connect_neuron_gap_junctions(start_row=next_row_set[0],
                                              end_row=next_row_set[1])

            # Find the next group of synapses
            next_row = next_row_set[1]  # 2nd number was not included in range
            next_row_set = self.find_next_synapse_group(next_row,
                                                        connection_type="gapjunctions")

    ############################################################################

    # Verify that the gap junctions conducts currents in both directions

    def connect_neuron_gap_junctions(self, start_row, end_row):

        source_id = self.gap_junctions[start_row:end_row, 0]
        dest_id = self.gap_junctions[start_row, 1]

        assert (self.gap_junctions[start_row:end_row, 1] == dest_id).all()

        # Double check mapping
        assert self.pc.gid2cell(dest_id) == self.neurons[dest_id].icell, \
            "GID mismatch: " + str(self.pc.gid2cell(dest_id)) \
            + " != " + str(self.neurons[dest_id].icell)

        source_sec_id = self.gap_junctions[start_row:end_row, 2]
        dest_sec_id = self.gap_junctions[start_row:end_row, 3]

        # !!! Double check we get number between 0.0 and 1.0
        source_sec_x = self.gap_junctions[start_row:end_row, 4] * 1e-4
        dest_sec_x = self.gap_junctions[start_row:end_row, 5] * 1e-4

        # conductances are stored in pS, Neuron wants it in microsiements??!
        # (reason for not storing SI units is that we use INTs)
        conductance = self.gap_junctions[start_row:end_row, 10] * 1e-6

        dest_loc = self.neurons[dest_id].map_id_to_compartment(dest_sec_id)

        for (srcID, srcSecID, srcSecX, dstLoc, dstSecX, rowIdx) \
                in zip(source_id, source_sec_id, source_sec_x,
                       dest_loc, dest_sec_x, range(start_row, end_row)):

            src_loc = self.neurons[srcID].map_id_to_compartment([srcSecID])[0]

            # Change the src and dest ID to be based on the row idx
            # to avoid overlaps
            gj_src_id = rowIdx * 2 + 10 * self.num_neurons
            gj_dest_id = gj_src_id + 1

            print("rowIdx = " + str(rowIdx))
            print("GJsrcID = " + str(gj_src_id))
            print("GJdestID = " + str(gj_dest_id))

            try:
                self.add_gap_junction(section=src_loc,
                                      section_dist=srcSecX,
                                      gid_source_gj=gj_src_id,
                                      gid_dest_gj=gj_dest_id)

                self.add_gap_junction(section=dstLoc,
                                      section_dist=dstSecX,
                                      gid_source_gj=gj_dest_id,
                                      gid_dest_gj=gj_src_id)
            except:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)
                import pdb
                pdb.set_trace()

    ############################################################################

    def find_gap_junction_compartments(self):

        all_loc = dict([])

        orig_gj_coords = self.network_info["origGJCoords"]

        for ID in self.neuron_id:

            idx_gj1 = np.where(np.logical_and(self.network_info["synapses"][:, 0] == ID,
                                              self.network_info["synapses"][:, 5] == 3))

            idx_gj2 = np.where(np.logical_and(self.network_info["synapses"][:, 2] == ID,
                                              self.network_info["synapses"][:, 5] == 3))

            # Coordinates on either side of the gap junction
            gj_coords1 = np.array([orig_gj_coords[x][0:3] for x in idx_gj1[0]])
            gj_coords2 = np.array([orig_gj_coords[x][3:6] for x in idx_gj2[0]])

            gj_id1 = np.array([orig_gj_coords[x][6:8] for x in idx_gj1[0]])
            gj_id2 = np.array([orig_gj_coords[x][6:8] for x in idx_gj2[0]])

            if gj_coords1.shape[0] == 0:
                gj_coords = gj_coords2
            elif gj_coords2.shape[0] == 0:
                gj_coords = gj_coords1
            else:
                gj_coords = np.concatenate([gj_coords1, gj_coords2], axis=0)

            gj_loc_type = 4 * np.ones(shape=(gj_coords.shape[0], 1))

            len_gj1 = len(gj_coords1)

            # import pdb
            # pdb.set_trace()

            if gj_coords.shape[0] > 0:

                self.write_log("Looking for " + str(gj_coords.shape[0]) + " gap junctions")

                # Get the compartment location of each coordinate
                gj_dend_loc = self.neurons[ID].find_dend_compartment(gj_coords,
                                                                     gj_loc_type,
                                                                     self.sim)
                gj_dend_loc1 = gj_dend_loc[:len_gj1]
                gj_dend_loc2 = gj_dend_loc[len_gj1:]

                assert (gj_coords1.shape[0] == len(gj_dend_loc1))
                assert (gj_coords2.shape[0] == len(gj_dend_loc2))

                for (idx, loc, id1) in zip(idx_gj1[0], gj_dend_loc1, gj_id1):
                    all_loc[(idx, 1)] = (loc, id1[0], id1[1])

                for (idx, loc, id2) in zip(idx_gj2[0], gj_dend_loc2, gj_id2):
                    all_loc[(idx, 2)] = (loc, id2[0], id2[1])

        return all_loc

    ############################################################################

    def add_synapse(self, cell_id_source, dend_compartment, section_dist, conductance,
                    parameter_id, synapse_type_id, axon_dist=None):

        # You can not locate a point process at
        # position 0 or 1 if it needs an ion
        if section_dist == 0.0:
            section_dist = 0.01
        if section_dist == 1.0:
            section_dist = 0.99

        (channel_module, par_data) = self.synapse_parameters[synapse_type_id]

        syn = channel_module(dend_compartment(section_dist))

        if par_data is not None:
            # Picking one of the parameter sets stored in parData
            par_id = parameter_id % len(par_data)

            par_set = par_data[par_id]
            for par in par_set:
                if par == "expdata" or par == "cond":
                    # expdata is not a parameter, and cond we take from synapse matrix
                    continue

                try:
                    # Can be value, or a tuple/list, if so second value is scale factor
                    # for SI -> natural units conversion
                    val = par_set[par]

                    # Do we need to convert from SI to natural units?
                    if type(val) == tuple or type(val) == list:
                        val = val[0] * val[1]

                    setattr(syn, par, val)

                    # evalStr = "syn." + par + "=" + str(parSet[par])
                    # self.writeLog("Updating synapse: " + evalStr)
                    # !!! Can we avoid an eval here, it is soooo SLOW
                    # exec(evalStr)
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)
                    import pdb
                    pdb.set_trace()

        # Just create a default expsyn for test, will need to create proper GABA
        # synapses later
        # if(synapseType == 'ExpSyn'):
        #  syn = self.sim.neuron.h.ExpSyn(dendCompartment(sectionDist))
        # elif(synapseType == 'GABA'):
        #  syn = self.sim.neuron.h.tmGabaA(dendCompartment(sectionDist))
        # elif(synapseType == "AMPA_NMDA"):
        #  syn = self.sim.neuron.h.tmGlut(dendCompartment(sectionDist))
        # else:
        #  self.writeLog("Synapse type not implemented: ", synapseType)
        #  import pdb
        #  pdb.set_trace()

        if axon_dist is not None:
            # axon dist is in micrometer, want delay in ms
            synapse_delay = (1e3 * 1e-6 * axon_dist) / self.axon_speed + self.synapse_delay
        else:
            synapse_delay = self.synapse_delay

        #    self.write_log("Synapse delay: " + str(synapse_delay) + " ms")

        # What do we do if the GID does not exist?
        # print("GID exists:" + str(self.pc.gid_exists(cellIDsource)))

        if self.is_virtual_neuron[cell_id_source]:
            # Source is a virtual neuron, need to read and connect input
            src_name = self.network_info["neurons"][cell_id_source]["name"]

            nc = self.pc.gid_connect(cell_id_source, syn)
            nc.weight[0] = conductance
            nc.delay = synapse_delay
            nc.threshold = self.spike_threshold

            # Prevent garbage collection in python
            self.net_con_list.append(nc)
            self.synapse_list.append(syn)

        else:

            nc = self.pc.gid_connect(cell_id_source, syn)
            nc.weight[0] = conductance
            nc.delay = synapse_delay
            nc.threshold = self.spike_threshold

            self.net_con_list.append(nc)
            self.synapse_list.append(syn)

        return syn

    ############################################################################

    # Add one gap junction to specific location

    def add_gap_junction(self,
                         section, section_dist,
                         gid_source_gj, gid_dest_gj,
                         g_gap_junction=0.5e-9,
                         gid=None):  # GID unused??

        # If neuron complains, make sure you have par_ggap.mod
        gj = h.gGapPar(section(section_dist))
        self.gap_junction_list.append(gj)

        self.pc.target_var(gj._ref_vgap, gid_dest_gj)

        self.pc.source_var(section(section_dist)._ref_v, gid_source_gj, sec=section)

        gj.g = g_gap_junction
        # print("Setting conductance: " + str(GJ.g))

    ############################################################################

    # Wilson 2007 - GABAergic inhibition in the neostriatum
    # 80% of synapses in Striatum are glutamatergic
    # Estimated 10000 glutamate and 2000 GABA synapses per MS,
    # 325 dopamine synapses per MS
    # Wilson 1996 - 10000 spines per MS = 10000 glutamatergic inputs

    # Ingham et al 1998, approx 1 glutamatergic synapse per 0.92 mum3
    # --> ~11000 glutamate synapses per MS
    # Kemp 1971 -- The synaptic organization of the caudate nucleus (85% glu)

    def add_external_input(self, input_file=None):

        if input_file is None:
            input_file = self.input_file

        self.write_log("Adding external (cortical, thalamic) input from " + input_file)

        self.input_data = h5py.File(input_file, 'r')

        for neuron_id, neuron in self.neurons.items():

            self.external_stim[neuron_id] = []
            name = neuron.name

            if str(neuron_id) not in self.input_data["input"]:
                self.write_log("Warning - No input specified for " + name)
                continue

            for inputType in self.input_data["input"][str(neuron_id)]:
                neuron_input = self.input_data["input"][str(neuron_id)][inputType]

                loc_type = 1 * np.ones((neuron_input["sectionID"].shape[0],))  # Axon-Dend

                sections = self.neurons[neuron_id].map_id_to_compartment(neuron_input["sectionID"])

                # Setting individual parameters for synapses
                mod_file = SnuddaLoad.to_str(neuron_input["modFile"][()])
                param_list = json.loads(neuron_input["parameterList"][()])

                eval_str = "self.sim.neuron.h." + mod_file
                channel_module = eval(eval_str)

                for inputID, (section, section_x, paramID, nSpikes) \
                        in enumerate(zip(sections,
                                         neuron_input["sectionX"],
                                         neuron_input["parameterID"],
                                         neuron_input["nSpikes"])):
                    # We need to find cellID (int) from neuronID (string, eg. MSD1_3)

                    idx = inputID
                    spikes = neuron_input["spikes"][inputID, :nSpikes] * 1e3  # Neuron uses ms
                    assert (spikes >= 0).all(), \
                        "Negative spike times for neuron " + str(neuron_id) + " " + inputType

                    # Creating NEURON VecStim and vector
                    # https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3125
                    # import pdb
                    # pdb.set_trace()
                    try:
                        vs = h.VecStim()
                        v = h.Vector(spikes.size)
                        v.from_python(spikes)
                        vs.play(v)
                    except:
                        import traceback
                        tstr = traceback.format_exc()
                        print(tstr)

                        assert False, "!!! If you see this, make sure that vecevent.mod is included in nrnivmodl compilation"

                    # NEURON: You can not locate a point process at position 0 or 1
                    # if it needs an ion
                    if section_x == 0.0:
                        section_x = 0.01
                    elif section_x == 1.0:
                        section_x = 0.99

                    # !!! Parameters for the tmGlut should be possible to set in the
                    # input specification !!!
                    # syn = self.sim.neuron.h.tmGlut(section(sectionX))
                    syn = channel_module(section(section_x))
                    nc = h.NetCon(vs, syn)

                    nc.delay = 0.0
                    # Should weight be between 0 and 1, or in microsiemens?
                    nc.weight[0] = neuron_input["conductance"][()] * 1e6  # !! what is unit? microsiemens?
                    nc.threshold = 0.1

                    # Get the modifications of synapse parameters, specific to
                    # this synapse
                    if param_list is not None and len(param_list) > 0:
                        syn_params = param_list[paramID % len(param_list)]["synapse"]

                        for par in syn_params:
                            if par == "expdata":
                                # Not a parameter
                                continue

                            if par == "cond":
                                # Ignoring cond value specified for synapse, using the
                                # one specified in the input information instead
                                continue

                            eval_str = "syn." + par + "=" + str(syn_params[par])
                            # self.writeLog("Updating synapse: " + evalStr)
                            # !!! Can we avoid an eval here, it is soooo SLOW
                            exec(eval_str)

                    # !!! Set parameters in synParams

                    # Need to save references, otherwise they will be freed
                    # So sorry, but that is how neuron is
                    self.external_stim[neuron_id].append((v, vs, nc, syn, spikes))

                    # ps = h.PatternStim()

                    # HOW DO WE USE PATTERNSTIM?

    ############################################################################

    def set_resting_voltage(self, neuron_id, rest_volt=None):

        if rest_volt is None:
            # If no resting voltage is given, extract it from parameters
            rest_volt = [x for x in self.neurons[neuron_id].parameters if x["param_name"] == "v_init"][0]["value"]
            self.write_log("Neuron " + self.neurons[neuron_id].name + " resting voltage = " + str(rest_volt))

        soma = [x for x in self.neurons[neuron_id].icell.soma]
        axon = [x for x in self.neurons[neuron_id].icell.axon]
        dend = [x for x in self.neurons[neuron_id].icell.dend]

        cell = soma + axon + dend

        for sec in cell:
            for seg in sec.allseg():
                seg.v = rest_volt

    ############################################################################

    def add_virtual_neuron_input(self):

        self.write_log("Adding inputs from virtual neurons")

        assert False, "addVirtualNeuronInput not implemented"

    ############################################################################

    def centre_neurons(self, side_len=None, neuron_id=None):
        if neuron_id is None:
            neuron_id = self.neuron_id

        if side_len is None:
            return neuron_id

        c_id = []

        positions = self.network_info["neuronPositions"]

        centre_pos = np.min(positions, axis=0)

        for nid in neuron_id:
            # pos = self.network_info["neurons"][nid]["position"]
            pos = positions[nid, :]

            if (abs(pos[0] - centre_pos[0]) <= side_len
                    and abs(pos[1] - centre_pos[1]) <= side_len
                    and abs(pos[2] - centre_pos[2]) <= side_len):
                c_id.append(nid)

        print("Centering: Keeping " + str(len(c_id)) + "/" + str(len(neuron_id)))

        return c_id

    ############################################################################

    def add_recording_of_type(self, neuron_type, num_neurons=None):

        cell_id = self.snudda_loader.get_cell_id_of_type(neuron_type=neuron_type,
                                                         num_neurons=num_neurons)

        self.add_recording(cell_id)

    ############################################################################

    def add_voltage_clamp(self, cell_id, voltage, duration, res=1e-9, save_i_flag=False):

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
                self.write_log(tstr)
                import pdb
                pdb.set_trace()

            self.write_log("Adding voltage clamp to " + str(cID))
            s = self.neurons[cID].icell.soma[0]
            vc = neuron.h.SEClamp(s(0.5))
            vc.rs = rs
            vc.amp1 = v * 1e3
            vc.dur1 = dur * 1e3

            self.write_log("Resistance: " + str(rs) + ", voltage: " + str(vc.amp1) + "mV")

            self.v_clamp_list.append(vc)

            if save_i_flag:
                cur = self.sim.neuron.h.Vector()
                cur.record(vc._ref_i)
                self.i_save.append(cur)
                self.i_key.append(cID)

    ############################################################################

    def add_recording(self, cell_id=None, side_len=None):
        self.write_log("Adding somatic recordings")

        if cell_id is None:
            cell_id = self.neuron_id

        # Does nothing if sideLen is not specified (otherwise, give neurons in
        # the centre)
        cell_id = self.centre_neurons(side_len=side_len, neuron_id=cell_id)

        # Only include neuron IDs on this worker, ie those in self.neuronID
        # (filtering in the if statement)
        cells = dict((k, self.neurons[k])
                     for k in cell_id if (not self.is_virtual_neuron[k]
                                          and k in self.neuron_id))

        if len(self.t_save) == 0 or self.t_save is None:
            self.t_save = self.sim.neuron.h.Vector()
            self.t_save.record(self.sim.neuron.h._ref_t)

        for cellKey in cells:
            cell = cells[cellKey]
            try:
                v = self.sim.neuron.h.Vector()
                # import pdb
                # pdb.set_trace()
                v.record(getattr(cell.icell.soma[0](0.5), '_ref_v'))
                self.v_save.append(v)
                self.v_key.append(cellKey)
            except Exception as e:
                self.write_log("Error: " + str(e))
                import pdb
                pdb.set_trace()

    ############################################################################

    def run(self, t=1000.0, hold_v=None):

        self.setup_print_sim_time(t)

        start_time = timeit.default_timer()

        # If we want to use a non-default initialisation voltage, we need to
        # explicitly set: h.v_init
        # self.sim.neuron.h.v_init = -78
        # self.sim.neuron.h.finitialize(-78)
        if hold_v is None:
            self.sim.neuron.h.finitialize()
        else:
            self.write_log("User override for holding voltage: " + str(hold_v * 1e3) + " mV")
            self.sim.neuron.h.finitialize(hold_v * 1e3)

        # Asked on neuron, check answer:
        # https://www.neuron.yale.edu/phpBB/viewtopic.php?f=2&t=4161&p=18021

        # Make sure all processes are synchronised
        self.pc.barrier()
        self.write_log("Running simulation for " + str(t / 1000) + " s")
        # self.sim.psolve(t)
        self.sim.run(t, dt=0.025)
        self.pc.barrier()
        self.write_log("Simulation done.")

        end_time = timeit.default_timer()
        self.write_log("Simulation run time: " + str(end_time - start_time) + " s")

    ############################################################################

    def plot(self):
        import matplotlib.pyplot as pyplot
        pyplot.figure()
        for v in self.v_save:
            pyplot.plot(self.t_save, v)

        pyplot.xlabel('Time (ms)')
        pyplot.ylabel('Voltage (mV)')
        pyplot.show()
        from os.path import basename
        name = basename(self.network_info["configFile"])
        pyplot.savefig('figures/Network-voltage-trace' + name + '.pdf')

    ############################################################################

    # def getSpikes(self):
    #  spiketrain.netconvecs_to_listoflists(self.tSpikes,self.idSpikes)

    ############################################################################

    def write_spikes(self, output_file=None):

        if output_file is None:
            output_file = self.get_spike_file_name()

        self.write_log("Writing spike times to " + output_file)

        for i in range(int(self.pc.nhost())):
            self.pc.barrier()  # sync all processes
            if i == int(self.pc.id()):
                if i == 0:
                    mode = 'w'
                else:
                    mode = 'a'
                with open(output_file, mode) as spikeFile:
                    for (t, id) in zip(self.t_spikes, self.id_spikes):
                        spikeFile.write('%.3f\t%d\n' % (t, id))
            self.pc.barrier()

    ############################################################################

    # secList is a list of sections
    # secXList is a list of X values 0 to 1.0
    # destID is the ID of the neuron receiving synapse (one value!)
    # voxel coords are the voxel that the synapse is in

    # We want to check that voxel coords transformed to local coordinate system
    # of neuron matches with where neuron places the synapse

    def verify_synapse_placement(self, sec_list, sec_x_list, dest_id, voxel_coords):

        # print("Running verify synapse placement")

        simulation_origo = self.network_info["simulationOrigo"]
        voxel_size = self.network_info["voxelSize"]
        neuron_position = self.network_info["neurons"][dest_id]["position"]
        neuron_rotation = self.network_info["neurons"][dest_id]["rotation"]

        # Transform voxel coordinates to local neuron coordinates to match neuron
        synapse_pos = (voxel_size * voxel_coords + simulation_origo - neuron_position) * 1e6

        syn_pos_nrn = np.zeros((len(sec_list), 3))

        for i, (sec, secX) in enumerate(zip(sec_list, sec_x_list)):
            num_points = h.n3d(sec=sec)
            arc_len = h.arc3d(num_points - 1, sec=sec)
            idx = int(np.round(secX * (num_points - 1)))
            arc_len_x = h.arc3d(idx, sec=sec)

            # print("X : " + str(secX) + " = " + str(arcLenX/arcLen) + " ???")

            syn_pos_nrn[i, 0] = h.x3d(idx, sec=sec)
            syn_pos_nrn[i, 1] = h.y3d(idx, sec=sec)
            syn_pos_nrn[i, 2] = h.z3d(idx, sec=sec)

        # We need to rotate the neuron to match the big simulation
        # !!! OBS, this assumes that some is in 0,0,0 local coordinates
        syn_pos_nrn_rot = np.transpose(np.matmul(neuron_rotation,
                                                 np.transpose(syn_pos_nrn)))

        syn_mismatch = np.sqrt(np.sum((syn_pos_nrn_rot - synapse_pos) ** 2, axis=1))

        bad_threshold = 50
        num_bad = np.sum(syn_mismatch > bad_threshold)

        if num_bad > 0:
            # If this happens, check that Neuron does not warn for removing sections
            # due to having only one point
            self.write_log("!!! Found " + str(num_bad) + " synapses on "
                           + self.network_info["neurons"][dest_id]["name"]
                           + "( " + str(dest_id) + ") "
                           + " that are further than " + str(bad_threshold) + "mum away."
                           + " morphology: "
                           + self.network_info["neurons"][dest_id]["morphology"])

            ### DEBUG PLOT!!!

            if True:
                import matplotlib.pyplot as plt
                plt.figure()

                soma_dist = np.sqrt(np.sum(synapse_pos ** 2, axis=1))
                plt.scatter(soma_dist * 1e6, syn_mismatch)
                plt.ion()
                plt.show()
                plt.title(self.network_info["neurons"][dest_id]["name"])

                from mpl_toolkits.mplot3d import Axes3D
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')

                ax.scatter(synapse_pos[:, 0],
                           synapse_pos[:, 1],
                           synapse_pos[:, 2], color="red")
                ax.scatter(syn_pos_nrn_rot[:, 0],
                           syn_pos_nrn_rot[:, 1],
                           syn_pos_nrn_rot[:, 2], color="black", s=50)

                if False:
                    # Draw neuron
                    all_sec = [x for x in neuron.h.allsec() if "axon" not in str(x)]
                    for x in np.linspace(0, 1, 10):
                        sec_pos = np.array([[h.x3d(x, sec=sec),
                                            h.y3d(x, sec=sec),
                                            h.z3d(x, sec=sec)]
                                           for sec in all_sec])

                        ax.scatter(sec_pos[:, 0], sec_pos[:, 1], sec_pos[:, 2], color="blue")

                import pdb
                pdb.set_trace()

        # voxelCoords *

        # Get local neuron position
        # self.neurons["position"]

        # for sec,secX
        # h.x3d(

    ############################################################################

    # File format for csv voltage file:
    # -1,t0,t1,t2,t3 ... (time)
    # cellID,v0,v1,v2,v3, ... (voltage for cell #ID)
    # repeat

    def write_voltage(self,
                      output_file="save/traces/network-voltage",
                      down_sampling=20):
        for i in range(int(self.pc.nhost())):
            self.pc.barrier()

            if i == int(self.pc.id()):
                if i == 0:
                    mode = 'w'
                else:
                    mode = 'a'

                with open(output_file, mode) as voltageFile:
                    if mode == 'w':
                        voltageFile.write('-1')  # Indiciate that first column is time

                        for tIdx in range(0, len(self.t_save), down_sampling):
                            voltageFile.write(',%.4f' % self.t_save[tIdx])

                    for vID, voltage in zip(self.v_key, self.v_save):
                        voltageFile.write('\n%d' % vID)

                        for vIdx in range(0, len(voltage), down_sampling):
                            voltageFile.write(',%.4f' % voltage[vIdx])

            self.pc.barrier()

    ############################################################################

    # File format for csv current file:
    # -1,t0,t1,t2,t3 ... (time)
    # cellID,i0,i1,i2,i3, ... (current for cell #ID)
    # repeat

    def write_current(self,
                      output_file="save/traces/network-current",
                      down_sampling=20):

        for i in range(int(self.pc.nhost())):
            self.pc.barrier()

            if i == int(self.pc.id()):
                if i == 0:
                    mode = 'w'
                else:
                    mode = 'a'

                with open(output_file, mode) as current_file:
                    if mode == 'w':
                        current_file.write('-1')  # Indiciate that first column is time

                        for tIdx in range(0, len(self.t_save), down_sampling):
                            current_file.write(',%.4f' % self.t_save[tIdx])

                    for iID, cur in zip(self.i_key, self.i_save):
                        current_file.write('\n%d' % iID)

                        for iIdx in range(0, len(cur), down_sampling):
                            current_file.write(',%.4f' % cur[iIdx])

            self.pc.barrier()

    ##############################################################################

    def write_log(self, text, flush=True):
        if self.log_file is not None:
            self.log_file.write(text + "\n")
            print(text)
            if flush:
                self.log_file.flush()
        else:
            if self.verbose:
                print(text)

    ############################################################################

    def create_dir(self, dir_name):
        if not os.path.isdir(dir_name):
            print("Creating " + str(dir_name))
            try:
                os.makedirs(dir_name)
            except:
                print("Failed to create dir. Already exists?")

    ############################################################################

    def add_current_injection(self, neuron_id, start_time, end_time, amplitude):

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

    def get_spike_file_name(self):

        spike_file = os.path.basename(self.network_file) + "/simulation/spike-data.txt"
        return spike_file

    ############################################################################

    def get_volt_file_name(self):

        volt_file = os.path.basename(self.network_file) + "/simulation/simulation-volt.txt"

        return volt_file

    ############################################################################

    # Use event handling

    def setup_print_sim_time(self, t_max):

        # Only have the first node print time estimates
        if self.pc.id() == 0:
            self.t_max = t_max
            self.sim_start_time = timeit.default_timer()
            self.fih_time = h.FInitializeHandler((self._setup_print_sim_time_helper, t_max))

    ############################################################################

    def _setup_print_sim_time_helper(self, t_max):
        update_points = np.arange(t_max / 100., t_max, t_max / 100.)
        for t in update_points:
            h.cvode.event(t, self.print_sim_time)

    ############################################################################

    def print_sim_time(self):
        cur_time = timeit.default_timer()
        elapsed_time = cur_time - self.sim_start_time
        fraction_done = h.t / self.t_max
        time_left = elapsed_time * ((self.t_max - h.t) / h.t)

        self.write_log("%.0f%% done. Elapsed: %.1f s, estimated time left: %.1f s"
                       % (fraction_done * 100, elapsed_time, time_left))

    ############################################################################

    def check_memory_status(self, threshold=0.1):

        mem_available, mem_total = snudda.utils.memory.memory_status()

        memory_ratio = mem_available / mem_total

        self.write_log(str(self.pc.id()) + ": Memory status: "
                       + str(int(memory_ratio * 100)) + "% free")

        return memory_ratio < threshold

    ############################################################################

    def set_dopamine_modulation(self, sec, transient_vector=None):

        if not transient_vector:
            transient_vector = []

        channel_list = {'spn': ['naf_ms', 'kas_ms', 'kaf_ms', 'kir_ms', 'cal12_ms', 'cal13_ms', 'can_ms', 'car_ms'],
                        'fs': ['kir_fs', 'kas_fs', 'kaf_fs', 'naf_fs'],
                        'chin': ['na_ch', 'na2_ch', 'kv4_ch', 'kir2_ch', 'hcn12_ch', 'cap_ch'],
                        'lts': ['na3_lts', 'hd_lts']}

        for cell_type in channel_list:
            for seg in sec:
                for mech in seg:
                    if mech.name in channel_list[cell_type]:
                        if len(transient_vector) == 0:
                            mech.damod = 1
                        else:
                            transient_vector.play(mech._ref_damod,
                                                  self.sim.neuron.h.dt)

    ############################################################################

    def apply_dopamine(self, cell_id=None, transient_vector=None):

        if not transient_vector:
            transient_vector = []

        if cell_id is None:
            cell_id = self.neuron_id

        cells = dict((k, self.neurons[k]) for k in cell_id if not self.is_virtual_neuron[k])

        for c in cells.values():
            for comp in [c.icell.dend, c.icell.axon, c.icell.soma]:
                for sec in comp:
                    self.set_dopamine_modulation(sec, transient_vector)

    ############################################################################

    def get_path(self, path_str):

        return path_str.replace("$DATA", os.path.dirname(__file__) + "/data")

    ############################################################################


def find_latest_file(file_mask):
    files = glob(file_mask)

    mod_time = [os.path.getmtime(f) for f in files]
    idx = np.argsort(mod_time)

    return files[idx[-1]]


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
    parser.add_argument("--spikesOut", "--spikesout",
                        default=None,
                        help="Name of spike output file (csv)")
    parser.add_argument("--voltOut", "--voltout",
                        default=None,
                        help="Name of voltage output file (csv)")
    parser.add_argument("--disableGJ", action="store_true",
                        help="Disable gap junctions")
    parser.add_argument("--time", type=float, default=1.5,
                        help="Duration of simulation in seconds")
    parser.add_argument("--verbose", action="store_true")

    # If called through "nrniv -python Network_simulate.py ..." then argparse
    # gets confused by -python flag, and we need to ignore it
    # parser.add_argument("-python",help=argparse.SUPPRESS,
    #                    action="store_true")

    args = parser.parse_args()
    network_data_file = args.networkFile
    input_file = args.inputFile
    log_file = os.path.dirname(args.networkFile) + "/network-simulation-log.txt"

    save_dir = os.path.dirname(args.networkFile) + "/simulation/"

    if not os.path.exists(save_dir):
        print("Creating directory " + save_dir)
        os.makedirs(save_dir, exist_ok=True)

    # Get the SlurmID, used in default file names
    slurm_id = os.getenv('SLURM_JOBID')

    if slurm_id is None:
        digits = re.findall(r'\d+', input_file)
        # Second to last digit is slurmID of old run, reuse that
        try:
            slurm_id = digits[-2]
        except:
            print("Failed to auto detect SlurmID, defaulting to 666")
            slurm_id = str(666)

    if args.voltOut is None:
        # Do not save neuron soma voltage
        volt_file = None
    else:
        # Save neuron voltage
        if args.voltOut == "default":
            volt_file = save_dir + 'network-voltage-' + slurm_id + '.csv'
        else:
            volt_file = args.voltOut

    if args.spikesOut is None or args.spikesOut == "default":
        spikes_file = save_dir + 'network-output-spikes-' + slurm_id + '.txt'
    else:
        spikes_file = args.spikesOut

    start = timeit.default_timer()

    disableGJ = args.disableGJ
    # assert disableGJ, "Please use --disableGJ for now, need to test code"
    if disableGJ:
        print("!!! WE HAVE DISABLED GAP JUNCTIONS !!!")

    pc = h.ParallelContext()

    sim = SnuddaSimulate(network_file=network_data_file,
                         input_file=input_file,
                         disable_gap_junctions=disableGJ,
                         log_file=log_file,
                         verbose=args.verbose)

    sim.add_external_input()
    sim.check_memory_status()

    if volt_file is not None:
        sim.add_recording(side_len=None)  # Side len let you record from a subset
        # sim.addRecordingOfType("dSPN",5) # Side len let you record from a subset
        # sim.addRecordingOfType("dSPN",2)
        # sim.addRecordingOfType("iSPN",2)
        # sim.addRecordingOfType("FSN",2)
        # sim.addRecordingOfType("LTS",2)
        # sim.addRecordingOfType("ChIN",2)

    tSim = args.time * 1000  # Convert from s to ms for Neuron simulator

    sim.check_memory_status()
    print("Running simulation for " + str(tSim) + " ms.")
    sim.run(tSim)  # In milliseconds

    print("Simulation done, saving output")
    if spikes_file is not None:
        sim.write_spikes(spikes_file)

    if volt_file is not None:
        sim.write_voltage(volt_file)

    stop = timeit.default_timer()
    if sim.pc.id() == 0:
        print("Program run time: " + str(stop - start))

    # sim.plot()
    exit(0)

# Check this code example
# Why are spikes not propagated from one neuron to another
# https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=188544&file=%2FLyttonEtAl2016%2FREADME.html#tabs-2
