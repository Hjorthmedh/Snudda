# snudda_prune.py
#
# Johannes Hjorth, Royal Institute of Technology (KTH)
# Human Brain Project 2019
#
# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union’s Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).
#
#
# Currently a n-way merge sort is done in serial using a heap
#
# New idea for how to speed up the merge:
#
# Divide neurons into groups depending on spatial location.
# We have a list with which hyper voxel each neuron is in.
# Each worker opens their files, and does an n-way merge
# sort ignoring all other neurons but their own
# Finally, the resulting files from the workers are merged together

import numpy as np
import scipy
import math
import os
from glob import glob
import collections

import heapq  # Priority queue

import time
import timeit

import h5py
import json
import pickle

# from snudda.Neuron_morphology import NeuronMorphology


# The x,y,z coordinates of the synapses in the original file are integer
# indexes referring to within the hyper voxel. Here they are transformed
# to refer to within the entire simulation, with the simulation origo as
# the origo (bin size unchanged)

# This code combines multiple voxel files
# !!! This also performs the pruning, since it requires knowledge of all
#     synapses between the neurons, something not possible within a hyper voxel
#     if the neuron crosses borders

class SnuddaPrune(object):

    ############################################################################

    def __init__(self,
                 network_path,
                 logfile=None, logfile_name=None,
                 d_view=None, lb_view=None, role="master", verbose=True,
                 config_file=None,  # Default is to use same config_file as for detect, but you can override it
                 scratch_path=None,
                 # pre_merge_only=False,
                 h5libver="latest",
                 random_seed=None,
                 clean_voxel_files=True):

        # Help with parallel debugging, when workers cant print to screen:
        # self.writeToRandomFile("WH = " + str(workHistoryFile) \
        #                       + "\nlog = " + str(logFileName) \
        #                       + "\nrole = " + role)

        self.work_history_file = os.path.join(network_path, "log", "network-detect-worklog.hdf5")    
        self.network_path = network_path

        self.logfile = logfile
        self.logfile_name = logfile_name
        self.verbose = verbose
        self.h5libver = h5libver

        if logfile is None and logfile_name is not None:
            self.logfile = open(logfile_name, 'w')
            self.write_log("Log file created.")

        self.write_log(f"Random seed: {random_seed}")
        self.random_seed = random_seed

        self.h5driver = "sec2"  # "core" # "stdio", "sec2"

        self.write_log(f"Using hdf5 driver {self.h5driver}, {self.h5libver} version")

        self.d_view = d_view
        self.lb_view = lb_view
        self.role = role
        self.workers_initialised = False

        self.synapse_type_lookup = {1: "GABA",
                                    2: "AMPA_NMDA",
                                    3: "GapJunction",
                                    4: "ACh",
                                    5: "NO"}

        self.synapse_type_reverse_lookup = \
            {v: k for k, v in self.synapse_type_lookup.items()}

        self.max_channel_type = None

        # Parameters for the HDF5 writing, this affects write speed
        self.synapse_chunk_size = 10000
        self.gap_junction_chunk_size = 10000
        self.h5compression = "lzf"

        # These are for the merge code
        self.synapse_write_buffer = None
        self.synapse_buffer_size = 100000
        self.next_buffer_write_pos = 0
        self.next_file_write_pos = 0
        self.buffer_out_file = None
        self.file_list = []
        self.file_buffers = []

        # These are for the pruning code
        self.merge_data_type = None  # "synapses" or "gapJunctions"
        self.synapse_read_buffer = None
        self.next_synapse_buffer_read_pos = 0
        self.next_synapse_file_read_pos = 0
        self.end_of_synapse_read_range = 0

        # Imported from HDF5 file
        self.slurm_id = None
        self.hyper_voxel_id_list = None
        self.all_hyper_id_list = None
        self.voxel_size = None
        self.hyper_voxel_size = None
        self.simulation_origo = None
        self.hyper_voxel_width = None
        self.axon_stump_id_flag = None
        self.hyper_voxel_offset = None
        self.num_synapses_total = None
        self.num_gap_junctions_total = None
        self.config_file = None
        self.position_file = None
        self.detect_config = None
        self.config = None

        self.population_unit_id = None
        self.connectivity_distributions = None
        self.type_id_lookup = None
        self.type_id_list = None

        self.next_file_read_pos = None
        self.next_buffer_read_pos = None

        self.hist_file = None
        self.out_file = None
        self.temp_file_list = []  # List of all temp files created
        self.scratch_path = None

        self.next_merge_file_id = 0

        self.voxel_overflow_counter = 0
        self.overflow_files = []

        self.open_work_history_file(work_history_file=self.work_history_file, config_file=config_file)

        self.set_scratch_path(scratch_path)
        self.load_pruning_information(config_file=config_file)

        # (locationOfMatrix,locationOfN,locationOfCoords)
        self.data_loc = {"synapses": ("network/synapses",
                                      "nSynapses",
                                      "network/synapseLookup"),  # range(2,5)),
                         "gapJunctions": ("network/gapJunctions",
                                          "nGapJunctions",
                                          "network/gapJunctionLookup")}  # range(6,9))}

        # TODO: Move this outside of init and into "prune"

    def prune(self, pre_merge_only=False):

        start_time = timeit.default_timer()

        if self.role == "master":

            # This bit is only done by the master, which then delegates job to
            # the worker nodes

            # For the larger networks the merge can take hours, so we
            # check if the file already exists
            (merge_done, synapses, synapse_file) = self.merge_file_exists()

            if not merge_done:

                if self.d_view:

                    self.write_log("Running parallel merge")
                    synapse_file = self.big_merge_parallel()

                    assert not (synapse_file["network/synapses"][-1, :] == 0).all(), \
                        "We are missing some synapses in the merge!"

                else:
                    self.write_log("Running merge in serial")

                    (synapses, synapse_file) = \
                        self.big_merge_lookup(merge_data_type="synapses",
                                              clean_voxel_files=True)

                    (gap_junctions, gap_junction_file) = \
                        self.big_merge_lookup(merge_data_type="gapJunctions",
                                              clean_voxel_files=True)

            # When running on a cluster, we might want to do the serial parts
            # in a separate run, hence this option that allows us to prepare
            # the data for parallel execution.
            if pre_merge_only:
                self.write_log(f"Pre-merge of synapses done. preMergeOnly = {pre_merge_only}, exiting.",
                               force_print=True)

                if self.hist_file is not None:
                    self.hist_file.close()
                    self.hist_file = None

                end_time = timeit.default_timer()

                self.write_log(f"Total duration: {end_time - start_time}")
                return

            # We need to close and reopen file with other driver, to allow swmr mode
            f_name = synapse_file.filename
            synapse_file.flush()
            synapse_file.close()
            synapse_file = h5py.File(f_name, "r", libver=self.h5libver, swmr=True)

            # Delegate pruning
            if self.d_view is not None:
                self.prune_synapses_parallel(synapse_file=synapse_file,
                                             merge_data_type="synapses",
                                             close_input_file=False,
                                             setup_out_file=True)

                self.prune_synapses_parallel(synapse_file=synapse_file,
                                             merge_data_type="gapJunctions",
                                             setup_out_file=False)
            else:
                # OBS, dont close the input file after first call
                self.prune_synapses(synapse_file=synapse_file,
                                    output_filename=None, row_range=None,
                                    close_out_file=False,
                                    close_input_file=False,
                                    merge_data_type="synapses")
                self.prune_synapses(synapse_file=synapse_file,
                                    output_filename=None, row_range=None,
                                    close_out_file=False,
                                    close_input_file=True,
                                    merge_data_type="gapJunctions")

            if self.out_file is None:
                self.write_log("No output file created, no synapses exist?", force_print=True)
                self.write_log("Creating symbolic link to MERGE file instead", force_print=True)

                f_dest = os.path.join(os.path.dirname(self.work_history_file).replace(f"{os.path.sep}log", os.path.sep),
                                      "network-pruned-synapses.hdf5")
                f_src = os.path.basename(f_name)
                self.write_log(f"{f_dest} -> {f_src}", force_print=True)

                if os.path.exists(f_dest):
                    os.remove(f_dest)

                os.symlink(f_src, f_dest)

                return

            n_syn_before = np.sum(self.hist_file["nSynapses"][()])
            n_syn_after = self.out_file["network/nSynapses"][0]

            n_overflow = np.sum(self.hist_file["voxelOverflowCounter"][()])
            n_gj_before = np.sum(self.hist_file["nGapJunctions"][()])
            n_gj_after = self.out_file["network/nGapJunctions"][0]

            self.write_log(f"Voxel overflows: {n_overflow} (should be zero)", is_error=(n_overflow > 0))

            if n_syn_before > 0:
                self.write_log(f"Synapses before pruning: {n_syn_before}", force_print=True)
                self.write_log(f"Synapses after pruning: {n_syn_after}"
                               f" ({round(100.0 * n_syn_after / n_syn_before, 2)} % kept)", force_print=True)
            else:
                self.write_log("No synapses to prune", force_print=True)

            if n_gj_before > 0:
                self.write_log(f"Gap junctions before pruning {n_gj_before}", force_print=True)
                self.write_log(f"Gap junctions after pruning {n_gj_after}"
                               f" ({round(100.0 * n_gj_after / n_gj_before, 2)} % kept)", force_print=True)
            else:
                self.write_log("No gap junctions to prune.", force_print=True)

            self.clean_up_merge_files()

            if self.hist_file is not None:
                self.hist_file.close()
                self.hist_file = None

            end_time = timeit.default_timer()

            self.write_log(f"Total duration: {end_time - start_time}")

            if self.voxel_overflow_counter > 0:
                self.write_log(f"Voxel overflows: {self.voxel_overflow_counter}\nIn files: ", is_error=True)
                for f in self.overflow_files:
                    self.write_log(f"Overflow in {f}", is_error=True)

    ############################################################################

    def set_scratch_path(self, scratch_path=None):

        assert self.work_history_file is not None and self.work_history_file is not "last", \
            "Need to call openWorkHistoryFile before setScratchPath"

        if scratch_path is None:
            self.scratch_path = os.path.join(self.network_path, "temp")

            if not os.path.exists(self.scratch_path):
                os.makedirs(self.scratch_path)

            self.write_log(f"Using default scratch path: {self.scratch_path}")
        else:
            self.scratch_path = scratch_path
            self.write_log(f"User selected scratch path: {self.scratch_path}")

    ############################################################################

    def __del__(self):

        try:
            if self.hist_file is not None:
                self.hist_file.close()
        except:
            print("Hist file already closed?")

        try:
            if self.out_file is not None:
                self.out_file.close()
                self.out_file = None
        except:
            print("Out file already closed?")

        # self.clean_up_merge_files()  # -- This caused old files to be cleaned up when aborting. Bad for debugging.

    ############################################################################

    def open_work_history_file(self, work_history_file=None, config_file=None):

        if work_history_file is None:
            work_history_file = self.work_history_file

        if work_history_file == "last":
            work_history_file = self.find_latest_file()

        self.write_log(f"Opening work history file: {work_history_file}")

        self.hist_file = h5py.File(work_history_file, 'r')
        self.work_history_file = work_history_file

        self.slurm_id = self.hist_file["meta/SlurmID"][()]
        self.hyper_voxel_id_list = self.hist_file["meta/hyperVoxelIDs"][()]
        self.all_hyper_id_list = self.hist_file["allHyperIDs"][()]
        self.voxel_size = self.hist_file["meta/voxelSize"][()]
        self.hyper_voxel_size = self.hist_file["meta/hyperVoxelSize"][()]  # num bins
        self.simulation_origo = self.hist_file["meta/simulationOrigo"][()]
        self.hyper_voxel_width = self.voxel_size * self.hyper_voxel_size

        # We need to make sure that detect finished correctly, that all hyper voxels are done
        all_id = set(self.hist_file["allHyperIDs"])
        n_completed = int(self.hist_file["nCompleted"][0])
        completed_id = set(self.hist_file["completed"][:n_completed])
        remaining = all_id - completed_id
        assert len(remaining) == 0, (f"Detection not done. There are {len(remaining)} hypervoxels "
                                     f"not completed: {', '.join([str(x) for x in remaining])}")

        # Network_simulate.py uses axonStumpIDFlag = True
        # Neurodamus uses axonStumpIDFlag = False
        self.axon_stump_id_flag = self.hist_file["meta/axonStumpIDFlag"]

        # We need a lookup table for offsets of hypervoxel origos
        self.hyper_voxel_offset = np.zeros((self.hyper_voxel_id_list.size, 3), dtype=int)
        for ix in range(0, self.hyper_voxel_id_list.shape[0]):
            for iy in range(0, self.hyper_voxel_id_list.shape[1]):
                for iz in range(0, self.hyper_voxel_id_list.shape[2]):
                    self.hyper_voxel_offset[self.hyper_voxel_id_list[ix, iy, iz], :] \
                        = np.array([ix, iy, iz]) * self.hyper_voxel_size

        # OBS the synapse and gap junction numbers are listed not in order of
        # hyperID but in order they were completed, to find out hyperID for
        # a particular one check self.histFile["completed"]
        self.num_synapses_total = np.sum(self.hist_file["nSynapses"][()])
        self.num_gap_junctions_total = np.sum(self.hist_file["nGapJunctions"][()])

        if config_file is None:
            self.config_file = self.hist_file["meta/configFile"][()]
        else:
            self.config_file = config_file
        self.position_file = self.hist_file["meta/positionFile"][()]

        self.detect_config = json.loads(self.hist_file["meta/config"][()])  # This was config data used for detection, might differ from pruning config
        with open(self.config_file, "r") as f:
            self.config = json.load(f)

        # This also loads random seed from config file while we have it open
        if self.random_seed is None:
            if "RandomSeed" in self.config and "prune" in self.config["RandomSeed"]:
                self.random_seed = self.config["RandomSeed"]["prune"]
                self.write_log(f"Reading random seed from config file: {self.random_seed}")
            else:
                # No random seed given, invent one
                self.random_seed = 1003
                self.write_log(f"No random seed provided, using: {self.random_seed}")
        else:
            self.write_log(f"Using random seed provided by command line: {self.random_seed}")

        ############################################################################

    def check_hyper_voxel_integrity(self, hypervoxel_file, hypervoxel_file_name, verbose=False):

        if verbose:
            self.write_log(f"Checking that {hypervoxel_file_name} matches circuit settings")

        check_list = ["voxelSize", "hyperVoxelSize", "simulationOrigo",
                      "configFile", "positionFile", "SlurmID", "axonStumpIDFlag"]

        # Just some sanity checks
        for c in check_list:
            test = self.hist_file["meta/" + c][()] == hypervoxel_file["meta/" + c][()]
            if type(test) == bool:
                assert test, \
                    "Mismatch of " + c + " in file " + hypervoxel_file_name
            else:
                assert test.all(), \
                    "Mismatch of " + c + " in file " + hypervoxel_file_name

                # Get xyz coordinates of hyper voxel
        xyz = np.where(self.hyper_voxel_id_list == hypervoxel_file["meta/hyperVoxelID"][()])
        xyz = np.array([x[0] for x in xyz])

        # Just do a sanity check that the hypervoxel origo matches stored value
        hypervoxel_origo = self.simulation_origo + self.hyper_voxel_width * xyz
        assert (hypervoxel_origo == hypervoxel_file["meta/hyperVoxelOrigo"][()]).all(), \
            "Hyper voxel origo mismatch in file " + hypervoxel_file_name

        ofc = hypervoxel_file["meta/voxelOverflowCounter"][()]

        if ofc > 0:
            self.voxel_overflow_counter += ofc
            self.overflow_files.append(hypervoxel_file_name)
            self.write_log(f"Overflow of {ofc} in {hypervoxel_file_name}", is_error=True)

    ############################################################################

    def open_hyper_voxel(self, hyper_voxel_id, verbose=False, verify=True):

        if verbose:
            self.write_log(f"Reading hypervoxel {hyper_voxel_id}")

        h_file_name = os.path.join(self.network_path, "voxels", f"network-putative-synapse-{hyper_voxel_id}.hdf5")
        h_file = h5py.File(h_file_name)

        # Just make sure the data we open is OK and match the other data
        if verify:
            self.check_hyper_voxel_integrity(h_file, h_file_name, verbose=verbose)

        return h_file

    ############################################################################

    # !!! Cache the results after creation, in case we want to rerun pruning

    def create_connection_matrix(self):

        if self.hist_file is None:
            self.open_work_history_file()

        # Do we have the matrix cached?
        c_mat = self.load_connection_matrix_cache()
        if c_mat is not None:
            return c_mat

        num_neurons = self.hist_file["network/neurons/neuronID"].shape[0]

        self.write_log(f"Creating the {num_neurons} x {num_neurons} sparse connection matrix")

        # After creation, convert to csr or csc matrix
        con_mat = scipy.sparse.lil_matrix((num_neurons, num_neurons), dtype=np.int16)

        self.write_log(f"Parsing {len(self.all_hyper_id_list)} hyper voxel files")

        for ctr, h_id in enumerate(self.all_hyper_id_list):
            if ctr % 1000 == 0 and ctr > 0:
                self.write_log(f"{ctr} files done")

            h_file = self.open_hyper_voxel(h_id)

            synapses = h_file["network/synapses"][()]

            for row in synapses:
                src_id = row[0]
                dest_id = row[1]

                con_mat[src_id, dest_id] += 1

            h_file.close()

        self.write_log("Converting to CSR sparse matrix format")
        c_mat = con_mat.tocsr()

        self.save_connection_matrix_cache(c_mat)

        return c_mat

    ############################################################################

    def get_con_mat_cache_filename(self):
        cache_file = os.path.join(self.network_path, "connection-matrix-cache.pickle")
        return cache_file

    ############################################################################

    def load_connection_matrix_cache(self):

        cache_file = self.get_con_mat_cache_filename()

        if os.path.isfile(cache_file):
            self.write_log(f"Loading connection matrix cache from {cache_file}")

            with open(cache_file, 'rb') as f:
                data = pickle.load(f)
                con_mat = data["conMat"]

                # Verify that the matrix matches the simulation
                check_list = [("meta/SlurmID", "SlurmID"),
                              ("meta/simulationOrigo", "simulationOrigo"),
                              ("meta/voxelSize", "voxelSize"),
                              ("meta/hyperVoxelSize", "hyperVoxelSize"),
                              ("meta/configFile", "configFile"),
                              ("meta/positionFile", "positionFile"),
                              ("meta/axonStumpIDFlag", "axonStumpIDFlag")]

                for checkNames in check_list:
                    test = self.hist_file[checkNames[0]][()] == data[checkNames[1]]

                    if type(test) == bool:
                        assert test, f"Connection matrix cache - mismatch for {checkNames[1]}"
                    else:
                        assert test.all(), f"Connection matrix cache - mismatch for {checkNames[1]}"

                assert self.num_synapses_total == data["nSynapsesTotal"] \
                    and self.num_gap_junctions_total == data["nGapJunctionsTotal"], \
                        " Synapse or gap junction count mismatch -- corrupt file?"

        else:
            con_mat = None

        return con_mat

    ############################################################################

    def save_connection_matrix_cache(self, con_mat):

        cache_file = self.get_con_mat_cache_filename()
        self.write_log(f"Saving connection matrix cache to {cache_file}")

        data = dict([])
        data["conMat"] = con_mat
        data["SlurmID"] = self.slurm_id
        data["simulationOrigo"] = self.simulation_origo
        data["voxelSize"] = self.voxel_size
        data["hyperVoxelSize"] = self.hyper_voxel_size
        data["nSynapsesTotal"] = self.num_synapses_total
        data["nGapJunctionsTotal"] = self.num_gap_junctions_total
        data["configFile"] = self.config_file
        data["positionFile"] = self.position_file
        data["axonStumpIDFlag"] = self.axon_stump_id_flag

        with open(cache_file, 'wb') as f:
            pickle.dump(data, f, -1)  # -1 latest version

    ############################################################################

    # This checks that all connections included in the pruning, were present
    # in the detection. If some are missing in the detection phase, then they
    # would incorrectly be missing after pruning.

    def check_network_config_integrity(self, config_file):

        detect_config = json.loads(self.hist_file["meta/config"][()])
        with open(config_file, "r") as f:
            prune_config = json.load(f)

        all_present = True

        for con in prune_config["Connectivity"]:
            if con not in detect_config["Connectivity"]:
                self.write_log(f"!!! Connection {con} is present in {config_file}, "
                               f"but was not included in config file used for detect. " 
                               f"Please rerun snudda detect", is_error=True)
                all_present = False

        assert all_present, "Please rerun snudda detect."

    ############################################################################

    # Parse the connection information in the config file
    def load_pruning_information(self, config_file=None):

        if config_file is None:
            config_file = self.hist_file["meta/configFile"][()]

        self.check_network_config_integrity(config_file=config_file)
        with open(config_file, "r") as f:
            self.config = json.load(f)

        self.population_unit_id = self.hist_file["network/neurons/populationUnitID"][()]

        # Normally we use type names as lookups, but since we will do this
        # many millions of times, we create an temporary typeID number
        self.make_type_numbering()

        orig_connectivity_distributions = \
            json.loads(self.hist_file["meta/connectivityDistributions"][()])

        config_connectivity_distributions = self.config["Connectivity"]

        self.connectivity_distributions = dict([])

        # For the pruning we merge the two into one
        for key in config_connectivity_distributions:
            (pre_type, post_type) = key.split(",")  # key.split("$$") -- if we instead loop over orig_connectivity_distribution
            orig_key = f"{pre_type}$${post_type}"

            # Need to handle if preType or postType don't exist, then skip this
            if pre_type not in self.type_id_lookup or post_type not in self.type_id_lookup:
                self.write_log(f"Skipping {pre_type} to {post_type} connection")
                continue

            pre_type_id = self.type_id_lookup[pre_type]
            post_type_id = self.type_id_lookup[post_type]

            for con_type in config_connectivity_distributions[key]:
                con_data = config_connectivity_distributions[key][con_type]

                pruning = self.complete_pruning_info(con_data["pruning"])

                if "pruningOther" in con_data:
                    pruning_other = self.complete_pruning_info(con_data["pruningOther"])
                else:
                    pruning_other = None

                # This data is added by detect, we need to take it from what was used during detection
                synapse_type_id = orig_connectivity_distributions[orig_key][con_type]["channelModelID"]

                self.connectivity_distributions[pre_type_id, post_type_id, synapse_type_id] \
                    = (pruning, pruning_other)

    ############################################################################

    # This makes sure all the variables exist, that way pruneSynapseHelper
    # does not have to check, but can assume that they will exist

    @staticmethod
    def complete_pruning_info(prune_info):

        if "distPruning" not in prune_info:
            prune_info["distPruning"] = None

        if "f1" not in prune_info or prune_info["f1"] is None:
            prune_info["f1"] = 1.0

        if "softMax" not in prune_info:
            prune_info["softMax"] = None

        if "mu2" not in prune_info:
            prune_info["mu2"] = None

        if "a3" not in prune_info:
            prune_info["a3"] = None

        return prune_info

    ############################################################################

    def make_type_numbering(self):

        name_list = [x if type(x) not in [bytes, np.bytes_] else x.decode()
                     for x in self.hist_file["network/neurons/name"]]

        type_name_list = [x.split("_")[0] for x in name_list]

        type_ctr = 1
        self.type_id_lookup = dict([])

        for x in type_name_list:
            if x not in self.type_id_lookup:
                self.type_id_lookup[x] = type_ctr
                type_ctr += 1

        self.type_id_list = [self.type_id_lookup[x] for x in type_name_list]

    ############################################################################

    def setup_output_file(self, output_file=None):

        if self.out_file is not None:
            self.write_log(f"Output file already set: {self.out_file.filename}")
            return

        if output_file is None:
            output_file = os.path.join(self.network_path, "network-pruned-synapses.hdf5")

        self.write_log(f"Writing to {output_file}")
        out_file = h5py.File(output_file, "w", libver=self.h5libver, driver=self.h5driver)
        out_file.create_dataset("config", data=json.dumps(self.config))

        # Copy over meta data
        self.hist_file.copy("meta", out_file)

        cfg = json.loads(self.hist_file["meta/config"][()])
        morph_group = out_file.create_group("morphologies")

        for name, definition in cfg["Neurons"].items():
            morph_file = definition["morphology"]

            with open(morph_file, "r") as f:
                swc_data = f.read()

            self.write_log(f"Saving morphology in HDF5 file: {morph_file}")
            swc_group = morph_group.create_group(name)
            swc_group.create_dataset("swc", data=swc_data)
            swc_group.create_dataset("location", data=morph_file)

        network_group = out_file.create_group("network")

        # Copy over neuron data
        self.hist_file.copy("network/neurons", network_group)

        network_group.create_dataset("synapses",
                                     dtype=np.int32,
                                     shape=(self.synapse_chunk_size, 13),
                                     chunks=(self.synapse_chunk_size, 13),
                                     maxshape=(None, 13),
                                     compression=self.h5compression)

        network_group.create_dataset("gapJunctions",
                                     dtype=np.int32,
                                     shape=(self.gap_junction_chunk_size, 11),
                                     chunks=(self.gap_junction_chunk_size, 11),
                                     maxshape=(None, 11),
                                     compression=self.h5compression)

        num_synapses = np.zeros((1,), dtype=np.uint64)
        num_gap_junctions = np.zeros((1,), dtype=np.uint64)

        network_group.create_dataset("nSynapses", data=num_synapses, dtype=np.uint64)
        network_group.create_dataset("nGapJunctions", data=num_gap_junctions,
                                     dtype=np.uint64)

        self.out_file = out_file

    ############################################################################

    def find_latest_file(self):

        # files = glob('save/network_connect_voxel_log-*-worklog.hdf5')
        files = glob(os.path.join('save","network-connect-synapse-voxel-file-*-worklog.hdf5'))

        mod_time = [os.path.getmtime(f) for f in files]
        idx = np.argsort(mod_time)

        self.write_log("Using the newest file: " + files[idx[-1]])

        return files[idx[-1]]

    ############################################################################

    def merge_file_exists(self):

        # check if merge file exists
        merge_file_name = os.path.join(self.network_path, "network-putative-synapses-MERGED.hdf5")

        merge_file_ok = False

        self.write_log(f"Checking for merge file {merge_file_name}")

        try:
            if os.path.isfile(merge_file_name):

                merge_file_ok = True

                f = h5py.File(merge_file_name, "r")

                # Check that SlurmID matches
                if f["meta/SlurmID"][()] != self.slurm_id:
                    merge_file_ok = False

                # Check that synapse matrix has right size
                if f["network/synapses"].shape[0] != self.num_synapses_total:
                    merge_file_ok = False

                # Check that last row is not empty
                if (f["network/synapses"][-1, :] == 0).all():
                    merge_file_ok = False

                if "gapJunctions" not in f["network"]:
                    merge_file_ok = False
                elif f["network/gapJunctions"].shape[0] != self.num_gap_junctions_total:
                    merge_file_ok = False
                elif self.num_gap_junctions_total > 0 and (f["network/gapJunctions"][-1, :] == 0).all():
                    # Last row not set, not complete
                    merge_file_ok = False

            else:
                f = None
        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)

            self.write_log("Something went wrong with reading old merge file, ignoring it")
            # Something went wrong
            merge_file_ok = False
            f = None

        if merge_file_ok:
            self.write_log(f"Found old merge file {merge_file_name}")
            return True, f["network/synapses"], f
        else:
            # Bad merge file, close it
            if f is not None:
                f.close()

            return False, None, None

    ############################################################################

    def setup_merge_file(self, verbose=False, big_cache=False,
                         outfile_name=None, save_morphologies=True,
                         num_synapses=None,
                         num_gap_junctions=None,
                         delete_after=True):

        if outfile_name is None:
            outfile_name = os.path.join(self.network_path, "network-putative-synapses-MERGED.hdf5")

        #  Make a list of all temporary files so we can remove them
        if delete_after:
            self.temp_file_list.append(outfile_name)

        self.write_log(f"Setting up out file {outfile_name}")
        if big_cache:
            # !!! Special code to increase h5py cache size
            propfaid = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
            settings = list(propfaid.get_cache())
            # print(settings)
            # [0, 521, 1048576, 0.75]

            settings[2] *= 20
            # propfaid.set_cache(*settings)
            # settings = propfaid.get_cache()
            # print(settings)
            # (0, 521, 5242880, 0.75)

            # Low level opening hdf5 file, to have greater cache size
            fid = h5py.h5f.create(outfile_name.encode(),
                                  flags=h5py.h5f.ACC_TRUNC,
                                  fapl=propfaid)
            out_file = h5py.File(fid, libver=self.h5libver, driver=self.h5driver)
        else:
            out_file = h5py.File(outfile_name, "w", libver=self.h5libver, driver=self.h5driver)

        network_group = out_file.create_group("network")

        # Copy over meta data
        self.hist_file.copy("meta", out_file)

        # Copy over neuron data
        # self.histFile.copy("neurons",outFile)
        self.hist_file.copy("network/neurons", network_group)

        cfg = json.loads(self.hist_file["meta/config"][()])

        # Save morphologies
        if save_morphologies:
            morph_group = out_file.create_group("morphologies")

            for name, definition in cfg["Neurons"].items():
                morph_file = definition["morphology"]

                with open(morph_file, "r") as f:
                    swc_data = f.read()

                self.write_log(f"Saving morphology in HDF5 file: {morph_file}")
                swc_group = morph_group.create_group(name)
                swc_group.create_dataset("swc", data=swc_data)
                swc_group.create_dataset("location", data=morph_file)

        chunk_size = self.synapse_chunk_size

        if num_synapses is None:
            num_synapses = self.num_synapses_total

        if num_gap_junctions is None:
            num_gap_junctions = self.num_gap_junctions_total

        # !!! We are allocating a lot more space than we might be needing
        # The matrices are resized when the file write buffer is flushed
        # in bufferMergeWrite (flush=True)

        if num_synapses > chunk_size:
            network_group.create_dataset("synapses",
                                         dtype=np.int32,
                                         shape=(num_synapses, 13),
                                         chunks=(chunk_size, 13),
                                         maxshape=(num_synapses, 13),
                                         compression=self.h5compression)
        elif num_synapses > 0:
            network_group.create_dataset("synapses",
                                         dtype=np.int32,
                                         shape=(num_synapses, 13),
                                         chunks=(num_synapses, 13))
        else:
            network_group.create_dataset("synapses",
                                         dtype=np.int32,
                                         shape=(num_synapses, 13))

        if num_gap_junctions > chunk_size:
            network_group.create_dataset("gapJunctions",
                                         dtype=np.int32,
                                         shape=(num_gap_junctions, 11),
                                         chunks=(chunk_size, 11),
                                         maxshape=(num_gap_junctions, 11),
                                         compression=self.h5compression)
        elif num_gap_junctions > 0:
            network_group.create_dataset("gapJunctions",
                                         dtype=np.int32,
                                         shape=(num_gap_junctions, 11),
                                         chunks=(num_gap_junctions, 11))
        else:
            network_group.create_dataset("gapJunctions",
                                         dtype=np.int32,
                                         shape=(num_gap_junctions, 11))

        self.next_merge_file_id += 1

        # Reset the position for the merge file
        self.next_file_write_pos = 0

        return out_file, outfile_name

    ############################################################################

    def prune_synapses_parallel(self, synapse_file,
                                output_file=None,
                                merge_data_type="synapses",
                                setup_out_file=True,
                                close_input_file=True):

        h5_syn_mat, h5_syn_n, h5_syn_loc = self.data_loc[merge_data_type]

        if synapse_file[h5_syn_mat].shape[0] == 0:
            self.write_log(f"prune_synapses_parallel: No {merge_data_type} skipping pruning", force_print=True)
            return
        else:
            self.write_log(f"prune_synapses_parallel, before pruning : "
                           f"{synapse_file[h5_syn_mat].shape[0]} {merge_data_type}")

        start_time = timeit.default_timer()

        if self.d_view is not None and self.role == "master":
            self.setup_parallel(d_view=self.d_view)

            # Make sure all synapse writes are on the disk
            synapse_file.flush()
            if not synapse_file.swmr_mode:
                synapse_file.swmr_mode = True  # Allow multiple readers from the file

        # 1. Pick names for the workers
        temp_output_file_name = [os.path.join(self.scratch_path, f"worker-temp-{merge_data_type}-file-{x}")
                                 for x in range(0, len(self.d_view))]
        self.d_view.scatter("output_filename", temp_output_file_name, block=True)

        # Add the files to a delete list, so we remove them after
        for f in temp_output_file_name:
            self.temp_file_list.append(f)

        # 2. Define what ranges each worker should do
        synapse_ranges = self.find_ranges(synapse_file[h5_syn_mat], len(self.d_view))

        if synapse_ranges is None or synapse_ranges[-1][-1] is None:
            self.write_log("There are few synapses, we will run it in serial instead", force_print=True)
            return self.prune_synapses(synapse_file=synapse_file,
                                       output_filename=None, row_range=None,
                                       close_out_file=False,
                                       close_input_file=close_input_file,
                                       merge_data_type=merge_data_type)

        # We got valid synapseRanges, continue

        self.d_view.scatter("synapse_range", synapse_ranges, block=True)
        self.write_log("synapse_ranges: " + str(synapse_ranges))

        self.d_view.push({"synapse_filename": synapse_file.filename}, block=True)
        self.d_view.push({"merge_data_type": merge_data_type}, block=True)

        # Close the synapse file on the master node
        if close_input_file:
            synapse_file.close()
            synapse_file = None

        # 3. Let the workers prune
        self.write_log("Sending pruning job to workers")

        cmd_str = ("nw.prune_synapses(synapse_file=synapse_filename,output_filename=output_filename[0]," 
                   "row_range=synapse_range[0],merge_data_type=merge_data_type)")

        self.d_view.execute(cmd_str, block=True)

        end_time = timeit.default_timer()
        self.write_log(f"Parallel pruning duration: {end_time - start_time}")

        # 4. Merge the resulting files -- this is easier, since synapses
        #    are already sorted in the file
        self.write_log("Merging parallel pruning results")

        if setup_out_file:
            assert self.out_file is None, "prune_synapses_parallel: Output file already setup"
            self.setup_output_file(output_file=output_file)
        else:
            # If there were no synapses, but there are gap junctions
            # then this assert could theoretically happen
            assert self.out_file is not None, "prune_synapses_parallel: Out file not set up"

        tmp_files = [h5py.File(f, 'r') for f in temp_output_file_name]

        num_syn = np.sum(f[h5_syn_mat].shape[0] for f in tmp_files)
        mat_width_all = [f[h5_syn_mat].shape[1] for f in tmp_files]

        assert (np.array(mat_width_all) == mat_width_all[0]).all(), \
            "prune_synapses_parallel: Internal error, width does not match"
        mat_width = mat_width_all[0]

        self.out_file[h5_syn_mat].resize((num_syn, mat_width))
        next_syn = 0

        for f in tmp_files:
            n = f[h5_syn_mat].shape[0]
            if n > 0:
                self.out_file[h5_syn_mat][next_syn:(next_syn + n), :] = f[h5_syn_mat][()]
                next_syn += n
            f.close()

        self.out_file["network/" + h5_syn_n][0] = next_syn

        end_time2 = timeit.default_timer()
        self.write_log(f"Parallel pruning + merging: {end_time2 - start_time}")

    ############################################################################

    # Find which ranges of the synapse matrix that each worker should take care of

    def find_ranges(self, synapses, num_workers, start_pos=0, num_syn=None):

        if num_syn is None:
            num_syn = synapses.shape[0] - start_pos

        block_size = max(1, int(math.floor(float(num_syn) / num_workers)))

        self.write_log(f"Find block ranges. From {start_pos} to {num_syn + start_pos} block size {block_size}")

        # We keep blockStart in a list, so we can remove if they overlap
        block_start = [x + start_pos for x in range(0, block_size * (num_workers + 1), block_size)]

        range_end = start_pos + num_syn

        # There can be multiple synapses between a cell pair, those must belong
        # to the same block since pruning is dependent on the number of synapses
        # between a pair of neurons
        for idx in range(1, len(block_start)):
            start_row = block_start[idx]

            # UPDATE: Now we make sure all synapses on a post synaptic cell is on the same worker
            while start_row < range_end and synapses[block_start[idx] - 1, 1] == synapses[start_row, 1]:
                start_row += 1

            block_start[idx] = start_row

        # If we have few synapses, and many workers, there might be some workers
        # assigned to the same range. Remove those, and pad with None

        block_start = [x for x in block_start if x <= range_end]
        block_start = list(collections.OrderedDict.fromkeys(block_start))

        while len(block_start) < num_workers + 1:
            block_start.append(None)

        synapse_range = [x for x in zip(block_start[:-1], block_start[1:])]

        self.write_log(f"synapse_range={synapse_range}")

        return synapse_range

    ############################################################################

    def clean_up_merge_files(self):

        if self.temp_file_list is None:
            # Nothing to do
            return

        self.write_log("Cleaning up old merge files")

        for f in self.temp_file_list:
            # self.writeLog("Removing old merge file: " + str(f))
            try:
                if os.path.exists(f):
                    os.remove(f)
            except:
                print(f"Closing of file failed: {f}")

        self.temp_file_list = None

    ############################################################################

    def write_log(self, text, flush=True, is_error=False, force_print=False):  # Change flush to False in future, debug
        try:
            if self.logfile is not None:
                self.logfile.write(f"{text}\n")
                print(text)
                if flush:
                    self.logfile.flush()
            else:
                if self.verbose or is_error or force_print:
                    print(text)
        except:
            print(text)
            print("Unable to write to log file. Is log file closed?")

    ############################################################################

    def setup_parallel(self, d_view):

        assert self.role == "master", "setupParallel: Should only be called by master node"

        if d_view is None:
            self.write_log("setupParallel called without dView, aborting.")
            return

        if self.workers_initialised:
            self.write_log("Workers already initialised.")
            return

        with d_view.sync_imports():
            from snudda.prune import SnuddaPrune  # Import on workers

        self.write_log(f"Setting up workers: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")

        # Create unique log file names for the workers
        if self.logfile_name is not None:
            engine_log_file = [f"{self.logfile_name}-{x}" for x in range(0, len(d_view))]
        else:
            engine_log_file = [[] for x in range(0, len(d_view))]

        d_view.scatter('logfile_name', engine_log_file, block=True)
        d_view.push({"network_path": self.network_path,
                     "random_seed": self.random_seed}, block=True)

        cmd_str = ("nw = SnuddaPrune(network_path=network_path, logfile_name=logfile_name[0]," 
                   "role='worker',random_seed=random_seed)")
        d_view.execute(cmd_str, block=True)

        self.write_log(f"Workers setup: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
        self.workers_initialised = True

    ########################################################################

    def clean_up_merge_read_buffers(self):

        self.merge_data_type = None

        for hid, f in enumerate(self.file_list):

            if f is not None:
                try:
                    f.close()
                except:
                    self.write_log(f"Problem closing file for HID: {hid}")

        if self.file_buffers is not None:
            for idx in range(0, len(self.file_buffers)):
                self.file_buffers[idx] = None

        self.file_buffers = None
        self.next_file_read_pos = None
        self.next_buffer_read_pos = None

    ############################################################################

    def buffer_merge_write(self, synapse_matrix_loc, synapses=None, flush=False):

        if self.synapse_write_buffer is None:
            # First time run
            if synapses is not None:
                buffer_width = synapses.shape[1]
            else:
                assert False, "bufferMergeWrite: Unable to guess width of matrix"

            self.synapse_write_buffer = np.zeros((self.synapse_buffer_size,
                                                  buffer_width),
                                                 dtype=np.int32)
            self.next_buffer_write_pos = 0

        if synapses is not None:

            # Is buffer almost full?
            if synapses.shape[0] + self.next_buffer_write_pos >= self.synapse_buffer_size:
                # Flush buffer first
                end_file_idx = self.next_file_write_pos + self.next_buffer_write_pos
                self.buffer_out_file[synapse_matrix_loc][self.next_file_write_pos:end_file_idx, :] \
                    = self.synapse_write_buffer[0:self.next_buffer_write_pos, :]

                self.next_file_write_pos += self.next_buffer_write_pos
                self.next_buffer_write_pos = 0
                # Ok, now we have clean buffer, back to normal program

            end_buf_idx = self.next_buffer_write_pos + synapses.shape[0]
            self.synapse_write_buffer[self.next_buffer_write_pos:end_buf_idx, :] = synapses
            self.next_buffer_write_pos += synapses.shape[0]

        if flush:
            self.write_log(f"Flushing {self.buffer_out_file.filename} data: {synapse_matrix_loc}")
            end_file_idx = self.next_file_write_pos + self.next_buffer_write_pos
            self.buffer_out_file[synapse_matrix_loc][self.next_file_write_pos:end_file_idx, :] \
                = self.synapse_write_buffer[0:self.next_buffer_write_pos, :]

            self.next_file_write_pos += self.next_buffer_write_pos
            self.next_buffer_write_pos = 0

            # Resize matrix to fit data
            w = self.buffer_out_file[synapse_matrix_loc].shape[1]
            self.buffer_out_file[synapse_matrix_loc].resize((end_file_idx, w))
            self.write_log(f"{synapse_matrix_loc} new size {(end_file_idx, w)}")

            # Remove buffer
            self.synapse_write_buffer = None

    ############################################################################

    def big_merge_parallel(self):

        if self.role != "master":
            self.write_log("big_merge_parallel is only run on master node, aborting")
            return

        self.write_log(f"big_merge_parallel, starting {self.role}")

        if self.d_view:
            self.setup_parallel(d_view=self.d_view)

        # Split neurons between nodes, we need the neurons to be in order
        num_neurons = self.hist_file["network/neurons/neuronID"].shape[0]
        assert num_neurons - 1 == self.hist_file["network/neurons/neuronID"][-1], \
            "neuronID should start from 0 and the end should be n-1"

        n_workers = len(self.d_view)

        neuron_ranges = []
        range_borders = np.linspace(0, num_neurons, n_workers + 1).astype(int)

        for idx in range(0, n_workers):
            neuron_ranges.append((range_borders[idx], range_borders[idx + 1]))

        assert neuron_ranges[-1][-1] == num_neurons, \
            "big_merge_parallel: Problem with neuron_ranges, last element incorrect"
        assert len(neuron_ranges) == n_workers, \
            "big_merge_parallel: Problem with neuron_ranges, bad length"

        # Send list of neurons to workers
        self.d_view.scatter("neuron_range", neuron_ranges, block=True)

        # Each worker sorts a subset of the neurons and write it to separate files
        cmd_str_syn = "merge_result_syn = nw.big_merge_helper(neuron_range=neuron_range[0],merge_data_type='synapses')"

        self.d_view.execute(cmd_str_syn, block=True)
        merge_results_syn = self.d_view["merge_result_syn"]

        # When we do scatter, it embeds the result in a list
        cmd_str_gj = "merge_result_gj = nw.big_merge_helper(neuron_range=neuron_range[0],merge_data_type='gapJunctions')"
        self.d_view.execute(cmd_str_gj, block=True)
        merge_results_gj = self.d_view["merge_result_gj"]

        # We need to sort the files in order, so we know how to add them
        self.write_log("Processing merge...")

        merge_start_syn = [x[1][0] for x in merge_results_syn]
        merge_start_gj = [x[1][0] for x in merge_results_gj]

        merge_order_syn = np.argsort(merge_start_syn)
        merge_order_gj = np.argsort(merge_start_gj)

        # We then need the file to merge them in
        (self.buffer_out_file, outFileName) = self.setup_merge_file(big_cache=False, delete_after=False)

        # Copy the data to the file
        for order, results, location in zip([merge_order_syn, merge_order_gj],
                                            [merge_results_syn, merge_results_gj],
                                            ["network/synapses", "network/gapJunctions"]):
            start_pos = 0

            for idx in order:
                data_file = results[idx][0]
                num_synapses = results[idx][2]   # This is synapses first iteration, gap junction 2nd iteration
                end_pos = start_pos + num_synapses

                if num_synapses == 0:
                    continue  # No synapses, skip this file.

                assert data_file is not None or num_synapses == 0, \
                    f"!!! Missing merge file {data_file}, internal problem"

                self.write_log(f"Extracting {location} from {data_file}")
                f_in = h5py.File(data_file, 'r')

                # Some idiot checks...
                assert f_in[location].shape[0] == num_synapses, "Internal inconsistency in number of rows stored"
                assert not (f_in[location][-1, :] == 0).all(), "Last row all zero, that should not happen"

                self.buffer_out_file[location][start_pos:end_pos, :] = f_in[location][:, :]

                f_in.close()
                start_pos = end_pos

        # We should check that the number of synapses in output file is correct
        assert self.buffer_out_file["network/synapses"].shape[0] == np.sum(self.hist_file["nSynapses"][:]), \
            "Not all synapses kept in merge, internal problem."
        assert self.buffer_out_file["network/gapJunctions"].shape[0] == np.sum(self.hist_file["nGapJunctions"][:]), \
            "Not all gap junctions kept in merge, internal problem."

        self.write_log("big_merge_parallel: done")

        return self.buffer_out_file

    ############################################################################

    def get_hyper_voxel_list(self, neuron_range):

        # We need to find all the hypervoxels that might contain synapses.
        # During touch detection we generated a list for each hypervoxel with
        # all the neurons that had any part of it within that hypervoxel.

        hv_list = []
        neuron_set = set(range(neuron_range[0], neuron_range[1]))

        for hid in self.hist_file["hyperVoxels"]:
            hv_neurons = self.hist_file["hyperVoxels"][hid]["neurons"]
            if len(neuron_set.intersection(hv_neurons)) > 0:
                hv_list.append(int(hid))

        return hv_list

    ############################################################################

    # Needs to handle both gap junctions and synapses
    # This is called using d_view for parallel execution, the code coverage algorithm does not understand that

    def big_merge_helper(self, neuron_range, merge_data_type):

        try:
            self.write_log(f"Neuron_range = {neuron_range}")

            output_filename = os.path.join(self.scratch_path,
                                           f"{merge_data_type}-for-neurons-"
                                           f"{neuron_range[0]}-to-{neuron_range[1]}-MERGE-ME.hdf5")

            self.merge_data_type = merge_data_type

            # Which hyper voxels are the neurons located in?
            hv_list = self.get_hyper_voxel_list(neuron_range)

            # Locations of data within the file
            h5_syn_mat, h5_syn_n, h5_syn_lookup = self.data_loc[merge_data_type]

            synapse_heap = []
            file_list = dict([])
            file_mat_iterator = dict([])

            # Next we need to open all the relevant files
            h_file_name_mask = os.path.join(self.network_path, "voxels", "network-putative-synapses-%s.hdf5")

            max_axon_voxel_ctr = 0
            max_dend_voxel_ctr = 0

            # !!! Special code to increase h5py cache size
            propfaid = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
            settings = list(propfaid.get_cache())
            settings[2] *= 20
            propfaid.set_cache(*settings)
            # !!! End special cache code

            n_hv = int(self.hist_file["nCompleted"][0])
            n_total = 0

            for h_id, n_syn, n_overflow in zip(self.hist_file["completed"][:n_hv],
                                               self.hist_file[h5_syn_n][:n_hv],
                                               self.hist_file["voxelOverflowCounter"][:n_hv]):
                n_total += n_syn

                if h_id not in hv_list:
                    # We only need a subset of the files
                    continue

                if n_syn > 0:
                    h_filename = h_file_name_mask % str(h_id)
                    self.write_log(f"Opening voxel file: {h_filename}")

                    # Low level opening hdf5 file, to have greater cache size #ACC_RDWR
                    fid = h5py.h5f.open(h_filename.encode(), flags=h5py.h5f.ACC_RDONLY, fapl=propfaid)
                    file_list[h_id] = h5py.File(fid, drive=self.h5driver)

                    # Verify hyper voxel
                    if merge_data_type == "synapses":
                        # No need to do it for both synapses and gap junctions, since same hyper voxels
                        self.check_hyper_voxel_integrity(file_list[h_id], h_filename.encode(), verbose=True)

                    if self.max_channel_type:
                        if self.max_channel_type != file_list[h_id]["network/maxChannelTypeID"][()]:
                            print("Investigate")
                            print(f"{self.max_channel_type} != {file_list[h_id]['network/maxChannelTypeID'][()]}")
                            # import pdb
                            # pdb.set_trace()

                        # These should be the same for all hypervoxels
                        assert self.max_channel_type == file_list[h_id]["network/maxChannelTypeID"][()], \
                            (f"max_channel_type = {self.max_channel_type} "
                             f"(differ with what is in file {file_list[h_id]['network/maxChannelTypeID'][()]})")
                    else:
                        self.max_channel_type = file_list[h_id]["network/maxChannelTypeID"][()]
                        print(f"Setting max_channel_type to {self.max_channel_type} from h_id={h_id}")

                    chunk_size = 10000
                    lookup_iterator = \
                        self.file_row_lookup_iterator_subset(h5mat_lookup=file_list[h_id][h5_syn_lookup],
                                                             min_dest_id=neuron_range[0],
                                                             max_dest_id=neuron_range[1],
                                                             chunk_size=chunk_size)

                    file_mat_iterator[h_id] \
                        = self.synapse_set_iterator(h5mat_lookup=file_list[h_id][h5_syn_lookup],
                                                    h5mat=file_list[h_id][h5_syn_mat],
                                                    chunk_size=chunk_size,
                                                    lookup_iterator=lookup_iterator)

                    syn_set, unique_id = next(file_mat_iterator[h_id], (None, None))

                    if syn_set is None:
                        # There were synapses in the hyper voxel, but none relevant to our
                        # selected files. Clear file List and fileMatIterator for this worker
                        del file_list[h_id]
                        del file_mat_iterator[h_id]
                        continue

                    # Create a heap containing the first subset of all files
                    heapq.heappush(synapse_heap, (unique_id, h_id, syn_set))

                    # This is so we can optimize the axon/dend voxelCtr and size
                    if "maxAxonVoxelCtr" in file_list[h_id]["meta"]:
                        max_axon_voxel_ctr = max(max_axon_voxel_ctr, file_list[h_id]["meta/maxAxonVoxelCtr"][()])
                    if "maxDendVoxelCtr" in file_list[h_id]["meta"]:
                        max_dend_voxel_ctr = max(max_dend_voxel_ctr, file_list[h_id]["meta/maxDendVoxelCtr"][()])

            if merge_data_type == "synapses":
                num_synapses = n_total
                num_gap_junctions = 0
            elif merge_data_type == "gapJunctions":
                num_synapses = 0
                num_gap_junctions = n_total
            else:
                assert False, f"Unknown mergeDataType {merge_data_type}"

            # Setup output file
            (self.buffer_out_file, outFileName) = self.setup_merge_file(big_cache=True,
                                                                        outfile_name=output_filename,
                                                                        save_morphologies=False,
                                                                        num_synapses=num_synapses,
                                                                        num_gap_junctions=num_gap_junctions)

            # Here we store the sorted connection matrix
            sorted_mat = self.buffer_out_file[h5_syn_mat]

            # Only save this meta data if doing the synapses call
            if max_axon_voxel_ctr > 0 and self.merge_data_type == "synapses":
                self.buffer_out_file["meta"].create_dataset("maxAxonVoxelCtr", data=max_axon_voxel_ctr)
                self.write_log(f"max_axon_voxel_ctr = {max_axon_voxel_ctr}")

            if max_dend_voxel_ctr > 0 and self.merge_data_type == "synapses":
                self.buffer_out_file["meta"].create_dataset("maxDendVoxelCtr", data=max_dend_voxel_ctr)
                self.write_log(f"max_dend_voxel_ctr = {max_dend_voxel_ctr}")

            # Take the first (smallest uniqueID) element from the heap
            if len(synapse_heap) > 0:
                # Get the first file to read synapses from
                (unique_id, h_id, syn_set) = heapq.heappop(synapse_heap)
            else:
                # No synapses at all, return
                self.clean_up_merge_read_buffers()
                return None, neuron_range, 0

            # Store synapse
            self.buffer_merge_write(h5_syn_mat, syn_set)
            syn_ctr = syn_set.shape[0]

            loop_ctr = 0
            done = False

            while not done:

                old_unique_id = unique_id

                if loop_ctr % 1000000 == 0:
                    self.write_log(f"Worker synapses: {syn_ctr}/{n_total} (heap size: {len(synapse_heap)})",
                                   force_print=True)

                # Get the next set of synapses from this file from the iterator
                next_row_set = next(file_mat_iterator[h_id], None)

                if next_row_set is not None:
                    # More synapses in file, push next pair to heap, and pop top pair
                    syn_set, unique_id = next_row_set
                    (unique_id, h_id, syn_set) = heapq.heappushpop(synapse_heap, (unique_id, h_id, syn_set))
                elif len(synapse_heap) > 0:
                    (unique_id, h_id, syn_set) = heapq.heappop(synapse_heap)

                else:
                    done = True
                    continue

                assert unique_id >= old_unique_id, "uniqueID should be increasing in file"

                self.buffer_merge_write(h5_syn_mat, syn_set)
                syn_ctr += syn_set.shape[0]
                loop_ctr += 1

            self.write_log(f"Worker synapses: {syn_ctr}/{n_total} (heap size: {len(synapse_heap)})", force_print=True)
            self.write_log(f"Read {syn_ctr} out of total {n_total} synapses", force_print=True)

            self.buffer_merge_write(h5_syn_mat, flush=True)
            self.write_log("bigMergeHelper: done")

            # Close the hyper voxel files
            for f in file_list:
                try:
                    file_list[f].close()
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)
                    self.write_log(f"Problems closing files {f}, {file_list[f]}")

            self.buffer_out_file.close()

            return output_filename, neuron_range, syn_ctr
        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr, is_error=True)
            os.sys.exit(-1)

        ############################################################################

    def big_merge_lookup(self, merge_data_type="synapses", clean_voxel_files=True):

        # Since we want the code to work for both synapses and gap junction
        # we need to know location of synapse matrix, eg "network/synapses",
        # number of synapses, eg "nSynapses", the lookup table to quickly
        # find which synapse rows belongs to each pair of connected neurons
        # eg "network/synapseLookup"
        h5_syn_mat, h5_syn_n, h5_syn_lookup = self.data_loc[merge_data_type]

        self.merge_data_type = merge_data_type
        self.write_log(f"Doing big_merge_loopup for {merge_data_type}")

        synapse_heap = []

        num_neurons = len(self.hist_file["network/neurons/neuronID"])
        assert np.max(self.hist_file["network/neurons/neuronID"]) + 1 == num_neurons, \
            "bigMerge (lookup): There are neuron IDs missing"

        max_hyper_id = np.max(self.all_hyper_id_list) + 1
        file_list = [None] * max_hyper_id
        file_mat_iterator = [None] * max_hyper_id

        # fileMat = [None] * maxHyperID # points to the synapse matrix in each file
        # fileMatLookup = [None] * maxHyperID # points to the matrix lookup in file

        num_synapses = np.zeros((max_hyper_id,), dtype=np.int64)

        # Open all files for reading
        h_file_name_mask = os.path.join(self.network_path, "voxels", "network-putative-synapses-%s.hdf5")

        max_axon_voxel_ctr = 0
        max_dend_voxel_ctr = 0

        num_syn_hist = self.hist_file[h5_syn_n]
        num_syn_total = np.sum(num_syn_hist)

        # !!! Special code to increase h5py cache size
        propfaid = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
        settings = list(propfaid.get_cache())
        print(settings)
        # [0, 521, 1048576, 0.75]

        settings[2] *= 20
        propfaid.set_cache(*settings)
        settings = propfaid.get_cache()
        print(settings)
        # (0, 521, 5242880, 0.75)

        for h_id, nSyn, nOverflow in zip(self.hist_file["completed"], num_syn_hist,
                                         self.hist_file["voxelOverflowCounter"]):

            h_file_name = h_file_name_mask % str(h_id)

            if clean_voxel_files:
                # This makes sure we remove the old voxel files afterwards
                self.temp_file_list.append(h_file_name)

            if nSyn > 0:
                # Open file, and add info about first pairs synapses to the heap
                self.write_log(f"Opening voxel file: {h_file_name}")

                # Low level opening hdf5 file, to have greater cache size #ACC_RDWR
                fid = h5py.h5f.open(h_file_name.encode(), flags=h5py.h5f.ACC_RDONLY, fapl=propfaid)

                # !!! Temp print to check cache size
                settings = list(fid.get_access_plist().get_cache())
                print(settings)

                # fileList[hID] = h5py.File(hFileName,'r')
                file_list[h_id] = h5py.File(fid, drive=self.h5driver)
                file_mat_iterator[h_id] = self.synapse_set_iterator(h5mat_lookup=file_list[h_id][h5_syn_lookup],
                                                                    h5mat=file_list[h_id][h5_syn_mat],
                                                                    chunk_size=10000)

                num_synapses[h_id] = nSyn


                if self.max_channel_type:
                    if self.max_channel_type != file_list[h_id]["network/maxChannelTypeID"][()]:
                        print("Investigate")
                        print(f"{self.max_channel_type} != {file_list[h_id]['network/maxChannelTypeID'][()]}")
                        # import pdb
                        # pdb.set_trace()

                    # These should be the same for all hypervoxels
                    assert self.max_channel_type == file_list[h_id]["network/maxChannelTypeID"][()], \
                        (f"max_channel_type = {self.max_channel_type} "
                         f"(differ with what is in file {file_list[h_id]['network/maxChannelTypeID'][()]})")
                else:
                    self.max_channel_type = file_list[h_id]["network/maxChannelTypeID"][()]
                    print(f"Setting max_channel_type to {self.max_channel_type} from h_id={h_id}")
                

                # There should be at least the first row, otherwise nSyn = 0
                syn_set, unique_id = next(file_mat_iterator[h_id], None)

                # Create a heap containing the first subset of all files
                heapq.heappush(synapse_heap, (unique_id, h_id, syn_set))

                # This is so we can optimize the axon/dend voxelCtr and size
                if "maxAxonVoxelCtr" in file_list[h_id]["meta"]:
                    max_axon_voxel_ctr = max(max_axon_voxel_ctr, file_list[h_id]["meta/maxAxonVoxelCtr"][()])
                if "maxDendVoxelCtr" in file_list[h_id]["meta"]:
                    max_dend_voxel_ctr = max(max_dend_voxel_ctr, file_list[h_id]["meta/maxDendVoxelCtr"][()])

        assert np.sum(num_synapses) == num_syn_total, \
            f"Mismatch between work log file and data files: {num_syn_total} vs {np.sum(num_synapses)} synapses"

        if self.buffer_out_file is None:
            # Create output file
            (self.buffer_out_file, outFileName) = self.setup_merge_file(big_cache=True, delete_after=False)
        else:
            # We need to reset the write pointer (GJ and synapses should start from 0)
            self.next_file_write_pos = 0

        # Here we store the sorted connection matrix
        sorted_mat = self.buffer_out_file[h5_syn_mat]

        # Only save this meta data if doing the synapses call
        if max_axon_voxel_ctr > 0 and self.merge_data_type == "synapses":
            self.buffer_out_file["meta"].create_dataset("maxAxonVoxelCtr", data=max_axon_voxel_ctr)
            self.write_log(f"max_axon_voxel_ctr = {max_axon_voxel_ctr}")

        if max_dend_voxel_ctr > 0 and self.merge_data_type == "synapses":
            self.buffer_out_file["meta"].create_dataset("maxDendVoxelCtr", data=max_dend_voxel_ctr)
            self.write_log(f"max_dend_voxel_ctr = {max_dend_voxel_ctr}")

        # 2. Pop the smallest element, and add it to the final file
        # -- check same file if there are more with the same source and dest
        # -- buffer the writes

        syn_ctr = 0

        if len(synapse_heap) > 0:
            # Get the first file to read synapses from
            (unique_id, h_id, syn_set) = heapq.heappop(synapse_heap)
        else:
            # No synapses at all, return
            self.clean_up_merge_read_buffers()
            return sorted_mat, self.buffer_out_file

        # Store synapses
        syn_ctr = syn_set.shape[0]
        self.buffer_merge_write(h5_syn_mat, syn_set)

        loop_ctr = 0
        done = False

        while not done:

            if loop_ctr % 1000000 == 0:
                self.write_log(f"Synapses: {syn_ctr}/{num_syn_total} (heap size: {len(synapse_heap)})",
                               force_print=True)

            # Get the next set of synapses from this file from the iterator
            next_row_set = next(file_mat_iterator[h_id], None)

            if next_row_set is not None:
                # More synapses in file, push next pair to heap, and pop top pair
                syn_set, unique_id = next_row_set
                (unique_id, h_id, syn_set) = heapq.heappushpop(synapse_heap, (unique_id, h_id, syn_set))
            elif len(synapse_heap) > 0:
                (unique_id, h_id, syn_set) = heapq.heappop(synapse_heap)
            else:
                done = True
                continue

            # Write synapses to file
            self.buffer_merge_write(h5_syn_mat, syn_set)
            syn_ctr += syn_set.shape[0]
            loop_ctr += 1

        # Flush the buffers to file
        self.buffer_merge_write(h5_syn_mat, flush=True)
        self.write_log("big_merge_lookup: done")

        return sorted_mat, self.buffer_out_file

    ############################################################################

    def prune_synapses(self, synapse_file, output_filename, row_range,
                       merge_data_type,
                       close_input_file=True,
                       close_out_file=True):

        h5_syn_mat, h5_syn_n, h5_syn_loc = self.data_loc[merge_data_type]

        if row_range is None:
            row_start = 0
            row_end = synapse_file[h5_syn_mat].shape[0]
        else:
            row_start = row_range[0]
            row_end = row_range[-1]

        if row_start is None or row_end is None:
            self.write_log("Nothing to do, empty row range")
            return

        self.write_log("pruneSynapses called.")

        if type(synapse_file) != str and synapse_file[h5_syn_mat].shape[0] == 0:
            self.write_log(f"pruneSynapses: No {merge_data_type} skipping pruning")
            return

        self.write_log(f"pruneSynapses: synapseFile={synapse_file}, outputFileName={output_filename}"
                       f", rowRange={row_range} ({merge_data_type})")

        if type(synapse_file) == str:
            self.write_log(f"Opening synapse file: {synapse_file}")
            if self.role != "master":
                # SWMR = one writer, multiple readers
                synapse_file = h5py.File(synapse_file, 'r', swmr=True)
            else:
                synapse_file = h5py.File(synapse_file, 'r')

        num_syn = row_end - row_start

        # We need to split the rowRange into smaller pieces that fit in memory
        chunk_size = 1000000

        # To avoid division by zero
        num_blocks = max(1, int(np.ceil(float(row_end - row_start) / chunk_size)))

        self.write_log(f"About to calculate block ranges ({num_blocks} blocks)")

        block_ranges = self.find_ranges(synapses=synapse_file[h5_syn_mat],
                                        num_workers=num_blocks,
                                        start_pos=row_start,
                                        num_syn=num_syn)

        self.write_log(f"blockRanges={block_ranges}")

        self.setup_output_file(output_filename)  # Sets self.outFile

        for synRange in block_ranges:
            self.write_log(f"Pruning range: {synRange}")

            synapses = synapse_file[h5_syn_mat][synRange[0]:synRange[-1]]
            self.prune_synapses_helper(synapses=synapses, output_file=self.out_file,
                                       merge_data_type=merge_data_type)

        # Close synapse input file
        if close_input_file:
            synapse_file.close()

        if close_out_file:
            self.out_file.close()
            self.out_file = None

    ############################################################################

    # This code prunes synapses, but you need to send it small chunks
    # so synapses are all in memory, also it buffers the write until the end

    # synapses -- subset of synapse matrix that fits in memory
    # rowRange -- which rows to read
    # outputFile -- where to write synapses, assumed to already exist
    # outFilePos -- which position to start writing from

    def prune_synapses_helper(self, synapses, output_file, merge_data_type):

        h5_syn_mat, h5_syn_n, h5_syn_loc = self.data_loc[merge_data_type]

        keep_row_flag = np.zeros((synapses.shape[0],), dtype=bool)

        next_read_pos = 0
        read_end_of_range = synapses.shape[0]

        # Random seeds for reproducability
        neuron_seeds = self.get_neuron_random_seeds()
        previous_post_synaptic_neuron_id = None

        # Init some stats
        n_all_removed = 0
        n_some_removed = 0
        n_too_few_removed = 0
        n_dist_dep_pruning = 0
        n_too_many_removed = 0
        n_not_connected = 0

        old_pos = -1

        while next_read_pos < read_end_of_range:

            assert old_pos != next_read_pos, "prune_synapses_helper: Stuck in a loop."
            old_pos = next_read_pos

            # How many lines contain synapses between this pair of neurons
            read_end_idx = next_read_pos + 1

            if merge_data_type == "gapJunctions":
                while (read_end_idx < read_end_of_range and
                       (synapses[next_read_pos, 0:2] == synapses[read_end_idx, 0:2]).all()):
                    read_end_idx += 1
            else:
                while (read_end_idx < read_end_of_range and
                       (synapses[next_read_pos, 0:2] == synapses[read_end_idx, 0:2]).all()  # Same neuron pair
                        and synapses[next_read_pos, 6] == synapses[read_end_idx, 6]):       # Same synapse type
                    read_end_idx += 1

            # Temp check
            assert ((synapses[next_read_pos:read_end_idx, 0] == synapses[next_read_pos, 0]).all()
                    and (synapses[next_read_pos:read_end_idx, 1] == synapses[next_read_pos, 1]).all()), \
                "prune_synapses_helper: Internal error, more than one neuron pair"

            # Stats
            n_pair_synapses = read_end_idx - next_read_pos

            src_id = synapses[next_read_pos, 0]
            dest_id = synapses[next_read_pos, 1]

            if dest_id != previous_post_synaptic_neuron_id:
                # New post synaptic cell, reseed random generator
                self.write_log(f"Random seed set for neuron {dest_id}: {neuron_seeds[dest_id]}")  # Temp logging
                post_rng = np.random.default_rng(neuron_seeds[dest_id])
                previous_post_synaptic_neuron_id = dest_id

            if merge_data_type == "gapJunctions":
                # All are gap junctions
                synapse_type = 3
            else:
                synapse_type = synapses[next_read_pos, 6]

                assert (synapses[next_read_pos:read_end_idx, 6] == synapse_type).all(), \
                    (f"More than one synapse type connecting "
                     f"{self.hist_file['network/neurons/name'][synapses[next_read_pos,0]]}  "
                     f"(ID {synapses[next_read_pos,0]}) and "
                     f"{self.hist_file['network/neurons/name'][synapses[next_read_pos,1]]} "
                     f"(ID {synapses[next_read_pos,1]})")

            con_id = (self.type_id_list[src_id], self.type_id_list[dest_id], synapse_type)

            if con_id in self.connectivity_distributions:

                # We have the option to separate between connections within a
                # population unit or not. If conInfo[1] != None then first
                # tuple is connection info within a population unit, and second item
                # is connection info between different population units
                con_info = self.connectivity_distributions[con_id]

                #
                if con_info[1] is None or self.population_unit_id[src_id] == self.population_unit_id[dest_id]:
                    # All or within population unit pruning parameters
                    c_info = con_info[0]
                else:
                    # Between population unit pruning parameters
                    c_info = con_info[1]

                # These will always exist thanks to completePruningInfo function

                dist_p = c_info["distPruning"]  # Dist dep pruning
                f1 = c_info["f1"]
                soft_max = c_info["softMax"]
                mu2 = c_info["mu2"]
                a3 = c_info["a3"]

            else:
                # Not listed in connectivityDistribution, skip neuron pair
                next_read_pos = read_end_idx
                # No need to update keepRowFlag since default set to 0

                # Stats
                n_not_connected += n_pair_synapses
                continue

            # 3. This is the last step of pruning, but we move it to the top
            # since there is no point doing the other steps if we going to
            # throw them away anyway
            if a3 is not None and post_rng.random() > a3:
                # Prune all synapses between pair, do not add to synapse file
                next_read_pos = read_end_idx
                # No need to update keepRowFlag since default set to 0

                # Stats
                n_all_removed += n_pair_synapses
                continue

            if dist_p is not None:
                assert synapse_type != 3, \
                    "Distance dependent pruning currently only supported for synapses, not gap junctions"
                # Distance dependent pruning, used for e.g. FS->MS connections

                # distP contains d (variable for distance to soma)
                d = synapses[next_read_pos:read_end_idx, 8] * 1e-6  # dendrite distance d, used in eval below
                p = eval(dist_p)
                frac_flag = post_rng.random(n_pair_synapses) < f1
                dist_flag = post_rng.random(n_pair_synapses) < p

                keep_row_flag[next_read_pos:read_end_idx] = np.logical_and(frac_flag, dist_flag)

                n_frac = sum(frac_flag)
                n_some_removed += n_pair_synapses - n_frac
                n_dist_dep_pruning += n_frac - sum(keep_row_flag[next_read_pos:read_end_idx])

            else:
                keep_row_flag[next_read_pos:read_end_idx] = post_rng.random(n_pair_synapses) < f1
                n_some_removed += n_pair_synapses - sum(keep_row_flag[next_read_pos:read_end_idx])

            # Check if too many synapses, trim it down a bit
            n_keep = np.sum(keep_row_flag[next_read_pos:read_end_idx])

            if soft_max is not None and n_keep > soft_max:
                # pKeep = float(softMax)/nKeep # OLD implementation
                soft_max = float(soft_max)
                # pKeep = 2*softMax*np.divide(1-np.exp(-nKeep/softMax),1+np.exp(-nKeep/softMax))/nKeep
                p_keep = np.divide(2 * soft_max, (1 + np.exp(-(n_keep - soft_max) / 5)) * n_keep)

                keep_row_flag[next_read_pos:read_end_idx] = \
                    np.logical_and(p_keep > post_rng.random(n_pair_synapses),
                                   keep_row_flag[next_read_pos:read_end_idx])

                # Stats
                n_too_many_removed += n_keep - sum(keep_row_flag[next_read_pos:read_end_idx])

                # Update count
                n_keep = np.sum(keep_row_flag[next_read_pos:read_end_idx])

            # If too few synapses, remove all synapses
            if mu2 is not None:
                p_mu = 1.0 / (1.0 + np.exp(-8.0 / mu2 * (n_keep - mu2)))

                if p_mu < post_rng.random():
                    # Too few synapses, remove all -- need to update keepRowFlag
                    keep_row_flag[next_read_pos:read_end_idx] = 0
                    next_read_pos = read_end_idx

                    # Stats
                    n_too_few_removed += n_keep
                    continue

            next_read_pos = read_end_idx

            # Time to write synapses to file
        n_keep_tot = sum(keep_row_flag)
        write_start_pos = int(output_file["network/" + h5_syn_n][0])
        write_end_pos = write_start_pos + n_keep_tot

        if n_keep_tot > 0:
            output_file[h5_syn_mat].resize((write_end_pos, output_file[h5_syn_mat].shape[1]))
            output_file[h5_syn_mat][write_start_pos:write_end_pos] = \
                synapses[keep_row_flag, :]

            # Update counters
            output_file["network/" + h5_syn_n][0] = write_end_pos

        else:
            self.write_log("No synapses kept, resizing")
            output_file[h5_syn_mat].resize((write_end_pos, output_file[h5_syn_mat].shape[1]))

        self.write_log(f"Number of synapses removed where synapse connection not allowed: {n_not_connected}"
                       f"\nNumber of synapses removed due to distance dependent pruning: {n_dist_dep_pruning}"
                       f"\nNumber of synapses removed randomly: {n_some_removed}"
                       f"\nNumber of synapses removed due to too many synapses between connected pair: "
                       f"{n_too_many_removed}"
                       f"\nNumber of synapses removed due to too few synapses between connected pairs: "
                       f"{n_too_few_removed}"
                       f"\nNumber of synapses removed where all synapses between pairs are removed: "
                       f"{n_all_removed}", force_print=True)

    ############################################################################

    @staticmethod
    def file_row_lookup_iterator(h5mat, chunk_size=10000):

        mat_size = h5mat.shape[0]

        # If matrix is smaller than chunk, buffer can be smaller than requested
        if mat_size < chunk_size:
            chunk_size = h5mat.shape[0]

        mat_buf = np.zeros((chunk_size, h5mat.shape[1]), dtype=h5mat.dtype)
        end_idx = 0

        while end_idx < mat_size:
            start_idx = end_idx
            end_idx = start_idx + chunk_size

            if end_idx < mat_size or mat_size == chunk_size:
                # Copy to existing buffer, to avoid memory allocation
                mat_buf[:, :] = h5mat[start_idx:end_idx, :]
            else:
                # Create a new buffer
                mat_buf = h5mat[start_idx:end_idx, :]

            for row in mat_buf:
                yield row

    ############################################################################

    # min_dest_id and max_dest_id are inclusive, only synapses with dest_id in that range are iterated over

    def file_row_lookup_iterator_subset(self, h5mat_lookup, min_dest_id, max_dest_id, chunk_size=10000):

        num_neurons = self.hist_file["network/neurons/neuronID"].shape[0]
        min_unique_id = min_dest_id * num_neurons * self.max_channel_type
        max_unique_id = max_dest_id * num_neurons * self.max_channel_type

        self.write_log(f"min_unique_id: {min_unique_id}, max_unique_id: {max_unique_id}")
        mat_size = h5mat_lookup.shape[0]

        # If matrix is smaller than chunk, buffer can be smaller than requested
        if mat_size < chunk_size:
            chunk_size = mat_size

        mat_buf = np.zeros((chunk_size, h5mat_lookup.shape[1]), dtype=h5mat_lookup.dtype)
        end_idx = 0

        while end_idx < mat_size:
            start_idx = end_idx
            end_idx = start_idx + chunk_size

            if end_idx < mat_size or mat_size == chunk_size:
                # Copy to existing buffer, to avoid memory allocation
                mat_buf[:, :] = h5mat_lookup[start_idx:end_idx, :]
            else:
                # Create a new buffer
                mat_buf = h5mat_lookup[start_idx:end_idx, :]

            for row in mat_buf:
                if min_unique_id <= row[0] < max_unique_id:
                    # Only return synapses that terminate on the neurons we are
                    # interested in here

                    yield row

    ############################################################################

    # Returns (subset of rows, uniqueID)

    def synapse_set_iterator(self, h5mat_lookup, h5mat, chunk_size=10000, lookup_iterator=None):

        # Allow the user to set an alternative lookupIterator if we only
        # want to iterate over a subset of the synapses
        if not lookup_iterator:
            lookup_iterator = self.file_row_lookup_iterator(h5mat_lookup, chunk_size=chunk_size)

        mat_size = h5mat.shape[0]
        if mat_size < chunk_size:
            chunk_size = mat_size

        old_synapses = None
        read_buffer = h5mat[:chunk_size, :].copy()
        buffer_start = 0  # What file pos does start of buffer correspond to
        buffer_end = chunk_size  # What file pos does end of buffer correspond to (+1)

        next_row_set = next(lookup_iterator, None)

        while next_row_set is not None:

            # start_idx and end_idx are the rows in the matrix we want to read between
            [unique_id, start_idx, end_idx] = next_row_set

            if start_idx >= buffer_end:
                # We need to jump forward...
                buffer_start = start_idx

                if start_idx + chunk_size <= h5mat.shape[0]:
                    buffer_end = start_idx + chunk_size
                    read_buffer[:, :] = h5mat[buffer_start:buffer_end, :]
                else:
                    buffer_end = h5mat.shape[0]
                    read_buffer = h5mat[buffer_start:buffer_end, :].copy()

                old_synapses = None

            if end_idx > buffer_end:
                assert old_synapses is None, "get_next_synapse_set: chunk_size too small"

                # Part of the synapse range requested is outside buffer
                old_synapses = read_buffer[(start_idx - buffer_start):, :].copy()

                buffer_start = buffer_end
                buffer_end = buffer_start + chunk_size

                if buffer_end > mat_size:
                    buffer_end = mat_size
                    read_buffer = h5mat[buffer_start:buffer_end, :].copy()
                else:
                    # Reuse old buffer storage
                    read_buffer[:, :] = h5mat[buffer_start:buffer_end, :]

                # Need to concatenate with old synapses
                syn_mat = np.concatenate([old_synapses,
                                         read_buffer[:(end_idx - buffer_start), :]],
                                         axis=0)

                assert end_idx == buffer_end \
                    or (read_buffer[start_idx - buffer_start, :2] != read_buffer[end_idx - buffer_start, :2]).any(), \
                    "We missed one synpase! (2)"

                assert (syn_mat[:, 0] == syn_mat[0, 0]).all() and (syn_mat[:, 1] == syn_mat[0, 1]).all(), \
                    f"Synapse matrix (2) contains more than one pair:\n{syn_mat}"

                assert syn_mat.shape[0] == end_idx - start_idx, "Synapse matrix has wrong size"

                yield (syn_mat,
                       unique_id)

                old_synapses = None
            else:
                syn_mat = read_buffer[(start_idx - buffer_start):(end_idx - buffer_start), :]

                assert end_idx == buffer_end \
                       or (read_buffer[start_idx - buffer_start, :2]
                           != read_buffer[end_idx - buffer_start, :2]).any(), \
                       "We missed one synpase! (1)"

                assert (syn_mat[:, 0] == syn_mat[0, 0]).all() and (syn_mat[:, 1] == syn_mat[0, 1]).all(), \
                       f"Synapse matrix (1) contains more than one pair:\n{syn_mat}"

                assert syn_mat.shape[0] == end_idx - start_idx, \
                    "Synapse matrix has wrong size"

                yield (syn_mat,
                       unique_id)

            next_row_set = next(lookup_iterator, None)

    def get_neuron_random_seeds(self):

        # Cache using ... self.old_seed = self.random_seed
        num_neurons = self.hist_file["network/neurons/neuronID"].shape[0]
        assert num_neurons - 1 == self.hist_file["network/neurons/neuronID"][-1], \
            "neuronID should start from 0 and the end should be n-1"

        # Need different seeds for each post synaptic neuron
        ss = np.random.SeedSequence(self.random_seed)
        neuron_seeds = ss.generate_state(num_neurons)
        return neuron_seeds


##############################################################################

if __name__ == "__main__":
    print("Please do not call this file directly, use snudda.py")
    os.sys.exit(-1)
