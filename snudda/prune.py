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

from snudda.Neuron_morphology import NeuronMorphology


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

    def __init__(self, work_history_file,
                 logfile=None, logfile_name=None,
                 d_view=None, lb_view=None, role="master", verbose=True,
                 scratch_path=None,
                 pre_merge_only=False,
                 h5libver="latest",
                 random_seed=None,
                 clean_voxel_files=True):

        # Help with parallel debugging, when workers cant print to screen:
        # self.writeToRandomFile("WH = " + str(workHistoryFile) \
        #                       + "\nlog = " + str(logFileName) \
        #                       + "\nrole = " + role)

        start_time = timeit.default_timer()
        self.work_history_file = work_history_file
        self.base_path = os.path.dirname(self.work_history_file) + "/"
        self.base_path = self.base_path.replace("/log/", "/")

        self.logfile = logfile
        self.logfile_name = logfile_name
        self.verbose = verbose
        self.h5libver = h5libver

        if logfile is None and logfile_name is not None:
            self.logfile = open(logfile_name, 'w')
            self.write_log("Log file created.")

        np.random.seed(random_seed)
        self.write_log("Setting random seed: " + str(random_seed))

        self.h5driver = "sec2"  # "core" # "stdio", "sec2"

        self.write_log("Using hdf5 driver " + str(self.h5driver) \
                       + ", " + str(self.h5libver) + " version")

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

        self.hist_file = None
        self.out_file = None
        self.temp_file_list = []  # List of all temp files created

        self.next_merge_file_id = 0

        self.voxel_overflow_counter = 0
        self.overflow_files = []

        self.open_work_history_file(work_history_file=work_history_file)

        self.set_scratch_path(scratch_path)
        self.load_pruning_information()

        # (locationOfMatrix,locationOfN,locationOfCoords)
        self.data_loc = {"synapses": ("network/synapses",
                                      "nSynapses",
                                      "network/synapseLookup"),  # range(2,5)),
                         "gapJunctions": ("network/gapJunctions",
                                          "nGapJunctions",
                                          "network/gapJunctionLookup")}  # range(6,9))}

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

                    # print("We do not currently get all the synapses, where are the missing ones going? About 1.4% of synapses are missed")
                    # import pdb
                    # pdb.set_trace()

                    assert not (synapse_file["network/synapses"][-1, :] == 0).all(), \
                        "We are missing some synapses in the merge!"

                else:
                    self.write_log("Running merge in serial")

                    (synapses, synapse_file) = \
                        self.big_merge_lookup(merge_data_type="synapses",
                                              clean_voxel_files=True)

                    (gapJunctions, gapJunctionFile) = \
                        self.big_merge_lookup(merge_data_type="gapJunctions",
                                              clean_voxel_files=True)

            # When running on a cluster, we might want to do the serial parts
            # in a separte run, hence this option that allows us to prepare
            # the data for parallel execution.
            if pre_merge_only:
                self.write_log("Pre-merge of synapses done. preMergeOnly = " \
                               + str(pre_merge_only) + ", exiting.")

                if self.hist_file is not None:
                    self.hist_file.close()
                    self.hist_file = None

                end_time = timeit.default_timer()

                self.write_log("Total duration: " + str(end_time - start_time))
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
                self.write_log("No output file created, no synapses exist?")
                self.write_log("Creating symbolic link to MERGE file instead")

                f_dest = (os.path.dirname(self.work_history_file) + "/").replace("/log/", "/") \
                         + "/network-pruned-synapses.hdf5"
                f_src = os.path.basename(f_name)

                print(str(f_dest) + "->" + str(f_src))
                os.symlink(f_src, f_dest)

                return

            n_syn_before = np.sum(self.hist_file["nSynapses"][()])
            n_syn_after = self.out_file["network/nSynapses"][0]

            n_overflow = np.sum(self.hist_file["voxelOverflowCounter"][()])
            n_gj_before = np.sum(self.hist_file["nGapJunctions"][()])
            n_gj_after = self.out_file["network/nGapJunctions"][0]

            self.write_log("Voxel overflows: " + str(n_overflow) \
                           + " (should be zero)")

            if n_syn_before > 0:
                self.write_log("Synapses before pruning: " + str(n_syn_before))
                self.write_log("Synapses after pruning: " + str(n_syn_after)
                               + " (" + str(round(100.0 * n_syn_after / n_syn_before, 2))
                               + " % kept)")
            else:
                self.write_log("No synapses to prune")

            if n_gj_before > 0:
                self.write_log("Gap junctions before pruning " + str(n_gj_before))
                self.write_log("Gap junctions after pruning " + str(n_gj_after)
                               + " (" + str(round(100.0 * n_gj_after / n_gj_before, 2))
                               + " % kept)")
            else:
                self.write_log("No gap junctions to prune.")

            self.clean_up_merge_files()

            if self.hist_file is not None:
                self.hist_file.close()
                self.hist_file = None

            end_time = timeit.default_timer()

            self.write_log("Total duration: " + str(end_time - start_time))

            if self.voxel_overflow_counter > 0:
                print("Voxel overflows: " + str(self.voxel_overflow_counter) + "\nIn files: ")
                for f in self.overflow_files:
                    print("Overflow in " + str(f))

    ############################################################################

    def set_scratch_path(self, scratch_path=None):

        assert self.work_history_file is not None \
               and self.work_history_file is not "last", \
            "Need to call openWorkHistoryFile before setScratchPath"

        if scratch_path is None:
            self.scratchPath = self.base_path + "/temp/"

            if not os.path.exists(self.scratchPath):
                os.makedirs(self.scratchPath)

            self.write_log("Using default scratch path: " + self.scratchPath)
        else:
            self.scratchPath = scratch_path
            self.write_log("User selected scratch path: " + self.scratchPath)

    ############################################################################

    def write_to_random_file(self, text):

        import uuid
        tmp = open("save/tmp-log-file-" + str(uuid.uuid4()), 'w')
        tmp.write(text)
        tmp.close()
        print(text)

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

        self.clean_up_merge_files()

    ############################################################################

    def open_work_history_file(self, work_history_file=None):

        if work_history_file is None:
            work_history_file = self.work_history_file

        if work_history_file == "last":
            work_history_file = self.find_latest_file()

        self.write_log("Opening work history file: " + work_history_file)

        self.hist_file = h5py.File(work_history_file, 'r')
        self.work_history_file = work_history_file

        self.slurm_id = self.hist_file["meta/SlurmID"][()]
        self.hyper_voxel_i_ds = self.hist_file["meta/hyperVoxelIDs"][()]
        self.all_hyper_i_ds = self.hist_file["allHyperIDs"][()]
        self.voxel_size = self.hist_file["meta/voxelSize"][()]
        self.hyper_voxel_size = self.hist_file["meta/hyperVoxelSize"][()]  # num bins
        self.simulation_origo = self.hist_file["meta/simulationOrigo"][()]
        self.hyper_voxel_width = self.voxel_size * self.hyper_voxel_size

        # Network_simulate.py uses axonStumpIDFlag = True
        # Neurodamus uses axonStumpIDFlag = False
        self.axon_stump_id_flag = self.hist_file["meta/axonStumpIDFlag"]

        # We need a lookup table for offsets of hypervoxel origos
        self.hyper_voxel_offset = np.zeros((self.hyper_voxel_i_ds.size, 3), dtype=int)
        for ix in range(0, self.hyper_voxel_i_ds.shape[0]):
            for iy in range(0, self.hyper_voxel_i_ds.shape[1]):
                for iz in range(0, self.hyper_voxel_i_ds.shape[2]):
                    self.hyper_voxel_offset[self.hyper_voxel_i_ds[ix, iy, iz], :] \
                        = np.array([ix, iy, iz]) * self.hyper_voxel_size

        # OBS the synapse and gap junction numbers are listed not in order of
        # hyperID but in order they were completed, to find out hyperID for
        # a particular one check self.histFile["completed"]
        self.num_synapses_total = np.sum(self.hist_file["nSynapses"][()])
        self.num_gap_junctions_total = np.sum(self.hist_file["nGapJunctions"][()])

        self.config_file = self.hist_file["meta/configFile"][()]
        self.position_file = self.hist_file["meta/positionFile"][()]

        self.detect_config = json.loads(self.hist_file["meta/config"][()])
        with open(self.hist_file["meta/configFile"][()], "r") as f:
            self.config = json.load(f)

            ############################################################################

    def check_hyper_voxel_integrity(self, hypervoxel_file, hypervoxel_file_name, verbose=False):

        if verbose:
            self.write_log("Checking that " + hypervoxel_file_name + " matches circuit settings")

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
        xyz = np.where(self.hyper_voxel_i_ds == hypervoxel_file["meta/hyperVoxelID"][()])
        xyz = np.array([x[0] for x in xyz])

        # Just do a sanity check that the hypervoxel origo matches stored value
        hypervoxel_origo = self.simulation_origo + self.hyper_voxel_width * xyz
        assert (hypervoxel_origo == hypervoxel_file["meta/hyperVoxelOrigo"][()]).all(), \
            "Hyper voxel origo mismatch in file " + hypervoxel_file_name

        ofc = hypervoxel_file["meta/voxelOverflowCounter"][()]

        if ofc > 0:
            self.voxel_overflow_counter += ofc
            self.overflow_files.append(hypervoxel_file_name)
            self.write_log("Overflow of " + str(ofc) + " in " + hypervoxel_file_name)

    ############################################################################

    def open_hyper_voxel(self, hyper_voxel_id, verbose=False, verify=True):

        if verbose:
            self.write_log("Reading hypervoxel " + str(hyper_voxel_id))

        h_file_name = self.base_path + "/voxels/network-putative-synapse-" \
                      + str(hyper_voxel_id) + ".hdf5"

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

        self.write_log("Creating the " + str(num_neurons) + "x" + str(num_neurons) \
                       + " sparse connection matrix")

        # After creation, convert to csr or csc matrix
        con_mat = scipy.sparse.lil_matrix((num_neurons, num_neurons), dtype=np.int16)

        self.write_log("Parsing " + str(len(self.all_hyper_i_ds)) + " hyper voxel files")

        for ctr, h_id in enumerate(self.all_hyper_i_ds):
            if ctr % 1000 == 0 and ctr > 0:
                self.write_log(str(ctr) + " files done")

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

        cache_file = self.base_path + "/connection-matrix-cache.pickle"

        return cache_file

    ############################################################################

    def load_connection_matrix_cache(self):

        cache_file = self.get_con_mat_cache_filename()

        if os.path.isfile(cache_file):
            self.write_log("Loading connection matrix cache from " + cache_file)

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
                        assert test, "Connection matrix cache - mismatch for " \
                                     + checkNames[1]
                    else:
                        assert test.all(), \
                            "Connection matrix cache - mismatch for " + checkNames[1]

                assert self.num_synapses_total == data["nSynapsesTotal"] \
                       and self.num_gap_junctions_total == data["nGapJunctionsTotal"], \
                    " Synapse or gap junction count mismatch -- corrupt file?"

        else:
            con_mat = None

        return con_mat

    ############################################################################

    def save_connection_matrix_cache(self, con_mat):

        cache_file = self.get_con_mat_cache_filename()
        self.write_log("Saving connection matrix cache to " + cache_file)

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

    def check_network_config_integrity(self):

        detect_config = json.loads(self.hist_file["meta/config"][()])
        with open(self.hist_file["meta/configFile"][()], "r") as f:
            prune_config = json.load(f)

        all_present = True

        for con in prune_config["Connectivity"]:
            if con not in detect_config["Connectivity"]:
                self.write_log("!!! Connection " + con + " has been added to " \
                               + self.hist_file["meta/configFile"][()] \
                               + " after detection, please rerun snudda detect")
                all_present = False

        assert all_present, "Please rerun snudda detect."

    ############################################################################

    # Parse the connection information in the config file, stored in the
    # the work history file

    def load_pruning_information(self):

        # self.config = json.loads(self.histFile["meta/config"][()])

        self.check_network_config_integrity()
        with open(self.hist_file["meta/configFile"][()], "r") as f:
            self.config = json.load(f)

        self.populationUnitID = self.hist_file["network/neurons/populationUnitID"][()]

        # Normally we use type names as lookups, but since we will do this
        # many millions of times, we create an temporary typeID number
        self.make_type_numbering()

        orig_connectivity_distributions = \
            json.loads(self.hist_file["meta/connectivityDistributions"][()])

        # origConnectivityDistributionsGJ = \
        #  json.loads(self.histFile["meta/connectivityDistributionsGJ"][()])

        self.connectivity_distributions = dict([])

        # For the pruning we merge the two into one
        for keys in orig_connectivity_distributions:
            (pre_type, post_type) = keys.split("$$")

            # Need to handle if preType or postType dont exist, then skip this
            if pre_type not in self.type_id_lookup or post_type not in self.type_id_lookup:
                print("Skipping " + pre_type + " to " + post_type + " connection")
                continue

            pre_type_id = self.type_id_lookup[pre_type]
            post_type_id = self.type_id_lookup[post_type]

            for con_type in orig_connectivity_distributions[keys]:
                con_data = orig_connectivity_distributions[keys][con_type]

                pruning = self.complete_pruning_info(con_data["pruning"])

                if "pruningOther" in con_data:
                    pruning_other = self.complete_pruning_info(con_data["pruningOther"])
                else:
                    pruning_other = None

                synapse_type_id = con_data["channelModelID"]

                self.connectivity_distributions[pre_type_id, post_type_id, synapse_type_id] \
                    = (pruning, pruning_other)

    ############################################################################

    # This makes sure all the variables exist, that way pruneSynapseHelper
    # does not have to check, but can assume that they will exist

    def complete_pruning_info(self, prune_info):

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

    def loadPruningInformationOLD(self):

        assert False, "loadPruningInformationOLD -- unused?? delete..."

        self.config = json.loads(self.hist_file["meta/config"][()])

        self.population_unit_id = self.hist_file["network/neurons/populationUnitID"][()]

        # Normally we use type names as lookups, but since we will do this
        # many millions of times, we create an temporary typeID number
        self.make_type_numbering()

        self.connectivity_distributions = dict([])

        for name, definition in self.config.items():

            if name in ["Volume", "PopulationUnits"]:
                # We are just loading the pruning information
                continue

            if "targets" in definition:
                conDist = definition["targets"]
            else:
                conDist = []

            for row in conDist:
                target = row[0]

                try:
                    synapseType = self.synapse_type_reverse_lookup[row[1][0]]
                    # Skipping loading conductance mean and std,
                    # and channelParamDictionary
                    distrib = row[2]
                except:
                    self.write_log("Something wrong with your network json file. row: " \
                                   + str(row))
                    import pdb
                    pdb.set_trace()

                if len(row) > 3:
                    # If distrib2 is specified, then distrib is within population Unit
                    # distribution and distrib2 is between population unit distribution
                    distrib2 = row[3]
                else:
                    distrib2 = None

                typeName = name.split("_")[0]

                try:
                    if target in self.type_id_lookup:
                        self.connectivity_distributions[self.type_id_lookup[typeName],
                                                        self.type_id_lookup[target],
                                                        synapseType] \
                            = (distrib, distrib2)
                    else:
                        self.write_log("WARNING: Target neuron " + str(target) + " not included in simulation")
                except:

                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)

                    import pdb
                    pdb.set_trace()

    ############################################################################

    def createConnectionMatrixOLD(self):

        assert False, "createConnectionMatrixOLD unused --- delete??"

        if self.hist_file is None:
            self.open_work_history_file()

        nNeurons = self.hist_file["network/neurons/neuronID"].shape[0]

        self.write_log("Creating the " + str(nNeurons) + "x" + str(nNeurons) \
                       + " sparse connection matrix")

        # After creation, convert to csr or csc matrix
        conMat = scipy.sparse.lil_matrix((nNeurons, nNeurons), dtype=np.int16)

        assert self.synapses.shape[0] > 0, "No synapses in synapse matrix"

        prevSrcID = self.synapses[0, 0]
        prevDestID = self.synapses[0, 1]
        prevCtr = 1

        # The reason we dont write the synapse directly is that each access
        # to the sparse matrix takes time, so better group synapses
        # between same pairs together to one insert.
        # This requires self.synapess to be sorted.

        for row in self.synapses[1:, :]:
            SrcID = row[0]
            DestID = row[1]

            if SrcID == prevSrcID and DestID == prevDestID:
                prevCtr += 1
            else:
                conMat[prevSrcID, prevDestID] += prevCtr
                prevSrcID = SrcID
                prevDestID = DestID
                prevCtr = 1

        # Dont forget last one
        conMat[prevSrcID, prevDestID] += prevCtr

        self.write_log("Converting to CSR sparse matrix format")
        return conMat.tocsr()

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
            self.write_log("Output file already set: " + str(self.out_file.filename))
            return

        if output_file is None:
            output_file = self.base_path + "/network-pruned-synapses.hdf5"

        self.write_log("Writing to " + output_file)

        out_file = h5py.File(output_file, "w",
                             libver=self.h5libver,
                             driver=self.h5driver)

        out_file.create_dataset("config", data=json.dumps(self.config))

        # Copy over meta data
        self.hist_file.copy("meta", out_file)

        # Morphologies not included at this stage -- maybe add them?
        # Same morphologies
        cfg = json.loads(self.hist_file["meta/config"][()])

        morph_group = out_file.create_group("morphologies")

        for name, definition in cfg["Neurons"].items():
            try:
                morph_file = definition["morphology"]
            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr)
                import pdb
                pdb.set_trace()

            with open(morph_file, "r") as f:
                swc_data = f.read()

            self.write_log("Saving morphology in HDF5 file: " + morph_file)
            swc_group = morph_group.create_group(name)
            swc_group.create_dataset("swc", data=swc_data)
            swc_group.create_dataset("location", data=morph_file)

        network_group = out_file.create_group("network")

        # Copy over neuron data
        # self.histFile.copy("neurons",outFile)
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
        # self.outFileSynapseCtr = np.int64(0)
        # self.outFileGapJunctionCtr = np.int64(0)

    ############################################################################

    def find_latest_file(self):

        # files = glob('save/network_connect_voxel_log-*-worklog.hdf5')
        files = glob('save/network-connect-synapse-voxel-file-*-worklog.hdf5')

        mod_time = [os.path.getmtime(f) for f in files]
        idx = np.argsort(mod_time)

        self.write_log("Using the newest file: " + files[idx[-1]])

        return files[idx[-1]]

    ############################################################################

    def merge_file_exists(self):

        # check if merge file exists
        merge_file_name = self.base_path + "/network-putative-synapses-MERGED.hdf5"

        merge_file_ok = False

        self.write_log("Checking for merge file " + str(merge_file_name))

        try:
            if os.path.isfile(merge_file_name):

                merge_file_ok = True

                f = h5py.File(merge_file_name)

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
                elif (f["network/gapJunctions"][-1, :] == 0).all():
                    # Last row not set, not complete
                    merge_file_ok = False

                # !!! Add tests for gap junctions also
                # import pdb
                # pdb.set_trace()
            else:
                f = None
        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)

            self.write_log("Something went wrong with reading old merge file, ignoring it")
            # Something went wrong
            merge_file_ok = False

        if merge_file_ok:
            self.write_log("Found old merge file " + str(merge_file_name))
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
            outfile_name = self.base_path + "/network-putative-synapses-MERGED.hdf5"

        #  Make a list of all temporary files so we can remove them
        if delete_after:
            self.temp_file_list.append(outfile_name)

        self.write_log("Setting up out file " + str(outfile_name))
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
                try:
                    morph_file = definition["morphology"]
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)
                    import pdb
                    pdb.set_trace()

                with open(morph_file, "r") as f:
                    swc_data = f.read()

                self.write_log("Saving morphology in HDF5 file: " + morph_file)
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

        # import pdb
        # pdb.set_trace()
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
            self.write_log("prune_synapses_parallel: No " + merge_data_type
                           + " skipping pruning")
            return
        else:
            self.write_log("prune_synapses_parallel, before pruning : "
                           + str(synapse_file[h5_syn_mat].shape[0]) + " "
                           + merge_data_type)

        start_time = timeit.default_timer()

        if self.d_view is not None and self.role == "master":
            self.setup_parallel(dView=self.d_view)

            try:
                # Make sure all synapse writes are on the disk
                synapse_file.flush()
                if not synapse_file.swmr_mode:
                    synapse_file.swmr_mode = True  # Allow multiple readers from the file
            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr)
                import pdb
                pdb.set_trace()

                # 1. Pick names for the workers
        temp_output_file_name = [self.scratchPath + "worker-temp-" + merge_data_type + "-file-" + str(x)
                                 for x in range(0, len(self.d_view))]
        self.d_view.scatter("output_filename", temp_output_file_name, block=True)

        # Add the files to a delete list, so we remove them after
        for f in temp_output_file_name:
            self.temp_file_list.append(f)

        # 2. Define what ranges each worker should do
        synapse_ranges = self.find_ranges(synapse_file[h5_syn_mat],
                                          len(self.d_view))

        if synapse_ranges is None or synapse_ranges[-1][-1] is None:
            self.write_log("There are few synapses, we will run it in serial instead")
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

        fn = synapse_file.filename

        # Close the synapse file on the master node
        if close_input_file:
            synapse_file.close()
            synapse_file = None

        # 3. Let the workers prune
        self.write_log("Sending pruning job to workers")

        cmd_str = "nw.prune_synapses(synapse_file=synapse_filename,output_filename=output_filename[0]," \
                  + "row_range=synapse_range[0],merge_data_type=merge_data_type)"

        self.d_view.execute(cmd_str, block=True)

        end_time = timeit.default_timer()
        self.write_log("Parallel pruning duration: " + str(end_time - start_time))

        # 4. Merge the resulting files -- this is easier, since synapses
        #    are already sorted in the file
        self.write_log("Merging parallel pruning results")

        if setup_out_file:
            assert self.out_file is None, \
                "pruneSynapsesParallel: Output file already setup"
            self.setup_output_file(output_file=output_file)
        else:
            # If there were no synapses, but there are gap junctions
            # then this assert could theoretically happen
            assert self.out_file is not None, \
                "pruneSynapsesParallel: Out file not set up"

        tmp_files = [h5py.File(f, 'r') for f in temp_output_file_name]

        num_syn = np.sum(f[h5_syn_mat].shape[0] for f in tmp_files)
        mat_width_all = [f[h5_syn_mat].shape[1] for f in tmp_files]

        try:

            assert (np.array(mat_width_all) == mat_width_all[0]).all(), \
                "Internal error, width does not match"
            mat_width = mat_width_all[0]

            self.out_file[h5_syn_mat].resize((num_syn, mat_width))

        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)
            import pdb
            pdb.set_trace()

        next_syn = 0

        for f in tmp_files:

            n = f[h5_syn_mat].shape[0]

            if n > 0:
                self.out_file[h5_syn_mat][next_syn:(next_syn + n), :] = \
                    f[h5_syn_mat][()]

                next_syn += n

            f.close()

        self.out_file["network/" + h5_syn_n][0] = next_syn

        end_time2 = timeit.default_timer()
        self.write_log("Parallel pruning + merging: " + str(end_time2 - start_time))

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

            while start_row < range_end and (synapses[block_start[idx] - 1, 0:2] == synapses[start_row, 0:2]).all():
                start_row += 1

            block_start[idx] = start_row

        # If we have few synapses, and many workers, there might be some workers
        # assigned to the same range. Remove those, and pad with None

        block_start = [x for x in block_start if x <= range_end]
        block_start = list(collections.OrderedDict.fromkeys(block_start))

        while len(block_start) < num_workers + 1:
            block_start.append(None)

        synapse_range = [x for x in zip(block_start[:-1], block_start[1:])]

        self.write_log("synapse_range=" + str(synapse_range))

        return synapse_range

    ############################################################################

    def findRangesOLD(self, synapses, nWorkers, startPos=0, nSyn=None):

        assert False, "findRangesOLD obsolete, remove?"

        if nSyn is None:
            nSyn = synapses.shape[0] - startPos

        blockSize = max(1, int(math.floor(float(nSyn) / nWorkers)))

        self.write_log("Find block ranges. From " + str(startPos) \
                       + " to " + str(nSyn + startPos) \
                       + " block size " + str(blockSize))

        try:
            blockStart = np.array([x + startPos for x \
                                   in range(0, blockSize * nWorkers, blockSize)])
        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)
            import pdb
            pdb.set_trace()

        self.write_log("blockStart=" + str(blockStart))

        for idx in range(1, len(blockStart)):

            assert blockStart[idx - 1] < blockStart[idx], \
                "Blocks seem to be quite small, run in serial instead. " \
                + "idx = " + str(idx) + " blockStart = " + str(blockStart)

            startRow = blockStart[idx]

            # Need to make sure that the synapses between pairs of neurons are
            # not split between workers
            try:
                while (synapses[blockStart[idx] - 1, 0:2] == synapses[startRow, 0:2]).all():
                    startRow += 1
            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr)

                self.write_log("Problem with setting ranges, fall back to serial")
                return None

            blockStart[idx] = startRow

        assert (np.diff(blockStart) > 0).all(), \
            "All workers should have different block starts. Try running in serial."

        # Calculate the corresponding block ends, special treatment for last pos
        blockEnd = np.zeros(blockStart.shape, dtype=int)
        blockEnd[0:-1] = blockStart[1:]
        blockEnd[-1] = startPos + nSyn + 1  # +1, python dont include element at end

        synapseRanges = [range(x, y) for x, y in zip(blockStart, blockEnd)]

        self.write_log("synapseRanges = " + str(synapseRanges))

        return synapseRanges

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
                print("Closing of file failed: " + str(f))

        self.temp_file_list = None

    ############################################################################

    def write_log(self, text, flush=True):  # Change flush to False in future, debug
        try:
            if self.logfile is not None:
                self.logfile.write(text + "\n")
                print(text)
                if flush:
                    self.logfile.flush()
            else:
                if self.verbose:
                    print(text)
        except:
            print(text)
            print("Unable to write to log file. Is log file closed?")

    ############################################################################

    def set_seed(self, random_seed):

        self.write_log("Setting random seed: " + str(random_seed))
        np.random.seed(random_seed)

    ############################################################################

    def new_worker_seeds(self, d_view):

        num_workers = len(self.d_view)
        worker_seeds = np.random.randint(0, np.iinfo(np.uint32).max,
                                         dtype=np.uint32,
                                         size=(num_workers,))
        self.d_view.scatter("worker_seed", worker_seeds, block=True)
        self.d_view.execute("nw.set_seed(worker_seed[0])", block=True)

        self.write_log("New worker seeds: " + str(worker_seeds))

    ############################################################################

    def setup_parallel(self, dView):

        assert self.role == "master", \
            "setupParallel: Should only be called by master node"

        if dView is None:
            self.write_log("setupParallel called without dView, aborting.")
            return

        if self.workers_initialised:
            self.write_log("Workers already initialised.")
            return

        with dView.sync_imports():
            from snudda.prune import SnuddaPrune  # Import on workers

        self.write_log("Setting up workers: " \
                       + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

        # Create unique log file names for the workers
        if self.logfile_name is not None:
            engine_log_file = [self.logfile_name + "-" + str(x) for x in range(0, len(dView))]
        else:
            engine_log_file = [[] for x in range(0, len(dView))]

        dView.scatter('logfile_name', engine_log_file, block=True)
        dView.push({"work_history_file": self.work_history_file}, block=True)

        cmd_str = "nw = SnuddaPrune(work_history_file=work_history_file, logfile_name=logfile_name[0],role='worker')"
        dView.execute(cmd_str, block=True)

        # Make sure we have different seeds for workers
        self.new_worker_seeds(dView)

        self.write_log("Workers setup: " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

        self.workers_initialised = True

    ############################################################################

    def low_memory(self, threshold=0.1) -> bool:
        with open('/proc/meminfo', 'r') as mem:
            ret = {}
            tmp = 0

            for i in mem:
                sline = i.split()
                if str(sline[0]) == 'MemTotal:':
                    ret['total'] = int(sline[1])
                elif str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
                    tmp += int(sline[1])
            ret['free'] = tmp
            ret['used'] = int(ret['total']) - int(ret['free'])

        # import pdb
        # pdb.set_trace()

        memory_ratio = ret['free'] / ret['total']

        # self.writeLog("Memory status: " + str(int(memoryRatio * 100)) + "% free")

        return memory_ratio < threshold

    ############################################################################

    # This code does a N-way merge for synapses

    # !!! We also need to merge gap junctions

    def big_merge(self, merge_data_type="synapses", clean_voxel_files=True):

        assert False, "big_merge is replaced by big_merge_lookup"

        # Since we want the code to work for both synapses and gap junction
        # we need to know location of synapse matrix, eg "network/synapses",
        # number of synapses, eg "nSynapses", and in some functions which
        # columns containe x,y,z voxel coords, eg. range(2,5)
        h5SynMat, h5SynN, h5SynLoc = self.data_loc[merge_data_type]

        self.merge_data_type = merge_data_type
        self.write_log("Doing bigMerge for " + merge_data_type)

        synapse_heap = []

        self.write_log("Starting big merge of all data")

        nNeurons = len(self.hist_file["network/neurons/neuronID"])
        assert np.max(self.hist_file["network/neurons/neuronID"]) + 1 == nNeurons, \
            "bigMerge: There are neuron IDs missing"

        maxHyperID = np.max(self.all_hyper_i_ds) + 1
        fileList = [None] * maxHyperID
        numSynapses = np.zeros((maxHyperID,), dtype=np.int)

        # Open all files for reading
        hFileNameMask = self.base_path + "/voxels/network-putative-synapses-%s.hdf5"

        maxAxonVoxelCtr = 0
        maxDendVoxelCtr = 0

        nSynHist = self.hist_file[h5SynN]
        nSynTotal = np.sum(nSynHist)

        for hID, nSyn, nOverflow in zip(self.hist_file["completed"], nSynHist,
                                        self.hist_file["voxelOverflowCounter"]):

            if nSyn > 0:
                # Open file, and add first row to the heap
                hFileName = hFileNameMask % str(hID)
                fileList[hID] = h5py.File(hFileName, 'r')
                numSynapses[hID] = nSyn

                srcID = fileList[hID][h5SynMat][0, 0]
                destID = fileList[hID][h5SynMat][0, 1]
                uniqueID = destID * nNeurons + srcID  # This is used for heap priority

                # Create a heap containing the first element of all files
                heapq.heappush(synapse_heap, (uniqueID, hID))

                # This is so we can optimize the axon/dend voxelCtr and size
                if "maxAxonVoxelCtr" in fileList[hID]["meta"]:
                    maxAxonVoxelCtr = max(maxAxonVoxelCtr,
                                          fileList[hID]["meta/maxAxonVoxelCtr"][()])
                if "maxDendVoxelCtr" in fileList[hID]["meta"]:
                    maxDendVoxelCtr = max(maxDendVoxelCtr,
                                          fileList[hID]["meta/maxDendVoxelCtr"][()])

                if clean_voxel_files:
                    # This makes sure we remove the old voxel files afterwards
                    self.temp_file_list.append(hFileName)

        # Setup read buffer
        self.setupBufferedMergeRead(fileList, h5SynMat)

        if self.buffer_out_file is None:
            # Create output file
            (self.buffer_out_file, outFileName) \
                = self.setup_merge_file(delete_after=False)
        else:
            # We need to reset the write pointer (GJ and synapses should start from 0)
            self.next_file_write_pos = 0

        # Only save this meta data if doing the synapses call
        if maxAxonVoxelCtr > 0 and self.merge_data_type == "synapses":
            self.buffer_out_file["meta"].create_dataset("maxAxonVoxelCtr",
                                                        data=maxAxonVoxelCtr)
            self.write_log("maxAxonVoxelCtr = " + str(maxAxonVoxelCtr))

        if maxDendVoxelCtr > 0 and self.merge_data_type == "synapses":
            self.buffer_out_file["meta"].create_dataset("maxDendVoxelCtr",
                                                        data=maxDendVoxelCtr)
            self.write_log("maxDendVoxelCtr = " + str(maxDendVoxelCtr))

        # 2. Pop the smallest element, and add it to the final file
        # -- check same file if there are more with the same source and dest
        # -- buffer the writes

        synCtr = 0

        if len(synapse_heap) > 0:
            # Get the first file to read synapses from
            (uniqueID, hID) = heapq.heappop(synapse_heap)
        else:
            # No synapses at all, return
            self.clean_up_merge_read_buffers()

            return self.buffer_out_file[h5SynMat], self.buffer_out_file

        # synapses coordinates are translated to simulation wide coordinates
        # from hyper voxel local coordinates
        (synapses, synapsesRemaining, nextSynapsePair) \
            = self.findPairSynapsesBufferedMerge(hID, h5SynMat, h5SynLoc)

        self.buffer_merge_write(h5SynMat, synapses)
        synCtr += synapses.shape[0]

        # !!! DEBUG
        prevPair = synapses[0, 0:2].copy()
        prevUID = uniqueID

        loopCtr = 0
        while synapsesRemaining or len(synapse_heap) > 0:

            if loopCtr % 1000000 == 0:
                self.write_log("Synapses: " + str(synCtr) \
                               + "/" + str(nSynTotal) \
                               + " (heap size: " + str(len(synapse_heap)) + ")")

            if synapsesRemaining:
                # Need to push the next synapse in line from that file onto the heap
                [srcID, destID] = nextSynapsePair

                uniqueID = destID * nNeurons + srcID  # This is used for heap priority
                (uniqueID, hID) = heapq.heappushpop(synapse_heap, (uniqueID, hID))
            else:
                (uniqueID, hID) = heapq.heappop(synapse_heap)

            (synapses, synapsesRemaining, nextSynapsePair) \
                = self.findPairSynapsesBufferedMerge(hID, h5SynMat, h5SynLoc)

            try:
                assert nextSynapsePair is None or \
                       synapses[0, 1] < nextSynapsePair[1] \
                       or (synapses[0, 1] == nextSynapsePair[1] \
                           and synapses[0, 0] < nextSynapsePair[0]), \
                    "They are not in order in file"

            except:
                print("This is realyl weird...")
                import pdb
                pdb.set_trace()

            assert uniqueID == synapses[0, 1] * nNeurons + synapses[0, 0], \
                "Oh no! Internal inconsistency"

            # !!! Just a check while writing code, to make sure it is ok
            assert (synapses[:, 0] == synapses[0, 0]).all() \
                   and (synapses[:, 1] == synapses[0, 1]).all(), \
                "Source and dest should all be the same for pair synapses"

            try:
                assert (synapses[0, 1] > prevPair[1]
                        or (synapses[0, 1] == prevPair[1] and synapses[0, 0] >= prevPair[0])), \
                    "Synapses not in order"
            except:
                print("We have a big problem! Synapses not in order!")
                import pdb
                pdb.set_trace()

            try:
                self.buffer_merge_write(h5SynMat, synapses)
                synCtr += synapses.shape[0]
                loopCtr += 1

            except:
                print("Why is this wrong?")
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr)
                import pdb
                pdb.set_trace()

                # DEBUG!!!
            prevPair = synapses[0, 0:2].copy()
            prevUID = uniqueID

        # Make sure all of the buffer is written to file
        self.buffer_merge_write(h5SynMat, flush=True)

        # The buffers should be cleared once the files are read from,
        # but just as a precaution
        self.clean_up_merge_read_buffers()

        return self.buffer_out_file[h5SynMat], self.buffer_out_file

    ############################################################################

    def findPairSynapses(self, synapses, startIdx):

        endIdx = startIdx + 1

        while (endIdx < synapses.shape[0]
               and synapses[startIdx, 0] == synapses[endIdx, 0]
               and synapses[startIdx, 1] == synapses[endIdx, 1]):
            endIdx += 1

        if endIdx < synapses.shape[0]:
            synapsesRemaining = True
        else:
            synapsesRemaining = False

        return synapses[startIdx:endIdx, :], endIdx, synapsesRemaining

    ############################################################################

    def setupBufferedMergeRead(self, fileList, synapseMatrixPath):

        self.file_list = fileList
        self.file_buffers = [None] * len(fileList)
        self.nextBufferReadPos = np.zeros((len(fileList),), dtype=int)
        self.nextFileReadPos = np.zeros((len(fileList),), dtype=int)

        for idx, f in enumerate(fileList):
            if f is not None:
                nSyn = f[synapseMatrixPath].shape[0]
                bufLen = min(nSyn, self.synapse_buffer_size)
                self.file_buffers[idx] = f[synapseMatrixPath][0:bufLen, :]
                self.nextFileReadPos[idx] = bufLen

    ############################################################################

    def updateMergeReadBuffer(self, hID, synapseMatrixPath):

        if self.file_buffers[hID] is None:
            # Nothing to do
            return

        if self.file_list[hID] is None:
            self.file_buffers[hID] = None
            return

        syn = self.file_list[hID][synapseMatrixPath]
        nSynToRead = min(syn.shape[0] - self.nextFileReadPos[hID],
                         self.synapse_buffer_size)

        if nSynToRead > 0:
            startPos = self.nextFileReadPos[hID]
            endPos = startPos + nSynToRead

            if endPos - startPos != self.file_buffers[hID].shape[0]:
                self.file_buffers[hID] = syn[startPos:endPos, :].copy()
            else:
                # Reuse memory
                self.file_buffers[hID][:, :] = syn[startPos:endPos, :]

            self.nextFileReadPos[hID] = endPos
            self.nextBufferReadPos[hID] = 0

        else:
            self.file_buffers[hID] = None
            self.nextBufferReadPos[hID] = 0

            ############################################################################

    def clean_up_merge_read_buffers(self):

        self.merge_data_type = None

        for hid, f in enumerate(self.file_list):

            if f is not None:
                try:
                    f.close()
                except:
                    self.write_log("Problem closing file for HID: " + str(hid))

        if self.file_buffers is not None:
            for idx in range(0, len(self.file_buffers)):
                self.file_buffers[idx] = None

        self.file_buffers = None
        self.nextFileReadPos = None
        self.nextBufferReadPos = None

    ############################################################################

    # Input: Hyper voxel ID, columns which contains location
    # (columns are different for synapses: range(2,5) and GJ range(6,9)
    # so we take them as parameter to reuse code for both cases

    # Returns: synapses between pair,
    #          synapses remaining flag,
    #          (srcID,destID) of next synapse
    #
    # OBS, voxel coordinates for synapses are translated from hyper voxel
    # coordiantes, to simulation wide coordinates

    def findPairSynapsesBufferedMerge(self, hID, synapseMatrixPath, locationColumns):

        startIdx = self.nextBufferReadPos[hID]
        endIdx = startIdx + 1

        synapses = self.file_buffers[hID]
        bufLen = synapses.shape[0]

        oldSynapses = None

        while True:

            if endIdx >= bufLen:

                # End of buffer, we need to read more synapses from file
                assert oldSynapses is None, \
                    "Buffer too small, should not overflow twice."  # self.synapseBufferSize

                oldSynapses = synapses[startIdx:, :].copy()

                # This resets self.nextBufferReadPos[hID]
                self.updateMergeReadBuffer(hID, synapseMatrixPath)

                startIdx = 0
                endIdx = 1

                if self.file_buffers[hID] is None:
                    # This was the end of the file
                    return oldSynapses, False, None

                synapses = self.file_buffers[hID]
                bufLen = synapses.shape[0]

                assert oldSynapses.shape[0] > 0, "This should never be empty"

                if (oldSynapses[0, 0] != synapses[0, 0]
                        or oldSynapses[0, 1] != synapses[0, 1]):
                    # New pair in new buffer, send old pair synapses to caller
                    return oldSynapses, True, synapses[0, :2].copy()

            # Check if the next synapse belongs to pair
            if (synapses[startIdx, 0] != synapses[endIdx, 0]
                    or synapses[startIdx, 1] != synapses[endIdx, 1]):
                # No more similar, stop
                break

            endIdx += 1

        # Update the buffer read pos
        self.nextBufferReadPos[hID] = endIdx

        if endIdx < bufLen:
            synapsesRemaining = True
            nextSynapsePair = synapses[endIdx, 0:2].copy()
        else:
            synapsesRemaining = False
            nextSynapsePair = None

        if oldSynapses is not None:
            syn = np.concatenate((oldSynapses, synapses[startIdx:endIdx, :]))
            # Need to add first half also from previous buffer read
        else:
            syn = synapses[startIdx:endIdx, :]

        # OBS, we need to convert to global voxel coordinates from the local
        # hyper voxel coordinates

        # !!! We need to convert from hyper voxel local coordinates for synapse
        # to simulation wide coordinates
        # !!! THIS IS NOW DONE AFTER TOUCH DETECTION
        # syn[:,locationColumns] += self.hyperVoxelOffset[hID]

        return syn, synapsesRemaining, nextSynapsePair

    ############################################################################

    def buffer_merge_write(self, synapseMatrixLoc, synapses=None, flush=False):

        #    if(synapseMatrixLoc == "network/gapJunctions" \
        #       and flush):
        #      import pdb # !!DEBUG
        #      pdb.set_trace()

        if self.synapse_write_buffer is None:
            # First time run

            if synapses is not None:
                bufferWidth = synapses.shape[1]
            else:
                self.write_log("bufferMergeWrite: Unable to guess width of matrix")
                import pdb
                pdb.set_trace()

            self.synapse_write_buffer = np.zeros((self.synapse_buffer_size,
                                                  bufferWidth),
                                                 dtype=np.int32)
            self.next_buffer_write_pos = 0
            # self.nextFileWritePos = 0

        if synapses is not None:

            # Is buffer almost full?
            if synapses.shape[0] + self.next_buffer_write_pos >= self.synapse_buffer_size:
                # Flush buffer first
                endFileIdx = self.next_file_write_pos + self.next_buffer_write_pos
                self.buffer_out_file[synapseMatrixLoc][self.next_file_write_pos:endFileIdx, :] \
                    = self.synapse_write_buffer[0:self.next_buffer_write_pos, :]

                self.next_file_write_pos += self.next_buffer_write_pos
                self.next_buffer_write_pos = 0
                # Ok, now we have clean buffer, back to normal program

            endBufIdx = self.next_buffer_write_pos + synapses.shape[0]
            self.synapse_write_buffer[self.next_buffer_write_pos:endBufIdx, :] = synapses
            self.next_buffer_write_pos += synapses.shape[0]

        if flush:
            self.write_log("Flushing " + str(self.buffer_out_file.filename) \
                           + " data: " + synapseMatrixLoc)
            endFileIdx = self.next_file_write_pos + self.next_buffer_write_pos
            self.buffer_out_file[synapseMatrixLoc][self.next_file_write_pos:endFileIdx, :] \
                = self.synapse_write_buffer[0:self.next_buffer_write_pos, :]

            self.next_file_write_pos += self.next_buffer_write_pos
            self.next_buffer_write_pos = 0

            # Resize matrix to fit data
            w = self.buffer_out_file[synapseMatrixLoc].shape[1]
            self.buffer_out_file[synapseMatrixLoc].resize((endFileIdx, w))
            self.write_log(synapseMatrixLoc + " new size " + str((endFileIdx, w)))

            # Remove buffer
            self.synapse_write_buffer = None

    ############################################################################

    def big_merge_parallel(self):

        if self.role != "master":
            self.write_log("bigMergeParallel is only run on master node, aborting")
            return

        self.write_log("bigMergeParallel, starting " + str(self.role))

        if self.d_view:
            self.setup_parallel(dView=self.d_view)

        # Split neurons between nodes, we need the neurons to be in order
        num_neurons = self.hist_file["network/neurons/neuronID"].shape[0]
        assert num_neurons - 1 == self.hist_file["network/neurons/neuronID"][-1], \
            "neuronID should start from 0 and the end should be n-1"

        nWorkers = len(self.d_view)

        neuron_ranges = []
        range_borders = np.linspace(0, num_neurons, nWorkers + 1).astype(int)

        for idx in range(0, nWorkers):
            neuron_ranges.append((range_borders[idx], range_borders[idx + 1]))

        assert neuron_ranges[-1][-1] == num_neurons, \
            "bigMergeParallel: Problem with neuronRanges, last element incorrect"
        assert len(neuron_ranges) == nWorkers, \
            "bigMergeParallel: Problem with neuronRanges, bad length"

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
        print("Processing merge...")

        merge_start_syn = [x[1][0] for x in merge_results_syn]
        merge_start_gj = [x[1][0] for x in merge_results_gj]

        merge_order_syn = np.argsort(merge_start_syn)
        merge_order_gj = np.argsort(merge_start_gj)

        # We then need the file to merge them in
        (self.buffer_out_file, outFileName) \
            = self.setup_merge_file(big_cache=False, delete_after=False)

        # import pdb
        # pdb.set_trace()

        # Copy the data to the file

        for order, results, location \
                in zip([merge_order_syn, merge_order_gj],
                       [merge_results_syn, merge_results_gj],
                       ["network/synapses", "network/gapJunctions"]):

            start_pos = 0

            for idx in order:

                data_file = results[idx][0]
                num_synapses = results[idx][2]
                end_pos = start_pos + num_synapses

                if data_file is None:
                    assert num_synapses == 0, "!!! Missing merge file " + str(data_file) \
                                              + ", internal problem"
                    continue

                self.write_log("Extracting " + location + " from " + data_file)

                f_in = h5py.File(data_file, 'r')

                # Some idiot checks...
                assert f_in[location].shape[0] == num_synapses, \
                    "Internal inconsistency in number of rows stored"
                assert not (f_in[location][-1, :] == 0).all(), \
                    "Last row all zero, that should not happen"

                self.buffer_out_file[location][start_pos:end_pos, :] = \
                    f_in[location][:, :]

                f_in.close()
                start_pos = end_pos

        # We should check that the number of synapses in output file is correct

        assert self.buffer_out_file["network/synapses"].shape[0] \
               == np.sum(self.hist_file["nSynapses"][:]), \
            "Not all synapses kept in merge, internal problem."
        assert self.buffer_out_file["network/gapJunctions"].shape[0] \
               == np.sum(self.hist_file["nGapJunctions"][:]), \
            "Not all gap junctions kept in merge, internal problem."

        self.write_log("bigMergeParallel: done")

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

    def big_merge_helper(self, neuron_range, merge_data_type):

        try:
            self.write_log("Neuronrange = " + str(neuron_range))

            output_filename = self.scratchPath + merge_data_type + "-for-neurons-" \
                              + str(neuron_range[0]) + "-to-" + str(neuron_range[1]) \
                              + "-MERGE-ME.hdf5"

            self.merge_data_type = merge_data_type

            # Which hyper voxels are the neurons located in?
            hv_list = self.get_hyper_voxel_list(neuron_range)

            # Locations of data within the file
            h5_syn_mat, h5_syn_n, h5_syn_lookup = self.data_loc[merge_data_type]

            synapse_heap = []
            file_list = dict([])
            file_mat_iterator = dict([])

            # Next we need to open all the relevant files
            h_file_name_mask = self.base_path + "/voxels/network-putative-synapses-%s.hdf5"

            max_axon_voxel_ctr = 0
            max_dend_voxel_ctr = 0

            # !!! Special code to increase h5py cache size
            propfaid = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
            settings = list(propfaid.get_cache())
            settings[2] *= 20
            propfaid.set_cache(*settings)
            # settings = propfaid.get_cache()
            # !!! End special cache code

            n_hv = int(self.hist_file["nCompleted"][0])
            n_total = 0

            # import pdb
            # pdb.set_trace()

            for h_id, n_syn, n_overflow in zip(self.hist_file["completed"][:n_hv],
                                               self.hist_file[h5_syn_n][:n_hv],
                                               self.hist_file["voxelOverflowCounter"][:n_hv]):
                n_total += n_syn

                if h_id not in hv_list:
                    # We only need a subset of the files
                    continue

                if n_syn > 0:
                    h_filename = h_file_name_mask % str(h_id)
                    self.write_log("Opening voxel file: " + h_filename)

                    # Low level opening hdf5 file, to have greater cache size #ACC_RDWR
                    fid = h5py.h5f.open(h_filename.encode(),
                                        flags=h5py.h5f.ACC_RDONLY,
                                        fapl=propfaid)

                    # !!! Temp print to check cache size
                    settings = list(fid.get_access_plist().get_cache())
                    print(settings)

                    # fileList[hID] = h5py.File(hFileName,'r')
                    file_list[h_id] = h5py.File(fid, drive=self.h5driver)

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
                        max_axon_voxel_ctr = max(max_axon_voxel_ctr,
                                              file_list[h_id]["meta/maxAxonVoxelCtr"][()])
                    if "maxDendVoxelCtr" in file_list[h_id]["meta"]:
                        max_dend_voxel_ctr = max(max_dend_voxel_ctr,
                                              file_list[h_id]["meta/maxDendVoxelCtr"][()])

            if merge_data_type == "synapses":
                num_synapses = n_total
                num_gap_junctions = 0
            elif merge_data_type == "gapJunctions":
                num_synapses = 0
                num_gap_junctions = n_total
            else:
                assert False, "Unknown mergeDataType " + str(merge_data_type)

            # Setup output file
            (self.buffer_out_file, outFileName) \
                = self.setup_merge_file(big_cache=True,
                                        outfile_name=output_filename,
                                        save_morphologies=False,
                                        num_synapses=num_synapses,
                                        num_gap_junctions=num_gap_junctions)

            # Here we store the sorted connection matrix
            sorted_mat = self.buffer_out_file[h5_syn_mat]

            # Only save this meta data if doing the synapses call
            if max_axon_voxel_ctr > 0 and self.merge_data_type == "synapses":
                self.buffer_out_file["meta"].create_dataset("maxAxonVoxelCtr",
                                                            data=max_axon_voxel_ctr)
                self.write_log("maxAxonVoxelCtr = " + str(max_axon_voxel_ctr))

            if max_dend_voxel_ctr > 0 and self.merge_data_type == "synapses":
                self.buffer_out_file["meta"].create_dataset("maxDendVoxelCtr",
                                                            data=max_dend_voxel_ctr)
                self.write_log("maxDendVoxelCtr = " + str(max_dend_voxel_ctr))

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
                    self.write_log("Worker synapses: " + str(syn_ctr) \
                                   + "/" + str(n_total) \
                                   + " (heap size: " + str(len(synapse_heap)) + ")")

                # Get the next set of synapses from this file from the iterator
                next_row_set = next(file_mat_iterator[h_id], None)

                if next_row_set is not None:
                    # More synapses in file, push next pair to heap, and pop top pair
                    syn_set, unique_id = next_row_set
                    (unique_id, h_id, syn_set) = heapq.heappushpop(synapse_heap,
                                                                (unique_id, h_id,
                                                                 syn_set))
                elif len(synapse_heap) > 0:
                    (unique_id, h_id, syn_set) = heapq.heappop(synapse_heap)

                else:
                    done = True
                    continue

                try:
                    assert unique_id >= old_unique_id, "uniqueID should be increasing in file"
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)

                    import pdb
                    pdb.set_trace()

                self.buffer_merge_write(h5_syn_mat, syn_set)
                syn_ctr += syn_set.shape[0]
                loop_ctr += 1

            self.write_log("Worker synapses: " + str(syn_ctr) \
                           + "/" + str(n_total) \
                           + " (heap size: " + str(len(synapse_heap)) + ")")

            self.write_log("Read " + str(syn_ctr) + " out of total " + str(n_total) \
                           + " synapses")

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

                    print("Problems closing files")
                    import pdb
                    pdb.set_trace()

            self.buffer_out_file.close()

            return output_filename, neuron_range, syn_ctr
        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)
            self.write_to_random_file(tstr)
            import pdb
            pdb.set_trace()

        ############################################################################

    def big_merge_lookup(self, merge_data_type="synapses", clean_voxel_files=True):

        # Since we want the code to work for both synapses and gap junction
        # we need to know location of synapse matrix, eg "network/synapses",
        # number of synapses, eg "nSynapses", the lookup table to quickly
        # find which synapse rows belongs to each pair of connected neurons
        # eg "network/synapseLookup"
        h5_syn_mat, h5_syn_n, h5_syn_lookup = self.data_loc[merge_data_type]

        self.merge_data_type = merge_data_type
        self.write_log("Doing bigMerge (lookup) for " + merge_data_type)

        synapse_heap = []

        num_neurons = len(self.hist_file["network/neurons/neuronID"])
        assert np.max(self.hist_file["network/neurons/neuronID"]) + 1 == num_neurons, \
            "bigMerge (lookup): There are neuron IDs missing"

        max_hyper_id = np.max(self.all_hyper_i_ds) + 1
        file_list = [None] * max_hyper_id
        file_mat_iterator = [None] * max_hyper_id

        # fileMat = [None] * maxHyperID # points to the synapse matrix in each file
        # fileMatLookup = [None] * maxHyperID # points to the matrix lookup in file

        num_synapses = np.zeros((max_hyper_id,), dtype=np.int)

        # Open all files for reading
        h_file_name_mask = self.base_path + "/voxels/network-putative-synapses-%s.hdf5"

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
                self.write_log("Opening voxel file: " + h_file_name)

                # Low level opening hdf5 file, to have greater cache size #ACC_RDWR
                fid = h5py.h5f.open(h_file_name.encode(),
                                    flags=h5py.h5f.ACC_RDONLY,
                                    fapl=propfaid)

                # !!! Temp print to check cache size
                settings = list(fid.get_access_plist().get_cache())
                print(settings)

                # fileList[hID] = h5py.File(hFileName,'r')
                try:
                    file_list[h_id] = h5py.File(fid, drive=self.h5driver)
                    file_mat_iterator[h_id] \
                        = self.synapse_set_iterator(h5mat_lookup=file_list[h_id][h5_syn_lookup],
                                                    h5mat=file_list[h_id][h5_syn_mat],
                                                    chunk_size=10000)
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)

                    print("This should not happen...")
                    import pdb
                    pdb.set_trace()

                num_synapses[h_id] = nSyn

                # There should be at least the first row, otherwise nSyn = 0
                syn_set, unique_id = next(file_mat_iterator[h_id], None)

                # Create a heap containing the first subset of all files
                heapq.heappush(synapse_heap, (unique_id, h_id, syn_set))

                # This is so we can optimize the axon/dend voxelCtr and size
                if "maxAxonVoxelCtr" in file_list[h_id]["meta"]:
                    max_axon_voxel_ctr = max(max_axon_voxel_ctr,
                                             file_list[h_id]["meta/maxAxonVoxelCtr"][()])
                if "maxDendVoxelCtr" in file_list[h_id]["meta"]:
                    max_dend_voxel_ctr = max(max_dend_voxel_ctr,
                                             file_list[h_id]["meta/maxDendVoxelCtr"][()])

        assert np.sum(num_synapses) == num_syn_total, \
            "Mismatch between work log file and data files: " \
            + str(num_syn_total) + " vs " + str(np.sum(num_synapses)) + " synapses"

        if self.buffer_out_file is None:
            # Create output file
            (self.buffer_out_file, outFileName) \
                = self.setup_merge_file(big_cache=True, delete_after=False)
        else:
            # We need to reset the write pointer (GJ and synapses should start from 0)
            self.next_file_write_pos = 0

        # Here we store the sorted connection matrix
        sorted_mat = self.buffer_out_file[h5_syn_mat]

        # Only save this meta data if doing the synapses call
        if max_axon_voxel_ctr > 0 and self.merge_data_type == "synapses":
            self.buffer_out_file["meta"].create_dataset("maxAxonVoxelCtr",
                                                        data=max_axon_voxel_ctr)
            self.write_log("maxAxonVoxelCtr = " + str(max_axon_voxel_ctr))

        if max_dend_voxel_ctr > 0 and self.merge_data_type == "synapses":
            self.buffer_out_file["meta"].create_dataset("maxDendVoxelCtr",
                                                        data=max_dend_voxel_ctr)
            self.write_log("maxDendVoxelCtr = " + str(max_dend_voxel_ctr))

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
        # synEndIdx = synCtr + synSet.shape[0]
        # sortedMat[synCtr:synEndIdx,:] = synSet
        # synCtr = synEndIdx
        syn_ctr = syn_set.shape[0]
        self.buffer_merge_write(h5_syn_mat, syn_set)

        loop_ctr = 0
        done = False

        while not done:

            if loop_ctr % 1000000 == 0:
                self.write_log("Synapses: " + str(syn_ctr) \
                               + "/" + str(num_syn_total) \
                               + " (heap size: " + str(len(synapse_heap)) + ")")

            # Get the next set of synapses from this file from the iterator
            next_row_set = next(file_mat_iterator[h_id], None)

            if next_row_set is not None:
                # More synapses in file, push next pair to heap, and pop top pair
                syn_set, unique_id = next_row_set
                (unique_id, h_id, syn_set) = heapq.heappushpop(synapse_heap,
                                                               (unique_id, h_id,
                                                                syn_set))
            elif len(synapse_heap) > 0:
                (unique_id, h_id, syn_set) = heapq.heappop(synapse_heap)

            else:
                done = True
                continue

            # Write synapses to file
            # synEndIdx = synCtr + synSet.shape[0]
            # sortedMat[synCtr:synEndIdx,:] = synSet
            # synCtr = synEndIdx
            self.buffer_merge_write(h5_syn_mat, syn_set)
            syn_ctr += syn_set.shape[0]
            loop_ctr += 1

            # assert uniqueID == synapses[0,1]*nNeurons+synapses[0,0], \
            #  "bigMergeLookup: Oh no! Internal inconsistency"

        # Flush the buffers to file
        # self.bufferOutFile.flush()
        self.buffer_merge_write(h5_syn_mat, flush=True)

        self.write_log("big_merge_lookup: done")

        return sorted_mat, self.buffer_out_file

    ############################################################################

    def bigMergeLookupNOCACHE(self, mergeDataType="synapses"):

        assert False, "Is this used anymore? Replaced by big_merge_lookup -- or used for big big files?"

        # Since we want the code to work for both synapses and gap junction
        # we need to know location of synapse matrix, eg "network/synapses",
        # number of synapses, eg "nSynapses", the lookup table to quickly
        # find which synapse rows belongs to each pair of connected neurons
        # eg "network/synapseLookup"
        h5SynMat, h5SynN, h5SynLookup = self.data_loc[mergeDataType]

        self.merge_data_type = mergeDataType
        self.write_log("Doing bigMerge (lookup) for " + mergeDataType)

        synapse_heap = []

        nNeurons = len(self.hist_file["network/neurons/neuronID"])
        assert np.max(self.hist_file["network/neurons/neuronID"]) + 1 == nNeurons, \
            "bigMerge (lookup): There are neuron IDs missing"

        maxHyperID = np.max(self.all_hyper_i_ds) + 1
        fileList = [None] * maxHyperID
        fileMat = [None] * maxHyperID  # points to the synapse matrix in each file
        fileMatLookup = [None] * maxHyperID  # points to the matrix lookup in file

        nextLookupIdx = np.zeros((maxHyperID,), dtype=np.int)
        lookupLength = np.zeros((maxHyperID,), dtype=np.int)
        numSynapses = np.zeros((maxHyperID,), dtype=np.int)

        # Open all files for reading
        hFileNameMask = self.base_path + "/voxels/network-putative-synapses-%s.hdf5"

        maxAxonVoxelCtr = 0
        maxDendVoxelCtr = 0

        nSynHist = self.hist_file[h5SynN]
        nSynTotal = np.sum(nSynHist)

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

        for h_id, nSyn, nOverflow in zip(self.hist_file["completed"], nSynHist,
                                         self.hist_file["voxelOverflowCounter"]):
            if nSyn > 0:
                # Open file, and add info about first pairs synapses to the heap
                hFileName = hFileNameMask % str(h_id)
                self.write_log("Opening voxel file: " + hFileName)

                # Low level opening hdf5 file, to have greater cache size #ACC_RDWR
                fid = h5py.h5f.open(hFileName.encode(),
                                    flags=h5py.h5f.ACC_RDONLY,
                                    fapl=propfaid)

                # !!! Temp print to check cache size
                settings = list(fid.get_access_plist().get_cache())
                print(settings)

                # fileList[hID] = h5py.File(hFileName,'r')
                fileList[h_id] = h5py.File(fid, drive=self.h5driver)
                fileMat[h_id] = fileList[h_id][h5SynMat]
                fileMatLookup[h_id] = fileList[h_id][h5SynLookup]
                numSynapses[h_id] = nSyn

                try:
                    unique_id, start_idx, end_idx = fileMatLookup[h_id][0, :]
                    nextLookupIdx[h_id] = 1
                    lookupLength[h_id] = fileMatLookup[h_id].shape[0]
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)

                    print("No what happened...")
                    import pdb
                    pdb.set_trace()

                # Create a heap containing the first element of all files
                heapq.heappush(synapse_heap, (unique_id, h_id, start_idx, end_idx))

                # This is so we can optimize the axon/dend voxelCtr and size
                if "maxAxonVoxelCtr" in fileList[h_id]["meta"]:
                    maxAxonVoxelCtr = max(maxAxonVoxelCtr,
                                          fileList[h_id]["meta/maxAxonVoxelCtr"][()])
                if "maxDendVoxelCtr" in fileList[h_id]["meta"]:
                    maxDendVoxelCtr = max(maxDendVoxelCtr,
                                          fileList[h_id]["meta/maxDendVoxelCtr"][()])

        assert np.sum(numSynapses) == nSynTotal, \
            "Mismatch between work log file and data files: " \
            + str(nSynTotal) + " vs " + str(np.sum(numSynapses)) + " synapses"

        if self.buffer_out_file is None:
            # Create output file
            (self.buffer_out_file, outFileName) \
                = self.setup_merge_file(big_cache=True)
        else:
            # We need to reset the write pointer (GJ and synapses should start from 0)
            self.next_file_write_pos = 0

        # Here we store the sorted connection matrix
        sortedMat = self.buffer_out_file[h5SynMat]

        # Only save this meta data if doing the synapses call
        if maxAxonVoxelCtr > 0 and self.merge_data_type == "synapses":
            self.buffer_out_file["meta"].create_dataset("maxAxonVoxelCtr",
                                                        data=maxAxonVoxelCtr)
            self.write_log("maxAxonVoxelCtr = " + str(maxAxonVoxelCtr))

        if maxDendVoxelCtr > 0 and self.merge_data_type == "synapses":
            self.buffer_out_file["meta"].create_dataset("maxDendVoxelCtr",
                                                        data=maxDendVoxelCtr)
            self.write_log("maxDendVoxelCtr = " + str(maxDendVoxelCtr))

        # 2. Pop the smallest element, and add it to the final file
        # -- check same file if there are more with the same source and dest
        # -- buffer the writes

        synCtr = 0

        if len(synapse_heap) > 0:
            # Get the first file to read synapses from
            (unique_id, h_id, start_idx, end_idx) = heapq.heappop(synapse_heap)
        else:
            # No synapses at all, return
            self.clean_up_merge_read_buffers()
            return sortedMat, self.buffer_out_file

        # Maybe store ref direct to synapses, avoiding one level of lookups
        synapses = fileMat[h_id][start_idx:end_idx, :]
        synEndIdx = synCtr + synapses.shape[0]

        sortedMat[synCtr:synEndIdx, :] = synapses
        synCtr = synEndIdx

        loopCtr = 0
        done = False

        while not done:

            if loopCtr % 10000 == 0:
                self.write_log("Synapses: " + str(synCtr) \
                               + "/" + str(nSynTotal) \
                               + " (heap size: " + str(len(synapse_heap)) + ")")

            if nextLookupIdx[h_id] < lookupLength[h_id]:
                # More synapses in file, push next pair to heap, and pop top pair
                [unique_id, start_idx, end_idx] = \
                    fileMatLookup[h_id][nextLookupIdx[h_id], :]
                nextLookupIdx[h_id] += 1

                (unique_id, h_id, start_idx, end_idx) = heapq.heappushpop(synapse_heap,
                                                                          (unique_id, h_id,
                                                                           start_idx, end_idx))
            elif len(synapse_heap) > 0:
                (unique_id, h_id, start_idx, end_idx) = heapq.heappop(synapse_heap)

            else:
                done = True
                continue

            # Write synapses to file
            synapses = fileMat[h_id][start_idx:end_idx, :]
            synEndIdx = synCtr + (end_idx - start_idx)  # synapses.shape[0]

            sortedMat[synCtr:synEndIdx, :] = synapses
            synCtr = synEndIdx
            loopCtr += 1

            # assert uniqueID == synapses[0,1]*nNeurons+synapses[0,0], \
            #  "bigMergeLookup: Oh no! Internal inconsistency"

        # Flush the buffers to file
        self.buffer_out_file.flush()

        self.write_log("bigMergeLookup NOCACHE: done")

        return sortedMat, self.buffer_out_file

    ############################################################################

    def prune_synapses(self, synapse_file, output_filename, row_range,
                       merge_data_type,
                       close_input_file=True,
                       close_out_file=True):

        try:
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
                self.write_log("pruneSynapses: No " + merge_data_type
                               + " skipping pruning")
                return

            self.write_log("pruneSynapses: synapseFile=" + str(synapse_file)
                           + ", outputFileName=" + str(output_filename)
                           + ", rowRange=" + str(row_range)
                           + " (" + merge_data_type + ")")

            if type(synapse_file) == str:
                self.write_log("Opening synapse file: " + synapse_file)
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

            self.write_log("blockRanges=" + str(block_ranges))

            self.setup_output_file(output_filename)  # Sets self.outFile

            for synRange in block_ranges:
                self.write_log("Pruning range : " + str(synRange))

                synapses = synapse_file[h5_syn_mat][synRange[0]:synRange[-1]]
                self.prune_synapses_helper(synapses=synapses,
                                           output_file=self.out_file,
                                           merge_data_type=merge_data_type)

        except:
            print("prune_synapses: Something went wrong... :/")

            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)

            import pdb
            pdb.set_trace()

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

        # Init some stats
        n_all_removed = 0
        n_some_removed = 0
        n_too_few_removed = 0
        n_dist_dep_pruning = 0
        n_too_many_removed = 0
        n_not_connected = 0

        old_pos = -1

        while next_read_pos < read_end_of_range:

            if old_pos == next_read_pos:
                print("pruneSynapsesHelper: Same again")
                import pdb
                pdb.set_trace()

            old_pos = next_read_pos

            # How many lines contain synapses between this pair of neurons
            read_end_idx = next_read_pos + 1
            while (read_end_idx < read_end_of_range and
                   (synapses[next_read_pos, 0:2] == synapses[read_end_idx, 0:2]).all()):
                read_end_idx += 1

            # Temp check
            assert ((synapses[next_read_pos:read_end_idx, 0] == synapses[next_read_pos, 0]).all()
                    and (synapses[next_read_pos:read_end_idx, 1] == synapses[next_read_pos, 1]).all()), \
                "pruneSynapsesHelper: Internal error, more than one neuron pair"

            # import pdb
            # pdb.set_trace()

            # Stats
            n_pair_synapses = read_end_idx - next_read_pos

            src_id = synapses[next_read_pos, 0]
            dest_id = synapses[next_read_pos, 1]

            if merge_data_type == "gapJunctions":
                # All are gap junctions
                synapse_type = 3
            else:
                synapse_type = synapses[next_read_pos, 6]

            con_id = (self.type_id_list[src_id], self.type_id_list[dest_id], synapse_type)

            if con_id in self.connectivity_distributions:

                # We have the option to separate between connections within a
                # population unit or not. If conInfo[1] != None then first
                # tuple is connection info within a population unit, and second item
                # is connection info between different population units
                con_info = self.connectivity_distributions[con_id]

                #
                if con_info[1] is None or self.populationUnitID[src_id] == self.populationUnitID[dest_id]:
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
            if a3 is not None and np.random.random() > a3:
                # Prune all synapses between pair, do not add to synapse file
                next_read_pos = read_end_idx
                # No need to update keepRowFlag since default set to 0

                # Stats
                n_all_removed += n_pair_synapses
                continue

            if dist_p is not None:
                # Distance dependent pruning, used for FS->MS connections

                # distP contains d (variable for distance to soma)
                d = synapses[next_read_pos:read_end_idx, 8] * 1e-6  # dendrite distance
                p = eval(dist_p)
                frac_flag = np.random.random(n_pair_synapses) < f1
                dist_flag = np.random.random(n_pair_synapses) < p

                keep_row_flag[next_read_pos:read_end_idx] = np.logical_and(frac_flag, dist_flag)

                n_frac = sum(frac_flag)
                n_some_removed += n_pair_synapses - n_frac
                n_dist_dep_pruning += n_frac - sum(keep_row_flag[next_read_pos:read_end_idx])

            else:
                keep_row_flag[next_read_pos:read_end_idx] \
                    = np.random.random(n_pair_synapses) < f1
                n_some_removed += n_pair_synapses - sum(keep_row_flag[next_read_pos:read_end_idx])

            # Check if too many synapses, trim it down a bit
            n_keep = np.sum(keep_row_flag[next_read_pos:read_end_idx])

            if soft_max is not None and n_keep > soft_max:
                # pKeep = float(softMax)/nKeep # OLD implementation
                soft_max = float(soft_max)
                # pKeep = 2*softMax*np.divide(1-np.exp(-nKeep/softMax),1+np.exp(-nKeep/softMax))/nKeep
                p_keep = np.divide(2 * soft_max, (1 + np.exp(-(n_keep - soft_max) / 5)) * n_keep)

                keep_row_flag[next_read_pos:read_end_idx] = \
                    np.logical_and(p_keep > np.random.random(n_pair_synapses),
                                   keep_row_flag[next_read_pos:read_end_idx])

                # Stats
                n_too_many_removed += n_keep - sum(keep_row_flag[next_read_pos:read_end_idx])

                # Update count
                n_keep = np.sum(keep_row_flag[next_read_pos:read_end_idx])

            # If too few synapses, remove all synapses
            if mu2 is not None:
                p_mu = 1.0 / (1.0 + np.exp(-8.0 / mu2 * (n_keep - mu2)))

                if p_mu < np.random.random():
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

        self.write_log("Number of synapses removed where synapse connection not allowed: " + str(n_not_connected) \
                       + "\nNumber of synapses removed due to distance dependent pruning: " + str(n_dist_dep_pruning) \
                       + "\nNumber of synapses removed randomly: " + str(n_some_removed) \
                       + "\nNumber of synapses removed due to too many synapses between connected pair: " + str(
            n_too_many_removed) \
                       + "\nNumber of synapses removed due to too few synapses between connected pairs: " + str(
            n_too_few_removed) \
                       + "\nNumber of synapses removed where all synapses between pairs are removed: " \
                       + str(n_all_removed))

    ############################################################################

    def file_row_lookup_iterator(self, h5mat, chunk_size=10000):

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

    # minDestID and maxDestID are inclusive, only synapses with destID in that
    # range are iterated over

    def file_row_lookup_iterator_subset(self, h5mat_lookup, min_dest_id, max_dest_id, chunk_size=10000):

        num_neurons = self.hist_file["network/neurons/neuronID"].shape[0]
        min_unique_id = min_dest_id * num_neurons
        max_unique_id = max_dest_id * num_neurons

        self.write_log("minUniqueID: " + str(min_unique_id) + ", maxUniqueID: " + str(max_unique_id))

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

    def synapse_set_iterator(self, h5mat_lookup, h5mat,
                             chunk_size=10000,
                             lookup_iterator=None):

        # Allow the user to set an alternative lookupIterator if we only
        # want to iterate over a subset of the synapses
        if not lookup_iterator:
            lookup_iterator = self.file_row_lookup_iterator(h5mat_lookup,
                                                            chunk_size=chunk_size)

        mat_size = h5mat.shape[0]
        if mat_size < chunk_size:
            chunk_size = mat_size

        old_synapses = None
        # readBuffer = np.zeros((chunkSize,h5mat.shape[1]),dtype=h5mat.dtype)
        read_buffer = h5mat[:chunk_size, :].copy()
        buffer_start = 0  # What file pos does start of buffer correspond to
        buffer_end = chunk_size  # What file pos does end of buffer correspond to (+1)

        next_row_set = next(lookup_iterator, None)

        while next_row_set is not None:

            # startIdx and endIdx are the rows in the matrix we want to read between
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
                synMat = np.concatenate([old_synapses,
                                         read_buffer[:(end_idx - buffer_start), :]],
                                        axis=0)

                try:
                    assert end_idx == buffer_end \
                           or (read_buffer[start_idx - buffer_start, :2]
                               != read_buffer[end_idx - buffer_start, :2]).any(), \
                        "We missed one synpase! (2)"

                    assert (synMat[:, 0] == synMat[0, 0]).all() \
                           and (synMat[:, 1] == synMat[0, 1]).all(), \
                        f"Synapse matrix (2) contains more than one pair:\n{synMat}"

                    assert synMat.shape[0] == end_idx - start_idx, \
                        "Synapse matrix has wrong size"
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)

                    import pdb
                    pdb.set_trace()

                yield (synMat,
                       unique_id)

                old_synapses = None

            else:

                synMat = read_buffer[(start_idx - buffer_start):(end_idx - buffer_start), :]

                try:
                    assert end_idx == buffer_end \
                           or (read_buffer[start_idx - buffer_start, :2]
                               != read_buffer[end_idx - buffer_start, :2]).any(), \
                        "We missed one synpase! (1)"

                    assert (synMat[:, 0] == synMat[0, 0]).all() \
                           and (synMat[:, 1] == synMat[0, 1]).all(), \
                        f"Synapse matrix (1) contains more than one pair:\n{synMat}"

                    assert synMat.shape[0] == end_idx - start_idx, \
                        "Synapse matrix has wrong size"
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)

                    import pdb
                    pdb.set_trace()

                yield (synMat,
                       unique_id)

            next_row_set = next(lookup_iterator, None)


##############################################################################

if __name__ == "__main__":
    print("Please do not call this file directly, use snudda.py")
    exit(-1)
