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

import os
import sys
import numpy as np
import scipy
from numba import jit
import math
import numexpr
import collections

import heapq  # Priority queue

import time
import timeit

import h5py
import json
import glob

from snudda.utils.numpy_encoder import NumpyEncoder

# from snudda.Neuron_morphology import NeuronMorphology


# The x,y,z coordinates of the synapses in the original file are integer
# indexes referring to within the hyper voxel. Here they are transformed
# to refer to within the entire simulation, with the simulation origo as
# the origo (bin size unchanged)

# This code combines multiple voxel files
# !!! This also performs the pruning, since it requires knowledge of all
#     synapses between the neurons, something not possible within a hyper voxel
#     if the neuron crosses borders
from snudda.utils.snudda_path import snudda_parse_path


class SnuddaPrune(object):

    """ Prunes detected network synapses to match pairwise connectivity statistics. """

    ############################################################################

    def __init__(self,
                 network_path,
                 logfile=None, logfile_name=None,
                 rc=None, d_view=None, role="master", verbose=False,
                 config_file=None,  # Default is to use same config_file as for detect, but you can override it
                 scratch_path=None,
                 # pre_merge_only=False,
                 h5libver="latest",
                 random_seed=None,
                 keep_files=False,   # If True then you can redo pruning multiple times without reruning detect
                 all_neuron_pair_synapses_share_parameter_id=True):

        """
        Constructor.

        Args:
            network_path (str): Network directory
            logfile : File pointer to logfile, derived from logfile_name if not given
            logfile_name (str): Path to logfile, if given will set logfile
            rc : ipyparallel RemoteClient, used for parallel execution
            d_view : ipyparallel direct view, derived from rc if not given
            role (str, optional) : "master" or "worker"
            verbose (bool) : Verbose output, default False
            config_file (str, optional): Path to config file, default is to use same config file as for detect,
                                         but the user can override it. Useful when testing different pruning values
                                         in combination with keep_files=True.
            scratch_path (str, optional): Path to scratch directory for temporary files
            h5libver (str, optional) : Default "latest"
            random_seed (int, optinoal): Random seed for pruning
            keep_files (bool, optional): If True then you can redo pruning multiple times without reruning detect
            all_neuron_pair_synapses_share_parameter_id (bool): Instead of each synapse having a unique parameter_id
                                                                all synapses between the same neuron pair will have
                                                                the same parameter id.
        """

        self.rc = rc
        self.d_view = d_view

        self.work_history_file = os.path.join(network_path, "log", "network-detect-worklog.hdf5")    
        self.network_path = network_path
        self.keep_files = keep_files
        self.merge_info_file = os.path.join(self.network_path, "pruning_merge_info.json")
        self.all_neuron_pair_synapses_share_parameter_id = all_neuron_pair_synapses_share_parameter_id

        self.logfile = logfile
        self.verbose = verbose
        self.h5libver = h5libver

        if logfile_name:
            self.logfile_name = logfile_name
        elif logfile is not None:
            self.logfile_name = logfile.name
        else:
            self.logfile_name = os.path.join(self.network_path, "log", "synapse-pruning.txt")

        if self.logfile is None and self.logfile_name is not None:
            self.logfile = open(self.logfile_name, 'w')
            self.write_log(f"Log file {self.logfile_name} created.")

        self.random_seed = random_seed

        self.h5driver = "sec2"  # "core" # "stdio", "sec2"

        self.write_log(f"Using hdf5 driver {self.h5driver}, {self.h5libver} version")

        if self.rc and not self.d_view:
            self.d_view = self.rc.direct_view(targets='all')

        self.role = role
        self.workers_initialised = False

        # TODO: Move this to external config file?
        self.synapse_type_lookup = {1: "GABA",
                                    2: "AMPA_NMDA",
                                    3: "GapJunction",
                                    4: "ACh",
                                    5: "NO"}

        self.synapse_type_reverse_lookup = \
            {v: k for k, v in self.synapse_type_lookup.items()}

        self.max_channel_type = None

        self.projection_synapse_file = os.path.join(self.network_path, "network-projection-synapses.hdf5")

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
        self.num_projection_synapses = None
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

        # Used by get_neuron_random_seeds
        self.cache_neuron_seeds = None
        self.cache_old_seed = None
        self.cache_num_neurons = None

        self.open_work_history_file(work_history_file=self.work_history_file, config_file=config_file)

        self.set_scratch_path(scratch_path)
        self.load_pruning_information(config_file=config_file)

        # (locationOfMatrix,locationOfN,locationOfCoords)
        self.data_loc = {"synapses": ("network/synapses",
                                      "nHypervoxelSynapses",
                                      "nSynapses",
                                      "network/synapseLookup"),  # range(2,5)),
                         "gapJunctions": ("network/gapJunctions",
                                          "nHypervoxelGapJunctions",
                                          "nGapJunctions",
                                          "network/gapJunctionLookup")}  # range(6,9))}

    ############################################################################

    def prune(self):

        """ Merges output files from detect.py and prunes the synapses. Writes network-synapses.hdf5 """

        start_time = timeit.default_timer()

        if self.role != "master":
            self.write_log("prune should only be called on master")
            return

        merge_info = self.get_merge_info()
        if merge_info:
            merge_files_syn, merge_neuron_range_syn, merge_syn_ctr, \
                merge_files_gj, merge_neuron_range_gj, merge_gj_ctr = merge_info
        else:
            # From the hyper voxels gather all synapses (and gap junctions) belonging to specific neurons
            merge_files_syn, merge_neuron_range_syn, merge_syn_ctr, \
                merge_files_gj, merge_neuron_range_gj, merge_gj_ctr = self.gather_neuron_synapses()

            self.save_merge_info(merge_files_syn=merge_files_syn,
                                 merge_neuron_range_syn=merge_neuron_range_syn,
                                 merge_syn_ctr=merge_syn_ctr,
                                 merge_files_gj=merge_files_gj,
                                 merge_neuron_range_gj=merge_neuron_range_gj,
                                 merge_gj_ctr=merge_gj_ctr)

        # Prune synapses and gap junctions
        self.prune_synapses_parallel(synapse_file=merge_files_syn,
                                     synapse_ctr=merge_syn_ctr,
                                     merge_data_type="synapses",
                                     close_input_file=False)

        self.prune_synapses_parallel(synapse_file=merge_files_gj,
                                     synapse_ctr=merge_gj_ctr,
                                     merge_data_type="gapJunctions",
                                     close_input_file=True)

        end_time = timeit.default_timer()

        if not self.keep_files:
            self.cleanup()

        self.write_log(f"prune synapses and gap junctions: {end_time - start_time:.1f}s")

        # Close output file
        try:
            if self.out_file:
                self.out_file.close()
                self.out_file = None
        except:
            self.write_log("Problem closing out file - already closed?")


    ############################################################################

    def save_merge_info(self,
                        merge_files_syn, merge_neuron_range_syn, merge_syn_ctr,
                        merge_files_gj, merge_neuron_range_gj, merge_gj_ctr):

        """
        Writes merge info to file (pruning_merge_info.json), so prune can be rerun without having
        to re-merge detection hyper voxel files.

        Args:
            merge_files_syn : List of files containing synapses
            merge_neuron_range_syn : List of neuron ranges in each file
            merge_syn_ctr : List of synapse counts for each file
            merge_files_gj : List of files containing gap junctions
            merge_neuron_range_gj : List of neuron ranges in each file
            merge_gj_ctr : List of gap junction count for each file

        """

        data = collections.OrderedDict()

        data["merge_files_syn"] = merge_files_syn
        data["merge_neuron_range_syn"] = merge_neuron_range_syn
        data["merge_syn_ctr"] = merge_syn_ctr
        data["merge_files_gj"] = merge_files_gj
        data["merge_neuron_range_gj"] = merge_neuron_range_gj
        data["merge_gj_ctr"] = merge_gj_ctr

        with open(self.merge_info_file, "w") as f:
            json.dump(data, f, indent=4, cls=NumpyEncoder)

    def get_merge_info_helper(self):

        """
        Helper function for get_merge_info. Reads merge info written by save_merge_info.
        Please use get_merge_info, which also includes sanity checks on the data returned.

        Returns:
            tuple of data:
                merge_files_syn : List of files containing synapses
                merge_neuron_range_syn : List of neuron ranges in each file
                merge_syn_ctr : List of synapse counts for each file
                merge_files_gj : List of files containing gap junctions
                merge_neuron_range_gj : List of neuron ranges in each file
                merge_gj_ctr : List of gap junction count for each file

        """

        with open(self.merge_info_file, "r") as f:
            data = json.load(f, object_pairs_hook=collections.OrderedDict)

        merge_files_syn = data["merge_files_syn"]
        merge_neuron_range_syn = data["merge_neuron_range_syn"]
        merge_syn_ctr = data["merge_syn_ctr"]
        merge_files_gj = data["merge_files_gj"]
        merge_neuron_range_gj = data["merge_neuron_range_gj"]
        merge_gj_ctr = data["merge_gj_ctr"]

        return merge_files_syn, merge_neuron_range_syn, merge_syn_ctr, \
            merge_files_gj, merge_neuron_range_gj, merge_gj_ctr

    def get_merge_info(self):

        """
        Reads merge info written by save_merge_info.

        Returns:
            tuple of data:
                merge_files_syn : List of files containing synapses
                merge_neuron_range_syn : List of neuron ranges in each file
                merge_syn_ctr : List of synapse counts for each file
                merge_files_gj : List of files containing gap junctions
                merge_neuron_range_gj : List of neuron ranges in each file
                merge_gj_ctr : List of gap junction count for each file

        """

        if not os.path.exists(self.merge_info_file):
            return None

        try:
            merge_files_syn, merge_neuron_range_syn, merge_syn_ctr, \
                merge_files_gj, merge_neuron_range_gj, merge_gj_ctr = self.get_merge_info_helper()
        except:
            self.write_log(f"Problem readin merge info from {self.merge_info_file}")
            return None

        # Check that the merge info file is more recent than all the files it refer to

        merge_info_time = os.path.getmtime(self.merge_info_file)
        for f_name in merge_files_syn + merge_files_gj:
            if f_name is not None and not os.path.exists(f_name):
                self.write_log(f"Bad merge info: {f_name} referenced by {self.merge_info_file} does not exist.")
                return None

            if f_name is not None and os.path.getmtime(f_name) > merge_info_time:
                self.write_log(f"Bad merge info: {f_name} is newer than {self.merge_info_file}")
                return None

        # Check that synapse and gap junction counts in file match those in merge_info_file
        for f_name, syn_ctr, neuron_range in zip(merge_files_syn, merge_syn_ctr, merge_neuron_range_syn):

            if f_name is None:
                if syn_ctr > 0:  # Just double check that it was not supposed to have synapses
                    self.write_log(f"Bad merge info: File {f_name} was expected to have {syn_ctr} synapses")
                    return None

                continue

            with h5py.File(f_name, "r") as f:
                if f["network/synapses"].shape[0] != syn_ctr:
                    self.write_log(f"Bad merge info: {f_name} expected {syn_ctr} synapses, " 
                                   f"found {f['network/synapses'].shape[0]}")
                    return None

                if syn_ctr > 0:
                    # While we are at it, check that neuron range is what we expect
                    if neuron_range[0] > f["network/synapses"][0, 1] \
                            or f["network/synapses"][-1, 1] > neuron_range[1]:
                        self.write_log(f"Bad merge info: {f_name} has post synaptic neuron outside range "
                                       f"{neuron_range[0]}-{neuron_range[1]}")
                        return None

        for f_name, gj_ctr, neuron_range in zip(merge_files_gj, merge_gj_ctr, merge_neuron_range_gj):
            if f_name is None:
                if gj_ctr > 0:
                    self.write_log(f"Bad merge info: {f_name} was supposed to have {gj_ctr} gap junctions")
                    return None
                continue

            with h5py.File(f_name, "r") as f:
                if f["network/gapJunctions"].shape[0] != gj_ctr:
                    self.write_log(f"Bad merge info: {f_name} expected {gj_ctr} gap junctions, " 
                                   f"found {f['network/gapJunctions'].shape[0]}")
                    return None

                if gj_ctr > 0:
                    if neuron_range[0] > f["network/gapJunctions"][0, 1] \
                            or f["network/gapJunctions"][-1, 1] > neuron_range[1]:
                        self.write_log(f"Bad merge info: {f_name} has post synaptic neuron outside range "
                                       f"{neuron_range[0]}-{neuron_range[1]}")
                        return None

        self.write_log(f"Found merge info in {self.merge_info_file}")

        return merge_files_syn, merge_neuron_range_syn, merge_syn_ctr, \
            merge_files_gj, merge_neuron_range_gj, merge_gj_ctr

    ############################################################################

    def set_scratch_path(self, scratch_path=None):

        """
        Sets scratch path. OBS, need to call open_work_history_file first.

        Args:
            scratch_path (str): Path to scratch directory
        """

        # Check why we required this?
        assert self.work_history_file is not None and self.work_history_file != "last", \
            "Need to call open_work_history_file before set_scratch_path"

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

        """ Destructor, cleans up temporary files after. """

        try:
            if self.hist_file is not None:
                self.hist_file.close()
        except Exception as e:
            print("Hist file already closed?")

        try:
            if self.out_file is not None:
                self.out_file.close()
                self.out_file = None
        except Exception as e:
            print("Out file already closed?")

        # self.clean_up_merge_files()  # -- This caused old files to be cleaned up when aborting. Bad for debugging.

        if self.rc:
            # Clean up memory on workers
            from snudda.utils import cleanup
            cleanup(self.rc, "prune")

    ############################################################################

    def open_work_history_file(self, work_history_file=None, config_file=None):

        """
        Opens work history file.

        Args:
            work_history_file (str, optional) : Path to work history file (if you want to override default,
                                                network-detect-worklog.hdf5)
            config_file (str, optional): Path to config file, if you want to override file used for detection.

        """

        if work_history_file is None:
            work_history_file = self.work_history_file

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

        # Also, we need to add any synapses from projections between volumes to the total
        if "nProjectionSynapses" in self.hist_file:
            self.num_projection_synapses = self.hist_file["nProjectionSynapses"][()]
        else:
            self.num_projection_synapses = 0

        self.num_synapses_total = np.sum(self.hist_file["nHypervoxelSynapses"][()]) + self.num_projection_synapses
        self.num_gap_junctions_total = np.sum(self.hist_file["nHypervoxelGapJunctions"][()])

        if config_file is None:
            self.config_file = self.hist_file["meta/configFile"][()]
        else:
            self.config_file = config_file

        self.position_file = self.hist_file["meta/positionFile"][()]

        # This was config data used for detection, might differ from pruning config
        self.detect_config = json.loads(self.hist_file["meta/config"][()], object_pairs_hook=collections.OrderedDict)
        with open(self.config_file, "r") as f:
            self.config = json.load(f, object_pairs_hook=collections.OrderedDict)

        # If connectivity is empty, then there was nothing to do touch detection on
        # But if it is non-empty, then there should be no remaining hyper voxels
        assert len(remaining) == 0 or len(self.config["Connectivity"]) == 0, \
            (f"Detection not done. There are {len(remaining)} hypervoxels "
             f"not completed: {', '.join([str(x) for x in remaining])}")

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

        """
        Sanity check on hyper voxels. Checks size of voxels, hyper voxels, simulation origo,
        hyper voxel origo, config files, position files, SlurmID and axon stump flags match.

        Args:
            hypervoxel_file : HDF5 file object
            hypervoxel_file_name : Name of HDF5 file
            verbose (bool) : Print extra information, default False
        """

        if verbose:
            self.write_log(f"Checking that {hypervoxel_file_name} matches circuit settings")

        check_list = ["voxelSize", "hyperVoxelSize", "simulationOrigo",
                      "configFile", "positionFile", "SlurmID", "axonStumpIDFlag"]

        # Just some sanity checks
        for c in check_list:
            test = self.hist_file[f"meta/{c}"][()] == hypervoxel_file[f"meta/{c}"][()]
            if type(test) == bool:
                assert test, f"Mismatch of {c} in file {hypervoxel_file_name}"
            else:
                assert test.all(), f"Mismatch of {c} in file {hypervoxel_file_name}"

                # Get xyz coordinates of hyper voxel
        xyz = np.where(self.hyper_voxel_id_list == hypervoxel_file["meta/hyperVoxelID"][()])
        xyz = np.array([x[0] for x in xyz])

        # Just do a sanity check that the hypervoxel origo matches stored value
        hypervoxel_origo = self.simulation_origo + self.hyper_voxel_width * xyz
        assert (hypervoxel_origo == hypervoxel_file["meta/hyperVoxelOrigo"][()]).all(), \
            f"Hyper voxel origo mismatch in file {hypervoxel_file_name}"

        ofc = hypervoxel_file["meta/voxelOverflowCounter"][()]

        if ofc > 0:
            self.voxel_overflow_counter += ofc
            self.overflow_files.append(hypervoxel_file_name)
            self.write_log(f"Overflow of {ofc} in {hypervoxel_file_name}", is_error=True)

    ############################################################################

    def open_hyper_voxel(self, hyper_voxel_id, verbose=False, verify=True):

        """
        Helper function, opens hyper voxel.

        Args:
            hyper_voxel_id (int) : ID of hyper voxel to open
            verbose (bool) : Print extra information
            verify (bool) : Verify hypervoxel integrity, default=True

        """

        if verbose:
            self.write_log(f"Reading hypervoxel {hyper_voxel_id}")

        h_file_name = os.path.join(self.network_path, "voxels", f"network-putative-synapse-{hyper_voxel_id}.hdf5")
        h_file = h5py.File(h_file_name)

        # Just make sure the data we open is OK and match the other data
        if verify:
            self.check_hyper_voxel_integrity(h_file, h_file_name, verbose=verbose)

        return h_file

    ############################################################################

    # This checks that all connections included in the pruning, were present
    # in the detection. If some are missing in the detection phase, then they
    # would incorrectly be missing after pruning.

    def check_network_config_integrity(self, config_file):

        """
        This checks that all connections included in the pruning, were present
        in the detection. If some are missing in the detection phase, then they
        would incorrectly be missing after pruning.

        Args:
            config_file (str) : Path to new config file
        """

        detect_config = json.loads(self.hist_file["meta/config"][()], object_pairs_hook=collections.OrderedDict)
        with open(config_file, "r") as f:
            prune_config = json.load(f, object_pairs_hook=collections.OrderedDict)

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

        """
        Parse the connection information in the config file

        Args:
            config_file (str): Path to config file
        """

        if config_file is None:
            config_file = self.hist_file["meta/configFile"][()]

        self.check_network_config_integrity(config_file=config_file)
        with open(config_file, "r") as f:
            self.config = json.load(f, object_pairs_hook=collections.OrderedDict)

        self.population_unit_id = self.hist_file["network/neurons/populationUnitID"][()]

        # Normally we use type names as lookups, but since we will do this
        # many millions of times, we create an temporary typeID number
        self.make_type_numbering()

        orig_connectivity_distributions = json.loads(self.hist_file["meta/connectivityDistributions"][()],
                                                     object_pairs_hook=collections.OrderedDict)

        config_connectivity_distributions = self.config["Connectivity"]

        self.connectivity_distributions = dict([])

        # For the pruning we merge the two into one
        for key in config_connectivity_distributions:
            (pre_type, post_type) = key.split(",")  # split on "$$" if we had looped over orig_connectivity_distribution
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

                self.connectivity_distributions[pre_type_id, post_type_id, synapse_type_id] = (pruning, pruning_other)

    ############################################################################

    # This makes sure all the variables exist, that way prune_synapses_helper
    # does not have to check, but can assume that they will exist

    @staticmethod
    def complete_pruning_info(prune_info):

        """
        This makes sure all the variables exist, that way prune_synapses_helper
        does not have to check, but can assume that they will exist

        Args:
            prune_info (dict) : Dictionary with pruning information to be complemented.

        """

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

        if "cluster" not in prune_info:
            prune_info["cluster"] = False

        return prune_info

    ############################################################################

    def make_type_numbering(self):

        """
        Normally we use type names as lookups in connection_distribution, but since we will do this
        many millions of times, we create an temporary type_id number.
        """

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

    # TODO: Remove save_morphologies flag, no longer valid since we added NeuronPrototypes
    def setup_output_file(self, output_file=None, save_morphologies=False):

        """
        Sets up output file (network-synapses.hdf5), copying over meta data.

        Args:
            output_file (str, optional): Path to output file, default network-synapses.hdf5
        """

        assert save_morphologies == False, "setup_output_file: save_morphologies currently disabled"

        if self.out_file is not None:
            self.write_log(f"Output file already set: {self.out_file.filename}")
            return

        if output_file is None:
            output_file = os.path.join(self.network_path, "network-synapses.hdf5")

        self.write_log(f"Writing to {output_file}")
        out_file = h5py.File(output_file, "w", libver=self.h5libver, driver=self.h5driver)
        out_file.create_dataset("config", data=json.dumps(self.config))

        # Copy over meta data
        self.hist_file.copy("meta", out_file)

        cfg = json.loads(self.hist_file["meta/config"][()], object_pairs_hook=collections.OrderedDict)
        morph_group = out_file.create_group("morphologies")

        for name, definition in cfg["Neurons"].items():
            morph_file = definition["morphology"]

            swc_group = morph_group.create_group(name)
            swc_group.create_dataset("location", data=morph_file)

            # We now allow multiple variations of each morphology, so we no longer save them in the HDF5 file
            #
            # if save_morphologies:
            #     self.write_log(f"Saving morphology in HDF5 file: {morph_file}")
            #     with open(snudda_parse_path(morph_file), "r") as f:
            #         swc_data = f.read()
            #
            #     swc_group.create_dataset("swc", data=swc_data)

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

    # TODO: Remove save_morphologies as option, no longer possible since we are using NeuronPrototypes
    #  with multiple morphologies

    def setup_merge_file(self, big_cache=False,
                         outfile_name=None,
                         save_morphologies=False,
                         num_synapses=None,
                         num_gap_junctions=None,
                         delete_after=True):

        """
        Sets up merge file, which is done after detection, but before synapse pruning.
        Writes network-putative-synapses-MERGED.hdf5

        Args:
            big_cache (bool): Increase HDF5 cache size
            outfile_name (str): Name of output file name (default network-putative-synapses-MERGED.hdf5)
            num_synapses (int): Number of expected synapses
            num_gap_junctions (int): Number of expected gap junctions
            delete_after (bool): Clean up hypervoxel files after

        """

        assert save_morphologies == False, "setup_merge_file: save_morphologies currently disabled"

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

        cfg = json.loads(self.hist_file["meta/config"][()], object_pairs_hook=collections.OrderedDict)

        # Save morphologies
        # -- We no longer allow morphologies to be saved in the HDF5 file, since now we can have
        #    multiple variations of each morphology. Maybe activate it again in the future

        # if save_morphologies:
        #     morph_group = out_file.create_group("morphologies")
        #
        #     for name, definition in cfg["Neurons"].items():
        #         morph_file = definition["morphology"]
        #
        #         with open(snudda_parse_path(morph_file), "r") as f:
        #             swc_data = f.read()
        #
        #         self.write_log(f"Saving morphology in HDF5 file: {morph_file}")
        #         swc_group = morph_group.create_group(name)
        #         swc_group.create_dataset("swc", data=swc_data)
        #         swc_group.create_dataset("location", data=morph_file)

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
                                synapse_ctr,
                                merge_data_type="synapses",
                                close_input_file=True):

        """
        Helper function, performs synapse (or gap junction) pruning in parallel.

        Args:
            synapse_file (list) : List with paths to synapse files
            synapse_ctr (int) : List or array with number of synapses
            merge_data_type (str) : "synapses" or "gapJunctions"
            close_input_file (bool) : Close input file after, or keep open (default True)

        """

        if type(synapse_ctr) == list:
            synapse_ctr = np.array(synapse_ctr)

        if type(synapse_file) != list:
            # Run in serial
            syn_before, syn_after = self.prune_synapses(synapse_file=synapse_file,
                                                        output_filename=None,
                                                        row_range=None,
                                                        close_out_file=False,
                                                        close_input_file=close_input_file,
                                                        merge_data_type=merge_data_type)

            # Just add a sanity check, that all synapses promised before pruning are accounted for
            assert syn_before == synapse_ctr, \
                f"prune_synapse_parallel: serial run, received {syn_before}, expected {synapse_ctr}"

            return

        if self.d_view:
            self.setup_parallel(d_view=self.d_view)

            assert len(synapse_file) == len(self.d_view), \
                (f"Internal mismatch, n_workers={len(self.d_view)}, n_synapse_files={len(synapse_file)}"
                 f" (first prune was run with {len(synapse_file)} workers), these must match when rerunning prune.")

            self.d_view.scatter("synapse_filename", synapse_file, block=True)

            # 1. Pick names for the workers
            temp_output_file_name = [os.path.join(self.scratch_path, f"worker-temp-{merge_data_type}-file-{x}")
                                     for x in range(0, len(self.d_view))]
            self.d_view.scatter("output_filename", temp_output_file_name, block=True)
            self.d_view.push({"merge_data_type": merge_data_type}, block=True)

            cmd_str = ("syn_before, syn_after = sp.prune_synapses(synapse_file=synapse_filename[0],"
                       "                                          output_filename=output_filename[0],"
                       "                                          merge_data_type=merge_data_type)")

            start_time = timeit.default_timer()

            self.d_view.execute(cmd_str, block=True)
            syn_before = np.array(self.d_view.gather("syn_before"))
            syn_after = np.array(self.d_view.gather("syn_after"))

            # Check all synapses prior to pruning existed
            assert (syn_before == synapse_ctr).all(), \
                f"prune_synapse_parallel: parallel run, received {syn_before.sum()}, expected {synapse_ctr.sum()}"

            # Add the files to a delete list, so we remove them after
            for f in temp_output_file_name:
                self.temp_file_list.append(f)

            syn_after_merge = self.combine_files(temp_output_file_name, merge_data_type)

            assert syn_after_merge == syn_after.sum(), \
                f"prune_synapses_parallel: parallel run, gathered {syn_after_merge}, expected {syn_after.sum()}"

            end_time2 = timeit.default_timer()
            self.write_log(f"prune_synapses_parallel "
                           f"({syn_after_merge}/{syn_before.sum()} {merge_data_type}, " 
                           f"{(100*syn_after_merge/np.maximum(1, syn_before.sum())):.1f}% kept)"
                           f": {end_time2 - start_time:.1f}s", force_print=True)

        else:
            # Multiple files but we are running in serial
            self.write_log(f"Warning, multiple_files but running {merge_data_type} in serial")

            syn_before_total = 0
            self.setup_output_file()

            for syn_file in synapse_file:
                syn_before, syn_after = self.prune_synapses(synapse_file=syn_file, output_filename=None,
                                                            merge_data_type=merge_data_type,
                                                            close_input_file=close_input_file, close_out_file=False)
                syn_before_total += syn_before
                assert syn_before == synapse_ctr.sum(), \
                    (f"prune_synapse_parallel: serial run (multi files), "
                     f"received {syn_before_total}, expected {synapse_ctr.sum()}")

            # Need to resize synapse matrix
            n_synapses = self.out_file["network/nSynapses"][0]
            n_gj = self.out_file["network/nGapJunctions"][0]
            self.out_file["network/synapses"].resize((n_synapses, self.out_file["network/synapses"].shape[1]))
            self.out_file["network/gapJunctions"].resize((n_gj, self.out_file["network/gapJunctions"].shape[1]))

    ############################################################################

    def save_putative_synapses(self):

        assert self.keep_files, "keep_files must be True to use this feature"

        putative_file_name = os.path.join(self.network_path, "network-putative-synapses.hdf5")

        self.write_log(f"Saving putative synapses to {putative_file_name}")

        old_out_file = self.out_file
        self.out_file = None

        # TODO: Do sanity check to make sure no files are duplicated etc.
        synapse_files = glob.glob(os.path.join(self.network_path, "temp", "synapses*MERGE-ME.hdf5"))
        gj_files = glob.glob(os.path.join(self.network_path, "temp", "synapses*MERGE-ME.hdf5"))

        self.combine_files(source_filenames=synapse_files, merge_data_type="synapses",
                           output_filename=putative_file_name)
        self.combine_files(source_filenames=gj_files, merge_data_type="gapJunctions",
                           output_filename=putative_file_name)

        # Restore old out_file
        self.out_file = old_out_file

    ############################################################################

    def combine_files(self, source_filenames, merge_data_type, output_filename=None):

        """
        Combines synapse files after pruning in network-synapses.hdf5.

        Args:
              source_filenames (list) : List with path to source files
              merge_data_type (str) : "synapses" or "gapJunctions"
              output_filename (str) : Output file name
        """

        start_time = timeit.default_timer()

        if not self.out_file:
            self.setup_output_file(output_file=output_filename)

        h5_syn_mat, h5_hyp_syn_n, h5_syn_n, h5_syn_loc = self.data_loc[merge_data_type]

        tmp_files = [h5py.File(f, 'r') for f in source_filenames if os.path.isfile(f)]
        num_syn = np.sum(np.fromiter(iter=(f[h5_syn_mat].shape[0] for f in tmp_files), dtype=int))
        mat_width_all = [f[h5_syn_mat].shape[1] for f in tmp_files]

        if mat_width_all:
            assert (np.array(mat_width_all) == mat_width_all[0]).all(), \
                "combine_files: Internal error, width does not match"
            mat_width = mat_width_all[0]
        else:
            mat_width = 0

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
        self.write_log(f"combine_files ({merge_data_type}): {end_time2 - start_time:.1f}s")

        # Number of synapses in total
        return next_syn

    ############################################################################

    # Find which ranges of the synapse matrix that each worker should take care of

    def find_ranges(self, synapses, num_workers, start_pos=0, num_syn=None):

        """
        Find which ranges of the synapse matrix that each worker should take care of

        Args:
            synapses : Synapse matrix
            num_workers (int): Number of workers
            start_pos (int) : Starting from position
            num_syn (int): Number of synapses

        """

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

        """ Cleans up merge files. """

        if self.role == "master" and self.d_view:
            # Make workers clean up their files also
            cmd_str = "sp.clean_up_merge_files()"
            self.d_view.execute(cmd_str, block=True)

        if self.temp_file_list is None:
            # Nothing to do
            return

        self.write_log("Cleaning up old merge files")

        for f in self.temp_file_list:
            # self.writeLog("Removing old merge file: " + str(f))
            try:
                if os.path.exists(f):
                    # TODO: Close file, Zahra fix
                    os.remove(f)
            except Exception as e:
                self.write_log(f"Closing of file {f} failed: {e}")

        self.temp_file_list = None

    ############################################################################

    # TODO: Switch to logger instead.

    def write_log(self, text, flush=True, is_error=False, force_print=False):  # Change flush to False in future, debug

        """
        Writes to log file. Use setup_log first. Text is only written to screen if self.verbose=True,
        or is_error = True, or force_print = True.

        test (str) : Text to write
        flush (bool) : Should all writes be flushed to disk directly?
        is_error (bool) : Is this an error, always written.
        force_print (bool) : Force printing, even if self.verbose=False.
        """

        try:
            if self.logfile is not None:
                self.logfile.write(f"{text}\n")
                if flush:
                    self.logfile.flush()

            if self.verbose or is_error or force_print:
                print(text, flush=True)
        except Exception as e:
            print(text)
            print("Unable to write to log file. Is log file closed?")

    ############################################################################

    def setup_parallel(self, d_view):

        """
        Setup workers for parallel execution.

        Args:
            d_view : ipyparallel direct view object
        """

        assert self.role == "master", "setup_parallel: Should only be called by master node"

        if d_view is None:
            self.write_log("setup_parallel called without dView, aborting.")
            return

        if self.workers_initialised:
            self.write_log("Workers already initialised.")
            return

        # This imports for the workers
        with d_view.sync_imports():
            from snudda.detect.prune import SnuddaPrune

        self.write_log(f"Setting up workers: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")

        # Create unique log file names for the workers
        if self.logfile_name is not None:
            engine_log_file = [f"{self.logfile_name}-{x}" for x in range(0, len(d_view))]
        else:
            engine_log_file = [[] for x in range(0, len(d_view))]

        d_view.scatter('logfile_name', engine_log_file, block=True)
        d_view.push({"network_path": self.network_path,
                     "random_seed": self.random_seed,
                     "config_file": self.config_file}, block=True)

        cmd_str = ("sp = SnuddaPrune(network_path=network_path, logfile_name=logfile_name[0]," 
                   "                 config_file=config_file,"
                   "                 role='worker',random_seed=random_seed)")
        d_view.execute(cmd_str, block=True)

        self.write_log(f"Workers setup: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
        self.workers_initialised = True

    ########################################################################

    def clean_up_merge_read_buffers(self):

        """ Housekeeping, clean up merge read buffers. """

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

        """
        Writes buffered synapses to file.

        Args:
            synapse_matrix_loc (str): location of synapses in HDF5 file
            synapses : Synapse to write
            flush (bool): Buffer needs to be flushed at last call, otherwise all is not written to file

        """

        if self.synapse_write_buffer is None:
            # First time run
            if synapses is not None:
                buffer_width = synapses.shape[1]
            else:
                assert False, "buffer_merge_write: Unable to guess width of matrix"

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

    # This goes through the hyper voxel synapse files, extracts a range of neurons,
    # and puts their synapses in its own file. Parallel execution, all neurons.

    def gather_neuron_synapses(self):

        """
        Goes through the hyper voxel synapse files, extracts a range of neurons,
        and puts their synapses in its own file. Parallel execution, all neurons.
        """

        if self.role != "master":
            self.write_log("gather_neuron_synapses is only run on master node, aborting")
            return

        start_time = timeit.default_timer()

        # Split neurons between nodes, we need the neurons to be in order
        num_neurons = self.hist_file["network/neurons/neuronID"].shape[0]
        assert num_neurons - 1 == self.hist_file["network/neurons/neuronID"][-1], \
            "neuronID should start from 0 and the end should be n-1"

        if not self.d_view:

            # Run in serial, save as a list to make result compatible with parallel version of code
            merge_results_syn = [self.big_merge_helper(neuron_range=np.array([0, num_neurons]),
                                                       merge_data_type='synapses')]
            merge_results_gj = [self.big_merge_helper(neuron_range=np.array([0, num_neurons]),
                                                      merge_data_type='gapJunctions')]

        else:
            self.setup_parallel(d_view=self.d_view)

            n_workers = len(self.d_view)
            neuron_ranges = []
            range_borders = np.linspace(0, num_neurons, n_workers + 1).astype(int)

            for idx in range(0, n_workers):
                neuron_ranges.append((range_borders[idx], range_borders[idx + 1]))

            assert neuron_ranges[-1][-1] == num_neurons, \
                "gather_neuron_synapses: Problem with neuron_ranges, last element incorrect"
            assert len(neuron_ranges) == n_workers, \
                "gather_neuron_synapses: Problem with neuron_ranges, bad length"

            # Send list of neurons to workers
            self.d_view.scatter("neuron_range", neuron_ranges, block=True)

            # Each worker sorts a subset of the neurons and write it to separate files
            cmd_str_syn = ("merge_result_syn = sp.big_merge_helper(neuron_range=neuron_range[0], "
                           "merge_data_type='synapses')")

            self.d_view.execute(cmd_str_syn, block=True)
            merge_results_syn = self.d_view["merge_result_syn"]

            # When we do scatter, it embeds the result in a list
            cmd_str_gj = ("merge_result_gj = sp.big_merge_helper(neuron_range=neuron_range[0], "
                          "merge_data_type='gapJunctions')")
            self.d_view.execute(cmd_str_gj, block=True)
            merge_results_gj = self.d_view["merge_result_gj"]

        # Sort the files in order
        merge_start_syn = [x[1][0] for x in merge_results_syn]
        merge_start_gj = [x[1][0] for x in merge_results_gj]

        merge_order_syn = np.argsort(merge_start_syn)
        merge_order_gj = np.argsort(merge_start_gj)

        # Extract data, and return it in more meaningful format
        merge_files_syn = [merge_results_syn[idx][0] for idx in merge_order_syn]
        merge_files_gj = [merge_results_gj[idx][0] for idx in merge_order_gj]

        merge_neuron_range_syn = [merge_results_syn[idx][1] for idx in merge_order_syn]
        merge_neuron_range_gj = [merge_results_gj[idx][1] for idx in merge_order_gj]

        merge_syn_ctr = [merge_results_syn[idx][2] for idx in merge_order_syn]
        merge_gj_ctr = [merge_results_gj[idx][2] for idx in merge_order_gj]

        end_time = timeit.default_timer()
        self.write_log(f"gather_neuron_synapses took {end_time - start_time:.1f} s")

        return merge_files_syn, merge_neuron_range_syn, merge_syn_ctr, \
            merge_files_gj, merge_neuron_range_gj, merge_gj_ctr

    ############################################################################

    #    We need to find all the hypervoxels that might contain synapses.
    #    During touch detection we generated a list for each hypervoxel with
    #    all the neurons that had any part of it within that hypervoxel.

    def get_hyper_voxel_list(self, neuron_range):

        """
        Gather a list of all hyper voxels which contains synapses belonging to neurons in neuron_range.

        Args:
            neuron_range : Range of neurons

        """

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

        """
        Gathers synapses (or gap junctions) belonging to neuron_range from multiple hyper voxels into one merge file.

        Args:
            neuron_range : Range of neurons
            merge_data_type : "synapses" or "gapJunctions"

        """

        try:
            self.write_log(f"big_merge_helper ({merge_data_type}): neuron_range = {neuron_range}")

            output_filename = os.path.join(self.scratch_path,
                                           f"{merge_data_type}-for-neurons-"
                                           f"{neuron_range[0]}-to-{neuron_range[1]}-MERGE-ME.hdf5")

            self.merge_data_type = merge_data_type

            # Which hyper voxels are the neurons located in?
            hv_list = self.get_hyper_voxel_list(neuron_range)

            # Locations of data within the file
            h5_syn_mat, h5_hyp_syn_n, h5_syn_n, h5_syn_lookup = self.data_loc[merge_data_type]

            synapse_heap = []
            file_list = dict([])
            file_mat_iterator = dict([])
            chunk_size = 10000

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
                                               self.hist_file[h5_hyp_syn_n][:n_hv],
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
                            self.write_log("Investigate", is_error=True)
                            self.write_log(f"{self.max_channel_type} != "
                                           f"{file_list[h_id]['network/maxChannelTypeID'][()]}", is_error=True)
                            # import pdb
                            # pdb.set_trace()

                        # These should be the same for all hypervoxels
                        assert self.max_channel_type == file_list[h_id]["network/maxChannelTypeID"][()], \
                            (f"max_channel_type = {self.max_channel_type} "
                             f"(differ with what is in file {file_list[h_id]['network/maxChannelTypeID'][()]})")
                    else:
                        self.max_channel_type = file_list[h_id]["network/maxChannelTypeID"][()]
                        self.write_log(f"Setting max_channel_type to {self.max_channel_type} from h_id={h_id}")

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

            # --- Start of special code for projection synapses from project.py

            if merge_data_type == "synapses" and self.num_projection_synapses > 0:

                assert os.path.exists(self.projection_synapse_file), \
                  f"Missing projection connection file: {self.projection_synapse_file}"

                self.write_log(f"Adding projection synapses from {self.projection_synapse_file}")
                # There is also a file with synapse connections, add it to the merge set
                # Since they were not created using hyper voxels, they get the special h_id = -1
                proj_connection = -1

                fid = h5py.h5f.open(self.projection_synapse_file.encode(),
                                    flags=h5py.h5f.ACC_RDONLY, fapl=propfaid)
                file_list[proj_connection] = h5py.File(fid, drive=self.h5driver)

                if self.max_channel_type:
                    assert self.max_channel_type == file_list[proj_connection]["network/maxChannelTypeID"][()], \
                        "max_channel_type does not match for projection file"
                else:
                    self.max_channel_type = file_list[proj_connection]["network/maxChannelTypeID"][()]

                assert file_list[proj_connection]["network/nSynapses"][()] == self.num_projection_synapses, \
                    (f"Mismatch between work history file and data file. "
                     f"nProjectionSynapses: {self.num_projection_synapses} vs {self.num_projection_synapses}")

                if file_list[proj_connection]["network/nSynapses"][()] > 0:
                    n_total += file_list[proj_connection]["network/nSynapses"][()]

                    lookup_iterator = \
                        self.file_row_lookup_iterator_subset(h5mat_lookup=file_list[proj_connection][h5_syn_lookup],
                                                             min_dest_id=neuron_range[0],
                                                             max_dest_id=neuron_range[1],
                                                             chunk_size=chunk_size)

                    file_mat_iterator[proj_connection] \
                        = self.synapse_set_iterator(h5mat_lookup=file_list[proj_connection][h5_syn_lookup],
                                                    h5mat=file_list[proj_connection][h5_syn_mat],
                                                    chunk_size=chunk_size,
                                                    lookup_iterator=lookup_iterator)

                    syn_set, unique_id = next(file_mat_iterator[proj_connection], (None, None))

                    if syn_set is None:
                        # No synapse in our range, let this worker skip the file
                        # Clear file List and fileMatIterator for this worker
                        del file_list[proj_connection]
                        del file_mat_iterator[proj_connection]
                    else:
                        # Create a heap containing the first subset of all files
                        heapq.heappush(synapse_heap, (unique_id, proj_connection, syn_set))

                else:
                    # No synapses in file, close it.
                    file_list[proj_connection].close()
                    del file_list[proj_connection]

            # --- End of code for projection synapses

            if merge_data_type == "synapses":
                num_synapses = n_total
                num_gap_junctions = 0
            elif merge_data_type == "gapJunctions":
                num_synapses = 0
                num_gap_junctions = n_total
            else:
                assert False, f"Unknown mergeDataType {merge_data_type}"

            # Setup output file
            (self.buffer_out_file, out_filename) = self.setup_merge_file(big_cache=True, outfile_name=output_filename,
                                                                         save_morphologies=False,
                                                                         num_synapses=num_synapses,
                                                                         num_gap_junctions=num_gap_junctions)

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

                if loop_ctr % 1000000 == 0 and n_total > 1000000:
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

            if n_total > 1000000:
                self.write_log(f"Worker {merge_data_type}: {syn_ctr}/{n_total} (heap size: {len(synapse_heap)})",
                               force_print=True)

            self.write_log(f"Read {syn_ctr} out of total {n_total} {merge_data_type}", force_print=True)

            self.buffer_merge_write(h5_syn_mat, flush=True)
            self.write_log("big_merge_helper: done")

            # Close the hyper voxel files
            for f in file_list:
                try:
                    file_list[f].close()
                except:
                    import traceback
                    t_str = traceback.format_exc()
                    self.write_log(t_str)
                    self.write_log(f"Problems closing files {f}, {file_list[f]}")

            self.buffer_out_file.close()

            return output_filename, neuron_range, syn_ctr
        except:
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str, is_error=True)
            sys.exit(-1)

    ############################################################################

    def prune_synapses(self, synapse_file, output_filename,
                       merge_data_type, row_range=None,
                       close_input_file=True,
                       close_out_file=True):

        """
        Prune synapses.

        Args:
            synapse_file : Path to file with putative synapses
            output_filename : Path to file to write pruned synapses to
            merge_data_type : "synapses" or "gapJunctions"
            row_range : Synapse row range to prune
            close_input_file (bool) : Close input files after
            close_out_file (bool) : Close output file after
        """

        try:
            h5_syn_mat, h5_hyp_syn_n, h5_syn_n, h5_syn_loc = self.data_loc[merge_data_type]

            if synapse_file is None:
                self.write_log(f"prune_synapses: No synapse_file specified for {merge_data_type} -- none detected?")
                return 0, 0

            if type(synapse_file) == str:
                self.write_log(f"Opening synapse file: {synapse_file}")
                synapse_file = h5py.File(synapse_file, 'r')
                # if self.role != "master":
                #     # SWMR = one writer, multiple readers
                #     synapse_file = h5py.File(synapse_file, 'r', swmr=True)
                # else:
                #    synapse_file = h5py.File(synapse_file, 'r')

            if row_range is None:
                row_start = 0
                row_end = synapse_file[h5_syn_mat].shape[0]
            else:
                row_start = row_range[0]
                row_end = row_range[-1]

            if row_start is None or row_end is None or row_start == row_end:
                self.write_log("prune_synapses: Nothing to do, empty row range")
                return 0, 0

            if synapse_file[h5_syn_mat].shape[0] == 0:
                self.write_log(f"prune_synapses: No {merge_data_type} skipping pruning")
                return 0, 0

            self.write_log(f"prune_synapses: synapseFile={synapse_file}, outputFileName={output_filename}"
                           f", rowRange={row_range} ({merge_data_type})")

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

            num_syn_kept = 0

            for synRange in block_ranges:
                self.write_log(f"Pruning range: {synRange}")

                synapses = synapse_file[h5_syn_mat][synRange[0]:synRange[-1]]
                num_syn_kept += self.prune_synapses_helper(synapses=synapses, output_file=self.out_file,
                                                           merge_data_type=merge_data_type)

            # Close synapse input file
            if close_input_file:
                synapse_file.close()

            if close_out_file:
                self.out_file.close()
                self.out_file = None

        except Exception as e:
            # Write error to log file to help trace it.
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str, is_error=True)

            sys.exit(-1)

        return num_syn, num_syn_kept

    ############################################################################

    def prune_synapses_helper(self, synapses, output_file, merge_data_type):

        """
        Helper function to prunes synapses. It takes a subset of the synapse matrix as input, as it needs to keep
        all synapses onto the same neuron in memory at the same time. It buffers writes until the end.

        Args:
            synapses: subset of synapse matrix that fits in memory
            output_file: where to write synapses, assumed to already exist
            merge_data_type : "synapses" or "gapJunctions"

        """
        # Tried to use numba, but dictionaries and h5py._hl.files.File not supported

        h5_syn_mat, h5_hyp_syn_n, h5_syn_n, h5_syn_loc = self.data_loc[merge_data_type]

        keep_row_flag = np.zeros((synapses.shape[0],), dtype=bool)

        next_read_pos = 0
        read_end_of_range = synapses.shape[0]

        # Random seeds for reproducability
        neuron_seeds = self.get_neuron_random_seeds()
        previous_post_synaptic_neuron_id = None
        post_rng = None

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

                if self.all_neuron_pair_synapses_share_parameter_id:
                    # Make all synapses between a particular pair of neurons have the same parameter ID
                    synapses[next_read_pos:read_end_idx, 12] = synapses[next_read_pos, 12]

            # Temp check
            assert ((synapses[next_read_pos:read_end_idx, 0] == synapses[next_read_pos, 0]).all()
                    and (synapses[next_read_pos:read_end_idx, 1] == synapses[next_read_pos, 1]).all()), \
                "prune_synapses_helper: Internal error, more than one neuron pair"

            n_pair_synapses = read_end_idx - next_read_pos

            src_id = synapses[next_read_pos, 0]
            dest_id = synapses[next_read_pos, 1]

            if dest_id != previous_post_synaptic_neuron_id:
                # New post synaptic cell, reseed random generator
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

                # These will always exist thanks to complete_pruning_info function

                dist_p = c_info["distPruning"]  # Dist dep pruning
                f1 = c_info["f1"]
                soft_max = c_info["softMax"]
                mu2 = c_info["mu2"]
                a3 = c_info["a3"]

                # If cluster_flag is set, then the synapses furthest from their companion synapses are removed first
                cluster_flag = c_info["cluster"]

            else:
                # Not listed in connectivityDistribution, skip neuron pair
                next_read_pos = read_end_idx
                # No need to update keepRowFlag since default set to 0

                continue

            # Lets get all cell pairs random numbers in one go, total: n_pair_synapses*3 + 2
            # f1        :n_pair_synapses
            # p         n_pair_synapses:2*n_pair_synapses
            # soft_max  2*n_pair_synapses:3*n_pair_synapses
            # a3       -1
            # p_mu     -2

            random_pool = post_rng.random(n_pair_synapses*3 + 2)

            # 3. This is the last step of pruning, but we move it to the top
            # since there is no point doing the other steps if we going to
            # throw them away anyway

            if a3 is not None and random_pool[-1] > a3:
                # Prune all synapses between pair, do not add to synapse file
                next_read_pos = read_end_idx
                # No need to update keepRowFlag since default set to 0

                continue

            if dist_p is not None:
                assert synapse_type != 3, \
                    "Distance dependent pruning currently only supported for synapses, not gap junctions"
                # Distance dependent pruning, used for e.g. FS->MS connections

                # distP contains d (variable for distance to soma)
                d = synapses[next_read_pos:read_end_idx, 8] * 1e-6  # dendrite distance d, used in eval below
                p = numexpr.evaluate(dist_p)

                frac_flag = random_pool[:n_pair_synapses] < f1
                dist_flag = random_pool[n_pair_synapses:2*n_pair_synapses] < p

                keep_row_flag[next_read_pos:read_end_idx] = np.logical_and(frac_flag, dist_flag)

            else:
                keep_row_flag[next_read_pos:read_end_idx] = random_pool[:n_pair_synapses] < f1
                dist_flag = None

            # Check if too many synapses, trim it down a bit
            n_keep = np.sum(keep_row_flag[next_read_pos:read_end_idx])

            if soft_max is not None and n_keep > soft_max:
                soft_max = float(soft_max)
                p_keep = np.divide(2 * soft_max, (1 + np.exp(-(n_keep - soft_max) / 5)) * n_keep)

                keep_row_flag[next_read_pos:read_end_idx] = \
                    np.logical_and(p_keep > random_pool[2*n_pair_synapses:3*n_pair_synapses],
                                   keep_row_flag[next_read_pos:read_end_idx])

                n_keep = np.sum(keep_row_flag[next_read_pos:read_end_idx])

            # If too few synapses, remove all synapses
            if mu2 is not None:
                # Markram et al, Cell 2015
                p_mu = 1.0 / (1.0 + np.exp(-8.0 / mu2 * (n_keep - mu2)))

                if p_mu < random_pool[-2]:

                    # Too few synapses, remove all -- need to update keepRowFlag
                    keep_row_flag[next_read_pos:read_end_idx] = 0
                    next_read_pos = read_end_idx

                    continue

            # This code remaps which synapses are kept, such that synapses in a cluster are more likely to be kept
            if cluster_flag and n_keep > 0:
                # The rows that passed distance dependent pruning are: dist_flag

                # 1. Calculate distance between all synapses, smallest total distance (sum to all neighbours) kept
                synapse_coords = synapses[next_read_pos:read_end_idx, 2:5]

                if dist_flag is not None:
                    # If dist_flag is set, we need to pick a subset from the ones that passed distance dependent pruning
                    synapse_coords = synapse_coords[dist_flag, :]
                    lookup_idx = np.where(dist_flag)[0]
                else:
                    lookup_idx = None

                # pdist faster, but does not give full distance matrix
                synapse_dist = scipy.spatial.distance.cdist(synapse_coords, synapse_coords)
                synapse_tot_dist = np.sum(synapse_dist, axis=0)
                synapse_priority = np.argsort(synapse_tot_dist)

                keep_idx = synapse_priority[:n_keep]
                if dist_flag is not None:
                    keep_idx = lookup_idx[keep_idx]

                keep_row_flag[next_read_pos:read_end_idx] = 0
                keep_row_flag[next_read_pos+keep_idx] = 1

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

        return n_keep_tot

    ############################################################################

    @staticmethod
    @jit(nopython=True)
    def file_row_lookup_iterator(h5mat, chunk_size=10000):

        """
        File row lookup, to quickly find range of synapses onto a neuron.

        Args:
            h5mat : Synapse matrix
            chunk_size : Chunk size to process each time
        """

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

        """
        File row lookup iterator, working on a subset of the synapse matrix.
        min_dest_id and max_dest_id are inclusive, only synapses with dest_id in that range are iterated over

        Args:
            h5mat_lookup : Synapse lookup
            min_dest_id : Minimum neuron destination ID for synapses  (inclusive)
            max_dest_id : Maximum neuron destination ID for synapses (inclusive)
            chunk_size : Chunk size for processing
        """

        num_neurons = self.hist_file["network/neurons/neuronID"].shape[0]

        assert self.max_channel_type is not None, "max_channel_type should not be None"

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

    @staticmethod
    def synapse_set_iterator(h5mat_lookup, h5mat, chunk_size=10000, lookup_iterator=None):

        """
        Iterator over synapses.

        Args:
            h5mat_lookup : Synapse lookup table
            h5mat : Synapse matrix
            chunk_size : Chunk size to process
            lookup_iterator : Synapse lookup iterator
        """

        # Allow the user to set an alternative lookupIterator if we only
        # want to iterate over a subset of the synapses
        if not lookup_iterator:
            lookup_iterator = SnuddaPrune.file_row_lookup_iterator(h5mat_lookup, chunk_size=chunk_size)

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
                    "We missed one synapse! (2)"

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
                       "We missed one synapse! (1)"

                assert (syn_mat[:, 0] == syn_mat[0, 0]).all() and (syn_mat[:, 1] == syn_mat[0, 1]).all(), \
                       "Synapse matrix (1) contains more than one pair:\n{syn_mat}"

                assert syn_mat.shape[0] == end_idx - start_idx, \
                    "Synapse matrix has wrong size"

                yield (syn_mat,
                       unique_id)

            next_row_set = next(lookup_iterator, None)

    def get_neuron_random_seeds(self):

        """ Derive neuron random seeds from pruning master seed. """

        num_neurons = self.hist_file["network/neurons/neuronID"].shape[0]

        if self.cache_neuron_seeds is not None \
            and self.cache_old_seed == self.random_seed \
                and self.cache_num_neurons == num_neurons:

            neuron_seeds = self.cache_neuron_seeds
        else:
            assert num_neurons - 1 == self.hist_file["network/neurons/neuronID"][-1], \
                "neuronID should start from 0 and the end should be n-1"

            # Need different seeds for each post synaptic neuron
            ss = np.random.SeedSequence(self.random_seed)
            neuron_seeds = ss.generate_state(num_neurons)

            # Cache results for next iteration
            self.cache_neuron_seeds = neuron_seeds
            self.cache_old_seed = self.random_seed
            self.cache_num_neurons = num_neurons

        return neuron_seeds

    # This removes the files in the temp and voxels directory, freeing up diskspace
    def cleanup(self):

        """ Removes files in temp and voxels directories, freeing up diskspace. """

        temp_path = os.path.join(self.network_path, "temp")
        voxel_path = os.path.join(self.network_path, "voxels")

        for path in [temp_path, voxel_path]:
            self.write_log(f"Removing temp files from {path}")
            files = glob.glob(os.path.join(path, "*"))
            for f in files:
                if os.path.isfile(f):
                    try:
                        os.remove(f)
                    except:
                        self.write_log(f"cleanup: Failed to close {f}, is it already closed?")

##############################################################################


if __name__ == "__main__":
    print("Please do not call this file directly, use snudda.py")
    sys.exit(-1)
