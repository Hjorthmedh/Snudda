# snudda_detect.py
#
# Johannes Hjorth, Royal Institute of Technology (KTH)
# Human Brain Project 2019
#
# Requires Python version 3+
#
# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union’s Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).
#
import functools
import itertools
import json
import os
import sys
import time
import timeit
import gc

import h5py
import numexpr
import numpy as np
from numba import jit
import copy

import snudda.utils.memory
from snudda.neurons import NeuronMorphologyExtended
from snudda.neurons.morphology_data import MorphologyData
from snudda.utils import NumpyEncoder
from snudda.utils.snudda_path import get_snudda_data, snudda_parse_path
from snudda.detect.projection_detection import ProjectionDetection
from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.utils.load import SnuddaLoad

# from memory_profiler import profile
# Put @profile decorator over function: https://pypi.org/project/memory-profiler/

# TODO: Exclude neurons without synapses or gap junctions from touch detection (ie if no pre/post connections possible)


class SnuddaDetect(object):
    """
    SnuddaDetect places synapses in the network based on touch detection.
    """

    def __init__(self,
                 config_file=None,
                 network_path=None,
                 snudda_data=None,
                 position_file=None,
                 voxel_size=3e-6,  # 2e-6,
                 hyper_voxel_size=100,  # 250, #100,
                 verbose=False,
                 logfile_name=None,
                 logfile=None,
                 save_file=None,
                 work_history_file=None,
                 slurm_id=0,
                 volume_id=None,
                 role=None,  # Default: "master"
                 rc=None,
                 simulation_origo=None,  # Auto detect
                 h5libver=None,  # Default: "latest"
                 random_seed=None,
                 debug_flag=False):

        """
        Constructor.

        Args:
            network_path (str): Network directory
            snudda_data (str, optional): Path to SNUDDA_DATA (if not specified will try to read from config file)
            config_file (str, optional): Network config file (default network-config.json in network_path)
            position_file (str, optional): Network position file (default network-neuron-positions in network_path)
            voxel_size (float, optional): Width of voxel (default 3e-6m)
            hyper_voxel_size (int, optional): Number of voxels per side (default: 100, ie 100x100x100 voxels total)
            verbose (bool, optional): Verbose mode (default False)
            logfile_name (str, optional): Name of log file
            logfile (_io.TextIOWrapper, optional): Pointer to already open log file
            save_file (str, optional): Name of output file (default voxels/network-putative-synapses.hdf5 in network_path)
            work_history_file (str, optional): Work log file (default network-detect-worklog.hdf5 in network_path)
            slurm_id (int, optional): slurm_id of job
            volume_id (str, optional): Volume ID to do touch detection on
            role (str, optional): Parallel role, i.e. "master" or "worker"
            rc (ipyparallel.Client, optional): iPyParallel client, if given program will run in parallel
            simulation_origo (np.array, optional): Origo for touch detection hypervoxels and voxels, voxel coordinates must always positive.
            h5libver (string, optional): h5py library version (default "latest")
            random_seed (int, optional): Random seed
            debug_flag (bool, optional): Save additional information for debugging (Default: False)

        """

        self.rc = rc

        if not role:
            self.role = "master"
        else:
            self.role = role

        assert self.role in ["master", "worker"], \
            "SnuddaDetect: Role must be master or worker"

        self.verbose = verbose

        if not h5libver:
            self.h5libver = "latest"
        else:
            self.h5libver = h5libver

        self.debug_flag = debug_flag

        self.random_seed = random_seed

        if config_file and not network_path:
            network_path = os.path.dirname(config_file)

        if network_path:
            self.network_path = network_path

            # Setting default config, position, save and log file if not user provided
            if not config_file:
                config_file = os.path.join(network_path, "network-config.json")

            if not position_file:
                position_file = os.path.join(network_path, "network-neuron-positions.hdf5")

            if not save_file:
                save_file = os.path.join(network_path, "voxels", "network-putative-synapses.hdf5")

            if not logfile and not logfile_name:
                log_filename = os.path.join(network_path, "log", "touch-detection.txt")
        else:
            self.network_path = None

        self.work_history_file = work_history_file  # Name of work history file
        self.work_history = None  # File pointer for actual file

        if logfile_name:
            self.logfile_name = logfile_name
        elif logfile is not None:
            self.logfile_name = logfile.name
        elif self.network_path is not None:
            self.logfile_name = os.path.join(self.network_path, "log", "touch-detection.txt")
        else:
            self.logfile_name = None

        self.logfile = logfile
        self.setup_log()

        self.config_file = config_file
        self.position_file = position_file
        self.save_file = save_file

        self.snudda_data = get_snudda_data(snudda_data=snudda_data,
                                           config_file=self.config_file,
                                           network_path=self.network_path)

        self.write_log(f"Using hdf5 driver version: {self.h5libver}")

        mem = self.memory()
        self.write_log(f"{mem}")

        self.slurm_id = int(slurm_id)  # Make sure integer
        self.workers_initialised = False

        self.voxel_size = voxel_size
        self.hyper_voxel_size = hyper_voxel_size  # = N,  N x N x N voxels in a hyper voxel
        self.hyper_voxel_origo = np.zeros((3,))
        self.voxel_overflow_counter = 0
        self.step_multiplier = 1.1  # 1.12  # 2.0

        self.hyper_voxel_offset = None
        self.hyper_voxel_id = 0
        self.hyper_voxel_rng = None

        self.num_bins = hyper_voxel_size * np.ones((3,), dtype=int)
        self.write_log("Each hyper voxel has %d x %d x %d voxels" % tuple(self.num_bins))

        # These are voxels that axons/dend occupy
        self.axon_voxels = None
        self.dend_voxels = None

        self.volume_id = volume_id
        if volume_id is not None:
            self.write_log(f"Touch detection only {volume_id}")
        else:
            self.write_log("Touch detecting all volumes")

        # These are counters, indicating how many (different) axons/dend in the voxel
        self.axon_voxel_ctr = None
        self.dend_voxel_ctr = None

        self.max_axon_voxel_ctr = None
        self.max_dend_voxel_ctr = None

        self.axon_soma_dist = None
        self.dend_sec_id = None
        self.dend_sec_x = None
        self.dend_soma_dist = None

        self.neurons = None
        self.neuron_positions = None
        self.population_unit = None

        self.hyper_voxels = None
        self.hyper_voxel_id_lookup = None
        self.num_hyper_voxels = None
        self.hyper_voxel_width = self.hyper_voxel_size * self.voxel_size
        self.simulation_origo = np.array(simulation_origo) if simulation_origo is not None else None

        self.config = None

        # Columns in hyperVoxelSynapses:
        # 0: source_cell_id, 1: dest_cell_id, 2: voxel_x, 3: voxel_y, 4: voxel_z,
        # 5: hyper_voxel_id, 6: channel_model_id,
        # 7: source_axon_soma_dist (not SI scaled 1e6, micrometers),
        # 8: dest_dend_soma_dist (not SI scaled 1e6, micrometers)
        # 9: dest_sec_id, 10: dest_sec_x (int 0 - 1000, SONATA wants float 0.0-1.0)
        # 11: conductance (int, not SI scaled 1e12, in pS)
        # 12: parameter_id
        #
        # Note on parameter_id:
        # If there are n parameter sets for the particular synapse type, then
        # the ID to use is parameterID % n, this way we can reuse connectivity
        # if we add more synapse parameter sets later.

        self.hyper_voxel_synapses = None

        # Columns in hyper_voxel_gap_junctions
        # 0: source_cell_id, 1: dest_cell_id, 2: source_sec_id, 3: dest_sec_id,
        # 4: source_seg_x, 5: dest_seg_x, 6: voxel_x, 7: voxel_y, 8: voxel_z,
        # 9: hyper_voxel_id, 10: conductance (integer, in pS)
        self.hyper_voxel_gap_junctions = None

        self.hyper_voxel_synapse_ctr = 0
        self.hyper_voxel_gap_junction_ctr = 0

        self.hyper_voxel_coords = dict([])

        # This is used by the heap sort, when merging hdf5 files for the different
        # hyper voxels
        self.hyper_voxel_synapse_lookup = None
        self.hyper_voxel_gap_junction_lookup = None

        # Parameters for the HDF5 writing, this affects write speed
        self.synapse_chunk_size = 10000
        self.gap_junction_chunk_size = 10000
        self.h5compression = "lzf"

        # This is an upper limit how many axon/dend we allow in each voxel max
        # 10 overflowed
        self.max_axon = 45
        self.max_dend = 20

        self.max_neurons = 10000
        self.max_synapses = 2000000
        self.max_gap_junctions = 100000

        self.connectivity_distributions = dict([])
        # self.connectivityDistributionsGJ = dict([])
        self.next_channel_model_id = 10

        self.prototype_neurons = dict([])
        self.neuron_cache = dict()
        self.extra_axon_cache = dict()

        self.axon_cum_density_cache = dict([])

        self.delete_old_merge()

        # Rather than load all neuron morphologies, we only load prototypes
        self.read_prototypes(config_file=config_file)

        # Read positions
        self.read_neuron_positions(position_file)

        self.projection_detection = None  # Helper class for handling projections between structures

    def detect(self, restart_detection_flag=True, rc=None):

        """
        Synapse placement based on touch detection. Space is divided into hyper voxels, containing 100x100x100 voxels.
        Each hyper voxel is processed separately.

        Args:
              restart_detection_flag (bool, optional): Restart detection or resume previous partial run.
              rc (ipyparallel.Client, optional): Remote Client, used for parallel execution

        """

        # Normally rc is assigned in init, but let's have option to get it here also
        if rc is not None:
            self.rc = rc

        # We need to set up the workers
        if self.rc is not None:
            d_view = self.rc.direct_view(targets='all')
        else:
            d_view = None

        if self.role == "master":

            # Make sure path exists
            if not os.path.exists(os.path.dirname(self.save_file)):
                self.write_log(f"Creating directory {os.path.dirname(self.save_file)}")
                os.mkdir(os.path.dirname(self.save_file))

            self.setup_parallel(d_view=d_view)

            if self.work_history_file is None:
                log_dir = os.path.join(self.network_path, "log")
                self.work_history_file = os.path.join(log_dir, "network-detect-worklog.hdf5")

            if restart_detection_flag:
                if os.path.isfile(self.work_history_file):
                    self.write_log("Removing old work history file")
                    os.remove(self.work_history_file)

                # Setup new work history
                self.setup_work_history(self.work_history_file)
            else:
                # Open old file with work history
                self.write_log(f"Reusing old work history file {self.work_history_file}")
                self.work_history = h5py.File(self.work_history_file, "r+", libver=self.h5libver)

            # For each neuron we need to find which hyper voxel it belongs to
            # (can be more than one)
            self.distribute_neurons_parallel(d_view=d_view)

            # We also need to start the projection code
            self.projection_detection = ProjectionDetection(snudda_detect=self, role=self.role, rc=self.rc)
            self.projection_detection.find_neurons_projections_in_hyper_voxels()

            if d_view is not None:
                self.parallel_process_hyper_voxels(rc=self.rc, d_view=d_view)

            else:
                # We are running it in serial

                (all_hyper_id, n_completed, remaining, self.voxel_overflow_counter) = \
                    self.setup_process_hyper_voxel_state_history()

                for hyper_id in remaining:  # self.hyperVoxels:
                    (hyper_id, n_syn, n_gj, exec_time, voxel_overflow_ctr) = \
                        self.process_hyper_voxel(hyper_id)

                    if voxel_overflow_ctr > 0:
                        self.write_log(f"!!! HyperID {hyper_id} OVERFLOWED {voxel_overflow_ctr} TIMES"
                                       f"{exec_time}s)")
                        self.voxel_overflow_counter += voxel_overflow_ctr
                    else:
                        self.write_log(f"HyperID {hyper_id} completed - {n_syn}  synapses and "
                                       f"{n_gj} gap junctions found (in {exec_time} s)")

                    self.update_process_hyper_voxel_state(hyper_id=hyper_id, num_syn=n_syn, num_gj=n_gj,
                                                          exec_time=exec_time,
                                                          voxel_overflow_counter=voxel_overflow_ctr)

        # We need to gather data from all the HDF5 files -- that is done in prune

    ############################################################################

    def __del__(self):

        if self.work_history is not None:
            try:
                self.work_history.close()
            except:
                print("Work history already closed")

        if self.logfile is not None:
            try:
                self.logfile.close()
            except:
                print("Log file already closed")

        if self.rc:
            # Clean up memory on workers
            from snudda.utils import cleanup
            cleanup(self.rc, "detect")

    ############################################################################

    def parallel_process_hyper_voxels(self, rc=None, d_view=None):

        """
        Distributes touch detection in hyper voxels to workers.

        Args:
            rc (ipyparallel.Client, optional): Remote client, for parallel execution

        """

        self.write_log("Starting parallel_process_hyper_voxels")

        start_time = timeit.default_timer()

        # Loads state if previously existed, otherwise creates new fresh history
        (all_hyper_id_list, num_completed, remaining, self.voxel_overflow_counter) = \
            self.setup_process_hyper_voxel_state_history()

        n_workers = len(rc.ids)
        worker_status = [None for x in rc.ids]
        worker_idx = 0
        job_idx = 0
        busy_ctr = 0
        no_change_ctr = 0
        num_syn = 1  # If nSyn is zero delay in loop is shorter

        self.write_log(f"parallel_process_hyper_voxels: Using {n_workers} worker")

        info_msg_written = False

        while job_idx < len(remaining) or busy_ctr > 0:

            if worker_status[worker_idx] is not None:

                # We have an async result, check status of it
                if worker_status[worker_idx].ready():

                    # Result is ready, get it
                    hyper_voxel_data = rc[worker_idx]["result"]

                    hyper_id = hyper_voxel_data[0]
                    num_syn = hyper_voxel_data[1]
                    n_gj = hyper_voxel_data[2]
                    exec_time = hyper_voxel_data[3]
                    voxel_overflow_ctr = hyper_voxel_data[4]

                    self.update_process_hyper_voxel_state(hyper_id=hyper_id,
                                                          num_syn=num_syn,
                                                          num_gj=n_gj,
                                                          exec_time=exec_time,
                                                          voxel_overflow_counter=voxel_overflow_ctr)
                    worker_status[worker_idx] = None
                    rc[worker_idx]["result"] = None  # Clear to be safe
                    busy_ctr -= 1

                    if voxel_overflow_ctr > 0:
                        self.write_log(f"!!! HyperID {hyper_id} OVERFLOWED {voxel_overflow_ctr} TIMES"
                                       f"(execution time {exec_time} s)", is_error=True)
                        self.voxel_overflow_counter += voxel_overflow_ctr
                    else:
                        if exec_time > 100 or self.verbose:
                            # Only print the hyper voxels taking a long time to complete
                            self.write_log(f"HyperID {hyper_id} completed "
                                           f"- {num_syn} synapses found ({np.around(exec_time, 1)} s)",
                                           force_print=True)
                        elif not info_msg_written:
                            self.write_log(f"Suppressing printouts for hyper voxels that complete in < 100 seconds.",
                                           force_print=True)
                            info_msg_written = True

                    if len(remaining) == 0:
                        # If there are no more hypervoxels to process, free the memory
                        # so that any workers still running can have more memory available
                        cmd_free_str = "sd.free_memory()"
                        rc["worker_idx"].execute(cmd_free_str, block=False)

            # Check that there are neurons in the hyper voxel, otherwise skip it.
            if worker_status[worker_idx] is None and job_idx < len(remaining):
                self.write_log(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}"
                               f" Starting hyper voxel {remaining[job_idx]} on worker {worker_idx}")

                cmd_str = f"result = sd.process_hyper_voxel({remaining[job_idx]})"
                worker_status[worker_idx] = rc[worker_idx].execute(cmd_str, block=False)

                job_idx += 1
                busy_ctr += 1
                no_change_ctr = 0

            else:
                no_change_ctr += 1

            worker_idx = (worker_idx + 1) % n_workers

            if no_change_ctr >= n_workers:
                # If all workers were working, sleep for a bit, then reset counter
                if num_syn > 0:
                    time.sleep(0.010)
                else:
                    # If previous hyper voxel had no synapses, be faster
                    time.sleep(1e-6)

                no_change_ctr = 0

        end_time = timeit.default_timer()

        self.write_log(f"Voxel overflows: {self.voxel_overflow_counter}", is_error=(self.voxel_overflow_counter > 0))
        self.write_log(f"Total number of synapses: {np.sum(self.work_history['num_hypervoxel_synapses'][:])}")
        self.write_log(f"parallel_process_hyper_voxels: {end_time - start_time:.1f} s")

        self.work_history.close()

    ############################################################################

    def generate_hyper_voxel_random_seeds(self):

        """ Generates a seed sequence for each hyper voxel based on the master seed for touch detection. """
        # https://albertcthomas.github.io/good-practices-random-number-generators/

        ss = np.random.SeedSequence(self.random_seed)
        all_seeds = ss.generate_state(len(self.hyper_voxels))
        all_hid = sorted(self.hyper_voxels.keys())

        for hi, s in zip(all_hid, all_seeds):
            self.hyper_voxels[hi]["random_seed"] = s

    def setup_hyper_voxel_id_lookup(self, max_coord, min_coord):

        self.num_hyper_voxels = np.ceil((max_coord - min_coord) / self.hyper_voxel_width).astype(int) + 1

        assert np.prod(self.num_hyper_voxels) < 3e9, \
            (f"Very large brain structure... Did you use SI units for neuron positions? "
             f"Number of hyper voxels: {self.num_hyper_voxels}")
        self.hyper_voxel_id_lookup = np.zeros(self.num_hyper_voxels, dtype=int)

        self.hyper_voxel_id_lookup[:] = \
            np.arange(0, self.hyper_voxel_id_lookup.size).reshape(self.hyper_voxel_id_lookup.shape)

        self.write_log(f"{self.hyper_voxel_id_lookup.size} hyper voxels in total")

    def get_hypervoxel_coords_and_section_id(self, neuron=None, neuron_info=None):

        """ Returns a tree info dictionary where each subtree is
            represented by an array with columns
            hyper_voxel_id, section_type, section_id"""

        # We do not check that neuron_info and neuron (if both give) are consistent.

        if neuron_info is None and neuron is None:
            raise KeyError("You must specify neuron_info or neuron")

        if neuron is None:
            neuron = self.load_neuron(neuron_info)

        tree_info = dict()
        hv_dim = np.array(self.hyper_voxel_id_lookup.shape)

        for subtree_name, subtree in neuron.morphology_data.items():
            hyper_voxel_coords = np.floor((subtree.geometry[:, :3] - self.simulation_origo[None, :]) / self.hyper_voxel_width).astype(int)

            # We need to do a range check, since virtual neurons and extra axons might be outside simulated region
            inside_idx = np.logical_and(0 <= hyper_voxel_coords,
                                        hyper_voxel_coords < hv_dim[None, :]).all(axis=1)

            hyper_voxel_id = self.hyper_voxel_id_lookup[tuple(hyper_voxel_coords[inside_idx, :].T)]
            section_type_id = subtree.section_data[inside_idx, :][:, [2, 0]]

            tree_info[subtree_name] = np.hstack([hyper_voxel_id.reshape([hyper_voxel_id.shape[0], 1]), section_type_id])

            # TODO: This should no longer be necessary! PARENT POINT should be included if
            #       parent section type is the same
            #
            # # We also need to add parent points with the child branch's section id
            # # This is to not miss the first bit between the parent point and the first real point of the branch
            # parent_rows = []
            # if inside_idx.all():
            #     # If inside_idx are all True then hyper_voxel_id is same length
            #     # as subtree.geometry (should be valid for all but possibly Virtual Axons)
            #
            #     for tree_type in subtree.sections:
            #         for section in subtree.sections[tree_type].values():
            #
            #             if section.section_type == 1:
            #                 # Soma has no parent, skip
            #                 continue
            #
            #             parent_idx = section.point_idx[0]
            #             parent_rows.append([hyper_voxel_id[parent_idx], section.section_type, section.section_id])
            #
            # if len(parent_rows) > 0:
            #     tree_info[subtree_name] = np.unique(np.vstack([hid_st_sid, parent_rows]), axis=0)
            # else:
            #     # This is the case if we only have a soma
            #     tree_info[subtree_name] = hid_st_sid

            # OBS, there is a rare case when a line segment starts in a hyper voxel, crosses a second hyper voxel
            # and ends up in a third hyper voxel. In this case the intermediate second hyper voxel will be missed if
            # there are no additional points from the section in that hyper voxel. This should be on the order of
            # 100 synapses per half a miljon synapses (if 40x40x40 hyper voxel size, ie very small hyper voxels)
            # For larger hyper voxels, the fraction of missed synapses should be even lower.

        # Check which hyper voxels the soma is in
        if "neuron" in neuron.morphology_data:
            morph_data = neuron.morphology_data["neuron"]
            if morph_data.section_data[0, 2] == 1:
                # First compartment is soma
                soma_pos = morph_data.geometry[0, :3]
                soma_radius = morph_data.geometry[0, 3]

                u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
                x = (soma_radius * np.cos(u) * np.sin(v) + soma_pos[0]).flatten()
                y = (soma_radius * np.sin(u) * np.sin(v) + soma_pos[1]).flatten()
                z = (soma_radius * np.cos(v) + soma_pos[2]).flatten()

                vx = np.floor((x - self.simulation_origo[0]) / self.hyper_voxel_width).astype(int)
                vy = np.floor((y - self.simulation_origo[1]) / self.hyper_voxel_width).astype(int)
                vz = np.floor((z - self.simulation_origo[2]) / self.hyper_voxel_width).astype(int)

                v_xyz = np.vstack([vx, vy, vz]).T

                inside_idx = np.logical_and(0 <= v_xyz, v_xyz < hv_dim[None, :]).all(axis=1)

                soma_hyper_voxels = np.unique(self.hyper_voxel_id_lookup[tuple(v_xyz[inside_idx, :].T)])

                soma_info = np.zeros((len(soma_hyper_voxels), 3), dtype=int)
                soma_info[:, 0] = soma_hyper_voxels
                soma_info[:, 1] = 1  # soma

                tree_info["neuron"] = np.unique(np.vstack([tree_info["neuron"], soma_info]), axis=0)

        # Columns in dictionary are: hyper_voxel_id, section_type, section_id

        return tree_info

    def group_section_info(self, tree_info):

        section_info = dict()

        for tree in tree_info:

            for hyper_id, section_type, section_id in tree_info[tree]:

                if hyper_id not in section_info:
                    section_info[hyper_id] = dict()

                if section_type not in section_info[hyper_id]:
                    section_info[hyper_id][section_type] = dict()

                if tree not in section_info[hyper_id][section_type]:
                    section_info[hyper_id][section_type][tree] = [section_id]
                else:
                    section_info[hyper_id][section_type][tree].append(section_id)

        # Sort the sections to make it reproducible
        for hyper_id in section_info:
            for section_type in section_info[hyper_id]:
                for tree in section_info[hyper_id][section_type]:
                    section_info[hyper_id][section_type][tree] = \
                        np.array(sorted(section_info[hyper_id][section_type][tree]))

        return section_info

    def get_density_location(self, neuron, seed):

        axon_loc = None
        hyper_voxel_id = []

        if neuron.axon_density_type == "r":

            rng = np.random.default_rng(seed)

            # We create a set of points corresponding approximately to the
            # extent of the axonal density, and check which hyper voxels
            # they occupy

            # Radius of sphere in hyper voxels, rounded up
            rad = np.ceil(neuron.max_axon_radius / (self.hyper_voxel_size * self.voxel_size))

            # Approximately how many hyper voxels will the dendritic tree occupy
            n_hv = (2 * rad) ** 3

            # Over sample
            num_points = int(30 * n_hv)

            # Randomly place these many points within a sphere of the given radius
            # and then check which hyper voxels these points belong to

            theta = 2 * np.pi * rng.random(num_points)
            phi = np.arccos(2 * rng.random(num_points) - 1)
            r = neuron.max_axon_radius * (rng.random(num_points) ** (1 / 3))

            x = np.multiply(r, np.multiply(np.sin(phi), np.cos(theta)))
            y = np.multiply(r, np.multiply(np.sin(phi), np.sin(theta)))
            z = np.multiply(r, np.cos(phi))

            axon_cloud = np.zeros((len(x), 3))
            axon_cloud[:, 0] = x + neuron.position[0]
            axon_cloud[:, 1] = y + neuron.position[1]
            axon_cloud[:, 2] = z + neuron.position[2]

            axon_loc = np.floor((axon_cloud[:, :3] - self.simulation_origo) / self.hyper_voxel_width).astype(int)

        elif neuron.axon_density_type == "xyz":

            # TODO: Maybe replace random points by a grid for this test step?

            rng = np.random.default_rng(seed)

            # Estimate how many points we need to randomly place
            num_points = 100 * np.prod(neuron.axon_density_bounds_xyz[1:6:2]
                                       - neuron.axon_density_bounds_xyz[0:6:2]) \
                         / ((self.hyper_voxel_size * self.voxel_size) ** 3)
            num_points = int(np.ceil(num_points))

            if num_points > 1e4:
                self.write_log(f"!!! Many many points placed for axon density of {neuron.name} : {num_points}")

            xmin = neuron.axon_density_bounds_xyz[0]
            xwidth = neuron.axon_density_bounds_xyz[1] - neuron.axon_density_bounds_xyz[0]
            ymin = neuron.axon_density_bounds_xyz[2]
            ywidth = neuron.axon_density_bounds_xyz[3] - neuron.axon_density_bounds_xyz[2]
            zmin = neuron.axon_density_bounds_xyz[4]
            zwidth = neuron.axon_density_bounds_xyz[5] - neuron.axon_density_bounds_xyz[4]

            # The purpose of this is to find out the range of the axon bounding box
            axon_cloud = rng.random((num_points, 3))
            axon_cloud[:, 0] = axon_cloud[:, 0] * xwidth + xmin
            axon_cloud[:, 1] = axon_cloud[:, 1] * ywidth + ymin
            axon_cloud[:, 2] = axon_cloud[:, 2] * zwidth + zmin

            # Don't forget to rotate
            axon_cloud = np.matmul(neuron.rotation,
                                   axon_cloud.transpose()).transpose() + neuron.position

            axon_loc = np.floor((axon_cloud[:, :3] - self.simulation_origo) / self.hyper_voxel_width).astype(int)

        if axon_loc is not None:
            inside_idx = np.sum(np.logical_and(0 <= axon_loc, axon_loc < self.hyper_voxel_id_lookup.shape), axis=1) == 3
            hyper_voxel_id = np.unique(self.hyper_voxel_id_lookup[tuple(axon_loc[inside_idx, :].T)])

        return hyper_voxel_id

    def preallocate_empty_hyper_voxel_dict(self):

        # First assign hyperVoxelID to the space
        self.hyper_voxels = dict([])

        for ix in range(0, self.num_hyper_voxels[0]):
            for iy in range(0, self.num_hyper_voxels[1]):
                for iz in range(0, self.num_hyper_voxels[2]):
                    hid = int(self.hyper_voxel_id_lookup[ix, iy, iz])

                    self.hyper_voxels[hid] = dict([])
                    self.hyper_voxels[hid]["origo"] = (self.simulation_origo
                                                       + self.hyper_voxel_width * np.array([ix, iy, iz]))

                    self.hyper_voxels[hid]["neurons"] = dict()

                    # Changed so we preallocate only empty, to preserve memory
                    # self.hyper_voxels[hid]["soma"] = []
                    # self.hyper_voxels[hid]["axon"] = []
                    # self.hyper_voxels[hid]["dend"] = []

                    self.hyper_voxels[hid]["axon_density"] = []

                    # self.hyper_voxels[hid]["neurons"] = []
                    self.hyper_voxels[hid]["neuron_ctr"] = 0

        self.write_log("Pre allocation done.")

    ############################################################################

    def generate_neuron_distribution_random_seeds(self):

        """ Generate seed sequence for neuron distribution from master seed. """
        # Need different master seed than hyper voxel seed sequence
        ss = np.random.SeedSequence(self.random_seed + 1337)
        distribution_seeds = ss.generate_state(len(self.neurons))
        return distribution_seeds

    ############################################################################

    def setup_work_history(self, work_history_file=None):

        """
        Sets up work history. By logging progress we are able to restart an earlier partial touch detection run.

        Args:
            work_history_file (str, optional): Path to work history file

        """

        if self.role != "master":
            return

        if work_history_file is None:
            work_history_file = self.work_history_file
        else:
            # Update internal state
            self.work_history_file = work_history_file

        dir_name = os.path.dirname(self.work_history_file)
        if not os.path.exists(dir_name):
            self.write_log(f"Creating directory {dir_name}")
            os.mkdir(dir_name)

        self.write_log(f"Work history file: {self.work_history_file}")

        self.work_history = h5py.File(work_history_file, "w", libver=self.h5libver)

        # We need to encode the connectivityDistributions tuple as a string
        # before saving it in json
        # we also need to parse this string and recreate a tuple afterwards

        tmp_con_dist = dict([])

        for keys in self.connectivity_distributions:
            tmp_con_dist["$$".join(keys)] = self.connectivity_distributions[keys]

        save_meta_data = [(self.snudda_data, "snudda_data"),
                          (self.slurm_id, "slurm_id"),
                          (self.config_file, "config_file"),
                          (self.position_file, "position_file"),
                          (self.voxel_size, "voxel_size"),
                          (self.hyper_voxel_size, "hyper_voxel_size"),
                          (self.hyper_voxel_width, "hyper_voxel_width"),
                          (json.dumps(self.config), "config"),
                          (json.dumps(tmp_con_dist), "connectivity_distributions")]

        if "meta" not in self.work_history:
            self.write_log("Writing metadata to work history file")
            meta = self.work_history.create_group("meta")

            for data, data_name in save_meta_data:
                meta.create_dataset(data_name, data=data)

        else:
            self.write_log("Work history file found, checking meta data")
            # Make sure config file etc match

            for data, data_name in save_meta_data:
                assert data == self.work_history["meta/" + data_name][()], \
                    (f"Mismatch in workhistory file: {data_name} mismatch {data} "
                     f"vs {self.work_history['meta/' + data_name][()]}"
                     f"\nPlease delete old work history file.")

        self.write_log("Write neuron data to file")

        network_group = self.work_history.create_group("network")

        # Finally the neuron information from the positon data file
        with h5py.File(self.position_file, "r") as f:
            f.copy("network/neurons", network_group)

    ############################################################################

    # Reading work

    def setup_process_hyper_voxel_state_history(self):

        """
        Initialises the variables for tracking progress of touch detection.

        Returns:
            progress_data (list, int, list, int) : The items returned are (all_hyper_id_list, num_completed,
                                                   remaining, voxel_overflow_counter). The remaining hypervoxel id list
                                                   is sorted descending by size.

        """

        if "completed" in self.work_history:
            self.write_log("setup_process_hyper_voxel_state_history: Resuming from old state")
            # We already have a run in progress, load the state
            all_hyper_id_list = set(self.work_history["all_hyper_ids"])
            num_completed = int(self.work_history["num_completed"][0])
            completed = set(self.work_history["completed"][:num_completed])
            remaining = self.sort_remaining_by_size(all_hyper_id_list - completed)
            voxel_overflow_counter = self.work_history["voxel_overflow_counter"][0]

        else:
            self.write_log("setup_process_hyper_voxel_state_history: Creating new work history.")
            # No history, add it to work history file
            num_hyper_voxels = len(self.hyper_voxels)
            minus_one = -1 * np.ones((num_hyper_voxels,), dtype=np.int32)
            self.work_history.create_dataset("completed", data=minus_one)

            num_completed = 0
            voxel_overflow_counter = 0

            # Could not rewrite scalars, so saving num_completed as a vector of length 1
            self.work_history.create_dataset("num_completed", data=np.zeros(1, ))
            all_hyper_id_list = np.array([x for x in self.hyper_voxels.keys()], dtype=np.int32)

            # Remove the empty hyper IDs
            (valid_hyper_id, empty_hyper_id) = self.remove_empty(all_hyper_id_list)
            all_hyper_id_list = valid_hyper_id
            remaining = self.sort_remaining_by_size(all_hyper_id_list)

            if len(self.connectivity_distributions) == 0:
                # We have no possible connections specified -- mark all voxels as done
                self.write_log("No connections specified in connectivity_distribution.", is_error=True)
                remaining = []
            else:

                assert (np.array([self.hyper_voxels[x]["neuron_ctr"] for x in
                                  empty_hyper_id]) == 0).all(), "All hyperIDs marked as empty are not empty!"

                self.write_log(f"Skipping {len(empty_hyper_id)} empty hyper voxels")

            self.work_history.create_dataset("all_hyper_ids", data=all_hyper_id_list)
            self.work_history.create_dataset("num_hypervoxel_synapses",
                                             data=np.zeros(num_hyper_voxels, ), dtype=np.int64)
            self.work_history.create_dataset("num_hypervoxel_gap_junctions",
                                             data=np.zeros(num_hyper_voxels, ), dtype=np.int64)
            self.work_history.create_dataset("voxel_overflow_counter", data=np.zeros(num_hyper_voxels, ), dtype=np.int64)

        return all_hyper_id_list, num_completed, remaining, voxel_overflow_counter

    ############################################################################

    # We want to do the hyper voxels with most neurons first, to minimize
    # the time waiting for lone cpu worker stragglers at the end.

    def sort_remaining_by_size(self, remaining):

        """
        Sorts the remaining hypervoxel ID list by descending size (number of neurons in hypervoxel)

        Args:
            remaining (list): List of hypervoxels

        Returns:
            sorted_remaining (list): Sorted list of hypervoxels
        """

        remaining = np.array(list(remaining), dtype=int)

        # Minus since we want them in descending order
        num_neurons = [-self.hyper_voxels[x]["neuron_ctr"] for x in remaining]
        sort_idx = np.argsort(num_neurons)

        return remaining[sort_idx]

    ############################################################################

    def remove_empty(self, hyper_id):
        """
        Removes empty hypervoxels from the list

        Args:
            hyper_id (list): List of hyper voxel IDs

        Returns:
            hyper_id_kept, hyper_id_removed (list, list) : List of remaining, and removed, hyper voxel ID

        """

        num_neurons = np.array([self.hyper_voxels[x]["neuron_ctr"] for x in hyper_id])
        keep_idx = np.where(num_neurons > 0)[0]
        remove_idx = np.where(num_neurons == 0)[0]

        return hyper_id[keep_idx], hyper_id[remove_idx]

    ############################################################################

    def get_neuron_distribution_history(self):

        """ Returns info about what neurons each hyper voxel contains etc.

        Returns:
            (tuple) : containing
                hyper_voxels (dictionary): dictionary with keys 'neurons', 'neuron_ctr', 'origo', 'random_seed'
                hyper_voxel_id_lookup (3D matrix with int): hypervoxel ID, spatially arranged
                n_hyper_voxels (int): number of hypervoxels
                simulation_origo (float, float, float): origo of entire simulation
        """

        if "hyper_voxels" in self.work_history:
            self.write_log("Using neuron distribution from work history.")

            hyper_voxels_str = json.loads(self.work_history["hyper_voxels"][()])
            hyper_voxels = dict()

            # When we load from the JSON file the integer keys have become strings, convert keys back to int
            for k, v in hyper_voxels_str.items():
                hyper_voxels[int(k)] = v

            hyper_voxel_id_lookup = self.work_history["meta/hyper_voxel_ids"][()]
            n_hyper_voxels = self.work_history["meta/num_hyper_voxels"][()]
            simulation_origo = self.work_history["meta/simulation_origo"][()]

            return hyper_voxels, hyper_voxel_id_lookup, n_hyper_voxels, simulation_origo
        else:
            # No information stored
            return None, None, None, None

    ############################################################################

    def save_neuron_distribution_history(self, hyper_voxels, min_coord, max_coord):

        """
        Save neuron distribution history to file.

        Args:
            hyper_voxels (3D matrix with int): Hypervoxel IDs, spatially arranged
            min_coord (float,float,float): minimal coordinate for all neurons/neurites in simulation
            max_coord (float,float,float): maximal coordinate for all neurons/neurites in simulation

        """

        # TODO: min_coord, max_coord not saved, no longer used? remove from input parameters?

        self.write_log("Writing neuron distribution history to file")

        assert "hyper_voxels" not in self.work_history, "save_neuron_distribution_history should only be called once"

        self.work_history.create_dataset("meta/hyper_voxel_ids", data=self.hyper_voxel_id_lookup)
        self.work_history.create_dataset("meta/num_hyper_voxels", data=self.num_hyper_voxels)
        self.work_history.create_dataset("meta/simulation_origo", data=self.simulation_origo)

        self.write_log("Generating JSON hyper voxel data")

        mem = self.memory()
        self.write_log(f"{mem}")

        try:
            hyper_voxels = json.dumps(self.hyper_voxels, cls=NumpyEncoder)
        except:
            import traceback
            self.write_log(traceback.format_exc(), is_error=True)
            import pdb
            pdb.set_trace()

        self.write_log("Writing JSON hyper voxel data to file")

        self.work_history.create_dataset("hyper_voxels", data=hyper_voxels)

        self.write_log("Save neuron distribution history done.")

    ############################################################################

    def update_process_hyper_voxel_state(self, hyper_id, num_syn, num_gj, exec_time, voxel_overflow_counter):

        """Updates the process log with new hypervoxel state

        Args:
            hyper_id (int) : Hypervoxel id completed
            num_syn (int) : Number of synapses detected in hyper voxel
            num_gj (int) : Number of gap junctions detected in hyper voxel
            exec_time : Execution time (currently not used!)
            voxel_overflow_counter : How many synapses/gap junctions did we miss due to memory overflow? (Should be 0)

        """

        num_completed = int(self.work_history["num_completed"][0])

        self.work_history["completed"][num_completed] = hyper_id
        self.work_history["num_hypervoxel_synapses"][num_completed] = num_syn
        self.work_history["num_hypervoxel_gap_junctions"][num_completed] = num_gj
        self.work_history["voxel_overflow_counter"][num_completed] = voxel_overflow_counter

        num_completed += 1
        self.work_history["num_completed"][0] = num_completed

    ############################################################################

    def setup_hyper_voxel(self, hyper_voxel_origo, hyper_voxel_id):

        """
        Initialise all variables for a hypervoxel (containing NxNxN voxels) before touch detection.

        Args:
            hyper_voxel_origo (float,float,float): Origo of new hyper voxel
            hyper_voxel_id (int): ID of hypervoxel to process
        """

        # Each hyper voxel has its own seed
        random_seed = self.hyper_voxels[hyper_voxel_id]["random_seed"]
        self.hyper_voxel_rng = np.random.default_rng(random_seed)

        self.hyper_voxel_coords[hyper_voxel_id] = hyper_voxel_origo

        self.hyper_voxel_origo = hyper_voxel_origo
        self.hyper_voxel_id = hyper_voxel_id

        # Clear lookup tables, just to be safe
        self.hyper_voxel_synapse_lookup = None
        self.hyper_voxel_gap_junction_lookup = None

        if SnuddaDetect.memory_fraction_free() < 0.2:
            # Clear some variables to free memory, reset max values to default, and perform garbage collection
            self.hyper_voxel_synapses = None
            self.hyper_voxel_gap_junctions = None
            self.max_synapses = 2000000
            self.max_gap_junctions = 100000
            self.neuron_cache = dict([])
            self.extra_axon_cache = dict([])
            gc.collect()

        if self.hyper_voxel_synapses is None:
            self.hyper_voxel_synapses = np.zeros((self.max_synapses, 13), dtype=np.int32)
            self.hyper_voxel_synapse_ctr = 0
        else:
            self.hyper_voxel_synapses[:] = 0
            self.hyper_voxel_synapse_ctr = 0

        if self.hyper_voxel_gap_junctions is None:
            self.hyper_voxel_gap_junctions = np.zeros((self.max_synapses, 11),
                                                      dtype=np.int32)
            self.hyper_voxel_gap_junction_ctr = 0
        else:
            self.hyper_voxel_gap_junctions[:] = 0
            self.hyper_voxel_gap_junction_ctr = 0

        # Used by plotHyperVoxel to make sure synapses are displayed correctly
        self.hyper_voxel_offset = None

        # Which axons populate the different voxels
        if self.axon_voxels is None:
            self.axon_voxels = np.zeros((self.num_bins[0],
                                         self.num_bins[1],
                                         self.num_bins[2],
                                         self.max_axon),
                                        dtype=np.int32)
            self.axon_voxel_ctr = np.zeros(self.num_bins, dtype=np.int32)

            # How far from the soma is this point
            self.axon_soma_dist = np.zeros((self.num_bins[0],
                                            self.num_bins[1],
                                            self.num_bins[2],
                                            self.max_axon),
                                           dtype=np.int16)
        else:
            # Already allocated, just clear it
            self.axon_voxels[:] = 0
            self.axon_voxel_ctr[:] = 0
            self.axon_soma_dist[:] = 0

        # Which dendrites populate the different voxels
        if self.dend_voxels is None:
            self.dend_voxels = np.zeros((self.num_bins[0],
                                         self.num_bins[1],
                                         self.num_bins[2],
                                         self.max_dend),
                                        dtype=np.int32)
            self.dend_voxel_ctr = np.zeros(self.num_bins, dtype=np.int32)

            # Which segment ID does the point belong to, and what segX
            self.dend_sec_id = np.zeros((self.num_bins[0],
                                         self.num_bins[1],
                                         self.num_bins[2],
                                         self.max_dend),
                                        dtype=np.int16)
            self.dend_sec_x = np.zeros((self.num_bins[0],
                                        self.num_bins[1],
                                        self.num_bins[2],
                                        self.max_dend),
                                       dtype=np.float64)  # 0 - 1.0, low pres   # float16 -> float64, numba requirement

            # How far from the soma is this point
            self.dend_soma_dist = np.zeros((self.num_bins[0],
                                            self.num_bins[1],
                                            self.num_bins[2],
                                            self.max_dend),
                                           dtype=np.int16)

        else:
            # Already allocated, just clear it
            self.dend_voxels[:] = 0
            self.dend_voxel_ctr[:] = 0
            self.dend_sec_id[:] = 0
            self.dend_sec_x[:] = 0
            self.dend_soma_dist[:] = 0

        self.voxel_overflow_counter = 0

        if self.role == "worker":
            # Let's clear the cache between hyper voxels if we are running in parallel
            # (less chance of cache hits between hyper voxels, and avoid too many copies cached)
            self.neuron_cache = dict()
            self.extra_axon_cache = dict()

    def free_memory(self):
        # Clear some variables to free memory, reset max values to default, and perform garbage collection
        self.hyper_voxel_synapses = None
        self.hyper_voxel_gap_junctions = None
        self.max_synapses = 2000000
        self.max_gap_junctions = 100000

        self.neuron_cache = dict([])
        self.extra_axon_cache = dict([])

        self.axon_voxels = None
        self.axon_voxel_ctr = None
        self.axon_soma_dist = None

        self.dend_voxels = None
        self.dend_voxel_ctr = None
        self.dend_sec_id = None
        self.dend_sec_x = None
        self.dend_soma_dist = None

        mem_available_before = self.memory_fraction_free()

        gc.collect()

        mem_available_after = self.memory_fraction_free()

        self.writelog(f"Hyper voxel memory freed. Free before {mem_available_before}%, free after {mem_available_after}%")

    ############################################################################

    # hyperID is only needed if we have neurons without axons, ie we use
    # axon density

    def detect_synapses(self):

        """ Helper function, triggers detection of synapses. Called by process_hyper_voxel. """

        start_time = timeit.default_timer()

        # assert self.hyperVoxelSynapseCtr == 0 \
        #   and self.hyperVoxelSynapses is not None, \
        #   "setupHyperVoxel must be called before detecting synapses"

        # Find all voxels that contain axon and dendrites
        [x_syn, y_syn, z_syn] = np.where(np.bitwise_and(self.dend_voxel_ctr > 0,
                                                        self.axon_voxel_ctr > 0))

        if True:
            # This gives us some statistics, turn off later for speed
            self.max_axon_voxel_ctr = np.amax(self.axon_voxel_ctr)
            self.max_dend_voxel_ctr = np.amax(self.dend_voxel_ctr)

        for x, y, z in zip(x_syn, y_syn, z_syn):
            axon_id_list = self.axon_voxels[x, y, z, :self.axon_voxel_ctr[x, y, z]]
            dend_id_list = self.dend_voxels[x, y, z, :self.dend_voxel_ctr[x, y, z]]

            axon_dist = self.axon_soma_dist[x, y, z, :self.axon_voxel_ctr[x, y, z]]
            dend_dist = self.dend_soma_dist[x, y, z, :self.dend_voxel_ctr[x, y, z]]

            dend_sec_id = self.dend_sec_id[x, y, z, :self.dend_voxel_ctr[x, y, z]]
            dend_sec_x = self.dend_sec_x[x, y, z, :self.dend_voxel_ctr[x, y, z]]

            # Maybe make dendrite loop outer, since it has more variables?
            # speedup??
            for (ax_id, ax_dist) in zip(axon_id_list, axon_dist):
                for (d_id, d_sec_id, d_sec_x, d_dist) \
                        in zip(dend_id_list, dend_sec_id, dend_sec_x, dend_dist):

                    if ax_id == d_id:
                        # Avoid self connections
                        continue

                    pre_type = self.neurons[ax_id]["type"]
                    post_type = self.neurons[d_id]["type"]

                    if (pre_type, post_type) in self.connectivity_distributions:

                        con_dict = self.connectivity_distributions[pre_type, post_type]

                        # We need to loop over conDict in case there are multiple
                        # types of synapses from this neuron
                        for con_type in con_dict:
                            if con_type == "gap_junction":
                                # This part detects only axon-dend synapses, skip gap junctions
                                continue

                            synapse_mu, synapse_sigma = con_dict[con_type]["lognormal_mu_sigma"]
                            mean_synapse_cond, std_synapse_cond = con_dict[con_type]["conductance"]
                            channel_model_id = con_dict[con_type]["channel_model_id"]

                            # Should we add just one synapse, or a cluster of synapses
                            cluster_size = con_dict[con_type]["cluster_size"]

                            if isinstance(cluster_size, (np.ndarray, list)):
                                cluster_size = round(self.hyper_voxel_rng.normal(loc=cluster_size[0],
                                                                                 scale=cluster_size[1]))
                                
                            if cluster_size > 1:
                                cluster_spread = con_dict[con_type]["cluster_spread"]

                                if cluster_spread is None:
                                    raise ValueError(f"cluster_spread ({cluster_spread}) must be a distance (in meters), e.g. 10e-6. ")
                                
                                if isinstance(cluster_spread, (np.ndarray, list)):
                                    cluster_spread = np.maximum(np.abs(self.hyper_voxel_rng.normal(loc=cluster_spread[0],
                                                                                                   scale=cluster_spread[1])),
                                                                5e-6)

                                # This uses clone in neuron_prototype which should be cached (not anymore, but will be cached for 2nd hit)
                                neuron = self.load_neuron(self.neurons[d_id])
                                
                                try:
                                    cluster_sec_x, syn_coords, soma_dist \
                                        = neuron.cluster_synapses(sec_id=d_sec_id, sec_x=d_sec_x,
                                                                  count=cluster_size, distance=cluster_spread,
                                                                  rng=self.hyper_voxel_rng)
                                except:
                                    import traceback
                                    tstr = traceback.format_exc()
                                    print(tstr)
                                    import pdb
                                    pdb.set_trace()

                                if self.hyper_voxel_synapse_ctr + cluster_size >= self.max_synapses:
                                    self.resize_hyper_voxel_synapses_matrix()

                                cluster_cond = self.hyper_voxel_rng.lognormal(synapse_mu, synapse_sigma, cluster_size)
                                cluster_cond = np.maximum(cluster_cond, mean_synapse_cond * 0.1)
                                cluster_param_id = self.hyper_voxel_rng.integers(1000000, size=cluster_size)

                                # We need to convert coords to hyper voxel coords, to fit with other coords
                                coords_all = np.floor((syn_coords - self.hyper_voxel_origo)/self.voxel_size).astype(int)

                                for d_sec_x, x, y, z, d_dist, cond, param_id \
                                        in zip(cluster_sec_x, coords_all[:, 0], coords_all[:, 1], coords_all[:, 2],
                                               soma_dist * 1e6, cluster_cond, cluster_param_id):
                                    assert cond > 0, f"Conductance should be larger than 0. cond = {cond}"

                                    self.hyper_voxel_synapses[self.hyper_voxel_synapse_ctr, :] = \
                                        [ax_id, d_id, x, y, z, self.hyper_voxel_id, channel_model_id,
                                         ax_dist, d_dist, d_sec_id, d_sec_x * 1000, cond * 1e12, param_id]

                                    # !!! OBS, dSegX is a value between 0 and 1, multiplied by 1000
                                    # need to divide by 1000 later

                                    self.hyper_voxel_synapse_ctr += 1

                            else:
                                # We can not do pruning at this stage, since we only see
                                # synapses within hyper voxel, and pruning depends on
                                # all synapses between two connected cells.

                                # Do we have enough space allocated?
                                if self.hyper_voxel_synapse_ctr >= self.max_synapses:
                                    self.resize_hyper_voxel_synapses_matrix()

                                # Synapse conductance varies between synapses
                                # cond = self.hyper_voxel_rng.normal(mean_synapse_cond, std_synapse_cond)

                                # lognormal distribution -- https://www.nature.com/articles/nrn3687
                                # https://en.wikipedia.org/wiki/Log-normal_distribution
                                cond = self.hyper_voxel_rng.lognormal(synapse_mu, synapse_sigma)

                                # Need to make sure the conductance is not negative,
                                # set lower cap at 10% of mean value
                                cond = np.maximum(cond, mean_synapse_cond * 0.1)
                                assert cond > 0, f"Conductance should be larger than 0. cond = {cond}"

                                param_id = self.hyper_voxel_rng.integers(1000000)

                                # Add synapse
                                self.hyper_voxel_synapses[self.hyper_voxel_synapse_ctr, :] = \
                                    [ax_id, d_id, x, y, z, self.hyper_voxel_id, channel_model_id,
                                     ax_dist, d_dist, d_sec_id, d_sec_x * 1000, cond * 1e12, param_id]

                                # !!! OBS, dSegX is a value between 0 and 1, multiplied by 1000
                                # need to divide by 1000 later

                                self.hyper_voxel_synapse_ctr += 1

        # Sort the synapses (note sortIdx will not contain the empty rows
        # at the end.

        self.sort_synapses()

        # Convert from hyper voxel local coordinates to simulation coordinates
        # basically how many voxel steps do we need to take to go from
        # simulationOrigo to hyperVoxelOrigo (those were not included, so add them)
        hyper_voxel_offset = np.round((self.hyper_voxel_origo - self.simulation_origo)
                                      / self.hyper_voxel_width).astype(int) * self.hyper_voxel_size

        # Just a double check...
        assert self.hyper_voxel_id_lookup[int(np.round(hyper_voxel_offset[0] / self.hyper_voxel_size))][
                   int(np.round(hyper_voxel_offset[1] / self.hyper_voxel_size))][
                   int(np.round(hyper_voxel_offset[2] / self.hyper_voxel_size))] == self.hyper_voxel_id, \
            "Internal inconsistency, have hyper voxel numbering or coordinates been changed?"

        self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr, :][:, range(2, 5)] \
            += hyper_voxel_offset

        # We need this in case plotHyperVoxel is called
        self.hyper_voxel_offset = hyper_voxel_offset

        # These are used when doing the heap sort of the hyper voxels
        self.hyper_voxel_synapse_lookup \
            = self.create_lookup_table(data=self.hyper_voxel_synapses,
                                       n_rows=self.hyper_voxel_synapse_ctr,
                                       data_type="synapses",
                                       num_neurons=len(self.neurons),
                                       max_synapse_type=self.next_channel_model_id)

        # if(self.hyperVoxelSynapseCtr > 0 and self.hyperVoxelSynapseCtr < 10):
        #  self.plotHyperVoxel()
        #  import pdb
        #  pdb.set_trace()

        end_time = timeit.default_timer()

        self.write_log(f"detect_synapses: {self.hyper_voxel_synapse_ctr} took {end_time - start_time:.1f} s")

        if False and self.hyper_voxel_synapse_ctr > 0:
            print("First plot shows dendrites, and the voxels that were marked")
            print("Second plot same, but for axons")
            self.plot_hyper_voxel(plot_neurons=True, draw_axons=False)
            self.plot_hyper_voxel(plot_neurons=True, draw_dendrites=False)
            # This is for debug purposes
            import pdb
            pdb.set_trace()

        return self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr, :]

    ############################################################################

    def place_synapses_no_axon(self, hyper_id, voxel_space, voxel_space_ctr,
                               voxel_axon_dist):

        """
        Places fake axon segments for neurons without axons.

        Args:
            hyper_id (int): Hypervoxel ID
            voxel_space: Axon voxel space (n_bins x n_bins x n_bins x max_syn, voxel space matrix)
            voxel_space_ctr: Synapse counter (int) for voxels (n_bins x n_bins x n_bins)
            voxel_axon_dist: Axonal distance from soma to synapses (n_bins x n_bins x n_bins)

        """

        start_time = timeit.default_timer()

        # 1. Find neurons within hyper voxel that have no axon
        no_axon_neurons = self.hyper_voxels[hyper_id]["axon_density"]

        if len(no_axon_neurons) == 0:
            # No neurons without axons
            return

        for na_neuron_id in no_axon_neurons:

            # There are two types of axon density specified
            # - Spherically symmetric
            # - f(x,y,z) in SWC coordinates

            na_neuron = self.neurons[na_neuron_id]

            if na_neuron["axon_density_type"] == "r":

                # 2. Check that we have cumulative probability distribution for
                #    radial distance, if not compute and cache

                if na_neuron["type"] in self.axon_cum_density_cache:
                    (na_cum_density, na_points) = self.axon_cum_density_cache[na_neuron["type"]]
                    self.write_log(f"Placing {na_points} random axon points for {na_neuron['name']} (cached)")

                else:
                    # r is used in numexpr.evaluate below
                    r = radius = np.arange(0, na_neuron["axon_density_radius"] + self.voxel_size, self.voxel_size)

                    # old way using eval, replaced with numpexpr.evaluate (safer, faster?)
                    # density_as_func = eval('lambda r: ' + na_neuron["axon_density"])
                    # na_p_density = np.array([density_as_func(r) for r in radius])

                    na_p_density = numexpr.evaluate(na_neuron["axon_density"])

                    # We need to scale by distance squared, since in the shell at distance
                    # d further from the soma has more voxels in it than a shell closer
                    # This cumulative distribution is only used to determine how far
                    # from the soma a synapse is located (not the direction)

                    # !!! Plot and verify this !!!
                    na_cum_density = np.cumsum(np.multiply(na_p_density, radius ** 2))
                    na_cum_density /= na_cum_density[-1]  # Normalise

                    # 3. Calculate how many points there should be within volume
                    #    based on (unscaled raw) probability density
                    # Volume at each distance is 4*pi*(r**2) * voxel_size
                    na_points = int(np.round(np.sum(4 * np.pi * self.voxel_size
                                                    * np.multiply(radius ** 2, na_p_density))))

                    self.write_log(f"Placing {na_points} random axon points for {na_neuron['name']}")

                    self.axon_cum_density_cache[na_neuron["type"]] = (na_cum_density, na_points)

                # 4. Randomize the points
                (na_voxel_coords, na_axon_dist) = self.no_axon_points_sphere(na_neuron["position"],
                                                                             na_cum_density,
                                                                             na_points)

            elif na_neuron["axon_density_type"] == "xyz":

                (na_voxel_coords, na_axon_dist) = self.no_axon_points_xyz(na_neuron["position"],
                                                                          na_neuron["rotation"],
                                                                          na_neuron["axon_density"],
                                                                          na_neuron["axon_density_bounds_xyz"])
            else:
                self.write_log(f"Unknown axon_density_type: {na_neuron['axon_density_type']}\n{na_neuron}", is_error=True)
                na_voxel_coords = np.zeros((0, 3))
                na_axon_dist = []

            neuron_id = na_neuron["neuron_id"]

            for idx in range(0, na_voxel_coords.shape[0]):
                x_idx = na_voxel_coords[idx, 0]
                y_idx = na_voxel_coords[idx, 1]
                z_idx = na_voxel_coords[idx, 2]
                axon_dist = na_axon_dist[idx]

                v_ctr = voxel_space_ctr[x_idx, y_idx, z_idx]
                if v_ctr > 0 and voxel_space[x_idx, y_idx, z_idx, v_ctr - 1] == neuron_id:
                    # Voxel already has neuron_id, skip
                    continue

                try:
                    voxel_space[x_idx, y_idx, z_idx, v_ctr] = neuron_id
                    voxel_axon_dist[x_idx, y_idx, z_idx, v_ctr] = axon_dist
                    voxel_space_ctr[x_idx, y_idx, z_idx] += 1
                except:
                    self.voxel_overflow_counter += 1
                    self.write_log(f"!!! Axon voxel space overflow: {voxel_space_ctr[x_idx, y_idx, z_idx]}",
                                   is_error=True)

            # if(True):
            #  # Debug plot
            #  self.plotHyperVoxel()
            #  import pdb
            #  pdb.set_trace()

        end_time = timeit.default_timer()

        self.write_log(f"place_synapses_no_axon: {end_time - start_time:.1f} s, hyper_id: {hyper_id}")

    ############################################################################

    def no_axon_points_sphere(self, soma_centre, r_cum_distribution, num_points):

        """
        Helper function placing axon segments with spherical probability distribution.
        This picks points around soma centre. num_points are randomized, points
        outside the hyper sphere are rejected, so fewer than nPoints might be returned.

        Args:
            soma_centre (float,float,float): x,y,z coordinates of soma centre
            r_cum_distribution: cumulative distribution
            num_points: number of points to place

        Returns:
            voxel_coordinates
            synapse_distance_to_soma
        """

        uvr = self.hyper_voxel_rng.random((num_points, 3))
        theta = 2 * np.pi * uvr[:, 0]
        phi = np.arccos(2 * uvr[:, 1] - 1)

        # Double check these are sorted
        # We want to sample from the supplied distance distribution
        r_p = np.sort(uvr[:, 2] * r_cum_distribution[-1], axis=0)
        next_idx = 0

        self.write_log(f"num_points = {num_points}")

        r = np.zeros((num_points,))

        for posIdx, rP1 in enumerate(r_p):
            while rP1 > r_cum_distribution[next_idx + 1]:
                next_idx += 1

            r[posIdx] = next_idx

        # Rescale index to radie
        r = r * self.voxel_size

        xyz = np.array([r * np.sin(theta) * np.cos(phi),
                        r * np.sin(theta) * np.sin(phi),
                        r * np.cos(theta)]).transpose() + soma_centre

        # Check which points are inside this hyper voxel
        vox_idx = np.floor((xyz - self.hyper_voxel_origo) / self.voxel_size).astype(int)

        inside_idx = np.where(np.sum(np.bitwise_and(0 <= vox_idx, vox_idx < self.hyper_voxel_size), axis=1) == 3)[0]

        return vox_idx[inside_idx, :], r[inside_idx]

    ############################################################################

    def get_hyper_voxel_axon_points(self,
                                    neuron_position,
                                    rotation,
                                    axon_density_bounds_xyz,
                                    num_points=1000):

        """
        Helper function to give points inside axon bounding box, that are inside hyper voxel

        Args:
            neuron_position (float,float,float): coordinates of neuron
            rotation: rotation matrix
            axon_density_bounds_xyz: boundary box for axon
            num_points: number of points to place
        """

        # Randomly place nPoints inside bounding box (SWC coordinates, soma (0,0,0))
        x_min = axon_density_bounds_xyz[0]
        x_width = axon_density_bounds_xyz[1] - axon_density_bounds_xyz[0]
        y_min = axon_density_bounds_xyz[2]
        y_width = axon_density_bounds_xyz[3] - axon_density_bounds_xyz[2]
        z_min = axon_density_bounds_xyz[4]
        z_width = axon_density_bounds_xyz[5] - axon_density_bounds_xyz[4]

        xyz = self.hyper_voxel_rng.random((num_points, 3))
        xyz[:, 0] = x_min + x_width * xyz[:, 0]
        xyz[:, 1] = y_min + y_width * xyz[:, 1]
        xyz[:, 2] = z_min + z_width * xyz[:, 2]

        # Check which of the points are inside hyper voxel (rotate+translate)
        vox_idx = np.floor((np.matmul(rotation, xyz.transpose()).transpose()
                           + neuron_position - self.hyper_voxel_origo)
                           / self.voxel_size).astype(int)

        inside_idx = np.where(np.sum(np.bitwise_and(0 <= vox_idx, vox_idx < self.hyper_voxel_size), axis=1) == 3)[0]

        return xyz[inside_idx, :], vox_idx[inside_idx, :]

    ############################################################################

    # somaCentre and rotation of neuron
    # axon_density_func should be written so that it can handle x,y,z (SWC coords)
    # as vectors
    # axon_density_bounds_xyz = [xmin,xmax,ymin,ymax,zmin,zmax] in SWC coordinates

    # axon_density_func = eval("lambda x,y,z: " + axonPstr)

    def no_axon_points_xyz(self, neuron_position, rotation,
                           axon_density_func, axon_density_bounds_xyz):

        """
        Placing fake axon segments based on probability distribution speicified with x,y,z.

        Args:
             neuron_position (float,float,float): location of neuron soma
             rotation (3x3 rotation matrix): rotation of neuron
             axon_density_func: axon density function in x,y,z coordinates (SWC coords), must handle x,y,z as vectors
                                e.g. axon_density_func = eval("lambda x,y,z: " + axonPstr)
             axon_density_bounds_xyz: [xmin,xmax,ymin,ymax,zmin,zmax] using the coordinates in the SWC file
        """

        # Points for initial sample
        n_points = 5000

        (xyz_inside, voxIdx) = self.get_hyper_voxel_axon_points(neuron_position,
                                                                rotation,
                                                                axon_density_bounds_xyz,
                                                                n_points)

        x_width = axon_density_bounds_xyz[1] - axon_density_bounds_xyz[0]
        y_width = axon_density_bounds_xyz[3] - axon_density_bounds_xyz[2]
        z_width = axon_density_bounds_xyz[5] - axon_density_bounds_xyz[4]

        point_volume = x_width * y_width * z_width / n_points
        voxel_volume = self.voxel_size ** 3

        # If no points were inside HV, then the intersection must be small
        # so we assume no voxels should be filled
        if xyz_inside.shape[0] == 0:
            # Double check correct data-types
            self.write_log("Bounding box appears to be outside hyper voxel")
            return np.zeros((0, 3), dtype=int), np.zeros((0, 1))

        # Calculate density at each of the points inside HV
        x, y, z = xyz_inside[:, 0], xyz_inside[:, 1], xyz_inside[:, 2]
        density_inside = numexpr.evaluate(axon_density_func)
        # OLD: density_inside = axon_density_func(xyz_inside[:, 0], xyz_inside[:, 1], xyz_inside[:, 2])

        # Estimate number of synapses from density, in this case we use a volume
        # equal to bounding box volume / nPoints for each point.
        # OBS that we only want to know how many synapses to place inside HV
        n_points_to_place = np.round(np.sum(density_inside * point_volume))

        if n_points_to_place <= 0:
            # To little probability mass inside
            self.write_log("Too little of axon appears to be inside hyper voxel")
            return np.zeros((0, 3), dtype=int), np.zeros((0, 1))

        # Calculate max density inside HV, divide all densities by that to
        # get Pkeep.
        max_density = np.max(density_inside)

        # We know that n out of N points placed were inside volume, so volume
        # acceptance rate is n/N.
        # In order to get about nPointsToPlace points placed we need to account
        # for how many outside volume, and also account for how many of those
        # inside volume gets accepted (Pkeep = density / maxDensity)
        n_tries = np.round(n_points_to_place * n_points
                           / np.sum(density_inside / max_density)).astype(int)

        if n_tries > 1e6:
            self.write_log(f"!!! no_axon_points_xyz: Warning trying to place {n_tries} points. "
                           "Bounds: {axon_density_bounds_xyz}")

        # Only print this in verbose mode
        if self.verbose:
            self.write_log(f"Trying to place {n_tries} points to get {n_points_to_place} axon voxels filled")

        if n_points >= n_tries:
            # We already have enough points, use a subset
            use_num = np.round(voxIdx.shape[0] * n_tries / n_points).astype(int)
            picked_idx = np.where(self.hyper_voxel_rng.random(use_num)
                                  < density_inside[:use_num] / max_density)[0]
            axon_dist = np.sqrt(np.sum((xyz_inside[picked_idx, :]) ** 2, axis=1))

            return voxIdx[picked_idx, :], axon_dist
        else:
            # Not enough points, use the ones we have, then get more
            picked_idx = np.where(self.hyper_voxel_rng.random(voxIdx.shape[0])
                                  < density_inside / max_density)[0]
            axon_dist = np.sqrt(np.sum((xyz_inside[picked_idx, :]) ** 2, axis=1))

            # Get more points

            (xyz_inside_b, voxIdxB) = \
                self.get_hyper_voxel_axon_points(neuron_position,
                                                 rotation,
                                                 axon_density_bounds_xyz,
                                                 n_tries - n_points)

            # OLD: density_inside_b = axon_density_func(xyz_inside_b[:, 0], xyz_inside_b[:, 1], xyz_inside_b[:, 2])

            # x,y,z used in axon_density_func below
            x, y, z = xyz_inside_b[:, 0], xyz_inside_b[:, 1], xyz_inside_b[:, 2]
            density_inside_b = numexpr.evaluate(axon_density_func)
            picked_idx_b = np.where(self.hyper_voxel_rng.random(voxIdxB.shape[0]) < density_inside_b / max_density)[0]
            axon_dist_b = np.sqrt(np.sum((xyz_inside_b[picked_idx_b, :]) ** 2, axis=1))

            return (np.concatenate([voxIdx[picked_idx, :],
                                    voxIdxB[picked_idx_b, :]]),
                    np.concatenate([axon_dist, axon_dist_b]))

    ############################################################################

    def resize_hyper_voxel_synapses_matrix(self, new_size=None):

        """
        Increase the maximal size of the synapse matrix used for the hypervoxel.

        Args:
            new_size (int): Number of rows in synapse matrix
        """

        if new_size is None:
            new_size = int(np.ceil(1.5 * self.max_synapses))

        assert new_size >= self.hyper_voxel_synapse_ctr, " Cannot shrink below existing number of synapses"

        # We need to increase matrix size
        old = self.hyper_voxel_synapses
        self.max_synapses = new_size
        self.write_log(f"Increasing max synapses to {self.max_synapses}")
        self.hyper_voxel_synapses = np.zeros((self.max_synapses, 13), dtype=np.int32)
        self.hyper_voxel_synapses[:old.shape[0], :] = old
        del old

    ############################################################################

    # This truncates and sorts the hyper voxel synapse matrix

    def sort_synapses(self):

        """
        Sort synapses stored in self.hyper_voxel_synapses.
        New sort order is columns 1 (dest), 0 (src), 6 (synapse type).
        """

        sort_idx = np.lexsort(self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr,
                              [6, 0, 1]].transpose())  # Sort order: columns 1 (dest), 0 (src), 6 (synapse type)

        self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr, :] = \
            self.hyper_voxel_synapses[sort_idx, :]

    ############################################################################

    def sort_gap_junctions(self):

        """
        Sort gap junctions in self.hyper_voxel_gap_junctions.
        New sort order is columns 1 (dest), 0 (src).
        """

        sort_idx = \
            np.lexsort(self.hyper_voxel_gap_junctions[:self.hyper_voxel_gap_junction_ctr, [0, 1]].transpose())

        self.hyper_voxel_gap_junctions[:self.hyper_voxel_gap_junction_ctr, :] = \
            self.hyper_voxel_gap_junctions[sort_idx, :]

    ############################################################################

    # First and second column specify the source and destination ID of a synapse
    # or a gap junction.
    #
    # This creates a lookup table where all synapses in the hyper voxel
    # between the same pair of neurons are grouped together
    # returns a matrix where first column is a UID = srcID*nNeurons + destID
    # and the following two columns are start row and end row (-1) in matrix

    # Turned the method into static so project.py can use it also

    @staticmethod
    def create_lookup_table(data, n_rows, data_type, num_neurons, max_synapse_type):

        """
        This creates a lookup table where all synapses in the hyper voxel
        between the same pair of neurons are grouped together in the synapse matrix.
        Returns a matrix where first column is a UID = srcID*nNeurons + destID
        and the following two columns are start row and end row (-1) in matrix

        Args:
            data : either synapse matrix, or gap junction matrix
            n_rows : number of rows in matrix that are used (matrix itself can be larger)
            data_type : "synapses" or "gap_junctions"
            num_neurons : number of neurons
            max_synapse_type : the synapse types are numbered, this number must not be too small.


        Returns a matrix where first column is a UID = src_ID*num_neurons + dest_ID
        and the following two columns are start row and end row (-1) in matrix
        """

        # self.write_log("Create lookup table")
        # nRows = data.shape[0] -- zero padded, cant use shape
        lookup_table = np.zeros((data.shape[0], 3), dtype=np.int64)  # Linux is 64-bit int, Windows was 32-bit... :-/

        next_idx = 0
        start_idx = 0

        lookup_idx = 0
        # num_neurons = len(self.neurons)

        if data_type == "synapses":
            hardcoded_synapse_type = None
        elif data_type == "gap_junctions":
            hardcoded_synapse_type = 3  # Hardcoded for gap junctions
        else:
            assert False, f"Unknown data_type {data_type}, should be 'synapses' or ' gap_junctions'"

        # max_synapse_type = self.next_channel_model_id   # This needs to be saved in HDF5 file
        max_synapse_type = np.int64(max_synapse_type)

        while next_idx < n_rows:
            src_id = data[next_idx, 0]
            dest_id = data[next_idx, 1]

            if hardcoded_synapse_type:
                synapse_type = hardcoded_synapse_type
            else:
                synapse_type = data[next_idx, 6]

            next_idx += 1

            while (next_idx < n_rows
                   and data[next_idx, 0] == src_id
                   and data[next_idx, 1] == dest_id):
                next_idx += 1

            lookup_table[lookup_idx, :] = [(dest_id * num_neurons + src_id) * max_synapse_type + synapse_type,
                                           start_idx, next_idx]

            start_idx = next_idx
            lookup_idx += 1

        return lookup_table[:lookup_idx, :]

    ############################################################################

    def includes_gap_junctions(self):

        """ Checks if any gap junctions are defined in self.connectivity_distribution. Returns True or False. """

        has_gap_junctions = False

        for key in self.connectivity_distributions:
            if "gap_junction" in self.connectivity_distributions[key]:
                has_gap_junctions = True

        return has_gap_junctions

    ############################################################################

    # Gap junctions are stored in self.hyperVoxelGapJunctions

    def detect_gap_junctions(self):

        """ Helper function, triggers detection of gap junctions. Called by process_hyper_voxel. """

        if not self.includes_gap_junctions():
            self.write_log("detect_gap_junctions: No gap junctions defined in connectivity rules")
            return

        start_time = timeit.default_timer()

        assert self.hyper_voxel_gap_junction_ctr == 0 and self.hyper_voxel_gap_junctions is not None, \
            "setup_hyper_voxel must be called before detecting gap junctions"

        [x_dv, y_dv, z_dv] = np.where(self.dend_voxel_ctr > 0)

        for x, y, z in zip(x_dv, y_dv, z_dv):

            neuron_id_list = self.dend_voxels[x, y, z, :self.dend_voxel_ctr[x, y, z]]
            seg_id = self.dend_sec_id[x, y, z, :self.dend_voxel_ctr[x, y, z]]
            seg_x = self.dend_sec_x[x, y, z, :self.dend_voxel_ctr[x, y, z]]

            # All possible pairs
            for pairs in itertools.combinations(np.arange(0, self.dend_voxel_ctr[x, y, z]), 2):
                neuron_id1 = self.dend_voxels[x, y, z, pairs[0]]
                neuron_id2 = self.dend_voxels[x, y, z, pairs[1]]

                # !!! Check no self connections??

                # Add type field, derived from name field MSD1_45 --> MSD1
                pre_type = self.neurons[neuron_id1]["type"]
                post_type = self.neurons[neuron_id2]["type"]

                if (pre_type, post_type) in self.connectivity_distributions:

                    if "gap_junction" in self.connectivity_distributions[pre_type, post_type]:
                        con_info = self.connectivity_distributions[pre_type, post_type]["gap_junction"]

                        seg_id1 = self.dend_sec_id[x, y, z, pairs[0]]
                        seg_id2 = self.dend_sec_id[x, y, z, pairs[1]]

                        seg_x1 = self.dend_sec_x[x, y, z, pairs[0]]
                        seg_x2 = self.dend_sec_x[x, y, z, pairs[1]]

                        mean_gj_cond, std_gj_cond = con_info["conductance"]
                        gj_mu, gj_sigma = con_info["lognormal_mu_sigma"]

                        # !!! Currently not using channelParamDict for GJ

                        # gj_cond = self.hyper_voxel_rng.normal(mean_gj_cond, std_gj_cond)
                        # lognormal distribution https://www.nature.com/articles/nrn3687
                        gj_cond = self.hyper_voxel_rng.lognormal(gj_mu, gj_sigma)

                        gj_cond = np.maximum(gj_cond, mean_gj_cond * 0.1)  # Avoid negative cond

                        self.hyper_voxel_gap_junctions[self.hyper_voxel_gap_junction_ctr, :] = \
                            [neuron_id1, neuron_id2, seg_id1, seg_id2, seg_x1 * 1e3, seg_x2 * 1e3,
                             x, y, z, self.hyper_voxel_id, gj_cond * 1e12]
                        self.hyper_voxel_gap_junction_ctr += 1

        self.sort_gap_junctions()

        # We also translate from local hyper voxel coordinates to simulation
        # voxel coordinates

        hyper_voxel_offset = np.round((self.hyper_voxel_origo - self.simulation_origo)
                                      / self.hyper_voxel_width).astype(int) * self.hyper_voxel_size

        self.hyper_voxel_gap_junctions[:self.hyper_voxel_gap_junction_ctr, :][:, range(6, 9)] += hyper_voxel_offset

        self.hyper_voxel_gap_junction_lookup = self.create_lookup_table(data=self.hyper_voxel_gap_junctions,
                                                                        n_rows=self.hyper_voxel_gap_junction_ctr,
                                                                        data_type="gap_junctions",
                                                                        num_neurons=len(self.neurons),
                                                                        max_synapse_type=self.next_channel_model_id)
        end_time = timeit.default_timer()

        self.write_log(f"detect_gap_junctions: {end_time - start_time:.1f} s")

        return self.hyper_voxel_gap_junctions[:self.hyper_voxel_gap_junction_ctr, :]

    ############################################################################

    def setup_log(self, logfile_name=None):

        """
        Initiates log file.

        Args:
            logfile_name (str) : Path to log file

        """

        if logfile_name is None:
            logfile_name = self.logfile_name

        if logfile_name is None or len(logfile_name) == 0:
            # Not a valid log file name
            return

        if self.logfile is not None:
            self.write_log("Already have a log file setup, ignoring")
            return

        # If directory does not exist, create
        dir_name = os.path.dirname(logfile_name)
        if not os.path.exists(dir_name):
            self.write_log(f"Creating directory {dir_name}")
            os.makedirs(dir_name)

        self.logfile = open(logfile_name, 'wt')

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
            self.logfile.write(f"{text}\n")
            if flush:
                self.logfile.flush()

        if self.verbose or is_error or force_print:
            print(text, flush=True)

    ############################################################################

    def read_prototypes(self, config_file=None):

        """
        Read in neuron prototypes. A neuron prototype can have multiple parameters, and morphology variations.

        Args:
            config_file (str): path to network config file
        """

        if config_file is None:
            config_file = self.config_file

        if config_file is None:
            self.write_log(f"No config file specified. Not reading any neuron prototypes.")
            return

        self.write_log(f"Loading from {config_file}")

        cfg_file = open(str(config_file), 'r')

        try:
            self.config = json.load(cfg_file)
        finally:
            cfg_file.close()

        # This also loads random seed from config file while we have it open
        if self.random_seed is None:
            if "random_seed" in self.config and "detect" in self.config["random_seed"]:
                self.random_seed = self.config["random_seed"]["detect"]
                self.write_log(f"Reading random seed from config file: {self.random_seed}")
            else:
                # No random seed given, invent one
                self.random_seed = 1002
                self.write_log(f"No random seed provided, using: {self.random_seed}")
        else:
            self.write_log(f"Using random seed provided by command line: {self.random_seed}")

        self.prototype_neurons = dict()

        for region_name, region_data in self.config["regions"].items():
            for name_type, definition in region_data["neurons"].items():
                for name, neuron_path in definition["neuron_path"].items():

                    self.write_log(f"Reading prototype for: {name}")

                    if "neuron_type" in definition:
                        neuron_type = definition["neuron_type"]
                    else:
                        neuron_type = "neuron"

                    if neuron_type == "virtual":
                        virtual_neuron = True
                    else:
                        virtual_neuron = False

                    if 'hoc' in definition:
                        hoc = definition["hoc"]
                        assert "hoc no longer passed to NeuronPrototype / NeuronMorphology -- need to add it later "
                    else:
                        hoc = None

                    self.prototype_neurons[name] = NeuronPrototype(neuron_name=name,
                                                                   neuron_path=neuron_path,
                                                                   snudda_data=self.snudda_data,
                                                                   # hoc=hoc,
                                                                   virtual_neuron=virtual_neuron)

                    if "axon_density" in definition:

                        # We need to do this so that we can apply the axon densities below
                        self.prototype_neurons[name].instantiate()

                        self.write_log("Setting axon density for neuron without axon")
                        axon_density_type = definition["axon_density"][0]

                        if axon_density_type == "r":
                            density = definition["axon_density"][1]
                            max_radius = definition["axon_density"][2]

                            self.prototype_neurons[name].apply("set_axon_voxel_radial_density", [density, max_radius])
                        elif axon_density_type == "xyz":
                            density = definition["axon_density"][1]
                            axon_density_bounds_xyz = np.array(definition["axon_density"][2])

                            self.prototype_neurons[name].apply("set_axon_voxel_xyz_density", [density, axon_density_bounds_xyz])

                        else:
                            self.write_log(f"{name}: Unknown axon density type : {axon_density_type}\n"
                                           f"{definition['axon_density']}", is_error=True)

                    else:
                        # If no axon density specified, then axon must be present in morphology
                        assert self.prototype_neurons[name].all_have_axon(), f"A morphlogy in {neuron_path} does not have an axon"

                    assert self.prototype_neurons[name].all_have_dend() or self.prototype_neurons[name].virtual_neuron, \
                        f"A morphology in {neuron_path} does not have a dendrite"

                    # Since we already have the config file open, let's read connectivity
                    # distributions also

        self.write_log("Loading connectivity information")
        self.next_channel_model_id = 10  # Reset counter

        for region_name, region_data in self.config["regions"].items():
            for name, connection_def in region_data["connectivity"].items():

                # This also enriches the self.config by adding channelModelID, lognormal_mu_sigma etc
                con_def = copy.deepcopy(connection_def)

                pre_type, post_type = name.split(",")

                for key in con_def:
                    if key == "gap_junction":
                        con_def[key]["channel_model_id"] = 3
                    else:
                        con_def[key]["channel_model_id"] = self.next_channel_model_id
                        self.next_channel_model_id += 1

                    # Also if conductance is just a number, add std 0
                    if type(con_def[key]["conductance"]) not in [list, tuple]:
                        con_def[key]["conductance"] = [con_def[key]["conductance"], 0]

                    # Precompute lognormal parameters
                    # https://en.wikipedia.org/wiki/Log-normal_distribution
                    mean_cond = con_def[key]["conductance"][0]
                    std_cond = con_def[key]["conductance"][1]
                    mu = np.log(mean_cond ** 2 / np.sqrt(mean_cond ** 2 + std_cond ** 2))
                    sigma = np.sqrt(np.log(1 + std_cond ** 2 / mean_cond ** 2))
                    con_def[key]["lognormal_mu_sigma"] = [mu, sigma]

                    if "cluster_size" not in con_def[key]:
                        con_def[key]["cluster_size"] = 1

                    if "cluster_spread" not in con_def[key]:
                        con_def[key]["cluster_spread"] = 20e-3

                    if "RxD" in con_def[key] and "weight_scale" not in con_def:
                        print(f"Connection {key} uses RxD, but does not specify weight_scale set, will use default scaling 1.")

                self.connectivity_distributions[pre_type, post_type] = con_def

    ############################################################################

    def read_neuron_positions(self, position_file):

        """
        Loads neuron positions from network's position_file.

        Args:
            position_file : path to network position file (network-neuron-positions.hdf5)
        """

        if position_file is None:
            position_file = self.position_file

        if position_file is None:
            self.write_log("No position file specified, not reading neuron positions")
            return

        mem = self.memory()
        self.write_log(f"{mem}")

        self.write_log(f"Reading positions from file: {position_file}")

        pos_info = SnuddaLoad(position_file).data

        mem = self.memory()
        self.write_log(f"{mem}")

        # Make sure we do not change config file unintentionally
        assert pos_info["config_file"] == self.config_file, \
            f"Not using original config file: {pos_info['config_file']} \nvs\n{self.config_file}"

        self.neurons = pos_info["neurons"]
        num_neurons = len(self.neurons)

        self.neuron_positions = np.zeros((num_neurons, 3))

        for ni, neuron in enumerate(pos_info["neurons"]):
            self.neuron_positions[ni, :] = neuron["position"]

            # Add a few sanity checks
            assert ni == neuron["neuron_id"], f"neuron_id={neuron['neuron_id']} and ni={ni} not equal, corruption?"
            assert neuron["name"] in self.prototype_neurons, \
                f"Neuron type {neuron['name']} not in prototype_neurons: {self.prototype_neurons}"

        # Also load population_unit data
        self.population_unit = pos_info["population_unit"]

        self.write_log("Position file read.")
        del pos_info

    ############################################################################

    # If the detect is rerun we need to make sure there are not old MERGE
    # files left that might remember old run accidentally

    def delete_old_merge(self):

        """ Cleans up data files from previous detection run. """

        if self.role == "master" and self.network_path is not None:
            del_files = [os.path.join(self.network_path, "network-putative-synapses-MERGED.hdf5"),
                         os.path.join(self.network_path, "network-putative-synapses-MERGED.hdf5-cache"),
                         os.path.join(self.network_path, "network-synapses.hdf5"),
                         os.path.join(self.network_path, "network-synapses.hdf5-cache"),
                         os.path.join(self.network_path, "network-projection-synapses.hdf5"),
                         os.path.join(self.network_path, "pruning_merge_info.json"),
                         os.path.join(self.network_path, "input-spikes.hdf5")
                         ]

            for f in del_files:
                if os.path.exists(f):
                    self.write_log(f"Removing old files {f}")
                    os.remove(f)

    ############################################################################

    def write_hyper_voxel_to_hdf5(self):

        """ Saves hyper voxel synapses to data file. """

        start_time = timeit.default_timer()

        output_name = self.save_file.replace(".hdf5", f"-{self.hyper_voxel_id}.hdf5")

        with h5py.File(output_name, "w", libver=self.h5libver) as out_file:

            out_file.create_dataset("config", data=json.dumps(self.config))

            meta_data = out_file.create_group("meta")
            meta_data.create_dataset("hyper_voxel_id", data=self.hyper_voxel_id)
            meta_data.create_dataset("hyper_voxel_origo", data=self.hyper_voxel_origo)
            meta_data.create_dataset("simulation_origo", data=self.simulation_origo)

            meta_data.create_dataset("snudda_data", data=self.snudda_data)
            meta_data.create_dataset("slurm_id", data=self.slurm_id)
            meta_data.create_dataset("voxel_size", data=self.voxel_size)
            meta_data.create_dataset("hyper_voxel_size", data=self.hyper_voxel_size)
            meta_data.create_dataset("num_hyper_voxels", data=self.num_hyper_voxels)

            meta_data.create_dataset("num_bins", data=self.num_bins)
            meta_data.create_dataset("voxel_overflow_counter", data=self.voxel_overflow_counter)

            meta_data.create_dataset("config_file", data=self.config_file)
            meta_data.create_dataset("position_file", data=self.position_file)

            # These may or may not exist, if they do, write them to file
            if self.max_axon_voxel_ctr is not None:
                meta_data.create_dataset("max_axon_voxel_ctr", data=self.max_axon_voxel_ctr)

            if self.max_dend_voxel_ctr is not None:
                meta_data.create_dataset("max_dend_voxel_ctr", data=self.max_dend_voxel_ctr)

            if self.voxel_overflow_counter > 0:
                self.write_log("!!! Voxel overflow detected, please increase max_axon and max_dend", is_error=True)

            network_group = out_file.create_group("network")
            network_group.create_dataset("synapses",
                                         data=self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr, :],
                                         dtype=np.int32,
                                         chunks=(self.synapse_chunk_size, 13),
                                         maxshape=(None, 13),
                                         compression=self.h5compression)
            network_group.create_dataset("gap_junctions",
                                         data=self.hyper_voxel_gap_junctions[:self.hyper_voxel_gap_junction_ctr, :],
                                         dtype=np.int32,
                                         chunks=(self.gap_junction_chunk_size, 11),
                                         maxshape=(None, 11),
                                         compression=self.h5compression)

            network_group.create_dataset("synapse_lookup",
                                         data=self.hyper_voxel_synapse_lookup,
                                         dtype=np.int64)   # Need to force 64-bit int for Windows

            network_group.create_dataset("gap_junction_lookup",
                                         data=self.hyper_voxel_gap_junction_lookup,
                                         dtype=np.int64)  # Need to force 64-bit int for Windows

            network_group.create_dataset("max_channel_type_id", data=self.next_channel_model_id, dtype=int)

            # Additional information useful for debugging
            if self.debug_flag:
                debug_group = out_file.create_group("debug")

                debug_group.create_dataset("dend_voxels", data=self.dend_voxels)
                debug_group.create_dataset("axon_voxels", data=self.axon_voxels)

                debug_group.create_dataset("dend_voxel_ctr", data=self.dend_voxel_ctr)
                debug_group.create_dataset("axon_voxel_ctr", data=self.axon_voxel_ctr)

            end_time = timeit.default_timer()

            out_file.close()

        self.write_log(f"Wrote hyper voxel {self.hyper_voxel_id}"
                       f" ({self.hyper_voxel_synapse_ctr} synapses, "
                       f"{self.hyper_voxel_gap_junction_ctr} gap junctions)")

    ############################################################################

    def load_neuron(self, neuron_info, use_cache=True):

        """
        Load neuron.

        Args:
            neuron_info : dictionary with neuron information, i.e. 'name', 'parameterID', 'morphologyID',
                          'modulationID', 'rotation', 'position'
        """

        neuron_id = neuron_info["neuron_id"]

        if use_cache and neuron_id in self.neuron_cache:
            return self.neuron_cache[neuron_id]

        morph_path = snudda_parse_path(neuron_info["morphology"], self.snudda_data)
        if os.path.isfile(morph_path):
            morphology_path = morph_path
        else:
            morphology_path = None  # Get morpholog automatically from morphology_key

        # Clone prototype neuron (it is centred, and not rotated)
        neuron = self.prototype_neurons[neuron_info["name"]].clone(parameter_key=neuron_info["parameter_key"],
                                                                   morphology_key=neuron_info["morphology_key"],
                                                                   modulation_key=neuron_info["modulation_key"],
                                                                   rotation=neuron_info["rotation"],
                                                                   position=neuron_info["position"],
                                                                   morphology_path=morphology_path)

        if "axon_density" in neuron_info and neuron_info["axon_density"] is not None:
            if neuron_info["axon_density_type"] == "r":
                neuron.set_axon_voxel_radial_density(neuron_info["axon_density"],
                                                     neuron_info["axon_density_radius"])
            elif neuron_info["axon_density_type"] == "xyz":
                neuron.set_axon_voxel_xyz_density(neuron_info["axon_density"],
                                                  neuron_info["axon_density_bounds_xyz"])

            else:
                raise ValueError(f"Unknown axon density: {neuron_info['axon_density_type']}")

        if "extra_axons" in neuron_info:
            for axon_name, axon_info in neuron_info["extra_axons"].items():
                # print(f"DEVELOPMENT: Adding extra axon to {neuron_info['name']}: {axon_name} morph {axon_info['morphology']}, pos {axon_info['position']}, rot {axon_info['rotation']}")

                if axon_info["morphology"] not in self.extra_axon_cache:
                    self.extra_axon_cache[axon_info["morphology"]] = MorphologyData(swc_file=axon_info["morphology"],
                                                                                    parent_tree_info=None,
                                                                                    snudda_data=self.snudda_data)

                neuron.add_morphology(swc_file=axon_info["morphology"],
                                      name=axon_name,
                                      position=axon_info["position"],
                                      rotation=axon_info["rotation"],
                                      morphology_data=self.extra_axon_cache[axon_info["morphology"]])

        if use_cache:
            self.neuron_cache[neuron_id] = neuron

        return neuron

    ############################################################################

    def distribute_neurons_parallel(self, d_view=None):

        """ Locates which hyper voxel each neuron is present in."""

        if self.role != "master":
            # Only run this as master
            return

        (hyper_voxels, hyper_voxel_id_lookup, num_hyper_voxels, simulation_origo) = \
            self.get_neuron_distribution_history()

        # Do we have old data that we can reuse?
        if hyper_voxels is not None:
            self.write_log("distribute_neurons_parallel: Reusing old neuron allocation")

            self.hyper_voxels = hyper_voxels
            self.hyper_voxel_id_lookup = hyper_voxel_id_lookup
            self.num_hyper_voxels = num_hyper_voxels
            self.simulation_origo = simulation_origo

            if d_view:
                # We need to push the data to the workers also
                d_view.push({"sd.simulation_origo": simulation_origo,
                             "sd.hyper_voxels": hyper_voxels,
                             "sd.hyper_voxel_id_lookup": hyper_voxel_id_lookup,
                             "sd.num_hyper_voxels": num_hyper_voxels}, block=True)
            return

        # No old data, we need to calculate it

        distribution_seeds = self.generate_neuron_distribution_random_seeds()

        if d_view is None:
            self.write_log("No d_view specified, running distribute neurons in serial", force_print=True)
            (min_coord, max_coord) = self.distribute_neurons(distribution_seeds=distribution_seeds)

            self.count_and_sort_neurons_in_hypervoxels()
            self.generate_hyper_voxel_random_seeds()

            self.save_neuron_distribution_history(hyper_voxels=self.hyper_voxels,
                                                  min_coord=min_coord,
                                                  max_coord=max_coord)

            return

        (min_coord, max_coord) = self.find_min_max_coord_parallel(d_view=d_view,
                                                                  volume_id=self.volume_id)

        # The order here should not affect reproducibility, each neuron has its own seed for distribution part
        # but only those with probabilistic axon clouds will use it.
        neuron_idx = np.random.permutation(np.arange(0, len(self.neurons), dtype=np.int32))

        # Split the neuron_idx between the workers
        d_view.scatter("neuron_idx", neuron_idx, block=True)
        d_view.scatter("distribution_seeds", distribution_seeds[neuron_idx], block=True)  # Need to preserve order
        d_view.push({"min_coord": min_coord, "max_coord": max_coord}, block=True)

        self.write_log("Distributing neurons, parallel.")

        # For the master node, run with empty list
        # This sets up internal state of master
        self.distribute_neurons(neuron_idx=[], min_coord=min_coord, max_coord=max_coord, distribution_seeds=[])

        cmd_str = ("sd.distribute_neurons(neuron_idx=neuron_idx, distribution_seeds=distribution_seeds, "
                   "min_coord=min_coord, max_coord=max_coord)")
        d_view.execute(cmd_str, block=True)

        self.write_log("Gathering neuron distribution from workers")

        # Collect all the neurons in the list from the workers
        # For each neuron we found out which hyper voxels it occupies,
        # now we want for each hyper voxel to know which neurons are in there
        hyper_voxel_list = d_view.gather("sd.hyper_voxels", block=True)

        self.write_log("Distributions received.")

        for hv in hyper_voxel_list:
            for hid in hv:

                assert (hv[hid]["origo"] == self.hyper_voxels[hid]["origo"]).all(), \
                    "Origo for hyper voxels do not match --- should never happen"

                for neuron_id in hv[hid]["neurons"]:
                    assert neuron_id not in self.hyper_voxels[hid]["neurons"], \
                        f"Internal error, neuron_id {neuron_id } already exists hyper_voxel {hid}"

                    self.hyper_voxels[hid]["neurons"][neuron_id] = hv[hid]["neurons"][neuron_id]

                for neuron_id in hv[hid]["axon_density"]:
                    assert neuron_id not in self.hyper_voxels[hid]["axon_density"], \
                        f"Internal error, neuron_id {neuron_id} already exists hyper_voxel {hid} (axon_density)"
                    self.hyper_voxels[hid]["axon_density"].append(neuron_id)

        # Sort for reproducibility
        self.count_and_sort_neurons_in_hypervoxels()
        self.generate_hyper_voxel_random_seeds()

        # Distribute the new list to all neurons
        d_view.push({"sd.hyper_voxels": self.hyper_voxels}, block=True)

        self.save_neuron_distribution_history(hyper_voxels=self.hyper_voxels,
                                              min_coord=min_coord,
                                              max_coord=max_coord)

    ############################################################################

    # This creates a list for each hyper voxel for the neurons that
    # has any neurites within its border (here defined as vertices inside region)

    def distribute_neurons(self, neuron_idx=None, distribution_seeds=None, min_coord=None, max_coord=None):

        """
        This creates a list for each hyper voxel of the neurons that
        has any neurites within its border (here defined as vertices inside region)

        Args:
            neuron_idx : neuron_id:s to process
            distribution_seeds : Random seed (used for neurons without axon)
            min_coord (float, float, float) : Minimum x,y,z coordinates
            max_coord (float, float, float) : Maximum x,y,z coordinates

        Updates self.hyper_voxels. Also returns min_coord, max_coord
        """

        try:

            if neuron_idx is None:
                neuron_idx = np.arange(0, len(self.neurons), dtype=np.int32)

            assert distribution_seeds is not None and len(neuron_idx) == len(distribution_seeds), \
                "distribute_neurons - distribution seeds needed for reproducability"

            self.write_log(f"distribute_neurons: neuron_idx = {neuron_idx} (n={len(neuron_idx)})")
            start_time = timeit.default_timer()

            if max_coord is None or min_coord is None:
                self.write_log("distribute_neurons: calculating min and max coords")
                (min_coord, max_coord) = self.find_min_max_coord()

            # Simulation origo is in meters
            if self.simulation_origo is None:
                # We align the simulation origo to the closest voxel (that is smaller)
                self.simulation_origo = np.floor(min_coord / self.voxel_size) * self.voxel_size
            elif (self.simulation_origo > min_coord).any():
                raise ValueError(f"Simulation origo ({self.simulation_origo}) must be smaller than {min_coord}. "
                                 f"This since all voxel and hyper voxel coordinates must be positive.")

            assert ((self.num_bins - self.num_bins[0]) == 0).all(), "Hyper voxels should be cubes"

            self.setup_hyper_voxel_id_lookup(max_coord=max_coord, min_coord=min_coord)
            self.preallocate_empty_hyper_voxel_dict()
            ctr = 0

            if neuron_idx is None:
                neurons = self.neurons
            elif len(neuron_idx) == 0:
                neurons = []
            else:
                neurons = [self.neurons[idx] for idx in neuron_idx]

            for n, d_seed in zip(neurons, distribution_seeds):

                ctr = ctr + 1
                if ctr % 10000 == 0:
                    self.write_log(f"Assignment counter: {ctr}")

                neuron = self.load_neuron(n, use_cache=False)
                neuron_id = n["neuron_id"]

                tree_info = self.get_hypervoxel_coords_and_section_id(neuron=neuron)
                section_info = self.group_section_info(tree_info=tree_info)
                density_hyper_voxel_id = self.get_density_location(neuron=neuron, seed=d_seed)

                # First loop over section info, add info
                for h_id in section_info:

                    neuron_data = dict()

                    if 1 in section_info[h_id] and "neuron" in section_info[h_id][1]:
                        neuron_data["soma"] = True

                    if 2 in section_info[h_id]:
                        for sub_tree_name in section_info[h_id][2]:
                            sec_id_list = section_info[h_id][2][sub_tree_name]

                            if "axon" not in neuron_data:
                                neuron_data["axon"] = [(sec_id_list, sub_tree_name)]
                            else:
                                neuron_data["axon"].append((sec_id_list, sub_tree_name))

                    if 3 in section_info[h_id] and "neuron" in section_info[h_id][3]:
                        sec_id_list = section_info[h_id][3]["neuron"]
                        neuron_data["dend"] = sec_id_list

                    if len(neuron_data) > 0:
                        self.hyper_voxels[h_id]["neurons"][int(neuron_id)] = neuron_data

                # Then loop over density info, add data
                for h_id in density_hyper_voxel_id:
                    self.hyper_voxels[h_id]["axon_density"].append(neuron_id)

            end_time = timeit.default_timer()

            if len(neurons) > 0:
                self.write_log(f"Calculated distribution of neurons: {end_time - start_time:.1f} seconds")

        except Exception as e:
            # Write error to log file to help trace it.
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str, is_error=True)

            sys.exit(-1)

        # For serial version of code, we need to return this, so we
        # can save work history
        return min_coord, max_coord

    def count_and_sort_neurons_in_hypervoxels(self):

        # Sorting the list of neurons (needed for reproducibility when axon is probability cloud and we sample them)
        # and converting to numpy arrays

        for hid in self.hyper_voxels:
            self.hyper_voxels[hid]["axon_density"] = np.array(sorted(self.hyper_voxels[hid]["axon_density"]))
            self.hyper_voxels[hid]["neuron_ctr"] = len(self.hyper_voxels[hid]["neurons"])

    ############################################################################

    def setup_parallel(self, d_view=None):

        """ Prepares workers for parallel execution if d_view is not None. """

        assert self.role == "master", \
            "setup_parallel: Should only be called by master node"

        if d_view is None:
            self.write_log("setup_parallel called without dView, aborting.")
            return

        if self.workers_initialised:
            self.write_log("Workers already initialised.")
            return

        self.write_log(f"setup_parallel: {d_view = }")

        with d_view.sync_imports():
            from snudda import SnuddaDetect

        self.write_log(f"Setting up workers: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")

        # Create unique log file names for the workers
        if self.logfile_name is not None:
            engine_logfile = [f"{self.logfile_name}-{x}" for x in range(0, len(d_view))]
        else:
            engine_logfile = [[] for x in range(0, len(d_view))]

        self.write_log(f"Scattering engine_logfile = {engine_logfile}")
        d_view.scatter('logfile_name', engine_logfile, block=True)
        self.write_log("Scatter done.")

        d_view.push({"position_file": self.position_file,
                     "config_file": self.config_file,
                     "snudda_data": self.snudda_data,
                     "voxel_size": self.voxel_size,
                     "hyper_voxel_size": self.hyper_voxel_size,
                     "verbose": self.verbose,
                     "slurm_id": self.slurm_id,
                     "save_file": self.save_file,
                     "random_seed": self.random_seed},
                    block=True)

        self.write_log("Init values pushed to workers")

        cmd_str = ("sd = SnuddaDetect(config_file=config_file, position_file=position_file,voxel_size=voxel_size,"
                   "snudda_data=snudda_data,"
                   "hyper_voxel_size=hyper_voxel_size,verbose=verbose,logfile_name=logfile_name[0],"
                   "save_file=save_file,slurm_id=slurm_id,role='worker', random_seed=random_seed)")

        cmd_str = """
try:
    sd = SnuddaDetect(config_file=config_file, position_file=position_file,voxel_size=voxel_size,
                      snudda_data=snudda_data,
                      hyper_voxel_size=hyper_voxel_size,verbose=verbose,logfile_name=logfile_name[0],
                      save_file=save_file,slurm_id=slurm_id,role='worker', random_seed=random_seed)
except Exception as e:
    import os
    import datetime
    import traceback
    engine_id = os.getpid()
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file_name = f"worker_error_{engine_id}_{timestamp}.log"
    with open(log_file_name, "w") as log_file:
        log_file.write(f"Error initializing SnuddaDetect:\\n{traceback.format_exc()}\\n")
        log_file.write(f"Engine ID: {engine_id}\\n")
        log_file.write(f"Timestamp: {timestamp}\\n")
        log_file.write(f"Parameters: config_file={config_file}, position_file={position_file}, voxel_size={voxel_size}\\n")
    
    raise  # Rethrow exception
        """

        d_view.execute(cmd_str, block=True)

        self.write_log(f"Workers setup: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")
        self.workers_initialised = True

    ############################################################################

    def find_min_max_coord_parallel(self, volume_id=None, d_view=None):

        """
        Finds the minimum and maximum coordinates in entire model for all neuron components.

        Args:
            volume_id : Volume ID to check, None means all
            d_view : Direct view object

        Returns:
            min_coord, max_coord
        """

        if d_view is None:
            self.write_log("find_min_max_coord_parallel: dView is None")
            return self.find_min_max_coord(volume_id=volume_id)

        self.write_log("Finding min/max coords parallel")

        neuron_idx = np.random.permutation(np.arange(0, len(self.neurons), dtype=np.int32))

        d_view.scatter("neuron_idx_find", neuron_idx, block=True)
        d_view.push({"volume_id": volume_id}, block=True)

        cmd_str = "min_max = sd.find_min_max_coord(volume_id=volume_id,neuron_idx=neuron_idx_find)"

        d_view.execute(cmd_str, block=True)

        self.write_log("Execution of min/max complete")
        all_min_max = d_view["min_max"]
        self.write_log("Gathered min/max - complete.")

        max_coord = -1e6 * np.ones((3,))
        min_coord = 1e6 * np.ones((3,))

        for (min_c, max_c) in all_min_max:
            max_coord = np.maximum(max_coord, max_c)
            min_coord = np.minimum(min_coord, min_c)

        return min_coord, max_coord

    ############################################################################

    def find_min_max_coord(self, volume_id=None, neuron_idx=None):

        """
        Finds the minimum and maximum coordinates in entire model for all neuron components.

        Args:
            volume_id : Volume ID to check, None means all

        Returns:
            min_coord, max_coord
        """

        try:
            if volume_id is None:
                volume_id = self.volume_id

            self.write_log(f"Finding minMax coord in volume_id = {volume_id}")

            max_coord = -1e6 * np.ones((3,))
            min_coord = 1e6 * np.ones((3,))

            if neuron_idx is None:
                neurons = self.neurons
            else:
                neurons = [self.neurons[idx] for idx in neuron_idx]

            for n in neurons:

                # By using "in" for comparison, we can pass a list of volume_id also
                if volume_id is not None and n["volume_id"] not in volume_id:
                    self.write_log(f"Skipping {n['name']} when calculating hyper voxel size")
                    # Only include neurons belonging to the volume ID
                    # we are looking at now
                    continue

                neuron = self.load_neuron(n, use_cache=False)
                for subtree in neuron.morphology_data.values():
                    try:
                        max_coord = np.maximum(max_coord, np.max(subtree.geometry[:, :3], axis=0))
                        min_coord = np.minimum(min_coord, np.min(subtree.geometry[:, :3], axis=0))

                        if np.isnan(min_coord).any() or np.isnan(max_coord).any():
                            print(f"Problem with morphology: {neuron}")
                            import pdb
                            pdb.set_trace()
                    except:
                        import traceback
                        print(traceback.format_exc())
                        import pdb
                        pdb.set_trace()

        except Exception as e:
            # Write error to log file to help trace it.
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str, is_error=True)

            sys.exit(-1)

        return min_coord, max_coord

    ############################################################################

    #TODO: We should probably NUMBA this method
    def fill_voxels_soma(self, voxel_space, voxel_space_ctr,
                         voxel_sec_id, voxel_sec_x,
                         neuron, neuron_id, verbose=False):

        """
        Marks all the dendrite voxels that all the somas in the hyper voxel occupy.

        voxel_space : n x n x n x k matrix holding the voxel content, normally self.dend_voxels (neuron IDs)
        voxel_space_ctr : n x n x n matrix holding count of how many items each voxel holds
        voxel_sec_id : n x n x n x k matrix, holding section ID of each item
        voxel_sec_x : n x n x n x k matrix, holding section X of each item
        soma_coord : (x,y,z,r) location of soma to place, and radius
        neuron_id : ID of the neurons
        verbose (bool) : how much to print out

        """

        if "neuron" not in  neuron.morphology_data:
            return  # No soma defined
        soma_coord = neuron.morphology_data["neuron"].geometry[0, :]

        if neuron.morphology_data["neuron"].section_data[0, 2] != 1:
            raise ValueError(f"First compartment should be soma. Neuron ID {neuron_id}")

        try:
            v_coords = np.floor((soma_coord[:3] - self.hyper_voxel_origo) / self.voxel_size).astype(int)
            radius2 = soma_coord[3] ** 2
            v_radius = np.ceil(soma_coord[3] / self.voxel_size).astype(int)
        except:
            import traceback
            traceback.format_exc()
            import pdb
            pdb.set_trace()

        assert v_radius < 1000, \
            f"fill_voxels_soma: v_radius={v_radius} soma coords = {soma_coord} (BIG SOMA, not SI units?)"

        # Range check, so we stay within hypervoxel
        vx_min = max(0, v_coords[0] - v_radius)
        vx_max = min(self.hyper_voxel_size, v_coords[0] + v_radius + 1)

        vy_min = max(0, v_coords[1] - v_radius)
        vy_max = min(self.hyper_voxel_size, v_coords[1] + v_radius + 1)

        vz_min = max(0, v_coords[2] - v_radius)
        vz_max = min(self.hyper_voxel_size, v_coords[2] + v_radius + 1)

        # self.write_log(f"Soma check x: {vx_min} - {vx_max} y: {vy_min} - {vy_max} z: {vz_min} - {vz_max}")

        for vx in range(vx_min, vx_max):
            for vy in range(vy_min, vy_max):
                for vz in range(vz_min, vz_max):

                    d2 = (((vx + 0.5) * self.voxel_size + self.hyper_voxel_origo[0] - soma_coord[0]) ** 2
                          + ((vy + 0.5) * self.voxel_size + self.hyper_voxel_origo[1] - soma_coord[1]) ** 2
                          + ((vz + 0.5) * self.voxel_size + self.hyper_voxel_origo[2] - soma_coord[2]) ** 2)

                    if d2 < radius2:
                        # Mark the point
                        try:
                            v_ctr = voxel_space_ctr[vx, vy, vz]

                            if (v_ctr > 0
                                    and voxel_space[vx, vy, vz, v_ctr - 1] == neuron_id):
                                # Voxel already has neuron_id, skip
                                continue

                            voxel_space[vx, vy, vz, v_ctr] = neuron_id
                            voxel_sec_id[vx, vy, vz, v_ctr] = -1  # Soma is section_id -1
                            voxel_sec_x[vx, vy, vz, v_ctr] = 0.5

                            voxel_space_ctr[vx, vy, vz] += 1
                        except:
                            self.voxel_overflow_counter += 1
                            self.write_log("!!! If you see this you need to increase max_dend above "
                                           f"{voxel_space_ctr[vx, vy, vz]}", is_error=True)
                            continue

    ############################################################################

    def fill_voxels_dend(self, voxel_space, voxel_space_ctr,
                         voxel_sec_id, voxel_sec_x, voxel_soma_dist,
                         neuron, neuron_id, section_id=None):

        """
        Mark all voxels containing dendrites.

        voxel_space : n x n x n x k matrix holding the voxel content, normally self.dend_voxels (neuron IDs)
        voxel_space_ctr : n x n x n matrix holding count of how many items each voxel holds
        voxel_sec_id : n x n x n x k matrix, holding section ID of each item
        voxel_sec_x : n x n x n x k matrix, holding section X of each item
        voxel_soma_dist : n x n x n x k matrix, holding distance to soma along dendrite
        neuron : NeuronMorphologyExtended object
        neuron_id : ID of the neurons
        section_id :

        """

        # Can we move the iterator into numba?
        for section in neuron.section_iterator_selective(section_type=3, section_id=section_id):

            voxel_overflow_ctr = SnuddaDetect.fill_voxels_dend_helper(voxel_space=voxel_space,
                                                                      voxel_space_ctr=voxel_space_ctr,
                                                                      voxel_sec_id=voxel_sec_id,
                                                                      voxel_sec_x=voxel_sec_x,
                                                                      voxel_soma_dist=voxel_soma_dist,
                                                                      point_idx=section.point_idx,
                                                                      geometry=section.morphology_data.geometry,
                                                                      section_data=section.morphology_data.section_data,
                                                                      neuron_id=neuron_id,
                                                                      self_hyper_voxel_origo=self.hyper_voxel_origo,
                                                                      self_voxel_size=self.voxel_size,
                                                                      self_num_bins=self.num_bins,
                                                                      self_max_dend=self.max_dend,
                                                                      self_step_multiplier=self.step_multiplier)
            self.voxel_overflow_counter += voxel_overflow_ctr

    @staticmethod
    @jit(nopython=True, fastmath=True, cache=True)
    def fill_voxels_dend_helper(voxel_space, voxel_space_ctr,
                                voxel_sec_id, voxel_sec_x,
                                voxel_soma_dist,
                                point_idx, geometry, section_data, neuron_id: int,
                                self_hyper_voxel_origo, self_voxel_size,
                                self_num_bins, self_max_dend,
                                self_step_multiplier):

        """ Helper function for fill_voxels_dend, static method needed for NUMBA. """

        self_voxel_overflow_counter = 0

        # We use padding to include line segments where both points are outside, but their
        # line intersect the hypervoxel.
        padding = 1
        lower_padding_bound = 0 - padding
        upper_padding_bound = self_num_bins + padding

        section_id = section_data[point_idx, 0]
        section_x = section_data[point_idx, 1] * 1e-3  # Stored as section_x*1000 (since int)
        section_x[0] = 0

        coords = geometry[point_idx, :3]
        voxel_coords = (coords - self_hyper_voxel_origo) / self_voxel_size
        point_inside = np.sum(np.logical_and(lower_padding_bound <= voxel_coords,
                                             voxel_coords < upper_padding_bound),
                              axis=1) == 3
        scaled_soma_dist = geometry[point_idx, 4] * 1e6  # Dist to soma

        # Numba does not support third argument axis of np.diff, so transpose it instead
        # num_steps = np.ceil(np.amax(np.abs(np.diff(voxel_coords, axis=0)) * self_step_multiplier, axis=1)).astype(int)
        # dv_step = np.diff(voxel_coords, axis=0) / num_steps[:, None]

        step_diff = np.abs(np.diff(voxel_coords.T).T)
        sec_length = np.zeros((step_diff.shape[0], ))
        for i in range(0, step_diff.shape[0]):
            sec_length[i] = np.sqrt(step_diff[i, 0] ** 2 + step_diff[i, 1] ** 2 + step_diff[i, 2]**2)

        # TODO: Should we change the code to do the steps using coords rather than in voxel coords, and then
        #       convert to voxel coords after? (slower, but would be more accurate)

        num_steps = np.ceil(sec_length*self_step_multiplier).astype(np.int64)

        # OLD VERSION
        # This is just does same as numpy.amax for axis=1, but inlined for NUMBA
        # max_val = step_diff[:, 0].copy()
        # for j in range(1, step_diff.shape[1]):
        #     for i in range(0, step_diff.shape[0]):
        #         if max_val[i] < step_diff[i, j]:
        #             max_val[i] = step_diff[i, j]
        #
        # num_steps = np.ceil(max_val * self_step_multiplier).astype(np.int64)

        # TODO: num_steps should perhaps instead depend on the total length (with a small oversampling?)
        #       the reference should be that a line should occupy the same number of voxels regardless
        #       of what orientation it has.

        # [:, None] is here used to divide first row by first element in num_step
        # second row divide by second element in num_step etc. Pretty clever.
        # https://stackoverflow.com/questions/19602187/numpy-divide-each-row-by-a-vector-element
        # Looks like NUMBA cant handle this though,
        # dv_step = np.diff(voxel_coords.T).T / num_steps[:, None]

        dv_step = np.diff(voxel_coords.T).T
        for i in range(0, dv_step.shape[0]):
            dv_step[i, :] /= num_steps[i]

        ds_step = np.divide(np.diff(section_x), num_steps)
        dd_step = np.divide(np.diff(scaled_soma_dist), num_steps)

        # Remove this check later... should be done in morphology_data
        if (num_steps <= 0).any(): 
            print(f"Found zero length dendrite segment in neuron_id {neuron_id}")
            # import pdb
            # pdb.set_trace()

            raise ValueError(f"Found zero length dendrite segment (please check morphologies).")


        # Loop through all point-pairs of the section
        for idx in range(0, len(scaled_soma_dist)-1):

            if point_inside[idx] or point_inside[idx+1]:
                # Either of the points are within the cube + padding zone

                steps = np.arange(0, num_steps[idx] + 1)
                # vp = np.floor(voxel_coords[idx, :] + dv_step[idx, :] * steps[:, None]).astype(np.int64)
                # OBS! np.floor below is crucial -- np.floor(-0.2) = -1, vs int(-0.2) = 0
                vp_x = np.floor(voxel_coords[idx, 0] + dv_step[idx, 0] * steps).astype(np.int64)
                vp_y = np.floor(voxel_coords[idx, 1] + dv_step[idx, 1] * steps).astype(np.int64)
                vp_z = np.floor(voxel_coords[idx, 2] + dv_step[idx, 2] * steps).astype(np.int64)

                # Dirty fix no longer needed, section_x for first element set to 0 now
                # if idx == 0:
                #     # Dirty fix to handle that parent point has section_x = 1, but
                #     # child's first point should have section_x = 0
                #     s_x = (section_x[1] / num_steps[0]) * steps
                # else:
                s_x = section_x[idx] + ds_step[idx] * steps

                soma_dist = np.floor(scaled_soma_dist[idx] + dd_step[idx]*steps).astype(np.int64)

                # p_inside = np.sum(np.logical_and(0 <= vp, vp < self_num_bins), axis=1) == 3
                p_inside = np.logical_and(np.logical_and(0 <= vp_x, vp_x < self_num_bins[0]),
                                          np.logical_and(np.logical_and(0 <= vp_y, vp_y < self_num_bins[1]),
                                                         np.logical_and(0 <= vp_z, vp_z < self_num_bins[2])))

                # For each point pair in a section, find the intermediate points and mark voxels
                for i in steps:
                    if p_inside[i]:
                        # v_idx = tuple(vp[i, :])
                        v_idx = (vp_x[i], vp_y[i], vp_z[i])
                        v_ctr = voxel_space_ctr[v_idx]

                        if v_ctr > 0 and voxel_space[v_idx][v_ctr-1] == neuron_id:
                            # Voxel already contains neuron_id, skip
                            continue

                        if v_ctr < self_max_dend:
                            voxel_space[v_idx][v_ctr] = neuron_id
                            # Use section id from last point for section as first point has parent's section id
                            voxel_sec_id[v_idx][v_ctr] = section_id[-1]
                            voxel_sec_x[v_idx][v_ctr] = s_x[i]
                            voxel_soma_dist[v_idx][v_ctr] = soma_dist[i]
                            voxel_space_ctr[v_idx] += 1
                        else:
                            # Overflow, not enough space to store info
                            self_voxel_overflow_counter += 1

        return self_voxel_overflow_counter

    ############################################################################

    def fill_voxels_axon(self, voxel_space, voxel_space_ctr,
                         voxel_axon_dist,
                         neuron, neuron_id, section_id=None, subtree="neuron"):

        """
        Mark all voxels containing axons.

        voxel_space : n x n x n x k matrix holding the voxel content, normally self.axon_voxels (neuron IDs)
        voxel_space_ctr : n x n x n matrix holding count of how many items each voxel holds
        neuron : NeuronMorphologyExtended object
        neuron_id : ID of the neurons
        section_id : section id to add
        subtree : which subtree

        """

        for section in neuron.section_iterator_selective(section_type=2, section_id=section_id, subtree=subtree):

            voxel_overflow_ctr = SnuddaDetect.fill_voxels_axon_helper(voxel_space=voxel_space,
                                                                      voxel_space_ctr=voxel_space_ctr,
                                                                      voxel_axon_dist=voxel_axon_dist,
                                                                      point_idx=section.point_idx,
                                                                      geometry=section.morphology_data.geometry,
                                                                      neuron_id=neuron_id,
                                                                      self_hyper_voxel_origo=self.hyper_voxel_origo,
                                                                      self_voxel_size=self.voxel_size,
                                                                      self_num_bins=self.num_bins,
                                                                      self_max_axon=self.max_axon,
                                                                      self_step_multiplier=self.step_multiplier)

            self.voxel_overflow_counter += voxel_overflow_ctr

    # Temporarily disabling NUMBA, since amax does not support axis in NUMBA
    @staticmethod
    @jit(nopython=True, fastmath=True, cache=True)
    def fill_voxels_axon_helper(voxel_space,
                                voxel_space_ctr,
                                voxel_axon_dist,
                                point_idx,
                                geometry,
                                neuron_id: int,
                                self_hyper_voxel_origo,
                                self_voxel_size,
                                self_num_bins,
                                self_max_axon,
                                self_step_multiplier):

        self_voxel_overflow_counter = 0

        # We use padding to include line segments where both points are outside, but their
        # line intersect the hypervoxel.
        padding = 1
        lower_padding_bound = 0 - padding
        upper_padding_bound = self_num_bins + padding

        coords = geometry[point_idx, :3]
        voxel_coords = (coords - self_hyper_voxel_origo) / self_voxel_size
        point_inside = np.sum(np.logical_and(lower_padding_bound <= voxel_coords,
                                             voxel_coords < upper_padding_bound),
                              axis=1) == 3
        scaled_soma_dist = geometry[point_idx, 4] * 1e6  # Dist to soma

        # Numba does not support third argument axis of np.diff, so transpose it instead
        # num_steps = np.ceil(np.amax(np.abs(np.diff(voxel_coords, axis=0)) * self_step_multiplier, axis=1)).astype(int)
        # dv_step = np.diff(voxel_coords, axis=0) / num_steps[:, None]

        step_diff = np.abs(np.diff(voxel_coords.T).T)
        sec_length = np.zeros((step_diff.shape[0], ))
        for i in range(0, step_diff.shape[0]):
            sec_length[i] = np.sqrt(step_diff[i, 0] ** 2 + step_diff[i, 1] ** 2 + step_diff[i, 2]**2)

        num_steps = np.ceil(sec_length*self_step_multiplier).astype(np.int64)

        # OLD VERSION
        # This is just does same as numpy.amax for axis=1, but inlined for NUMBA
        # max_val = step_diff[:, 0].copy()
        # for j in range(1, step_diff.shape[1]):
        #     for i in range(0, step_diff.shape[0]):
        #         if max_val[i] < step_diff[i, j]:
        #             max_val[i] = step_diff[i, j]
        #
        # num_steps = np.ceil(max_val * self_step_multiplier).astype(np.int64)

        # [:, None] is here used to divide first row by first element in num_step
        # second row divide by second element in num_step etc. Pretty clever.
        # https://stackoverflow.com/questions/19602187/numpy-divide-each-row-by-a-vector-element
        # Looks like NUMBA cant handle this though,
        # dv_step = np.diff(voxel_coords.T).T / num_steps[:, None]
        dv_step = np.diff(voxel_coords.T).T
        for i in range(0, dv_step.shape[0]):
            dv_step[i, :] /= num_steps[i]

        dd_step = np.divide(np.diff(scaled_soma_dist), num_steps)

        # Remove this check later... should be done in morphology_data
        if (num_steps <= 0).any():
            print(f"Found zero length axon segment in neuron_id {neuron_id}")
            # Numba does not allow variables in exceptions...
            # raise ValueError(f"Found zero length axon segment in neuron_id {neuron_id}")
            raise ValueError(f"Found zero length axon segment (please check morphologies")

        # Loop through all point-pairs of the section
        for idx in range(0, len(scaled_soma_dist)-1):

            if point_inside[idx] or point_inside[idx+1]:
                # Either of the points are within the cube or padding zone

                steps = np.arange(0, num_steps[idx] + 1)
                # vp = np.floor(voxel_coords[idx, :] + dv_step[idx, :] * steps[:, None]).astype(np.int64)
                # OBS, we must have np.floor here, since np.floor(-0.2) = -1, and int(-0.2) = 0... we want the former
                vp_x = np.floor(voxel_coords[idx, 0] + dv_step[idx, 0] * steps).astype(np.int64)
                vp_y = np.floor(voxel_coords[idx, 1] + dv_step[idx, 1] * steps).astype(np.int64)
                vp_z = np.floor(voxel_coords[idx, 2] + dv_step[idx, 2] * steps).astype(np.int64)

                soma_dist = np.floor(scaled_soma_dist[idx] + dd_step[idx]*steps).astype(np.int64)

                # p_inside = np.sum(np.logical_and(0 <= vp, vp < self_num_bins), axis=1) == 3
                p_inside = np.logical_and(np.logical_and(0 <= vp_x, vp_x < self_num_bins[0]),
                                          np.logical_and(np.logical_and(0 <= vp_y, vp_y < self_num_bins[1]),
                                                         np.logical_and(0 <= vp_z, vp_z < self_num_bins[2])))

                # For each point pair in a section, find the intermediate points and mark voxels
                for i in steps:
                    if p_inside[i]:
                        # v_idx = tuple(vp[i, :])
                        v_idx = (vp_x[i], vp_y[i], vp_z[i])
                        v_ctr = voxel_space_ctr[v_idx]

                        if v_ctr > 0 and voxel_space[v_idx][v_ctr-1] == neuron_id:
                            # Voxel already contains neuron_id, skip
                            continue

                        if v_ctr < self_max_axon:
                            voxel_space[v_idx][v_ctr] = neuron_id
                            voxel_axon_dist[v_idx][v_ctr] = soma_dist[i]
                            voxel_space_ctr[v_idx] += 1
                        else:
                            # Overflow, not enough space to store info
                            self_voxel_overflow_counter += 1

        return self_voxel_overflow_counter

    ############################################################################

    # TODO: Add a filter, neurons that are not included in connectivity definition as either
    #       source or target are excluded. Also, if none of their sources/targets are in hypervoxel
    #       then the neurons are also excluded.

    def process_hyper_voxel(self, hyper_id):

        """
        Process hyper voxel, ie do touch detection, and save results.

        Args:
            hyper_id : ID of hyper voxel to process
        """

        start_time = timeit.default_timer()
        end_time = None

        try:
            if self.hyper_voxels[hyper_id]["neuron_ctr"] == 0:
                # No neurons, return quickly - do not write hdf5 file
                end_time = timeit.default_timer()
                return hyper_id, 0, 0, end_time - start_time, 0

            hyp_origo = self.hyper_voxels[hyper_id]["origo"]
            self.setup_hyper_voxel(hyp_origo, hyper_id)

            num_neurons = self.hyper_voxels[hyper_id]["neuron_ctr"]

            # TODO: Should we not force print this?
            self.write_log(f"Processing hyper voxel : {hyper_id}/{self.hyper_voxel_id_lookup.size}"
                           f" ({num_neurons} neurons)", force_print=True)

            # !!! Suggestion for optimisation. Place neurons with GJ first, then do
            # GJ touch detection, after that add rest of neurons (to get complete set)
            # and then do axon-dend synapse touch detection

            for neuron_id in sorted(self.hyper_voxels[hyper_id]["neurons"].keys()):

                neuron_info = self.hyper_voxels[hyper_id]["neurons"][neuron_id]
                neuron = self.load_neuron(self.neurons[neuron_id], use_cache=False)  # !!! Cached objects get huge

                if "soma" in neuron_info:
                    self.fill_voxels_soma(self.dend_voxels,
                                          self.dend_voxel_ctr,
                                          self.dend_sec_id,
                                          self.dend_sec_x,
                                          neuron,
                                          neuron_id)

                if "dend" in neuron_info:
                    section_id = neuron_info["dend"]

                    self.fill_voxels_dend(voxel_space=self.dend_voxels,
                                          voxel_space_ctr=self.dend_voxel_ctr,
                                          voxel_sec_id=self.dend_sec_id,
                                          voxel_sec_x=self.dend_sec_x,
                                          voxel_soma_dist=self.dend_soma_dist,
                                          neuron=neuron,
                                          neuron_id=neuron_id,
                                          section_id=section_id)

                if "axon" in neuron_info:
                    for section_id, subtree in neuron_info["axon"]:

                        self.fill_voxels_axon(voxel_space=self.axon_voxels,
                                              voxel_space_ctr=self.axon_voxel_ctr,
                                              voxel_axon_dist=self.axon_soma_dist,
                                              neuron=neuron,
                                              neuron_id=neuron_id,
                                              section_id=section_id,
                                              subtree=subtree)

            # This should be outside the neuron loop
            # This places axon voxels for neurons without axon morphologies
            self.place_synapses_no_axon(hyper_id,
                                        self.axon_voxels,
                                        self.axon_voxel_ctr,
                                        self.axon_soma_dist)

            # Finally this adds axon voxels for projections coming from other structures using projection maps
            try:
                self.projection_detection.voxelise_projections()
            except:
                # !!! TODO remove this bit of logging code
                import traceback
                t_str = traceback.format_exc()
                self.write_log(t_str, is_error=True)
                print(t_str)
                import pdb
                pdb.set_trace()

            # The normal voxel synapse detection
            self.detect_synapses()

            self.detect_gap_junctions()

            self.write_hyper_voxel_to_hdf5()

            end_time = timeit.default_timer()

            self.write_log(f"process_hyper_voxel: {hyper_id} took {end_time - start_time:.1f} s")

        except Exception as e:
            # Write error to log file to help trace it.
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str, is_error=True)

            sys.exit(-1)

        return (hyper_id, self.hyper_voxel_synapse_ctr,
                self.hyper_voxel_gap_junction_ctr, end_time - start_time,
                self.voxel_overflow_counter)

    ############################################################################

    # hyperID is just needed if we want to plotNeurons also

    def plot_hyper_voxel(self, plot_neurons=False, draw_axons=True, draw_dendrites=True,
                         draw_axon_voxels=True, draw_dendrite_voxels=True,
                         detect_done=True, elev_azim=None, show_axis=True, title=None,
                         fig_file_name=None, dpi=300, plot_neuron_id=None):

        """
        Plot hyper voxel.

        Args:
            plot_neurons : Should neurons be plotted
            draw_axons : Draw axons
            draw_dendrites : Draw dendrites
            draw_axon_voxels : Draw axon voxels marked
            draw_dendrite_voxels : Draw dendrite voxels marked
            detect_done :
            elev_azim : View angle
            show_axis : Show x,y,z axis?
            title : Title of plot
            fig_file_name : Fig name to save figure to
            dpi : Resolution

        """

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        colors = np.zeros((self.dend_voxel_ctr.shape[0],
                           self.dend_voxel_ctr.shape[1],
                           self.dend_voxel_ctr.shape[2], 4))
        colors[:, :, :, 3] = 0.3

        voxel_data = np.zeros((self.dend_voxel_ctr.shape[0],
                               self.dend_voxel_ctr.shape[1],
                               self.dend_voxel_ctr.shape[2]))

        if draw_axon_voxels:
            colors[:, :, :, 0] = self.axon_voxel_ctr / max(np.max(self.axon_voxel_ctr), 1)
            voxel_data += self.axon_voxel_ctr

        if draw_dendrite_voxels:
            colors[:, :, :, 2] = self.dend_voxel_ctr / max(np.max(self.dend_voxel_ctr), 1)
            voxel_data += self.dend_voxel_ctr

        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.add_subplot(projection='3d')
        ax.voxels(voxel_data > 0, facecolors=colors, edgecolor=None)

        if self.hyper_voxel_synapse_ctr > 0:
            syn_coord = self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr, 2:5]

            # In case hyperVoxelOffset has been applied, we need to subtract it
            # to draw within the hyper voxel
            if self.hyper_voxel_offset is not None:
                syn_coord -= self.hyper_voxel_offset
            ax.scatter(syn_coord[:, 0] + 0.5, syn_coord[:, 1] + 0.5, syn_coord[:, 2] + 0.5, c="green", s=64)

        if self.hyper_voxel_gap_junction_ctr > 0:
            gj_coord = self.hyper_voxel_gap_junctions[:self.hyper_voxel_gap_junction_ctr, 6:9]
            if self.hyper_voxel_offset is not None:
                gj_coord -= self.hyper_voxel_offset
            ax.scatter(gj_coord[:, 0] + 0.5, gj_coord[:, 1] + 0.5, gj_coord[:, 2] + 0.5, c="yellow", s=64)

        if elev_azim:
            ax.view_init(elev_azim[0], elev_azim[1])

        if not show_axis:
            plt.axis("off")

        # If plot is empty then this might be problem (ie axis scaled wrong)
        ax.set_xlim3d(0, self.hyper_voxel_size)
        ax.set_ylim3d(0, self.hyper_voxel_size)
        ax.set_zlim3d(0, self.hyper_voxel_size)

        ax.xaxis.set_tick_params(labelsize=18)
        ax.yaxis.set_tick_params(labelsize=18)
        ax.zaxis.set_tick_params(labelsize=18)

        plt.tight_layout()
        plt.ion()

        if plot_neurons:
            # Also plot the neurons overlayed, to verify

            num_neurons = self.hyper_voxels[self.hyper_voxel_id]["neuron_ctr"]

            if plot_neuron_id is None:
                try:
                    plot_neuron_id = np.array(list(self.hyper_voxels[self.hyper_voxel_id]["neurons"].keys()))

                    # plot_neuron_id = self.hyper_voxels[self.hyper_voxel_id]["neurons"][:num_neurons]
                except:
                    import traceback
                    print(traceback.format_exc())
                    import pdb
                    pdb.set_trace()

            for neuron_id in plot_neuron_id:
                neuron = self.load_neuron(self.neurons[neuron_id], use_cache=False)

                neuron.plot_neuron(axis=ax,
                                   plot_axon=draw_axons,
                                   plot_dendrite=draw_dendrites,
                                   plot_origo=self.hyper_voxel_origo, plot_scale=1 / self.voxel_size,
                                   soma_colour=(0, 0, 0),
                                   axon_colour=(1, 0, 0),
                                   dend_colour=(0, 0, 0))

        if title is None:
            plt.title("Scale is in voxels, not micrometers")
        else:
            plt.title(title)

        if fig_file_name is None:
            fig_file_name = f"Hypervoxel-{self.slurm_id}-{self.hyper_voxel_id}.png"

        if self.network_path is not None:
            fig_name = os.path.join(self.network_path, "figures", fig_file_name)
        else:
            fig_name = fig_file_name

        if not os.path.exists(os.path.dirname(fig_name)):
            if len(os.path.dirname(fig_name)) > 0:
                print(f"plot_hyper_voxel: Creating directory : {os.path.dirname(fig_name)}")
                os.mkdir(os.path.dirname(fig_name))

        plt.savefig(fig_name, dpi=dpi)

        plt.show()
        plt.ion()
        plt.pause(0.001)

        return plt, ax

    ############################################################################

    def export_voxel_visualisation_csv(self, neuron_id):

        """ Export CSV file with voxel data, used for visualisation."""

        # x,y,z = coords
        # shape = "cube" or "sphere"
        # type = "axon", "dendrite", "synapse"
        # id = neuron_id
        # x,y,z,shape,type,id

        header_str = "# x,y,z,shape,type,id\n"
        axon_str = ""
        dend_str = ""
        synapse_str = ""

        for x in range(0, self.axon_voxel_ctr.shape[0]):
            for y in range(0, self.axon_voxel_ctr.shape[1]):
                for z in range(0, self.axon_voxel_ctr.shape[2]):
                    for c in range(0, self.axon_voxel_ctr[x, y, z]):
                        n_id = self.axon_voxels[x, y, z, c]
                        if n_id in neuron_id:
                            axon_str += str(x) + "," + str(y) + "," + str(z) \
                                        + ",cube,axon," + str(n_id) + "\n"

        for x in range(0, self.dend_voxel_ctr.shape[0]):
            for y in range(0, self.dend_voxel_ctr.shape[1]):
                for z in range(0, self.dend_voxel_ctr.shape[2]):
                    for c in range(0, self.dend_voxel_ctr[x, y, z]):
                        n_id = self.dend_voxels[x, y, z, c]
                        if n_id in neuron_id:
                            dend_str += str(x) + "," + str(y) + "," + str(z) \
                                        + ",cube,dend," + str(n_id) + "\n"

        syn_list = []
        for ir, row in enumerate(self.hyper_voxel_synapses):
            if row[0] in neuron_id and row[1] in neuron_id:
                syn_list.append(ir)

        for i in syn_list:
            xyz = self.hyper_voxel_synapses[i, 2:5]
            synapse_str += str(xyz[0]) + "," + str(xyz[1]) + "," + str(xyz[2]) \
                           + ",sphere,synapse," + str(self.hyper_voxel_synapses[i, 1]) + "\n"

        f_name = "hypervoxel-" + str(self.slurm_id) + ".csv"
        with open(f_name, 'w') as f:
            f.write(header_str)
            f.write(axon_str)
            f.write(dend_str)
            f.write(synapse_str)

    ############################################################################

    # Example usage:
    # sd.plot_neurons_in_hyperVoxel(neuron_id=[1,20],
    #                               neuron_colour=np.array([[0,0,1],[0,1,0]]),
    #                               axon_alpha=[1,0.3], dend_alpha=[0.3,1])

    # each row in neuronColour is a colour for a neuron

    def plot_neurons_in_hyper_voxel(self, neuron_id, neuron_colour,
                                    axon_alpha=None, dend_alpha=None,
                                    show_plot=True, dpi=300):

        """
        Plot neurons in hyper voltage

        Args:
            neuron_id : ID of neurons to plot
            neuron_colour : Colur of neurons to plot
            axon_alpha : Alpha value of neuron axons
            dend_alpha : Alpha value of neuron dendrites
            show_plot : Should we dispaly the plot or keep it hidden
            dpi : Resolution of output file
        """

        if axon_alpha is None:
            axon_alpha = np.ones((len(neuron_id),))

        if dend_alpha is None:
            dend_alpha = np.ones((len(neuron_id),))

        alpha_axon_lookup = dict([])
        alpha_dend_lookup = dict([])
        neuron_colour_lookup = dict([])

        for ni, aa, da, nc in zip(neuron_id, axon_alpha, dend_alpha, neuron_colour):
            alpha_axon_lookup[ni] = aa
            alpha_dend_lookup[ni] = da
            neuron_colour_lookup[ni] = nc

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        colours = np.zeros((self.dend_voxel_ctr.shape[0],
                            self.dend_voxel_ctr.shape[1],
                            self.dend_voxel_ctr.shape[2], 4))

        voxel_data = np.zeros((self.dend_voxel_ctr.shape[0],
                               self.dend_voxel_ctr.shape[1],
                               self.dend_voxel_ctr.shape[2]))

        for ix in range(0, self.axon_voxel_ctr.shape[0]):
            for iy in range(0, self.axon_voxel_ctr.shape[1]):
                for iz in range(0, self.axon_voxel_ctr.shape[2]):
                    for ic in range(0, self.axon_voxel_ctr[ix, iy, iz]):
                        n_id = self.axon_voxels[ix, iy, iz, ic]
                        if n_id in neuron_id:
                            colours[ix, iy, iz, 0:3] = neuron_colour_lookup[n_id]
                            colours[ix, iy, iz, 3] = alpha_axon_lookup[n_id]
                            voxel_data[ix, iy, iz] = 1

        for ix in range(0, self.dend_voxel_ctr.shape[0]):
            for iy in range(0, self.dend_voxel_ctr.shape[1]):
                for iz in range(0, self.dend_voxel_ctr.shape[2]):
                    for ic in range(0, self.dend_voxel_ctr[ix, iy, iz]):
                        n_id = self.dend_voxels[ix, iy, iz, ic]
                        if n_id in neuron_id:
                            colours[ix, iy, iz, 0:3] = neuron_colour_lookup[n_id]
                            colours[ix, iy, iz, 3] = alpha_dend_lookup[n_id]
                            voxel_data[ix, iy, iz] = 1

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.voxels(voxel_data > 0,
                  facecolors=colours, edgecolor=None)

        print_list = []
        for ir, row in enumerate(self.hyper_voxel_synapses):
            if row[0] in neuron_id and row[1] in neuron_id:
                print_list.append(ir)

        # This should really only plot those between the neurons indicated
        if len(print_list) > 0:
            s_coord = self.hyper_voxel_synapses[print_list, 2:5]
            ax.scatter(s_coord[:, 0] + 0.5, s_coord[:, 1] + 0.5, s_coord[:, 2] + 0.5, c="red", s=100)

        plt.axis("off")

        fig_name = os.path.join(self.network_path, "figures",
                                f"Hypervoxel-{self.slurm_id}-{self.hyper_voxel_id}-someNeurons.png")

        if not os.path.exists(os.path.dirname(fig_name)):
            os.mkdir(os.path.dirname(fig_name))

        plt.savefig(fig_name, dpi=dpi)  # dpi = 900

        if show_plot:
            plt.ion()
            plt.show()
            plt.pause(0.001)

    ############################################################################

    # Memory check code taken from
    # https://stackoverflow.com/questions/17718449/determine-free-ram-in-python/17718729#17718729
    #

    @staticmethod
    def memory():

        memory_available, memory_total = snudda.utils.memory.memory_status()
        res = f"Memory: {memory_available} free, {memory_total} total"

        return res

    @staticmethod
    def memory_fraction_free():
        memory_available, memory_total = snudda.utils.memory.memory_status()

        return memory_available / memory_total

@staticmethod
def amax_helper(matrix) -> np.array:
    # Same as numpy.amax with axis=1

    max_val = matrix[:, 0].copy()
    for j in range(1, matrix.shape[1]):
        for i in range(0, matrix.shape[0]):
            if max_val[i] < matrix[i, j]:
                max_val[i] = matrix[i, j]

    return max_val

############################################################################


if __name__ == "__main__":
    print("Please do not call this file directly, use snudda.py")
    sys.exit(-1)
