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


import numpy as np
import os
import itertools

import time
import timeit

import h5py
import json
import pickle

from .Neuron_morphology import NeuronMorphology
from .load import SnuddaLoad

status = None
hyperVoxelData = None


class SnuddaDetect(object):

    def __init__(self,
                 config_file=None,
                 position_file=None,
                 voxel_size=3e-6,  # 2e-6,
                 hyper_voxel_size=100,  # 250, #100,
                 verbose=True,
                 logfile_name=None,
                 logfile=None,
                 save_file=None,
                 work_history_file=None,
                 restart_detection_flag=True,  # False = continue old detection
                 slurm_id=0,
                 volume_id=None,
                 role="master",
                 d_view=None,
                 lb_view=None,
                 rc=None,
                 axon_stump_id_flag=False,
                 h5libver="latest",
                 debug_flag=False):

        if rc is not None:
            d_view = rc.direct_view(targets='all')
            lb_view = rc.load_balanced_view(targets='all')

        assert role in ["master", "worker"], \
            "SnuddaDetect: Role must be master or worker"
        self.role = role

        self.verbose = verbose
        self.h5libver = h5libver
        self.debug_flag = debug_flag

        self.logfile = logfile
        self.config_file = config_file
        self.position_file = position_file
        self.save_file = save_file

        self.work_history_file = work_history_file  # Name of work history file
        self.work_history = None  # File pointer for actual file

        if logfile_name is None and logfile is not None:
            self.logfile_name = logfile.name
        else:
            self.logfile_name = logfile_name

        self.setup_log()

        self.write_log("Using hdf5 driver version: " + str(self.h5libver))

        mem = self.memory()
        self.write_log(str(mem))

        self.slurm_id = int(slurm_id)  # Make sure integer
        self.workers_initialised = False

        self.voxel_size = voxel_size
        self.hyper_voxel_size = hyper_voxel_size  # = N,  N x N x N voxels in a hyper voxel
        self.hyper_voxel_origo = np.zeros((3,))
        self.voxel_overflow_counter = 0

        self.num_bins = hyper_voxel_size * np.ones((3,), dtype=int)
        self.write_log("Each hyper voxel has %d x %d x %d voxels" \
                       % tuple(self.num_bins))

        # These are voxels that axons/dend occupy
        self.axon_voxels = None
        self.dend_voxels = None

        self.volume_id = volume_id
        if volume_id is not None:
            self.write_log("Touch detection only " + str(volume_id))
        else:
            self.write_log("Touch detecting all volumes")

        # These are counters, how many different axons/dend in the voxel
        self.axon_voxel_ctr = None
        self.dend_voxel_ctr = None

        self.max_axon_voxel_ctr = None
        self.max_dend_voxel_ctr = None

        # Columns in hyperVoxelSynapses:
        # 0: sourceCellID, 1: destCellID, 2: voxelX, 3: voxelY, 4: voxelZ,
        # 5: hyperVoxelID, 6: channelModelID,
        # 7: sourceAxonSomaDist (not SI scaled 1e6, micrometers),
        # 8: destDendSomaDist (not SI scalled 1e6, micrometers)
        # 9: destSegID, 10: destSegX (int 0 - 1000, SONATA wants float 0.0-1.0)
        # 11: conductance (int, not SI scaled 1e12, in pS)
        # 12: parameterID
        #
        # Note on parameterID:
        # If there are n parameter sets for the particular synapse type, then
        # the ID to use is parameterID % n, this way we can reuse connectivity
        # if we add more synapse parameter sets later.

        self.hyper_voxel_synapses = None

        # Columns in hyperVoxelGapJunctions
        # 0: sourceCellID, 1: destCellID, 2: sourceSegID, 3: destSegID,
        # 4: sourceSegX, 5: destSegX, 6: voxelX, 7: voxelY, 8: voxelZ,
        # 9: hyperVoxelID, 10: conductance (integer, in pS)
        self.hyper_voxel_gap_junctions = None

        self.hyper_voxel_synapse_ctr = 0
        self.hyper_voxel_gap_junction_ctr = 0

        self.hyper_voxel_coords = dict([])

        # This is used by the heap sort, when merging hdf5 files for the different
        # hyper voxels
        self.hyper_voxel_synapse_lookup = None
        self.hyper_voxel_gap_gunction_lookup = None

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

        # We have to dynamically create this lookup
        # self.synapseTypeLookup = { 1 : "GABA",
        #                           2 : "AMPA_NMDA",
        #                           3 : "GapJunction",
        #                           4 : "ACh",
        #                           5 : "NO"}
        #
        # self.synapseTypeReverseLookup = \
        #    {v: k for k, v in self.synapseTypeLookup.items()}

        self.connectivity_distributions = dict([])
        # self.connectivityDistributionsGJ = dict([])
        self.next_channel_model_id = 10

        self.prototype_neurons = dict([])

        self.axon_cum_density_cache = dict([])

        self.delete_old_merge()

        # Rather than load all neuron morphologies, we only load prototypes
        self.readPrototypes(configFile=config_file,
                            axonStumpIDFlag=axon_stump_id_flag)

        # Read positions
        self.read_neuron_positions(position_file)

        # Then we need to setup the workers

        if self.role == "master":

            self.setup_parallel(dView=d_view)

            if self.work_history_file is None:
                work_dir = os.path.dirname(self.save_file)
                work_dir = work_dir.replace("/voxels", "/")
                log_dir = work_dir + "/log/"
                # self.workHistoryFile = self.saveFile.replace(".hdf5","-worklog.hdf5")
                # self.workHistoryFile = self.workHistoryFile.replace("/voxels/","/")
                self.work_history_file = log_dir + "network-detect-worklog.hdf5"

                # workHistoryFile = re.sub("/voxels-\d+/","/",workHistoryFile)

                if self.work_history_file == self.save_file:
                    self.write_log("Unable to set a good worklog name")
                    self.work_history_file = "worklog.hdf5"

            if restart_detection_flag:
                if os.path.isfile(self.work_history_file):
                    self.write_log("Removing old work history file")
                    os.remove(self.work_history_file)

                # Setup new work history
                self.setup_work_history(self.work_history_file)
            else:
                # Open old file with work history
                print("Reusing old work history file " + str(self.work_history_file))
                self.work_history = h5py.File(self.work_history_file, "r+",
                                              libver=self.h5libver)

            # For each neuron we need to find which hyper voxel it belongs to
            # (can be more than one)
            self.distribute_neurons_parallel(dView=d_view)

            if d_view is not None:
                self.parallel_process_hyper_voxels(rc=rc, dView=d_view)

            else:
                # We are running it in serial

                (all_hyper_id, n_completed, remaining, self.voxel_overflow_counter) = \
                    self.setup_process_hyper_voxel_state_history()

                for hyper_id in remaining:  # self.hyperVoxels:
                    (hyper_id, n_syn, n_gj, exec_time, voxel_overflow_ctr) = \
                        self.processHyperVoxel(hyper_id)

                    if voxel_overflow_ctr > 0:
                        self.write_log("!!! HyperID " + str(hyper_id) + " OVERFLOWED " \
                                       + str(voxel_overflow_ctr) + " TIMES (" + \
                                       str(exec_time) + "s)")
                        self.voxel_overflow_counter += voxel_overflow_ctr
                    else:
                        self.write_log("HyperID " + str(hyper_id) + " completed - " \
                                       + str(n_syn) + " synapses and " \
                                       + str(n_gj) + " gap junctions found (" \
                                       + str(exec_time) + "s)")

                    self.update_process_hyper_voxel_state(hyperID=hyper_id, nSyn=n_syn, nGJ=n_gj,
                                                          execTime=exec_time,
                                                          voxelOverflowCounter=voxel_overflow_ctr)

        # We need to gather data from all the HDF5 files

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

    ############################################################################

    # The original code had problems with not being able to access the nc
    # object on the worker from the map call.

    def parallel_process_hyper_voxels(self, rc=None, dView=None):

        self.write_log("Starting parallelProcessHyperVoxels")

        start_time = timeit.default_timer()

        # Loads state if previously existed, otherwise creates new fresh history
        (allHyperIDs, nCompleted, remaining, self.voxel_overflow_counter) = \
            self.setup_process_hyper_voxel_state_history()

        n_workers = len(rc.ids)
        worker_status = [None for x in rc.ids]
        worker_idx = 0
        job_idx = 0
        busy_ctr = 0
        no_change_ctr = 0
        num_syn = 1  # If nSyn is zero delay in loop is shorter

        self.write_log("parallelProcessHyperVoxels: Using " \
                       + str(n_workers) + " worker")

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

                    self.update_process_hyper_voxel_state(hyperID=hyper_id,
                                                          nSyn=num_syn,
                                                          nGJ=n_gj,
                                                          execTime=exec_time,
                                                          voxelOverflowCounter=voxel_overflow_ctr)
                    worker_status[worker_idx] = None
                    rc[worker_idx]["result"] = None  # Clear to be safe
                    busy_ctr -= 1

                    if voxel_overflow_ctr > 0:
                        self.write_log("!!! HyperID " + str(hyper_id) + " OVERFLOWED " \
                                       + str(voxel_overflow_ctr) + " TIMES (" + \
                                       str(exec_time) + "s)")
                        self.voxel_overflow_counter += voxel_overflow_ctr
                    else:
                        self.write_log("HyperID " + str(hyper_id) + " completed - " \
                                       + str(num_syn) + " synapses found (" \
                                       + str(exec_time) + "s)")

            # Check that there are neurons in the hyper voxel, otherwise skip it.
            if worker_status[worker_idx] is None and job_idx < len(remaining):
                try:
                    self.write_log(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) \
                                   + " Starting hyper voxel " + str(remaining[job_idx]) \
                                   + " on worker " + str(worker_idx))
                    cmd_str = "result = nc.processHyperVoxel(" + str(remaining[job_idx]) + ")"

                    worker_status[worker_idx] = rc[worker_idx].execute(cmd_str, block=False)

                    job_idx += 1
                    busy_ctr += 1
                    no_change_ctr = 0

                except:
                    print("!!! Problem with worker : " + str(worker_idx) + " (" \
                          + str(n_workers) + " total workers)")
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)
                    import pdb
                    pdb.set_trace()

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

        self.write_log("Voxel overflows: " + str(self.voxel_overflow_counter))
        self.write_log("Total number of synapses: " \
                       + str(np.sum(self.work_history["nSynapses"][:])))
        self.write_log("parallelProcessHyperVoxels: " \
                       + str(end_time - start_time) + "s")

        self.work_history.close()

    ############################################################################

    def setup_work_history(self, work_history_file=None):

        if self.role != "master":
            return

        if work_history_file is None:
            work_history_file = self.work_history_file
        else:
            # Update internal state
            self.work_history_file = work_history_file

        self.write_log("Work history file: " + str(self.work_history_file))

        self.work_history = h5py.File(work_history_file, "w",
                                      libver=self.h5libver)

        # We need to encode the connectivityDistributions tuple as a string
        # before saving it in json
        # we also need to parse this string and recreate a tuple afterwards

        tmp_con_dist = dict([])
        tmp_con_dist_gj = dict([])

        for keys in self.connectivity_distributions:
            tmp_con_dist["$$".join(keys)] = self.connectivity_distributions[keys]

        # for keys in self.connectivityDistributionsGJ:
        #  tmpConDistGJ["$$".join(keys)] = self.connectivityDistributionsGJ[keys]

        save_meta_data = [(self.slurm_id, "SlurmID"),
                          (self.config_file, "configFile"),
                          (self.position_file, "positionFile"),
                          (self.voxel_size, "voxelSize"),
                          (self.hyper_voxel_size, "hyperVoxelSize"),
                          (self.axon_stump_id_flag, "axonStumpIDFlag"),
                          (json.dumps(self.config), "config"),
                          (json.dumps(tmp_con_dist),
                           "connectivityDistributions")]
        #    ,
        #                    (json.dumps(tmpConDistGJ),
        #                     "connectivityDistributionsGJ")]

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
                    data_name + " mismatch " + str(data) + " vs " \
                    + str(self.work_history["meta/" + data_name][()])

        print("Write neuron data to file")

        network_group = self.work_history.create_group("network")

        # Finally the neuron information
        neuron_group = network_group.create_group("neurons")

        # If the name list is longer than 20 chars, increase S20
        name_list = [n["name"].encode("ascii", "ignore") for n in self.neurons]
        str_type = 'S' + str(max(1, max([len(x) for x in name_list])))
        neuron_group.create_dataset("name", (len(name_list),), str_type, name_list,
                                    compression=self.h5compression)

        neuron_id_list = [n["neuronID"] for n in self.neurons]
        neuron_group.create_dataset("neuronID", (len(neuron_id_list),),
                                    'int', neuron_id_list)

        # Just make sure there is at least one neuron in volumeIDlist
        # that i inside volumeID

        volume_set = set([n["volumeID"] for n in self.neurons])
        assert self.volume_id is None or self.volume_id in volume_set, "VolumeID contains no neurons: " + str(
            self.volume_id)

        volume_id_list = [n["volumeID"].encode("ascii", "ignore") \
                          for n in self.neurons]
        str_type_vid = 'S' + str(max(1, max([len(x) for x in volume_id_list])))

        neuron_group.create_dataset("volumeID",
                                    (len(volume_id_list),), str_type_vid, volume_id_list,
                                    compression=self.h5compression)

        hoc_list = [n["hoc"].encode("ascii", "ignore") for n in self.neurons]
        neuron_group.create_dataset("hoc", (len(hoc_list),), 'S100', hoc_list,
                                    compression=self.h5compression)

        virtual_neuron_list = np.array([n["virtualNeuron"] for n in self.neurons],
                                       dtype=bool)
        virtual_neuron = neuron_group.create_dataset("virtualNeuron",
                                                     data=virtual_neuron_list,
                                                     compression=self.h5compression)

        swc_list = [n["morphology"].encode("ascii", "ignore") for n in self.neurons]
        max_swc_len = max([len(x) for x in swc_list])
        neuron_group.create_dataset("morphology", (len(swc_list),),
                                    'S' + str(max_swc_len), swc_list,
                                    compression=self.h5compression)

        neuron_position = neuron_group.create_dataset("position",
                                                      (len(self.neurons), 3),
                                                      "float",
                                                      compression=self.h5compression)

        neuron_rotation = neuron_group.create_dataset("rotation",
                                                      (len(self.neurons), 9),
                                                      "float",
                                                      compression=self.h5compression)

        neuron_dend_radius = neuron_group.create_dataset("maxDendRadius",
                                                         (len(self.neurons),),
                                                         "float",
                                                         compression=self.h5compression)

        neuron_axon_radius = neuron_group.create_dataset("maxAxonRadius",
                                                         (len(self.neurons),),
                                                         "float",
                                                         compression=self.h5compression)

        neuron_param_id = neuron_group.create_dataset("parameterID",
                                                      (len(self.neurons),),
                                                      "int",
                                                      compression=self.h5compression)

        neuron_modulation_id = neuron_group.create_dataset("modulationID",
                                                           (len(self.neurons),),
                                                           "int",
                                                           compression=self.h5compression)

        for (i, n) in enumerate(self.neurons):
            neuron_position[i] = n["position"]
            neuron_rotation[i] = n["rotation"].reshape(1, 9)
            neuron_dend_radius[i] = n["maxDendRadius"]
            neuron_axon_radius[i] = n["maxAxonRadius"]
            neuron_param_id[i] = n["parameterID"]
            neuron_modulation_id[i] = n["modulationID"]

        # Store input information
        neuron_group.create_dataset("populationUnitID", data=self.populationUnit,
                                    compression=self.h5compression, dtype=int)

        neuron_group.create_dataset("nPopulationUnits", data=self.nPopulationUnits)
        neuron_group.create_dataset("populationUnitPlacementMethod", data=self.populationUnitPlacementMethod)

        # Variable for axon density "r", "xyz" or "" (No axon density)
        axon_density_type = [n["axonDensityType"].encode("ascii", "ignore") \
                                 if n["axonDensityType"] is not None \
                                 else b"" \
                             for n in self.neurons]

        ad_str_type2 = "S" + str(max(1, max([len(x) if x is not None else 1 \
                                             for x in axon_density_type])))
        neuron_group.create_dataset("axonDensityType", (len(axon_density_type),),
                                    ad_str_type2, data=axon_density_type,
                                    compression=self.h5compression)

        axon_density = [n["axonDensity"].encode("ascii", "ignore") \
                            if n["axonDensity"] is not None \
                            else b"" \
                        for n in self.neurons]
        ad_str_type = "S" + str(max(1, max([len(x) if x is not None else 1 \
                                            for x in axon_density])))

        neuron_group.create_dataset("axonDensity", (len(axon_density),),
                                    ad_str_type, data=axon_density,
                                    compression=self.h5compression)

        axon_density_radius = [n["axonDensityRadius"] \
                                   if n["axonDensity"] is not None \
                                      and n["axonDensityType"] == "r" \
                                   else np.nan for n in self.neurons]

        neuron_group.create_dataset("axonDensityRadius", data=axon_density_radius)

        # This is for the density function where it uses x,y,z
        axon_density_bounds_xyz = np.nan * np.zeros((len(self.neurons), 6))

        for ni, n in enumerate(self.neurons):

            if n["axonDensity"] is None:
                # No axon density specified, skip
                continue

            if n["axonDensityType"] == "xyz":
                axon_density_bounds_xyz[ni, :] = n["axonDensityBoundsXYZ"]

        neuron_group.create_dataset("axonDensityBoundsXYZ", data=axon_density_bounds_xyz)

    ############################################################################

    # Reading work

    def setup_process_hyper_voxel_state_history(self):

        if "completed" in self.work_history:
            self.write_log("setupProcessHyperVoxelStateHistory: " \
                           + "Resuming from old state")
            # We already have a run in progress, load the state
            all_hyper_i_ds = set(self.work_history["allHyperIDs"])
            num_completed = int(self.work_history["nCompleted"][0])
            completed = set(self.work_history["completed"][:num_completed])
            remaining = self.sort_remaining_by_size(all_hyper_i_ds - completed)
            num_synapses = int(self.work_history["nSynapses"][0])
            num_gap_junctions = int(self.work_history["nGapJunctions"][0])
            voxel_overflow_counter = self.work_history["voxelOverflowCounter"][0]

        else:
            self.write_log("setupProcessHyperVoxelStateHistory: " \
                           + "Creating new work history.")
            # No history, add it to work history file
            num_hyper_voxels = len(self.hyper_voxels)
            minus_one = -1 * np.ones((num_hyper_voxels,), dtype=np.int32)
            self.work_history.create_dataset("completed", data=minus_one)

            num_completed = 0
            voxel_overflow_counter = 0

            # Could not rewrite scalars, so saving nCompleted as a vector of length 1
            self.work_history.create_dataset("nCompleted", data=np.zeros(1, ))
            all_hyper_i_ds = np.array([x for x in self.hyper_voxels.keys()], dtype=np.int)

            # Remove the empty hyper IDs
            (valid_hyper_id, empty_hyper_id) = self.remove_empty(all_hyper_i_ds)
            all_hyper_i_ds = valid_hyper_id
            remaining = self.sort_remaining_by_size(all_hyper_i_ds)

            # This should never happen,
            try:
                assert (np.array([self.hyper_voxels[x]["neuronCtr"] for x in
                                  empty_hyper_id]) == 0).all(), "All hyperIDs marked as empty are not empty!"
            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr)
                import pdb
                pdb.set_trace()

            self.write_log("Skipping " + str(len(empty_hyper_id)) + " empty hyper voxels")

            self.work_history.create_dataset("allHyperIDs", data=all_hyper_i_ds)
            self.work_history.create_dataset("nSynapses",
                                             data=np.zeros(num_hyper_voxels, ), dtype=np.int)
            self.work_history.create_dataset("nGapJunctions",
                                             data=np.zeros(num_hyper_voxels, ), dtype=np.int)
            self.work_history.create_dataset("voxelOverflowCounter",
                                             data=np.zeros(num_hyper_voxels, ), dtype=np.int)

        return all_hyper_i_ds, num_completed, remaining, voxel_overflow_counter

    ############################################################################

    # We want to do the hyper voxels with most neurons first, to minimize
    # the time waiting for lone cpu worker stragglers at the end.

    def sort_remaining_by_size(self, remaining):

        remaining = np.array(list(remaining), dtype=int)

        # Minus since we want them in descending order
        num_neurons = [-self.hyper_voxels[x]["neuronCtr"] for x in remaining]
        sort_idx = np.argsort(num_neurons)

        return remaining[sort_idx]

    ############################################################################

    def remove_empty(self, hyper_id):

        num_neurons = np.array([self.hyper_voxels[x]["neuronCtr"] for x in hyper_id])
        keep_idx = np.where(num_neurons > 0)[0]
        remove_idx = np.where(num_neurons == 0)[0]

        return hyper_id[keep_idx], hyper_id[remove_idx]

    ############################################################################

    def get_neuron_distribution_history(self):

        if "hyperVoxels" in self.work_history:
            self.write_log("Using neuron distribution from work history.")

            # We have hyper voxel information, load it
            hyper_voxels = dict([])

            for h_id_str in self.work_history["hyperVoxels"]:
                h_id = int(h_id_str)

                hyper_voxels[h_id] = dict([])

                hyper_voxels[h_id]["neurons"] = \
                    self.work_history["hyperVoxels"][h_id_str]["neurons"][()]

                hyper_voxels[h_id]["neuronCtr"] = \
                    self.work_history["hyperVoxels"][h_id_str]["neuronCtr"][()]

                hyper_voxels[h_id]["origo"] = \
                    self.work_history["hyperVoxels"][h_id_str]["origo"][()]

            hyper_voxel_id = self.work_history["meta/hyperVoxelIDs"][()]
            n_hyper_voxels = self.work_history["meta/nHyperVoxels"][()]
            simulation_origo = self.work_history["meta/simulationOrigo"][()]

            return hyper_voxels, hyper_voxel_id, n_hyper_voxels, simulation_origo
        else:
            # No information stored
            return None, None, None, None

    ############################################################################

    def save_neuron_distribution_history(self, hyperVoxels, minCoord, maxCoord):

        self.write_log("Writing neuron distribution history to file")

        assert "hyperVoxels" not in self.work_history, \
            "saveNeuronDistributionHistory should only be called once"

        self.work_history.create_dataset("meta/hyperVoxelIDs",
                                         data=self.hyperVoxelIDs)
        self.work_history.create_dataset("meta/nHyperVoxels",
                                         data=self.nHyperVoxels)
        self.work_history.create_dataset("meta/simulationOrigo",
                                         data=self.simulationOrigo)

        hv = self.work_history.create_group("hyperVoxels")

        for hID in hyperVoxels:
            hData = hv.create_group(str(hID))
            neurons = hyperVoxels[hID]["neurons"]
            neuronCtr = hyperVoxels[hID]["neuronCtr"]
            origo = hyperVoxels[hID]["origo"]
            hData.create_dataset("neurons", data=neurons[:neuronCtr])
            hData.create_dataset("neuronCtr", data=neuronCtr)
            hData.create_dataset("origo", data=origo)

    ############################################################################

    def update_process_hyper_voxel_state(self, hyperID, nSyn, nGJ, execTime,
                                         voxelOverflowCounter):

        nCompleted = self.work_history["nCompleted"][0]

        self.work_history["completed"][nCompleted] = hyperID
        self.work_history["nSynapses"][nCompleted] = nSyn
        self.work_history["nGapJunctions"][nCompleted] = nGJ
        self.work_history["voxelOverflowCounter"][nCompleted] = voxelOverflowCounter

        nCompleted += 1
        self.work_history["nCompleted"][0] = nCompleted

    ############################################################################

    def setupHyperVoxel(self, hyperVoxelOrigo, hyperVoxelID):

        # hypervoxel = a set of NxNxN voxels
        # hyperVoxelSynapses = list of all synapses detected in the hypervoxel

        self.hyper_voxel_coords[hyperVoxelID] = hyperVoxelOrigo  # Used???

        self.hyper_voxel_origo = hyperVoxelOrigo
        self.hyperVoxelID = hyperVoxelID

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

        # Clear lookup tables, just to be safe
        self.hyper_voxel_synapse_lookup = None
        self.hyper_voxel_gap_gunction_lookup = None

        # Used by plotHyperVoxel to make sure synapses are displayed correctly
        self.hyperVoxelOffset = None

        # Which axons populate the different voxels
        if self.axon_voxels is None:
            self.axon_voxels = np.zeros((self.num_bins[0],
                                         self.num_bins[1],
                                         self.num_bins[2],
                                         self.max_axon),
                                        dtype=np.int32)
            self.axon_voxel_ctr = np.zeros(self.num_bins, dtype=np.int32)

            # How far from the soma is this point
            self.axonSomaDist = np.zeros((self.num_bins[0],
                                          self.num_bins[1],
                                          self.num_bins[2],
                                          self.max_axon),
                                         dtype=np.int16)
        else:
            # Already allocated, just clear it
            self.axon_voxels[:] = 0
            self.axon_voxel_ctr[:] = 0
            self.axonSomaDist[:] = 0

        # Which dendrites populate the different voxels
        if self.dend_voxels is None:
            self.dend_voxels = np.zeros((self.num_bins[0],
                                         self.num_bins[1],
                                         self.num_bins[2],
                                         self.max_dend),
                                        dtype=np.int32)
            self.dend_voxel_ctr = np.zeros(self.num_bins, dtype=np.int32)

            # Which segment ID does the point belong to, and what segX
            self.dendSecID = np.zeros((self.num_bins[0],
                                       self.num_bins[1],
                                       self.num_bins[2],
                                       self.max_dend),
                                      dtype=np.int16)
            self.dendSecX = np.zeros((self.num_bins[0],
                                      self.num_bins[1],
                                      self.num_bins[2],
                                      self.max_dend),
                                     dtype=np.float16)  # 0 - 1.0, low pres

            # How far from the soma is this point
            self.dendSomaDist = np.zeros((self.num_bins[0],
                                          self.num_bins[1],
                                          self.num_bins[2],
                                          self.max_dend),
                                         dtype=np.int16)

        else:
            # Already allocated, just clear it
            self.dend_voxels[:] = 0
            self.dend_voxel_ctr[:] = 0
            self.dendSecID[:] = 0
            self.dendSecX[:] = 0
            self.dendSomaDist[:] = 0

        self.voxel_overflow_counter = 0

    ############################################################################

    # hyperID is only needed if we have neurons without axons, ie we use
    # axon density

    def detectSynapses(self):

        startTime = timeit.default_timer()

        # assert self.hyperVoxelSynapseCtr == 0 \
        #   and self.hyperVoxelSynapses is not None, \
        #   "setupHyperVoxel must be called before detecting synapses"

        # Find all voxels that contain axon and dendrites
        [xSyn, ySyn, zSyn] = np.where(np.bitwise_and(self.dend_voxel_ctr > 0,
                                                     self.axon_voxel_ctr > 0))

        if True:
            # This gives us some statistics, turn off later for speed
            self.max_axon_voxel_ctr = np.amax(self.axon_voxel_ctr)
            self.max_dend_voxel_ctr = np.amax(self.dend_voxel_ctr)

        for x, y, z in zip(xSyn, ySyn, zSyn):
            axonIDs = self.axon_voxels[x, y, z, :self.axon_voxel_ctr[x, y, z]]
            dendIDs = self.dend_voxels[x, y, z, :self.dend_voxel_ctr[x, y, z]]

            axonDist = self.axonSomaDist[x, y, z, :self.axon_voxel_ctr[x, y, z]]
            dendDist = self.dendSomaDist[x, y, z, :self.dend_voxel_ctr[x, y, z]]

            dendSecID = self.dendSecID[x, y, z, :self.dend_voxel_ctr[x, y, z]]
            dendSecX = self.dendSecX[x, y, z, :self.dend_voxel_ctr[x, y, z]]

            # Maybe make dendrite loop outer, since it has more variables?
            # speedup??
            for (axID, axDist) in zip(axonIDs, axonDist):
                for (dID, dSegID, dSegX, dDist) \
                        in zip(dendIDs, dendSecID, dendSecX, dendDist):

                    if axID == dID:
                        # Avoid self connections
                        continue

                    preType = self.neurons[axID]["type"]
                    postType = self.neurons[dID]["type"]

                    if (preType, postType) in self.connectivity_distributions:

                        conDict = self.connectivity_distributions[preType, postType]

                        # We need to loop over conDict in case there are multiple
                        # types of synapses from this neuron
                        for conType in conDict:
                            if conType == "GapJunction":
                                # This part detects only axon-dend synapses, skip gap junctions
                                continue

                            meanSynapseCond, stdSynapseCond = conDict[conType]["conductance"]
                            channelModelID = conDict[conType]["channelModelID"]

                            # We can not do pruning at this stage, since we only see
                            # synapses within hyper voxel, and pruning depends on
                            # all synapses between two connected cells.

                            # Do we have enough space allocated?
                            if self.hyper_voxel_synapse_ctr >= self.max_synapses:
                                self.resizeHyperVoxelSynapsesMatrix()

                            try:

                                # Synapse conductance varies between synapses
                                cond = np.random.normal(meanSynapseCond,
                                                        stdSynapseCond)

                                # Need to make sure the conductance is not negative,
                                # set lower cap at 10% of mean value
                                cond = np.maximum(cond, meanSynapseCond * 0.1)

                                paramID = np.random.randint(1000000)

                                # Add synapse
                                self.hyper_voxel_synapses[self.hyper_voxel_synapse_ctr, :] = \
                                    [axID, dID, x, y, z, self.hyperVoxelID, channelModelID,
                                     axDist, dDist, dSegID, dSegX * 1000, cond * 1e12, paramID]
                            except:
                                import traceback
                                tstr = traceback.format_exc()
                                self.write_log(tstr)
                                import pdb
                                pdb.set_trace()

                            # !!! OBS, dSegX is a value between 0 and 1, multiplied by 1000
                            # need to divide by 1000 later

                            self.hyper_voxel_synapse_ctr += 1

        # Sort the synapses (note sortIdx will not contain the empty rows
        # at the end.

        self.sortSynapses()

        # Convert from hyper voxel local coordinates to simulation coordinates
        # basically how many voxel steps do we need to take to go from
        # simulationOrigo to hyperVoxelOrigo (those were not included, so add them)
        try:
            hyperVoxelOffset = np.round((self.hyper_voxel_origo - self.simulationOrigo) \
                                        / self.hyperVoxelWidth).astype(int) \
                               * self.hyper_voxel_size

            # Just a double check...
            assert self.hyperVoxelIDs[int(np.round(hyperVoxelOffset[0] / self.hyper_voxel_size))][
                       int(np.round(hyperVoxelOffset[1] / self.hyper_voxel_size))][
                       int(np.round(hyperVoxelOffset[2] / self.hyper_voxel_size))] == self.hyperVoxelID, \
                "Internal inconsistency, have hyper voxel numbering or coordinates been changed?"

            self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr, :][:, range(2, 5)] \
                += hyperVoxelOffset

            # We need this in case plotHyperVoxel is called
            self.hyperVoxelOffset = hyperVoxelOffset

        except:
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)

            print("ARgh, what happened...")
            import pdb
            pdb.set_trace()

        # These are used when doing the heap sort of the hyper voxels
        self.hyper_voxel_synapse_lookup \
            = self.createLookupTable(data=self.hyper_voxel_synapses,
                                     nRows=self.hyper_voxel_synapse_ctr)

        # if(self.hyperVoxelSynapseCtr > 0 and self.hyperVoxelSynapseCtr < 10):
        #  self.plotHyperVoxel()
        #  import pdb
        #  pdb.set_trace()

        endTime = timeit.default_timer()

        self.write_log("detectSynapses: " + str(self.hyper_voxel_synapse_ctr) \
                       + " took " + str(endTime - startTime) + "s")

        if False and self.hyper_voxel_synapse_ctr > 0:
            print("First plot shows dendrites, and the voxels that were marked")
            print("Second plot same, but for axons")
            self.plotHyperVoxel(plotNeurons=True, drawAxons=False)
            self.plotHyperVoxel(plotNeurons=True, drawDendrites=False)

            import pdb
            pdb.set_trace()

        return self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr, :]

    ############################################################################

    def placeSynapsesNoAxon(self, hyperID, voxelSpace, voxelSpaceCtr,
                            voxelAxonDist):

        startTime = timeit.default_timer()

        # 1. Find neurons within hyper voxel that have no axon

        nNeurons = self.hyper_voxels[hyperID]["neuronCtr"]
        hNeurons = self.hyper_voxels[hyperID]["neurons"][:nNeurons]

        noAxonNeurons = [self.neurons[x] for x in hNeurons
                         if self.neurons[x]["axonDensity"] is not None]

        if len(noAxonNeurons) == 0:
            # No neurons without axons
            return

        for naNeuron in noAxonNeurons:

            # There are two types of axon density specified
            # - Spherically symmetric
            # - f(x,y,z) in SWC coordinates

            if naNeuron["axonDensityType"] == "r":

                # 2. Check that we have cumulative probability distribution for
                #    radial distance, if not compute and cache

                if naNeuron["type"] in self.axon_cum_density_cache:
                    (naCumDensity, naPoints) = self.axon_cum_density_cache[naNeuron["type"]]

                    self.write_log("Placing " + str(naPoints) + " random axon points for " \
                                   + str(naNeuron["name"]) + "(cached)")

                else:
                    radius = np.arange(0,
                                       naNeuron["axonDensityRadius"] + self.voxel_size,
                                       self.voxel_size)

                    density_as_func = eval('lambda r: ' + naNeuron["axonDensity"])
                    naPDensity = np.array([density_as_func(r) for r in radius])

                    # We need to scale by distance squared, since in the shell at distance
                    # d further from the soma has more voxels in it than a shell closer
                    # This cumulative distribution is only used to determine how far
                    # from the soma a synapse is located (not the direction)

                    # !!! Plot and verify this !!!
                    naCumDensity = np.cumsum(np.multiply(naPDensity, radius ** 2))
                    naCumDensity /= naCumDensity[-1]  # Normalise

                    # 3. Calculate how many points there should be within volume
                    #    based on (unscaled raw) probability density
                    # Volume at each distance is 4*pi*(r**2) * voxelSize
                    naPoints = int(np.round(np.sum(4 * np.pi * self.voxel_size \
                                                   * np.multiply(radius ** 2, naPDensity))))

                    self.write_log("Placing " + str(naPoints) + " random axon points for " \
                                   + str(naNeuron["name"]))

                    self.axon_cum_density_cache[naNeuron["type"]] = (naCumDensity, naPoints)

                # print("Check naPoints")
                # import pdb
                # pdb.set_trace()

                # 4. Randomize the points
                (naVoxelCoords, naAxonDist) = \
                    self.noAxonPointsSphere(naNeuron["position"],
                                            naCumDensity,
                                            naPoints)

            elif naNeuron["axonDensityType"] == "xyz":

                axonDensityFunc = eval("lambda x,y,z: " + naNeuron["axonDensity"])

                (naVoxelCoords, naAxonDist) = \
                    self.noAxonPointsXYZ(naNeuron["position"],
                                         naNeuron["rotation"],
                                         axonDensityFunc,
                                         naNeuron["axonDensityBoundsXYZ"])
            else:
                self.write_log("Unknown axonDensityType: " \
                               + str(naNeuron["axonDensityType"]) \
                               + "\n" + str(naNeuron))

            neuronID = naNeuron["neuronID"]

            for idx in range(0, naVoxelCoords.shape[0]):
                try:
                    xIdx = naVoxelCoords[idx, 0]
                    yIdx = naVoxelCoords[idx, 1]
                    zIdx = naVoxelCoords[idx, 2]
                    axonDist = naAxonDist[idx]
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)

                    import pdb
                    pdb.set_trace()

                vCtr = voxelSpaceCtr[xIdx, yIdx, zIdx]
                if vCtr > 0 and voxelSpace[xIdx, yIdx, zIdx, vCtr - 1] == neuronID:
                    # Voxel already has neuronID, skip
                    continue

                try:
                    voxelSpace[xIdx, yIdx, zIdx, vCtr] = neuronID
                    voxelAxonDist[xIdx, yIdx, zIdx, vCtr] = axonDist
                    voxelSpaceCtr[xIdx, yIdx, zIdx] += 1
                except:
                    self.voxel_overflow_counter += 1
                    self.write_log("!!! Axon voxel space overflow: " \
                                   + str(voxelSpaceCtr[xIdx, yIdx, zIdx]))

            # if(True):
            #  # Debug plot
            #  self.plotHyperVoxel()
            #  import pdb
            #  pdb.set_trace()

        endTime = timeit.default_timer()

        self.write_log("placeSynapsesNoAxonSphere: " + str(endTime - startTime) + "s" \
                       + ", hyperID: " + str(hyperID))

    ############################################################################

    # This picks points around soma centre. nPoints are randomized, points
    # outside the hyper sphere are rejected, so fewer than nPoints might be
    # returned.

    def noAxonPointsSphere(self, somaCentre, rCumDistribution, nPoints):

        uvr = np.random.rand(nPoints, 3)
        theta = 2 * np.pi * uvr[:, 0]
        phi = np.arccos(2 * uvr[:, 1] - 1)

        # Double check these are sorted
        # We want to sample from the supplied distance distribution
        rP = np.sort(uvr[:, 2] * rCumDistribution[-1], axis=0)
        nextIdx = 0

        print("nPoints = " + str(nPoints))

        r = np.zeros((nPoints,))

        for posIdx, rP1 in enumerate(rP):
            try:
                while rP1 > rCumDistribution[nextIdx + 1]:
                    nextIdx += 1

                r[posIdx] = nextIdx

            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr)
                import pdb
                pdb.set_trace()

        # Rescale index to radie
        r = r * self.voxel_size

        xyz = np.array([r * np.sin(theta) * np.cos(phi),
                        r * np.sin(theta) * np.sin(phi),
                        r * np.cos(theta)]).transpose() + somaCentre

        # Check which points are inside this hyper voxel
        voxIdx = np.floor((xyz - self.hyper_voxel_origo) / self.voxel_size).astype(int)

        # print("Verify that the points are on a sphere, and that inside check is ok")
        # import pdb
        # pdb.set_trace()

        insideIdx = np.where(np.sum(np.bitwise_and(0 <= voxIdx,
                                                   voxIdx < self.hyper_voxel_size),
                                    axis=1) == 3)[0]

        return voxIdx[insideIdx, :], r[insideIdx]

    ############################################################################

    # Helper function to give points inside axon bounding box, that are
    # inside hyper voxel

    def getHVaxonPoints(self,
                        neuronPosition,
                        rotation,
                        axonDensityBoundsXYZ,
                        nPoints=1000):

        # Randomly place nPoints inside bounding box (SWC coordintes, soma (0,0,0))
        xMin = axonDensityBoundsXYZ[0]
        xWidth = axonDensityBoundsXYZ[1] - axonDensityBoundsXYZ[0]
        yMin = axonDensityBoundsXYZ[2]
        yWidth = axonDensityBoundsXYZ[3] - axonDensityBoundsXYZ[2]
        zMin = axonDensityBoundsXYZ[4]
        zWidth = axonDensityBoundsXYZ[5] - axonDensityBoundsXYZ[4]

        xyz = np.random.rand(nPoints, 3)
        xyz[:, 0] = xMin + xWidth * xyz[:, 0]
        xyz[:, 1] = yMin + yWidth * xyz[:, 1]
        xyz[:, 2] = zMin + zWidth * xyz[:, 2]

        # Check which of the points are inside hyper voxel (rotate+translate)
        voxIdx = ((np.matmul(rotation, xyz.transpose()).transpose()
                   + neuronPosition - self.hyper_voxel_origo)
                  / self.voxel_size).astype(int)

        insideIdx = np.where(np.sum(np.bitwise_and(0 <= voxIdx,
                                                   voxIdx < self.hyper_voxel_size),
                                    axis=1) == 3)[0]

        return xyz[insideIdx, :], voxIdx[insideIdx, :]

    ############################################################################

    # somaCentre and rotation of neuron
    # axonDensityFunc should be written so that it can handle x,y,z (SWC coords)
    # as vectors
    # axonDensityBoundsXYZ = [xmin,xmax,ymin,ymax,zmin,zmax] in SWC coordinates

    # axonDensityFunc = eval("lambda x,y,z: " + axonPstr)

    def noAxonPointsXYZ(self, neuronPosition, rotation,
                        axonDensityFunc, axonDensityBoundsXYZ):

        # Points for initial sample
        nPoints = 5000

        (xyzInside, voxIdx) = self.getHVaxonPoints(neuronPosition,
                                                   rotation,
                                                   axonDensityBoundsXYZ,
                                                   nPoints)

        xWidth = axonDensityBoundsXYZ[1] - axonDensityBoundsXYZ[0]
        yWidth = axonDensityBoundsXYZ[3] - axonDensityBoundsXYZ[2]
        zWidth = axonDensityBoundsXYZ[5] - axonDensityBoundsXYZ[4]

        pointVolume = xWidth * yWidth * zWidth / nPoints
        voxelVolume = self.voxel_size ** 3

        # If no points were inside HV, then the intersection must be small
        # so we assume no voxels should be filled
        if xyzInside.shape[0] == 0:
            # Double check correct data-types
            self.write_log("Bounding box appears to be outside hyper voxel")
            return np.zeros((0, 3), dtype=int), np.zeros((0, 1))

            # Calculate density at each of the points inside HV
        densityInside = axonDensityFunc(xyzInside[:, 0],
                                        xyzInside[:, 1],
                                        xyzInside[:, 2])

        # Estimate number of synapses from density, in this case we use a volume
        # equal to bounding box volume / nPoints for each point.
        # OBS that we only want to know how many synapses to place inside HV
        nPointsToPlace = np.round(np.sum(densityInside * pointVolume))

        if nPointsToPlace <= 0:
            # To little probability mass inside
            self.write_log("Too little of axon appears to be inside hyper voxel")
            return np.zeros((0, 3), dtype=int), np.zeros((0, 1))

        # Calculate max density inside HV, divide all densities by that to
        # get Pkeep.
        maxDensity = np.max(densityInside)

        # We know that n out of N points placed were inside volume, so volume
        # acceptance rate is n/N.
        # In order to get about nPointsToPlace points placed we need to account
        # for how many outside volume, and also account for how many of those
        # inside volume gets accepted (Pkeep = density / maxDensity)
        nTries = np.round(nPointsToPlace * nPoints \
                          / np.sum(densityInside / maxDensity)).astype(int)

        if nTries > 1e6:
            self.write_log("!!! noAxonPointsXYZ: Warning trying to place " \
                           + str(nTries) + " points. Bounds: " \
                           + str(axonDensityBoundsXYZ))

        # Only print this in verbose mode
        if self.verbose:
            self.write_log("Trying to place " + str(nTries) + " points to get " \
                           + str(nPointsToPlace) + " axon voxels filled")

        if nPoints >= nTries:
            # We already have enough points, use a subset
            useNum = np.round(voxIdx.shape[0] * nTries / nPoints).astype(int)
            pickedIdx = np.where(np.random.rand(useNum) \
                                 < densityInside[:useNum] / maxDensity)[0]
            axonDist = np.sqrt(np.sum((xyzInside[pickedIdx, :]) ** 2, axis=1))

            return voxIdx[pickedIdx, :], axonDist
        else:
            # Not enough points, use the ones we have, then get more
            pickedIdx = np.where(np.random.rand(voxIdx.shape[0])
                                 < densityInside / maxDensity)[0]
            axonDist = np.sqrt(np.sum((xyzInside[pickedIdx, :]) ** 2, axis=1))

            # Get more points

            (xyzInsideB, voxIdxB) = \
                self.getHVaxonPoints(neuronPosition,
                                     rotation,
                                     axonDensityBoundsXYZ,
                                     nTries - nPoints)

            densityInsideB = axonDensityFunc(xyzInsideB[:, 0],
                                             xyzInsideB[:, 1],
                                             xyzInsideB[:, 2])

            pickedIdxB = np.where(np.random.rand(voxIdxB.shape[0]) \
                                  < densityInsideB / maxDensity)[0]

            axonDistB = np.sqrt(np.sum((xyzInsideB[pickedIdxB, :]) ** 2, axis=1))

            return (np.concatenate([voxIdx[pickedIdx, :],
                                    voxIdxB[pickedIdxB, :]]),
                    np.concatenate([axonDist, axonDistB]))

    ############################################################################

    def resizeHyperVoxelSynapsesMatrix(self, newSize=None):

        if newSize is None:
            newSize = int(np.ceil(1.5 * self.max_synapses))

        assert newSize >= self.hyper_voxel_synapse_ctr, \
            " Can not shrink below existing number of synapses"

        # We need to increase matrix size
        old = self.hyper_voxel_synapses
        self.max_synapses = newSize
        self.write_log("!!! Increasing max synapses to " + str(self.max_synapses))
        self.hyper_voxel_synapses = np.zeros((self.max_synapses, 13),
                                             dtype=np.int32)
        self.hyper_voxel_synapses[:old.shape[0], :] = old
        del old

    ############################################################################

    # This truncates and sorts the hyper voxel synapse matrix

    def sortSynapses(self):

        sortIdx = np.lexsort(self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr,
                             [0, 1]].transpose())

        self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr, :] = \
            self.hyper_voxel_synapses[sortIdx, :]

    ############################################################################

    def sortGapJunctions(self):

        sortIdx = \
            np.lexsort(self.hyper_voxel_gap_junctions[:self.hyper_voxel_gap_junction_ctr,
                       [0, 1]].transpose())

        self.hyper_voxel_gap_junctions[:self.hyper_voxel_gap_junction_ctr, :] = \
            self.hyper_voxel_gap_junctions[sortIdx, :]

    ############################################################################

    # First and second column specify the source and destination ID of a synapse
    # or a gap junction.
    #
    # This creates a lookup table where all synapses in the hyper voxel
    # between the same pair of neurons are grouped together
    # returns a matrix where first column is a UID = srcID*nNeurons + destID
    # and the following two columns are start row and end row (-1) in matrix

    def createLookupTable(self, data, nRows):

        self.write_log("Create lookup table")
        # nRows = data.shape[0] -- zero padded, cant use shape
        lookupTable = np.zeros((data.shape[0], 3), dtype=int)

        nextIdx = 0
        startIdx = 0

        lookupIdx = 0
        nNeurons = len(self.neurons)

        while nextIdx < nRows:
            srcID = data[nextIdx, 0]
            destID = data[nextIdx, 1]
            nextIdx += 1

            while (nextIdx < nRows
                   and data[nextIdx, 0] == srcID
                   and data[nextIdx, 1] == destID):
                nextIdx += 1

            lookupTable[lookupIdx, :] = [destID * nNeurons + srcID, startIdx, nextIdx]

            startIdx = nextIdx
            lookupIdx += 1

        return lookupTable[:lookupIdx, :]

    ############################################################################

    def includesGapJunctions(self):

        hasGapJunctions = False

        for key in self.connectivity_distributions:
            if "GapJunction" in self.connectivity_distributions[key]:
                hasGapJunctions = True

        return hasGapJunctions

    ############################################################################

    # Gap junctions are stored in self.hyperVoxelGapJunctions

    def detectGapJunctions(self):

        if not self.includesGapJunctions():
            self.write_log("detectGapJunctions: No gap junctions defined in connectivity rules")
            return

        startTime = timeit.default_timer()

        assert self.hyper_voxel_gap_junction_ctr == 0 \
               and self.hyper_voxel_gap_junctions is not None, \
            "setupHyperVoxel must be called before detecting gap junctions"

        [xDV, yDV, zDV] = np.where(self.dend_voxel_ctr > 0)

        for x, y, z in zip(xDV, yDV, zDV):

            neuronIDs = self.dend_voxels[x, y, z, :self.dend_voxel_ctr[x, y, z]]
            segID = self.dendSecID[x, y, z, :self.dend_voxel_ctr[x, y, z]]
            segX = self.dendSecX[x, y, z, :self.dend_voxel_ctr[x, y, z]]

            # All possible pairs
            for pairs in itertools.combinations(np.arange(0, self.dend_voxel_ctr[x, y, z]), 2):
                neuronID1 = self.dend_voxels[x, y, z, pairs[0]]
                neuronID2 = self.dend_voxels[x, y, z, pairs[1]]

                # !!! Check no self connections??

                # Add type field, derived from name field MSD1_45 --> MSD1
                preType = self.neurons[neuronID1]["type"]
                postType = self.neurons[neuronID2]["type"]

                if (preType, postType) in self.connectivity_distributions:

                    if ("GapJunction"
                            in self.connectivity_distributions[preType, postType]):
                        conInfo \
                            = self.connectivity_distributions[preType, postType]["GapJunction"]

                        segID1 = self.dendSecID[x, y, z, pairs[0]]
                        segID2 = self.dendSecID[x, y, z, pairs[1]]

                        segX1 = self.dendSecX[x, y, z, pairs[0]]
                        segX2 = self.dendSecX[x, y, z, pairs[1]]

                        meanGJcond, stdGJcond = conInfo["conductance"]

                        # !!! Currently not using channelParamDict for GJ

                        GJcond = np.random.normal(meanGJcond, stdGJcond)
                        GJcond = np.maximum(GJcond, meanGJcond * 0.1)  # Avoid negative cond

                        self.hyper_voxel_gap_junctions[self.hyper_voxel_gap_junction_ctr, :] = \
                            [neuronID1, neuronID2, segID1, segID2, segX1 * 1e3, segX2 * 1e3,
                             x, y, z, self.hyperVoxelID, GJcond * 1e12]
                        self.hyper_voxel_gap_junction_ctr += 1

        self.sortGapJunctions()

        # We also translate from local hyper voxel coordinates to simulation
        # voxel coordinates

        hyperVoxelOffset = np.round((self.hyper_voxel_origo
                                     - self.simulationOrigo) \
                                    / self.hyperVoxelWidth).astype(int) \
                           * self.hyper_voxel_size

        self.hyper_voxel_gap_junctions[:self.hyper_voxel_gap_junction_ctr, :][:, range(6, 9)] \
            += hyperVoxelOffset

        self.hyper_voxel_gap_gunction_lookup \
            = self.createLookupTable(data=self.hyper_voxel_gap_junctions,
                                     nRows=self.hyper_voxel_gap_junction_ctr)

        endTime = timeit.default_timer()

        self.write_log("detectGapJunctions: " + str(endTime - startTime) + "s")

        return self.hyper_voxel_gap_junctions[:self.hyper_voxel_gap_junction_ctr, :]

    ############################################################################

    def setup_log(self, logFileName=None):

        if logFileName is None:
            logFileName = self.logfile_name

        if logFileName is None or len(logFileName) == 0:
            # Not a valid log file name
            return

        if self.logfile is not None:
            self.write_log("Already have a log file setup, ignoring")
            return

        self.logfile = open(logFileName, 'wt')

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

    def readPrototypes(self, configFile=None, axonStumpIDFlag=False):

        if configFile is None:
            configFile = self.config_file

        configFile = self.getPath(configFile)

        self.axon_stump_id_flag = axonStumpIDFlag

        print("Loading from " + configFile)

        cfgFile = open(str(configFile), 'r')

        try:
            self.config = json.load(cfgFile)
        finally:
            cfgFile.close()

        self.prototype_neurons = dict()

        for name, definition in self.config["Neurons"].items():

            self.write_log("Reading prototype for: " + name)

            morph = self.getPath(definition["morphology"])
            param = self.getPath(definition["parameters"])
            mech = self.getPath(definition["mechanisms"])

            if "neuronType" in definition:
                neuronType = definition["neuronType"]
            else:
                neuronType = "neuron"

            if neuronType == "virtual":
                virtualNeuron = True
            else:
                virtualNeuron = False

            if 'hoc' in definition:
                hoc = definition["hoc"]
            else:
                hoc = None

            self.prototype_neurons[name] \
                = NeuronMorphology(name=name,
                                   swc_filename=morph,
                                   param_filename=param,
                                   mech_filename=mech,
                                   hoc=hoc,
                                   virtualNeuron=virtualNeuron,
                                   axonStumpIDFlag=axonStumpIDFlag)

            if "axonDensity" in definition:
                self.write_log("Setting axon density for neuron without axon")
                axonDensityType = definition["axonDensity"][0]

                if axonDensityType == "r":
                    density = definition["axonDensity"][1]
                    maxRadius = definition["axonDensity"][2]

                    self.prototype_neurons[name].setAxonVoxelRadialDensity(density,
                                                                           maxRadius)
                elif axonDensityType == "xyz":
                    density = definition["axonDensity"][1]
                    axonDensityBoundsXYZ = np.array(definition["axonDensity"][2])

                    self.prototype_neurons[name].setAxonVoxelXYZDensity(density,
                                                                        axonDensityBoundsXYZ)

                else:
                    self.write_log(str(name) + ": Unknown axon density type : " \
                                   + str(axonDensityType) \
                                   + "\n" + str(definition["axonDensity"]))

            else:
                # If no axon density specified, then axon must be present in morphology
                assert (len(self.prototype_neurons[name].axon) > 0), \
                    "File: " + morph + " does not have an axon"

            assert len(self.prototype_neurons[name].dend) > 0 \
                   or self.prototype_neurons[name].virtualNeuron, \
                "File: " + morph + " does not have a dendrite"

            # Since we already have the config file open, let's read connectivity
            # distributions also

        self.write_log("Loading connectivity information")

        for name, definition in self.config["Connectivity"].items():

            preType, postType = name.split(",")

            conDef = definition.copy()

            for key in conDef:
                if key == "GapJunction":
                    conDef[key]["channelModelID"] = 3
                else:
                    conDef[key]["channelModelID"] = self.next_channel_model_id
                    self.next_channel_model_id += 1

                # Also if conductance is just a number, add std 0
                if type(conDef[key]["conductance"]) not in [list, tuple]:
                    conDef[key]["conductance"] = [conDef[key]["conductance"], 0]

            self.connectivity_distributions[preType, postType] = conDef

    ############################################################################

    def read_neuron_positions(self, positionFile):

        if positionFile is None:
            positionFile = self.position_file

        mem = self.memory()
        self.write_log(str(mem))

        self.write_log("Reading positions from file: " + positionFile)

        posInfo = SnuddaLoad(positionFile).data

        mem = self.memory()
        self.write_log(str(mem))

        # Make sure we do not change config file unintentionally
        assert posInfo["configFile"] == self.config_file, \
            "Not using original config file: " \
            + str(posInfo["configFile"]) + " vs " + self.config_file

        self.neurons = posInfo["neurons"]
        nNeurons = len(self.neurons)

        self.neuronPositions = np.zeros((nNeurons, 3))

        for ni, neuron in enumerate(posInfo["neurons"]):
            self.neuronPositions[ni, :] = neuron["position"]

            # Add a few sanity checks
            assert ni == neuron["neuronID"], \
                "NeuronID=" + str(neuron["neuronID"]) + "and ni=" + str(ni) \
                + " not equal, corruption?"
            assert neuron["name"] in self.prototype_neurons, \
                "Neuron type " + neuron["name"] + " not in prototypeNeurons: " \
                + str(self.prototype_neurons)

        # Also load the channel data
        self.nPopulationUnits = posInfo["nPopulationUnits"]
        self.populationUnit = posInfo["populationUnit"]
        self.populationUnitPlacementMethod = posInfo["populationUnitPlacementMethod"]

        self.populationUnits = dict([])
        for i in range(0, self.nPopulationUnits):
            self.populationUnits[i] = np.where(self.populationUnit == i)[0]

        self.write_log("Position file read.")
        del posInfo

    ############################################################################

    # If the detect is rerun we need to make sure there are not old MERGE
    # files left that might remember old run accidentally

    def delete_old_merge(self):

        if self.role == "master":

            workDir = os.path.dirname(self.save_file)
            workDir = (workDir + "/").replace("/voxels/", "/")

            delFiles = [workDir + "network-putative-synapses-MERGED.hdf5",
                        workDir + "network-putative-synapses-MERGED.hdf5-cache",
                        workDir + "network-pruned-synapses.hdf5",
                        workDir + "network-pruned-synapses.hdf5-cache"]

            for f in delFiles:
                if os.path.exists(f):
                    self.write_log("Removing old files " + str(f))
                    os.remove(f)

    ############################################################################

    def writeHyperVoxelToHDF5(self):

        startTime = timeit.default_timer()

        outputName = self.save_file.replace(".hdf5", "-" + str(self.hyperVoxelID) + ".hdf5")

        with h5py.File(outputName, "w", libver=self.h5libver) as outFile:

            configData = outFile.create_dataset("config",
                                                data=json.dumps(self.config))

            metaData = outFile.create_group("meta")

            metaData.create_dataset("hyperVoxelID", data=self.hyperVoxelID)
            metaData.create_dataset("hyperVoxelOrigo", data=self.hyper_voxel_origo)
            metaData.create_dataset("simulationOrigo", data=self.simulationOrigo)

            metaData.create_dataset("SlurmID", data=self.slurm_id)
            metaData.create_dataset("voxelSize", data=self.voxel_size)
            metaData.create_dataset("hyperVoxelSize", data=self.hyper_voxel_size)
            metaData.create_dataset("nBins", data=self.num_bins)
            metaData.create_dataset("voxelOverflowCounter",
                                    data=self.voxel_overflow_counter)

            metaData.create_dataset("configFile", data=self.config_file)
            metaData.create_dataset("positionFile", data=self.position_file)

            metaData.create_dataset("axonStumpIDFlag", data=self.axon_stump_id_flag)

            # These may or may not exist, if they do, write them to file
            if self.max_axon_voxel_ctr is not None:
                metaData.create_dataset("maxAxonVoxelCtr", data=self.max_axon_voxel_ctr)

            if self.max_dend_voxel_ctr is not None:
                metaData.create_dataset("maxDendVoxelCtr", data=self.max_dend_voxel_ctr)

            if self.voxel_overflow_counter > 0:
                self.write_log("!!! Voxel overflow detected, please increase maxAxon and maxDend")

            networkGroup = outFile.create_group("network")
            networkGroup.create_dataset("synapses",
                                        data=self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr, :],
                                        dtype=np.int32,
                                        chunks=(self.synapse_chunk_size, 13),
                                        maxshape=(None, 13),
                                        compression=self.h5compression)
            networkGroup.create_dataset("gapJunctions",
                                        data=self.hyper_voxel_gap_junctions[:self.hyper_voxel_gap_junction_ctr, :],
                                        dtype=np.int32,
                                        chunks=(self.gap_junction_chunk_size, 11),
                                        maxshape=(None, 11),
                                        compression=self.h5compression)

            networkGroup.create_dataset("synapseLookup",
                                        data=self.hyper_voxel_synapse_lookup,
                                        dtype=int)

            networkGroup.create_dataset("gapJunctionLookup",
                                        data=self.hyper_voxel_gap_gunction_lookup,
                                        dtype=int)

            # Additional information useful for debugging
            if self.debug_flag:
                debugGroup = outFile.create_group("debug")

                debugGroup.create_dataset("dendVoxels", data=self.dend_voxels)
                debugGroup.create_dataset("axonVoxels", data=self.axon_voxels)

                debugGroup.create_dataset("dendVoxelCtr", data=self.dend_voxel_ctr)
                debugGroup.create_dataset("axonVoxelCtr", data=self.axon_voxel_ctr)

            endTime = timeit.default_timer()

            outFile.close()

        self.write_log("Wrote hyper voxel " + str(self.hyperVoxelID) \
                       + " (" + str(self.hyper_voxel_synapse_ctr) + " synapses, "
                       + str(self.hyper_voxel_gap_junction_ctr) + " gap junctions)")

    ############################################################################

    def loadNeuron(self, neuronInfo):

        # Clone prototype neuron (it is centred, and not rotated)
        neuron = self.prototype_neurons[neuronInfo["name"]].clone()

        # Rotate and place neuron in correct location
        neuron.place(rotation=neuronInfo["rotation"],
                     position=neuronInfo["position"])

        return neuron

    ############################################################################

    def distribute_neurons_parallel(self, dView=None):

        if self.role != "master":
            # Only run this as master
            return

        (hyperVoxels, hyperVoxelIDs, nHyperVoxels, simulationOrigo) = \
            self.get_neuron_distribution_history()

        # Do we have old data that we can reuse?
        if hyperVoxels is not None:
            self.write_log("distributeNeuronsParallel: Reusing old neuron allocation")

            self.hyper_voxels = hyperVoxels
            self.hyperVoxelIDs = hyperVoxelIDs
            self.nHyperVoxels = nHyperVoxels
            self.simulationOrigo = simulationOrigo

            # We need to push the data to the workers also
            dView.push({"simulationOrigo": simulationOrigo,
                        "nc.hyperVoxels": hyperVoxels,
                        "nc.hyperVoxelIDs": hyperVoxelIDs,
                        "nc.nHyperVoxels": nHyperVoxels}, block=True)
            return

        # No old data, we need to calculate it

        if dView is None:
            self.write_log("No dView specified, running distribute neurons in serial")
            (minCoord, maxCoord) = self.distributeNeurons()

            self.save_neuron_distribution_history(hyperVoxels=self.hyper_voxels,
                                                  minCoord=minCoord,
                                                  maxCoord=maxCoord)

            return

        (minCoord, maxCoord) = self.findMinMaxCoordParallell(dView=dView,
                                                             volumeID=self.volume_id)

        neuronIdx = np.random.permutation(np.arange(0, len(self.neurons),
                                                    dtype=np.int32))

        # Split the neuronIdx between the workers
        dView.scatter("neuronIdx", neuronIdx, block=True)
        dView.push({"minCoord": minCoord,
                    "maxCoord": maxCoord}, block=True)

        self.write_log("Distributing neurons, parallell.")

        # For the master node, run with empty list
        # This sets up internal state of master
        self.distributeNeurons(neuronIdx=[], minCoord=minCoord, maxCoord=maxCoord)

        cmdStr = "nc.distributeNeurons(neuronIdx=neuronIdx," \
                 + "minCoord=minCoord," \
                 + "maxCoord=maxCoord)"
        dView.execute(cmdStr, block=True)

        self.write_log("Gathering neuron distribution from workers")

        # Collect all the neurons in the list from the workers
        # For each neuron we found out which hyper voxels it occupies,
        # now we want for each hyper voxel to know which neurons are in there
        hyperVoxelList = dView.gather("nc.hyperVoxels", block=True)

        self.write_log("Distributions received.")

        for hv in hyperVoxelList:
            for hID in hv:

                assert (hv[hID]["origo"] == self.hyper_voxels[hID]["origo"]).all(), \
                    "Origo for hyper voxels do not match --- should never happen"

                nNeurons = int(hv[hID]["neuronCtr"])
                startIdx = int(self.hyper_voxels[hID]["neuronCtr"])
                endIdx = startIdx + nNeurons

                if endIdx >= len(self.hyper_voxels[hID]["neurons"]):
                    # Not enough space, reallocating

                    old = self.hyper_voxels[hID]["neurons"]
                    newMax = endIdx + self.max_neurons

                    self.hyper_voxels[hID]["neurons"] = np.zeros((newMax,), dtype=np.int32)

                    # Copying back the old data to new vector
                    if len(old) > 0:
                        self.hyper_voxels[hID]["neurons"][:len(old)] = old

                    del old

                # Adding the new neurons
                self.hyper_voxels[hID]["neurons"][startIdx:endIdx] = \
                    hv[hID]["neurons"][:nNeurons]

                # Increment counter
                self.hyper_voxels[hID]["neuronCtr"] += nNeurons

        # Sorting the list of neurons.
        # -- check why this order matters to number of synapses detected,
        #    it should not matter (except in case of voxel overflows).
        if False:
            for hID in self.hyper_voxels:
                nCtr = self.hyper_voxels[hID]["neuronCtr"]

                self.hyper_voxels[hID]["neurons"] = \
                    np.sort(self.hyper_voxels[hID]["neurons"][:nCtr])

        # Distribute the new list to all neurons
        dView.push({"nc.hyperVoxels": self.hyper_voxels}, block=True)

        self.save_neuron_distribution_history(hyperVoxels=self.hyper_voxels,
                                              minCoord=minCoord,
                                              maxCoord=maxCoord)

    ############################################################################

    # This creates a list for each hyper voxel for the neurons that
    # has any neurites within its border (here defined as vertices inside region)

    def distributeNeurons(self, neuronIdx=None, minCoord=None, maxCoord=None):

        if neuronIdx is None:
            neuronIdx = np.arange(0, len(self.neurons), dtype=np.int32)

        self.write_log("distributeNeurons: neuronIdx = " + str(neuronIdx) \
                       + " (n=" + str(len(neuronIdx)) + ")")

        startTime = timeit.default_timer()

        if maxCoord is None or minCoord is None:
            self.write_log("distributeNeurons: calculating min and max coords")
            (minCoord, maxCoord) = self.findMinMaxCoord()

        # Simulation origo is in meters
        self.simulationOrigo = minCoord

        assert ((self.num_bins - self.num_bins[0]) == 0).all(), \
            "Hyper voxels should be cubes"

        self.hyperVoxelWidth = self.num_bins[0] * self.voxel_size

        self.nHyperVoxels = np.ceil((maxCoord - minCoord)
                                    / self.hyperVoxelWidth).astype(int) + 1

        self.hyperVoxelIDs = np.zeros(self.nHyperVoxels, dtype=int)

        self.hyperVoxelIDs[:] = \
            np.arange(0, self.hyperVoxelIDs.size).reshape(self.hyperVoxelIDs.shape)

        self.write_log(str(self.hyperVoxelIDs.size) + " hyper voxels in total")

        # First assign hyperVoxelID to the space

        self.hyper_voxels = dict([])

        for ix in range(0, self.nHyperVoxels[0]):
            for iy in range(0, self.nHyperVoxels[1]):
                for iz in range(0, self.nHyperVoxels[2]):
                    hID = self.hyperVoxelIDs[ix, iy, iz]

                    self.hyper_voxels[hID] = dict([])
                    self.hyper_voxels[hID]["origo"] = self.simulationOrigo \
                                                      + self.hyperVoxelWidth \
                                                      * np.array([ix, iy, iz])

                    # Changed so we preallocate only empty, to preserve memory
                    self.hyper_voxels[hID]["neurons"] = np.zeros((0,), dtype=np.int32)
                    self.hyper_voxels[hID]["neuronCtr"] = 0

        self.write_log("Pre allocation done.")

        ctr = 0

        if neuronIdx is None:
            neurons = self.neurons
        elif len(neuronIdx) == 0:
            neurons = []
        else:
            neurons = [self.neurons[idx] for idx in neuronIdx]

        for n in neurons:

            ctr = ctr + 1
            if ctr % 10000 == 0:
                print("Assignment counter: " + str(ctr))

            neuron = self.loadNeuron(n)
            neuronID = n["neuronID"]

            if neuron.dend.shape[0] > 0:
                dendLoc = np.floor((neuron.dend[:, :3] - self.simulationOrigo) \
                                   / self.hyperVoxelWidth).astype(int)
            else:
                dendLoc = np.zeros((0, 3))

            if neuron.axon.shape[0] > 0:
                # We have an axon, use it
                axonLoc = np.floor((neuron.axon[:, :3] - self.simulationOrigo) \
                                   / self.hyperVoxelWidth).astype(int)

            elif neuron.axonDensityType == "r":
                # axonLoc = np.zeros((0,3))

                # We create a set of points corresponding approximately to the
                # extent of the axonal density, and check which hyper voxels
                # they occupy

                # Radius of sphere in hyper voxels, rounded up
                rad = np.ceil(neuron.maxAxonRadius \
                              / (self.hyper_voxel_size * self.voxel_size))

                # Approximately how many hyper voxels will the dendritic tree occupy
                nHV = (2 * rad) ** 3

                # Over sample
                nPoints = int(30 * nHV)

                # Randomly place these many points within a sphere of the given radius
                # and then check which hyper voxels these points belong to

                theta = 2 * np.pi * np.random.rand(nPoints)
                phi = np.arccos(2 * np.random.rand(nPoints) - 1)
                r = neuron.maxAxonRadius * (np.random.rand(nPoints) ** (1 / 3))

                x = np.multiply(r, np.multiply(np.sin(phi), np.cos(theta)))
                y = np.multiply(r, np.multiply(np.sin(phi), np.sin(theta)))
                z = np.multiply(r, np.cos(phi))

                axonCloud = np.zeros((len(x), 3))
                axonCloud[:, 0] = x + neuron.soma[0, 0]
                axonCloud[:, 1] = y + neuron.soma[0, 1]
                axonCloud[:, 2] = z + neuron.soma[0, 2]

                axonLoc = np.floor((axonCloud[:, :3] - self.simulationOrigo) \
                                   / self.hyperVoxelWidth).astype(int)

                axonInsideFlag = [xa >= 0 and xa < self.hyperVoxelIDs.shape[0] \
                                  and ya >= 0 and ya < self.hyperVoxelIDs.shape[1] \
                                  and za >= 0 and za < self.hyperVoxelIDs.shape[2] \
                                  for xa, ya, za in axonLoc]

                axonLoc = axonLoc[axonInsideFlag, :]

                # We need to remove the axon volumes outside the modelled volume

                if False:
                    # Verify
                    import matplotlib.pyplot as plt
                    from mpl_toolkits.mplot3d import Axes3D

                    fig = plt.figure()
                    ax = fig.gca(projection='3d')
                    ax.scatter(axonCloud[:, 0], axonCloud[:, 1], axonCloud[:, 2])
                    plt.ion()
                    plt.show()

                    import pdb
                    pdb.set_trace()

                # !!! If there is no axon, and we use probability density,
                #     then we need to include the neuron in hyper voxels
                #     within the axon volume specified

            elif neuron.axonDensityType == "xyz":

                # Estimate how many points we need to randomly place
                nPoints = 100 * np.prod(neuron.axonDensityBoundsXYZ[1:6:2] \
                                        - neuron.axonDensityBoundsXYZ[0:6:2]) \
                          / ((self.hyper_voxel_size * self.voxel_size) ** 3)
                nPoints = int(np.ceil(nPoints))

                if nPoints > 1e4:
                    self.write_log("!!! Many many points placed for axon density of" \
                                   + str(neuron.name) + ": " + str(nPoints))

                xmin = neuron.axonDensityBoundsXYZ[0]
                xwidth = neuron.axonDensityBoundsXYZ[1] - neuron.axonDensityBoundsXYZ[0]
                ymin = neuron.axonDensityBoundsXYZ[2]
                ywidth = neuron.axonDensityBoundsXYZ[3] - neuron.axonDensityBoundsXYZ[2]
                zmin = neuron.axonDensityBoundsXYZ[4]
                zwidth = neuron.axonDensityBoundsXYZ[5] - neuron.axonDensityBoundsXYZ[4]

                # The purpose of this is to find out the range of the axon bounding box
                axonCloud = np.random.rand(nPoints, 3)
                axonCloud[:, 0] = axonCloud[:, 0] * xwidth + xmin
                axonCloud[:, 1] = axonCloud[:, 1] * ywidth + ymin
                axonCloud[:, 2] = axonCloud[:, 2] * zwidth + zmin

                # Dont forget to rotate
                axonCloud = np.matmul(neuron.rotation,
                                      axonCloud.transpose()).transpose() \
                            + neuron.position

                axonLoc = np.floor((axonCloud[:, :3] - self.simulationOrigo) \
                                   / self.hyperVoxelWidth).astype(int)

                axonInsideFlag = [x >= 0 and x < self.hyperVoxelIDs.shape[0] \
                                  and y >= 0 and y < self.hyperVoxelIDs.shape[1] \
                                  and z >= 0 and z < self.hyperVoxelIDs.shape[2] \
                                  for x, y, z in axonLoc]

                axonLoc = axonLoc[axonInsideFlag, :]

                if False:
                    # Verify
                    import matplotlib.pyplot as plt
                    from mpl_toolkits.mplot3d import Axes3D

                    fig = plt.figure()
                    ax = fig.gca(projection='3d')
                    ax.scatter(axonCloud[:, 0], axonCloud[:, 1], axonCloud[:, 2])
                    plt.ion()
                    plt.show()

                    import pdb
                    pdb.set_trace()


            else:
                self.write_log(str(neuron.name) \
                               + ": No axon and unknown axon density type: " \
                               + str(neuron.axonDensityType))
                assert False, "No axon for " + str(neuron.name)

            # Find unique hyper voxel coordinates
            hLoc = np.unique(np.concatenate([axonLoc, dendLoc]), axis=0).astype(int)

            if n["virtualNeuron"]:
                # Range check since we have neurons coming in from outside the volume
                # the parts outside should be ignored
                try:
                    hyperID = [self.hyperVoxelIDs[x, y, z] for x, y, z in hLoc \
                               if x >= 0 and x < self.hyperVoxelIDs.shape[0] \
                               and y >= 0 and y < self.hyperVoxelIDs.shape[1] \
                               and z >= 0 and z < self.hyperVoxelIDs.shape[2]]
                except:
                    self.write_log("Hyper ID problem")
                    assert False, "Hyper ID problem. x=" + str(x) + " y=" \
                                  + str(y) + " z=" + str(z)
                    import pdb
                    pdb.set_trace()
            else:
                # Not a virtual neuron, should all be inside volume
                try:
                    hyperID = [self.hyperVoxelIDs[x, y, z] for x, y, z in hLoc]
                except Exception as e:
                    import traceback
                    tstr = traceback.format_exc()
                    self.write_log(tstr)
                    self.write_log("Affected neuron: " + str(n))
                    self.write_log("Range check fuckup : x = " + str(x) \
                                   + " y = " + str(y) + " z = " + str(z))
                    assert False, "Range check fuckup : x = " + str(x) + " y=" \
                                  + str(y) + " z=" + str(z)
                    import pdb
                    pdb.set_trace()

            # Add the neuron to the hyper voxel's list over neurons
            for hID in hyperID:

                nextPos = self.hyper_voxels[hID]["neuronCtr"]

                if nextPos >= len(self.hyper_voxels[hID]["neurons"]):
                    old = self.hyper_voxels[hID]["neurons"]
                    newMax = nextPos + self.max_neurons
                    self.hyper_voxels[hID]["neurons"] = np.zeros((newMax,), dtype=np.int32)

                    if nextPos > 0:
                        self.hyper_voxels[hID]["neurons"][:len(old)] = old

                    del old

                self.hyper_voxels[hID]["neurons"][nextPos] = neuronID
                self.hyper_voxels[hID]["neuronCtr"] += 1

        endTime = timeit.default_timer()

        if len(neurons) > 0:
            self.write_log("Calculated distribution of neurons: " \
                           + str(endTime - startTime) + " seconds")

        # For serial version of code, we need to return this, so we
        # can save work history
        return minCoord, maxCoord

    ############################################################################

    def setup_parallel(self, dView=None):

        assert self.role == "master", \
            "setupParallel: Should only be called by master node"

        if dView is None:
            self.write_log("setupParallel called without dView, aborting.")
            return

        if self.workers_initialised:
            self.write_log("Workers already initialised.")
            return

        with dView.sync_imports():
            from snudda.detect import SnuddaDetect

        self.write_log("Setting up workers: " \
                       + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

        # Create unique log file names for the workers
        if self.logfile_name is not None:
            engineLogFile = [self.logfile_name + "-" \
                             + str(x) for x in range(0, len(dView))]
        else:
            engineLogFile = [[] for x in range(0, len(dView))]

        self.write_log("Scattering " + str(engineLogFile))

        dView.scatter('logFileName', engineLogFile, block=True)

        self.write_log("Scatter done.")

        dView.push({"positionFile": self.position_file,
                    "configFile": self.config_file,
                    "voxelSize": self.voxel_size,
                    "hyperVoxelSize": self.hyper_voxel_size,
                    "verbose": self.verbose,
                    "SlurmID": self.slurm_id,
                    "saveFile": self.save_file},
                   block=True)

        self.write_log("Init values pushed to workers")

        cmdStr = "nc = SnuddaDetect(configFile=configFile, positionFile=positionFile,voxelSize=voxelSize,hyperVoxelSize=hyperVoxelSize,verbose=verbose,logFileName=logFileName[0],saveFile=saveFile,SlurmID=SlurmID,role='worker')"

        # import pdb
        # pdb.set_trace()

        dView.execute(cmdStr, block=True)

        self.write_log("Workers setup: " \
                       + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

        self.workers_initialised = True

    ############################################################################

    def findMinMaxCoordParallell(self, volumeID=None, dView=None):

        if dView is None:
            self.write_log("findMinMaxCoordParallell: dView is None")
            return self.findMinMaxCoord(volumeID=volumeID)

        self.write_log("Finding min/max coords parallel")

        neuronIdx = np.random.permutation(np.arange(0, len(self.neurons),
                                                    dtype=np.int32))

        dView.scatter("neuronIdxFind", neuronIdx, block=True)
        dView.push({"volumeID": volumeID}, block=True)

        cmdStr = "minMax = nc.findMinMaxCoord(volumeID=volumeID," \
                 + "neuronIdx=neuronIdxFind)"

        dView.execute(cmdStr, block=True)

        self.write_log("Execution of min/max complete")
        # allMinMax = dView.gather("minMax",block=True)
        allMinMax = dView["minMax"]
        self.write_log("Gathered min/max - complete.")

        maxCoord = -1e6 * np.ones((3,))
        minCoord = 1e6 * np.ones((3,))

        try:
            for (minC, maxC) in allMinMax:
                maxCoord = np.maximum(maxCoord, maxC)
                minCoord = np.minimum(minCoord, minC)
        except:
            import traceback
            tstr = traceback.format_exc()
            print(tstr)
            import pdb
            pdb.set_trace()

        return minCoord, maxCoord

    ############################################################################

    def findMinMaxCoord(self, volumeID=None, neuronIdx=None):

        if volumeID is None:
            volumeID = self.volume_id

        print("Finding minMax coord in volumeID = " + str(volumeID))

        maxCoord = -1e6 * np.ones((3,))
        minCoord = 1e6 * np.ones((3,))

        if neuronIdx is None:
            neurons = self.neurons
        else:
            neurons = [self.neurons[idx] for idx in neuronIdx]

        for n in neurons:

            # By using "in" for comparison, we can pass a list of volumeID also
            if volumeID is not None and n["volumeID"] not in volumeID:
                self.write_log("Skipping " + n["name"] \
                               + " when calculating hyper voxel size")
                # Only include neurons belonging to the volume ID
                # we are looking at now
                # import pdb
                # pdb.set_trace()
                continue

            neuron = self.loadNeuron(n)

            if len(neuron.dend) > 0:
                maxCoord = np.maximum(maxCoord, np.max(neuron.dend[:, :3], axis=0))
                minCoord = np.minimum(minCoord, np.min(neuron.dend[:, :3], axis=0))

            if len(neuron.axon) > 0:
                maxCoord = np.maximum(maxCoord, np.max(neuron.axon[:, :3], axis=0))
                minCoord = np.minimum(minCoord, np.min(neuron.axon[:, :3], axis=0))

            maxCoord = np.maximum(maxCoord, np.max(neuron.soma[:, :3], axis=0))
            minCoord = np.minimum(minCoord, np.min(neuron.soma[:, :3], axis=0))

        return minCoord, maxCoord

    ############################################################################

    def fillVoxelsSoma(self, voxelSpace, voxelSpaceCtr,
                       voxelSecID, voxelSecX,
                       somaCoord, neuronID, verbose=False):

        vCoords = np.floor((somaCoord[0, :3] - self.hyper_voxel_origo) \
                           / self.voxel_size).astype(int)
        radius2 = somaCoord[0, 3] ** 2
        vRadius = np.ceil(somaCoord[0, 3] / self.voxel_size).astype(int)

        assert vRadius < 1000, "fillVoxelsSoma: vRadius=" + str(vRadius) \
                               + " soma coords = " + str(somaCoord) + " (BIG SOMA, not SI units?)"

        # Range check, so we stay within hypervoxel
        vxMin = max(0, vCoords[0] - vRadius)
        vxMax = min(self.hyper_voxel_size, vCoords[0] + vRadius + 1)

        vyMin = max(0, vCoords[1] - vRadius)
        vyMax = min(self.hyper_voxel_size, vCoords[1] + vRadius + 1)

        vzMin = max(0, vCoords[2] - vRadius)
        vzMax = min(self.hyper_voxel_size, vCoords[2] + vRadius + 1)

        if verbose:
            print("Soma check x: " + str(vxMin) + " - " + str(vxMax) \
                  + " y: " + str(vyMin) + " - " + str(vyMax) \
                  + " z: " + str(vzMin) + " - " + str(vzMax))

        for vx in range(vxMin, vxMax):
            for vy in range(vyMin, vyMax):
                for vz in range(vzMin, vzMax):

                    d2 = ((vx + 0.5) * self.voxel_size
                          + self.hyper_voxel_origo[0] - somaCoord[0, 0]) ** 2 \
                         + ((vy + 0.5) * self.voxel_size
                            + self.hyper_voxel_origo[1] - somaCoord[0, 1]) ** 2 \
                         + ((vz + 0.5) * self.voxel_size
                            + self.hyper_voxel_origo[2] - somaCoord[0, 2]) ** 2

                    if d2 < radius2:
                        # Mark the point
                        try:
                            vCtr = voxelSpaceCtr[vx, vy, vz]

                            if (vCtr > 0
                                    and voxelSpace[vx, vy, vz, vCtr - 1] == neuronID):
                                # Voxel already has neuronID, skip
                                continue

                            voxelSpace[vx, vy, vz, vCtr] = neuronID
                            voxelSecID[vx, vy, vz, vCtr] = 0  # Soma is 0
                            voxelSecX[vx, vy, vz, vCtr] = 0.5

                            voxelSpaceCtr[vx, vy, vz] += 1
                        except:
                            self.voxel_overflow_counter += 1
                            self.write_log("!!! If you see this you need to increase " \
                                           + "maxDend above " \
                                           + str(voxelSpaceCtr[vx, vy, vz]))
                            continue

    ############################################################################

    # This uses self.hyperVoxelOrigo, self.voxelSize, self.nBins

    # !!! OBS segX must be an integer here, so to get true segX divide by 10000

    def fillVoxelsDend(self, voxelSpace, voxelSpaceCtr,
                       voxelSecID, voxelSecX,
                       voxelSomaDist,
                       coords, links,
                       segID, segX, neuronID):

        # segID gives segment ID for each link
        # segX gives segmentX for each link

        for line, segmentID, segmentX in zip(links, segID, segX):
            p1 = coords[line[0], :3]
            p2 = coords[line[1], :3]
            p1Dist = coords[line[0], 4] * 1e6  # Dist to soma
            p2Dist = coords[line[1], 4] * 1e6

            vp1 = np.floor((p1 - self.hyper_voxel_origo) / self.voxel_size).astype(int)
            vp2 = np.floor((p2 - self.hyper_voxel_origo) / self.voxel_size).astype(int)

            vp1Inside = ((vp1 >= 0).all() and (vp1 < self.num_bins).all())
            vp2Inside = ((vp2 >= 0).all() and (vp2 < self.num_bins).all())

            # Four cases, if neither inside, skip line
            # If one inside but not the other, start at inside point and
            # continue until outside
            # If both inside, add all points without checking if they are inside

            if not vp1Inside and not vp2Inside:
                # No points inside, skip
                continue

            if (vp1 == vp2).all():
                # Line is only one voxel, steps will be 0, so treat it separately
                # We know it is inside, since they are same and both not outside

                vCtr = voxelSpaceCtr[vp1[0], vp1[1], vp1[2]]
                if vCtr > 0 and voxelSpace[vp1[0], vp1[1], vp1[2], vCtr - 1] == neuronID:
                    # Voxel already has neuronID, skip
                    continue

                try:
                    voxelSpace[vp1[0], vp1[1], vp1[2], vCtr] = neuronID
                    voxelSecID[vp1[0], vp1[1], vp1[2], vCtr] = segmentID
                    voxelSecX[vp1[0], vp1[1], vp1[2], vCtr] = segmentX[0]
                    voxelSomaDist[vp1[0], vp1[1], vp1[2], vCtr] = p1Dist

                    voxelSpaceCtr[vp1[0], vp1[1], vp1[2]] += 1

                except:
                    self.voxel_overflow_counter += 1
                    self.write_log("!!! If you see this you need to increase " \
                                   + "maxDend above " \
                                   + str(voxelSpaceCtr[vp1[0], vp1[1], vp1[2]]))
                    continue

                # Done, next voxel
                continue

            if not vp1Inside:
                if not vp2Inside:
                    # No point inside, skip
                    continue
                else:
                    # Start with vp2 continue until outside cube
                    steps = max(np.abs(vp2 - vp1))
                    dv = (vp1 - vp2) / steps
                    ds = (segmentX[0] - segmentX[1]) / steps
                    dd = (p1Dist - p2Dist) / steps

                    # We want the end element "steps" also, hence +1
                    for i in range(0, steps + 1):
                        vp = (vp2 + dv * i).astype(int)
                        sX = segmentX[1] + ds * i  # float
                        somaDist = (p2Dist + dd * i).astype(int)

                        if (vp < 0).any() or (vp >= self.num_bins).any():
                            # Rest of line outside
                            break

                        try:
                            vCtr = voxelSpaceCtr[vp[0], vp[1], vp[2]]
                            if vCtr > 0 and voxelSpace[vp[0], vp[1], vp[2], vCtr - 1] == neuronID:
                                # Voxel already contains neuronID, skip
                                continue

                            voxelSpace[vp[0], vp[1], vp[2], vCtr] = neuronID
                            voxelSecID[vp[0], vp[1], vp[2], vCtr] = segmentID
                            voxelSecX[vp[0], vp[1], vp[2], vCtr] = sX
                            voxelSomaDist[vp[0], vp[1], vp[2], vCtr] = somaDist

                            voxelSpaceCtr[vp[0], vp[1], vp[2]] += 1
                        except:
                            # Increase maxAxon and maxDend
                            self.write_log("!!! If you see this you need to increase " \
                                           + "maxDend above " \
                                           + str(voxelSpaceCtr[vp[0], vp[1], vp[2]]))
                            self.voxel_overflow_counter += 1
                            continue

            elif not vp2Inside:
                # Start with vp1 continue until outside cube
                steps = max(np.abs(vp2 - vp1))
                dv = (vp2 - vp1) / steps
                ds = (segmentX[1] - segmentX[0]) / steps
                dd = (p2Dist - p1Dist) / steps

                # We want the end element "steps" also, hence +1
                for i in range(0, steps + 1):
                    vp = (vp1 + dv * i).astype(int)
                    sX = segmentX[0] + ds * i  # float
                    somaDist = (p1Dist + dd * i).astype(int)

                    if (vp < 0).any() or (vp >= self.num_bins).any():
                        # Rest of line outside
                        break

                    try:
                        vCtr = voxelSpaceCtr[vp[0], vp[1], vp[2]]

                        if vCtr > 0 and voxelSpace[vp[0], vp[1], vp[2], vCtr - 1] == neuronID:
                            # Voxel already contains neuronID, skip
                            continue

                        voxelSpace[vp[0], vp[1], vp[2], vCtr] = neuronID
                        voxelSecID[vp[0], vp[1], vp[2], vCtr] = segmentID
                        voxelSecX[vp[0], vp[1], vp[2], vCtr] = sX
                        voxelSomaDist[vp[0], vp[1], vp[2], vCtr] = somaDist

                        voxelSpaceCtr[vp[0], vp[1], vp[2]] += 1
                    except:
                        self.write_log("!!! If you see this you need to increase " \
                                       + "maxDend above " \
                                       + str(voxelSpaceCtr[vp[0], vp[1], vp[2]]))
                        self.voxel_overflow_counter += 1
                        continue

            else:
                # Entire line inside
                steps = max(np.abs(vp2 - vp1))
                dv = (vp2 - vp1) / steps
                ds = (segmentX[1] - segmentX[0]) / steps
                dd = (p2Dist - p1Dist) / steps

                for i in range(0, steps + 1):
                    vp = (vp1 + dv * i).astype(int)
                    sX = segmentX[0] + ds * i  # float
                    somaDist = (p1Dist + dd * i).astype(int)

                    try:
                        vCtr = voxelSpaceCtr[vp[0], vp[1], vp[2]]

                        if vCtr > 0 and voxelSpace[vp[0], vp[1], vp[2], vCtr - 1] == neuronID:
                            # Voxel already has neuronID, skip
                            continue

                        voxelSpace[vp[0], vp[1], vp[2], vCtr] = neuronID
                        voxelSecID[vp[0], vp[1], vp[2], vCtr] = segmentID
                        voxelSecX[vp[0], vp[1], vp[2], vCtr] = sX
                        voxelSomaDist[vp[0], vp[1], vp[2], vCtr] = somaDist

                        voxelSpaceCtr[vp[0], vp[1], vp[2]] += 1
                    except:
                        self.write_log("!!! If you see this you need to increase " \
                                       + "maxDend above " \
                                       + str(voxelSpaceCtr[vp[0], vp[1], vp[2]]))
                        self.voxel_overflow_counter += 1
                        continue

            # Potentially faster?
            # http://code.activestate.com/recipes/578112-bresenhams-line-algorithm-in-n-dimensions/

    ############################################################################

    def fillVoxelsAxon(self, voxelSpace, voxelSpaceCtr,
                       voxelAxonDist,
                       coords, links,
                       neuronID):

        # segID gives segment ID for each link
        # segX gives segmentX for each link

        for line in links:
            p1 = coords[line[0], :3]
            p2 = coords[line[1], :3]
            p1Dist = coords[line[0], 4] * 1e6  # Dist to soma
            p2Dist = coords[line[1], 4] * 1e6

            vp1 = np.floor((p1 - self.hyper_voxel_origo) / self.voxel_size).astype(int)
            vp2 = np.floor((p2 - self.hyper_voxel_origo) / self.voxel_size).astype(int)

            vp1Inside = ((vp1 >= 0).all() and (vp1 < self.num_bins).all())
            vp2Inside = ((vp2 >= 0).all() and (vp2 < self.num_bins).all())

            # Four cases, if neither inside, skip line
            # If one inside but not the other, start at inside point and
            # continue until outside
            # If both inside, add all points without checking if they are inside

            if not vp1Inside and not vp2Inside:
                # No points inside, skip
                continue

            if (vp1 == vp2).all():
                # Line is only one voxel, steps will be 0, so treat it separately
                # We know it is inside, since they are same and both not outside
                try:
                    vCtr = voxelSpaceCtr[vp1[0], vp1[1], vp1[2]]
                    if vCtr > 0 and voxelSpace[vp1[0], vp1[1], vp1[2], vCtr - 1] == neuronID:
                        # Voxel already has neuronID, skip
                        continue

                    voxelSpace[vp1[0], vp1[1], vp1[2], vCtr] = neuronID
                    voxelAxonDist[vp1[0], vp1[1], vp1[2], vCtr] = p1Dist

                    voxelSpaceCtr[vp1[0], vp1[1], vp1[2]] += 1

                except Exception as e:

                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)

                    self.voxel_overflow_counter += 1
                    self.write_log("!!! If you see this you need to increase " \
                                   + "maxAxon above " \
                                   + str(voxelSpaceCtr[vp1[0], vp1[1], vp1[2]]))
                    continue

                # Done, next voxel
                continue

            if not vp1Inside:
                if not vp2Inside:
                    # No point inside, skip
                    continue
                else:
                    # Start with vp2 continue until outside cube
                    steps = max(np.abs(vp2 - vp1))
                    dv = (vp1 - vp2) / steps
                    dd = (p1Dist - p2Dist) / steps

                    # We want the end element "steps" also, hence +1
                    for i in range(0, steps + 1):
                        vp = (vp2 + dv * i).astype(int)
                        axDist = (p2Dist + dd * i).astype(int)

                        if (vp < 0).any() or (vp >= self.num_bins).any():
                            # Rest of line outside
                            break

                        try:
                            vCtr = voxelSpaceCtr[vp[0], vp[1], vp[2]]
                            if vCtr > 0 and voxelSpace[vp[0], vp[1], vp[2], vCtr - 1] == neuronID:
                                # Voxel already has neuronID, skip
                                continue

                            voxelSpace[vp[0], vp[1], vp[2], vCtr] = neuronID
                            voxelAxonDist[vp[0], vp[1], vp[2], vCtr] = axDist

                            voxelSpaceCtr[vp[0], vp[1], vp[2]] += 1
                        except Exception as e:
                            import traceback
                            tstr = traceback.format_exc()
                            print(tstr)

                            # Increase maxAxon and maxDend
                            self.write_log("!!! If you see this you need to increase " \
                                           + "maxAxon above " \
                                           + str(voxelSpaceCtr[vp[0], vp[1], vp[2]]))
                            self.voxel_overflow_counter += 1
                            continue

            elif not vp2Inside:
                # Start with vp1 continue until outside cube
                steps = max(np.abs(vp2 - vp1))
                dv = (vp2 - vp1) / steps
                dd = (p2Dist - p1Dist) / steps

                # We want the end element "steps" also, hence +1
                for i in range(0, steps + 1):
                    vp = (vp1 + dv * i).astype(int)
                    axDist = (p1Dist + dd * i).astype(int)

                    if (vp < 0).any() or (vp >= self.num_bins).any():
                        # Rest of line outside
                        break

                    try:
                        vCtr = voxelSpaceCtr[vp[0], vp[1], vp[2]]
                        if vCtr > 0 and voxelSpace[vp[0], vp[1], vp[2], vCtr - 1] == neuronID:
                            # Voxel already has neuronID, skip
                            continue

                        voxelSpace[vp[0], vp[1], vp[2], vCtr] = neuronID
                        voxelAxonDist[vp[0], vp[1], vp[2], vCtr] = axDist

                        voxelSpaceCtr[vp[0], vp[1], vp[2]] += 1
                    except Exception as e:

                        import traceback
                        tstr = traceback.format_exc()
                        print(tstr)

                        self.write_log("!!! If you see this you need to increase " \
                                       + "maxAxon above " \
                                       + str(voxelSpaceCtr[vp[0], vp[1], vp[2]]))
                        self.voxel_overflow_counter += 1
                        continue

            else:
                # Entire line inside
                steps = max(np.abs(vp2 - vp1))
                dv = (vp2 - vp1) / steps
                dd = (p2Dist - p1Dist) / steps

                for i in range(0, steps + 1):
                    vp = (vp1 + dv * i).astype(int)
                    axDist = (p1Dist + dd * i).astype(int)

                    try:
                        vCtr = voxelSpaceCtr[vp[0], vp[1], vp[2]]
                        if vCtr > 0 and voxelSpace[vp[0], vp[1], vp[2], vCtr - 1] == neuronID:
                            # Voxel already has neuronID, skip
                            continue

                        voxelSpace[vp[0], vp[1], vp[2], vCtr] = neuronID
                        voxelAxonDist[vp[0], vp[1], vp[2], vCtr] = axDist

                        voxelSpaceCtr[vp[0], vp[1], vp[2]] += 1
                    except Exception as e:
                        import traceback
                        tstr = traceback.format_exc()
                        print(tstr)

                        self.write_log("!!! If you see this you need to increase " \
                                       + "maxAxon above " \
                                       + str(voxelSpaceCtr[vp[0], vp[1], vp[2]]))
                        self.voxel_overflow_counter += 1
                        continue

            # Potentially faster?
            # http://code.activestate.com/recipes/578112-bresenhams-line-algorithm-in-n-dimensions/

    ############################################################################

    def getPath(self, pathStr):

        return pathStr.replace("$DATA", os.path.dirname(__file__) + "/data")

    ############################################################################

    def processHyperVoxel(self, hyperID):

        startTime = timeit.default_timer()

        if self.hyper_voxels[hyperID]["neuronCtr"] == 0:
            # No neurons, return quickly - do not write hdf5 file
            endTime = timeit.default_timer()
            return hyperID, 0, 0, endTime - startTime, 0

        hOrigo = self.hyper_voxels[hyperID]["origo"]
        self.setupHyperVoxel(hOrigo, hyperID)

        nNeurons = self.hyper_voxels[hyperID]["neuronCtr"]

        self.write_log("Processing hyper voxel : " + str(hyperID) \
                       + "/" + str(self.hyperVoxelIDs.size) \
                       + " (" + str(nNeurons) + " neurons)")

        # !!! Suggestion for optimisation. Place neurons with GJ first, then do
        # GJ touch detection, after that add rest of neurons (to get comlete set)
        # and then do axon-dend synapse touch detection

        for neuronID in self.hyper_voxels[hyperID]["neurons"][:nNeurons]:

            neuron = self.loadNeuron(self.neurons[neuronID])

            try:
                self.fillVoxelsAxon(self.axon_voxels,
                                    self.axon_voxel_ctr,
                                    self.axonSomaDist,
                                    neuron.axon,
                                    neuron.axonLinks,
                                    neuronID)

                self.fillVoxelsSoma(self.dend_voxels,
                                    self.dend_voxel_ctr,
                                    self.dendSecID,
                                    self.dendSecX,
                                    neuron.soma,
                                    neuronID)

                self.fillVoxelsDend(self.dend_voxels,
                                    self.dend_voxel_ctr,
                                    self.dendSecID,
                                    self.dendSecX,
                                    self.dendSomaDist,
                                    neuron.dend,
                                    neuron.dendLinks,
                                    neuron.dendSecID,
                                    neuron.dendSecX,
                                    neuronID)

            except Exception as e:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)

                self.write_log("Something went wrong: " + str(tstr))
                import pdb
                pdb.set_trace()

        # This should be outside the neuron loop
        try:
            # This places axon voxels for neurons without axon morphologies
            self.placeSynapsesNoAxon(hyperID,
                                     self.axon_voxels,
                                     self.axon_voxel_ctr,
                                     self.axonSomaDist)
        except Exception as e:
            import traceback
            tstr = traceback.format_exc()
            print(tstr)

            self.write_log("Something wrong: " + str(tstr))
            import pdb
            pdb.set_trace()

        # This detects the synapses where we use a density distribution for axons
        # self.detectSynapsesNoAxonSLOW (hyperID) # --replaced by placeSynapseNoAxon

        # The normal voxel synapse detection
        self.detectSynapses()

        self.detectGapJunctions()

        self.writeHyperVoxelToHDF5()

        endTime = timeit.default_timer()

        return (hyperID, self.hyper_voxel_synapse_ctr,
                self.hyper_voxel_gap_junction_ctr, endTime - startTime,
                self.voxel_overflow_counter)

    ############################################################################

    # hyperID is just needed if we want to plotNeurons also

    def plotHyperVoxel(self, plotNeurons=False, drawAxons=True, drawDendrites=True,
                       detectDone=True):

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        colors = np.zeros((self.dend_voxel_ctr.shape[0],
                           self.dend_voxel_ctr.shape[1],
                           self.dend_voxel_ctr.shape[2], 4))
        colors[:, :, :, 3] = 0.3

        voxelData = np.zeros((self.dend_voxel_ctr.shape[0],
                              self.dend_voxel_ctr.shape[1],
                              self.dend_voxel_ctr.shape[2]))

        if drawAxons:
            colors[:, :, :, 0] = self.axon_voxel_ctr / max(np.max(self.axon_voxel_ctr), 1)
            voxelData += self.axon_voxel_ctr

        if drawDendrites:
            colors[:, :, :, 2] = self.dend_voxel_ctr / max(np.max(self.dend_voxel_ctr), 1)
            voxelData += self.dend_voxel_ctr

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.voxels(voxelData > 0,
                  facecolors=colors, edgecolor=None)

        if self.hyper_voxel_synapse_ctr > 0:
            sCoord = self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr, 2:5]

            # In case hyperVoxelOffset has been applied, we need to subtract it
            # to draw within the hyper voxel
            if self.hyperVoxelOffset is not None:
                sCoord - self.hyperVoxelOffset
            ax.scatter(sCoord[:, 0], sCoord[:, 1], sCoord[:, 2], c="green")

        plt.ion()
        plt.show()

        if plotNeurons:
            # Also plot the neurons overlayed, to verify

            nNeurons = self.hyper_voxels[self.hyperVoxelID]["neuronCtr"]

            for neuronID in self.hyper_voxels[self.hyperVoxelID]["neurons"][:nNeurons]:
                neuron = self.loadNeuron(self.neurons[neuronID])

                neuron.plotNeuron(axis=ax,
                                  plotAxon=drawAxons,
                                  plotDendrite=drawDendrites,
                                  plotOrigo=self.hyper_voxel_origo, plotScale=1 / self.voxel_size)

        plt.show()
        plt.ion()

        plt.pause(0.001)
        figName = "figures/Hypervoxel-" + str(self.slurm_id) \
                  + "-" + str(self.hyperVoxelID) + ".png"
        plt.savefig(figName)

    ############################################################################

    def exportVoxelVisualisationCSV(self, neuronID):

        # x,y,z = coords
        # shape = "cube" or "sphere"
        # type = "axon", "dendrite", "synapse"
        # id = neuronID
        # x,y,z,shape,type,id

        headerStr = "# x,y,z,shape,type,id\n"
        axonStr = ""
        dendStr = ""
        synapseStr = ""

        for x in range(0, self.axon_voxel_ctr.shape[0]):
            for y in range(0, self.axon_voxel_ctr.shape[1]):
                for z in range(0, self.axon_voxel_ctr.shape[2]):
                    for c in range(0, self.axon_voxel_ctr[x, y, z]):
                        nID = self.axon_voxels[x, y, z, c]
                        if nID in neuronID:
                            axonStr += str(x) + "," + str(y) + "," + str(z) \
                                       + ",cube,axon," + str(nID) + "\n"

        for x in range(0, self.dend_voxel_ctr.shape[0]):
            for y in range(0, self.dend_voxel_ctr.shape[1]):
                for z in range(0, self.dend_voxel_ctr.shape[2]):
                    for c in range(0, self.dend_voxel_ctr[x, y, z]):
                        nID = self.dend_voxels[x, y, z, c]
                        if nID in neuronID:
                            dendStr += str(x) + "," + str(y) + "," + str(z) \
                                       + ",cube,dend," + str(nID) + "\n"

        synList = []
        for ir, row in enumerate(self.hyper_voxel_synapses):
            if row[0] in neuronID and row[1] in neuronID:
                synList.append(ir)

        for i in synList:
            xyz = self.hyper_voxel_synapses[i, 2:5]
            synapseStr += str(xyz[0]) + "," + str(xyz[1]) + "," + str(xyz[2]) \
                          + ",sphere,synapse," + str(self.hyper_voxel_synapses[i, 1]) + "\n"

        fName = "hypervoxel-" + str(self.slurm_id) + ".csv"
        with open(fName, 'w') as f:
            f.write(headerStr)
            f.write(axonStr)
            f.write(dendStr)
            f.write(synapseStr)

    ############################################################################

    # Example usage:
    # nc.plotNeuronsInHyperVoxel(neuronID=[1,20],
    #                            neuronColour=np.array([[0,0,1],[0,1,0]]),
    #                            axonAlpha=[1,0.3],dendAlpha=[0.3,1])

    # each row in neuronColour is a colour for a neuron

    def plotNeuronsInHyperVoxel(self, neuronID, neuronColour,
                                axonAlpha=None, dendAlpha=None):

        if axonAlpha is None:
            axonAlpha = np.ones((len(neuronID),))

        if dendAlpha is None:
            dendAlpha = np.ones((len(neuronID),))

        alphaAxonLookup = dict([])
        alphaDendLookup = dict([])
        neuronColourLookup = dict([])

        for ni, aa, da, nc in zip(neuronID, axonAlpha, dendAlpha, neuronColour):
            alphaAxonLookup[ni] = aa
            alphaDendLookup[ni] = da
            neuronColourLookup[ni] = nc

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        colors = np.zeros((self.dend_voxel_ctr.shape[0],
                           self.dend_voxel_ctr.shape[1],
                           self.dend_voxel_ctr.shape[2], 4))

        voxelData = np.zeros((self.dend_voxel_ctr.shape[0],
                              self.dend_voxel_ctr.shape[1],
                              self.dend_voxel_ctr.shape[2]))

        for ix in range(0, self.axon_voxel_ctr.shape[0]):
            for iy in range(0, self.axon_voxel_ctr.shape[1]):
                for iz in range(0, self.axon_voxel_ctr.shape[2]):
                    for ic in range(0, self.axon_voxel_ctr[ix, iy, iz]):
                        nID = self.axon_voxels[ix, iy, iz, ic]
                        if nID in neuronID:
                            colors[ix, iy, iz, 0:3] = neuronColourLookup[nID]
                            colors[ix, iy, iz, 3] = alphaAxonLookup[nID]
                            voxelData[ix, iy, iz] = 1

        for ix in range(0, self.dend_voxel_ctr.shape[0]):
            for iy in range(0, self.dend_voxel_ctr.shape[1]):
                for iz in range(0, self.dend_voxel_ctr.shape[2]):
                    for ic in range(0, self.dend_voxel_ctr[ix, iy, iz]):
                        nID = self.dend_voxels[ix, iy, iz, ic]
                        if nID in neuronID:
                            colors[ix, iy, iz, 0:3] = neuronColourLookup[nID]
                            colors[ix, iy, iz, 3] = alphaDendLookup[nID]
                            voxelData[ix, iy, iz] = 1

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.voxels(voxelData > 0,
                  facecolors=colors, edgecolor=None)

        printList = []
        for ir, row in enumerate(self.hyper_voxel_synapses):
            if row[0] in neuronID and row[1] in neuronID:
                printList.append(ir)

        # This should really only plot those between the neurons indicated
        if len(printList) > 0:
            sCoord = self.hyper_voxel_synapses[printList, 2:5]
            ax.scatter(sCoord[:, 0] + 0.5, sCoord[:, 1] + 0.5, sCoord[:, 2] + 0.5, c="red", s=100)

        plt.axis("off")
        plt.ion()
        plt.show()

        # import pdb
        # pdb.set_trace()

        plt.pause(0.001)
        figName = "figures/Hypervoxel-" + str(self.slurm_id) \
                  + "-" + str(self.hyperVoxelID) + "-someNeurons.png"
        plt.savefig(figName, dpi=900)

    ############################################################################

    def trivialExample(self):

        # This places two neurons in a hyper voxel and tests it
        # --- It destroys the current state of the program.

        self.hyper_voxels[0]["neuronCtr"] = 2
        self.hyper_voxels[0]["origo"] = np.array([0, 0, 0])
        self.neurons[50]["position"] = np.array([-10e-6, 120e-6, -10e-6])
        self.neurons[51]["position"] = np.array([120e-6, 120e-6, -10e-6])
        self.hyper_voxels[0]["neurons"] = np.array([50, 51])

        self.processHyperVoxel(0)

        # self.plotHyperVoxel(plotNeurons=True)
        self.plotHyperVoxel(plotNeurons=False)

        print("Synapses: " + str(self.hyper_voxel_synapse_ctr))
        print(str(self.hyper_voxel_synapses[:self.hyper_voxel_synapse_ctr, :]))

        import pdb
        pdb.set_trace()

    ############################################################################

    def testVoxelDraw(self):

        print("This changes internal state of the object, restart after test run.")

        self.hyperVoxelID = -1
        self.hyper_voxel_origo = np.zeros((3,))
        self.voxel_size = 2
        self.num_bins = np.ones((3, 1)) * 10

        voxels = np.zeros((10, 10, 10, 10), dtype=int)
        voxelCtr = np.zeros((10, 10, 10), dtype=int)
        voxelSecID = np.zeros((10, 10, 10, 10), dtype=int)
        voxelSecX = np.zeros((10, 10, 10, 10), dtype=float)
        voxelSomaDist = np.zeros((10, 10, 10, 10), dtype=int)

        coords = np.array([[2, 2, 2, 1.1, 40], [8, 10, 8, 1.2, 50], [0, 23, 22, 1.3, 60]])
        links = np.array([[0, 1], [0, 2], [2, 1]], dtype=int)
        segID = np.array([1, 2, 3], dtype=int)
        segX = np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]], dtype=float)

        if False:
            self.fillVoxelsDend(voxelSpace=voxels,
                                voxelSpaceCtr=voxelCtr,
                                voxelSecID=voxelSecID,
                                voxelSecX=voxelSecX,
                                voxelSomaDist=voxelSomaDist,
                                coords=coords,
                                links=links,
                                segID=segID,
                                segX=segX,
                                neuronID=13)

        if False:
            self.fillVoxelsSoma(voxelSpace=voxels,
                                voxelSpaceCtr=voxelCtr,
                                voxelSecID=voxelSecID,
                                voxelSecX=voxelSecX,
                                somaCoord=np.array([[10, 10, 10, 8]]),
                                neuronID=14)

        # import pdb
        # pdb.set_trace()

        # We also need to check axon filling

        voxels[:] = 0
        voxelCtr[:] = 0
        voxelSomaDist[:] = 0

        self.fillVoxelsAxon(voxelSpace=voxels,
                            voxelSpaceCtr=voxelCtr,
                            voxelAxonDist=voxelSomaDist,
                            coords=coords,
                            links=links,
                            neuronID=13)

        import pdb
        pdb.set_trace()

    ############################################################################

    # Memory check code taken from
    # https://stackoverflow.com/questions/17718449/determine-free-ram-in-python/17718729#17718729
    #

    def memory(self):
        """
      Get node total memory and memory usage
      """
        try:
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
            return ret
        except:
            return "Non-linux system, /proc/meminfo unavailable"


############################################################################


# @staticmethod
# def _processHyperVoxelHelper(hyperID):
#
#   mem = nc.memory()
#   nc.writeLog("Memory status, before processing " + str(hyperID) \
#               + ": "+ str(mem))
#
#   return nc.processHyperVoxel(hyperID)


############################################################################

def nextRunID():
    runIDfile = ".runID.pickle"

    try:
        if os.path.isfile(runIDfile):

            with open(runIDfile, 'rb') as f:
                runID = pickle.load(f)
                nextID = int(np.ceil(np.max(runID)) + 1)

            runID.append(nextID)

            with open(runIDfile, 'wb') as f:
                pickle.dump(runID, f, -1)

        else:

            with open(runIDfile, 'wb') as f:
                nextID = 1
                runID = [1]
                pickle.dump(runID, f, -1)

    except Exception as e:
        import traceback
        tstr = traceback.format_exc()
        print(tstr)

        print("Problem reading .runID.pickle file, setting runID to 0")
        import pdb
        pdb.set_trace()
        return 0

    print("Using runID = " + str(nextID))

    return nextID


############################################################################

if __name__ == "__main__":
    print("Please do not call this file directly, use snudda.py")
    exit(-1)
