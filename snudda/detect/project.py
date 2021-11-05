# In addition to adding synapses using touch detection we also want the ability to add connections
# when we do not have an axon, for example over long range between structures, to connect them together.
# This is what project.py is responsible for.
from collections import OrderedDict

import numpy as np
import json
import os
import h5py

from scipy.interpolate import griddata

from snudda.detect.detect import SnuddaDetect
from snudda.neurons.neuron_morphology import NeuronMorphology
from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.utils.load import SnuddaLoad


class SnuddaProject(object):

    """ Adds projections between neurons, useful for connecting different regions with long range connections. """

    # TODO: Add support for log files!!
    # TODO: Add support for parallel execution
    def __init__(self, network_path, rng=None, random_seed=None, h5libver=None):

        """
        Constructor.

        Args:
            network_path: Path to network directory
            rng: Numpy random stream
            random_seed: Random seed
            h5libver: Version of hdf5 driver to use
        """

        self.network_path = network_path
        self.network_info = None
        self.work_history_file = os.path.join(self.network_path, "log", "network-detect-worklog.hdf5")
        self.output_file_name = os.path.join(self.network_path, "network-projection-synapses.hdf5")

        max_synapses = 100000
        self.synapses = np.zeros((max_synapses, 13), dtype=np.int32)
        self.synapse_ctr = 0
        self.connectivity_distributions = dict()
        self.prototype_neurons = dict()
        self.next_channel_model_id = 10

        self.simulation_origo = None
        self.voxel_size = None

        # Parameters for the HDF5 writing, this affects write speed
        self.synapse_chunk_size = 10000
        self.h5compression = "lzf"

        config_file = os.path.join(self.network_path, "network-config.json")

        with open(config_file, "r") as f:
            self.config = json.load(f, object_pairs_hook=OrderedDict)

        if not h5libver:
            self.h5libver = "latest"
        else:
            self.h5libver = h5libver

        # Setup random generator,
        # TODO: this assumes serial execution. Update for parallel version
        if rng:
            self.rng = rng
        elif random_seed:
            self.rng = np.random.default_rng(random_seed)
        else:
            random_seed = self.config["RandomSeed"]["project"]
            self.rng = np.random.default_rng(random_seed)

        self.read_neuron_positions()
        self.read_prototypes()

    def read_neuron_positions(self):

        """ Reads in neuron positions from network-neuron-positions.hdf5 """

        position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        self.network_info = SnuddaLoad(position_file)

        # We also need simulation origo and voxel size
        work_history_file = os.path.join(self.network_path, "log", "network-detect-worklog.hdf5")
        with h5py.File(work_history_file, "r") as work_hist:
            self.simulation_origo = work_hist["meta/simulationOrigo"][()]
            self.voxel_size = work_hist["meta/voxelSize"][()]

    # This is a simplified version of the prototype load in detect
    def read_prototypes(self):

        """ Reads in neuron prototypes. Simplified version of what same function in detect.py does. """

        for name, definition in self.config["Neurons"].items():

            morph = definition["morphology"]
            param = definition["parameters"]

            if "modulation" in definition:
                modulation = definition["modulation"]
            else:
                modulation = None

            mechanisms = definition["mechanisms"]

            # TODO: Need to update to use NeuronPrototype !!!
            self.prototype_neurons[name] = NeuronPrototype(neuron_name=name,
                                                           neuron_path=None,
                                                           morphology_path=morph,
                                                           parameter_path=param,
                                                           modulation_path=modulation,
                                                           mechanism_path=mechanisms)

        # TODO: The code below is duplicate from detect.py, update so both use same code base
        for name, definition in self.config["Connectivity"].items():

            pre_type, post_type = name.split(",")

            con_def = definition.copy()

            for key in con_def:
                if key == "GapJunction":
                    con_def[key]["channelModelID"] = 3
                else:
                    con_def[key]["channelModelID"] = self.next_channel_model_id
                    self.next_channel_model_id += 1

                # Also if conductance is just a number, add std 0
                if type(con_def[key]["conductance"]) not in [list, tuple]:
                    con_def[key]["conductance"] = [con_def[key]["conductance"], 0]

            self.connectivity_distributions[pre_type, post_type] = con_def

    def project(self, write=True):

        """ Create projections between neurons. """

        for (pre_type, post_type), connection_info in self.connectivity_distributions.items():
            self.connect_projection_helper(pre_type, post_type, connection_info)

        if write:
            self.write()

    def connect_projection_helper(self, pre_neuron_type, post_neuron_type, connection_info):

        """
        Helper function for project.

        Args:
            pre_neuron_type: Type of presynaptic neuron
            post_neuron_type: Type of postsynaptic neuron
            connection_info: Connection info

        """

        for connection_type, con_info in connection_info.items():

            if "projectionFile" not in con_info:
                # Not a projection, skipping
                continue

            projection_file = con_info["projectionFile"]
            with open(projection_file, "r") as f:
                projection_data = json.load(f, object_pairs_hook=OrderedDict)

            if "projectionName" in con_info:
                proj_name = con_info["projectionName"]
                projection_source = np.array(projection_data[proj_name]["source"])*1e-6
                projection_destination = np.array(projection_data[proj_name]["destination"])*1e-6
            else:
                projection_source = np.array(projection_data["source"])*1e-6
                projection_destination = np.array(projection_data["destination"])*1e-6

            if "projectionRadius" in con_info:
                projection_radius = con_info["projectionRadius"]
            else:
                projection_radius = None  # Find the closest neurons

            # TODO: Add projectionDensity later
            # if "projectionDensity" in con_info:
            #     projection_density = con_info["projectionDensity"]
            # else:
            #    projection_density = None  # All neurons within projection radius equally likely

            if "numberOfTargets" in con_info:
                if type(con_info["numberOfTargets"]) == list:
                    number_of_targets = np.array(con_info["numberOfTargets"])  # mean, std
                else:
                    number_of_targets = np.array([con_info["numberOfTargets"], 0])

            if "numberOfSynapses" in con_info:
                if type(con_info["numberOfSynapses"]) == list:
                    number_of_synapses = np.array(con_info["numberOfSynapses"])  # mean, std
                else:
                    number_of_synapses = np.array([con_info["numberOfSynapses"], 0])
            else:
                number_of_synapses = np.array([1, 0])

            if "dendriteSynapseDensity" in con_info:
                dendrite_synapse_density = con_info["dendriteSynapseDensity"]

            if type(con_info["conductance"]) == list:
                conductance_mean, conductance_std = con_info["conductance"]
            else:
                conductance_mean, conductance_std = con_info["conductance"], 0

            # The channelModelID is added to config information
            channel_model_id = con_info["channelModelID"]

            # Find all the presynaptic neurons in the network
            pre_id_list = self.network_info.get_neuron_id_of_type(pre_neuron_type)
            pre_positions = self.network_info.data["neuronPositions"][pre_id_list, :]

            # Find all the postsynaptic neurons in the network
            post_id_list = self.network_info.get_neuron_id_of_type(post_neuron_type)
            post_name_list = [self.network_info.data["name"][x] for x in post_id_list]
            post_positions = self.network_info.data["neuronPositions"][post_id_list, :]

            # For each presynaptic neuron, find their target regions.
            # -- if you want two distinct target regions, you have to create two separate maps
            target_centres = griddata(points=projection_source,
                                      values=projection_destination,
                                      xi=pre_positions, method="linear")

            num_targets = self.rng.normal(number_of_targets[0],
                                          number_of_targets[1],
                                          len(pre_id_list)).astype(int)

            # For each presynaptic neuron, using the supplied map, find the potential post synaptic targets
            for pre_id, centre_pos, n_targets in zip(pre_id_list, target_centres, num_targets):

                d = np.linalg.norm(centre_pos - post_positions, axis=1)  # !! Double check right axis
                d_idx = np.argsort(d)

                if projection_radius:
                    d_idx = d_idx[np.where(d[d_idx] <= projection_radius)[0]]
                    if len(d_idx) > n_targets:
                        d_idx = self.rng.permutation(d_idx)[:n_targets]

                elif len(d_idx) > n_targets:
                    d_idx = d_idx[:n_targets]

                target_id = [post_id_list[x] for x in d_idx]
                target_name = [post_name_list[x] for x in d_idx]
                axon_dist = d[d_idx]

                n_synapses = self.rng.normal(number_of_synapses[0],
                                             number_of_synapses[1],
                                             len(target_id)).astype(int)

                for t_id, t_name, n_syn, ax_dist in zip(target_id, target_name, n_synapses, axon_dist):

                    # We need to place neuron correctly in space (work on clone),
                    # so that synapse coordinates are correct
                    morph_prototype = self.prototype_neurons[t_name]
                    position = self.network_info.data["neurons"][t_id]["position"]
                    rotation = self.network_info.data["neurons"][t_id]["rotation"]
                    parameter_id = self.network_info.data["neurons"][t_id]["parameterID"]
                    morphology_id = self.network_info.data["neurons"][t_id]["morphologyID"]
                    modulation_id = self.network_info.data["neurons"][t_id]["modulationID"]

                    morph = morph_prototype.clone(parameter_id=parameter_id,
                                                  morphology_id=morphology_id,
                                                  modulation_id=modulation_id,
                                                  position=position,
                                                  rotation=rotation)

                    # We are not guaranteed to get n_syn positions, so use len(sec_x) to get how many after
                    # TODO: Fix so dendrite_input_locations always returns  n_syn synapses
                    xyz, sec_id, sec_x, dist_to_soma = morph.dendrite_input_locations(dendrite_synapse_density,
                                                                                      self.rng,
                                                                                      num_locations=n_syn)

                    # We need to convert xyz into voxel coordinates to match data format of synapse matrix
                    xyz = np.round((xyz - self.simulation_origo) / self.voxel_size)

                    # https://www.nature.com/articles/nrn3687 -- lognormal
                    # TODO: Precompute these
                    mu = np.log(conductance_mean ** 2 / np.sqrt(conductance_mean ** 2 + conductance_std ** 2))
                    sigma = np.sqrt(np.log(1 + conductance_std ** 2 / conductance_mean ** 2))

                    cond = self.rng.lognormal(mu, sigma, len(sec_x))
                    cond = np.maximum(cond, conductance_mean*0.1)  # Lower bound, prevent negative.
                    param_id = self.rng.integers(1000000, size=len(sec_x))

                    # TODO: Add code to extend synapses matrix if it is full
                    for i in range(len(sec_id)):
                        self.synapses[self.synapse_ctr, :] = \
                            [pre_id, t_id,
                             xyz[i, 0], xyz[i, 1], xyz[i, 2],
                             -1,  # Hypervoxelid
                             channel_model_id,
                             ax_dist*1e6, dist_to_soma[i]*1e6,
                             sec_id[i], sec_x[i] * 1000,
                             cond[i] * 1e12, param_id[i]]
                        self.synapse_ctr += 1

    def write(self):

        """ Writes projection data to file. """

        # Before writing synapses, lets make sure they are sorted.
        # Sort order: columns 1 (dest), 0 (src), 6 (synapse type)
        sort_idx = np.lexsort(self.synapses[:self.synapse_ctr, [6, 0, 1]].transpose())
        self.synapses[:self.synapse_ctr, :] = self.synapses[sort_idx, :]

        # Write synapses to file
        with h5py.File(self.output_file_name, "w", libver=self.h5libver) as out_file:

            out_file.create_dataset("config", data=json.dumps(self.config))

            network_group = out_file.create_group("network")
            network_group.create_dataset("synapses",
                                         data=self.synapses[:self.synapse_ctr, :],
                                         dtype=np.int32,
                                         chunks=(self.synapse_chunk_size, 13),
                                         maxshape=(None, 13),
                                         compression=self.h5compression)

            network_group.create_dataset("nSynapses", data=self.synapse_ctr, dtype=int)
            network_group.create_dataset("nNeurons", data=self.network_info.data["nNeurons"], dtype=int)

            # This is useful so the merge_helper knows if they need to search this file for synapses
            all_target_id = np.unique(self.synapses[:self.synapse_ctr, 1])
            network_group.create_dataset("allTargetId", data=all_target_id)

            # This creates a lookup that is used for merging later
            synapse_lookup = SnuddaDetect.create_lookup_table(data=self.synapses,
                                                              n_rows=self.synapse_ctr,
                                                              data_type="synapses",
                                                              num_neurons=self.network_info.data["nNeurons"],
                                                              max_synapse_type=self.next_channel_model_id)

            network_group.create_dataset("synapseLookup", data=synapse_lookup)
            network_group.create_dataset("maxChannelTypeID", data=self.next_channel_model_id, dtype=int)

        # We also need to update the work history file with how many synapses we created
        # for the projections between volumes

        with h5py.File(self.work_history_file, "a", libver=self.h5libver) as hist_file:
            if "nProjectionSynapses" in hist_file:
                hist_file["nProjectionSynapses"][()] = self.synapse_ctr
            else:
                hist_file.create_dataset("nProjectionSynapses", data=self.synapse_ctr, dtype=int)
