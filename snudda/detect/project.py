# In addition to adding synapses using touch detection we also want the ability to add connections
# when we do not have an axon, for example over long range between structures, to connect them together.
# This is what project.py is responsible for.
import copy
import json
import os
import traceback
from collections import OrderedDict

import h5py
import numpy as np
from scipy.interpolate import griddata

from snudda.utils.snudda_path import get_snudda_data
from snudda.detect.detect import SnuddaDetect
from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.utils.load import SnuddaLoad


class SnuddaProject(object):
    """ Adds projections between neurons, useful for connecting different regions with long range connections. """

    # TODO: Add support for log files!!
    # TODO: Add support for parallel execution
    def __init__(self, network_path, snudda_data=None, rng=None, random_seed=None, h5libver=None):

        """
        Constructor.

        Args:
            network_path: Path to network directory
            rng: Numpy random stream
            random_seed: Random seed
            h5libver: Version of hdf5 driver to use
        """

        self.network_path = network_path
        self.snudda_data = get_snudda_data(snudda_data=snudda_data,
                                           network_path=self.network_path)

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
            random_seed = self.config["random_seed"]["project"]
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
            self.simulation_origo = work_hist["meta/simulation_origo"][()]
            self.voxel_size = work_hist["meta/voxel_size"][()]

    # This is a simplified version of the prototype load in detect
    def read_prototypes(self):

        """ Reads in neuron prototypes. Simplified version of what same function in detect.py does. """

        for region_name, region_data in self.config["regions"].items():
            for name_type, definition in region_data["neurons"].items():
                for name, neuron_path in definition["neuron_path"].items():

                    self.prototype_neurons[name] = NeuronPrototype(neuron_name=name,
                                                                   neuron_path=neuron_path,
                                                                   snudda_data=self.snudda_data)

        # TODO: The code below is duplicate from detect.py, update so both use same code base
        for region_name, region_data in self.config["regions"].items():
            for name, connection_def in region_data["connectivity"].items():

                pre_type, post_type = name.split(",")

                con_def = copy.deepcopy(connection_def)

                for key in con_def:
                    if key == "gap_junction":
                        con_def[key]["channel_model_id"] = 3
                    else:
                        con_def[key]["channel_model_id"] = self.next_channel_model_id
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

            if "projection_file" not in con_info or "projection_name" not in con_info:
                # Not a projection, skipping
                continue

            projection_file = con_info["projection_file"]
            with open(projection_file, "r") as f:
                projection_data = json.load(f, object_pairs_hook=OrderedDict)

            if "projection_name" in con_info:
                proj_name = con_info["projection_name"]
                projection_source = np.array(projection_data[proj_name]["source"]) * 1e-6
                projection_destination = np.array(projection_data[proj_name]["destination"]) * 1e-6
            else:
                projection_source = np.array(projection_data["source"]) * 1e-6
                projection_destination = np.array(projection_data["destination"]) * 1e-6

            if "projection_radius" in con_info:
                projection_radius = con_info["projection_radius"]
            else:
                projection_radius = None  # Find the closest neurons

            # TODO: Add projection_density later
            # if "projection_density" in con_info:
            #     projection_density = con_info["projection_density"]
            # else:
            #    projection_density = None  # All neurons within projection radius equally likely

            if "number_of_targets" in con_info:
                if type(con_info["number_of_targets"]) == list:
                    number_of_targets = np.array(con_info["number_of_targets"])  # mean, std
                else:
                    number_of_targets = np.array([con_info["number_of_targets"], 0])

            if "number_of_synapses" in con_info:
                if type(con_info["number_of_synapses"]) == list:
                    number_of_synapses = np.array(con_info["number_of_synapses"])  # mean, std
                else:
                    number_of_synapses = np.array([con_info["number_of_synapses"], 0])
            else:
                number_of_synapses = np.array([1, 0])

            if "dendrite_synapse_density" in con_info:
                dendrite_synapse_density = con_info["dendrite_synapse_density"]

            if type(con_info["conductance"]) == list:
                conductance_mean, conductance_std = con_info["conductance"]
            else:
                conductance_mean, conductance_std = con_info["conductance"], 0

            # The channel_model_id is added to config information
            channel_model_id = con_info["channel_model_id"]

            # Find all the presynaptic neurons in the network
            pre_id_list = self.network_info.get_neuron_id_of_type(pre_neuron_type)
            pre_positions = self.network_info.data["neuron_positions"][pre_id_list, :]

            # Find all the postsynaptic neurons in the network
            post_id_list = self.network_info.get_neuron_id_of_type(post_neuron_type)
            post_name_list = [self.network_info.data["name"][x] for x in post_id_list]
            post_positions = self.network_info.data["neuron_positions"][post_id_list, :]

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

                n_synapses = np.maximum(0, self.rng.normal(number_of_synapses[0],
                                                           number_of_synapses[1],
                                                           len(target_id))).astype(int)


                for t_id, t_name, n_syn, ax_dist in zip(target_id, target_name, n_synapses, axon_dist):

                    # We need to place neuron correctly in space (work on clone),
                    # so that synapse coordinates are correct
                    morph_prototype = self.prototype_neurons[t_name]
                    position = self.network_info.data["neurons"][t_id]["position"]
                    rotation = self.network_info.data["neurons"][t_id]["rotation"]
                    parameter_key = self.network_info.data["neurons"][t_id]["parameter_key"]
                    morphology_key = self.network_info.data["neurons"][t_id]["morphology_key"]
                    modulation_key = self.network_info.data["neurons"][t_id]["modulation_key"]

                    morph = morph_prototype.clone(parameter_key=parameter_key,
                                                  morphology_key=morphology_key,
                                                  modulation_key=modulation_key,
                                                  position=position,
                                                  rotation=rotation)

                    # We are not guaranteed to get n_syn positions, so use len(sec_x) to get how many after
                    # TODO: Fix so dendrite_input_locations always returns  n_syn synapses
                    xyz, sec_id, sec_x, dist_to_soma = morph.dendrite_input_locations(synapse_density_str=dendrite_synapse_density,
                                                                                      rng=self.rng,
                                                                                      num_locations=n_syn)

                    # We need to convert xyz into voxel coordinates to match data format of synapse matrix
                    xyz = np.round((xyz - self.simulation_origo) / self.voxel_size)

                    # https://www.nature.com/articles/nrn3687 -- lognormal
                    # TODO: Precompute these
                    mu = np.log(conductance_mean ** 2 / np.sqrt(conductance_mean ** 2 + conductance_std ** 2))
                    sigma = np.sqrt(np.log(1 + conductance_std ** 2 / conductance_mean ** 2))

                    cond = self.rng.lognormal(mu, sigma, len(sec_x))
                    cond = np.maximum(cond, conductance_mean * 0.1)  # Lower bound, prevent negative.
                    param_id = self.rng.integers(1000000, size=len(sec_x))

                    # TODO: Add code to extend synapses matrix if it is full
                    for i in range(len(sec_id)):
                        self.synapses[self.synapse_ctr, :] = \
                            [pre_id, t_id,
                             xyz[i, 0], xyz[i, 1], xyz[i, 2],
                             -1,  # Hypervoxelid
                             channel_model_id,
                             ax_dist * 1e6, dist_to_soma[i] * 1e6,
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

            network_group.create_dataset("num_synapses", data=self.synapse_ctr, dtype=int)
            network_group.create_dataset("num_neurons", data=self.network_info.data["num_neurons"], dtype=int)

            # This is useful so the merge_helper knows if they need to search this file for synapses
            all_target_id = np.unique(self.synapses[:self.synapse_ctr, 1])
            network_group.create_dataset("all_target_id", data=all_target_id)

            # This creates a lookup that is used for merging later
            synapse_lookup = SnuddaDetect.create_lookup_table(data=self.synapses,
                                                              n_rows=self.synapse_ctr,
                                                              data_type="synapses",
                                                              num_neurons=self.network_info.data["num_neurons"],
                                                              max_synapse_type=self.next_channel_model_id)

            network_group.create_dataset("synapse_lookup", data=synapse_lookup)
            network_group.create_dataset("max_channel_type_id", data=self.next_channel_model_id, dtype=int)

        # We also need to update the work history file with how many synapses we created
        # for the projections between volumes

        with h5py.File(self.work_history_file, "a", libver=self.h5libver) as hist_file:
            if "num_projection_synapses" in hist_file:
                hist_file["num_projection_synapses"][()] = self.synapse_ctr
            else:
                hist_file.create_dataset("num_projection_synapses", data=self.synapse_ctr, dtype=int)
