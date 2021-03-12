# In addition to adding synapses using touch detection we also want the ability to add connections
# when we do not have an axon, for example over long range between structures, to connect them together.
# This is what connect.py is responsible for.

import numpy as np
import scipy
import json
import os
import h5py

from scipy.interpolate import griddata
from snudda.neuron_morphology import NeuronMorphology
from snudda.load import SnuddaLoad


class SnuddaConnect(object):

    def __init__(self, network_path, rng=None, random_seed=None, h5libver=None):

        self.network_path = network_path
        self.network_info = None
        max_synapses = 100000
        self.synapses = np.zeros((max_synapses, 13), dtype=np.int32)
        self.synapse_ctr = 0
        self.connectivity_distributions = dict()
        self.prototype_neurons = dict()
        self.next_channel_model_id = 10

        # Parameters for the HDF5 writing, this affects write speed
        self.synapse_chunk_size = 10000
        self.h5compression = "lzf"

        config_file = os.path.join(self.network_path, "network-config.json")

        with open(config_file, "r") as f:
            self.config = json.load(f)

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
            self.rng = np.random.default_rng()

        self.read_neuron_positions()
        self.read_prototypes()

    def read_neuron_positions(self):

        position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        self.network_info = SnuddaLoad(position_file)

    # This is a simplified version of the prototype load in detect
    def read_prototypes(self):

        for name, definition in self.config["Neurons"].items():

            morph = definition["morphology"]

            self.prototype_neurons[name] = NeuronMorphology(name=name, swc_filename=morph)

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

    def connect(self):

        for (pre_type, post_type), connection_info in self.connectivity_distributions.items():
            print(f"pre {pre_type}, post {post_type}")
            self.connect_projection_helper(pre_type, post_type, connection_info)

    def connect_projection_helper(self, pre_neuron_type, post_neuron_type, connection_info):

        for connection_type, con_info in connection_info.items():

            if "projectionFile" not in con_info:
                # Not a projection, skipping
                continue

            projection_file = con_info["projectionFile"]
            with open(projection_file, "r") as f:
                projection_data = json.load(f)

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
            pre_id_list = self.network_info.get_cell_id_of_type(pre_neuron_type)
            pre_positions = self.network_info.data["neuronPositions"][pre_id_list, :]

            # Find all the postsynaptic neurons in the network
            post_id_list = self.network_info.get_cell_id_of_type(post_neuron_type)
            post_name_list = [self.network_info.data["name"][x] for x in post_id_list]
            post_positions = self.network_info.data["neuronPositions"][post_id_list, :]

            # For each presynaptic neuron, find their target regions.
            # -- if you want two distinct target regions, you have to create two separate maps
            target_centres = griddata(points=projection_source,
                                      values=projection_destination,
                                      xi=pre_positions, method="linear")

            num_targets = self.rng.normal(number_of_targets[0], number_of_targets[1], len(pre_id_list))

            # For each presynaptic neuron, using the supplied map, find the potential post synaptic targets
            for pre_id, centre_pos, n_targets in zip(pre_id_list, target_centres, num_targets):

                d = np.linalg.norm(centre_pos - post_positions, axis=1)  # !! Double check right axis
                d_idx = np.argsort(d)

                if projection_radius:
                    d_idx = d_idx[np.where(d[d_idx] <= projection_radius)[0]]
                    if len(d_idx) > n_targets:
                        d_idx = np.randperm(d_idx)[:n_targets]

                elif len(d_idx) > n_targets:
                    d_idx = d_idx[:n_targets]

                target_id = post_id_list[d_idx]
                target_name = post_name_list[d_idx]
                axon_dist = d[d_idx]

                n_synapses = self.rng.normal(number_of_synapses[0], number_of_synapses[1], len(target_id))

                for t_id, t_name, n_syn, ax_dist in zip(target_id, target_name, n_synapses, axon_dist):
                    morph = self.read_prototypes()[t_name]

                    # We are not guaranteed to get n_syn positions, so use len(sec_x) to get how many after
                    xyz, sec_id, sec_x, dist_to_soma = morph.dendrite_input_locations(dendrite_synapse_density,
                                                                                      self.rng,
                                                                                      num_locations=n_syn)

                    cond = self.rng.normal(conductance_mean, conductance_std, len(sec_x))
                    cond = np.maximum(cond, conductance_mean*0.1)  # Lower bound, prevent negative.
                    param_id = self.hyper_voxel_rng.integers(1000000, size=len(sec_x))

                    # TODO: Add code to extend synapses matrix if it is full
                    for i in range(len(sec_id)):
                        self.synapses[self.synapse_ctr, :] = \
                            [pre_id, t_id,
                             xyz[i, 0], xyz[i, 1], xyz[i, 2],
                             np.nan,  # Hypervoxelid
                             channel_model_id,
                             ax_dist, dist_to_soma[i],
                             sec_id[i], sec_x[i] * 1000,
                             cond[i] * 1e12, param_id[i]]
                        self.synapse_ctr += 1

    def write(self):

        output_file_name = os.path.join(self.network_path, "network-connection-synapses.hdf5")
        with h5py.File(output_file_name, "w", libver=self.h5libver) as out_file:

            out_file.create_dataset("config", data=json.dumps(self.config))

            network_group = out_file.create_group("network")
            network_group.create_dataset("synapses",
                                         data=self.synapses[:self.synapse_ctr, :],
                                         dtype=np.int32,
                                         chunks=(self.synapse_chunk_size, 13),
                                         maxshape=(None, 13),
                                         compression=self.h5compression)


if __name__ == "__main__":

    print("This main function is only here for testing that things work.")

    sc = SnuddaConnect(network_path=os.path.join("networks", "con_test"))
