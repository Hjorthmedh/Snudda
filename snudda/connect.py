# In addition to adding synapses using touch detection we also want the ability to add connections
# when we do not have an axon, for example over long range between structures, to connect them together.
# This is what connect.py is responsible for.

import numpy as np
import scipy
import json
import os


from snudda.neuron_morphology import NeuronMorphology
from snudda.load import SnuddaLoad

class SnuddaConnect(object):

    def __init__(self, network_path, rng=None, random_seed=None):

        self.network_path = network_path
        self.network_info = None
        config_file = os.path.join(self.network_path, "network-config.json")

        with open(config_file, "r") as f:
            self.config = json.load(f)

        # Setup random generator,
        # TODO: this assumes serial execution. Update for parallel version
        if rng:
            self.rng = rng
        elif random_seed:
            self.rng = np.random.default_rng(random_seed)
        else:
            self.rng = np.random.default_rng()

    def read_neuron_positions(self):

        position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        self.network_info = SnuddaLoad(position_file)

    # This is a simplified version of the prototype load in detect
    def read_prototypes(self):

        for name, definition in self.config["Neurons"].items():

            morph = definition["morphology"]

            self.prototype_neurons[name] = NeuronMorphology(name=name, swc_filename=morph)

    #   We want to use: dendrite_input_locations(self, synapse_density, rng, num_locations=None, return_density=False):

    def connect(self):

        # This looks through

        pass

    def connect_projection_helper(self, pre_neuron_type, post_neuron_type, connection_info):

        if "projectionFile" not in connection_info:
            # Not a projection, skipping
            return

        projection_file = connection_info["projectionFile"]
        with open(projection_file, "r") as f:
            projection_data = json.load(f)

        if "projectionName" in connection_info:
            proj_name = connection_info["projectionName"]
            projection_source = np.array(projection_data[proj_name]["source"])*1e-6
            projection_destination = np.array(projection_data[proj_name]["destination"])*1e-6
        else:
            projection_source = np.array(projection_data["source"])*1e-6
            projection_destination = np.array(projection_data["destination"])*1e-6

        if "projectionRadius" in connection_info:
            projection_radius = connection_info["projectionRadius"]
        else:
            projection_radius = None  # Find the closest neurons

        # TODO: Add projectionDensity later
        # if "projectionDensity" in connection_info:
        #     projection_density = connection_info["projectionDensity"]
        # else:
        #    projection_density = None  # All neurons within projection radius equally likely

        if "numberOfTargets" in connection_info:
            if type(connection_info["numberOfTargets"]) == list:
                number_of_targets = np.array(connection_info["numberOfTargets"])  # mean, std
            else:
                number_of_targets = np.array([connection_info["numberOfTargets"], 0])

        if "numberOfSynapses" in connection_info:
            if type(connection_info["numberOfSynapses"]) == list:
                number_of_synapses = np.array(connection_info["numberOfSynapses"])  # mean, std
            else:
                number_of_synapses = np.array([connection_info["numberOfSynapses"], 0])
        else:
            number_of_synapses = np.array([1, 0])

        if "dendriteSynapseDensity" in connection_info:
            dendrite_synapse_density = connection_info["dendriteSynapseDensity"]

        # Find all the presynaptic neurons in the network
        pre_id_list = self.network_info.get_cell_id_of_type(pre_neuron_type)
        pre_positions = self.network_info.data["neuronPositions"][pre_id_list, :]

        # Find all the postsynaptic neurons in the network
        post_id_list = self.network_info.get_cell_id_of_type(post_neuron_type)
        post_name_list = [self.network_info.data["name"][x] for x in post_id_list]
        post_positions = self.network_info.data["neuronPositions"][post_id_list, :]

        # For each presynaptic neuron, find their target regions.
        # -- if you want two distinct target regions, you have to create two separate maps
        target_centres = scipy.interpolate.griddata(points=projection_source,
                                                    values=projection_destination,
                                                    xi=pre_positions, method="linear")

        # For each presynaptic neuron, using the supplied map, find the potential post synaptic targets
        for pre_id, centre_pos in zip(pre_id_list, target_centres):

            d = np.linalg.norm(centre_pos - post_positions, axis=1)  # !! Double check right axis
            d_idx = np.argsort(d)

            if projection_radius:
                d_idx = d_idx[np.where(d[d_idx] <= projection_radius)[0]]
                if len(d_idx) > number_of_targets:
                    d_idx = np.randperm(d_idx)[:nuber_of_targets]

            elif len(d_idx) > number_of_targets:
                d_idx = d_idx[:number_of_targets]

            target_id = post_id_list[d_idx]
            target_name = post_name_list[d_idx]

            n_synapses = self.rng.normal(number_of_synapses[0], number_of_synapses[1], len(target_id))

            for t_id, t_name, n_syn in zip(target_id, target_name, n_synapses):
                morph = self.read_prototypes()[t_name]
                xyz, sec_id, sec_x, dist_to_soma = morph.dendrite_input_locations(dendrite_synapse_density,
                                                                                  self.rng,
                                                                                  num_locations=n_syn)
                # !!! TODO save this synapse

            # Next we need to


            pass

