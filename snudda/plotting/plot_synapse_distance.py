# We want to investigate the distribution of glutamate and GABA synapses on the dendrites
#
# How many GABA synapses are there close to each glutamate synapse?
# Does it differ within the dendritic tree?
# Do we see clustering of synapses? -- This can be even more interesting to investigate after we have applied plasticity
#
#

import numpy as np
import os
from snudda import SnuddaLoad
import h5py

from snudda.neurons.morphology_data import MorphologyData


class PlotSynapseDistance:

    def __init__(self, network_path, network_file=None, input_file=None):

        self.network_path = network_path

        if network_file:
            self.network_file = network_file
        else:
            self.network_file = os.path.join(network_path, "network-synapses.hdf5")

        if input_file:
            self.input_file = input_file
        else:
            self.input_file = os.path.join(network_path, "input-spikes.hdf5")

        self.snudda_load = SnuddaLoad(network_file=network_file)
        self.input = h5py.File(self.input_file, "r")

    def gather_synapses(self, neuron_id):

        # Find the synapses from the connection matrix
        synapses = self.snudda_load.find_synapses(post_id=neuron_id)

        # Find input
        external_input_data = self.input["input"][str(neuron_id)]
        external_input = dict()

        for input_name in external_input_data.keys():
            external_input[input_name] = dict()
            external_input[input_name]["section_id"] = external_input_data[input_name].attrs["SectionID"].copy()
            external_input[input_name]["section_x"] = external_input_data[input_name].attrs["SectionX"].copy()

            # Do we also want to load the input spikes?  They are in neuron_input_data[input_name]["spikes"][()].copy()

        return synapses, external_input

    def load_morphology(self, neuron_id):

        morphology_path = self.snudda_load.data["neurons"][neuron_id]["morphology"]
        morphology_data = MorphologyData(swc_file=morphology_path)

        return morphology_data

    def get_synapse_neighbours(self, neuron_id, input_name,
                               soma_distance_range=None,
                               max_neighbour_distance=30e-6):

        # This function find the distances to the surrounding synapses from different external
        # input and synapse types

        # 1. Loop through all the sections (and note the distance to soma)
        # 2. Find the synapses belonging to input_name
        # 3. Find all parent and child comparments.
        # 4. Find all external input synapses (include of input_name) in the current section_id
        #    and in all the parent and child compartments
        # 5. Calculate distance to the input_name synapses in that section id.

        morphology_data = self.load_morphology(neuron_id=neuron_id)

        synapses, external_input = self.gather_synapses(neuron_id=neuron_id)

        neighbour_synapses = dict()

        for section in morphology_data.section_iterator(section_type=3):

            # Find the input locations
            idx = np.where(external_input[input_name]["section_id"] == section.section_id)[0]
            input_section_x = external_input[input_name]["section_x"][idx]
            input_soma_distance = section.soma_distance_at(input_section_x)
            input_soma_distance_all = input_soma_distance

            if soma_distance_range:
                keep_idx = np.where(np.logical_and(soma_distance_range[0] <= input_soma_distance,
                                                   input_soma_distance < soma_distance_range[1]))[0]
                input_section_x = input_section_x[keep_idx]
                input_soma_distance = input_soma_distance[keep_idx]

            # Now find all the neighbouring synapses of different types!
            for neighbour_input_name in external_input:
                neighbour_synapses[neighbour_input_name] = []

                # First we compare with synapses in the section itself
                neigh_dist = np.repeat(input_soma_distance.reshape(len(input_soma_distance_all, 1)),
                                       repeats=len(input_soma_distance_all),
                                       axis=1) \
                            - np.repeat(input_soma_distance_all.reshape(1, len(input_soma_distance)),
                                        repeats=len(input_soma_distance),
                                        axis=0)

                # !!! TODO: Work in progress

                # Then in the parent section

                # The in the children sections

            neighbour_section_id = [section.section_id, section.parent_section_id] + list(section.child_section_id)

            neighbour_idx = np.where(external_input[input_name]["section_id"] in neighbour_section_id)[0]
            neighbour_section_x = external_input[input_name]["section_x"][neighbour_idx]
            neighbour_soma_distance = se

            parent_section = morphology_data.section_data[3][section.parent_section_id]
            child_sections = [morphology_data.section_data[3][s_id] for s_id in section.child_section_id]



        pass
