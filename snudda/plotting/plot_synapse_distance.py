# We want to investigate the distribution of glutamate and GABA synapses on the dendrites
#
# How many GABA synapses are there close to each glutamate synapse?
# Does it differ within the dendritic tree?
# Do we see clustering of synapses?
#
import os
from snudda import SnuddaLoad
import h5py


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
        neuron_input_data = self.input["input"][str(neuron_id)]
        neuron_input = dict()

        for input_name in neuron_input_data.keys():
            neuron_input[input_name] = dict()
            neuron_input[input_name]["section_id"] = neuron_input_data[input_name].attrs["SectionID"].copy()
            neuron_input[input_name]["section_x"] = neuron_input_data[input_name].attrs["SectionX"].copy()

            # Do we also want to load the input spikes?  They are in neuron_input_data[input_name]["spikes"][()].copy()




        pass
