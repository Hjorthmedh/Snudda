import os
import numpy as np

from snudda.utils.load import SnuddaLoad
from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation
from snudda.utils.export_connection_matrix import SnuddaExportConnectionMatrix
from collections import OrderedDict


class SnuddaAnalyseTopologyActivity:

    def __init__(self):

        self.simulation_data = dict()

    def load_simulation_data(self, data_key, simulation_output=None):
        self.simulation_data[data_key] = SnuddaLoadNetworkSimulation(network_simulation_output_file=simulation_output)

    def check_same_neurons(self, data_key_a, data_key_b):

        sim_a = self.simulation_data[data_key_a]
        sim_b = self.simulation_data[data_key_b]

        is_same = len(list(sim_a.iter_neuron_id())) == len(list(sim_b.iter_neuron_id()))

        for nid_A, nid_B in zip(sim_a.iter_neuron_id(), sim_b.iter_neuron_id()):
            is_same = is_same and nid_A == nid_B
            is_same = is_same and sim_a.get_neuron_keys(nid_A) == sim_b.get_neuron_keys(nid_B)

        return is_same

    def get_spike_deltas(self, data_key_a, data_key_b):

        # Check that the neurons compared are the same (by verifying parameter key, morphology key, modulation key)
        assert self.check_same_neurons(data_key_a, data_key_b), f"data_keys have different neurons in the network"

        # Match spikes against each other, compute change...

        sim_a = self.simulation_data[data_key_a]
        sim_b = self.simulation_data[data_key_b]

        spikes_a = sim_a.get_spikes()
        spikes_b = sim_b.get_spikes()

        spike_time_difference = dict()

        for neuron_id in spikes_a.keys():
            s_a = spikes_a[neuron_id]
            s_b = spikes_b[neuron_id]

            t_diff = np.kron(np.ones(s_a.shape), s_b.T) - np.kron(s_a, np.ones(s_b.T.shape))
            min_pos_a = np.argmin(np.abs(t_diff), axis=0)
            min_pos_b = np.argmin(np.abs(t_diff), axis=1)

            t_min_diff_a = [t_diff[m[0], m[1]] for m in zip(min_pos_a, range(len(min_pos_a)))]
            t_min_diff_b = [-t_diff[m[0], m[1]] for m in zip(min_pos_b, range(len(min_pos_b)))]

            spike_time_difference[neuron_id] = np.array(t_min_diff_a), np.array(t_min_diff_b)

        return spike_time_difference
