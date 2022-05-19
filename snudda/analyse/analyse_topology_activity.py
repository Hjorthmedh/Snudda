import os
import json
import numpy as np

from snudda.utils.load import SnuddaLoad
from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation
from snudda.utils.export_connection_matrix import SnuddaExportConnectionMatrix
from collections import OrderedDict

import matplotlib.pyplot as plt


class SnuddaAnalyseTopologyActivity:

    def __init__(self):

        self.simulation_data = dict()
        self.mapping_list = dict()
        self.mapping_dictionary = dict()

    def load_simulation_data(self, data_key, simulation_output=None):
        self.simulation_data[data_key] = SnuddaLoadNetworkSimulation(network_simulation_output_file=simulation_output)
        self.load_mapping_file(data_key)

    def load_mapping_file(self, data_key):

        network_file = SnuddaLoad.to_str(self.simulation_data[data_key].network_simulation_file["metaData"]["networkFile"][()])
        mapping_file = f"{network_file}-remapping.txt"

        self.mapping_list[data_key] = np.genfromtxt(mapping_file, delimiter=',', dtype=int)
        self.mapping_dictionary[data_key] = OrderedDict()

        for row in self.mapping_list[data_key]:
            self.mapping_dictionary[data_key][row[0]] = row[1]

    def check_same_neurons(self, data_key_a, data_key_b):

        sim_a = self.simulation_data[data_key_a]
        sim_b = self.simulation_data[data_key_b]

        is_same = len(list(sim_a.iter_neuron_id())) == len(list(sim_b.iter_neuron_id()))

        for nid_A, nid_B in zip(sim_a.iter_neuron_id(), sim_b.iter_neuron_id()):
            is_same = is_same and nid_A == nid_B
            is_same = is_same and sim_a.get_neuron_keys(nid_A) == sim_b.get_neuron_keys(nid_B)

        return is_same

    def get_spike_deltas(self, data_key_a, data_key_b, match_closest=True):

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

            if s_a.size > 0 and s_b.size > 0:

                if match_closest:
                    t_diff = np.kron(np.ones(s_a.shape), s_b.T) - np.kron(s_a, np.ones(s_b.T.shape))
                    min_pos_a = np.argmin(np.abs(t_diff), axis=1)
                    min_pos_b = np.argmin(np.abs(t_diff), axis=0)

                    t_min_diff_a = [t_diff[m[0], m[1]] for m in zip(range(len(min_pos_b)), min_pos_a)]
                    t_min_diff_b = [-t_diff[m[0], m[1]] for m in zip(min_pos_b, range(len(min_pos_a)))]

                    spike_time_difference[neuron_id] = np.array(t_min_diff_a), np.array(t_min_diff_b)
                else:
                    n_compare = min(s_a.size, s_b.size)
                    spike_time_difference[neuron_id] = s_b[0, :n_compare] - s_a[0, :n_compare]
            else:
                # At least one of the spike trains does not have any spikes
                spike_time_difference[neuron_id] = np.array([]), np.array([])

        return spike_time_difference

    def plot_spike_delta_histogram(self, data_key_a=None, data_key_b=None,
                                   plot_title=None, direction=0,
                                   match_closest=True,
                                   range_min=-10e-3, range_max=2e-3, bin_size=0.5e-3):

        spike_time_difference = self.get_spike_deltas(data_key_a=data_key_a, data_key_b=data_key_b,
                                                      match_closest=match_closest)
        n_bins = int(np.ceil((range_max-range_min) / bin_size)) + 1

        fig = plt.figure()
        neuron_names = self.simulation_data[data_key_a].get_neuron_name()
        plt.set_cmap("autumn")

        hist_data = []
        hist_label = []

        for nid in spike_time_difference.keys():
            hist_data.append(spike_time_difference[nid][direction])
            hist_label.append(f"{neuron_names[nid]} ({nid})")

        plt.hist(hist_data, range=(range_min, range_max), bins=n_bins,
                 label=hist_label, histtype="barstacked")

            #plt.ion()
            #plt.show()
            #import pdb
            #pdb.set_trace()

        plt.xlabel("Time difference (s)")
        plt.ylabel("Count")
        plt.title(plot_title)
        plt.legend()
        plt.ion()
        plt.show()