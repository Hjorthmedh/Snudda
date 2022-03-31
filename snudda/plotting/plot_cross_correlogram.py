import os

import numpy as np
import matplotlib.pyplot as plt
from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation


class PlotCrossCorrelogram:

    def __init__(self, simulation_file):

        self.sim_data = SnuddaLoadNetworkSimulation(network_simulation_output_file=simulation_file)

    def calculate_all_pair_cross_correlogram(self, neuron_id):

        bin_count_total = None
        bin_edges = None

        spike_data = self.sim_data.get_spikes(neuron_id=neuron_id)

        for na in spike_data.keys():

            if spike_data[na].size == 0:
                continue

            for nb in spike_data.keys():
                if na == nb:
                    continue

                if spike_data[nb].size == 0:
                    continue

                bin_count, edges = self.calculate_cross_correlogram(spike_data[na], spike_data[nb])

                if bin_edges is None:
                    bin_edges = edges
                    bin_count_total = bin_count
                else:
                    assert (bin_edges == edges).all()
                    bin_count_total += bin_count

        return bin_count_total, bin_edges

    @staticmethod
    def calculate_cross_correlogram(spike_times_a, spike_times_b, n_bins=101, width=50e-3):

        t_diff = (np.kron(spike_times_a, np.ones(spike_times_b.shape).T)
                  - np.kron(np.ones(spike_times_a.shape), spike_times_b.T)).flatten()

        t_diff = t_diff[np.where(abs(t_diff) <= width)[0]]

        bin_count, bin_edges = np.histogram(t_diff, bins=n_bins, range=[-width, width])
        return bin_count, bin_edges

    def plot_cross_correlogram(self, spike_times_a, spike_times_b, fig_file_name=None):

        bin_count, bin_edges = self.calculate_cross_correlogram(spike_times_a=spike_times_a,
                                                                spike_times_b=spike_times_b)
        plt.figure()
        plt.stairs(values=bin_count, edges=bin_edges)
        plt.xlabel("Time (s)")
        plt.ylabel("Count")
        plt.show()
        if fig_file_name:
            plt.savefig(fig_file_name, dpi=300)

    def plot_all_pair_cross_correlogram(self, neuron_id, fig_file_name=None):

        bin_count, bin_edges = self.calculate_all_pair_cross_correlogram(neuron_id=neuron_id)

        plt.figure()
        plt.stairs(values=bin_count, edges=bin_edges)
        plt.xlabel("Time (s)")
        plt.ylabel("Count")
        plt.show()

        if fig_file_name:
            if not os.path.isdir(os.path.dirname(fig_file_name)):
                print(f"Creating directory {os.path.dirname(fig_file_name)}")
                os.mkdir(os.path.dirname(fig_file_name))

            plt.savefig(fig_file_name, dpi=300)

