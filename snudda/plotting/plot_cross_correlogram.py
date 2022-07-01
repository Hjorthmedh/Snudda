import os

import numpy as np
import matplotlib.pyplot as plt
from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation


class PlotCrossCorrelogram:

    def __init__(self, simulation_file):

        self.sim_data = SnuddaLoadNetworkSimulation(network_simulation_output_file=simulation_file)

    def calculate_all_pair_cross_correlogram(self, neuron_id, time_range=None, shuffle_correct=True):

        bin_count_total = None
        bin_edges = None

        shuffle_count_total = None

        spike_data = self.sim_data.get_spikes(neuron_id=neuron_id)

        for na in spike_data.keys():

            if spike_data[na].size == 0:
                continue

            for nb in spike_data.keys():
                if na == nb:
                    continue

                if spike_data[nb].size == 0:
                    continue

                bin_count, edges = self.calculate_cross_correlogram(spike_data[na],
                                                                    spike_data[nb],
                                                                    time_range=time_range)

                
                shuffle_count, shuffle_edges = self.calculate_cross_correlogram(self.shuffle_spikes(spike_data[na], time_range=time_range),
                                                                                self.shuffle_spikes(spike_data[nb], time_range=time_range),
                                                                                time_range=time_range)

                
                assert (edges == shuffle_edges).all()
                
                if bin_edges is None:
                    bin_edges = edges
                    bin_count_total = bin_count
                else:
                    assert (bin_edges == edges).all()
                    bin_count_total += bin_count

                if shuffle_count_total is None:
                    shuffle_count_total = shuffle_count
                else:
                    shuffle_count_total += shuffle_count

                    
        if shuffle_correct:
            bin_count_total -= shuffle_count_total

        return bin_count_total, bin_edges

    @staticmethod
    def calculate_cross_correlogram(spike_times_a, spike_times_b, n_bins=101, width=50e-3, time_range=None):

        if time_range is not None:
            idx_a = np.where(np.logical_and(time_range[0] <= spike_times_a,
                                            spike_times_a <= time_range[1]))[1]

            idx_b = np.where(np.logical_and(time_range[0] <= spike_times_b,
                                            spike_times_b <= time_range[1]))[1]

            if len(idx_a) == 0 or len(idx_b) == 0:
                t_diff = np.array([])
            else:
                try:
                    t_diff = (np.kron(spike_times_a[:, idx_a], np.ones(spike_times_b[:, idx_b].shape).T)
                              - np.kron(np.ones(spike_times_a[:, idx_a].shape), spike_times_b[:, idx_b].T)).flatten()
                except:
                    import traceback
                    print(traceback.format_exc())
                    import pdb
                    pdb.set_trace()
                
        else:

            t_diff = (np.kron(spike_times_a, np.ones(spike_times_b.shape).T)
                      - np.kron(np.ones(spike_times_a.shape), spike_times_b.T)).flatten()

        t_diff = t_diff[np.where(abs(t_diff) <= width)[0]]

        bin_count, bin_edges = np.histogram(t_diff, bins=n_bins, range=[-width, width])
        return bin_count, bin_edges

    def get_range(self, spike_times, time_range):

        idx = np.where(np.logical_and(time_range[0] <= spike_times.flatten(),
                                      spike_times.flatten() <= time_range[1]))[0]

        return spike_times[:, idx]
        
    
    def shuffle_spikes(self, spike_times, time_range):

        if time_range:
            
            idx = np.where(np.logical_and(time_range[0] <= spike_times.flatten(),
                                          spike_times.flatten() <= time_range[1]))[0]

            if np.size(idx) == 0:
                return np.array([[]])
            
            isi = np.diff(spike_times.flatten()[idx])
            isi = np.append(isi, spike_times.flatten()[idx[0]] - time_range[0])
            st = np.cumsum(np.random.permutation(isi)).reshape(1, (len(idx)))
            st += np.random.uniform(low=time_range[0], high=time_range[0])
            st = np.sort(st % (time_range[1] - time_range[0]) + time_range[0])
            
            return st
            
        isi = np.diff(spike_times)
        isi = np.append(isi, spike_times.flatten()[0])
        
        return np.cumsum(np.random.permutation(isi)).reshape(spike_times.shape)

    def plot_cross_correlogram(self, spike_times_a, spike_times_b, fig_file_name=None):

        bin_count, bin_edges = self.calculate_cross_correlogram(spike_times_a=spike_times_a,
                                                                spike_times_b=spike_times_b)
        plt.rcParams.update({'font.size': 24,
                             'xtick.labelsize': 20,
                             'ytick.labelsize': 20,
                             'legend.loc': 'best'})

        plt.figure()
        plt.stairs(values=bin_count, edges=bin_edges)
        plt.xlabel("Time (s)", fontsize=20)
        plt.ylabel("Count", fontsize=20)
        plt.ion()
        plt.show()

        if fig_file_name:
            plt.savefig(fig_file_name, dpi=300)

    def plot_all_pair_cross_correlogram(self, neuron_id, fig_file_name=None, time_range=None, shuffle_correct=True):

        bin_count, bin_edges = self.calculate_all_pair_cross_correlogram(neuron_id=neuron_id, time_range=time_range,
                                                                         shuffle_correct=shuffle_correct)

        plt.rcParams.update({'font.size': 24,
                             'xtick.labelsize': 20,
                             'ytick.labelsize': 20,
                             'legend.loc': 'best'})

        if bin_count is not None:

            plt.figure()
            plt.stairs(values=bin_count, edges=bin_edges)
            plt.xlabel("Time (s)", fontsize=20)
            plt.ylabel("Count", fontsize=20)
            if time_range is not None:
                plt.title(f"From {time_range[0]}s to {time_range[1]}s")
            plt.show()

            if fig_file_name:
                if not os.path.isdir(os.path.dirname(fig_file_name)):
                    print(f"Creating directory {os.path.dirname(fig_file_name)}")
                    os.mkdir(os.path.dirname(fig_file_name))

                plt.savefig(fig_file_name, dpi=300)

