import os
import json
import numpy as np

from snudda.utils.load import SnuddaLoad
from snudda.utils.load_network_simulation import SnuddaLoadSimulation
from snudda.utils.export_connection_matrix import SnuddaExportConnectionMatrix
from collections import OrderedDict

import matplotlib.pyplot as plt


class SnuddaAnalyseTopologyActivity:

    def __init__(self):

        self.simulation_data = dict()
        self.mapping_list = dict()
        self.mapping_dictionary = dict()

    def load_simulation_data(self, data_key, simulation_output=None):
        self.simulation_data[data_key] = SnuddaLoadSimulation(network_simulation_output_file=simulation_output)
        self.load_mapping_file(data_key)

    def load_mapping_file(self, data_key):

        network_file = SnuddaLoad.to_str(self.simulation_data[data_key].network_simulation_file["meta_data"]["network_file"][()])
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

    def match_closest_spikes(self, spike_train_a, spike_train_b):

        # This function assumes there are spike in both spike trains
        t_diff = (np.kron(np.ones((spike_train_a.size, 1)), spike_train_b)
                  - np.kron(spike_train_a.reshape(spike_train_a.size, 1), np.ones(spike_train_b.shape)))
        min_pos_a = np.argmin(np.abs(t_diff), axis=1)
        min_pos_b = np.argmin(np.abs(t_diff), axis=0)

        t_min_diff_a = [t_diff[m[0], m[1]] for m in zip(range(len(min_pos_b)), min_pos_a)]
        t_min_diff_b = [-t_diff[m[0], m[1]] for m in zip(min_pos_b, range(len(min_pos_a)))]

        return np.array(t_min_diff_a), np.array(t_min_diff_b)

    def match_order_spikes(self, spike_train_a, spike_train_b):

        n_compare = min(spike_train_a.size, spike_train_b.size)
        return spike_train_b[:n_compare] - spike_train_a[:n_compare], \
            spike_train_a[:n_compare] - spike_train_b[:n_compare]

    def match_closest_unique(self, spike_train_a, spike_train_b, delta_t):

        dt_list_a = np.full(spike_train_a.shape, np.nan, dtype=float)
        dt_list_b = np.full(spike_train_b.shape, np.nan, dtype=float)

        i_b = 0

        for i_a, t_a in enumerate(spike_train_a):

            while i_b < len(spike_train_b) and spike_train_b[i_b] - t_a < delta_t:
                if spike_train_b[i_b] - t_a < -delta_t:
                    i_b += 1
                else:
                    dt_list_a[i_a] = spike_train_b[i_b] - t_a
                    i_b += 1
                    break

        i_a = 0

        for i_b, t_b in enumerate(spike_train_b):

            while i_a < len(spike_train_a) and spike_train_a[i_a] - t_b < delta_t:
                if spike_train_a[i_a] - t_b < -delta_t:
                    i_a += 1
                else:
                    dt_list_b[i_b] = spike_train_a[i_a] - t_b
                    i_a += 1
                    break

        return dt_list_a, dt_list_b

    def get_spike_triggered_deltas(self, spike_times_a, spike_times_b):

        # The idea here is that for each spike in spike_train_a, we want to find the delta_t to all the spikes
        # following in spike_train_b (but before the next spike in spike_train_a). Sort of like a JPSTH.

        spike_dt = np.full(spike_times_b.shape, np.nan)

        idx_a = len(spike_times_a) - 1

        for idx_b, t_b in zip(np.flip(np.arange(0, len(spike_times_b))), np.flip(spike_times_b)):
            while idx_a >= 0 and spike_times_a[idx_a] > t_b:
                idx_a -= 1

            spike_dt[idx_b] = t_b - spike_times_a[idx_a] if idx_a >= 0 else np.nan

        return spike_dt

    def get_triggered_spikes(self, spike_times, trigger_times, duration, subtract_trigger=True):

        start_idx = 0

        for trigger_t in trigger_times:
            while start_idx < len(spike_times) and spike_times[start_idx] < trigger_t:
                start_idx += 1

            end_idx = start_idx

            while end_idx < len(spike_times) and spike_times[end_idx] <= trigger_t + duration:
                end_idx += 1

            if start_idx < len(spike_times):
                if subtract_trigger:
                    yield spike_times[start_idx:end_idx] - trigger_t
                else:
                    yield spike_times[start_idx:end_idx]

    def get_jpsth(self, spike_times_a, spike_times_b, trigger_times, duration=None, bin_size=2e-3):

        if duration is None and len(trigger_times) > 1:
            duration = trigger_times[1] - trigger_times[0]

        jpsth = np.zeros((bin_size, bin_size), dtype=int)

        for spikes_a, spikes_b in zip(self.get_triggered_spikes(spike_times=spike_times_a,
                                                                trigger_times=trigger_times,
                                                                duration=duration),
                                      self.get_triggered_spikes(spike_times=spike_times_b,
                                                                trigger_times=trigger_times,
                                                                duration=duration)):

            for idx_a in np.floor(spikes_a / bin_size):
                for idx_b in np.floor(spikes_b / bin_size):
                    jpsth[idx_a, idx_b] += 1

        return jpsth

    def plot_jpsth(self, spike_times_a, spike_times_b, trigger_times, duration=None, bin_size=2e-3, title=None):

        if duration is None and len(trigger_times) > 1:
            duration = trigger_times[1] - trigger_times[0]

        spike_combos = []

        for spikes_a, spikes_b in zip(self.get_triggered_spikes(spike_times=spike_times_a,
                                                                trigger_times=trigger_times,
                                                                duration=duration),
                                      self.get_triggered_spikes(spike_times=spike_times_b,
                                                                trigger_times=trigger_times,
                                                                duration=duration)):
            for t_a in spikes_a:
                for t_b in spikes_b:
                    spike_combos.append((t_a, t_b))

        x, y = zip(*spike_combos)

        plt.figure()
        plt.hist2d(x=np.array(x), y=np.array(y),
                   bins=int(np.ceil(duration/bin_size)+1), range=[[0, duration], [0, duration]],
                   cmap=plt.get_cmap("Reds"))
        plt.plot([0, duration], [0, duration], color="k", linestyle="dotted")
        plt.colorbar()
        plt.title(title)
        plt.xlabel("Ablated: time (s)")
        plt.ylabel("Fully ablated: time (s)")
        plt.ion()
        plt.show()

    def plot_jpsth_all(self, data_key_a, data_key_b, trigger_times, duration=None, bin_size=2e-3, time_range=None):

        if duration is None and len(trigger_times) > 1:
            duration = trigger_times[1] - trigger_times[0]

        assert self.check_same_neurons(data_key_a, data_key_b), f"data_keys have different neurons in the network"

        neuron_names = self.simulation_data[data_key_a].get_neuron_name()

        # Match spikes against each other, compute change...
        sim_a = self.simulation_data[data_key_a]
        sim_b = self.simulation_data[data_key_b]

        spikes_a = sim_a.get_spikes()
        spikes_b = sim_b.get_spikes()

        for neuron_id in spikes_a.keys():
            s_a = self.filter_times(spikes_a[neuron_id].flatten(), time_range)
            s_b = self.filter_times(spikes_b[neuron_id].flatten(), time_range)

            if s_a.size > 0 and s_b.size > 0:
                self.plot_jpsth(spike_times_a=s_a, spike_times_b=s_b,
                                trigger_times=trigger_times, duration=duration, bin_size=bin_size,
                                title=f"{neuron_names[neuron_id]} ({neuron_id})")

    def get_spike_deltas(self, data_key_a, data_key_b, matching_method, delta_t=5e-3, time_range=None):

        """
        Args:
            data_key_a : Data key for first dataset
            data_key_b : Data key for second dataset
            matching_method : Method to use for spike matching "closest", "order", "closestunique", "spike_triggered"
            delta_t : Optional, used by "closestunique" method (default 5e-3s)

        """

        # Check that the neurons compared are the same (by verifying parameter key, morphology key, modulation key)
        assert self.check_same_neurons(data_key_a, data_key_b), f"data_keys have different neurons in the network"

        # Match spikes against each other, compute change...

        sim_a = self.simulation_data[data_key_a]
        sim_b = self.simulation_data[data_key_b]

        spikes_a = sim_a.get_spikes()
        spikes_b = sim_b.get_spikes()

        spike_time_difference = dict()

        for neuron_id in spikes_a.keys():
            s_a = self.filter_times(spikes_a[neuron_id].flatten(), time_range)
            s_b = self.filter_times(spikes_b[neuron_id].flatten(), time_range)

            if s_a.size > 0 and s_b.size > 0:

                if matching_method == "closest":
                    spike_time_difference[neuron_id] = self.match_closest_spikes(s_a, s_b)
                elif matching_method == "order":
                    n_compare = min(s_a.size, s_b.size)
                    spike_time_difference[neuron_id] = self.match_order_spikes(s_a, s_b)
                elif matching_method == "closestunique":
                    spike_time_difference[neuron_id] = self.match_closest_unique(s_a, s_b, delta_t=delta_t)
                elif matching_method == "spike_triggered":
                    spike_time_difference[neuron_id] = self.get_spike_triggered_deltas(s_a, s_b), np.array([])
                else:
                    assert False, f"Unknown matching_method={matching_method}, " \
                                  f"please use ('closest', 'order', 'closestunique', 'spike_triggered')"
            else:
                # At least one of the spike trains does not have any spikes
                spike_time_difference[neuron_id] = np.array([]), np.array([])

        return spike_time_difference

    def filter_times(self, spike_times, time_range=None):

        if time_range is None:
            return spike_times
        else:
            return spike_times[time_range[0] <= spike_times <= time_range[1]]

    def plot_spike_delta_histogram(self, data_key_a=None, data_key_b=None,
                                   plot_title=None, direction=0,
                                   matching_method=None,
                                   range_min=-10e-3, range_max=2e-3, bin_size=0.5e-3):

        spike_time_difference = self.get_spike_deltas(data_key_a=data_key_a, data_key_b=data_key_b,
                                                      matching_method=matching_method)
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