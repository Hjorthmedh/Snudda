import os
import numpy as np
from numba import jit
import matplotlib.pyplot as plt

import h5py
# from elephant.spike_train_correlation import spike_time_tiling_coefficient
from snudda.analyse.spike_time_tiling_coefficient import spike_time_tiling_coefficient
# from neo import SpikeTrain as NeoSpikeTrain
# import quantities as pq

from snudda.utils.load_network_simulation import SnuddaLoadSimulation
from snudda.utils.load import SnuddaLoad


class AnalyseSpikeTrains:

    def __init__(self, network_file=None, output_file=None, input_file=None, network_path=None):

        if network_path is None:
            if network_file is not None:
                network_path = os.path.dirname(network_file)
            elif output_file is not None:
                if "simulation" in output_file:
                    network_path = os.path.dirname(os.path.dirname(output_file))
            elif input_file is not None:
                network_path = os.path.dirname(input_file)

        if network_path is not None:
            if network_file is None:
                putative_network_file = os.path.join(network_path, "network-synapses.hdf5")
                if os.path.exists(putative_network_file):
                    network_file = putative_network_file

            if output_file is None:
                putative_output_file = os.path.join(network_path, "simulation", "output.hdf5")
                if os.path.exists(putative_output_file):
                    output_file = putative_output_file

            if input_file is None:
                putative_input_file = os.path.join(network_path, "input-spikes.hdf5")
                if os.path.exists(putative_input_file):
                    input_file = putative_input_file

        self.network_path = network_path
        self.network_file = network_file
        self.output_file = output_file
        self.input_file = input_file

        if output_file is not None:
            self.output_data = SnuddaLoadSimulation(network_simulation_output_file=self.output_file)
        else:
            self.output_data = None

        if self.network_file is not None:
            self.network_data = SnuddaLoad(network_file=self.network_file)
        else:
            self.network_data = None

        if self.input_file is not None:
            self.input_data = h5py.File(self.input_file, "r")

    def calculate_sttc(self, spike_train_a, spike_train_b, dt, start_time=0, end_time=None):

        if end_time is None:
            end_time = self.get_end_time()

        idx_a = np.where(np.logical_and(start_time <= spike_train_a, spike_train_a <= end_time))[0]
        idx_b = np.where(np.logical_and(start_time <= spike_train_b, spike_train_b <= end_time))[0]

        return spike_time_tiling_coefficient(spiketrain_i=spike_train_a[idx_a], spiketrain_j=spike_train_b[idx_b],
                                             dt=dt, end_time=end_time, start_time=start_time)

    def get_end_time(self):
        return np.max(self.output_data.get_time())

    def calculate_sttc_all_to_all(self, spike_trains, n_spikes, dt, start_time=0, end_time=None):
        n_spike_trains = len(n_spikes)
        assert spike_trains.shape[0] == n_spike_trains

        corr = []

        if end_time is None:
            end_time = self.get_end_time()

        pruned_spike_trains = []
        for st in spike_trains:
            idx = np.where(np.logical_and(start_time <= st, st <= end_time))[0]
            pruned_spike_trains.append(st[idx].flatten())

        for i in range(0, n_spike_trains):
            if i % 50 == 0:
                print(f'{i} / {n_spike_trains}')
            for j in range(i+1, n_spike_trains):
                corr.append(spike_time_tiling_coefficient(spiketrain_a=pruned_spike_trains[i],
                                                          spiketrain_b=pruned_spike_trains[j],
                                                          dt=dt, start_time=start_time, end_time=end_time))
        print(f'{i} / {n_spike_trains}')

        return np.array(corr)

    def calculate_sttc_one_to_all(self, spike_train, spike_trains, n_spikes, dt, start_time=0, end_time=None):
        n_spike_trains = len(n_spikes)
        assert spike_trains.shape[0] == n_spike_trains
        corr = []

        if end_time is None:
            end_time = self.get_end_time()

        st = spike_train.flatten()
        idx = np.where(np.logical_and(start_time <= st, st <= end_time))[0]
        spike_train_b = st[idx]

        for i in range(1, n_spike_trains):
            idx_i = np.where(np.logical_and(start_time <= spike_trains[i, :], spike_trains[i, :] <= end_time))[0]
            spike_train_i = spike_trains[i, idx_i].flatten()

            corr.append(spike_time_tiling_coefficient(spiketrain_a=spike_train_i,
                                                      spiketrain_b=spike_train_b,
                                                      dt=dt, end_time=end_time))

        return np.array(corr)

    def input_correlation(self, neuron_id, input_type, dt):
        input_spikes = self.input_data[f"input/{neuron_id}/{input_type}/spikes"][()].copy()
        n_spikes = self.input_data[f"input/{neuron_id}/{input_type}/spikes"].attrs["num_spikes"].copy()

        corr = self.calculate_sttc_all_to_all(spike_trains=input_spikes, n_spikes=n_spikes, dt=dt)
        return corr

    def input_output_correlation(self, neuron_id, dt, start_time=0, end_time=None):

        input_data = dict()
        corr = dict()

        # First gather all input spikes
        for input_type, input_spikes in self.input_data[f"input/{neuron_id}"].items():
            input_data[input_type] = input_spikes["spikes"][()].copy(), input_spikes["spikes"].attrs["num_spikes"].copy()

        output_data = self.output_data.get_spikes(neuron_id=neuron_id)

        for input_type, (input_spikes, n_spikes) in input_data.items():
            corr[input_type] = self.calculate_sttc_one_to_all(spike_train=output_data,
                                                              spike_trains=input_spikes,
                                                              n_spikes=n_spikes,
                                                              dt=dt, start_time=start_time, end_time=end_time)
        return corr

    def plot_spike_multiplicity(self, neuron_id, input_type, jitter=0, start_time=None, end_time=None):

        mult_list, mult = self.calculate_spike_multiplicity(neuron_id=neuron_id, input_type=input_type, jitter=jitter,
                                                            start_time=None, end_time=None)

        plt.figure()
        plt.title(f"{self.network_path} {neuron_id} {input_type}")
        # plt.stairs(edges=np.arange(0, len(mult)+1), values=mult)
        plt.hist(mult_list, bins=50)
        plt.yscale('log')

    def calculate_spike_multiplicity(self, neuron_id, input_type, jitter=0, start_time=None, end_time=None):

        input_spikes = self.input_data[f"input/{neuron_id}/{input_type}/spikes"][()].flatten()
        return self.calculate_multiplicity_helper(input_spikes=input_spikes, jitter=jitter,
                                                  start_time=start_time, end_time=end_time)

    @staticmethod
    # @jit(nopython=True, fastmath=True, cache=True)
    def calculate_multiplicity_helper(input_spikes, jitter=0, start_time=None, end_time=None):

        assert (start_time is None) ^ (end_time is None) == False, "Either both start_time and end_time is set, or None"

        max_mult = input_spikes.shape[0]
        multiplicity = np.zeros((max_mult+1, ), dtype=int)

        input_spikes = input_spikes.flatten()
        if end_time is None:
            input_spikes = np.sort(input_spikes[np.where(input_spikes >= 0)[0]])
        else:
            input_spikes = np.sort(input_spikes[np.where(np.logical_and(start_time <= input_spikes,
                                                                        input_spikes <= end_time))[0]])

        mul_ctr = 1
        first_spike = input_spikes[0]
        max_mult = 1

        mult_list = []

        for idx in range(1, len(input_spikes)):
            if input_spikes[idx] <= first_spike + jitter:
                mul_ctr += 1
            else:
                mult_list.append(mul_ctr)
                multiplicity[mul_ctr] += 1
                max_mult = max(max_mult, mul_ctr)
                mul_ctr = 1
                first_spike = input_spikes[idx]

        # Need to add the last spike also
        mult_list.append(mul_ctr)
        multiplicity[mul_ctr] += 1
        max_mult = max(max_mult, mul_ctr)

        return mult_list, multiplicity[:max_mult+1]


