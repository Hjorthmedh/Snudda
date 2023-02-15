import os
import numpy as np

import h5py
from elephant.spike_train_correlation import spike_time_tiling_coefficient
from neo import SpikeTrain as NeoSpikeTrain
import quantities as pq

from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation
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
            self.output_data = SnuddaLoadNetworkSimulation(network_simulation_output_file=self.output_file)
        else:
            self.output_data = None

        if self.network_file is not None:
            self.network_data = SnuddaLoad(network_file=self.network_file)
        else:
            self.network_data = None

        if self.input_file is not None:
            self.input_data = h5py.File(self.input_file, "r")

    def calculate_sttc(self, spike_train_a, spike_train_b, dt):
        t_end = self.get_end_time() * pq.s
        sa = NeoSpikeTrain(spike_train_a.flatten(), t_stop=t_end, units="s")
        sb = NeoSpikeTrain(spike_train_b.flatten(), t_stop=t_end, units="s")

        return spike_time_tiling_coefficient(spiketrain_i=sa, spiketrain_j=sb, dt=dt)

    def get_end_time(self):
        return np.max(self.output_data.get_time()) * pq.s

    def calculate_sttc_all_to_all(self, spike_trains, n_spikes, dt):
        n_spike_trains = len(n_spikes)
        assert spike_trains.shape[0] == n_spike_trains

        corr = []
        t_end = self.get_end_time()

        neo_spike_trains = []
        for i in range(0, n_spike_trains):
            neo_spike_trains.append(NeoSpikeTrain(spike_trains[i, :n_spikes[i]].flatten(), t_stop=t_end, units="s"))

        for i in range(0, n_spike_trains):
            if i % 50 == 0:
                print(f'{i} / {n_spike_trains}')
            for j in range(i+1, n_spike_trains):
                corr.append(spike_time_tiling_coefficient(spiketrain_i=neo_spike_trains[i],
                                                          spiketrain_j=neo_spike_trains[j],
                                                          dt=dt))
        print(f'{i} / {n_spike_trains}')

        return np.array(corr)

    def calculate_sttc_one_to_all(self, spike_train, spike_trains, n_spikes, dt):
        n_spike_trains = len(n_spikes)
        assert spike_trains.shape[0] == n_spike_trains
        corr = []
        t_end = self.get_end_time()

        spike_train_b = NeoSpikeTrain(spike_train.flatten(), t_stop=t_end, units="s")

        for i in range(1, n_spike_trains):
            spike_train_a = NeoSpikeTrain(spike_trains[i, :n_spikes[i]].flatten(), t_stop=t_end, units="s")
            corr.append(spike_time_tiling_coefficient(spiketrain_i=spike_train_a,
                                                      spiketrain_j=spike_train_b,
                                                      dt=dt))

        return np.array(corr)

    def input_correlation(self, neuron_id, input_type, dt):
        input_spikes = self.input_data[f"input/{neuron_id}/{input_type}/spikes"][()].copy()
        n_spikes = self.input_data[f"input/{neuron_id}/{input_type}/spikes"].attrs["nSpikes"].copy()

        corr = self.calculate_sttc_all_to_all(spike_trains=input_spikes, n_spikes=n_spikes, dt=dt)
        return corr

    def input_output_correlation(self, neuron_id, dt):

        input_data = dict()
        corr = dict()

        # First gather all input spikes
        for input_type, input_spikes in self.input_data[f"input/{neuron_id}"].items():
            input_data[input_type] = input_spikes["spikes"][()].copy(), input_spikes["spikes"].attrs["nSpikes"].copy()

        output_data = self.output_data.get_spikes(neuron_id=neuron_id)

        for input_type, (input_spikes, n_spikes) in input_data.items():
            corr[input_type] = self.calculate_sttc_one_to_all(spike_train=output_data,
                                                              spike_trains=input_spikes,
                                                              n_spikes=n_spikes,
                                                              dt=dt)
        return corr

