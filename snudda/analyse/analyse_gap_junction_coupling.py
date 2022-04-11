import os
import numpy as np
import matplotlib.pyplot as plt

from snudda.simulate.pair_recording import PairRecording
from snudda.utils.load import SnuddaLoad
from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation


class AnalyseGapJunctionCoupling:

    def __init__(self, network_path, experiment_config_file):

        self.network_path = network_path
        self.experiment_file = experiment_config_file
        self.experiment_config = PairRecording.read_experiment_config(experiment_config_file)
        self.network_simulation = None
        self.network_info = SnuddaLoad(network_file=os.path.join(network_path, "network-synapses.hdf5"))
        self.neighbour_id_cache = dict()

        self.load_simulation_data()

    def load_simulation_data(self):

        if "meta" in self.experiment_config and "pairRecordingOutputFile" in self.experiment_config["meta"]:
            output_file = self.experiment_config["meta"]["pairRecordingOutputFile"]
        else:
            output_file = None

        self.network_simulation = \
            SnuddaLoadNetworkSimulation(network_simulation_output_file=os.path.join(self.network_path,
                                                                                    "simulation",
                                                                                    output_file))

    def find_gap_junction_neighbours(self, neuron_id):

        if neuron_id not in self.neighbour_id_cache:
            gj_matrix, gj_coords = self.network_info.find_gap_junctions(neuron_id=neuron_id)

            # First and second columns are id of neurons coupled by gap junctions, neuron_id is in one or the other
            # we need to find the id of the coupled neuron
            coupled_neuron_id = set(gj_matrix[:, 0]).union(gj_matrix[:, 1])
            coupled_neuron_id.remove(neuron_id)
            self.neighbour_id_cache[neuron_id] = coupled_neuron_id

        return self.neighbour_id_cache[neuron_id]

    def get_intervals(self, amplitude, duration):

        for cur_inj in self.experiment_config["currentInjection"]:
            neuron_id = cur_inj["neuronID"]
            cur_start_time = cur_inj["start"]
            cur_end_time = cur_inj["end"]
            cur_amplitude = cur_inj["amplitude"]

            if type(cur_start_time) != list:
                cur_start_time = [cur_start_time]

            if type(cur_end_time) != list:
                cur_end_time = [cur_end_time]

            if type(cur_amplitude) != list:
                cur_amplitude = [cur_amplitude] * len(cur_start_time)

            for st, et, amp in zip(cur_start_time, cur_end_time, cur_amplitude):
                dur = et - st
                if (duration is None or dur == duration) and (amplitude is None or amp == amplitude):
                    yield neuron_id, st, et

    def extract_amplitude(self, neuron_id, start_time, end_time):

        baseline_time = [0.1, 0.2]

        time = self.network_simulation.get_time()
        idx = np.where(np.logical_and(start_time <= time, time <= end_time))[0]

        volt = self.network_simulation.get_voltage(neuron_id=neuron_id)
        baseline_volt = np.mean(volt[np.where(np.logical_and(baseline_time[0] <= time, time <= baseline_time[1]))[0]])

        peak_volt = np.max(np.abs(volt[idx] - baseline_volt))

        return peak_volt

    def extract_coupling(self, duration=None, amplitude=None):

        coupling_factors = []

        for neuron_id, start_time, end_time in self.get_intervals(duration=duration, amplitude=amplitude):

            coupled_neuron_id = self.find_gap_junction_neighbours(neuron_id=neuron_id)
            for cid in coupled_neuron_id:
                coupling_factors.append(self.extract_amplitude(cid, start_time, end_time))

        return coupling_factors
