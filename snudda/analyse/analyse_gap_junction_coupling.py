import os
import numpy as np
import matplotlib.pyplot as plt

from snudda.simulate.pair_recording import PairRecording
from snudda.utils.load import SnuddaLoad
from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation


class AnalyseGapJunctionCoupling:

    def __init__(self, network_path, experiment_config_file):

        self.network_path = network_path
        self.figure_path = os.path.join(network_path, "figures")
        self.experiment_file = experiment_config_file
        self.experiment_config = PairRecording.read_experiment_config(experiment_config_file)
        self.network_simulation = None
        self.network_info = SnuddaLoad(network_file=os.path.join(network_path, "network-synapses.hdf5"))
        self.neighbour_id_cache = dict()

        self.load_simulation_data()

        if not os.path.exists(self.figure_path):
            os.makedirs(self.figure_path)

    def load_simulation_data(self):

        if "meta" in self.experiment_config and "pairRecordingOutputFile" in self.experiment_config["meta"]:
            output_file = self.experiment_config["meta"]["pairRecordingOutputFile"]
        else:
            output_file = None

        self.network_simulation = \
            SnuddaLoadNetworkSimulation(network_simulation_output_file=os.path.join(self.network_path,
                                                                                    "simulation",
                                                                                    output_file))

    def get_gap_junction_connection_matrix(self):

        num_neurons = len(self.network_info.data["neuronID"])
        gj_connection_matrix = np.zeros((num_neurons, num_neurons), dtype=int)

        for gj_list in self.network_info.gap_junction_iterator():
            for gj in gj_list:
                gj_connection_matrix[gj[0], gj[1]] += 1
                gj_connection_matrix[gj[1], gj[0]] += 1

        return gj_connection_matrix

    def get_gj_list(self, n_pairs=10):

        gj_matrix = self.get_gap_junction_connection_matrix()
        # We sort the items in the upper triangular part of the matrix
        gj_list = np.triu(gj_matrix).flatten()
        idx = np.argsort(-gj_list)

        # Get back original 2-d indexing
        row, col = np.unravel_index(idx, gj_matrix.shape)

        n_max = np.sum(gj_list[idx] > 0)

        if n_pairs is not None:
            n_shown = min(n_pairs, n_max)
        else:
            n_shown = n_max

        return row[:n_shown], col[:n_shown], gj_list[idx[:n_shown]]

    def print_gj_list(self, n_pairs=10):

        gj_a_list, gj_b_list, num_gj_list = self.get_gj_list(n_pairs=n_pairs)

        for gj_a, gj_b, num_gj in zip(gj_a_list, gj_b_list, num_gj_list):
            print(f"{gj_a} -- {gj_b} : {num_gj} gap junctions")

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

                if ((duration is None or np.abs(dur - duration) < 1e-5)
                    and (amplitude is None or np.abs(amp - amplitude) < 1e-13)):
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
        trace_info = dict()

        for neuron_id, start_time, end_time in self.get_intervals(duration=duration, amplitude=amplitude):

            info = dict()
            info["pre_neuron_id"] = neuron_id
            info["start_time"] = start_time
            info["end_time"] = end_time

            info["post_neuron_id"] = []
            info["coupling"] = []

            pre_amp = self.extract_amplitude(neuron_id, start_time, end_time)

            coupled_neuron_id = self.find_gap_junction_neighbours(neuron_id=neuron_id)
            for cid in coupled_neuron_id:
                post_amp = self.extract_amplitude(cid, start_time, end_time)
                coupling = post_amp / pre_amp
                coupling_factors.append(coupling)
                info["post_neuron_id"].append(cid)
                info["coupling"].append(coupling)

            if neuron_id not in trace_info:
                trace_info[neuron_id] = []

            trace_info[neuron_id].append(info)

        return coupling_factors, trace_info

    def get_trace(self, neuron_id_list, start_time, end_time, time_padding=0):

        time = self.network_simulation.get_time()
        idx = np.where(np.logical_and(start_time - time_padding <= time, time <= end_time + time_padding))[0]

        traces = np.zeros((len(idx), len(neuron_id_list)))

        for i, neuron_id in enumerate(neuron_id_list):
            volt = self.network_simulation.get_voltage(neuron_id=neuron_id)
            traces[:, i] = volt[idx][:, 0]  # 0 column is usually soma

        return traces, time[idx]

    def plot_coupling(self, duration, amplitude):

        coupling_factors, trace_info = self.extract_coupling(duration=duration, amplitude=amplitude)

        for neuron_id, trace_data_list in trace_info.items():

            fig = plt.figure()
            coupling = []

            for trace_data in trace_data_list:
                assert neuron_id == trace_data["pre_neuron_id"]
                post_neuron_id = trace_data["post_neuron_id"]
                coupling = coupling + trace_data["coupling"]

                get_id = [neuron_id, *post_neuron_id]
                volt, time = self.get_trace(neuron_id_list=get_id,
                                            start_time=trace_data["start_time"],
                                            end_time=trace_data["end_time"],
                                            time_padding=0.05)
                plt.plot(time, volt[:, 0], 'k', linewidth=2)
                plt.plot(time, volt[:, 1:], 'k', linewidth=1)

            plt.title(f"Coupling: {np.mean(coupling):.4f} ({np.min(coupling):.4f} - {np.max(coupling):.4f}), "
                      f"duration {duration*1e3:.0f} ms (neuron {neuron_id})")

            fig_name = os.path.join(self.figure_path,
                                    f"gap-junction-coupling-neuron-{neuron_id}-duration-{duration}.png")

            plt.savefig(fig_name, dpi=300)

        plt.show()
        plt.pause(1)
