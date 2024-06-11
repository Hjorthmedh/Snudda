import os
import numpy as np
import h5py
from snudda.utils import SnuddaLoad


class AnalyseInput:

    def __init__(self, network_path=None, input_file=None, network_file=None):

        if network_path is None:
            if network_file is not None:
                network_path = os.path.dirname(network_file)
            elif input_file is not None:
                network_path = os.path.dirname(input_file)
        self.network_path = network_path

        if input_file is None and network_path is not None:
            input_file = os.path.join(self.network_path, "network-spikes.hdf5")
        self.input_file = input_file

        if network_file is None and network_path is not None:
            network_file = os.path.join(network_path, "network-synapses.hdf5")
        self.network_file = network_file

        self.snudda_load = SnuddaLoad(network_file=network_file)
        self.input_data = h5py.File(input_file, "r")

        self.figure_path = os.path.join(network_path, "figures")

    def count_inputs_helper(self):

        input_num = dict()

        for neuron_id_str in self.input_data["input"].keys():
            neuron_id = int(neuron_id_str)
            input_num[neuron_id] = dict()

            for input_name in self.input_data["input"][neuron_id_str].keys():
                input_num[neuron_id][input_name] = \
                    self.input_data["input"][neuron_id_str][input_name]["spikes"].shape[0]

        return input_num

    def count_inputs(self):

        input_num = self.count_inputs_helper()
        input_summary = dict()

        neuron_types = self.snudda_load.get_neuron_types(return_set=True)

        for neuron_type in neuron_types:
            input_summary[neuron_type] = dict()

            neuron_id_list = self.snudda_load.get_neuron_id_of_type(neuron_type=neuron_type)

            for neuron_id in neuron_id_list:
                for input_type in input_num[neuron_id].keys():
                    if input_type not in input_summary[neuron_type]:
                        input_summary[neuron_type][input_type] = []

                    input_summary[neuron_type][input_type].append(input_num[neuron_id][input_type])

        return input_summary

    def plot_input_count(self, fig_name=None):

        if not os.path.exists(self.figure_path):
            os.mkdir(self.figure_path)

        import matplotlib.pyplot as plt

        input_summary = self.count_inputs()

        input_list = []
        label_list = []

        for neuron_type in sorted(list(input_summary.keys())):
            for input_type in input_summary[neuron_type].keys():
                if len(input_summary[neuron_type][input_type]) > 0:
                    input_list.append(input_summary[neuron_type][input_type])
                    label_list.append(f"{neuron_type}\n{input_type}")

        plt.figure()
        ax = plt.subplot()
        ax.violinplot(dataset=input_list)
        ax.set_xticks(np.arange(1, len(label_list)+1))
        ax.set_xticklabels(label_list)
        ax.set_ylabel("Number of input synapses")

        if fig_name is not None:
            fig_path = os.path.join(self.figure_path, fig_name)
            plt.savefig(fig_path)
            
        plt.ion()
        plt.show()
        
