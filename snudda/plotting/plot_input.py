import numpy as np
import h5py
import os
from collections import OrderedDict
import matplotlib.pyplot as plt
from matplotlib import cm

from snudda.utils.load import SnuddaLoad


class PlotInput(object):

    def __init__(self, input_file, network_path=None):

        self.input_data = None

        if input_file:
            self.load_input(input_file)

        if not network_path:
            network_path = os.path.dirname(input_file)

        network_file = os.path.join(network_path, "network-synapses.hdf5")

        if os.path.exists(network_file):
            self.network_info = SnuddaLoad(network_file)
        else:
            print(f"Specify a network_path with a network file, to get neuron type in figure")
            self.network_info = None

    def load_input(self, input_file):
        self.input_data = h5py.File(input_file, "r")

    def extract_input(self, input_target):

        data = OrderedDict()

        if input_target in self.input_data["input"]:
            for input_type in self.input_data["input"][input_target]:
                input_info = self.input_data["input"][input_target][input_type]

                data[input_type] = input_info["spikes"][()]

        return data
    
    def get_neuron_name(self, neuron_id):

        neuron_id = int(neuron_id)

        if self.network_info:
            neuron_name = self.network_info.data["neurons"][neuron_id]["name"]
        else:
            neuron_name = ""
            
        return neuron_name

    def plot_input(self, neuron_type, num_neurons, fig_size=None):

        neuron_id = self.network_info.get_neuron_id_of_type(neuron_type=neuron_type,
                                                            num_neurons=num_neurons,
                                                            random_permute=True)
        target_id = [str(x) for x in np.sort(neuron_id)]

        if len(target_id) == 0:
            print(f"No neurons of type {neuron_type}")
            return

        self.plot_input_to_target(target_id, fig_size=fig_size)

    def plot_input_population_unit(self, population_unit_id, num_neurons, neuron_type=None, fig_size=None):

        if not population_unit_id:
            population_unit_id = 0  # 0 = no population

        assert type(population_unit_id) == int

        neuron_id = self.network_info.get_population_unit_members(population_unit_id)

        assert np.array([self.network_info.data["populationUnit"][x] == population_unit_id for x in neuron_id]).all()

        if neuron_type:
            neuron_id2 = self.network_info.get_neuron_id_of_type(neuron_type)
            neuron_id = list(set(neuron_id).intersection(set(neuron_id2)))

        if num_neurons:
            num_neurons = min(num_neurons, len(neuron_id))
            neuron_id = np.random.permutation(neuron_id)[:num_neurons]

        target_id = [str(x) for x in np.sort(neuron_id)]

        if len(target_id) == 0:
            print(f"No neurons with population id {population_unit_id}")
            return

        assert np.array([self.network_info.data["populationUnit"][int(x)] == population_unit_id
                         for x in target_id]).all()

        self.plot_input_to_target(target_id, fig_size=fig_size)

    def plot_input_to_target(self, input_target, fig_size=None):

        if not fig_size:
            fig_size = (10, 5)

        if type(input_target) != list:
            input_target = [input_target]

        # Make sure each target is a str
        input_target = [str(x) for x in input_target]
        colours = cm.get_cmap('tab20', len(input_target) * 2)

        y_pos = 0
        input_ctr = 0
        plt.figure(figsize=fig_size)

        ytick_pos = []
        ytick_label = []

        for it in input_target:

            data = self.extract_input(it)

            for input_type in data:

                y_pos_start = y_pos
                for spike_train in data[input_type]:
                    idx = np.where(spike_train > 0)[0]
                    plt.scatter(spike_train[idx], y_pos * np.ones((len(idx),)),
                                color=colours(input_ctr), marker='.', s=7)
                    y_pos += 1

                y_pos_avg = (y_pos + y_pos_start)/2
                ytick_pos.append(y_pos_avg)
                ytick_label.append(f"{input_type}→{self.get_neuron_name(it)} ({it})")

                input_ctr += 1
                y_pos += 5

            y_pos += 5

        # Add yticks
        ax = plt.gca()
        ax.invert_yaxis()
        ax.set_yticks(ytick_pos)
        ax.set_yticklabels(ytick_label)
        ax.set_xlabel("Time (s)")
        plt.ion()
        plt.show()

