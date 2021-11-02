import os
import numpy as np
from snudda.utils.load import SnuddaLoad
import matplotlib.pyplot as plt


class PlotDensity(object):

    def __init__(self, network):

        if os.path.isdir(network):
            network_file = os.path.join(network, "network-synapses.hdf5")
        else:
            network_file = network

        self.network_file = network_file
        self.network_path = os.path.dirname(self.network_file)

        self.sl = SnuddaLoad(self.network_file)

    def close(self):
        self.sl.close()

    def plot(self, neuron_type, plot_axis, n_bins=5):

        p_axis = {"x": 0, "y": 1, "z": 2}
        assert plot_axis in p_axis, f"plot_axis must be one of {', '.join(p_axis.keys())}"
        neuron_pos = self.sl.data["neuronPositions"]
        cell_id = self.sl.get_neuron_id_of_type(neuron_type)

        fig = plt.figure()
        plt.hist(neuron_pos[cell_id, p_axis[plot_axis]], bins=n_bins)
        plt.title(f"Density of {neuron_type}")
        plt.xlabel(plot_axis)
        plt.ylabel("Count")

        fig_path = os.path.join(self.network_path, "figures")
        if not os.path.exists(fig_path):
            os.mkdir(fig_path)

        plt.savefig(os.path.join(fig_path, f"density-hist-{plot_axis}.png"))
        plt.ion()
        plt.show()