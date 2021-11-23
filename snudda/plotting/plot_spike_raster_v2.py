import os
import numpy as np

from snudda.utils.load import SnuddaLoad

import matplotlib.pyplot as plt

from snudda.utils.load_network_simulation import SnuddaLoadNetworkSimulation


class SnuddaPlotSpikeRaster2:

    def __init__(self, network_path, network_file=None, simulation_file=None, figure_path=None):

        self.network_path = network_path

        if network_file:
            self.network_file = network_file
        else:
            self.network_file = os.path.join(self.network_path, "network-synapses.hdf5")

        if simulation_file:
            self.simulation_file = simulation_file
        else:
            self.simulation_file = os.path.join(self.network_path, "simulation", "network-output.hdf5")

        if figure_path:
            self.figure_path = figure_path
        else:
            self.figure_path = os.path.join(self.network_path, "figures", "network-spike-raster.png")

        self.snudda_load = SnuddaLoad(network_file=self.network_file)

        self.snudda_simulation_load = SnuddaLoadNetworkSimulation(network_simulation_output_file=self.simulation_file)
        spike_data = self.snudda_simulation_load.merge_spikes()

        self.spike_time = spike_data[:, 0]
        self.spike_neuron_id = spike_data[:, 1].astype(int)

    def make_figures_directory(self):

        fig_dir = os.path.join(self.network_path, "figures")

        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)

    @staticmethod
    def get_colours(neuron_type):

        # TODO: Read colours from a JSON file

        colours = {"dSPN".lower(): (77. / 255, 151. / 255, 1.0),
                   "iSPN".lower(): (67. / 255, 55. / 255, 181. / 255),
                   "FS".lower(): (6. / 255, 31. / 255, 85. / 255),
                   "FSN".lower(): (6. / 255, 31. / 255, 85. / 255),
                   "ChIN".lower(): (252. / 255, 102. / 255, 0.0),
                   "LTS".lower(): (150. / 255, 63. / 255, 212. / 255)}

        return colours[neuron_type.lower()]

    def get_all_colours(self):

        neuron_type_list = self.snudda_load.get_neuron_types(return_set=False)
        neuron_colours = np.zeros((len(neuron_type_list), 3))

        for idx, nt in enumerate(neuron_type_list):
            neuron_colours[idx, :] = self.get_colours(nt)

        return neuron_colours

    def plot_spike_raster(self, type_order=None, skip_time=0, end_time=None, fig_size=None):

        self.make_figures_directory()

        fig = plt.figure(figsize=fig_size)
        ax = fig.add_subplot()

        # Gets a list of all the neurons' types
        neuron_type_list = self.snudda_load.get_neuron_types(return_set=False)
        neuron_type_map = dict()

        if type_order is None:
            unique_neuron_types = set(neuron_type_list)
        else:
            unique_neuron_types = type_order + list(set(neuron_type_list) - set(type_order))

        for nt_idx, nt in enumerate(unique_neuron_types):
            neuron_type_map[nt] = nt_idx

        # For each neuron, associate the number of the type it is
        neuron_type_idx = np.array([neuron_type_map[x] for x in neuron_type_list])
        neuron_order = np.argsort(neuron_type_idx)
        neuron_order_lookup = np.zeros(neuron_order.shape)

        for idx, no in enumerate(neuron_order):
            neuron_order_lookup[no] = idx

        spike_y = np.take(neuron_order_lookup, self.spike_neuron_id)

        colour_lookup = self.get_all_colours()
        sc = np.zeros((len(spike_y), 3))

        for i in range(0, 3):
            sc[:, i] = np.take(colour_lookup[:, i], self.spike_neuron_id)

        ax.scatter(self.spike_time-skip_time, spike_y, color=sc, s=5, linewidths=0.1)

        # Get position of labels
        unique_neuron_types = set(neuron_type_list)
        y_tick = []
        y_tick_label = []
        for nt in unique_neuron_types:
            y_tick_label.append(nt)
            y_tick.append(np.mean(neuron_order_lookup[np.where([x == nt for x in neuron_type_list])[0]]))

        ax.invert_yaxis()
        ax.set_xlabel('Time (s)')
        ax.set_yticks(y_tick)
        ax.set_yticklabels(y_tick_label)

        if skip_time or end_time:
            x_lim = ax.get_xlim()
            x_lim[0] = 0
            if end_time:
                x_lim[1] = end_time
            ax.set_xlim(x_lim)

        if not os.path.isdir(os.path.basename(self.figure_path)):
            os.makedirs(os.path.basename(self.figure_path))

        plt.savefig(self.figure_path, dpi=300)

        plt.ion()
        plt.show()


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser("Scatter plot")
    parser.add_argument("network_path", type=str, help="Network path")

    args = parser.parse_args()
    ps = SnuddaPlotSpikeRaster2(network_path=args.network_path)

    type_order = ["FS", "dSPN", "LTS", "iSPN", "ChIN"]
    ps.plot_spike_raster(type_order)
