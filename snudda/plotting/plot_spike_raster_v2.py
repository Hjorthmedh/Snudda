import os
import numpy as np

from snudda import SnuddaLoad

import matplotlib.pyplot as plt


class SnuddaPlotSpikeRaster2:

    def __init__(self, network_path, network_file=None, spike_file=None, figure_path=None):

        self.network_path = network_path

        if network_file:
            self.network_file = network_file
        else:
            self.network_file = os.path.join(self.network_path, "network-synapses.hdf5")

        if spike_file:
            self.spike_file = spike_file
        else:
            self.spike_file = os.path.join(self.network_path, "simulation", "network-output-spikes-666.txt")

        if figure_path:
            self.figure_path = figure_path
        else:
            self.figure_path = os.path.join(self.network_path, "figures", "network-spike-raster.pdf")

        assert self.spike_file and os.path.isfile(self.spike_file), f"Input spike file {self.spike_file} does not exist"
        data = np.loadtxt(self.spike_file, delimiter="\t")
        self.spike_time = data[:, 0] / 1e3
        self.spike_neuron_id = data[:, 1].astype(int)

        self.snudda_load = SnuddaLoad(network_file=self.network_file)

    def get_colours(self, neuron_type):

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

    def plot_spike_raster(self, type_order=None):

        fig = plt.figure()
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

        ax.scatter(self.spike_time, spike_y, color=sc, s=1, linewidths=0.1)

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
