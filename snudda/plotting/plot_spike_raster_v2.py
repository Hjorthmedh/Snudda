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
        self.spike_time = data[:, 0]
        self.spike_neuron_id = data[:, 1]

        self.snudda_load = SnuddaLoad(network_file=self.network_file)

    def get_colours(self, neuron_type):

        colours = {"dSPN".lower(): (77. / 255, 151. / 255, 1.0),
                   "iSPN".lower(): (67. / 255, 55. / 255, 181. / 255),
                   "FS".lower(): (6. / 255, 31. / 255, 85. / 255),
                   "FSN".lower(): (6. / 255, 31. / 255, 85. / 255),
                   "ChIN".lower(): (252. / 255, 102. / 255, 0.0),
                   "LTS".lower(): (150. / 255, 63. / 255, 212. / 255)}

        return colours[neuron_type.lower()]

    def plot_spike_raster(self):

        fig = plt.figure()
        ax = fig.add_subplot(), fig.add_subplot()

        neuron_type_list = self.snudda_load.get_neuron_types(return_set=False)
        neuron_type_map = dict()

        for nt_idx, nt in enumerate(neuron_type_list):
            neuron_type_map[nt] = nt_idx

        # For each neuron, associate the number of the type it is
        neuron_type_idx = np.array([neuron_type_map[x] for x in neuron_type_list])
        neuron_order = np.argsort(neuron_type_idx)

        spike_y = np.take(neuron_order, self.spike_neuron_id)
        sc = [self.get_colours(x) for x in neuron_type_list]

        ax[0].scatter(self.time, spike_y, color=sc, s=1, linewidths=0.1)

        # Get position of labels
        unique_neuron_types = set(neuron_type_list)
        y_tick = []
        y_tick_label = []
        for nt in unique_neuron_types:
            y_tick_label.append(nt)
            y_tick = np.mean(neuron_order[np.where(neuron_type_list == nt)[0]])

        ax[0].set_xlabel('Time (s)')
        ax[0].set_yticks(y_tick)
        ax[0].set_yticklabels(y_tick_label)

        plt.savefig(self.figure_path, dpi=300)

        plt.ion()
        plt.show()


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser("Scatter plot")
    parser.add_argument("network_path", type=str, help="Network path")

    args = parser.parse_args()
    ps = SnuddaPlotSpikeRaster2(network_path=args.network_path)

