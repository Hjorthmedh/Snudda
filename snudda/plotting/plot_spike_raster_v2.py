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
            self.simulation_file = os.path.join(self.network_path, "simulation", "output.hdf5")

        if figure_path:
            self.figure_path = figure_path
        else:
            self.figure_path = os.path.join(self.network_path, "figures", "network-spike-raster.png")

        self.snudda_load = SnuddaLoad(network_file=self.network_file)

        self.snudda_simulation_load = SnuddaLoadNetworkSimulation(network_simulation_output_file=self.simulation_file)
        spike_data = self.snudda_simulation_load.merge_spikes()

        self.spike_time = spike_data[:, 0]
        self.spike_neuron_id = spike_data[:, 1].astype(int)
        self.traces = self.snudda_simulation_load.get_voltage()
        self.time = self.snudda_simulation_load.get_time()

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

    def plot_hist_raster(self, type_order=None, skip_time=0, end_time=None, fig_size=None, type_division=None):
        # type_division: divides plot in two parts based on neuron type.
        #   Example 1: [["dspn","ispn"],["chin", "fsn","lts"]]
        #   Example 2: [["dspn","ispn","chin", "fsn","lts"],[]]

        self.make_figures_directory()
        # Gets a list of all the neurons' types
        neuron_type_list = self.snudda_load.get_neuron_types(return_set=False)
        neuron_type_map = dict()
        if end_time is None:
            end_time = 1.02 * max(self.time)
        if type_order is None:
            unique_neuron_types = np.unique(neuron_type_list)
            type_order = unique_neuron_types
        else:
            unique_neuron_types = type_order + list(set(neuron_type_list) - set(type_order))

        if type_division is None:
            type_division = [unique_neuron_types, []]

        for nt_idx, nt in enumerate(unique_neuron_types):
            neuron_type_map[nt] = nt_idx

        # For each neuron, associate the number of the type it is
        neuron_type_idx = np.array([neuron_type_map[x] for x in neuron_type_list])
        neuron_order = np.argsort(neuron_type_idx)
        neuron_order_lookup = np.zeros(neuron_order.shape)

        for idx, no in enumerate(neuron_order):
            neuron_order_lookup[no] = idx

        # new dict with cell type specific spikes
        type_dict = {'order': type_order}
        for t in type_order:
            type_dict[t] = [i for i, x in enumerate(neuron_type_list) if x.lower() == t.lower()]

        # Prepare figure
        if not fig_size:
            fig_size = (10, 10)
        fig = plt.figure(figsize=fig_size)
        r = 4
        grid = plt.GridSpec(r, r, hspace=0, wspace=0)
        ax = fig.add_subplot(grid[2:, :])
        if (len(type_division[0]) > 0) & (len(type_division[1]) > 0):
            atop = fig.add_subplot(grid[0, :])
            atop2 = fig.add_subplot(grid[1, :])
        elif (len(type_division[0]) == 0) & (len(type_division[1]) == 0):
            print("No neuron type to show in histogram due to empty variable type_division.")
        elif (len(type_division[0]) == 0) | (len(type_division[1]) == 0):
            atop = fig.add_subplot(grid[0:2, :])
            atop2 = atop

        # Plot raster plot
        spike_y = np.take(neuron_order_lookup, self.spike_neuron_id)
        colour_lookup = self.get_all_colours()
        sc = np.zeros((len(spike_y), 3))
        for i in range(0, 3):
            sc[:, i] = np.take(colour_lookup[:, i], self.spike_neuron_id)
        ax.scatter(self.spike_time - skip_time, spike_y, color=sc, s=5, linewidths=0.1)

        max0 = 0
        max1 = 0
        spike2type = np.take(neuron_type_list, self.spike_neuron_id)

        # Spike rates
        for t in type_order:
            cell_spike_rates = np.array([])
            spikes_of_type = [i for i, x in enumerate(spike2type) if x.lower() == t.lower()]
            pruned_spikes = self.spike_time[spikes_of_type] - skip_time
            pruned_spikes = pruned_spikes[pruned_spikes > 0]
            spike_y_t = spike_y[spikes_of_type]
            spike_y_t = spike_y_t[pruned_spikes > 0]
            uspike_y_t = np.unique(spike_y_t)

            for cellid in uspike_y_t:
                cell_spike_rate = len(pruned_spikes[cellid == spike_y_t]) / (end_time - skip_time)
                cell_spike_rates = np.append(cell_spike_rates, cell_spike_rate)
            num_of_type = len(type_dict[t])
            cell_spike_rates = np.append(cell_spike_rates, np.zeros(num_of_type - len(cell_spike_rates)))

            print("NeuronType: {0}, Mean: {1}Hz Stdev: {2}Hz".format(t, np.mean(cell_spike_rates),
                                                                     np.std(cell_spike_rates)))

        # histogram
        for t in type_order:
            spikes_of_type = [i for i, x in enumerate(spike2type) if x.lower() == t.lower()]
            pruned_spikes = self.spike_time[spikes_of_type] - skip_time
            pruned_spikes = pruned_spikes[pruned_spikes > 0]
            num_of_type = len(type_dict[t])  #
            nspikes = len(pruned_spikes)
            bin_width = 0.05  # 10  # ms
            bin_range = np.arange(0, end_time + skip_time + bin_width, bin_width)
            if (nspikes > 0) & (t.lower() in type_division[0]):
                counts0, bins0, bars0 = atop.hist(pruned_spikes,
                                                  bins=bin_range,
                                                  range=(skip_time, end_time),
                                                  density=0,
                                                  color=sc[spikes_of_type][0],
                                                  alpha=1.0,
                                                  histtype='step',
                                                  weights=np.ones_like(pruned_spikes) / bin_width / num_of_type)

                max0 = max(np.append(counts0, max0))
            elif (nspikes > 0) & (t.lower() in type_division[1]):
                counts1, bins1, bars1 = atop2.hist(pruned_spikes,
                                                   bins=bin_range,
                                                   range=(skip_time, end_time),
                                                   density=0,
                                                   color=sc[spike2type == t][0],
                                                   alpha=1.0,
                                                   histtype='step',
                                                   weights=np.ones_like(pruned_spikes) / bin_width / num_of_type)

                max1 = max(np.append(counts1, max1))
            spike_rate = len(pruned_spikes) / num_of_type / (end_time - skip_time)
            print("{0}: {1} Hz".format(t, spike_rate))
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
            x_lim = (0, x_lim[1])
            if end_time:
                x_lim = (x_lim[0], end_time)
            ax.set_xlim(x_lim)
        atop.set_xlim(x_lim)
        atop2.set_xlim(x_lim)
        atop.set_xticklabels([])
        atop.set_ylabel('Mean spikes/s')
        if not os.path.isdir(os.path.basename(self.figure_path)):
            os.makedirs(os.path.basename(self.figure_path))

        print(f"Writing figure to {self.figure_path}")
        plt.savefig(self.figure_path, dpi=300)

        plt.ion()
        plt.show()

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

        ax.scatter(self.spike_time - skip_time, spike_y, color=sc, s=5, linewidths=0.1)

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
            x_lim = (0, x_lim[1])
            if end_time:
                x_lim = (x_lim[0], end_time)
            ax.set_xlim(x_lim)

        if not os.path.isdir(os.path.basename(self.figure_path)):
            os.makedirs(os.path.basename(self.figure_path))

        print(f"Saving figure to {self.figure_path}")
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
