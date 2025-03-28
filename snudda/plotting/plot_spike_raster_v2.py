import os
from collections import OrderedDict

import numpy as np

from snudda.utils.load import SnuddaLoad

import matplotlib.pyplot as plt

from snudda.utils.load_network_simulation import SnuddaLoadSimulation


class SnuddaPlotSpikeRaster2:

    def __init__(self, network_path, network_file=None, simulation_file=None, figure_path=None,
                 snudda_load=None, snudda_simulation_load=None):

        """ If you pass snudda_load and snudda_simulation_load those will be used, so you do not
            have to reload the data if you do multiple types of plots. """

        if network_path is not None:
            self.network_path = network_path
        elif snudda_load is not None:
            self.network_path = os.path.dirname(snudda_load.network_file)
        else:
            self.network_path = None

        if network_file:
            self.network_file = network_file
        elif network_path is not None:
            self.network_file = os.path.join(self.network_path, "network-synapses.hdf5")
        elif snudda_load:
            self.network_file = snudda_load.network_file
        else:
            self.network_file = None

        if simulation_file:
            self.simulation_file = simulation_file
        else:
            self.simulation_file = os.path.join(self.network_path, "simulation", "output.hdf5")

        if figure_path:
            self.figure_path = figure_path
        else:
            self.figure_path = os.path.join(self.network_path, "figures")

        if snudda_load:
            self.snudda_load = snudda_load
            assert network_file is None or snudda_load.network_file == self.network_file, \
                f"snudda_load refers to {snudda_load.network_file}, but user passed network_file={self.network_file}"
        else:
            self.snudda_load = SnuddaLoad(network_file=self.network_file, load_synapses=False)

        if snudda_simulation_load:
            self.snudda_simulation_load = snudda_simulation_load
            assert simulation_file is None \
                   or snudda_simulation_load.network_simulation_output_file_name == simulation_file, \
                f"snudda_simulation_load refers to {snudda_simulation_load.self.network_simulation_output_file_name}," \
                f" but user passed simulation_file={simulation_file}"
        else:
            self.snudda_simulation_load = SnuddaLoadSimulation(network_simulation_output_file=self.simulation_file)

        spike_data = self.snudda_simulation_load.merge_spikes()

        self.spike_time = spike_data[:, 0]
        self.spike_neuron_id = spike_data[:, 1].astype(int)
        self.traces = self.snudda_simulation_load.get_voltage()
        self.time = self.snudda_simulation_load.get_time()

    def make_figures_directory(self):

        if not os.path.isdir(self.figure_path):
            os.makedirs(self.figure_path)

    @staticmethod
    def get_colours(neuron_type):

        # TODO: Read colours from a JSON file

        colours = {"dSPN".lower(): (77. / 255, 151. / 255, 1.0),
                   "iSPN".lower(): (67. / 255, 55. / 255, 181. / 255),
                   "FS".lower(): (6. / 255, 31. / 255, 85. / 255),
                   "FSN".lower(): (6. / 255, 31. / 255, 85. / 255),
                   "ChIN".lower(): (252. / 255, 102. / 255, 0.0),
                   "LTS".lower(): (150. / 255, 63. / 255, 212. / 255)}

        if neuron_type.lower() in colours:
            return colours[neuron_type.lower()]
        else:
            return (0, 0, 0)

    def get_all_colours(self):

        neuron_type_list = self.snudda_load.get_neuron_types(return_set=False)
        neuron_colours = np.zeros((len(neuron_type_list), 3))

        for idx, nt in enumerate(neuron_type_list):
            neuron_colours[idx, :] = self.get_colours(nt)

        return neuron_colours

    # TODO: Add background colour to population units
    def plot_hist_raster(self, type_order=None, skip_time=0, end_time=None, fig_size=None, type_division=None, fig_file=None):
        # type_division: divides plot in two parts based on neuron type.
        #   Example 1: [["dspn","ispn"],["chin", "fsn","lts"]]
        #   Example 2: [["dspn","ispn","chin", "fsn","lts"],[]]

        self.make_figures_directory()
        # Gets a list of all the neurons' types
        neuron_type_list = self.snudda_load.get_neuron_types(return_set=False)
        neuron_population_unit_list = self.snudda_load.get_neuron_population_units(return_set=False)

        neuron_type_map = dict()

        if end_time is None:
            end_time = 1.02 * max(self.time)
        if type_order is None:
            unique_neuron_types = sorted(list(set(neuron_type_list)))
            type_order = unique_neuron_types
        else:
            unique_neuron_types = type_order + sorted(list(set(neuron_type_list) - set(type_order)))

        if type_division is None:
            type_division = [unique_neuron_types, []]

        for nt_idx, nt in enumerate(unique_neuron_types):
            neuron_type_map[nt] = nt_idx

        # For each neuron, associate the number of the type it is
        neuron_type_idx = np.array([neuron_type_map[x] for x in neuron_type_list])

        neuron_order = np.lexsort((neuron_population_unit_list, neuron_type_idx))

        # neuron_order = np.argsort(neuron_type_idx)
        neuron_order_lookup = np.zeros(neuron_order.shape)

        for idx, no in enumerate(neuron_order):
            neuron_order_lookup[no] = idx

        # new dict with cell type specific spikes
        type_dict = {'order': type_order}
        for t in type_order:
            type_dict[t] = [i for i, x in enumerate(neuron_type_list) if x.lower() == t.lower()]

        plt.rcParams.update({'font.size': 24,
                             'xtick.labelsize': 20,
                             'ytick.labelsize': 20,
                             'legend.loc': 'best'})
            
        # Prepare figure
        if not fig_size:
            fig_size = (10, 10)
        fig = plt.figure(figsize=fig_size)
        plt.rcParams.update({'font.size': 22})

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
            bin_range = np.arange(0, end_time - skip_time + bin_width, bin_width)
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
        ax.set_xlabel('Time (s)', fontsize=20)
        ax.set_yticks(y_tick, fontsize=20)
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

        if fig_file is None:
            fig_file = os.path.join(self.figure_path, "spike-histogram-raster.pdf")
        else:
            fig_file = os.path.join(self.figure_path, fig_file)

        print(f"Writing figure to {fig_file}")
        plt.tight_layout()
        self.make_figures_directory()
        plt.savefig(fig_file, dpi=300)

        plt.ion()
        plt.show()

        return ax

    def calculate_period_synchrony(self, period, neuron_id=None, time_range=None):

        spikes = self.snudda_simulation_load.get_spikes(neuron_id=neuron_id)
        all_spikes = []

        for s in spikes.values():

            sf = s.flatten()

            if time_range is not None:
                idx = np.where(np.logical_and(time_range[0] <= sf, sf <= time_range[1]))[0]
                all_spikes = all_spikes + list(sf[idx] % period)
            else:
                all_spikes = all_spikes + list(sf % period)

        phases = 2*np.pi/period * np.array(all_spikes)
        x = np.sum(np.cos(phases))
        y = np.sum(np.sin(phases))
        vs = np.sqrt(x**2 + y**2) / phases.size

        # Verify this is correct

        return vs

    def plot_spike_histogram_type(self, neuron_type, time_range=None, bin_size=50e-3, fig_size=None,
                                  fig_file=None, label_text=None, show_figure=True, n_core=None,
                                  linestyle="-", line_colours=None, linewidth=3,
                                  legend_loc="best", bbox_anchor=None, ax=None):

        self.make_figures_directory()

        plt.rcParams.update({'font.size': 24,
                             'xtick.labelsize': 20,
                             'ytick.labelsize': 20})

        assert type(neuron_type) == list, "neuron_type should be a list of neuron types"

        all_spikes = OrderedDict()
        neurons_of_type = OrderedDict()

        if time_range is None:
            time_range = (0, self.snudda_simulation_load.get_time()[-1])

        for nt in neuron_type:
            if n_core:
                neuron_id = [x for x, y
                             in self.snudda_load.get_centre_neurons_iterator(neuron_type=nt, n_neurons=n_core)]
            else:
                neuron_id = self.snudda_load.get_neuron_id_of_type(nt)

            spikes = self.snudda_simulation_load.get_spikes(neuron_id=neuron_id)
            neurons_of_type[nt] = neuron_id
            all_spikes[nt] = self.snudda_simulation_load.merge_spikes(spikes)[:, 0]

        bins = np.arange(time_range[0], time_range[1]+bin_size/2, bin_size)
        weights = [np.full(y.shape, 1/(len(x)*bin_size)) for x, y in zip(neurons_of_type.values(), all_spikes.values())]

        if label_text is None:
            label_text = ""

        if ax is None:
            fig = plt.figure(figsize=fig_size)
            ax = fig.add_subplot()

        if len(all_spikes.keys()) > 1:
            all_labels = [f"{label_text}{x}" for x in all_spikes.keys()]
        else:
            all_labels = [label_text]

        if line_colours is None:
            line_colours = [self.get_colours(x) for x in all_spikes.keys()]

        ax.hist(x=all_spikes.values(), bins=bins, weights=weights, linewidth=linewidth, linestyle=linestyle,
                histtype="step", color=line_colours,
                label=all_labels)

        plt.xlabel("Time (s)", fontsize=20)
        plt.ylabel("Frequency (Hz)", fontsize=20)
        ax.legend(loc=legend_loc, bbox_to_anchor=bbox_anchor)

        plt.tight_layout()

        if fig_file:
            fig_path = os.path.join(self.figure_path, fig_file)
            print(f"Writing figure {fig_path}")
            self.make_figures_directory()

            plt.savefig(fig_path)

        if show_figure:
            plt.ion()
            plt.show()

        return ax

    def plot_spike_histogram(self, population_id=None, neuron_type=None,
                             skip_time=0, end_time=None, fig_size=None, bin_size=50e-3,
                             fig_file=None, ax=None, label_text=None, show_figure=True, save_figure=True, colour=None,
                             linestyle="-", legend_loc="best", title=None, bbox_anchor=None):

        if population_id is None:
            population_id = self.snudda_load.get_neuron_population_units(return_set=True)

        if neuron_type is not None:

            if type(neuron_type) != list:
                neuron_type = [neuron_type]

            keep_neuron_id = set()
            for nt in neuron_type:
                keep_neuron_id |= set(self.snudda_load.get_neuron_id_of_type(neuron_type=nt))

            print(f"Processing {len(keep_neuron_id)} neurons of type {neuron_type}")
        else:
            keep_neuron_id = None

        self.make_figures_directory()

        plt.rcParams.update({'font.size': 24,
                             'xtick.labelsize': 20,
                             'ytick.labelsize': 20},
                            )

        if ax is None:
            fig = plt.figure(figsize=fig_size)
            ax = fig.add_subplot()
        
        pop_members = dict()
        pop_spikes = dict()

        if np.issubdtype(type(population_id), np.integer):
            population_id = np.array([population_id])

        for pid in population_id:
            members = self.snudda_load.get_population_unit_members(pid)

            if keep_neuron_id is not None:
                members = np.array(list(set(members) & keep_neuron_id))

            if len(members) == 0:
                # No members, skip it
                continue

            pop_members[pid] = members

            spikes = self.snudda_simulation_load.get_spikes(pop_members[pid])
            pop_spikes[pid] = self.snudda_simulation_load.merge_spikes(spikes)[:, 0]

        if end_time is None:
            end_time = self.snudda_simulation_load.get_time()[-1]

        try:
            bins = np.arange(skip_time, end_time+bin_size/2, bin_size)
            weights = [np.full(y.shape, 1/(len(x)*bin_size)) for x, y in zip(pop_members.values(), pop_spikes.values())]
        except:
            import traceback
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()

        if label_text is None:
            label_text = ""

        N, bins, patches = ax.hist(x=pop_spikes.values(), bins=bins, weights=weights, linewidth=3, linestyle=linestyle,
                                   histtype="step", color=colour,
                                   label=[f"{label_text}{x}" for x in pop_spikes.keys()])

        if type(colour) == list:
            for patch, col in zip(patches, colour):
                patch[0].set_facecolor(col)

        plt.xlabel("Time (s)", fontsize=20)
        plt.ylabel("Frequency (Hz)", fontsize=20)
        ax.legend(loc=legend_loc, bbox_to_anchor=bbox_anchor)

        if title:
            plt.title(title)

        if fig_file is None:
            fig_file = os.path.join(self.figure_path,
                                    f"spike-frequency-pop-units{'-'.join([f'{x}' for x in pop_members.keys()])}.png")
        else:
            fig_file = os.path.join(self.figure_path, fig_file)

        if save_figure:
            print(f"Saving figure {fig_file}")
            plt.tight_layout()
            self.make_figures_directory()

            plt.savefig(fig_file, dpi=300)

        if show_figure:
            plt.ion()
            plt.show()

        return ax

    # Use this to plot a histogram for an arbitrary group of neurons specified with neuron_id
    def plot_group_spike_histogram(self, neuron_id=None,
                                    skip_time=0, end_time=None, fig_size=None, bin_size=50e-3,
                                    fig_file=None, ax=None, label_text=None, show_figure=True, save_figure=True, colour=None,
                                    linestyle="-", legend_loc="best", title=None):

        self.make_figures_directory()

        plt.rcParams.update({'font.size': 24,
                             'xtick.labelsize': 20,
                             'ytick.labelsize': 20,
                             'legend.loc': legend_loc})

        if ax is None:
            fig = plt.figure(figsize=fig_size)
            ax = fig.add_subplot()

        spikes = self.snudda_simulation_load.get_spikes(neuron_id)
        spikes = self.snudda_simulation_load.merge_spikes(spikes)[:, 0]

        if end_time is None:
            end_time = self.snudda_simulation_load.get_time()[-1]

        bins = np.arange(skip_time, end_time + bin_size / 2, bin_size)
        weights = 1 / (len(neuron_id) * bin_size)

        if label_text is None:
            label_text = ""

        try:
            N, bins, patches = ax.hist(x=spikes, bins=bins, weights=np.full(spikes.shape, weights), linewidth=3, linestyle=linestyle,
                                       histtype="step", color=colour,
                                       label=label_text)
        except:
            import traceback
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()

        if type(colour) == list:
            for patch, col in zip(patches, colour):
                patch[0].set_facecolor(col)

        plt.xlabel("Time (s)", fontsize=20)
        plt.ylabel("Frequency (Hz)", fontsize=20)
        ax.legend()

        if title:
            plt.title(title)

        if fig_file is None:
            fig_file = os.path.join(self.figure_path,
                                    f"spike-frequency-for-group.pdf")
        else:
            fig_file = os.path.join(self.figure_path, fig_file)

        if save_figure:
            self.make_figures_directory()

            print(f"Saving figure {fig_file}")
            plt.tight_layout()
            plt.savefig(fig_file, dpi=300)

        if show_figure:
            plt.ion()
            plt.show()

        return ax

    def plot_spike_raster(self, type_order=None, skip_time=0, end_time=None, fig_size=None, fig_file=None,
                          time_range=None, title=None):

        self.make_figures_directory()

        plt.rcParams.update({'font.size': 24,
                             'xtick.labelsize': 20,
                             'ytick.labelsize': 20,
                             'legend.loc': 'best'})
        
        fig = plt.figure(figsize=fig_size)
        ax = fig.add_subplot()

        # Gets a list of all the neurons' types
        neuron_type_list = self.snudda_load.get_neuron_types(return_set=False)
        neuron_population_unit_list = self.snudda_load.get_neuron_population_units(return_set=False)

        neuron_type_map = dict()

        if type_order is None:
            unique_neuron_types = sorted(list(set(neuron_type_list)))
        else:
            unique_neuron_types = type_order + sorted(list(set(neuron_type_list) - set(type_order)))

        for nt_idx, nt in enumerate(unique_neuron_types):
            neuron_type_map[nt] = nt_idx

        # For each neuron, associate the number of the type it is
        neuron_type_idx = np.array([neuron_type_map[x] for x in neuron_type_list])
        # neuron_order = np.argsort(neuron_type_idx)
        neuron_order = np.lexsort((neuron_population_unit_list, neuron_type_idx))

        neuron_order_lookup = np.zeros(neuron_order.shape)

        # for idx, no in enumerate(neuron_order):
        #     neuron_order_lookup[no] = idx

        idx = 0
        for no in neuron_order:
            # Skip the virtual neurons
            if not self.snudda_load.data["neurons"][no]["virtual_neuron"]:
                neuron_order_lookup[no] = idx
                idx += 1

        spike_y = np.take(neuron_order_lookup, self.spike_neuron_id)

        colour_lookup = self.get_all_colours()
        sc = np.zeros((len(spike_y), 3))

        for i in range(0, 3):
            sc[:, i] = np.take(colour_lookup[:, i], self.spike_neuron_id)

        ax.scatter(self.spike_time - skip_time, spike_y, color=sc, s=5, linewidths=0.1)

        # Optionally we should also show the virtual neuron spikes...


        # Get position of labels
        unique_neuron_types = set(neuron_type_list)
        y_tick = []
        y_tick_label = []
        for nt in unique_neuron_types:
            y_tick_label.append(nt)
            y_tick.append(np.mean(neuron_order_lookup[np.where([x == nt for x in neuron_type_list])[0]]))

        ax.invert_yaxis()
        ax.set_xlabel('Time (s)', fontsize=20)
        ax.set_yticks(y_tick)
        ax.set_yticklabels(y_tick_label, fontsize=20)

        if skip_time or end_time:
            x_lim = ax.get_xlim()
            x_lim = (0, self.time[1] - skip_time)
            if end_time:
                x_lim = (x_lim[0], end_time)
            ax.set_xlim(x_lim)
        else:
            x_lim = ax.get_xlim()
            x_lim = (self.time[0], self.time[-1])
            ax.set_xlim(x_lim)

        assert not ((skip_time or end_time) and time_range), \
            f"time_range can not be specified with skip_time and end_time"

        if time_range:
            ax.set_xlim(time_range)

        if fig_file is None:
            fig_file = "spike_raster.pdf"

        if title is not None:
            ax.set_title(title)

        fig_file = os.path.join(self.figure_path, fig_file)

        print(f"Saving figure to {fig_file}")
        plt.tight_layout()
        self.make_figures_directory()
        plt.savefig(fig_file, dpi=300)

        plt.ion()
        plt.show()

    def plot_firing_frequency_distribution(self, time_range=None, figure_name=None, bins=20, title=None):

        neuron_types = self.snudda_load.get_neuron_types(return_set=True)

        plt.figure()

        time = self.snudda_simulation_load.get_time()

        for nt in neuron_types:

            neuron_id = self.snudda_load.get_neuron_id_of_type(neuron_type=nt, include_virtual=False)
            spikes = self.snudda_simulation_load.get_spikes(neuron_id=neuron_id)

            if time_range is None:
                freq = [s.size/(time[-1] - time[0]) for s in spikes.values()]
            else:
                freq = [len(np.where(np.logical_and(time_range[0] <= s, s <= time_range[1]))[0])
                        / (time_range[1] - time_range[0])
                        for s in spikes.values()]

            if not isinstance(bins, int):
                if np.max(freq) > max(bins):
                    raise ValueError(f"Max frequency {np.max(freq)} larger than bin range specified {max(bins)}.")

            colour = SnuddaPlotSpikeRaster2.get_colours(nt)
            count, bin = np.histogram(freq, bins=bins)
            plt.stairs(count, bin, label=nt, color=colour, linewidth=3)

        plt.legend()
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Count")

        if title is not None:
            plt.title(title)

        plt.tight_layout()

        if figure_name is not None:
            self.make_figures_directory()

            plt.savefig(os.path.join(self.figure_path, figure_name))

        plt.ion()
        plt.show()

    def plot_population_frequency(self, population_id, time_ranges=None):

        raise NotImplementedError()

        if population_id is None:
            population_id = sorted(list(self.snudda_load.get_neuron_population_units(return_set=True)))
            label = [str(x) for x in population_id]

        if time_ranges is None:
            time_ranges = [(0, np.max(self.snudda_simulation_load.get_time()))]

        x_labels = [f"{x[0]} -- {x[1]}" for x in time_ranges]

        data = dict()

        for pop_id in population_id:
            neuron_id = self.snudda_load.get_population_unit_members(population_unit=pop_id)
            freq_table = self.snudda_simulation_load.get_frequency(neuron_id=neuron_id, time_ranges=time_ranges)
            avg_freq = np.sum(freq_table, axis=0).flatten()

            data[str(pop_id)] = avg_freq

        self.plot_grouped_bars(legend_labels_and_data=data)

    def plot_grouped_bars(self, legend_labels_and_data, x_labels, y_unit_label, title):

        raise NotImplementedError()

        # Derived from matplotlib example

        # x_labels = ("Adelie", "Chinstrap", "Gentoo")
        # legend_labels_and_data = {
        #     'Bill Depth': (18.35, 18.43, 14.98),
        #     'Bill Length': (38.79, 48.83, 47.50),
        #     'Flipper Length': (189.95, 195.82, 217.19),
        # }

        x = np.arange(len(x_labels))  # the label locations
        width = 0.25  # the width of the bars
        multiplier = 0

        fig, ax = plt.subplots(layout='constrained')

        for attribute, measurement in legend_labels_and_data.items():
            offset = width * multiplier
            rects = ax.bar(x + offset, measurement, width, label=attribute)
            ax.bar_label(rects, padding=3)
            multiplier += 1

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_ylabel(y_unit_label)
        ax.set_title(title)
        ax.set_xticks(x + width, x_labels)
        ax.legend(loc='upper left', ncols=3)
        ax.set_ylim(0, 250)

        plt.show()


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("Scatter plot")
    parser.add_argument("network_path", type=str, help="Network path")

    args = parser.parse_args()
    ps = SnuddaPlotSpikeRaster2(network_path=args.network_path)

    type_order = ["FS", "dSPN", "LTS", "iSPN", "ChIN"]
    ps.plot_spike_raster(type_order)
    ps.plot_spike_histogram(type_order)
