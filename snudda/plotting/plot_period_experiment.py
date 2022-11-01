import os
import numpy as np
import matplotlib.pyplot as plt

from snudda.plotting.plot_spike_raster_v2 import SnuddaPlotSpikeRaster2


class PlotPeriodExperiment(SnuddaPlotSpikeRaster2):

    def calculate_period_histogram_mod(self, period, neuron_id, time_range,
                                       exclude_depolarisation_blocked_neurons=False,
                                       n_bins=10):

        if exclude_depolarisation_blocked_neurons:
            if neuron_id is None:
                neuron_id = list(self.snudda_simulation_load.network_simulation_file["neurons"].keys())

            if self.snudda_simulation_load.depolarisation_block is None:
                self.snudda_simulation_load.depolarisation_block = self.snudda_simulation_load.check_depolarisation_block()

            bad_cells = sorted(list(set([x for x, ts, te in self.snudda_simulation_load.depolarisation_block])))
            neuron_id = np.array([x for x in neuron_id if x not in bad_cells])

        spikes = self.snudda_simulation_load.get_spikes(neuron_id=neuron_id)
        all_spikes = []

        n_spike_trains = len(spikes)

        for s in spikes.values():

            sf = s.flatten()

            if time_range is not None:
                idx = np.where(np.logical_and(time_range[0] <= sf, sf <= time_range[1]))[0]
                all_spikes = all_spikes + list(sf[idx] % period)
            else:
                all_spikes = all_spikes + list(sf % period)

        counts, bins = np.histogram(all_spikes, bins=n_bins)

        if time_range:
            t = time_range[1] - time_range[0]
        else:
            t = self.time[-1] - self.time[0]

        freq = counts/(n_spike_trains * (bins[1] - bins[0]) * t/period)

        return freq, bins

    def plot_period_histogram_mod(self, period, neuron_id=None, time_range=None, n_bins=20,
                                  fig_file=None, ax=None, fig_size=None, label=None, color=None,
                                  show_figure=True, exclude_depolarisation_blocked_neurons=False,
                                  save_figure=True, linestyle="-", legend_loc="best"):

        self.make_figures_directory()

        plt.rcParams.update({'font.size': 24,
                             'xtick.labelsize': 20,
                             'ytick.labelsize': 20,
                             'legend.loc': legend_loc})

        if ax is None:
            fig = plt.figure(figsize=fig_size)
            ax = fig.add_subplot()

        freq, bins = self.calculate_period_histogram_mod(period=period, neuron_id=neuron_id, time_range=time_range,
                                                         n_bins=n_bins,
                                                         exclude_depolarisation_blocked_neurons=exclude_depolarisation_blocked_neurons)

        ax.stairs(freq, bins, label=label, color=color, linewidth=3, linestyle=linestyle)

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Frequency (Hz)")
        ax.legend()

        if fig_file is None:
            fig_file = os.path.join(self.figure_path, "spike-period-histogram.pdf")
        else:
            fig_file = os.path.join(self.figure_path, fig_file)

        plt.tight_layout()

        if save_figure:
            print(f"Writing figure to {fig_file}")
            plt.savefig(fig_file, dpi=300)

        if show_figure:
            plt.ion()
            plt.show()

        # ax
        return ax, freq, bins

    def plot_period_voltage(self, period, neuron_type, ax=None, fig_size=None, figure_name=None, show_plot=True):

        """ This plots the average voltage for the neurons during the spikes, split into two groups. Those that
            had a presynaptic neuron spike before, and those that did not. """

        post_id_list = self.snudda_load.get_neuron_id_of_type(neuron_type=neuron_type)

        pre_id_lookup = dict()

        spike_with_pre_spike = []
        spike_no_pre_spike = []
        no_spike_with_pre_spike = []
        no_spike_no_pre_spike = []

        time = self.snudda_simulation_load.get_time()
        period_start = np.arange(0, time[-1], period)
        period_end = period_start + period

        for post_id in post_id_list:
            synapses, _ = self.snudda_load.find_synapses(post_id=post_id)
            pre_id_lookup[post_id] = set(synapses[:, 0])

        for post_id in post_id_list:
            spikes = self.snudda_simulation_load.get_spikes(neuron_id=post_id)
            pre_spikes = self.snudda_simulation_load.get_spikes(neuron_id=pre_id_lookup[post_id])

            for p_start, p_end in zip(period_start, period_end):

                # TODO: TO BE CONTINUED...
                # Extract all the traces where a presynaptic neuron spikes.
                pass

        assert False, "Not implemented yet"




        if ax is None:
            fig = plt.figure(figsize=fig_size)
            ax = fig.add_subplot()