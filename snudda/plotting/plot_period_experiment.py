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

        ax.stairs(freq, bins * 1e3, label=label, color=color, linewidth=3, linestyle=linestyle)

        ax.set_xlabel("Time (ms)")
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

    def plot_period_voltage(self, period, neuron_type, ax=None, fig_size=None, figure_name=None, show_plot=True,
                            time_range=None, show_pre_histogram=True):

        """ This plots the average voltage for the neurons during the spikes, split into two groups. Those that
            had a presynaptic neuron spike before, and those that did not. """

        if ax is None:
            fig = plt.figure(figsize=fig_size)
            ax = fig.add_subplot()

        post_id_list = self.snudda_load.get_neuron_id_of_type(neuron_type=neuron_type)

        pre_id_lookup = dict()

        spike_with_pre_spike = []
        spike_no_pre_spike = []
        no_spike_with_pre_spike = []
        no_spike_no_pre_spike = []

        time = np.around(self.snudda_simulation_load.get_time(), 8)
        period_start = np.around(np.arange(0, time[-1], period), 8)  # Use around to avoid rounding errors
        period_end = np.around(period_start + period, 8)

        n_periods = int(np.ceil(time[-1] / period))
        spike_in_period = np.zeros((self.snudda_load.data["nNeurons"], n_periods))

        for neuron_id in self.snudda_load.get_neuron_id():
            spikes = self.snudda_simulation_load.get_spikes(neuron_id=neuron_id)

            for spike in spikes.flatten():
                p_idx = int(spike / period)
                spike_in_period[neuron_id, p_idx] += 1

        for post_id in post_id_list:
            synapses, _ = self.snudda_load.find_synapses(post_id=post_id)
            pre_id_lookup[post_id] = np.array(sorted(list(set(synapses[:, 0]))))

        for post_id in post_id_list:

            pre_id = pre_id_lookup[post_id]
            has_pre_spike = np.sum(spike_in_period[pre_id, :], axis=0)
            has_post_spike = spike_in_period[post_id, :]

            voltage = self.snudda_simulation_load.get_voltage(neuron_id=post_id).flatten()

            for t_start, t_end, pre_spikes, post_spikes \
                    in zip(period_start, period_end, has_pre_spike, has_post_spike):

                if time_range and not (time_range[0] <= t_start and t_end <= time_range[1]):
                    continue

                t_idx = np.where(np.logical_and(t_start <= time, time < t_end))[0]

                if len(t_idx) > 2000:

                    print("Tell me why?!")
                    import pdb
                    pdb.set_trace()

                if pre_spikes:
                    if post_spikes:
                        spike_with_pre_spike.append(voltage[t_idx])
                    else:
                        no_spike_with_pre_spike.append(voltage[t_idx])
                else:
                    if post_spikes:
                        spike_no_pre_spike.append(voltage[t_idx])
                    else:
                        no_spike_no_pre_spike.append(voltage[t_idx])

        t_idx = np.where(np.logical_and(period_start[0] <= time, time < period_end[0]))[0]
        t = time[t_idx] - time[t_idx[0]]

        if len(spike_with_pre_spike) > 0:
            spike_with_pre_spike_mat = np.stack(spike_with_pre_spike, axis=1)
            ax.plot(t * 1e3, 1e3 * np.mean(spike_with_pre_spike_mat, axis=1),
                    label=f"{neuron_type} spike, pre-spike (n={spike_with_pre_spike_mat.shape[1]})")

        if len(no_spike_with_pre_spike) > 0:
            no_spike_with_pre_spike_mat = np.stack(no_spike_with_pre_spike, axis=1)
            ax.plot(t * 1e3, 1e3 * np.mean(no_spike_with_pre_spike_mat, axis=1),
                    label=f"No {neuron_type} spike, pre-spike (n={no_spike_with_pre_spike_mat.shape[1]})")

        if len(spike_no_pre_spike) > 0:
            spike_no_pre_spike_mat = np.stack(spike_no_pre_spike, axis=1)
            ax.plot(t * 1e3, 1e3 * np.mean(spike_no_pre_spike_mat, axis=1),
                    label=f"{neuron_type} spike, no pre-spike (n={spike_no_pre_spike_mat.shape[1]})")

        if len(no_spike_no_pre_spike) > 0:
            no_spike_no_pre_spike_mat = np.stack(no_spike_no_pre_spike, axis=1)
            ax.plot(t*1e3, 1e3*np.mean(no_spike_no_pre_spike_mat, axis=1),
                    label=f"No {neuron_type} spike, no pre-spike (n={no_spike_no_pre_spike_mat.shape[1]})")

        ax.legend()
        plt.xlabel("Time (ms)")
        plt.ylabel("Voltage (mV)")

        if figure_name:
            fig_path = os.path.join(self.figure_path, figure_name)
            plt.savefig(fig_path)

        if show_plot:
            plt.ion()
            plt.show()

        # if show_pre_histogram:
        #     plt.figure(figsize=fig_size)

        return ax
