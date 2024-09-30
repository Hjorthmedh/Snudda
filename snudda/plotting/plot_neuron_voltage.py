import os
import numpy as np

import matplotlib.pyplot as plt

from snudda.neurons import NeuronMorphologyExtended
from snudda.utils import SnuddaLoad, SnuddaLoadSimulation, snudda_parse_path


class PlotNeuronVoltage:

    def __init__(self, network_path, network_file=None, simulation_file=None, snudda_data=None):

        self.network_path = network_path

        if network_file is None:
            self.network_file = os.path.join(network_path, "network-synapses.hdf5")
        else:
            self.network_file = network_file

        if simulation_file is None:
            self.simulation_file = os.path.join(network_path, "simulation", "output.hdf5")
        else:
            self.simulation_file = simulation_file

        print(f"Loading network data {self.network_file}")
        self.network_loader = SnuddaLoad(network_file=self.network_file)
        self.network_data = self.network_loader.data

        if snudda_data is None:
            self.snudda_data = self.network_data["SnuddaData"]
        else:
            self.snudda_data = snudda_data

        print(f"Loading simulation data {self.simulation_file}")
        self.simulation_data = SnuddaLoadSimulation(network_simulation_output_file=self.simulation_file,
                                                    network_path=self.network_path,
                                                    do_test=False)

    def load_morphology(self, neuron_id):
        morphology_file = snudda_parse_path(self.network_data["neurons"][neuron_id]["morphology"],
                                            snudda_data=self.snudda_data)
        print(f"Loading morphology: {morphology_file}")

        morphology = NeuronMorphologyExtended(swc_filename=morphology_file)
        return morphology

    def get_voltage_sec(self, neuron_id, section_id=None):

        voltage, sec_id_x, _ = self.simulation_data.get_data("voltage", neuron_id=neuron_id)

        volt = voltage[neuron_id]
        sec_id = sec_id_x[neuron_id][0]
        sec_x = sec_id_x[neuron_id][1]

        if section_id is None:
            use_idx = np.arange(len(sec_id))
        else:
            use_id = [section_id is None or s in section_id for s in sec_id]
            use_idx = np.where(use_id)[0]

        time = self.simulation_data.get_time()

        return time, volt[:, use_idx], sec_id[use_idx], sec_x[use_idx]

    def plot_neuron_voltage(self, neuron_id, section_id=None,
                            sliding_window_size=False,
                            title=None, axis=None,
                            show_plot=True, fig_name=None, fig_size=None):

        if axis is None:
            fig = plt.figure(figsize=fig_size)
            ax = fig.add_subplot()
        else:
            ax = axis

        time, volt, sec_id, sec_x = self.get_voltage_sec(neuron_id=neuron_id, section_id=section_id)

        soma_dist = self.get_soma_dist(neuron_id=neuron_id, section_id=sec_id, section_x=sec_x)
        label = [f"Neuron {neuron_id}, sec {sid}:{sx}, dist {sd*1e6:.1f} um" for (sid, sx, sd) in zip(sec_id, sec_x, soma_dist)]

        if sliding_window_size:

            print(f"Volt size {volt.shape}")

            mean_volt_list = []
            for i in range(0, volt.shape[1]):
                mean_volt_list.append(np.convolve(volt[:, i],
                                                  np.ones(sliding_window_size) / sliding_window_size, mode='valid'))

            mean_volt = np.vstack(mean_volt_list).T

            mean_time = np.convolve(time, np.ones(sliding_window_size) / sliding_window_size, mode='valid')

            # I do the same convolution for time, to make sure edges get correct time
            ax.plot(mean_time, mean_volt, label=label)

            title_info = f"{self.network_data['neurons'][neuron_id]['name']} ({neuron_id}) - sliding window {sliding_window_size}"

        else:
            ax.plot(time, volt, label=label)

            title_info = f"{self.network_data['neurons'][neuron_id]['name']} ({neuron_id})"

        if title:
            plt.title(f"{title} - {title_info}")
        else:
            plt.title(title_info)

        plt.legend()

        if fig_name:
            self.create_fig_dir(fig_name)
            plt.savefig(fig_name, dpi=300)

        if show_plot:
            plt.ion()
            plt.show()

        return ax

    def list_all_section_dist(self, neuron_id):

        morphology = self.load_morphology(neuron_id=neuron_id)

        for section in morphology.morphology_data["neuron"].sections[3].values():
            print(f"Section {section.section_id}, center dist {section.soma_distance_at(0.5)*1e6:.0f} um, length {section.section_length()*1e6:.0f} um")

    def create_fig_dir(self, fig_path):

        fig_dir = os.path.dirname(fig_path)
        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)

    def extract_mean_values(self, time, volt, time_bins):
        # This assumes that only the volt traces relevant are passed to the function

        volt_mean = []
        volt_std = []

        for (t_low, t_high) in time_bins:
            t_idx = np.where(np.logical_and(t_low <= time, time < t_high))[0]
            volt_mean.append(np.mean(volt[t_idx, :]))
            volt_std.append(np.std(volt[t_idx, :]))

        return np.array(volt_mean), np.array(volt_std)

    def plot_binned_neuron_voltage(self, neuron_id, time_bins, section_id=None, dist_bin_size=50e-6,
                                   title=None, axis=None, show_plot=True, fig_name=None, fig_size=None,
                                   y_range=None):

        if axis is None:
            fig = plt.figure(figsize=fig_size)
            ax = fig.add_subplot()
        else:
            ax = axis

        time, volt, sec_id, sec_x = self.get_voltage_sec(neuron_id=neuron_id, section_id=section_id)
        soma_dist = self.get_soma_dist(neuron_id=neuron_id, section_id=sec_id, section_x=sec_x)

        # First bin is soma only, bins after are for dendrites (with binwidth dist_bin_size)
        bin_idx = np.floor(soma_dist / dist_bin_size).astype(int) + 1
        bin_idx[np.where(soma_dist == 0)[0]] = 0

        bin_members = dict()
        for idx, bin_pos in enumerate(bin_idx):
            if bin_pos not in bin_members:
                bin_members[bin_pos] = [idx]
            else:
                bin_members[bin_pos].append(idx)

        for bin_pos in bin_members.keys():
            bin_members[bin_pos] = np.array(bin_members[bin_pos])

        n_bins = np.max(list(bin_members.keys()))+1
        volt_data = np.zeros((len(time_bins), n_bins))
        volt_std = np.zeros((len(time_bins), n_bins))

        for bin_pos, bin_idx in bin_members.items():
            volt_data[:, bin_pos], volt_std[:, bin_pos] = self.extract_mean_values(time=time, volt=volt[:, bin_idx], time_bins=time_bins)

        time_center = np.array([np.mean(t) for t in time_bins])

        for bin_pos, bin_idx in sorted(bin_members.items()):
            if bin_pos == 0:
                legend_text = "soma"
            else:
                legend_text = f"dist {np.mean(soma_dist[bin_idx])*1e6:.1f} \mum"

            if len(bin_idx) > 0:
                # t_jitter = 0.1*np.random.uniform(size=time_center.shape)
                t_jitter = 0.01*bin_pos
                ax.errorbar(time_center + t_jitter, volt_data[:, bin_pos], volt_std[:, bin_pos],
                            label=legend_text, marker='o')

        if y_range:
            ax.set_ylim(y_range)

        plt.legend(loc="upper left")

        title_info = f"{self.network_data['neurons'][neuron_id]['name']} ({neuron_id})"

        if title:
            plt.title(f"{title} - {title_info}")
        else:
            plt.title(title_info)

        fig.tight_layout()

        if fig_name:
            self.create_fig_dir(fig_name)
            plt.savefig(fig_name, dpi=300)

        if show_plot:
            plt.ion()
            plt.show()

        return ax

    def get_soma_dist(self, neuron_id, section_id, section_x):

        morphology = self.load_morphology(neuron_id=neuron_id)

        soma_dist = np.zeros(section_id.shape)

        for idx, (sec_id, sec_x) in enumerate(zip(section_id, section_x)):

            if sec_id < 0:
                sec_type = 1  # soma
                sec_id = 0
            else:
                sec_type = 3  # dend

            soma_dist[idx] = morphology.morphology_data["neuron"].sections[sec_type][sec_id].soma_distance_at(sec_x)

        return soma_dist


def cli():

    import argparse
    parser = argparse.ArgumentParser("Plot Neuron Voltage")
    parser.add_argument("network_path", type=str)
    parser.add_argument("neuron_id", type=int)
    parser.add_argument("--network_file", default=None, type=str)
    parser.add_argument("--simulation_file", default=None, type=str)
    parser.add_argument("--snudda_data", default=None, type=str)
    parser.add_argument("--pause", action="store_true", dest="pause")
    parser.add_argument("--time_bins", type=str, default="0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5")
    args = parser.parse_args()

    pnv = PlotNeuronVoltage(network_path=args.network_path,
                            network_file=args.network_file,
                            simulation_file=args.simulation_file,
                            snudda_data=args.snudda_data)

    pnv.plot_neuron_voltage(neuron_id=args.neuron_id)

    time_edges = np.array([float(x) for x in args.time_bins.split(",")])

    time_bins = list(zip(time_edges[:-1], time_edges[1:]))
    print(f"Using time bins {time_bins}")
    pnv.plot_binned_neuron_voltage(neuron_id=args.neuron_id, time_bins=time_bins)

    if args.pause:
        input("Press a key to continue")


if __name__ == "__main__":
    cli()
