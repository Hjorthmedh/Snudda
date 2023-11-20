import os
import numpy as np

import matplotlib.pyplot as plt

from snudda.neurons import NeuronMorphologyExtended
from snudda.utils import SnuddaLoad, SnuddaLoadNetworkSimulation, snudda_parse_path


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
        self.simulation_data = SnuddaLoadNetworkSimulation(network_simulation_output_file=self.simulation_file,
                                                           network_path=self.network_path,
                                                           do_test=False)

    def load_morphology(self, neuron_id):
        morphology_file = snudda_parse_path(self.network_data["neurons"][neuron_id]["morphology"],
                                            snudda_data=self.snudda_data)
        print(f"Loading morphology: {morphology_file}")

        morphology = NeuronMorphologyExtended(swc_filename=morphology_file)
        return morphology

    def plot_neuron_voltage(self, neuron_id, time_range=None, section_id=None, axis=None,
                            show_plot=True, fig_name=None):

        if axis is None:
            fig = plt.figure()
            ax = fig.add_subplot()
        else:
            ax = axis

        morphology = self.load_morphology(neuron_id=neuron_id)

        time = self.simulation_data.get_time()
        voltage, sec_id_x, _ = self.simulation_data.get_data("voltage", neuron_id=neuron_id)

        volt = voltage[neuron_id]
        sec_id = sec_id_x[neuron_id][0]
        sec_x =  sec_id_x[neuron_id][1]

        if section_id is None:
            use_idx = np.arange(len(sec_id))
        else:
            use_mask = [s in section_id for s in sec_id]
            use_idx = np.where(sec_id)[0]

        ax.plot(time, volt[:, use_idx])

        soma_dist = np.zeros(sec_id.shape)

        for idx, (s_id, s_x) in enumerate(zip(sec_id, sec_x)):

            if s_id < 0:
                sec_type = 1  # soma
            else:
                sec_type = 3  # dend

            soma_dist[idx] = morphology.morphology_data["neuron"].sections[sec_type][s_id].soma_distanc_at(s_x)

        # TODO: We need to bin the data somehow...

        import pdb
        pdb.set_trace()

        if fig_name:
            plt.savefig(fig_name, dpi=300)

        if show_plot:
            plt.ion()
            plt.show()

        return ax


def cli():

    import argparse
    parser = argparse.ArgumentParser("Plot Neuron Voltage")
    parser.add_argument("network_path", type=str)
    parser.add_argument("neuron_id", type=int)
    parser.add_argument("--network_file", default=None, type=str)
    parser.add_argument("--simulation_file", default=None, type=str)
    parser.add_argument("--snudda_data", default=None, type=str)
    parser.add_argument("--pause", action="store_true", dest="pause")
    args = parser.parse_args()

    pnv = PlotNeuronVoltage(network_path=args.network_path,
                            network_file=args.network_file,
                            simulation_file=args.simulation_file,
                            snudda_data=args.snudda_data)

    pnv.plot_neuron_voltage(neuron_id=args.neuron_id)

    if args.pause:
        input("Press a key to continue")


if __name__ == "__main__":
    cli()
