import os
import numpy as np
from snudda.load import SnuddaLoad
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from snudda.neuron_morphology import NeuronMorphology

class PlotNetwork(object):

    def __init__(self, network_file):

        if os.path.isdir(network_file):
            network_file = os.path.join(network_file, "network-pruned-synapses.hdf5")

        self.network_file = network_file
        self.sl = SnuddaLoad(self.network_file)
        self.prototype_neurons = dict()

    def close(self):
        self.sl.close()

    def plot(self, plot_axon=True, plot_dendrite=True, plot_synapses=True,
             title=None, title_pad=None, show_axis=True,
             elev_azim=None, fig_name=None, dpi=600):

        if type(plot_axon) == bool:
            plot_axon = np.ones((self.sl.data["nNeurons"],), dtype=bool) * plot_axon

        if type(plot_dendrite) == bool:
            plot_dendrite = np.ones((self.sl.data["nNeurons"],), dtype=bool) * plot_dendrite

        assert len(plot_axon) == len(plot_dendrite) == len(self.sl.data["neurons"])

        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.gca(projection='3d')

        if "simulationOrigo" in self.sl.data:
            simulation_origo = self.sl.data["simulationOrigo"]
        else:
            simulation_origo = np.array([0, 0, 0])

        # Plot neurons
        for neuron_info, pa, pd in zip(self.sl.data["neurons"], plot_axon, plot_dendrite):
            neuron = self.load_neuron(neuron_info)
            neuron.plot_neuron(axis=ax,
                               plot_axon=pa,
                               plot_dendrite=pd,
                               soma_colour="black",
                               axon_colour="darksalmon",  #"maroon",
                               dend_colour="silver")   # Can also write colours as (0, 0, 0) -- rgb

        # Plot synapses
        if plot_synapses and "synapseCoords" in self.sl.data:
            ax.scatter(self.sl.data["synapseCoords"][:, 0],
                       self.sl.data["synapseCoords"][:, 1],
                       self.sl.data["synapseCoords"][:, 2], color=(0.1, 0.1, 0.1))

            plt.figtext(0.5, 0.20, f"{self.sl.data['nSynapses']} synapses", ha="center", fontsize=18)
            
        if elev_azim:
            ax.view_init(elev_azim[0], elev_azim[1])

        if not show_axis:
            plt.axis("off")

        plt.tight_layout()

        if title is None:
            title = ""

        if title_pad is not None:
            plt.rcParams['axes.titley'] = 0.95  # y is in axes-relative co-ordinates.
            plt.rcParams['axes.titlepad'] = title_pad  # pad is in points...

        plt.title(title, fontsize=18)

        # ax.dist = 8

        if fig_name is not None:
            fig_path = os.path.join(os.path.dirname(self.network_file), "figures", fig_name)
            if not os.path.exists(os.path.dirname(fig_path)):
                os.mkdir(os.path.dirname(fig_path))
            plt.savefig(fig_path, dpi=dpi, bbox_inches="tight")

        plt.ion()
        plt.show()

        return plt, ax


    def load_neuron(self, neuron_info=None, neuron_id=None):

        assert (neuron_info is None) + (neuron_id is None) == 1, "Specify neuron_info or neuron_id"

        if neuron_id is not None:
            print(f"Using id {neuron_id}")
            neuron_info = self.sl.data["neurons"][neuron_id]

        neuron_name = neuron_info["name"]

        if neuron_name not in self.prototype_neurons:
            self.prototype_neurons[neuron_name] = NeuronMorphology(name=neuron_name,
                                                                   swc_filename=neuron_info["morphology"])

        neuron = self.prototype_neurons[neuron_name].clone()
        neuron.place(rotation=neuron_info["rotation"],
                     position=neuron_info["position"])

        return neuron


if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser(description="Plot snudda network from file (hdf5)")
    parser.add_argument("networkFile", help="Network file (hdf5)", type=str)
    args = parser.parse_args()

    pn = PlotNetwork(args.networkFile)
    pn.plot()
