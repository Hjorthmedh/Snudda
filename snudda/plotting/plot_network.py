import os
import numpy as np
from snudda.utils.load import SnuddaLoad
import matplotlib.pyplot as plt

from snudda.neurons.neuron_morphology import NeuronMorphology


class PlotNetwork(object):

    def __init__(self, network):

        if os.path.isdir(network):
            network_file = os.path.join(network, "network-synapses.hdf5")
        else:
            network_file = network

        self.network_file = network_file
        self.network_path = os.path.dirname(self.network_file)

        self.sl = SnuddaLoad(self.network_file)
        self.prototype_neurons = dict()

    def close(self):
        self.sl.close()

    def plot(self, plot_axon=True, plot_dendrite=True, plot_synapses=True,
             title=None, title_pad=None, show_axis=True,
             elev_azim=None, fig_name=None, dpi=600,
             colour_population_unit=False):

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

        if "populationUnit" in self.sl.data and colour_population_unit:
            population_unit = self.sl.data["populationUnit"]
            pop_units = sorted(list(set(population_unit)))
            cmap = plt.get_cmap('tab20', len(pop_units))
            colour_lookup_helper = dict()

            for idx, pu in enumerate(population_unit):
                if pu > 0:
                    colour_lookup_helper[idx] = cmap(pu)
                else:
                    colour_lookup_helper[idx] = 'lightgrey'

            colour_lookup = lambda x: colour_lookup_helper[x]
        else:
            colour_lookup = lambda x: 'black'

        # Plot neurons
        for neuron_info, pa, pd in zip(self.sl.data["neurons"], plot_axon, plot_dendrite):

            soma_colour = colour_lookup(neuron_info["neuronID"])
            neuron = self.load_neuron(neuron_info)
            neuron.plot_neuron(axis=ax,
                               plot_axon=pa,
                               plot_dendrite=pd,
                               soma_colour=soma_colour,
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

        self.equal_axis(ax)

        if fig_name is not None:
            fig_path = os.path.join(self.network_path, "figures", fig_name)
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

    def plot_populations(self):

        fig = plt.figure(figsize=(6, 6.5))
        ax = fig.gca(projection='3d')

        assert "populationUnit" in self.sl.data

        population_unit = self.sl.data["populationUnit"]
        pop_units = sorted(list(set(population_unit)))
        cmap = plt.get_cmap('tab20', len(pop_units))
        neuron_colours = []

        for idx, pu in enumerate(population_unit):
            if pu > 0:
                neuron_colours.append(list(cmap(pu)))
            else:
                neuron_colours.append([0.7, 0.7, 0.7, 1.0])

        neuron_colours = np.array(neuron_colours)
        positions = self.sl.data["neuronPositions"]

        ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c=neuron_colours, marker='o', s=20)

        self.equal_axis(ax)

    def equal_axis(self, ax):

        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        z_min, z_max = ax.get_zlim()

        x_mean = (x_min + x_max)/2
        y_mean = (y_min + y_max)/2
        z_mean = (z_min + z_max)/2

        max_half_width = np.max([x_max-x_min, y_max-y_min, z_max-z_min])/2

        ax.set_xlim(x_mean-max_half_width, x_mean+max_half_width)
        ax.set_ylim(y_mean-max_half_width, y_mean+max_half_width)
        ax.set_zlim(z_mean-max_half_width, z_mean+max_half_width)

if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser(description="Plot snudda network from file (hdf5)")
    parser.add_argument("networkFile", help="Network file (hdf5)", type=str)
    args = parser.parse_args()

    pn = PlotNetwork(args.networkFile)
    pn.plot()
