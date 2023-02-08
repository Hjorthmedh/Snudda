import os
import numpy as np

from snudda.neurons.morphology_data import MorphologyData
from snudda.plotting.plot_spike_raster_v2 import SnuddaPlotSpikeRaster2
from snudda.utils.snudda_path import get_snudda_data
from snudda.neurons.neuron_prototype import NeuronPrototype
from snudda.utils.load import SnuddaLoad
import matplotlib.pyplot as plt


class PlotNetwork(object):

    def __init__(self, network, snudda_data=None):

        if os.path.isdir(network):
            network_file = os.path.join(network, "network-synapses.hdf5")
        else:
            network_file = network

        self.network_file = network_file
        self.network_path = os.path.dirname(self.network_file)

        self.sl = SnuddaLoad(self.network_file)
        self.snudda_data = get_snudda_data(snudda_data=snudda_data,
                                           network_path=self.network_path)
        self.prototype_neurons = dict()
        self.extra_axon_cache = dict()

    def close(self):
        self.sl.close()

    def plot(self, plot_axon=True, plot_dendrite=True, plot_synapses=True,
             neuron_id_list=None, filter_synapses_pre_id_list=None,
             title=None, title_pad=None, show_axis=True,
             elev_azim=None, fig_name=None, dpi=600,
             colour_population_unit=False,
             colour_neuron_type=True):

        if type(plot_axon) == bool:
            plot_axon = np.ones((self.sl.data["nNeurons"],), dtype=bool) * plot_axon

        if type(plot_dendrite) == bool:
            plot_dendrite = np.ones((self.sl.data["nNeurons"],), dtype=bool) * plot_dendrite

        assert len(plot_axon) == len(plot_dendrite) == len(self.sl.data["neurons"])

        fig = plt.figure(figsize=(6, 6.5))
        ax = plt.axes(projection='3d')

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

        if colour_neuron_type:
            colour_lookup_helper = dict()

            for nid in self.sl.get_neuron_id():
                colour_lookup_helper[nid] = SnuddaPlotSpikeRaster2.get_colours(self.sl.data["neurons"][nid]["type"])

            colour_lookup = lambda x: colour_lookup_helper[x]

        # Plot neurons
        for neuron_info, pa, pd in zip(self.sl.data["neurons"], plot_axon, plot_dendrite):

            if neuron_id_list:
                # We only plot certain neuronID
                if neuron_info["neuronID"] not in neuron_id_list:
                    continue

            soma_colour = colour_lookup(neuron_info["neuronID"])
            neuron = self.load_neuron(neuron_info)

            neuron.plot_neuron(axis=ax,
                               plot_axon=pa,
                               plot_dendrite=pd,
                               soma_colour=soma_colour,
                               axon_colour=None,  #"darksalmon",  #"maroon",
                               dend_colour=None)  #"silver")   # Can also write colours as (0, 0, 0) -- rgb

        # Plot synapses
        if plot_synapses and "synapseCoords" in self.sl.data and self.sl.data["synapseCoords"].size > 0:

            if neuron_id_list:
                post_id = self.sl.data["synapses"][:, 1]

                keep_flag = np.zeros(post_id.shape, dtype=bool)
                for nid in neuron_id_list:
                    keep_idx = np.where(post_id == nid)[0]
                    keep_flag[keep_idx] = True

                if filter_synapses_pre_id_list:
                    pre_id = self.sl.data["synapses"][:, 0]

                    pre_keep_flag = np.zeros(post_id.shape, dtype=bool)
                    for nid in filter_synapses_pre_id_list:
                        keep_idx = np.where(pre_id == nid)[0]
                        pre_keep_flag[keep_idx] = True

                    keep_flag = np.logical_and(keep_flag, pre_keep_flag)

                keep_idx = np.where(keep_flag)[0]

                x = self.sl.data["synapseCoords"][:, 0][keep_idx]
                y = self.sl.data["synapseCoords"][:, 1][keep_idx]
                z = self.sl.data["synapseCoords"][:, 2][keep_idx]

                plt.figtext(0.5, 0.20, f"{keep_idx.size} synapses", ha="center", fontsize=18)

            else:
                x = self.sl.data["synapseCoords"][:, 0]
                y = self.sl.data["synapseCoords"][:, 1]
                z = self.sl.data["synapseCoords"][:, 2]

            ax.scatter(x, y, z, color=(1, 0, 0))

            if neuron_id_list is None:
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
            print(f"Writing figure: {fig_path}")

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
            self.prototype_neurons[neuron_name] = NeuronPrototype(neuron_name=neuron_info["name"],
                                                                  neuron_path=neuron_info["neuronPath"],
                                                                  snudda_data=self.snudda_data,
                                                                  load_morphology=True,
                                                                  virtual_neuron=False)

        neuron = self.prototype_neurons[neuron_name].clone(position=neuron_info["position"],
                                                           rotation=neuron_info["rotation"],
                                                           parameter_key=neuron_info["parameterKey"],
                                                           morphology_key=neuron_info["morphologyKey"],
                                                           modulation_key=neuron_info["modulationKey"])

        if "extraAxons" in neuron_info:
            for axon_name, axon_info in neuron_info["extraAxons"].items():

                if axon_info["morphology"] not in self.extra_axon_cache:
                    self.extra_axon_cache[axon_info["morphology"]] = MorphologyData(swc_file=axon_info["morphology"],
                                                                                    parent_tree_info=None,
                                                                                    snudda_data=self.snudda_data)

                neuron.add_morphology(swc_file=axon_info["morphology"],
                                      name=axon_name,
                                      position=axon_info["position"],
                                      rotation=axon_info["rotation"],
                                      morphology_data=self.extra_axon_cache[axon_info["morphology"]])

        return neuron

    def plot_populations(self, unmarked_alpha=0.3):

        fig = plt.figure(figsize=(6, 6.5))
        ax = plt.axes(projection='3d')

        assert "populationUnit" in self.sl.data

        population_unit = self.sl.data["populationUnit"]
        pop_units = sorted(list(set(population_unit)))
        cmap = plt.get_cmap('tab20', len(pop_units))
        neuron_colours = []

        for idx, pu in enumerate(population_unit):
            if pu > 0:
                neuron_colours.append(list(cmap(pu)))
            else:
                neuron_colours.append([0.6, 0.6, 0.6, 1.0])

        neuron_colours = np.array(neuron_colours)
        positions = self.sl.data["neuronPositions"]

        pop_unit_idx = np.where(population_unit > 0)[0]
        unmarked_idx = np.where(population_unit == 0)[0]

        ax.scatter(positions[pop_unit_idx, 0], positions[pop_unit_idx, 1], positions[pop_unit_idx, 2],
                   c=neuron_colours[pop_unit_idx, :], marker='o', s=20, alpha=1)

        ax.scatter(positions[unmarked_idx, 0], positions[unmarked_idx, 1], positions[unmarked_idx, 2],
                   c=neuron_colours[unmarked_idx, :], marker='o', s=20, alpha=unmarked_alpha)

        # ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c=neuron_colours, marker='o', s=20,alpha=alphas)

        self.equal_axis(ax)

        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")

        for pop in pop_units:
            print(f"Population unit {pop} has {len(np.where(population_unit == pop)[0])} neurons")

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


def snudda_plot_network_cli():

    from argparse import ArgumentParser
    parser = ArgumentParser(description="Plot snudda network from file (hdf5)")
    parser.add_argument("networkFile", help="Network file (hdf5)", type=str)
    parser.add_argument("--neuronID", help="List of Neuron ID to show", type=int, default=None, nargs="+",
                        dest="neuron_id_list")
    parser.add_argument("--showAxons", help="Show Axons of neurons", action="store_true")
    parser.add_argument("--showDendrites", help="Show dendrites of neurons", action="store_true")
    parser.add_argument("--showSynapses", help="Show synapses of neurons", action="store_true")
    parser.add_argument("--preID", help="Only show synapses from pre ID neurons",
                        dest="pre_id_list",
                        nargs="+", type=int, default=None)
    parser.add_argument("--wait", action="store_true")
    args = parser.parse_args()

    pn = PlotNetwork(args.networkFile)
    pn.plot(fig_name="network-plot.png", neuron_id_list=args.neuron_id_list,
            plot_axon=args.showAxons, plot_dendrite=args.showDendrites, plot_synapses=args.showSynapses,
            filter_synapses_pre_id_list=args.pre_id_list)

    if args.wait:
        input("Press any key to exit")


if __name__ == "__main__":
    snudda_plot_network_cli()
