import numpy as np
from snudda.load import SnuddaLoad
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from snudda.neuron_morphology import NeuronMorphology


class PlotNetwork(object):

    def __init__(self, network_file):

        self.network_file = network_file
        self.sl = SnuddaLoad(self.network_file)
        self.prototype_neurons = dict()

    def plot(self, plot_axon=True, plot_dendrite=True, plot_synapses=True):

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        if "simulationOrigo" in self.sl.data:
            simulation_origo = self.sl.data["simulationOrigo"]
        else:
            simulation_origo = np.array([0, 0, 0])

        # Plot neurons
        for neuron_info in self.sl.data["neurons"]:
            neuron = self.load_neuron(neuron_info)
            neuron.plot_neuron(axis=ax,
                               plot_axon=plot_axon,
                               plot_dendrite=plot_dendrite,
                               soma_colour=(0, 0, 0),
                               axon_colour=(1, 0, 0),
                               dend_colour=(0, 0, 0))

        if plot_synapses and "synapseCoords" in self.sl.data:
            ax.scatter(self.sl.data["synapseCoords"][:, 0],
                       self.sl.data["synapseCoords"][:, 1],
                       self.sl.data["synapseCoords"][:, 2], c="royalblue")

        # Plot synapses

        plt.ion()
        plt.show()

        import pdb
        pdb.set_trace()

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