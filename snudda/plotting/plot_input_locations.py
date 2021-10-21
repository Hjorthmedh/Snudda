import os
import h5py
import numpy as np

from snudda.utils import SnuddaLoad
from snudda.neurons.neuron_morphology import NeuronMorphology
from snudda.neurons.neuron_prototype import NeuronPrototype
import matplotlib.pyplot as plt


class SnuddaPlotInputLocations:

    def __init__(self, network_path):

        self.network_path = network_path
        self.network_file = os.path.join(network_path, "network-synapses.hdf5")
        self.input_file = os.path.join(network_path, "input-spikes.hdf5")
        self.neuron_cache = dict()

        self.snudda_load = SnuddaLoad(self.network_file)

        if os.path.exists(self.input_file):
            self.input_data = h5py.File(self.input_file, "r")
        else:
            self.input_data = None

    def plot_neuron_inputs(self, neuron_id, input_type=None, colour=None):

        # TODO: Add ability to plot touch detected inputs also (use blue colour for them)

        coords = self.get_input_coords(neuron_id=neuron_id, input_type=input_type)

        if not colour:
            colour = "r"

        nm = self.load_neuron(neuron_id=neuron_id)
        ax = nm.plot_neuron(soma_colour=[0, 0, 0], dend_colour=[0, 0, 0], plot_axon=False, plot_dendrite=True)
        ax.scatter(xs=coords[:, 0], ys=coords[:, 1], zs=coords[:, 2], c=colour, marker=".")

        neuron_name = self.snudda_load.data["neurons"][neuron_id]["name"]

        if input_type:
            plt.title(f"{input_type} input to {neuron_name} ({neuron_id})")
            f_name = f"input-to-{neuron_id}-{neuron_name}-{input_type}.png"
        else:
            plt.title(f"Input to {neuron_name} ({neuron_id})")
            f_name = f"input-to-{neuron_id}-{neuron_name}.png"

        fig_name = os.path.join(self.network_path, "figures", f_name)
        plt.savefig(fig_name, dpi=300)

    def get_input_locations(self, neuron_id, input_type=None):

        if not self.input_data:
            # No input data available
            return None

        section_id = []
        section_x = []

        for input_name in self.input_data["input"][str(neuron_id)]:
            if input_type and input_name != input_type:
                continue

            section_id = section_id + list(self.input_data["input"][str(neuron_id)][input_name]["sectionID"])
            section_x = section_x + list(self.input_data["input"][str(neuron_id)][input_name]["sectionX"])

        return np.array(section_id), np.array(section_x)

    def get_input_coords(self, neuron_id, input_type=None):

        """ Returns input coordinates for all external input of input_type

            Args:
                neuron_id (int): Neuron ID
                input_type (str): Input type to show, eg. "Cortical" or None for all.

            Returns:
                coords: n x 3 dimensional array with x,y,z coordinates of external input synapses

        """

        section_id, section_x = self.get_input_locations(neuron_id=neuron_id, input_type=input_type)
        nm = self.load_neuron(neuron_id=neuron_id)
        coords = np.zeros((len(section_id), 3))

        for idx, (sec_id, sec_x) in enumerate(zip(section_id, section_x)):
            coords[idx, :] = nm.get_section_coordinates(section_id=sec_id, section_x=sec_x)

        return coords

    def load_neuron(self, neuron_id):

        if neuron_id not in self.neuron_cache:
            neuron_info = self.snudda_load.data["neurons"][neuron_id]
            nm = NeuronMorphology(name=neuron_info["name"],
                                  swc_filename=neuron_info["morphology"],
                                  position=neuron_info["position"],
                                  rotation=neuron_info["rotation"])
            self.neuron_cache[neuron_id] = nm

        return self.neuron_cache[neuron_id]


def plot_input_location_cli():

    import argparse
    parser = argparse.ArgumentParser("Plot input locations")
    parser.add_argument("networkPath", help="Path to network directory")
    parser.add_argument("neuronID", help="NeuronID to inspect", type=int)
    parser.add_argument("--inputType", help="Input type to show (default all)", default=None, type=str)

    args = parser.parse_args()

    pl = SnuddaPlotInputLocations(network_path=args.networkPath)
    pl.plot_neuron_inputs(neuron_id=args.neuronID, input_type=args.inputType)


if __name__ == "__main__":
    plot_input_location_cli()