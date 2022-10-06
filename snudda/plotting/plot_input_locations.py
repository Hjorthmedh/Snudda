import os
import h5py
import numpy as np

from snudda.utils.snudda_path import get_snudda_data
from snudda.utils import SnuddaLoad
from snudda.neurons.neuron_morphology import NeuronMorphology
from snudda.neurons.neuron_prototype import NeuronPrototype
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class SnuddaPlotInputLocations:

    def __init__(self, network_path=None, network_file=None, input_file=None):

        if not network_path:
            if network_file:
                network_path = os.path.dirname(network_file)
            elif input_file:
                network_path = os.path.dirname(input_file)
            else:
                network_path = "./"

        self.network_path = network_path

        self.snudda_data = get_snudda_data(network_path=self.network_path)

        if network_file is None:
            self.network_file = os.path.join(network_path, "network-synapses.hdf5")
        else:
            self.network_file = network_file

        if input_file is None:
            self.input_file = os.path.join(network_path, "input-spikes.hdf5")
        else:
            self.input_file = os.path.join(network_path, input_file)
            
        self.neuron_cache = dict()

        self.snudda_load = SnuddaLoad(self.network_file)

        if os.path.exists(self.input_file):
            self.input_data = h5py.File(self.input_file, "r")
        else:
            self.input_data = None

    def plot_neuron_inputs(self, neuron_id,
                           input_type=None,
                           show_internal_synapses=True,
                           external_colour=None, internal_colour=None,
                           ax=None, neuron_colour=None,
                           size=10,
                           save_fig=True,
                           dpi=300,
                           show_figure=True):

        coords = self.get_input_coords(neuron_id=neuron_id, input_type=input_type)

        if external_colour is None:
            external_colour = "r"

        if internal_colour is None:
            internal_colour = "b"

        nm = self.load_neuron(neuron_id=neuron_id)

        if ax is None:
            fig = plt.figure(visible=show_figure)
            ax = fig.add_subplot(111, projection='3d')

        if neuron_colour is None:
            neuron_colour = np.array([0, 0, 0])

        ax = nm.plot_neuron(soma_colour=neuron_colour, dend_colour=neuron_colour,
                            plot_axon=False, plot_dendrite=True,
                            show_plot=False, axis=ax)

        if len(coords) > 0:
            print(f"Plotting {len(coords)} external synapses")
           
            ax.scatter(xs=coords[:, 0], ys=coords[:, 1], zs=coords[:, 2],
                       c=external_colour, marker=".", s=size)

        if show_internal_synapses:
            syn_coords = self.get_synapse_coords(neuron_id=neuron_id)

            if syn_coords is not None:
                ax.scatter(xs=syn_coords[:, 0], ys=syn_coords[:, 1], zs=syn_coords[:, 2],
                           c=internal_colour, marker=".", s=size)

        neuron_name = self.snudda_load.data["neurons"][neuron_id]["name"]

        if show_internal_synapses:
            syn_txt = "-and-internal-synapses"
            syn_title = "and internal "
        else:
            syn_txt = ""
            syn_title = ""

        if input_type:
            plt.title(f"{input_type} {syn_title}input to {neuron_name} ({neuron_id})")
            f_name = f"input-to-{neuron_id}-{neuron_name}-{input_type}{syn_txt}.png"
        else:
            plt.title(f"Input to {neuron_name} ({neuron_id})")
            f_name = f"input-to-{neuron_id}-{neuron_name}{syn_txt}.png"

        fig_path = os.path.join(self.network_path, "figures")
        if not os.path.exists(fig_path):
            os.mkdir(fig_path)

        if save_fig:
            fig_name = os.path.join(fig_path, f_name)
            plt.savefig(fig_name, dpi=dpi)
            print(f"Figure written: {fig_name}")

        if show_figure:
            fig = ax.get_figure()
            fig.set_visible(show_figure)

        return ax

    def get_input_locations(self, neuron_id, input_type=None):

        if not self.input_data:
            # No input data available
            return None, None

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

        if section_id is None:
            return np.zeros((0, 3))

        nm = self.load_neuron(neuron_id=neuron_id)
        coords = np.zeros((len(section_id), 3))

        for idx, (sec_id, sec_x) in enumerate(zip(section_id, section_x)):
            coords[idx, :] = nm.get_section_coordinates(section_id=sec_id, section_x=sec_x)

        return coords

    def load_neuron(self, neuron_id):

        if neuron_id not in self.neuron_cache:
            neuron_info = self.snudda_load.data["neurons"][neuron_id]

            prot = NeuronPrototype(neuron_path=neuron_info["neuronPath"], neuron_name=neuron_info["name"],
                                 snudda_data=self.snudda_data)
            nm = prot.clone(parameter_key=neuron_info["parameterKey"], morphology_key=neuron_info["morphologyKey"],
                            position=neuron_info["position"], rotation=neuron_info["rotation"])
            self.neuron_cache[neuron_id] = nm

        return self.neuron_cache[neuron_id]

    def get_synapse_coords(self, neuron_id, pre_type=None):

        synapses, synapse_coords = self.snudda_load.find_synapses(post_id=neuron_id)

        if synapses is None:
            return None

        assert (neuron_id == synapses[:, 1]).all(), \
            f"Internal error, post_id should be {neuron_id} for all"

        # Also get synapse coords from section id and section x
        section_id = synapses[:, 9]
        section_x = synapses[:, 10] * 1e-3

        nm = self.load_neuron(neuron_id=neuron_id)
        coords = np.zeros((len(section_id), 3))

        for idx, (sec_id, sec_x) in enumerate(zip(section_id, section_x)):
            coords[idx, :] = nm.get_section_coordinates(section_id=sec_id, section_x=sec_x)

        dist = np.linalg.norm(coords - synapse_coords, axis=-1)
        max_dist = 10e-6  # np.sqrt(3*(5e-6 ** 2))
        if (dist > max_dist).any():
            print(f"Synapse coordinates mismatch {synapse_coords[np.where(dist > max_dist)[0], :]} "
                  f"vs {coords[np.where(dist > max_dist)[0], :]}"
                  f" (distances {dist[np.where(dist > max_dist)[0]]} with allowed max_dist = {max_dist})")

        pre_id = synapses[:, 0]

        if pre_type:
            pre_id_list = self.snudda_load.get_neuron_id_of_type(pre_type=pre_type)

            keep_idx = []
            for idx, pid in enumerate(pre_id):
                if pid in pre_id_list:
                    keep_idx.append(idx)

            coords = coords[idx, :]

        return coords


def plot_input_location_cli():

    import argparse
    parser = argparse.ArgumentParser("plot_input_locations")
    parser.add_argument("networkPath", help="Path to network directory")
    parser.add_argument("neuronID", help="NeuronID to inspect", type=int)
    parser.add_argument("--inputType", help="Input type to show (default all)", default=None, type=str)
    parser.add_argument("--showSynapses", help="Show internal synapses", action="store_true")

    args = parser.parse_args()

    pl = SnuddaPlotInputLocations(network_path=args.networkPath)
    pl.plot_neuron_inputs(neuron_id=args.neuronID, input_type=args.inputType,
                          show_internal_synapses=args.showSynapses)


if __name__ == "__main__":
    plot_input_location_cli()
