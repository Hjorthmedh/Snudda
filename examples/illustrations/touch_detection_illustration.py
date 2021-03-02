import os

import numpy as np
from snudda.create_cube_mesh import create_cube_mesh
from snudda.detect import SnuddaDetect
from snudda.place import SnuddaPlace
from snudda.plotting.plot_network import PlotNetwork
from snudda.prune import SnuddaPrune
from snudda.utils.reposition_neurons import RepositionNeurons


class TouchDetectionIllustration(object):

    def __init__(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        self.network_path = "touch_detection_illustration_network"
        self.config_file = os.path.join(self.network_path, "network-config.json")
        self.position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        self.save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        create_cube_mesh(file_name=os.path.join(self.network_path, "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=500e-6)

        sp = SnuddaPlace(config_file=self.config_file, d_view=None)
        sp.parse_config()
        sp.write_data(self.position_file)

        self.sd = SnuddaDetect(config_file=self.config_file, position_file=self.position_file,
                               save_file=self.save_file, rc=None,
                               hyper_voxel_size=150)

        # Reposition the neurons so we know how many synapses and where they will be located before pruning
        neuron_positions = np.array([[0, 59, 0],  # Postsynaptiska
                                     [0, 89, 0],
                                     [0, 119, 0],
                                     [0, 149, 0],
                                     [0, 179, 0],
                                     [0, 209, 0],
                                     [0, 239, 0],
                                     [0, 269, 0],
                                     [0, 299, 0],
                                     [0, 329, 0],
                                     [59, 0, 0],  # Presynaptiska
                                     [89, 0, 0],
                                     [119, 0, 0],
                                     [149, 0, 0],
                                     [179, 0, 0],
                                     [209, 0, 0],
                                     [239, 0, 0],
                                     [269, 0, 0],
                                     [299, 0, 0],
                                     [329, 0, 0],
                                     ]) * 1e-6

        for idx, pos in enumerate(neuron_positions):
            self.sd.neurons[idx]["position"] = pos

        ang = -np.pi / 2
        R_x = np.array([[1, 0, 0],
                        [0, np.cos(ang), -np.sin(ang)],
                        [0, np.sin(ang), np.cos(ang)]])

        ang = np.pi / 2
        R_y = np.array([[np.cos(ang), 0, np.sin(ang)],
                        [0, 1, 0],
                        [-np.sin(ang), 0, np.cos(ang)]])

        for idx in range(0, 10):  # Post synaptic neurons
            self.sd.neurons[idx]["rotation"] = R_x

        for idx in range(10, 20):  # Presynaptic neurons
            self.sd.neurons[idx]["rotation"] = R_y

        self.sd.detect(restart_detection_flag=True)

        # Also update so that the new positions are saved in the place file
        rn = RepositionNeurons(self.position_file)
        for neuron_info in self.sd.neurons:
            rn.place(neuron_info["neuronID"], position=neuron_info["position"], rotation=neuron_info["rotation"])
        rn.close()

        sp = SnuddaPrune(network_path=self.network_path)  # Use default config file
        sp.prune()
        sp = []

    def plot(self, fig_name):

        plot_axon = True
        plot_dendrite = True

        pn = PlotNetwork(self.network_path)
        plt, ax = pn.plot(fig_name=fig_name, show_axis=False,
                          plot_axon=plot_axon, plot_dendrite=plot_dendrite,
                          title="", title_pad=-14,
                          elev_azim=(-71, -35))

if __name__ == "__main__":

    tdi = TouchDetectionIllustration()
    tdi.plot("touch_detection_synapses.pdf")

    pass