import os

import numpy as np
from snudda.create_cube_mesh import create_cube_mesh
from snudda.detect import SnuddaDetect
from snudda.place import SnuddaPlace
from snudda.plotting.plot_network import PlotNetwork
from snudda.prune import SnuddaPrune
from snudda.utils.reposition_neurons import RepositionNeurons


class TouchDetectionHypervoxelIllustration(object):

    def __init__(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        self.network_path = "touch_detection_hypervoxel_illustration_network"
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
                               hyper_voxel_size=60)

        neuron_positions = np.array([[10, 30, 70],  # Postsynaptiska
                                     [50, 60, 70],  # Presynaptiska
                                     ]) * 1e-6

        for idx, pos in enumerate(neuron_positions):
            self.sd.neurons[idx]["position"] = pos

        ang = -np.pi / 2
        R_x = np.array([[1, 0, 0],
                        [0, np.cos(ang), -np.sin(ang)],
                        [0, np.sin(ang), np.cos(ang)]])

        ang = np.pi * 0.2
        R_y = np.array([[np.cos(ang), 0, np.sin(ang)],
                        [0, 1, 0],
                        [-np.sin(ang), 0, np.cos(ang)]])


        ang = np.pi * (0.5 + 0.2)
        R_z0 = np.array([[np.cos(ang), -np.sin(ang), 0],
                         [np.sin(ang), np.cos(ang), 0],
                         [0, 0, 1]])

        ang = np.pi*0.4
        R_z1 = np.array([[np.cos(ang), -np.sin(ang), 0],
                         [np.sin(ang), np.cos(ang), 0],
                         [0, 0, 1]])

        # Post synaptic
        self.sd.neurons[0]["rotation"] = R_z0

        # Presynaptic neuron
        self.sd.neurons[1]["rotation"] = np.matmul(R_z1, R_y)

        self.sd.detect(restart_detection_flag=True)

        self.sd.process_hyper_voxel(0)
        plt, ax = self.sd.plot_hyper_voxel(plot_neurons=False,
                                           draw_axon_voxels=True,
                                           draw_dendrite_voxels=True,
                                           elev_azim=(50, -22),
                                           title="",
                                           fig_file_name="touch_detection_illustration-voxels.pdf",
                                           dpi=300)

        plt, ax = self.sd.plot_hyper_voxel(plot_neurons=True,
                                           draw_axon_voxels=False,
                                           draw_dendrite_voxels=False,
                                           elev_azim=(50, -22),
                                           title="",
                                           fig_file_name="touch_detection_illustration-morph.pdf",
                                           dpi=300)

        print(f"\n--> Figures written to {self.sd.network_path}/figures")



if __name__ == "__main__":

    tdhi = TouchDetectionHypervoxelIllustration()
