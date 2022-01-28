import unittest
import os
import sys
import numpy as np

from snudda.place.create_cube_mesh import create_cube_mesh
from snudda.detect.detect import SnuddaDetect
from snudda.place.place import SnuddaPlace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


class TestDetectSynapseCluster(unittest.TestCase):

    def setUp(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        network_path = os.path.join(os.path.dirname(__file__), "networks", "network_testing_detect_cluster")

        create_cube_mesh(file_name=os.path.join(network_path, "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=500e-6)

        config_file = os.path.join(network_path, "network-config.json")
        position_file = os.path.join(network_path, "network-neuron-positions.hdf5")
        save_file = os.path.join(network_path, "voxels", "network-putative-synapses.hdf5")

        #  TODO: If d_view is None code run sin serial, add test parallel
        sp = SnuddaPlace(config_file=config_file, d_view=None, verbose=True)
        sp.place()

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 100 micrometer distance from centre of soma.

        self.sd = SnuddaDetect(config_file=config_file, position_file=position_file,
                               save_file=save_file, rc=None,
                               hyper_voxel_size=150, verbose=True)


    def test_clustering(self):

        neuron_positions = np.array([[0, 20, 0],    # Postsynaptiska
                                     [0, 80, 0],
                                     [0, 140, 0],
                                     [0, 200, 0],
                                     [20, 0, 0],    # Presynaptiska
                                     [80, 0, 0],
                                     [140, 0, 0],
                                     [200, 0, 0]
                                     ])*1e-6

        for idx, pos in enumerate(neuron_positions):
            self.sd.neurons[idx]["position"] = pos

        ang = -np.pi/2
        R_x = np.array([[1, 0, 0],
                        [0, np.cos(ang), -np.sin(ang)],
                        [0, np.sin(ang), np.cos(ang)]])

        ang = np.pi/2
        R_y = np.array([[np.cos(ang), 0, np.sin(ang)],
                        [0, 1, 0],
                        [-np.sin(ang), 0, np.cos(ang)]])

        for idx in range(0, 4):    # Post-synaptic neurons
            self.sd.neurons[idx]["rotation"] = R_x

        for idx in range(4, 8):   # Pre-synaptic neurons
            self.sd.neurons[idx]["rotation"] = R_y

        self.sd.detect(restart_detection_flag=True)

        synapse_voxel_loc = self.sd.hyper_voxel_synapses[:self.sd.hyper_voxel_synapse_ctr, 2:5]
        synapse_coords = synapse_voxel_loc * self.sd.voxel_size + self.sd.hyper_voxel_origo

        fig_path = os.path.join("networks", "network_testing_detect_cluster", "figures")

        if not os.path.exists(fig_path):
            os.mkdir(fig_path)

        if True:   # Set to True to include plot
            self.sd.plot_hyper_voxel(plot_neurons=True, fig_file_name="touch-detection-clusters-validation")

        # Lägg till tester som kollar positionerna på synapserna
        # import pdb
        # pdb.set_trace()

if __name__ == '__main__':
    unittest.main()
