import unittest
import os
import sys
import json
import numpy as np

from snudda.create_cube_mesh import create_cube_mesh
from snudda.detect import SnuddaDetect
from snudda.place import SnuddaPlace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


class TestDetect(unittest.TestCase):

    def setUp(self):

        network_path = os.path.join("tests", "network_testing_detect")

        create_cube_mesh(file_name=os.path.join(network_path, "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=500e-6)

        config_file = os.path.join(network_path, "network-config.json")
        position_file = os.path.join(network_path, "network-neuron-positions.hdf5")
        save_file = os.path.join(network_path, "voxels", "network-putative-synapses.hdf5")

        #  TODO: If d_view is None code run sin serial, add test parallel
        sp = SnuddaPlace(config_file=config_file, d_view=None)

        sp.read_config()
        sp.write_data_HDF5(position_file)

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 100 micrometer distance from centre of soma.

        self.sd = SnuddaDetect(config_file=config_file, position_file=position_file,
                               save_file=save_file, rc=None,
                               hyper_voxel_size=120)

        # Reposition the neurons for the


    def test_detect(self):

        neuron_positions = np.array([[0, 20, 0],    # Postsynaptiska
                                     [0, 40, 0],
                                     [0, 60, 0],
                                     [0, 80, 0],
                                     [0, 100, 0],
                                     [0, 120, 0],
                                     [0, 140, 0],
                                     [0, 160, 0],
                                     [0, 180, 0],
                                     [0, 200, 0],
                                     [20, 0, 0],    # Presynaptiska
                                     [40, 0, 0],
                                     [60, 0, 0],
                                     [80, 0, 0],
                                     [100, 0, 0],
                                     [120, 0, 0],
                                     [140, 0, 0],
                                     [160, 0, 0],
                                     [180, 0, 0],
                                     [200, 0, 0],
                                     [100, 100, -40],  # To get a gap junction
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

        ang = -np.pi/2
        R_gj = np.array([[np.cos(ang), 0, np.sin(ang)],
                         [0, 1, 0],
                         [-np.sin(ang), 0, np.cos(ang)]])

        for idx in range(0, 10):    # Post synaptic neurons
            self.sd.neurons[idx]["rotation"] = R_x

        for idx in range(10, 20):   # Presynaptic neurons
            self.sd.neurons[idx]["rotation"] = R_y

        self.sd.neurons[20]["rotation"] = R_gj

        self.sd.detect(restart_detection_flag=True)

        synapse_voxel_loc = self.sd.hyper_voxel_synapses[:self.sd.hyper_voxel_synapse_ctr, 2:5]
        synapse_coords = synapse_voxel_loc * self.sd.voxel_size + self.sd.hyper_voxel_origo

        fig_path = os.path.join("tests", "network_testing_detect", "figures")

        if not os.path.exists(fig_path):
            os.mkdir(fig_path)

        self.sd.plot_hyper_voxel(plot_neurons=True)

        # TODO: Also add tests for soma-axon synapses
        self.assertEqual(self.sd.hyper_voxel_synapse_ctr, 101)
        self.assertEqual(self.sd.hyper_voxel_gap_junction_ctr, 1)


#        import pdb
#        pdb.set_trace()


if __name__ == '__main__':
    unittest.main()

# python3 -m unittest test_detect
    