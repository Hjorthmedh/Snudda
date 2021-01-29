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

        if not os.path.exists(network_path):
            os.mkdir(network_path)

        print(f"Current directory: {os.path.dirname(os.path.realpath(__file__))}")
        print(f"network_path: {network_path}")

        create_cube_mesh(file_name=os.path.join(network_path, "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=500e-6)

        config_file = os.path.join(network_path, "network-config.json")
        position_file = os.path.join(network_path, "network-neuron-positions.hdf5")
        save_file = os.path.join(network_path, "voxels", "network-putative-synapses.hdf5")

        #  TODO: If d_view is None code run sin serial, add test parallel
        sp = SnuddaPlace(config_file=config_file, d_view=None)

        sp.read_config()
        sp.write_data(position_file)

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 100 micrometer distance from centre of soma.

        self.sd = SnuddaDetect(config_file=config_file, position_file=position_file,
                               save_file=save_file, rc=None,
                               hyper_voxel_size=120)



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

        # Check location and source, targets.
        with self.subTest(stage="gap_junction_check"):
            self.assertEqual(self.sd.hyper_voxel_gap_junction_ctr, 1)
            self.assertTrue((self.sd.hyper_voxel_gap_junctions[0, 0:2] == [4, 20]).all())
            self.assertTrue(self.compare_voxels_to_coordinates(self.sd.hyper_voxel_gap_junctions[0, 6:9],
                                                               np.array([100, 100, 0])*1e-6))
        with self.subTest(stage="synapses_check"):
            self.assertEqual(self.sd.hyper_voxel_synapse_ctr, 101)

            for pre_id in range(0, 10):
                for post_id in range(0, 10):
                    self.assertEqual(self.check_neuron_pair_has_synapses(pre_id, post_id), 1)

        # We should probably store the matrix as unsigned.
        self.assertTrue((self.sd.hyper_voxel_synapses > 0).all())
        self.assertTrue((self.sd.hyper_voxel_gap_junctions > 0).all())

    def check_neuron_pair_has_synapse(self, pre_neuron, post_neuron):

        connections = dict()

        for synapse_row in self.sd.hyper_voxel_synapses:

            loc = (synapse_row[1], synapse_row[0])
            if loc in connections:
                connections[loc] += 1
            else:
                connections[loc] = 1

        if (pre_neuron, post_neuron) in connections:
            return connections[(pre_neuron, post_neuron)]
        else:
            return False

    def convert_to_coordinates(self, voxel_xyz):
        return np.array(voxel_xyz) * self.sd.voxel_size + self.sd.hyper_voxel_origo

    def compare_voxels_to_coordinates(self, voxel_index, coordinates):
        return (np.abs(self.convert_to_coordinates(voxel_index) - np.array(coordinates)) < self.sd.voxel_size).all()

if __name__ == '__main__':
    unittest.main()

# python3 -m unittest test_detect
    