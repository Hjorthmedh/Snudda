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

        self.sd = SnuddaDetect(config_file=config_file, position_file=position_file, save_file=save_file, rc=None)

        # Reposition the neurons for the


    def test_detect(self):

        neuron_positions = np.array([[0, 20, 0],
                                     [0, 40, 0],
                                     [0, 60, 0],
                                     [0, 80, 0],
                                     [-20, 100, 0],  # This one is intentionally pulled back
                                     [20, 0, 0],
                                     [40, 0, 0],
                                     [60, 0, 0],
                                     [80, 0, 0],
                                     [100, 0, 0]])*1e-6

        for idx, pos in enumerate(neuron_positions):
            self.sd.neurons[idx]["position"] = pos

        ang = np.pi/2
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(ang), -np.sin(ang)],
                       [0, np.sin(ang), np.cos(ang)]])

        for idx in range(0, 5):
            self.sd.neurons[idx]["rotation"] = Rx

        self.sd.detect(restart_detection_flag=True)
        self.assertEqual(self.sd.hyper_voxel_synapse_ctr, 24)  # 5x5 - 1 (one pulled back slightly)

        synapse_voxel_loc = self.sd.hyper_voxel_synapses[:self.sd.hyper_voxel_synapse_ctr, 2:5]
        synapse_coords = synapse_voxel_loc * self.sd.voxel_size + self.sd.hyper_voxel_origo

        fig_path = os.path.join("tests", "network_testing_detect", "figures")

        if not os.path.exists(fig_path):
            os.mkdir(fig_path)

        self.sd.plot_hyper_voxel(plot_neurons=True)

        # TODO: Also add tests for gap junctions and for soma-axon synapses


        import pdb
        pdb.set_trace()


if __name__ == '__main__':
    unittest.main()

# python3 -m unittest test_detect
    