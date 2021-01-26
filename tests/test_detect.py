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
                         side_len=100)

        position_file = os.path.join(network_path, "network-neuron-positions.hdf5")

        #  TODO: If d_view is None code run sin serial, add test parallel
        sp = SnuddaPlace(config_file=self.config_file, d_view=None)  

        sp.read_config()
        sp.write_data_HDF5(position_file)

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 100 micrometer distance from centre of soma.

        config_file = os.path.join(network_path, "network-config.json")
        position_file = os.path.join(network_path, "network-neuron-positions.hdf5")
        save_file = os.path.join(network_path, "voxels", "network-putative-synapses.hdf5")

        self.sd = SnuddaDetect(config_file=config_file, position_file=position_file, save_file=save_file, rc=None)

        # Reposition the neurons for the


    def test_detect(self):

        import pdb
        pdb.set_trace()

        self.sd.detect(restart_detection_flag=True)

if __name__ == '__main__':
    unittest.main()

# python3 -m unittest test_detect
    