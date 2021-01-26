import unittest
import os
import sys
import json
import numpy as np

from snudda.create_cube_mesh import create_cube_mesh
from snudda.place import SnuddaPlace

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


class TestDetect(unittest.TestCase):

    def setUp(self):

        create_cube_mesh(file_name=os.path.join("tests", "network_testing_detect", "mesh", "simple_mesh.obj"),
                         centre_point=(0, 0, 0),
                         side_len=100)

        position_file = os.path.join(self.sim_name, "network-neuron-positions.hdf5")

        sp = SnuddaPlace(config_file=self.config_file,
                         log_file=None,
                         verbose=True,
                         d_view=None,  # TODO: If d_view is None code run sin serial, add test parallel
                         h5libver="latest")

        sp.read_config()
        sp.write_data_HDF5(position_file)

        # We want to load in the ball and stick neuron that has 20 micrometer soma diameter, and axon (along y-axis),
        # and dendrite along (x-axis) out to 100 micrometer distance from centre of soma.

        config_file = self.network_path + "/network-config.json"
        position_file = self.network_path + "/network-neuron-positions.hdf5"
        log_filename = self.network_path + "/log/logFile-touch-detection.txt"
        save_file = self.network_path + "/voxels/network-putative-synapses.hdf5"


        pass

    