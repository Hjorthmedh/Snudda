import unittest
import os
import sys
import json
import numpy as np

from snudda.utils import SnuddaLoad

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from snudda.init.init import SnuddaInit
from snudda.place.place import SnuddaPlace
from snudda.place import create_slice_mesh


class SnuddaDensityTest(unittest.TestCase):

    def setUp(self):

        # We want to setup a volume with density variations
        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))
        print(f"Current directory (detect): {os.path.dirname(os.path.realpath(__file__))}")

        neuron_dir = os.path.join(os.path.dirname(__file__), "validation")
        self.network_path = os.path.join("networks", "network_density")
        self.config_file = os.path.join(self.network_path, "network-config.json")
        cnc = SnuddaInit(struct_def={}, config_file=self.config_file, random_seed=1234)

        mesh_file = os.path.join(self.network_path, "mesh", "slice.obj")

        cnc.define_striatum(num_dSPN=0, num_iSPN=0, num_FS=1000, num_LTS=0, num_ChIN=0,
                            mesh_file=mesh_file, mesh_bin_width=1e-4, neurons_dir=neuron_dir)

        create_slice_mesh(file_name=mesh_file,
                          centre_point=[0.5e-3, 0.5e-3, 0.5e-3],
                          x_len=2e-3,
                          y_len=2e-3,
                          z_len=2e-3,
                          description="Test slice")

        # Linear density = x coordinate, obs we give a relative density profile
        # (real density is scaled based on number of neurons)
        density_function = "x"
        cnc.add_neuron_density("Striatum", "FSN", density_func=density_function)

        cnc.write_json(self.config_file)

        self.position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")

        npn = SnuddaPlace(config_file=self.config_file,
                          log_file=None,
                          verbose=True,
                          d_view=None,
                          h5libver="latest")

        npn.parse_config()
        npn.write_data(self.position_file)

    def test_density(self):

        # Load the positions
        sl = SnuddaLoad(self.position_file)
        neuron_pos = sl.data["neuronPositions"]

        if not os.path.isdir("figures"):
            os.mkdir("figures")

        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(xs=neuron_pos[:, 0], ys=neuron_pos[:, 1], zs=neuron_pos[:, 2])
        plt.savefig(os.path.join("figures", "density-scatter-plot.png"))

        for coord_str, coord in zip(['x', 'y', 'z'], range(0, 3)):
            plt.figure()
            plt.hist(neuron_pos[:, coord])
            plt.savefig(os.path.join("figures", f"density-hist-{coord_str}.png"))

        # import pdb
        # pdb.set_trace()

        # self.assertEqual(True, False)

        # TODO: Add a check for density -- currently it appears to not be working
        

if __name__ == '__main__':
    unittest.main()
