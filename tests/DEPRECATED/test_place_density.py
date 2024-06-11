import unittest
import os
import sys
import json
import numpy as np

from snudda.place.region_mesh_redux import RegionMeshRedux
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

        cnc.define_striatum(num_dSPN=0, num_iSPN=0, num_FS=2000, num_LTS=0, num_ChIN=0,
                            mesh_file=mesh_file, mesh_bin_width=5e-4, neurons_dir=neuron_dir)

        create_slice_mesh(file_name=mesh_file,
                          centre_point=[1e-3, 1e-3, 1e-3],
                          x_len=2e-3,
                          y_len=2e-3,
                          z_len=2e-3,
                          description="Test slice")

        # Linear density = x coordinate, obs we give a relative density profile
        # (real density is scaled based on number of neurons)
        density_function = "abs(x)"
        cnc.add_neuron_density("Striatum", "FS", density_func=density_function)

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
        neuron_pos = sl.data["neuron_positions"]

        if not os.path.isdir("figures"):
            os.mkdir("figures")

        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(xs=neuron_pos[:, 0], ys=neuron_pos[:, 1], zs=neuron_pos[:, 2])
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        plt.savefig(os.path.join("figures", "density-scatter-plot.png"))

        for coord_str, coord in zip(['x', 'y', 'z'], range(0, 3)):
            plt.figure()
            plt.hist(neuron_pos[:, coord])
            plt.savefig(os.path.join("figures", f"density-hist-{coord_str}.png"))

        binned, bin_edges = np.histogram(neuron_pos[:, 0], 5)
        # Density should increase with x-coordinate
        self.assertTrue((np.diff(binned) > 0).all())
        self.assertTrue(binned[0]*3 < binned[-1])

        # import pdb
        # pdb.set_trace()

        # self.assertEqual(True, False)

        # TODO: Add a check for density -- currently it appears to not be working

    @staticmethod
    def is_inside(point):
        if (0 <= point).all() and (point <= 2e-3).all():
            return True
        else:
            return False

    def test_ray_casting(self):

        print("Testing ray casting algorithm")

        # This sets up a small cube, then tries the ray tracing function with different exterior points,
        # to see if they all give the same result

        mesh_file = os.path.join(self.network_path, "mesh-test.obj")

        create_slice_mesh(file_name=mesh_file,
                          centre_point=[1e-3, 1e-3, 1e-3],
                          x_len=2e-3,
                          y_len=2e-3,
                          z_len=2e-3,
                          description="Test slice")

        rm = RegionMeshRedux(mesh_file, use_cache=False, bin_width=5e-4)

        inside_point = rm.random_generator.uniform(low=0, high=2e-3, size=3)
        a_point = rm.random_generator.uniform(low=-1e-3, high=3e-3, size=3)

        # [0.00092975 0.00117448 0.00070073] -- check this point

        for test_ctr in range(0, 1000):
            a_point = rm.random_generator.uniform(low=-1e-3, high=3e-3, size=3)
            self.assertEqual(rm.ray_casting(a_point), self.is_inside(a_point))

            # assert rm.ray_casting(a_point) == self.is_inside(a_point), \
            # f"{a_point} ray_casting={rm.ray_casting(a_point)} should be {self.is_inside(a_point)} for point {a_point}"
            # print(f"Testing {a_point} ray_casting={rm.ray_casting(a_point)} ({self.is_inside(a_point)})")


if __name__ == '__main__':
    unittest.main()
