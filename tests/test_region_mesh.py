import unittest
import numpy as np
from snudda.place.create_cube_mesh import create_cube_mesh
from snudda.place.region_mesh_vedo import RegionMesh


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:

        mesh_path = "my_test_mesh.obj"
        centre = (100e-6, 110e-6, 120e-6)
        side_len = 100e-6
        create_cube_mesh(mesh_path, centre, side_len)
        self.rm = RegionMesh(mesh_path=mesh_path)

    def test_inside(self):

        rng = np.random.default_rng(1)
        test_points = rng.uniform(low=-100e-6, high=300e-6, size=(1000, 3))

        inside_flag = self.rm.check_inside(test_points)

        ref_inside_flag = (50e-6 <= test_points[:, 0]) & (test_points[:, 0] <= 150e-6) \
            & (60e-6 <= test_points[:, 1]) & (test_points[:, 1] <= 160e-6) \
            & (70e-6 <= test_points[:, 2]) & (test_points[:, 2] <= 170e-6)

        self.assertTrue((inside_flag == ref_inside_flag).all())

    def test_distance(self):

        test_points = np.array([[40e-6, 110e-6, 120e-6],  # distance 10e-6
                                [100e-6, 90e-6, 120e-6],  # distane -20e-6
                                [160e-6, 170e-6, 180e-6],  # distance np.sqrt(3*100)
                                ])

        ref_distance = np.array([10e-6, -30e-6, np.sqrt(3*100e-12)])

        dist = self.rm.distance_to_border(test_points)

        self.assertTrue(np.allclose(dist, ref_distance, rtol=1e-7))


if __name__ == '__main__':
    unittest.main()
