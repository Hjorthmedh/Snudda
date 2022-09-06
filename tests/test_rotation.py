import os
import unittest
import numpy as np
from snudda.place.rotation import SnuddaRotate


class TestRotation(unittest.TestCase):

    def setUp(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        self.sr = SnuddaRotate()
        self.sr.parse_config_file("data/RotationTest.json")

    def test_rotations(self):

        rng = np.random.default_rng()

        rot = self.sr.get_rotations(volume_name="Striatum", neuron_type="LTS",
                                    neuron_positions=[[0, 0, 0]],
                                    rng=rng)
        self.assertEqual(np.linalg.det(rot[0]), 1)

        vec_test = np.array([0, 0, 1])
        vec_rot = rot[0].dot(vec_test)
        vec_facit = np.array([-1, -1, -1])
        self.assertTrue(np.allclose(vec_rot / np.linalg.norm(vec_rot), vec_facit / np.linalg.norm(vec_facit)))

        rot2 = self.sr.get_rotations(volume_name="Striatum", neuron_type="LTS",
                                     neuron_positions=[[0, 0, 0.5e-3]],
                                     rng=rng)
        self.assertTrue(np.abs(np.linalg.det(rot2[0]) - 1) < 1e-6)

        vec_rot2 = rot2[0].dot(vec_test)
        vec_facit2 = np.array([-1, -1, 0])

        self.assertTrue(np.allclose(vec_rot2 / np.linalg.norm(vec_rot2),
                                    vec_facit2 / np.linalg.norm(vec_facit2),
                            atol=1e-6))

        # TODO: Also need to test the random rotations


if __name__ == '__main__':
    unittest.main()
