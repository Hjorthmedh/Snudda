import unittest
import numpy as np
from snudda.rotation import SnuddaRotate


class TestRotation(unittest.TestCase):

    def setUp(self):

        self.sr = SnuddaRotate()
        self.sr.parse_config_file("data/RotationTest.json")

    def test_rotations(self):

        rot = self.sr.get_rotation(volume_name="Striatum", neuron_type="LTS", neuron_position=[0, 0, 0])
        self.assertEqual(np.linalg.det(rot), 1)

        vec_test = np.array([0, 0, 1])
        vec_rot = rot.dot(vec_test)
        vec_facit = np.array([-1, -1, -1])
        self.assertTrue(np.allclose(vec_rot / np.linalg.norm(vec_rot), vec_facit / np.linalg.norm(vec_facit)))

        rot2 = self.sr.get_rotation(volume_name="Striatum", neuron_type="LTS", neuron_position=[0, 0, 0.5e-3])
        self.assertEqual(np.linalg.det(rot), 1)

        vec_rot2 = rot2.dot(vec_test)
        vec_facit2 = np.array([-1, -1, 0])
        self.assertTrue(np.allclose(vec_rot2 / np.linalg.norm(vec_rot2), vec_facit2 / np.linalg.norm(vec_facit2)))

        # TODO: Also need to test the random rotations


if __name__ == '__main__':
    unittest.main()
