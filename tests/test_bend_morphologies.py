import os
import numpy as np
import unittest

from snudda.neurons import NeuronMorphologyExtended
from snudda.place.bend_morphologies import BendMorphologies
import snudda.place.rotation

class TestBendMorphologies(unittest.TestCase):
    def test_something(self):

        base_path = os.path.dirname(__file__)

        morph_path = os.path.join(base_path, "validation", "striatum", "dspn",
                                  "str-dspn-e150602_c1_D1-mWT-0728MSN01-v20190508",
                                  "WT-0728MSN01-cor-rep-ax.swc")
        # mesh_path = os.path.join(base_path, "validation", "mesh", "Striatum-d-right.obj")
        mesh_path = os.path.join(base_path, "validation", "mesh", "bend_test_mesh.obj")

        nm = NeuronMorphologyExtended(swc_filename=morph_path)
        bm = BendMorphologies(mesh_path, rng=np.random.default_rng(12345))

        # pos = np.array([0.006, 0.004, 0.00205])
        pos = np.array([400e-6, 0, 0])


        rng = np.random.default_rng(1234)

        rot_mat = snudda.place.rotation.SnuddaRotate.rand_rotation_matrix(rng=rng)

        before = nm.clone(position=pos, rotation=rot_mat) # np.eye(3)
        after = nm.clone(position=pos, rotation=rot_mat) # np.eye(3)

        before_morph = before.get_morphology()
        old_rot_rep = bm.get_full_rotation_representation(morphology=before_morph)
        coords = bm.apply_rotation(morphology=before_morph, rotation_representation=old_rot_rep)

        # Verify that rotation representation works
        self.assertTrue((np.abs(before_morph.geometry[:, :3] - coords) < 1e-6).all())

        new_rot_rep, _ = bm.bend_morphology(after.get_morphology())

        new_coord = bm.apply_rotation(after.get_morphology(), new_rot_rep)
        after.get_morphology().geometry[:, :3] = new_coord

        change = np.sum(np.abs(before.get_morphology().geometry[:, :3] - after.get_morphology().geometry[:, :3]))
        print(f"Change = {change}")

        before_inside = bm.check_if_inside(before)
        after_inside = bm.check_if_inside(after)

        n_before = np.sum(before_inside)
        n_after = np.sum(after_inside)
        n_all = len(before_inside)

        # The bending is statistical, so we some parts of neuron might go a little outside
        self.assertTrue(n_after > n_before)

        self.assertTrue(n_after - n_before > 0.8 * (n_all - n_before))



if __name__ == '__main__':
    unittest.main()
