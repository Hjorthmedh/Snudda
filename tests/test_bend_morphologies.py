import numpy as np
import unittest

from snudda.neurons import NeuronMorphologyExtended
from snudda.place.bend_morphologies import BendMorphologies


class TestBendMorphologies(unittest.TestCase):
    def test_something(self):

        morph_path = "validation/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20190508/WT-0728MSN01-cor-rep-ax.swc"
        mesh_path = "validation/mesh/Striatum-d-right.obj"
        nm = NeuronMorphologyExtended(swc_filename=morph_path)

        bm = BendMorphologies(mesh_path, rng=np.random.default_rng(1))

        pos = np.array([0.006, 0.004, 0.00205])

        before = nm.clone(position=pos, rotation=np.eye(3))
        after = nm.clone(position=pos, rotation=np.eye(3))

        new_rot_rep, _ = bm.bend_morphology(after.get_morphology())
        new_coord = bm.apply_rotation(after.get_morphology(), new_rot_rep)
        after.get_morphology().geometry[:, :3] = new_coord

        change = np.sum(np.abs(before.get_morphology().geometry[:, :3] - after.get_morphology().geometry[:, :3]))
        print(f"Change = {change}")

        return

        import pdb
        pdb.set_trace()

        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
