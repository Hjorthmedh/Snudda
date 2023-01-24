import unittest

import numpy as np
from snudda.detect import SnuddaDetect
from snudda.neurons.neuron_morphology_extended import NeuronMorphologyExtended
from snudda.place.rotation import SnuddaRotate

# We want to create a hyper voxel, then place simple ball and stick neurons to see that
# regardless of the direction they are rotated we get (approximately) the same number of
# voxels marked during touch detection


class TestVoxelCount(unittest.TestCase):

    def setUp(self) -> None:
        n_iter = 100

        self.sd = SnuddaDetect()
        self.hyper_voxels = {0: {"randomSeed": 123}}

        self.sd.setup_hyper_voxel(hyper_voxel_id=0, hyper_voxel_origo=np.zeros((3,)))

        self.neuron = NeuronMorphologyExtended(name="simple",
                                               swc_filename="validation/shortballandstick/simpleshort.swc")

        neuron_pos = np.full((3, 1), 50e-6)
        neuron_clones = []

        # Clone neuron and add multiple copies with different rotation centred in center
        for neuron_id in range(n_iter):
            clone_rotation = SnuddaRotate.rand_rotation_matrix()
            neuron_clones.append(self.neuron.clone(position=neuron_pos, rotation=clone_rotation))

            self.sd.fill_voxels_dend(voxel_space=self.sd.dend_voxels,
                                     voxel_space_ctr=self.sd.dend_voxel_ctr,
                                     voxel_sec_id=self.sd.dend_sec_id,
                                     voxel_sec_x=self.sd.dend_sec_x,
                                     voxel_soma_dist=self.sd.dend_soma_dist,
                                     neuron=neuron_clones[neuron_id],
                                     neuron_id=neuron_id)

            self.sd.fill_voxels_axon(voxel_space=self.sd.axon_voxels,
                                     voxel_space_ctr=self.sd.axon_voxel_ctr,
                                     voxel_axon_dist=self.sd.axon_soma_dist,
                                     neuron=neuron_clones[neuron_id],
                                     neuron_id=neuron_id)

    def test_something(self):

        import pdb
        pdb.set_trace()

        pass
        # self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
