import unittest

import numpy as np
from snudda.detect import SnuddaDetect
from snudda.neurons.neuron_morphology_extended import NeuronMorphologyExtended

# We want to create a hyper voxel, then place simple ball and stick neurons to see that
# regardless of the direction they are rotated we get (approximately) the same number of
# voxels marked during touch detection


class TestVoxelCount(unittest.TestCase):

    def setUp(self) -> None:
        self.sd = SnuddaDetect()
        self.hyper_voxels = {0: {"randomSeed": 123}}

        self.sd.setup_hyper_voxel(hyper_voxel_id=0, hyper_voxel_origo=np.zeros((3,)))

        self.neuron = NeuronMorphologyExtended(name="simple",
                                               swc_filename="validation/shortballandstick/simpleshort.swc")

        # Clone neuron and add multiple copies with different rotation centred in center

    def test_something(self):

        pass
        # self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
