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
        n_iter = 10000

        self.sd = SnuddaDetect()
        self.sd.hyper_voxels = {0: {"randomSeed": 123}}

        self.sd.setup_hyper_voxel(hyper_voxel_id=0, hyper_voxel_origo=np.zeros((3,)))

        self.neuron = NeuronMorphologyExtended(name="simple",
                                               swc_filename="validation/shortballandstick/simpleshort.swc")

        neuron_pos = np.full((1, 3), 150e-6)
        neuron_clones = [None]

        # Clone neuron and add multiple copies with different rotation centred in center
        for neuron_id in range(1, n_iter+1):
            print(f"Adding neuron {neuron_id}")

            if neuron_id == 1:
                clone_rotation = np.eye(3)
            else:
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
        print("TESTING")
        self.assertEqual(True, True)

    def test_voxel_counts(self):

        dend_neuron_id, dend_count = np.unique(self.sd.dend_voxels.flatten(), return_counts=True)
        axon_neuron_id, axon_count = np.unique(self.sd.axon_voxels.flatten(), return_counts=True)

        if True:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(2, 1)
            ax[0].hist(dend_count[1:], align="right", color="black")  # 0 is empty elements
            ax[1].hist(axon_count[1:], align="right", color="blue")

            ax[0].plot(np.full((2,), dend_count[1]), ax[0].get_ylim(), color="red")
            ax[1].plot(np.full((2,), axon_count[1]), ax[1].get_ylim(), color="red")

            ax[1].set_xlabel("Number of voxels")
            ax[0].set_ylabel("Count")
            ax[1].set_ylabel("Count")
            ax[0].set_title(f"Unrotated dendrite voxel count: {dend_count[1]}")
            ax[1].set_title(f"Unrotated axon voxel count: {axon_count[1]}")

            fig.suptitle(f"Step multiplier: {self.sd.step_multiplier}")

            plt.tight_layout()

            fig_name = f"figures/detect_voxel_count-step-multiplier-{self.sd.step_multiplier}.png"
            plt.savefig(fig_name)

            plt.ion()
            plt.show()



            # self.sd.plot_hyper_voxel()

        print("TESTME")

        # import pdb
        # pdb.set_trace()

        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
