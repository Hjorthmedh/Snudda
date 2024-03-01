import os.path
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
        self.n_iter = 500

        self.sd = SnuddaDetect(hyper_voxel_size=160)  # Use 250 for real neurons
        self.sd.hyper_voxels = {0: {"random_seed": 1234}}
        self.sd.max_axon = 35  # 27 -- for real morphologies
        self.sd.max_dend = 30

        centre_pos = self.sd.num_bins * self.sd.voxel_size / 2

        self.sd.setup_hyper_voxel(hyper_voxel_id=0, hyper_voxel_origo=np.zeros((3,)))

        # neuron_file = "validation/shortballandstick/simpleshort.swc"
        # neuron_file = "validation/shortballandstick/simpleshortshort.swc"

        neuron_file = "validation/striatum-var/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026/morphology/WT-0728MSN01-cor-rep-ax-res3-var8.swc"

        # neuron_file = "../../BasalGangliaData/data/neurons/striatum/dspn/str-dspn-e150917_c6_D1-m21-6-DE-v20211028/morphology/21-6-DE-cor-rep-ax-res3.swc"
        # neuron_file = "../../BasalGangliaData/data/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026/morphology/WT-0728MSN01-cor-rep-ax-res3.swc"
        # neuron_file = "../../BasalGangliaData/data/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026/morphology/WT-0728MSN01-cor-rep-ax-res3-var2.swc"

        self.neuron = NeuronMorphologyExtended(name="simple",
                                               swc_filename=neuron_file)

        rng = self.sd.hyper_voxel_rng
        neuron_pos = rng.uniform(low=centre_pos[0]-40e-6,
                                 high=centre_pos[0]+40e-6, size=(self.n_iter, 3))
        neuron_clones = [None]

        # Clone neuron and add multiple copies with different rotation centred in center
        for neuron_id in range(1, self.n_iter+1):
            print(f"Adding neuron {neuron_id}")

            if neuron_id == 1:
                clone_rotation = np.eye(3)
            else:
                clone_rotation = SnuddaRotate.rand_rotation_matrix(rng=rng)
            neuron_clones.append(self.neuron.clone(position=neuron_pos[neuron_id-1, :],
                                                   rotation=clone_rotation))

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

    def test_voxel_counts(self):

        if self.sd.voxel_overflow_counter > 0:
            print(f"Warning voxel overflow counter: {self.sd.voxel_overflow_counter}")

        dend_neuron_id, dend_count = np.unique(self.sd.dend_voxels.flatten(), return_counts=True)
        axon_neuron_id, axon_count = np.unique(self.sd.axon_voxels.flatten(), return_counts=True)

        min_dend_count = np.min(dend_count[1:])
        min_axon_count = np.min(axon_count[1:])
        max_dend_count = np.max(dend_count[1:])
        max_axon_count = np.max(axon_count[1:])

        dend_bins = np.arange(min_dend_count-1.5, max_dend_count+1.5)
        axon_bins = np.arange(min_axon_count-1.5, max_axon_count+1.5)

        if True:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(2, 1)
            ax[0].hist(dend_count[1:], bins=dend_bins, align="mid", color="black")  # 0 is empty elements
            ax[1].hist(axon_count[1:], bins=axon_bins, align="mid", color="blue")

            ax[0].plot(np.full((2,), dend_count[1]), ax[0].get_ylim(), color="red")
            ax[1].plot(np.full((2,), axon_count[1]), ax[1].get_ylim(), color="red")

            ax[0].plot(np.full((2,), np.mean(dend_count[1:])), ax[0].get_ylim(), color="black")
            ax[1].plot(np.full((2,), np.mean(axon_count[1:])), ax[1].get_ylim(), color="blue")

            ax[1].set_xlabel("Number of voxels")
            ax[0].set_ylabel("Count")
            ax[1].set_ylabel("Count")
            ax[0].set_title(f"Unrotated dendrite voxel count: {dend_count[1]} (spread {np.std(dend_count[1:]):.1f} relspread {np.std(dend_count[1:])/np.mean(dend_count[1:]):.5f} )")
            ax[1].set_title(f"Unrotated axon voxel count: {axon_count[1]} (spread {np.std(axon_count[1:]):.1f}, relspread {np.std(axon_count[1:])/np.mean(axon_count[1:]):.5f}))")

            fig.suptitle(f"Step multiplier: {self.sd.step_multiplier} (overflow: {self.sd.voxel_overflow_counter})")

            plt.tight_layout()

            if not os.path.isdir("figures"):
                os.mkdir("figures")

            fig_name = f"figures/detect_voxel_count-step-multiplier-{self.sd.step_multiplier}-num-{self.n_iter}.png"
            plt.savefig(fig_name)

            # plt.ion()
            #   plt.show()

            if False:
                print("Plotting hyper voxel (this is slow)")
                self.sd.plot_hyper_voxel()
                import pdb
                pdb.set_trace()

        self.assertTrue(np.std(dend_count[1:])/np.mean(dend_count[1:]) < 0.05)
        self.assertTrue(np.std(axon_count[1:])/np.mean(axon_count[1:]) < 0.05)

        # import pdb
        # pdb.set_trace()


if __name__ == '__main__':
    unittest.main()
