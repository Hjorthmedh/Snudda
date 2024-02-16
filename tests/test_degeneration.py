import os
import numpy as np
import h5py

import unittest

from snudda import SnuddaDetect, SnuddaPrune
from snudda.plotting import PlotNetwork
from snudda.utils.swap_to_degenerated_morphologies_extended import SwapToDegeneratedMorphologiesExtended


class MyTestCase(unittest.TestCase):

    def setUp(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        network_A = os.path.join("networks", "network_testing_degeneration", "A")
        network_B = os.path.join("networks", "network_testing_degeneration", "B")
        network_C = os.path.join("networks", "network_testing_degeneration", "C")

        self.network_A = network_A
        self.network_C = network_C

        from snudda.place import SnuddaPlace

        # Setup original network
        # https://www.americanscientist.org/article/the-quest-for-randomness
        sp = SnuddaPlace(network_path=network_A, random_seed=9)
        sp.place()

        network_A_pos_file = os.path.join(network_A, "network-neuron-positions.hdf5")
        self.set_neuron_positions(network_file=network_A_pos_file)

        sd = SnuddaDetect(network_path=network_A, random_seed=11)
        sd.detect()

        spr = SnuddaPrune(network_path=network_A, random_seed=13)
        spr.prune()

        # Setup degenerated network
        sp2 = SnuddaPlace(network_path=network_B, random_seed=9)  # same as original
        sp2.place()

        network_B_pos_file = os.path.join(network_B, "network-neuron-positions.hdf5")
        self.set_neuron_positions(network_file=network_B_pos_file)

        sd2 = SnuddaDetect(network_path=network_B, random_seed=11)
        sd2.detect()

        spr2 = SnuddaPrune(network_path=network_B, random_seed=13)
        spr2.prune()

        network_A_file = os.path.join(network_A, "network-synapses.hdf5")
        network_B_file = os.path.join(network_B, "network-synapses.hdf5")
        network_C_file = os.path.join(network_C, "network-synapses.hdf5")

        # These paths should normally be to the SNUDDA_DATA directory, and not directly to the neuron folders
        # but because we only have one neuron type in our model in this test, this will work... for now.
        original_data_dir = os.path.join("validation", "ballanddoublestick_original")
        updated_data_dir = os.path.join("validation", "ballanddoublestick_degenerated")

        swap = SwapToDegeneratedMorphologiesExtended(original_network_file=network_A_file,
                                                     updated_network_file=network_B_file,
                                                     output_network_file=network_C_file,
                                                     original_snudda_data_dir=original_data_dir,
                                                     updated_snudda_data_dir=updated_data_dir,
                                                     random_seed=15)

        swap.write_new_network_file()

        # Create config files for the network using ballanddoublestick, and ballanddoublestick_degenerated
        # (Here axons are elongated from 220 to 300 micrometers, and dendrites are shrunk from 200 to 150 micrometers)

        # Create placement files for the WT and degenerated networks, important that they both have the same neurons.

        # Do touch detection both networks

        # Do swap_to_degenerated_morphologies_extended

        pn = PlotNetwork(network_A)
        pn.plot(plot_axon=True, plot_dendrite=True, plot_synapses=True, fig_name="original-network.png")
        pn.close()

        pn2 = PlotNetwork(network_B)
        pn2.plot(plot_axon=True, plot_dendrite=True, plot_synapses=True, fig_name="degenerated-network.png")
        pn2.close()

        pn3 = PlotNetwork(network_C)
        pn3.plot(plot_axon=True, plot_dendrite=True, plot_synapses=True, fig_name="degenerated-network.png")
        pn3.close()

    def set_neuron_positions(self, network_file):

        with h5py.File(network_file, "r+") as hdf5_file:

            neuron_positions = np.array([[0, 20, 0],  # Postsynaptiska
                                         [0, 40, 0],
                                         [0, 60, 0],
                                         [0, 80, 0],
                                         [0, 100, 0],
                                         [0, 120, 0],
                                         [0, 140, 0],
                                         [0, 160, 0],
                                         [0, 180, 0],
                                         [0, 200, 0],
                                         [20, 0, 0],  # Presynaptiska
                                         [40, 0, 0],
                                         [60, 0, 0],
                                         [80, 0, 0],
                                         [100, 0, 0],
                                         [120, 0, 0],
                                         [140, 0, 0],
                                         [160, 0, 0],
                                         [180, 0, 0],
                                         [200, 0, 0],
                                         [100, 100, -240],  # To get a gap junction
                                         ]) * 1e-6

            hdf5_file["network/neurons/position"][:, :] = neuron_positions

            ang = -np.pi / 2
            R_x = np.array([[1, 0, 0],
                            [0, np.cos(ang), -np.sin(ang)],
                            [0, np.sin(ang), np.cos(ang)]])

            ang = np.pi / 2
            R_y = np.array([[np.cos(ang), 0, np.sin(ang)],
                            [0, 1, 0],
                            [-np.sin(ang), 0, np.cos(ang)]])

            ang = -np.pi / 2
            R_gj = np.array([[np.cos(ang), 0, np.sin(ang)],
                             [0, 1, 0],
                             [-np.sin(ang), 0, np.cos(ang)]])

            for idx in range(0, 10):  # Post synaptic neurons
                hdf5_file["network/neurons/rotation"][idx, :] = R_x.flatten()

            for idx in range(10, 20):  # Presynaptic neurons
                hdf5_file["network/neurons/rotation"][idx, :] = R_y.flatten()

                hdf5_file["network/neurons/rotation"][20, :] = R_gj.flatten()

            hdf5_file.close()

    def test_something(self):

        # Load the networks
        from snudda.utils import SnuddaLoad
        orig_load = SnuddaLoad(self.network_A)
        degen_load = SnuddaLoad(self.network_C)

        tmp = [(x["neuron_id"], x["morphology_key"]) for x in orig_load.data["neurons"]]
        print(f"Morphologies: {tmp}")

        #import pdb
        #pdb.set_trace()

        # self.assertEqual(orig_load.data["num_synapses"], 165)

        # Verify that it should be 99 synapses -- now it is just a regression test
        # Old version gave 99, new gives 165 --- CHECK WHY!

        # TODO: CHECK WHY NOT 99 SYNAPSES NOW
        # self.assertEqual(degen_load.data["num_synapses"], 155)  # -- Ilaria, we need to check what the true value should be?


if __name__ == '__main__':
    unittest.main()
