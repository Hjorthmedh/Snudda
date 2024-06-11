import unittest
import os
import numpy as np

from snudda.utils.ablate_network import SnuddaAblateNetwork
from snudda.utils.load import SnuddaLoad
from snudda.analyse import SnuddaAnalyse

#
# This unit test uses a pre-generated network, and checks that neurons and connections in it are correctly removed
#
# To regenerate the network, stand in the Snudda/tests directory and do:
#
#    snudda init networks/modify_network --size 200
#    snudda place networks/modify_network
#    snudda detect networks/modify_network
#    snudda prune networks/modify_network


class SnuddaModifyNetworkTestCase(unittest.TestCase):

    def setUp(self) -> None:

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        self.network_path = os.path.join("networks", "modify_network")
        self.original_file = os.path.join(self.network_path, "network-synapses.hdf5")

        if not os.path.exists(self.original_file):
            print("Missing test network, generating it first.")
            from snudda.init import SnuddaInit
            SnuddaInit(network_path=self.network_path, struct_def={"Striatum": 200}, random_seed=999)

            from snudda.place import SnuddaPlace
            sp = SnuddaPlace(network_path=self.network_path)
            sp.place()

            from snudda.detect import SnuddaDetect, SnuddaPrune
            sd = SnuddaDetect(network_path=self.network_path)
            sd.detect()

            spp = SnuddaPrune(network_path=self.network_path)
            spp.prune()

        self.snudda_load_original = SnuddaLoad(network_file=self.original_file, load_synapses=True)
        self.original_type_count = self.get_number_of_type(self.snudda_load_original.data)

    def test_remove_neurons(self):

        """ Verify that neuron type is removed correctly """

        profile_run = False

        new_file_dspn = os.path.join(self.network_path, "dspn-neurons-removed.hdf5")
        new_file_ispn = os.path.join(self.network_path, "ispn-neurons-removed.hdf5")
        new_file_ispn_dspn = os.path.join(self.network_path, "ispn-dspn-connections-removed.hdf5")

        mod_network = SnuddaAblateNetwork(network_file=self.original_file)
        sa_orig = SnuddaAnalyse(hdf5_file=self.original_file)
        dspn_id_orig = self.snudda_load_original.get_neuron_id_of_type(neuron_type="dSPN")
        ispn_id_orig = self.snudda_load_original.get_neuron_id_of_type(neuron_type="iSPN")

        with self.subTest(msg="Testing removal of iSPN to dSPN synapses"):
            mod_network.reset_network()
            mod_network.remove_connection(pre_neuron_type="iSPN", post_neuron_type="dSPN")
            mod_network.write_network(out_file_name=new_file_ispn_dspn)

            sa_new = SnuddaAnalyse(hdf5_file=new_file_ispn_dspn)

            print("Comparing connection matrices")

            ref_mat = sa_orig.connection_matrix.copy()
            for i_id in ispn_id_orig:
                ref_mat[i_id, dspn_id_orig] = 0

            # When all ispn to dspn connections are removed in ref matrix, it should match the new connection matrix
            print("Checking that new connection matrix match reference")
            self.assertTrue((sa_new.connection_matrix != ref_mat).nnz == 0)

        with self.subTest(msg="Testing removal of neuron type dSPN"):
            mod_network.reset_network()
            mod_network.remove_neuron_type(neuron_type="dSPN")

            if profile_run:
                import cProfile
                prof_file = "profiling-info.prof"
                cProfile.runctx("mod_network.write_network(out_file_name=new_file_dspn)",
                                None, locals(), filename=prof_file)

                # To analyse profile data:
                import pstats
                from pstats import SortKey
                p = pstats.Stats(prof_file)
                p.strip_dirs().sort_stats(SortKey.CUMULATIVE).print_stats(100)

            else:
                mod_network.write_network(out_file_name=new_file_dspn)

            snudda_load_new = SnuddaLoad(network_file=new_file_dspn, load_synapses=True)
            new_type_count = self.get_number_of_type(snudda_load_new.data)

            for neuron_type in self.original_type_count:
                if neuron_type == "dSPN":
                    self.assertTrue(neuron_type not in new_type_count)
                else:
                    self.assertEqual(new_type_count[neuron_type], self.original_type_count[neuron_type])

            sa_new = SnuddaAnalyse(hdf5_file=new_file_dspn)
            # Remove dSPN from original connection matrix and check it matches with new connection matrix
            keep_idx = list(set(np.arange(0, len(self.snudda_load_original.data["neurons"]))) - set(dspn_id_orig))
            ref_mat = sa_orig.connection_matrix.copy()[keep_idx, :][:, keep_idx]

            self.assertTrue((sa_new.connection_matrix != ref_mat).nnz == 0)

        with self.subTest(msg="Testing removal of neuron type iSPN, p=0.3"):
            mod_network.reset_network()
            mod_network.remove_neuron_type(neuron_type="iSPN", p_remove=0.3)
            mod_network.write_network(out_file_name=new_file_ispn)

            snudda_load_new = SnuddaLoad(network_file=new_file_ispn, load_synapses=True)
            new_type_count = self.get_number_of_type(snudda_load_new.data)

            for neuron_type in self.original_type_count:
                if neuron_type == "iSPN":
                    # Check that we are within a range of 30% removed
                    self.assertTrue(0.5 * self.original_type_count[neuron_type]
                                    < new_type_count[neuron_type]
                                    < 0.9 * self.original_type_count[neuron_type])
                else:
                    self.assertEqual(new_type_count[neuron_type], self.original_type_count[neuron_type])

    def get_number_of_type(self, snudda_load_data):

        type_counting = dict()

        neuron_types = [n["type"] for n in snudda_load_data["neurons"]]
        unique_neuron_types = set(neuron_types)

        for unt in unique_neuron_types:
            type_counting[unt] = np.sum([nt == unt for nt in neuron_types])

        return type_counting

        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
