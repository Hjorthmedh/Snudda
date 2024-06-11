import os
import unittest
import numpy as np

from snudda.init.init import SnuddaInit
from snudda import SnuddaPlace
from snudda import SnuddaDetect
from snudda import SnuddaPrune
from snudda.utils import SnuddaLoad
from snudda.utils.ablate_network import SnuddaAblateNetwork


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:

        self.network_path = os.path.join("networks", "ablation")
        self.modified_network = os.path.join(self.network_path, "network-synapses-ablated.hdf5")

        rc = None  # Run in serial for now

        bg_data = os.path.join("..", "..", "BasalGangliaData", "data")
        if os.path.exists(bg_data):
            snudda_data = bg_data
        else:
            snudda_data = None

        cnc = SnuddaInit(struct_def={}, network_path=self.network_path, random_seed=1337,
                         snudda_data=snudda_data)

        self.num_fs = 50
        self.num_dspn_orig = 100  # 1000
        self.num_ispn_orig = 100  # 1000
        self.num_dspn_final = 20
        self.num_ispn_final = 20

        cnc.define_striatum(num_dSPN=self.num_dspn_orig, num_iSPN=self.num_ispn_orig,
                            num_FS=self.num_fs, num_LTS=0, num_ChIN=0,
                            volume_type="cube",
                            neuron_density=0.013 * 80500 * (1 + (self.num_dspn_orig + self.num_ispn_orig) / self.num_fs))
        cnc.add_population_unit_random("Striatum", "FS", 1.0, 1)

        # Overwrite adult GJ connectivity, with higher juvenile connectivity
        cnc.add_neuron_target(neuron_name="FS",
                              target_name="FS",
                              region_name="Striatum",
                              connection_type="gap_junction",
                              dist_pruning=None,
                              f1=0.7, soft_max=8, mu2=2, a3=1.0,
                              conductance=[0.5e-9, 0.1e-9],
                              cluster_synapses=False,
                              channel_param_dictionary=None)

        cnc.write_json()
        del cnc

        sp = SnuddaPlace(network_path=self.network_path, verbose=False)
        sp.place()
        del sp

        sd = SnuddaDetect(network_path=self.network_path, verbose=False, rc=rc)
        sd.detect()
        del sd

        sp = SnuddaPrune(network_path=self.network_path, rc=rc)
        sp.prune()
        del sp

        # Ablate the network
        mod_network = SnuddaAblateNetwork(network_file=os.path.join(self.network_path, "network-synapses.hdf5"))
        orig_sl = SnuddaLoad(self.network_path)
        orig_fs_id = orig_sl.get_neuron_id_of_type("FS")
        orig_dspn_id = [x for x, y in orig_sl.get_centre_neurons_iterator(neuron_type="dSPN",
                                                                          n_neurons=self.num_dspn_final)]
        orig_ispn_id = [x for x, y in orig_sl.get_centre_neurons_iterator(neuron_type="iSPN",
                                                                          n_neurons=self.num_ispn_final)]

        self.keep_id = set(list(orig_fs_id) + list(orig_dspn_id) + list(orig_ispn_id))
        mod_network.only_keep_neuron_id(self.keep_id)
        mod_network.write_network(self.modified_network)

        if False:
            from snudda.plotting.plot_network import PlotNetwork
            pn = PlotNetwork(self.modified_network)
            pn.plot(plot_axon=False, plot_dendrite=False, plot_synapses=False)

    def test_mapping(self):

        # We started with 50 + 2x 1000 = 2050 neurons, and kept only self.keep_id
        new_id_lookup = dict()
        old_id_lookup = dict()
        for new_id, old_id in enumerate(sorted(list(self.keep_id))):
            new_id_lookup[old_id] = new_id
            old_id_lookup[new_id] = old_id

        keep_mask = np.zeros((self.num_fs + self.num_dspn_orig + self.num_ispn_orig,), dtype=bool)
        keep_mask[np.array(sorted(list(self.keep_id)))] = True

        sl_orig = SnuddaLoad(self.network_path)
        sl_ablated = SnuddaLoad(self.modified_network)

        for new_id, old_id in old_id_lookup.items():
            # Check synapses ok
            new_synapses, _ = sl_ablated.find_synapses(post_id=new_id)
            orig_synapses_all, _ = sl_orig.find_synapses(post_id=old_id)

            if orig_synapses_all is not None:
                # Remove synapses coming from neurons removed
                orig_synapses = orig_synapses_all[keep_mask[orig_synapses_all[:, 0]], :]

                if new_synapses is not None or orig_synapses is not None:
                    # Pre-synaptic neurons should match
                    self.assertTrue((new_synapses[:, 0] == [new_id_lookup[x] for x in orig_synapses[:, 0]]).all())
                    self.assertTrue((new_synapses[:, 1] == new_id).all())
                    self.assertTrue((orig_synapses[:, 1] == old_id).all())
                    self.assertTrue((new_synapses[:, 2:] == orig_synapses[:, 2:]).all())

            # Check gap junctions
            new_gj, _ = sl_ablated.find_gap_junctions(neuron_id=new_id)
            old_gj_all, _ = sl_orig.find_gap_junctions(neuron_id=old_id)

            # Remove GJ coming from neurons removed
            gj_mask = np.logical_and(keep_mask[old_gj_all[:, 0]], keep_mask[old_gj_all[:, 1]])
            try:
                old_gj = old_gj_all[np.where(gj_mask)[0], :]
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

            if new_gj is not None or old_gj is not None:
                self.assertTrue((new_gj[:, 0] == [new_id_lookup[x] for x in old_gj[:, 0]]).all())
                self.assertTrue((new_gj[:, 1] == [new_id_lookup[x] for x in old_gj[:, 1]]).all())
                self.assertTrue((new_gj[:, 2:] == old_gj[:, 2:]).all())

        # Check all the neuron information also
        for new_id, old_id in old_id_lookup.items():

            for var in sl_orig.data["neurons"][old_id]:
                if var == "neuron_id":
                    self.assertTrue(sl_ablated.data["neurons"][new_id]["neuron_id"] == new_id)
                    self.assertTrue(sl_orig.data["neurons"][old_id]["neuron_id"] == old_id)
                elif var == "axon_density_radius":
                    self.assertTrue((np.isnan(sl_ablated.data["neurons"][new_id][var]) and np.isnan(sl_orig.data["neurons"][old_id][var]))
                                    or sl_ablated.data["neurons"][new_id][var] == sl_orig.data["neurons"][old_id][var])

                elif type(sl_orig.data["neurons"][old_id][var]) in [str, type(None), np.float64]:
                    try:
                        self.assertTrue(sl_ablated.data["neurons"][new_id][var] == sl_orig.data["neurons"][old_id][var])
                    except:
                        print("Tell me why?")
                        import traceback
                        print(traceback.format_exc())
                        import pdb
                        pdb.set_trace()
                else:
                    try:
                        self.assertTrue((sl_ablated.data["neurons"][new_id][var] \
                                         == sl_orig.data["neurons"][old_id][var]).all())
                    except:
                        print("Tell me why?")
                        import traceback
                        print(traceback.format_exc())
                        import pdb
                        pdb.set_trace()


if __name__ == '__main__':
    unittest.main()
