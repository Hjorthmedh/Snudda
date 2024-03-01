import unittest

import os
import numpy as np

from snudda import SnuddaInit
from snudda import SnuddaPlace
from snudda import SnuddaDetect
from snudda import SnuddaPrune
from snudda import SnuddaLoad
from snudda.input import SnuddaInput
from snudda.utils.export_sonata import ExportSonata

from sonata.circuit import File as SonataFile


class TestSonata(unittest.TestCase):

    def setUp(self, create_network=True):

        self.network_path = os.path.join("networks", "sonata_example")

        # !!! TEMP SKIP setup while developing
        # print("SKIPPING NETWORK CREATIONG DURING DEVELOPMENT")
        # return

        if create_network:
            si = SnuddaInit(network_path=self.network_path, random_seed=12345)
            # Not the correct proportions, this is a unit test
            si.define_striatum(num_dSPN=10, num_iSPN=10, num_FS=10, num_LTS=2, num_ChIN=2,
                               volume_type="cube")
            si.write_json()

            sp = SnuddaPlace(network_path=self.network_path)
            sp.place()

            sd = SnuddaDetect(network_path=self.network_path)
            sd.detect()

            sp = SnuddaPrune(network_path=self.network_path)
            sp.prune()

            input_config_file = os.path.join("networks", "network_testing_input","input-test-1.json")
            si = SnuddaInput(network_path=self.network_path, input_config_file=input_config_file)
            si.generate()

    def test_sonata_export(self):

        network_file = os.path.join(self.network_path, "network-synapses.hdf5")
        input_file = os.path.join(self.network_path, "input-spikes.hdf5")
        out_dir = os.path.join(self.network_path, "SONATA")
        se = ExportSonata(network_file=network_file, input_file=input_file, out_dir=out_dir)

        data_files = [os.path.join(out_dir, "networks", "Striatum_nodes.hdf5"),
                      os.path.join(out_dir, "networks", "Striatum_edges.hdf5")]

        data_type_files = [os.path.join(out_dir, "networks", "Striatum_node_types.csv"),
                           os.path.join(out_dir, "networks", "Striatum_edge_types.csv")]

        # Load Sonata network
        sf = SonataFile(data_files=data_files, data_type_files=data_type_files)

        # Load Snudda reference network
        sl = SnuddaLoad(network_file=network_file)

        with self.subTest("Check nodes"):

            neuron_populations = set(['ChIN', 'FS', 'LTS', 'dSPN', 'iSPN'])

            self.assertEqual(set(sf.nodes.population_names), neuron_populations)
            nodes = dict()

            for pop_name in neuron_populations:

                nodes[pop_name] = sf.nodes.get_population(pop_name)

            # SONATA seems to have separate node_id for each population
            neuron_types = np.array([x["type"] for x in sl.data["neurons"]])
            within_type_idx = np.full(shape=(len(neuron_types),), fill_value=-1, dtype=int)
            for nt in set(neuron_types):
                idx = np.where(neuron_types == nt)[0]
                within_type_idx[idx] = np.arange(0, len(idx))

            assert (within_type_idx >= 0).all()

            for neuron in sl.data["neurons"]:
                neuron_type = neuron["type"]

                try:
                    sonata_node = nodes[neuron_type].get_node_id(within_type_idx[neuron["neuron_id"]])
                except:
                    import traceback
                    print(traceback.format_exc())
                    import pdb
                    pdb.set_trace()

                # Remember to convert from natural units to SI units
                self.assertAlmostEqual(neuron["position"][0], sonata_node["x"]*1e-6, 8)
                self.assertAlmostEqual(neuron["position"][1], sonata_node["y"]*1e-6, 8)
                self.assertAlmostEqual(neuron["position"][2], sonata_node["z"]*1e-6, 8)

                # TODO: We need to also check the rotation

                self.assertEqual(os.path.basename(neuron["morphology"]).replace(".swc", ""), sonata_node["morphology"])

                if se.target_simulator != "NEST":
                    self.assertEqual(sonata_node["model_type"], "biophysical")

                self.assertEqual(neuron["type"], sonata_node["model_name"].split("_")[0])

        with self.subTest("Check edges"):

            edge_count = 0

            con_mat = sl.create_connection_matrix()
            new_con_mat = np.zeros(shape=con_mat.shape, dtype=int)

            neuron_types = np.array([x["type"] for x in sl.data["neurons"]])
            neuron_idx = dict()
            for nt in set(neuron_types):
                neuron_idx[nt] = np.where(neuron_types == nt)[0]

            for edge_pop_name in sf.edges.population_names:
                edges_pop = sf.edges[edge_pop_name]

                pre_pop, post_pop = edge_pop_name.split("_")

                for edge in edges_pop:
                    # print(f"{edge}")
                    try:
                        src_id = neuron_idx[pre_pop][edge.source_node_id]
                        target_id = neuron_idx[post_pop][edge.target_node_id]
                        new_con_mat[src_id, target_id] += 1

                        edge_count += 1
                    except:
                        import traceback
                        print(traceback.format_exc())
                        import pdb
                        pdb.set_trace()

            # Check that all edges are accounted for

            self.assertEqual(edge_count, sl.data["num_synapses"])
            self.assertTrue((con_mat == new_con_mat).all())

        # !!! TODO: Add test for input, also write dedicated input.json file for this test


if __name__ == '__main__':
    unittest.main()
