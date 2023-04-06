import unittest

import os
import numpy as np

from snudda import SnuddaInit
from snudda import SnuddaPlace
from snudda import SnuddaDetect
from snudda import SnuddaPrune
from snudda import SnuddaLoad
from snudda.utils.export_sonata import ExportSonata

from sonata.circuit import File as SonataFile


class TestSonata(unittest.TestCase):

    def setUp(self, create_network=True):

        self.network_path = os.path.join("networks", "sonata_example")

        # !!! TEMP SKIP setup while developing
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
            self.assertEqual(sf.nodes.population_names, ["Striatum"])
            nodes = sf.nodes.get_population("Striatum")

            for neuron in sl.data["neurons"]:
                sonata_node = nodes.get_node_id(neuron["neuronID"])

                # Remember to convert from natural units to SI units
                self.assertAlmostEqual(neuron["position"][0], sonata_node["x"]*1e-6, 8)
                self.assertAlmostEqual(neuron["position"][1], sonata_node["y"]*1e-6, 8)
                self.assertAlmostEqual(neuron["position"][2], sonata_node["z"]*1e-6, 8)

                # TODO: We need to also check the rotation

                self.assertEqual(os.path.basename(neuron["morphology"]), sonata_node["morphology"])
                self.assertEqual(sonata_node["model_type"], "biophysical")
                self.assertEqual(neuron["type"], sonata_node["model_name"].split("_")[0])

        with self.subTest("Check edges"):

            edge_count = 0

            con_mat = sl.create_connection_matrix()
            new_con_mat = np.zeros(shape=con_mat.shape, dtype=int)

            for edge_pop_name in sf.edges.population_names:
                edges_pop = sf.edges[edge_pop_name]

                for edge in edges_pop:
                    # print(f"{edge}")
                    src_id = edge.source_node_id
                    target_id = edge.target_node_id
                    new_con_mat[src_id, target_id] += 1

                    edge_count += 1

            # Check that all edges are accounted for

            self.assertEqual(edge_count, sl.data["nSynapses"])
            self.assertTrue((con_mat == new_con_mat).all())


if __name__ == '__main__':
    unittest.main()