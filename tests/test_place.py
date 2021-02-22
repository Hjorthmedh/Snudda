import unittest
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


class TestPlace(unittest.TestCase):

    def setUp(self):
        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))
        print(f"Current directory (detect): {os.path.dirname(os.path.realpath(__file__))}")

        from snudda.init import SnuddaInit
        cell_spec = os.path.join(os.path.dirname(__file__), "validation")
        self.sim_name = os.path.join("networks", "network_testing_place")
        self.config_file = os.path.join(self.sim_name, "network-config.json")
        cnc = SnuddaInit(struct_def={}, config_file=self.config_file, num_population_units=1, random_seed=1234)
        cnc.define_striatum(num_dSPN=10, num_iSPN=0, num_FS=10, num_LTS=0, num_ChIN=0,
                            volume_type="cube", cell_spec_dir=cell_spec)
        cnc.write_json(self.config_file)

    def tearDown(self):
        pass

    def test_init(self):
        with open(self.config_file, "r") as f:
            config_data = json.load(f)

        self.assertTrue("RandomSeed" in config_data)
        self.assertTrue(config_data["RandomSeed"]["masterseed"] == 1234)
        self.assertTrue(config_data["RandomSeed"]["init"] == 2906597030)
        self.assertTrue(config_data["RandomSeed"]["place"] == 1602421836)
        self.assertTrue(config_data["RandomSeed"]["detect"] == 216676975)
        self.assertTrue(config_data["RandomSeed"]["prune"] == 2698621808)
        self.assertTrue(config_data["RandomSeed"]["input"] == 507409703)
        self.assertTrue(config_data["RandomSeed"]["simulate"] == 2825158027)

        self.assertTrue("Volume" in config_data)
        self.assertTrue("PopulationUnits" in config_data)
        self.assertTrue("Connectivity" in config_data)
        self.assertTrue("Neurons" in config_data)

        for neuron in config_data["Neurons"]:
            self.assertTrue("morphology" in config_data["Neurons"][neuron])
            self.assertTrue("parameters" in config_data["Neurons"][neuron])
            self.assertTrue("mechanisms" in config_data["Neurons"][neuron])
            self.assertTrue("modulation" in config_data["Neurons"][neuron])
            self.assertTrue("num" in config_data["Neurons"][neuron])

    def test_place(self):

        # Place neurons
        from snudda.place import SnuddaPlace

        position_file = os.path.join(self.sim_name, "network-neuron-positions.hdf5")

        npn = SnuddaPlace(config_file=self.config_file,
                          log_file=None,
                          verbose=True,
                          d_view=None,          # TODO: If d_view is None code run sin serial, add test parallel
                          h5libver="latest")

        npn.read_config()
        npn.write_data(position_file)

        with self.subTest(stage="neuron_count"):
            num_cells = 20
            self.assertEqual(npn.all_neuron_positions().shape[0], num_cells)
            self.assertEqual(len(npn.neurons), num_cells)

        with self.subTest(stage="d_min"):
            # Check that minimum distance between neurons, d_min is fulfilled
            all_pos = npn.all_neuron_positions()
            for pos in all_pos:
                d_min = 1e-5
                # Not too close (self comparison allowed to be equal, hence -1)
                self.assertTrue(np.sum(np.linalg.norm(all_pos - pos, axis=1) >= d_min) == all_pos.shape[0]-1)

                # Not too far apart
                d_max = 1e4  # This is just for this simulation, to catch if there are some spurious neurons
                self.assertTrue(np.sum(np.linalg.norm(all_pos - pos, axis=1) <= d_max) == all_pos.shape[0])

        # TODO: Load hdf5 file and check that data is what we expect


if __name__ == '__main__':
    unittest.main()

# python3 -m unittest test_place