import unittest
import os
import sys
import json
import numpy as np

from snudda import SnuddaLoad

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from snudda.init.init import SnuddaInit
from snudda.place.place import SnuddaPlace


class TestPlace(unittest.TestCase):

    def setUp(self):
        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        print(f"Current directory (detect): {os.path.dirname(os.path.realpath(__file__))}")

        neuron_dir = os.path.join(os.path.dirname(__file__), "validation")
        self.sim_name = os.path.join("networks", "network_testing_place")
        self.config_file = os.path.join(self.sim_name, "network-config.json")
        cnc = SnuddaInit(struct_def={}, config_file=self.config_file, random_seed=1234)
        cnc.define_striatum(num_dSPN=10, num_iSPN=0, num_FS=10, num_LTS=0, num_ChIN=0,
                            volume_type="cube", neurons_dir=neuron_dir)
        cnc.write_json(self.config_file)

    def tearDown(self):
        pass

    def test_init(self):
        with open(self.config_file, "r") as f:
            config_data = json.load(f)

        self.assertTrue("random_seed" in config_data)
        self.assertTrue(config_data["random_seed"]["master_seed"] == 1234)
        self.assertTrue(config_data["random_seed"]["init"] == 2906597030)
        self.assertTrue(config_data["random_seed"]["place"] == 1602421836)
        self.assertTrue(config_data["random_seed"]["detect"] == 216676975)
        self.assertTrue(config_data["random_seed"]["project"] == 2698621808)
        self.assertTrue(config_data["random_seed"]["prune"] == 507409703)
        self.assertTrue(config_data["random_seed"]["input"] == 2825158027)
        self.assertTrue(config_data["random_seed"]["simulate"] == 3613074)

        self.assertTrue("regions" in config_data)
        self.assertTrue("Striatum" in config_data["regions"])

        self.assertTrue("connectivity" in config_data["regions"]["Striatum"])
        self.assertTrue("neurons" in config_data["regions"]["Striatum"])

        for neuron in config_data["regions"]["Striatum"]["neurons"]:
            self.assertTrue("neuron_path" in config_data["regions"]["Striatum"]["neurons"][neuron])
            self.assertTrue("num_neurons" in config_data["regions"]["Striatum"]["neurons"][neuron])
            self.assertTrue("neuron_type" in config_data["regions"]["Striatum"]["neurons"][neuron])
            self.assertTrue("volume_id" in config_data["regions"]["Striatum"]["neurons"][neuron])

    def test_place(self):

        # Place neurons

        position_file = os.path.join(self.sim_name, "network-neuron-positions.hdf5")

        npn = SnuddaPlace(config_file=self.config_file,
                          log_file=None,
                          verbose=True,
                          d_view=None,          # TODO: If d_view is None code run sin serial, add test parallel
                          h5libver="latest")

        npn.parse_config()
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

    def test_population_units(self, stage="place-pop-unit-random"):

        network_path = os.path.join(os.path.dirname(__file__), "networks", "network_place_pop_unit_random")
        cnc = SnuddaInit(struct_def={}, network_path=network_path, random_seed=123457)
        cnc.define_striatum(num_dSPN=1000, num_iSPN=1000, num_FS=20, num_LTS=0, num_ChIN=0,
                            volume_type="cube")
        cnc.add_population_unit_random(structure_name="Striatum", neuron_types=["dSPN", "iSPN"],
                                       fraction_of_neurons=0.5)
        cnc.add_population_unit_random(structure_name="Striatum", neuron_types=["dSPN", "iSPN"],
                                       fraction_of_neurons=0.2)
        cnc.add_population_unit_random(structure_name="Striatum", neuron_types=["dSPN"],
                                       fraction_of_neurons=0.3)
        cnc.add_population_unit_random(structure_name="Striatum", neuron_types=["iSPN"],
                                       fraction_of_neurons=0.15)
        cnc.add_population_unit_random(structure_name="Striatum", neuron_types=["iSPN"],
                                       fraction_of_neurons=0.15, unit_id=10)
        cnc.write_json()

        npn = SnuddaPlace(network_path=network_path, h5libver="latest", verbose=True)
        npn.parse_config()
        npn.write_data()

        network_file = os.path.join(network_path, "network-neuron-positions.hdf5")
        sl = SnuddaLoad(network_file)

        pop_units = sl.get_neuron_population_units()
        neuron_types = sl.get_neuron_types()

        self.assertTrue(np.abs(np.sum(pop_units == 1) - 2000*0.5) < 50)  # 50% of dSPN, iSPN should be pop unit 1
        self.assertTrue(np.abs(np.sum(pop_units == 2) - 2000*0.2) < 50)  # 20% should be pop unit 1
        self.assertTrue(np.abs(np.sum(pop_units == 3) - 1000*0.3) < 35)  # 30% of dSPN should be pop unit 1
        self.assertTrue(np.abs(np.sum(pop_units == 4) - 1000*0.15) < 35)  # 15% of iSPN should be pop unit 1
        self.assertTrue(np.abs(np.sum(pop_units == 10) - 1000*0.15) < 35)  # 15% of iSPN should be pop unit 1


if __name__ == '__main__':
    unittest.main()

# python3 -m unittest test_place
