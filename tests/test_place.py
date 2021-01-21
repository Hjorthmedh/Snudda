import unittest
import os
import sys
import json

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

class TestPlace(unittest.TestCase):

    def setUp(self):
        from snudda.init import SnuddaInit
        cell_spec = os.path.join(os.path.dirname(__file__), "..", "snudda", "data", "validation")
        self.sim_name = os.path.join("tests", "network_testing_place")
        self.config_name = os.path.join(self.sim_name, "network-config.json")
        cnc = SnuddaInit(struct_def={}, config_name=self.config_name, num_population_units=1, random_seed=1234)
        cnc.define_striatum(num_dSPN=10, num_iSPN=0, num_FS=10, num_LTS=0, num_ChIN=0,
                            volume_type="cube", cell_spec_dir=cell_spec)
        cnc.write_json(self.config_name)

    def tearDown(self):
        pass

    def test_init(self):
        with open(self.config_name, "r") as f:
            config_data = json.load(f)

        self.assertTrue("RandomSeed" in config_data)
        self.assertTrue(config_data["RandomSeed"]["masterseed"] == 1234)
        self.assertTrue(config_data["RandomSeed"]["init"] == 2906597030)
        self.assertTrue(config_data["RandomSeed"]["place"] == 1602421836)
        self.assertTrue(config_data["RandomSeed"]["detect"] == 216676975)
        self.assertTrue(config_data["RandomSeed"]["prune"] == 2698621808)
        self.assertTrue(config_data["RandomSeed"]["input"] == 507409703)
        self.assertTrue(config_data["RandomSeed"]["simulate"] == 2825158027)


def test_place(self):


        pass


if __name__ == '__main__':
    unittest.main()

# python3 -m unittest test_place