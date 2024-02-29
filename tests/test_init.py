import unittest
import os
import json

from snudda.init.init import SnuddaInit


class TestInit(unittest.TestCase):

    def setUp(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

    def test_init(self):
        
        with self.subTest(stage="init_cube"):

            network_path = os.path.join(os.path.dirname(__file__), "networks", "network_testing_init_cube")
    
            config_name = os.path.join(network_path, "network-config.json")
            snudda_data = os.path.join(os.path.dirname(__file__), "..", "snudda", "data")
            cnc = SnuddaInit(struct_def={}, config_file=config_name, snudda_data=snudda_data)
            cnc.define_striatum(num_dSPN=47500, num_iSPN=47500, num_FS=1300, num_LTS=0, num_ChIN=0,
                                volume_type="cube")
            cnc.write_json()

            with open(config_name, "r") as f:
                cfg_data = json.load(f)

            self.assertTrue("regions" in cfg_data)
            self.assertTrue("Striatum" in cfg_data["regions"])
            self.assertTrue("neurons" in cfg_data["regions"]["Striatum"])

            self.assertTrue("dSPN" in cfg_data["regions"]["Striatum"]["neurons"])
            self.assertEqual(cfg_data["regions"]["Striatum"]["neurons"]["dSPN"]["num_neurons"], 47500)

            self.assertTrue("iSPN" in cfg_data["regions"]["Striatum"]["neurons"])
            self.assertEqual(cfg_data["regions"]["Striatum"]["neurons"]["iSPN"]["num_neurons"], 47500)

            self.assertTrue("FS" in cfg_data["regions"]["Striatum"]["neurons"])
            self.assertEqual(cfg_data["regions"]["Striatum"]["neurons"]["FS"]["num_neurons"], 1300)

            self.assertTrue("LTS" not in cfg_data["regions"]["Striatum"]["neurons"])
            self.assertTrue("ChIN" not in cfg_data["regions"]["Striatum"]["neurons"])

            self.assertTrue(os.path.isfile(cfg_data["regions"]["Striatum"]["volume"]["mesh_file"]))
            self.assertEqual(cfg_data["regions"]["Striatum"]["volume"]["d_min"], 1.5e-5)

            self.assertTrue("connectivity" in cfg_data["regions"]["Striatum"])
            self.assertTrue("FS,FS" in cfg_data["regions"]["Striatum"]["connectivity"])
            self.assertTrue("FS,dSPN" in cfg_data["regions"]["Striatum"]["connectivity"])
            self.assertTrue("FS,iSPN" in cfg_data["regions"]["Striatum"]["connectivity"])
            self.assertTrue("dSPN,dSPN" in cfg_data["regions"]["Striatum"]["connectivity"])
            self.assertTrue("iSPN,iSPN" in cfg_data["regions"]["Striatum"]["connectivity"])
            self.assertTrue("iSPN,dSPN" in cfg_data["regions"]["Striatum"]["connectivity"])
            self.assertTrue("dSPN,iSPN" in cfg_data["regions"]["Striatum"]["connectivity"])

            self.assertTrue("GABA" in cfg_data["regions"]["Striatum"]["connectivity"]["FS,FS"])
            self.assertTrue("gap_junction" in cfg_data["regions"]["Striatum"]["connectivity"]["FS,FS"])

            self.assertTrue("conductance" in cfg_data["regions"]["Striatum"]["connectivity"]["FS,FS"]["GABA"])
            self.assertTrue("channel_parameters" in cfg_data["regions"]["Striatum"]["connectivity"]["FS,FS"]["GABA"])
            self.assertTrue("pruning" in cfg_data["regions"]["Striatum"]["connectivity"]["FS,FS"]["GABA"])

        with self.subTest(stage="init_slice"):
            network_path = os.path.join(os.path.dirname(__file__), "networks", "network_testing_init_slice")

            config_name = os.path.join(network_path)
            cnc = SnuddaInit(struct_def={}, network_path=network_path)
            cnc.define_striatum(num_dSPN=47500, num_iSPN=47500, num_FS=1300, num_LTS=0, num_ChIN=0,
                                volume_type="slice")
            cnc.write_json()
            
        with self.subTest(stage="init_full"):
            network_path = os.path.join(os.path.dirname(__file__), "networks", "network_testing_init_full")

            config_name = os.path.join(network_path, "network-config.json")
            cnc = SnuddaInit(struct_def={}, config_file=config_name, random_seed=1237)
            cnc.define_striatum(num_neurons=1670000)
            cnc.write_json()

            with open(config_name, "r") as f:
                cfg_data = json.load(f)

            self.assertEqual(cfg_data["random_seed"]["master_seed"], 1237)

        with self.subTest(stage="population-unit-random"):
            network_path = os.path.join(os.path.dirname(__file__), "networks", "network_testing_init_pop_unit_random")
            cnc = SnuddaInit(struct_def={}, network_path=network_path)
            cnc.define_striatum(num_dSPN=1000, num_iSPN=1000, num_FS=20, num_LTS=0, num_ChIN=0,
                                volume_type="slice")
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

            config_name = os.path.join(network_path, "network-config.json")
            with open(config_name, "r") as f:
                cfg_data = json.load(f)

            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["unit_id"] == [1, 2, 3, 4, 10])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["method"] == "random")
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["fraction_of_neurons"] == [0.5, 0.2, 0.3, 0.15, 0.15])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["unit_id"] == [1, 2, 3, 4, 10])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["neuron_types"][0] == ["dSPN", "iSPN"])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["neuron_types"][1] == ["dSPN", "iSPN"])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["neuron_types"][2] == ["dSPN"])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["neuron_types"][3] == ["iSPN"])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["neuron_types"][4] == ["iSPN"])




        with self.subTest(stage="population-unit-sphere"):
            network_path = os.path.join(os.path.dirname(__file__), "networks", "network_testing_init_pop_unit_sphere")
            cnc = SnuddaInit(struct_def={}, network_path=network_path)
            cnc.define_striatum(num_dSPN=10000, num_iSPN=10000, num_FS=200, num_LTS=0, num_ChIN=0,
                                volume_type="slice")

            cnc.add_population_unit_density(structure_name="Striatum",
                                            neuron_types=["dSPN", "iSPN"],
                                            unit_centre=[0, 0, 0],
                                            probability_function="(d < 100e-6)*1")

            cnc.add_population_unit_density(structure_name="Striatum",
                                            neuron_types=["dSPN", "iSPN"],
                                            unit_centre=[300e-6, 0, 0],
                                            probability_function="exp(-d/200e-6)")

            cnc.add_population_unit_density(structure_name="Striatum",
                                            neuron_types=["dSPN"],
                                            unit_centre=[0, 300e-6, 0],
                                            probability_function="exp(-d/150e-6)",
                                            num_neurons=100)

            cnc.write_json()

            config_name = os.path.join(network_path, "network-config.json")
            with open(config_name, "r") as f:
                cfg_data = json.load(f)

            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["unit_id"] == [1, 2, 3])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["method"] == "radial_density")
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["centres"][0] == [0, 0, 0])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["centres"][1] == [0.000300, 0, 0])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["centres"][2] == [0, 0.000300, 0])

            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["neuron_types"][0] == ["dSPN", "iSPN"])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["neuron_types"][1] == ["dSPN", "iSPN"])
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["neuron_types"][2] == ["dSPN"])

            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["probability_functions"][0] == "(d < 100e-6)*1")
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["probability_functions"][1] == "exp(-d/200e-6)")
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["probability_functions"][2] == "exp(-d/150e-6)")

            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["num_neurons"][0] is None)
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["num_neurons"][1] is None)
            self.assertTrue(cfg_data["regions"]["Striatum"]["population_units"]["num_neurons"][2] == 100)

