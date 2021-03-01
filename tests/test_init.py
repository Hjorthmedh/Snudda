import unittest
import os
import numpy as np

from snudda.init import SnuddaInit


class TestInit(unittest.TestCase):

    def setUp(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

    def test_init(self):
        
        with self.subTest(stage="init_cube"):

            network_path = os.path.join(os.path.dirname(__file__), "networks", "network_testing_init_cube")
    
            config_name = os.path.join(network_path, "network-config.json")
            cnc = SnuddaInit(struct_def={}, config_file=config_name)
            cnc.define_striatum(num_dSPN=47500, num_iSPN=47500, num_FS=1300, num_LTS=0, num_ChIN=0,
                                volume_type="cube")
            cnc.write_json()

            # TODO: This only checks that the code runs. Add check for output

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
            cnc = SnuddaInit(struct_def={}, config_file=config_name, random_seed=123)
            cnc.define_striatum(num_neurons=1670000)
            cnc.write_json()

        with self.subTest(stage="population-unit-random"):
            network_path = os.path.join(os.path.dirname(__file__), "networks", "network_pop_unit_random")
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

        with self.subTest(stage="population-unit-sphere"):
            network_path = os.path.join(os.path.dirname(__file__), "networks", "network_pop_unit_sphere")
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
                                            probability_function="exp(-d/150e-6)")

            cnc.write_json()

