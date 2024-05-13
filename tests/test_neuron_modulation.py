import unittest
import os
from snudda.simulate import SnuddaSimulate
from snudda import Snudda


class NeuromodulationTestCase(unittest.TestCase):

    def setUp(self):

        self.neuron_path = "validation/ballanddoublestick_rxd"

        self.snudda = Snudda(network_path="networks/network_rxd")
        self.snudda.init_tiny(neuron_paths=self.neuron_path,
                              neuron_names="neuron_1",
                              number_of_neurons=[10])
        self.snudda.create_network()

        # Check why file size is so large, and why it is so slow to generate!

        self.sim = self.snudda.simulate()

    def test_something(self):

        n = self.sim.neurons[0]
        from snudda.neurons.neuron_modulation import NeuronModulation

        nm = NeuronModulation(neuron=n)
        config = os.path.join(self.neuron_path, "reaction_diffusion.json")
        nm.load_json(config_path=config)

        import pdb
        pdb.set_trace()

        self.assertEqual(True, True)  # add assertion here

    def tearDown(self):
        # Remember to clear old neuron, for next unit test!
        self.sim.clear_neuron()


if __name__ == '__main__':
    unittest.main()
