import unittest

from snudda.simulate import SnuddaSimulate
from snudda import Snudda


class NeuromodulationTestCase(unittest.TestCase):

    def setUp(self):

        self.snudda = Snudda(network_path="networks/network_rxd")
        self.snudda.init_tiny(neuron_paths="validation/ballanddoublestick_rxd",
                              neuron_names="neuron_1",
                              number_of_neurons=[10])
        self.snudda.create_network()

        # Check why file size is so large, and why it is so slow to generate!

        self.snudda.simulate()

    def test_something(self):
        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
