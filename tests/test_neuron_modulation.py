import unittest
import os
import numpy as np
from snudda.simulate import SnuddaSimulate
from snudda import Snudda
from snudda.utils import SnuddaLoadNetworkSimulation


class NeuromodulationTestCase(unittest.TestCase):

    def setUp(self):

        test_path = "test_project"
        if os.path.isdir(test_path):
            import shutil
            shutil.rmtree(test_path)
        os.mkdir(test_path)
        os.chdir(test_path)

        self.neuron_path = "../validation/dspn_rxd"
        self.network_path = "networks/network_rxd"

        self.snudda = Snudda(network_path=self.network_path)
        self.snudda.init_tiny(neuron_paths=self.neuron_path,
                              neuron_names="neuron_1",
                              number_of_neurons=[10])
        self.snudda.create_network()

        # Check why file size is so large, and why it is so slow to generate!

        mech_dir = "../validation/mechanisms_rxd"
        # self.snudda.compile_mechanisms(mech_dir=mech_dir)
        self.sim = self.snudda.simulate(time=0, mech_dir=mech_dir)

    def test_reaction(self):

        n = self.sim.neurons[0]

        self.sim.add_rxd_concentration_recording(species="A", neuron_id=0,
                                                 region="soma_internal",
                                                 sec_type="soma",
                                                 sec_id=0,
                                                 sec_x=0.5)

        self.sim.add_rxd_concentration_recording(species="B", neuron_id=0,
                                                 region="soma_internal",
                                                 sec_type="soma",
                                                 sec_id=0,
                                                 sec_x=0.5)

        self.sim.add_rxd_concentration_recording(species="AB", neuron_id=0,
                                                 region="soma_internal",
                                                 sec_type="soma",
                                                 sec_id=0,
                                                 sec_x=0.5)

        self.sim.run(t=1000)

        output_file = os.path.join(self.network_path, "simulation", "output-2.hdf5")
        self.sim.record.set_new_output_file(output_file)
        self.sim.record.write()

        nd = SnuddaLoadNetworkSimulation(output_file)
        time = nd.get_time()
        data_a = nd.get_data("A", 0)
        data_b = nd.get_data("B", 0)
        data_ab = nd.get_data("AB", 0)

        self.assertAlmostEqual(data_a[0][0][0], 0.5, 10)
        self.assertAlmostEqual(data_b[0][0][0], 0.7, 10)
        self.assertAlmostEqual(data_ab[0][0][0], 0.1, 10)

        self.assertTrue(data_a[0][0][-1] < 0.5)
        self.assertTrue(data_b[0][0][-1] < 0.7)
        self.assertTrue(data_ab[0][0][-1] > 0.1)

        da = data_a[0][0]
        db = data_b[0][0]
        dab = data_ab[0][0]

        self.plot_ractants(time, np.hstack([da, db, dab]), legend=["A", "B", "AB"])

    def plot_ractants(self, time, data, legend, filename=None):
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(time, data, label=legend)
        plt.xlabel("Time (s)")
        plt.ylabel("Concentration")
        plt.legend()

        if filename is None:
            filename = "concentration.pdf"

        plt.savefig(filename, dpi=300)
        plt.ion()
        plt.show()



    def tearDown(self):
        # Remember to clear old neuron, for next unit test!
        self.sim.clear_neuron()
        os.chdir("..")


if __name__ == '__main__':
    unittest.main()
