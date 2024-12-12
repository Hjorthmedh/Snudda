import unittest
import os
import numpy as np
from snudda.simulate import SnuddaSimulate
from snudda import Snudda
from snudda.utils import SnuddaLoadSimulation

# TODO: Write example that uses Anu's SBML Dopamine cascade
# https://www.ebi.ac.uk/biomodels/BIOMD0000000636#Files
# Skriv SBML --> vårt json format, converter
# https://modeldb.science/237653?tab=2&file=Beskow/speedy_reduced2.mod
#
# Use: https://github.com/sbmlteam/libsbml


class NeuromodulationTestCase(unittest.TestCase):

    def setUp(self):

        test_path = "test_project"
        if os.path.isdir(test_path):
            import shutil
            shutil.rmtree(test_path)
        os.mkdir(test_path)
        os.chdir(test_path)

        self.neuron_path = "../validation/dspn_neurons_rxd"
        self.network_path = "networks/network_rxd"

        self.snudda = Snudda(network_path=self.network_path)
        self.snudda.init_tiny(neuron_paths=self.neuron_path,
                              neuron_names="neuron",
                              number_of_neurons=[10],
                              random_seed=123456)
        self.snudda.create_network()

        # Check why file size is so large, and why it is so slow to generate!

        # mech_dir = "../validation/mechanisms_rxd"
        mech_dir = "../validation/mechanisms"   # Added the kirrxd and DASyn as symbolic links to mechanisms

        # self.snudda.compile_mechanisms(mech_dir=mech_dir)
        self.sim = self.snudda.simulate(time=0, mech_dir=mech_dir, use_rxd_neuromodulation=True)

    def test_reaction(self):

        n = self.sim.neurons[0]

        self.sim.add_rxd_concentration_recording(species="DA", neuron_id=0,
                                                 region="soma_internal",
                                                 sec_id=-1,
                                                 sec_x=0.5)

        self.sim.add_rxd_concentration_recording(species="B", neuron_id=0,
                                                 region="soma_internal",
                                                 sec_id=-1,
                                                 sec_x=0.5)

        self.sim.add_rxd_concentration_recording(species="PKA", neuron_id=0,
                                                 region="soma_internal",
                                                 sec_id=-1,
                                                 sec_x=0.5)

        self.sim.add_density_mechanism_recording(density_mechanism="kirrxd",
                                                 variable="modulation_factor",
                                                 neuron_id=0,
                                                 sec_id=-1,
                                                 sec_x=0.5)

        self.sim.add_density_mechanism_recording(density_mechanism="kirrxd",
                                                 variable="m",
                                                 neuron_id=0,
                                                 sec_id=-1,
                                                 sec_x=0.5)

        self.sim.add_membrane_recording(variable="PKAi",
                                        neuron_id=0,
                                        sec_id=-1,
                                        sec_x=0.5)

        # Add DA synapse

        #da_syn = h.DASyn(soma(0.5))
        mod_file = "DASyn"
        eval_str = f"self.sim.sim.neuron.h.{mod_file}"
        channel_module = eval(eval_str)

        da_syn = self.sim.get_external_input_synapse(channel_module=channel_module,
                                                     section=self.sim.neurons[0].icell.soma[0],
                                                     section_x=0.5)
        da_syn.tau = 1

        net_stim = self.sim.sim.neuron.h.NetStim()
        net_stim.number = 100
        net_stim.start = 300
        net_stim.interval = 5

        nc = self.sim.sim.neuron.h.NetCon(net_stim, da_syn)
        nc.weight[0] = 100_000_000.0   #units : molecules/ms

        self.sim.neurons[0].modulation.link_synapse(species_name="DA",
                                                    region="soma_internal",
                                                    synapse=da_syn,
                                                    flux_variable="open")

        self.sim.run(t=1000)

        output_file = os.path.join(self.network_path, "simulation", "output-2.hdf5")
        self.sim.record.set_new_output_file(output_file)
        self.sim.record.write()

        nd = SnuddaLoadSimulation(output_file)
        time = nd.get_time()
        data_a = nd.get_data("DA", 0)
        data_b = nd.get_data("B", 0)
        data_ab = nd.get_data("PKA", 0)

        data_kir_modulation_factor = nd.get_data("kirrxd.modulation_factor", 0)[0][0]
        data_kir_m = nd.get_data("kirrxd.m", 0)[0][0]
        data_voltage = nd.get_data("voltage", 0)[0][0]
        data_pka = nd.get_data("membrane.PKAi", 0)[0][0]

        self.assertTrue(np.max(np.abs(data_a[0][0][0] - 0)) < 1e-7)
        self.assertTrue(np.max(np.abs(data_b[0][0][0] - 0.7e-3)) < 1e-7)
        self.assertTrue(np.max(np.abs(data_ab[0][0][0] - 0.1e-3)) < 1e-7)

        #self.assertTrue(data_a[0][0][-1] < data_a[0][0][0])
        #self.assertTrue(data_b[0][0][-1] < data_b[0][0][0])
        #self.assertTrue(data_ab[0][0][-1] > data_ab[0][0][0])

        # Plot A, B, PKA activity
        da = data_a[0][0]
        db = data_b[0][0]
        dab = data_ab[0][0]

        self.plot_data(time, np.hstack([da, db, dab]), legend=["DA", "B", "PKA"],
                       ylabel="Concentration", filename="concentration.png")

        # Plot the voltage and KIR modulation
        self.plot_data(time, np.hstack([data_kir_modulation_factor,
                                        data_kir_m, data_pka]),
                       legend=["kir_mod_factor", "kir_m", "membrane.pka"],
                       filename="kir_activation.png")

        self.plot_data(time, data_voltage,
                       legend=["voltage"],
                       filename="voltage.png")

    def plot_data(self, time, data, legend, filename=None, ylabel=None):
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(time, data, label=legend)
        plt.xlabel("Time (s)")
        plt.ylabel(ylabel)
        plt.legend()

        if filename is not None:
            plt.savefig(filename, dpi=300)

        plt.ion()
        plt.show()

    def plot_kir_data(self):
        pass

    def tearDown(self):
        # Remember to clear old neuron, for next unit test!
        self.sim.clear_neuron()
        os.chdir("..")


if __name__ == '__main__':
    unittest.main()
