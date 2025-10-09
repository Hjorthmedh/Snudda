import unittest
import os
import numpy as np
import multiprocessing

try:
    multiprocessing.set_start_method("spawn", force=True)
except RuntimeError:
    # Already set, fine
    pass

from snudda import Snudda
from snudda.utils import SnuddaLoadSimulation

# TODO: Write example that uses Anu's SBML Dopamine cascade
# https://www.ebi.ac.uk/biomodels/BIOMD0000000636#Files
# Skriv SBML --> vårt json format, converter
# https://modeldb.science/237653?tab=2&file=Beskow/speedy_reduced2.mod
#
# Use: https://github.com/sbmlteam/libsbml

def run_sim(network_path, output_file):
    print(f"Running simulation: {output_file}")

    import os
    print("CWD in run_sim:", os.getcwd())

    mech_dir = "../validation/mechanisms"  # Added the kirrxd and DASyn as symbolic links to mechanisms

    # os.system(f"nrnivmodl {mech_dir}")

    snudda = Snudda(network_path=network_path)
    # snudda.compile_mechanisms(mech_dir=mech_dir)

    from neuron import h, nrn

    so_path = os.path.join("aarch64", ".libs", "libnrnmech.so")
    h.nrn_load_dll(so_path)

    print("Mechanisms loaded:")
    print("kirrxd visible?", hasattr(h, "kirrxd"))

    sim = snudda.simulate(time=0, mech_dir=mech_dir, use_rxd_neuromodulation=True)

    n = sim.neurons[0]

    sim.add_rxd_concentration_recording(species="DA", neuron_id=0,
                                        region="soma_internal",
                                        sec_id=-1,
                                        sec_x=0.5)

    sim.add_rxd_concentration_recording(species="B", neuron_id=0,
                                        region="soma_internal",
                                        sec_id=-1,
                                        sec_x=0.5)

    sim.add_rxd_concentration_recording(species="PKA", neuron_id=0,
                                        region="soma_internal",
                                        sec_id=-1,
                                        sec_x=0.5)

    sim.add_density_mechanism_recording(density_mechanism="kirrxd",
                                        variable="modulation_factor",
                                        neuron_id=0,
                                        sec_id=-1,
                                        sec_x=0.5)

    sim.add_density_mechanism_recording(density_mechanism="kirrxd",
                                        variable="m",
                                        neuron_id=0,
                                        sec_id=-1,
                                        sec_x=0.5)

    sim.add_membrane_recording(variable="PKAi",
                               neuron_id=0,
                               sec_id=-1,
                               sec_x=0.5)

    # Add DA synapse
    # da_syn = h.DASyn(soma(0.5))
    mod_file = "DASyn"
    eval_str = f"sim.sim.neuron.h.{mod_file}"
    channel_module = eval(eval_str)

    da_syn = sim.get_external_input_synapse(channel_module=channel_module,
                                            section=sim.neurons[0].icell.soma[0],
                                            section_x=0.5)
    da_syn.tau = 1

    net_stim = sim.sim.neuron.h.NetStim()
    net_stim.number = 100
    net_stim.start = 300
    net_stim.interval = 5

    nc = sim.sim.neuron.h.NetCon(net_stim, da_syn)
    nc.weight[0] = 100_000_000.0  # units : molecules/ms

    sim.neurons[0].modulation.link_synapse(species_name="DA",
                                           region="soma_internal",
                                           synapse=da_syn,
                                           flux_variable="open")

    sim.run(t=1000)

    # output_file = os.path.join(self.network_path, "simulation", "output-2.hdf5")
    sim.record.set_new_output_file(output_file)
    sim.record.write()

    # Important we need to delete the old da_syn
    del da_syn
    del net_stim
    del nc

    sim.clear_neuron()


class NeuromodulationTestCase(unittest.TestCase):

    def setUp(self):

        test_path = "test_project2"
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

        mech_path = os.path.join(os.path.dirname(__file__), "validation/mechanisms")
        os.system(f"nrnivmodl {mech_path}")


    def test_reaction(self):

        output_file = os.path.join(self.network_path, "simulation", "output-2.hdf5")

        p = multiprocessing.Process(target=run_sim, args=(self.network_path, output_file,))
        p.start()
        p.join()

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
        self.assertTrue(np.max(np.abs(data_b[0][0][0] - 0.7)) < 1e-7)
        self.assertTrue(np.max(np.abs(data_ab[0][0][0] - 0.1)) < 1e-7)

        # self.assertTrue(data_a[0][0][-1] < data_a[0][0][0])
        # self.assertTrue(data_b[0][0][-1] < data_b[0][0][0])
        # self.assertTrue(data_ab[0][0][-1] > data_ab[0][0][0])

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
        os.chdir("..")

        return

        print("\n[DEBUG] Running tearDown() for", self.id())

        # Remember to clear old neuron, for next unit test!
        self.sim.clear_neuron()

        # Really trying to make sure noting remains in NEURON
        # --- Hard NEURON reset to avoid stale rxd or mod pointer references ---
        from neuron import h
        import gc

        # Clear all rxd species, regions, and reactions
        try:
            from neuron import rxd
            rxd.region._all_regions.clear()
            rxd.species._all_species.clear()
            rxd.reaction._all_reactions.clear()
            rxd.rxdmath._all_reactions.clear()
            rxd.node._all_nodes.clear()
            rxd.include_flux._all_include_fluxes.clear()
        except Exception:
            pass

        # Reset ParallelContext and clear gids
        try:
            pc = h.ParallelContext()
            pc.gid_clear()
            del pc
        except Exception:
            pass

        # Delete any remaining sections and point processes
        for sec in list(h.allsec()):
            h.delete_section(sec=sec)
        h('forall delete_section()')

        # Force NEURON to release internal handles
        h('objref nil')
        h('forall delete_section()')
        h('nrnpython("import gc; gc.collect()")')

        gc.collect()

        import pdb
        pdb.set_trace()


if __name__ == '__main__':
    unittest.main()
