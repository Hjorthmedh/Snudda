import unittest
import os
import argparse
from snudda.neuromodulation.modulation_synapse import NeuromodulationSynapse
from snudda.core import Snudda
import json
import numpy as np
from collections import OrderedDict


class TestNeuromodulation(unittest.TestCase):

    def setUp(self):

        args = argparse.Namespace()
        os.environ["SNUDDA_DATA"] = os.path.join(os.path.dirname(__file__), "neuromodulation", "data")
        args.action = "init"
        args.size = 10
        args.overwrite = True
        args.path = os.path.join(os.path.dirname(__file__), "test_network_neuromodulation")
        args.randomseed = 12345
        args.neurons_dir = os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "neurons")
        args.input = os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "input", "input.json")

        os.system(f"snudda init {args.path} -size {args.size} -overwrite -randomseed {args.randomseed}")

        os.system(f"snudda place {args.path}")

        os.system(f"snudda detect {args.path}")

        os.system(f"snudda prune {args.path}")

        os.system(f"cp -a {args.input} {args.path}/input.json")

        os.system(f"snudda input {args.path} --time 5")

    def test_dump(self):

        args = argparse.Namespace()
        os.environ["SNUDDA_DATA"] = os.path.join(os.path.dirname(__file__), "neuromodulation", "data")
        args.neuromodulation = os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "modulation",
                                            "dopamine_modulation.json")
        args.path = os.path.join(os.path.dirname(__file__), "test_network_neuromodulation")
        args.output_file = os.path.join(os.path.dirname(__file__), "simulation", "test.hdf5")
        args.time = 0.01
        args.nrnivmodl = os.path.join(os.environ["SNUDDA_DATA"], "neurons", "mechanisms")
        args.network_file = None

        args.disable_gj = False
        args.exportCoreNeuron = False
        args.input_file = None
        args.mech_dir = None
        args.network_file = None
        args.profile = False
        args.randomseed = None
        args.record_all = None
        args.record_volt = True
        args.verbose = False

        network_path = args.path

        if args.network_file:
            network_file = args.network_file
        else:
            network_file = os.path.join(network_path, "network-synapses.hdf5")

        if args.input_file:
            input_file = args.input_file
        else:
            input_file = os.path.join(network_path, "input-spikes.hdf5")

        if args.output_file:
            output_file = args.output_file
        else:
            output_file = os.path.join(network_path, "simulation", "output.hdf5")

        os.makedirs(os.path.join(network_path, "simulation"), exist_ok=True)

        print(f"Using input file {input_file}")

        save_dir = os.path.join(os.path.dirname(network_file), "simulation")

        if not os.path.exists(save_dir):
            print(f"Creating directory {save_dir}")
            os.makedirs(save_dir, exist_ok=True)

        print(f"args: {args}")

        disable_gj = args.disable_gj
        if disable_gj:
            print("!!! WE HAVE DISABLED GAP JUNCTIONS !!!")

        log_file = os.path.join(os.path.dirname(network_file), "log", "network-simulation-log.txt")
        log_dir = os.path.join(os.path.dirname(network_file), "log")
        if not os.path.exists(log_dir):
            print(f"Creating directory {log_dir}")
            os.makedirs(log_dir, exist_ok=True)

        with open(args.neuromodulation, "r") as f:
            neuromod_dict = json.load(f, object_pairs_hook=OrderedDict)

        from snudda.neuromodulation.neuromodulation import SnuddaSimulateNeuromodulation

        sim = SnuddaSimulateNeuromodulation(network_file=network_file,
                                            input_file=input_file,
                                            output_file=output_file,
                                            disable_gap_junctions=disable_gj,
                                            log_file=log_file,
                                            verbose=args.verbose)

        sim.setup()
        sim.add_external_input()
        sim.apply_neuromodulation(neuromod_dict)
        sim.neuromodulation_network_wide()

        sim.check_memory_status()

        sim.add_volt_recording_soma()
        record_cell_id = np.array([0, 1])
        sim.add_volt_recording_all(cell_id=record_cell_id)
        sim.add_synapse_current_recording_all(record_cell_id)

        t_sim = args.time * 1000  # Convert from s to ms for Neuron simulator

        sim.check_memory_status()
        sim.run(t_sim)  # In milliseconds
        sim.write_output()

        from snudda.utils.dump import dump

        dump(snudda_simulate_obj=sim)

    def test_modulation(self):
        pass

    def test_modulation_synapse(self):

        sw = NeuromodulationSynapse()
        sw.set_weight(weight=1e-2)

        # Acetylcholine

        sw.set_connection_type(connector="concACh", neuromodulation_key="ACh")

        sw.add_cell_modulation(neuromodulation_key="ACh",
                               cell="dSPN",
                               ion_channels={
                                   "soma": ["kir_ms", "cal12_ms", "cal13_ms", "can_ms", "Im_ms"],
                                   "dendrite": ["kir_ms", "cal12_ms", "cal13_ms"]},
                               type_connection="spiking-concentration")

        # Dopamine

        sw.set_connection_type(connector="concDA", neuromodulation_key="DA")

        sw.add_cell_modulation(neuromodulation_key="DA",
                               cell="dSPN",
                               ion_channels={
                                   "soma": ["cal12_ms", "cal13_ms", "kir_ms", "kas_ms", "naf_ms", "can_ms"],
                                   "dendrite": ["kir_ms", "kas_ms", "kaf_ms", "cal12_ms", "cal13_ms"],
                                   "axon": []},
                               receptors={"tmGabaA": {"maxMod": 0.8},
                                          "tmGlut": {"maxMod_AMPA": 1.2,
                                                     "maxMod_NMDA": 1.3,
                                                     "failRate": 0.7}},
                               extrinsic=["CorticalBase", "CorticalSignal", "Thalamic"],
                               type_connection="spiking-concentration")

        sw.save(dir_path=os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "modulation"),
                name="dopamine_acetylcholine_modulation.json")

    def test_neuromodulation_adaptive(self):
        args = argparse.Namespace()
        os.environ["SNUDDA_DATA"] = os.path.join(os.path.dirname(__file__), "neuromodulation", "data")
        args.neuromodulation = os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "modulation",
                                            "dopamine_acetylcholine_modulation.json")
        args.path = os.path.join(os.path.dirname(__file__), "test_network_neuromodulation_adaptive")
        args.output_file = os.path.join(os.path.dirname(__file__), "simulation", "test.hdf5")
        args.time = 0.1
        args.nrnivmodl = os.path.join(os.environ["SNUDDA_DATA"], "neurons", "mechanisms")
        args.network_file = None

        args.disable_gj = False
        args.exportCoreNeuron = False
        args.input_file = None
        args.mech_dir = None
        args.network_file = None
        args.profile = False
        args.randomseed = None
        args.record_all = None
        args.record_volt = True
        args.verbose = True

        s = Snudda(network_path=args.path)
        s.simulate(args=args)

    def test_neuromodulation_replay(self):
        args = argparse.Namespace()
        os.environ["SNUDDA_DATA"] = os.path.join(os.path.dirname(__file__), "neuromodulation", "data")
        args.neuromodulation = os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "modulation",
                                            "dopamine_modulation.json")
        args.path = os.path.join(os.path.dirname(__file__), "test_network_neuromodulation")
        args.output_file = os.path.join(os.path.dirname(__file__), "simulation", "test.hdf5")
        args.time = 0.1
        args.nrnivmodl = os.path.join(os.environ["SNUDDA_DATA"], "neurons", "mechanisms")
        args.network_file = None

        args.disable_gj = False
        args.exportCoreNeuron = False
        args.input_file = None
        args.mech_dir = None
        args.network_file = None
        args.profile = False
        args.randomseed = None
        args.record_all = None
        args.record_volt = True
        args.verbose = False

        if os.path.exists("mechanisms"):
            pass
        else:
            os.system(f"ln -s {args.nrnivmodl}")
            os.system("nrnivmodl mechanisms")

        from snudda.core import Snudda

        s = Snudda(network_path=args.path)
        s.simulate(args=args)

        self.assertTrue(os.path.exists(args.output_file))
        self.assertTrue(os.path.exists(args.path))
        self.assertTrue(os.path.exists(os.path.join(args.path, "input-spikes.hdf5")))


if __name__ == "__main__":
    unittest.main()
