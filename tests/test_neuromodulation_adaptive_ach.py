import unittest
import os
import argparse
from snudda.neuromodulation.modulation_synapse import NeuromodulationSynapse
from snudda.core import Snudda
import json
import numpy as np
from collections import OrderedDict
import os


class TestNeuromodulationAdaptiveDA(unittest.TestCase):

    def setUp(self):
        import os
        args = argparse.Namespace()
        os.environ["SNUDDA_DATA"] = os.path.join(os.path.dirname(__file__), "neuromodulation", "data")
        args.action = "init"
        args.size = 10
        args.overwrite = True
        args.path = os.path.join(os.path.dirname(__file__), "networks", "test_network_neuromodulation_adaptive_ach")
        args.randomseed = 12345
        args.neurons_dir = os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "neurons")
        args.input = os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "input", "input.json")

        from neuromodulation.neuromodulationInitSNc import neuromodulationInit

        config_name = os.path.join(args.path, "network-config.json")
        cnc = neuromodulationInit(config_file=config_name, random_seed=12345)

        cnc.define_striatum_neuromodulation(num_dSPN=2, num_iSPN=2, num_ChIN=5, volume_type="cube", neurons_dir=args.neurons_dir)
        dirName = os.path.dirname(config_name)

        if not os.path.exists(dirName):
            os.makedirs(dirName)

        cnc.write_json(config_name)

        os.system(f"snudda place {args.path}")

        os.system(f"snudda detect {args.path}")

        os.system(f"snudda prune {args.path}")

        os.system(f"cp -a {args.input} {args.path}/input.json")

        os.system(f"snudda input {args.path} --time 5")

    def test_modulation_synapse(self):

        sw = NeuromodulationSynapse()
        sw.set_weight(weight=1e-2)

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

        sw.add_cell_modulation(neuromodulation_key="ACh",
                               cell="iSPN",
                               ion_channels={
                                   "soma": ["kir_ms", "cal12_ms", "cal13_ms", "can_ms"],
                                   "dendrite": ["kir_ms", "cal12_ms", "cal13_ms"]},
                               type_connection="spiking-concentration")

        sw.save(dir_path=os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "modulation"),
                name="acetylcholine_modulation.json")

    def test_neuromodulation_adaptive_ach(self):
        args = argparse.Namespace()
        os.environ["SNUDDA_DATA"] = os.path.join(os.path.dirname(__file__), "neuromodulation", "data")
        args.neuromodulation = os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "modulation",
                                            "acetylcholine_modulation.json")
        args.path = os.path.join(os.path.dirname(__file__), "networks", "test_network_neuromodulation_adaptive_ach")
        args.output_file = os.path.join(os.path.dirname(__file__), "simulation", "test.hdf5")
        args.time = 0.01
        args.nrnivmodl = os.path.join(os.environ["SNUDDA_DATA"], "mechanisms-ptr", "ach")

        if os.path.exists("ach"):
            pass
        else:
            os.system(f"ln -s {args.nrnivmodl}")
            os.system("nrnivmodl ach")

        args.network_file = None

        args.disable_gj = False
        args.disable_synapses = False
        args.exportCoreNeuron = False
        args.input_file = None
        args.mech_dir = os.path.join(os.environ["SNUDDA_DATA"], "mechanisms-ptr", "ach")
        args.network_file = None
        args.profile = False
        args.randomseed = None
        args.record_all = None
        args.record_volt = True
        args.verbose = True

        s = Snudda(network_path=args.path)
        s.simulate(args=args)

        self.assertTrue(os.path.exists(args.output_file))
        self.assertTrue(os.path.exists(args.path))
        self.assertTrue(os.path.exists(os.path.join(args.path, "input-spikes.hdf5")))

if __name__ == "__main__":

    unittest.main()