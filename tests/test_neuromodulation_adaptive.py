import unittest
import os
import argparse
from snudda.neuromodulation.modulation_synapse import NeuromodulationSynapse
from snudda.core import Snudda
import json
import numpy as np
from collections import OrderedDict
import os


class TestNeuromodulationAdaptive(unittest.TestCase):

    def setUp(self):
        import os
        args = argparse.Namespace()
        os.environ["SNUDDA_DATA"] = os.path.join(os.path.dirname(__file__), "neuromodulation", "data")
        args.action = "init"
        args.size = 10
        args.overwrite = True
        args.path = os.path.join(os.path.dirname(__file__), "test_network_neuromodulation_adaptive")
        args.randomseed = 12345
        args.neurons_dir = os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "neurons")
        args.input = os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "input", "input.json")

        from neuromodulation.neuromodulationInitSNc import neuromodulationInit

        config_name = os.path.join(args.path, "network-config.json")
        cnc = neuromodulationInit(config_file=config_name, random_seed=12345)

        cnc.define_striatum_neuromodulation(num_dSPN=25, volume_type="cube", neurons_dir=args.neurons_dir)
        cnc.define_snc(nNeurons=4, neuron_dir=args.neurons_dir)
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
                                   "soma": ["kas_ms", "kaf_ms", "can_ms"],
                                   "dendrite": ["kaf_ms", "kas_ms"],
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
        args.nrnivmodl = os.path.join(os.environ["SNUDDA_DATA"], "neurons", "mechanisms_ptr", "da_ch")
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


if __name__ == "__main__":
    unittest.main()