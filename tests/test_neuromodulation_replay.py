import unittest
import os
import argparse
from snudda.neuromodulation.modulation_synapse import NeuromodulationSynapse
from snudda.core import Snudda
import json
import numpy as np
from collections import OrderedDict


class TestNeuromodulationReplay(unittest.TestCase):

    def setUp(self):

        args = argparse.Namespace()
        args.action = "init"
        args.size = 10
        args.overwrite = True
        args.path = os.path.join(os.path.dirname(__file__), "networks", "test_network_neuromodulation")
        args.randomseed = 12345

        os.system(f"snudda init {args.path} -size {args.size} -overwrite -randomseed {args.randomseed}")

        os.system(f"snudda place {args.path}")

        os.system(f"snudda detect {args.path}")

        os.system(f"snudda prune {args.path}")

        os.system(f"snudda input {args.path} --time 5")

    def test_modulation(self):
        from snudda.neuromodulation.modulation_network import Neuromodulation

        nl = Neuromodulation()
        nl.set_timestep(dt=0.025)
        nl.set_modulation(neurotransmitter="dopamine", neurotransmitter_key="DA")
        nl.transient(neurotransmitter="dopamine", method="bath_application", duration=3000,
                     parameters={"gmax": 0})

        nl.ion_channel_modulation(neurotransmitter="dopamine",
                                  cell_type="dSPN",
                                  section="soma",
                                  ion_channels=["kas_ms", "kaf_ms", "can_ms"])
        nl.ion_channel_modulation(neurotransmitter="dopamine",
                                  cell_type="dSPN",
                                  section="dendrite",
                                  ion_channels=["kas_ms", "kaf_ms"])

        nl.ion_channel_modulation(neurotransmitter="dopamine",
                                  cell_type="iSPN",
                                  section="soma",
                                  ion_channels=["kir_ms", "kas_ms", "kaf_ms", "naf_ms", "cal12_ms", "cal13_ms",
                                                "can_ms", "car_ms"])
        nl.ion_channel_modulation(neurotransmitter="dopamine",
                                  cell_type="iSPN",
                                  section="dendrite",
                                  ion_channels=["kir_ms", "kas_ms", "kaf_ms", "naf_ms", "cal12_ms", "cal13_ms",
                                                "can_ms", "car_ms"])
        nl.ion_channel_modulation(neurotransmitter="dopamine",
                                  cell_type="iSPN",
                                  section="axon",
                                  ion_channels=["kir_ms", "kas_ms", "kaf_ms", "naf_ms", "cal12_ms", "cal13_ms",
                                                "can_ms", "car_ms"])

        nl.ion_channel_modulation(neurotransmitter="dopamine",
                                  cell_type="FSN",
                                  section="soma",
                                  ion_channels=["kir_fs", "kas_fs", "kaf_fs", "naf_fs"])
        nl.ion_channel_modulation(neurotransmitter="dopamine",
                                  cell_type="FSN",
                                  section="dendrite",
                                  ion_channels=["kir_fs"])
        
        nl.save(dir_path=os.path.join(os.path.dirname(__file__), "networks", "test_network_neuromodulation"), name="dopamine_modulation.json")

    def test_neuromodulation_replay_receptors(self):
        args = argparse.Namespace()
        args.neuromodulation = os.path.join(os.path.dirname(__file__), "networks", "test_network_neuromodulation",
                                            "dopamine_modulation.json")
        args.path = os.path.join(os.path.dirname(__file__), "networks", "test_network_neuromodulation")
        args.output_file = os.path.join(os.path.dirname(__file__), "simulation", "test.hdf5")
        args.time = 2.0
        args.nrnivmodl = os.path.join(os.path.dirname(__file__), os.path.pardir,
                                    "snudda", "data", "neurons", "mechanisms")
        args.network_file = None

        args.disable_gj = False
        args.exportCoreNeuron = False
        args.input_file = None
        args.mech_dir = None
        args.network_file = None
        args.disable_synapses = False
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
