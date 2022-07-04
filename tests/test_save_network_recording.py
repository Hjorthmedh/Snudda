import unittest
import os
import argparse
import h5py
import neuron
import numpy as np

from snudda.simulate.save_network_recording import *


class TestSaveNetworkRecording(unittest.TestCase):

    def testSavingTimestep(self):

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

        from snudda.neuromodulation.modulation_network import Neuromodulation

        nl = Neuromodulation()
        nl.set_timestep(dt=0.025)
        nl.set_modulation(neurotransmitter="dopamine", neurotransmitter_key="DA")
        nl.transient(neurotransmitter="dopamine", method="alpha", duration=3000,
                     parameters={"tstart": 700, "gmax": 1, "tau": 300})

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

        nl.save(dir_path=os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "modulation"),
                name="dopamine_modulation.json")

        args = argparse.Namespace()
        os.environ["SNUDDA_DATA"] = os.path.join(os.path.dirname(__file__), "neuromodulation", "data")
        args.neuromodulation = os.path.join(os.path.dirname(__file__), "neuromodulation", "data", "modulation",
                                            "dopamine_modulation.json")
        args.path = os.path.join(os.path.dirname(__file__), "test_network_neuromodulation")
        args.output_file = os.path.join(os.path.dirname(__file__), "simulation", "test.hdf5")
        args.time = 0.1
        args.nrnivmodl = os.path.join(os.environ["SNUDDA_DATA"], "mechanisms")

        args.network_file = None

        args.disable_gj = False
        args.disable_synapses = False
        args.exportCoreNeuron = False
        args.input_file = None
        args.mech_dir = None
        args.network_file = None
        args.profile = False
        args.randomseed = None
        args.record_all = "0,1"
        args.record_volt = True
        args.verbose = False

        if os.path.exists("mechanisms"):
            pass
        else:
            os.system(f"ln -s {args.nrnivmodl}")
            os.system("nrnivmodl mechanisms")

        from mpi4py import MPI  # This must be imported before neuron, to run parallel
        from neuron import h  # , gui
        import neuron

        if os.path.exists("x86_64/.libs/libnrnmech.so"):
            print("!!! Manually loading libraries")
            try:
                h.nrn_load_dll("x86_64/.libs/libnrnmech.so")
            except:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)

        from snudda.core import Snudda

        s = Snudda(network_path=args.path)
        s.simulate(args=args)

        self.assertTrue(os.path.exists(args.output_file))
        self.assertTrue(os.path.exists(args.path))
        self.assertTrue(os.path.exists(os.path.join(args.path, "input-spikes.hdf5")))

        f = h5py.File(os.path.join("simulation", "test.hdf5"), "r")

        test_array = f["neurons"]["0"]["voltage"]["sec_x"][()]
        should_be = np.ones_like(test_array) * 0.5

        dt_step = f["time"][()][1] - f["time"][()][0]
        d = np.round(dt_step, 4)

        self.assertTrue((should_be == test_array).all())
        self.assertTrue(d == 0.0005)

    def testNeuronRecordings(self):

        test_vector = neuron.h.Vector([1, 2, 3, 4])
        test_recordings = NeuronRecordings(neuron_id=12345)
        test_recordings.register_compartment_data(data=test_vector,
                                                  data_type="voltage",
                                                  sec_id=0,
                                                  sec_x=0.5)
        test_recordings.register_synapse_data(data=test_vector,
                                              data_type="synaptic_current",
                                              synapse_type="tmGlut",
                                              presynaptic_id=0,
                                              sec_id=1,
                                              sec_x=0.1,
                                              cond=1e-3)
        test_recordings.register_spike_data(data=test_vector,
                                            sec_id=0,
                                            sec_x=0.5)

        self.assertTrue("voltage" in test_recordings.data)
        self.assertTrue("synaptic_current" in test_recordings.data)
        self.assertTrue("spikes" in test_recordings.data)

        self.assertTrue(isinstance(test_recordings.data["voltage"], CompartmentData))
        self.assertTrue(isinstance(test_recordings.data["synaptic_current"], SynapseData))
        self.assertTrue(isinstance(test_recordings.data["spikes"], SpikeData))
