import multiprocessing
import os
import unittest
import json
import numpy as np
from snudda.utils.load import SnuddaLoad
from snudda.utils.load_network_simulation import SnuddaLoadSimulation

import traceback
import io
from contextlib import redirect_stdout, redirect_stderr

try:
    multiprocessing.set_start_method("spawn", force=True)
except RuntimeError:
    # Already set, fine
    pass

def run_pair_recording(network_path, log_file=None):
    from neuron import h
    from snudda.simulate.pair_recording import PairRecording
    experiment_config_file = os.path.join(network_path, "experiment.json")
    network_file = os.path.join(network_path, "network-synapses.hdf5")

    pr = PairRecording(network_path=network_path,
                       network_file=network_file,
                       log_file=log_file,
                       experiment_config_file=experiment_config_file)

    # Explicitly load mechanisms compiled by nrnivmodl
    # pr.load_mechanisms()

    pr.run()
    if pr.log_file is not None:
        pr.log_file.close()

def run_pair_recording_safe(network_path, q):

    log_file = "pair_recording_test_log.txt"

    f = io.StringIO()
    try:
        with redirect_stdout(f), redirect_stderr(f):
            run_pair_recording(network_path=network_path, log_file=log_file)
        q.put(("ok", f.getvalue()))
    except Exception:
        q.put(("error", f.getvalue() + "\n" + traceback.format_exc()))
    finally:
        f.close()

    if os.path.exists(log_file):
        print("\n--- PairRecording log file contents ---")
        with open(log_file) as f:
            print(f.read())
        print("--- End of log ---\n")
    else:
        print(f"Log file not found: {log_file}")

class PairRecordingTestCase(unittest.TestCase):

    def setUp(self):

        if os.path.dirname(__file__):
            os.chdir(os.path.dirname(__file__))

        # Temporarily disable the creation of network while testing...
        # print("test_pair_recoring.py -- Network creation disabled -- we need to retune or verify the if_info.json files before running")
        # return

        self.network_path = os.path.join("networks", "pair_recording_test")
        rc = None

        if not os.path.isdir("x86_64") and not os.path.isdir("aarch64"):
            print(f"Mech files not compiled --- for some reason {os.getcwd()}, compiling!")
            mech_path = os.path.join(os.path.dirname(__file__), "validation/mechanisms")
            os.system(f"nrnivmodl {mech_path}")

        from snudda import SnuddaPlace
        from snudda import SnuddaDetect
        from snudda import SnuddaPrune

        sp = SnuddaPlace(network_path=self.network_path, rc=rc)
        sp.place()
        sp.close_log_file()

        sd = SnuddaDetect(network_path=self.network_path, rc=rc)
        sd.detect()
        sd.close_log_file()

        sp = SnuddaPrune(network_path=self.network_path, rc=rc)
        sp.prune()
        sp.close_log_file()

        self.experiment_config_file = os.path.join(self.network_path, "experiment.json")
        network_file = os.path.join(self.network_path, "network-synapses.hdf5")

        q = multiprocessing.Queue()
        process = multiprocessing.Process(
            target=run_pair_recording_safe,
            args=(self.network_path, q)
        )
        process.start()
        process.join(timeout=300)

        if process.is_alive():
            process.terminate()
            raise RuntimeError("PairRecording timed out after 120s")

        if not q.empty():
            status, output = q.get()
            print("Child process output:\n", output)
            if status == "error":
                raise RuntimeError("PairRecording failed:\n" + output)
        else:
            raise RuntimeError("PairRecording failed: no output received from process")

        # self.pr = PairRecording(network_path=self.network_path, network_file=network_file,
        #                         experiment_config_file=self.experiment_config_file)
        # self.pr.run()


    def test_frequency(self):

        print("test_pair_recording.py -- This UNIT test is not currently active.")

        # TODO: TEMP line, this is normally done in setUp... remove later
        self.network_path = os.path.join("networks", "pair_recording_test")
        self.experiment_config_file = os.path.join(self.network_path, "experiment.json")

        sim_file = os.path.join(self.network_path, "simulation", "pair-recording-simulation.hdf5")

        sl = SnuddaLoad(self.network_path)
        sns = SnuddaLoadSimulation(network_path=self.network_path,
                                   network_simulation_output_file=sim_file)

        with open(self.experiment_config_file, "r") as f:
            experiment_config = json.load(f)

        for inj_info in experiment_config["current_injection"]:

            neuron_id = inj_info["neuron_id"]
            start = inj_info["start"]
            end = inj_info["end"]
            requested_freq = inj_info["requested_frequency"]

            with self.subTest(f"Testing neuron {neuron_id}"):

                spikes = sns.get_spikes(neuron_id)
                for s, e in zip(start, end):
                    n_spikes = len(np.where(np.logical_and(s <= spikes, spikes <= e))[0])
                    freq = n_spikes / (e - s)

                    print(f"neuron_id: {neuron_id}, requested_freq: {requested_freq} Hz, actual freq: {freq}")
                    self.assertTrue(requested_freq * 0.7 < freq < requested_freq * 1.3,
                                    f"neuron_id: {neuron_id}, requested_freq: {requested_freq} Hz, actual freq: {freq}")

#        import pdb
#        pdb.set_trace()

        # TODO: Read in the frequency of each neuron. Compare the simulated frequency to the requested frequency

        self.assertEqual(True, True)  # add assertion here

    # def tearDown(self):
    #     self.pr.clear_neuron()


if __name__ == '__main__':
    unittest.main()
