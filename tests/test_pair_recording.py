import os
import unittest
from snudda.simulate.pair_recording import PairRecording


class PairRecordingTestCase(unittest.TestCase):

    def setUp(self):

        # Temporarily disable the creation of network while testing...
        return

        network_path = os.path.join("networks", "pair_recording_test")
        rc = None

        if os.path.isdir("x86_64"):
            import shutil
            shutil.rmtree("x86_64")
        os.system(f"nrnivmodl {os.path.join('validation', 'mechanisms')}")

        from snudda import SnuddaPlace
        from snudda import SnuddaDetect
        from snudda import SnuddaPrune

        sp = SnuddaPlace(network_path=network_path, rc=rc)
        sp.place()

        sd = SnuddaDetect(network_path=network_path, rc=rc)
        sd.detect()

        sp = SnuddaPrune(network_path=network_path, rc=rc)
        sp.prune()

        experiment_config_file = os.path.join(network_path, "experiment.json")
        network_file = os.path.join(network_path, "network-synapses.hdf5")
        self.pr = PairRecording(network_path=network_path, network_file=network_file,
                                experiment_config_file=experiment_config_file)
        self.pr.run()

    def test_frequency(self):

        # TODO: Read in the frequency of each neuron. Compare the simulated frequency to the requested frequency

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
