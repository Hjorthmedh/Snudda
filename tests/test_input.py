import unittest
import os
import h5py
import json

from snudda.detect import SnuddaDetect
from snudda.input import SnuddaInput
from snudda.prune import SnuddaPrune


class MyTestCase(unittest.TestCase):

    def setUp(self):

        os.chdir(os.path.dirname(__file__))

        self.network_path = os.path.join("tests", "network_testing_input")
        self.config_file = os.path.join(self.network_path, "network-config.json")
        self.position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        self.save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        # Setup network so we can test input generation
        from snudda.init import SnuddaInit
        cell_spec = os.path.join(os.path.dirname(__file__), "validation")
        cnc = SnuddaInit(struct_def={}, config_name=self.config_file, num_population_units=1, random_seed=1234)
        cnc.define_striatum(num_dSPN=10, num_iSPN=0, num_FS=10, num_LTS=0, num_ChIN=0,
                            volume_type="cube", cell_spec_dir=cell_spec)
        cnc.write_json(self.config_file)

        # Place neurons
        from snudda.place import SnuddaPlace
        npn = SnuddaPlace(config_file=self.config_file,
                          log_file=None,
                          verbose=True,
                          d_view=None,          # TODO: If d_view is None code run sin serial, add test parallel
                          h5libver="latest")
        npn.read_config()
        npn.write_data(self.position_file)

        # Detect
        self.sd = SnuddaDetect(config_file=self.config_file, position_file=self.position_file,
                               save_file=self.save_file, rc=None,
                               hyper_voxel_size=120)

        self.sd.detect(restart_detection_flag=True)

        # Prune
        self.work_log = os.path.join(self.network_path, "log", "network-detect-worklog.hdf5")
        self.network_file = os.path.join(self.network_path, "network-pruned-synapses.hdf5")

        sp = SnuddaPrune(work_history_file=self.work_log, config_file=None)  # Use default config file
        sp.prune(pre_merge_only=False)

    def test_input_1(self):

        input_time = 10
        input_config = os.path.join(self.network_path, "input-test-1.json")
        spike_file = os.path.join(self.network_path, "input-spikes.hdf5")

        si = SnuddaInput(input_config_file=input_config,
                         hdf5_network_file=self.network_file,
                         spike_data_filename=spike_file,
                         time=input_time)
        si.generate()

        input_data = h5py.File(spike_file, 'r')
        config_data = json.loads(input_data["config"][()])

        #TODO: Add checks
        #import pdb
        #pdb.set_trace()

        pass

    # TODO: Test the disabling of input by adding a ! mark in config file before input name

if __name__ == '__main__':
    unittest.main()
