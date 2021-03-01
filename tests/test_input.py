import unittest
import os
import h5py
import json
import numpy as np

from snudda.detect import SnuddaDetect
from snudda.input import SnuddaInput
from snudda.prune import SnuddaPrune


class MyTestCase(unittest.TestCase):

    def setUp(self):

        os.chdir(os.path.dirname(__file__))

        self.network_path = os.path.join("networks", "network_testing_input")
        self.config_file = os.path.join(self.network_path, "network-config.json")
        self.position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        self.save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        # Setup network so we can test input generation
        from snudda.init import SnuddaInit
        cell_spec = os.path.join(os.path.dirname(__file__), "validation")
        cnc = SnuddaInit(struct_def={}, config_file=self.config_file, random_seed=1234)
        cnc.define_striatum(num_dSPN=10, num_iSPN=0, num_FS=10, num_LTS=0, num_ChIN=0,
                            volume_type="cube", neurons_dir=cell_spec)
        cnc.write_json(self.config_file)

        # Place neurons
        from snudda.place import SnuddaPlace
        npn = SnuddaPlace(config_file=self.config_file,
                          log_file=None,
                          verbose=True,
                          d_view=None,          # TODO: If d_view is None code run sin serial, add test parallel
                          h5libver="latest")
        npn.parse_config()
        npn.write_data(self.position_file)

        # Detect
        self.sd = SnuddaDetect(config_file=self.config_file, position_file=self.position_file,
                               save_file=self.save_file, rc=None,
                               hyper_voxel_size=120, verbose=True)

        self.sd.detect(restart_detection_flag=True)

        # Prune
        self.network_file = os.path.join(self.network_path, "network-pruned-synapses.hdf5")

        sp = SnuddaPrune(network_path=self.network_path, config_file=None)  # Use default config file
        sp.prune(pre_merge_only=False)

    def test_input_1(self):

        input_time = 10
        input_config = os.path.join(self.network_path, "input-test-1.json")
        spike_file = os.path.join(self.network_path, "input-spikes.hdf5")

        si = SnuddaInput(input_config_file=input_config,
                         hdf5_network_file=self.network_file,
                         spike_data_filename=spike_file,
                         time=input_time, verbose=True)
        si.generate()

        input_data = h5py.File(spike_file, 'r')
        config_data = json.loads(input_data["config"][()])

        # TODO: Add checks

        # Loop through all inputs, and verify them

        for neuron_id_str in input_data["input"].keys():
            neuron_id = int(neuron_id_str)
            neuron_name = si.network_info["neurons"][neuron_id]["name"]
            neuron_type = neuron_name.split("_")[0]

            # Check frequency is as advertised...
            for input_type in input_data["input"][neuron_id_str]:
                input_info = input_data["input"][neuron_id_str][input_type]

                start_time = input_info["start"][()].copy()
                end_time = input_info["end"][()].copy()
                freq = input_info["freq"][()].copy()
                spikes = input_info["spikes"][()]
                n_traces = spikes.shape[0]

                max_len = 1
                if type(start_time) is np.ndarray:
                    max_len = np.maximum(max_len, len(start_time))

                if type(end_time) is np.ndarray:
                    max_len = np.maximum(max_len, len(end_time))

                if type(freq) is np.ndarray:
                    max_len = np.maximum(max_len, len(freq))

                if type(start_time) != np.ndarray:
                    start_time = np.array([start_time]*max_len)

                if type(end_time) != np.ndarray:
                    end_time = np.array([end_time]*max_len)

                if type(freq) != np.ndarray:
                    freq = np.array([freq]*max_len)

                for st, et, f in zip(start_time, end_time, freq):
                    idx_x, idx_y = np.where(np.logical_and(st <= spikes, spikes <= et))

                    f_gen = len(idx_x)/(n_traces * (et-st))
                    print(f"ID {neuron_id_str} {neuron_name} {input_type} f={f}, f_gen={f_gen}")

                    self.assertTrue(f_gen > f - 4*np.sqrt(f)/np.sqrt(n_traces))
                    self.assertTrue(f_gen < f + 4*np.sqrt(f)/np.sqrt(n_traces))


if __name__ == '__main__':
    unittest.main()
