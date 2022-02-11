import unittest
import os
import h5py
import json
import numpy as np

from snudda.detect.detect import SnuddaDetect
from snudda.input.input import SnuddaInput
from snudda.detect.prune import SnuddaPrune


class InputTestCase(unittest.TestCase):

    def setUp(self):

        os.chdir(os.path.dirname(__file__))

        self.network_path = os.path.join("networks", "network_testing_input")
        self.config_file = os.path.join(self.network_path, "network-config.json")
        self.position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        self.save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")

        # Setup network so we can test input generation
        from snudda.init.init import SnuddaInit
        cell_spec = os.path.join(os.path.dirname(__file__), "validation")
        cnc = SnuddaInit(struct_def={}, config_file=self.config_file, random_seed=1234)
        cnc.define_striatum(num_dSPN=5, num_iSPN=0, num_FS=5, num_LTS=0, num_ChIN=0,
                            volume_type="cube", neurons_dir=cell_spec)
        cnc.write_json(self.config_file)

        # Place neurons
        from snudda.place.place import SnuddaPlace
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
        self.network_file = os.path.join(self.network_path, "network-synapses.hdf5")

        sp = SnuddaPrune(network_path=self.network_path, config_file=None)  # Use default config file
        sp.prune()

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
            neuron_name = si.network_data["neurons"][neuron_id]["name"]
            neuron_type = neuron_name.split("_")[0]

            # Check frequency is as advertised...
            for input_type in input_data["input"][neuron_id_str]:
                input_info = input_data["input"][neuron_id_str][input_type]

                start_time = input_info["start"][()].copy()
                end_time = input_info["end"][()].copy()
                freq = input_info["freq"][()].copy()
                spikes = input_info["spikes"][()]
                n_traces = spikes.shape[0]

                if "nInputs" in config_data[neuron_type][input_type]:
                    if "clusterSize" in config_data[neuron_type][input_type]:
                        cluster_size = config_data[neuron_type][input_type]["clusterSize"]
                    else:
                        cluster_size = 1

                    if isinstance(config_data[neuron_type][input_type]['nInputs'], dict):
                        config_n_inputs = config_data[neuron_type][input_type]['nInputs'][neuron_name]
                    else:
                        config_n_inputs = config_data[neuron_type][input_type]['nInputs']
                    print(f"Checking number of inputs is {config_n_inputs} * {cluster_size}")
                    self.assertEqual(config_n_inputs * cluster_size, n_traces)

                    if cluster_size > 1:
                        # Verify that all the clusters have the right size
                        for ctr in range(0, cluster_size-1):
                            self.assertTrue(np.all(np.diff(input_info["sectionID"])[ctr::cluster_size] == 0))

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
                    t_idx = np.where(np.logical_and(st <= spikes, spikes <= et))[0]

                    f_gen = len(t_idx) / (n_traces * (et - st))
                    print(f"ID {neuron_id_str} {neuron_name} {input_type} f={f}, f_gen={f_gen}")

                    self.assertTrue(f_gen > f - 4*np.sqrt(f)/np.sqrt(n_traces))
                    self.assertTrue(f_gen < f + 4*np.sqrt(f)/np.sqrt(n_traces))

                if "populationUnitCorrelation" in config_data[neuron_type][input_type]:
                    correlation = config_data[neuron_type][input_type]["populationUnitCorrelation"]

                    if "jitter" in config_data[neuron_type][input_type]:
                        jitter = config_data[neuron_type][input_type]["jitter"]
                    else:
                        jitter = 0

                    p_keep = 1 / (n_traces - np.sqrt(correlation) * (n_traces - 1))

                    # Is correlation what we expect?
                    bin_size = 2*jitter + 1e-3
                    n_bins = int(np.ceil(input_time / bin_size)) + 1
                    binned_data = np.zeros((n_bins,))

                    for t_idx in (spikes.flatten() / bin_size).astype(int):
                        if t_idx >= 0:
                            binned_data[t_idx] += 1

                    readout = np.zeros((spikes.size, ))
                    ctr = 0
                    for t_idx in (spikes.flatten() / bin_size).astype(int):
                        try:
                            if t_idx > 0:
                                readout[ctr] = binned_data[t_idx]
                                ctr += 1
                        except:
                            import traceback
                            t_str = traceback.format_exc()
                            print(t_str)
                            import pdb
                            pdb.set_trace()

                    readout = readout[:ctr]

                    if len(freq) == 1:
                        mean_freq = freq[0]
                    else:
                        # Note this is the mean freq during period of spiking (since we dont sample silent periods)
                        mean_freq = np.sum(np.multiply(end_time - start_time, freq)) / np.sum(end_time - start_time)

                    # If we look at a spike in a spike train, then with P=p_keep it is a mother spike,
                    # and then there should be (N-1) * p_keep + 1 spikes in that bin.
                    # With P=(1-p_keep) it is just a normal spike, and then there should be 1 + f*dt*(N-1) spikes
                    # in the bin

                    if len(freq) == 1:
                        expected_mean = (p_keep * ((n_traces - 1) * p_keep + 1 + freq[0] * bin_size * n_traces)
                                         + (1 - p_keep) * (1 + freq[0] * bin_size * (n_traces - 1)))

                    else:
                        # When calculating expected mean number of simultaneous spikes for a bin with a spike
                        # we need to take into account that high freq periods are more likely, and they also have
                        # higher freq during that period
                        picked_ctr = 0
                        spike_cnt = 0
                        for st, et, f in zip(start_time, end_time, freq):
                            picked_ctr += f*(et-st)  # Number of readouts in this time interval
                            spike_cnt += f*(et-st) * (p_keep * ((n_traces - 1) * p_keep + 1 + f * bin_size * n_traces)
                                                      + (1 - p_keep) * (1 + f * bin_size * (n_traces - 1)))

                        expected_mean = spike_cnt / picked_ctr

#                        expected_mean = (np.sum(np.multiply((p_keep * ((n_traces - 1) * p_keep + 1)
#                                                            + (1 - p_keep) * (1 + freq * bin_size * (n_traces - 1))),
#                                                            np.multiply(end_time - start_time, freq)))
#                                         / (np.sum(np.multiply(end_time - start_time, freq))))

                    print(f"Simultaneous spikes: {np.mean(readout):.2f} (expected {expected_mean:.2f}) "
                          f"- correlation {correlation}")
                    self.assertTrue(expected_mean * 0.9 < np.mean(readout) < expected_mean * 1.1)


if __name__ == '__main__':
    unittest.main()
