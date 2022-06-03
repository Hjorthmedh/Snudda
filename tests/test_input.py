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

    def test_generate(self):

        print("Checking generate_spikes.")

        rng = np.random.default_rng(777)
        si2 = SnuddaInput(verbose=True)
        freq = [100, 200, 300]
        start_times = [2, 5, 7]
        end_times = [3, 6, 9]
        spike_times = si2.generate_spikes(freq=freq, time_range=(start_times, end_times), rng=rng)

        for st, et, f in zip(start_times, end_times, freq):
            t_idx = np.where(np.logical_and(st <= spike_times, spike_times <= et))[0]
            f_gen = len(t_idx) / (et - st)
            print(f"test_generate: Expected frequency {f}, generated frequency {f_gen}")
            self.assertTrue(f - np.sqrt(f)*2 < f_gen < f + np.sqrt(f)*2,
                            f"test_generate: Expected frequency {f}, generated frequency {f_gen}")

        # Also check the in between ranges, that they are empty
        print("Checking in between ranges empty.")
        for st, et in zip(end_times[:-2], start_times[1:]):
            t_idx = np.where(np.logical_and(st <= spike_times, spike_times <= et))[0]
            self.assertEqual(len(t_idx), 0, f"Time range {st} to {et} should be empty but contains spikes.")

        p_keep = 0.3
        culled_spike_times = si2.cull_spikes(spikes=spike_times, p_keep=p_keep, rng=rng)
        s = np.sqrt(len(spike_times) * p_keep * (1 - p_keep))

        self.assertTrue(p_keep*len(spike_times) - 2*s < len(culled_spike_times) < p_keep * len(spike_times) + 2*s,
                        f"Problem with culling of spikes. n_before={len(spike_times)}, "
                        f"n_after={len(culled_spike_times)} (expected {len(spike_times)*p_keep} +/- {2*s}), "
                        f"p_keep={p_keep}")

        spike_times2 = si2.generate_spikes(freq=10, time_range=[0, 10], rng=rng)
        mixed_spike_times = si2.mix_spikes([spike_times, spike_times2])

        self.assertTrue((np.diff(mixed_spike_times) >= 0).all())
        self.assertTrue(len(spike_times) + len(spike_times2) == len(mixed_spike_times))

        jitter_dt = 10e-3
        jittered_spikes = si2.jitter_spikes(spike_trains=[spike_times], dt=jitter_dt, rng=rng)

        self.assertTrue((np.abs(spike_times - jittered_spikes[0]) < 4*jitter_dt).all())

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

                    try:
                        self.assertTrue(f_gen > f - 5*np.sqrt(f)/np.sqrt(n_traces))
                        self.assertTrue(f_gen < f + 5*np.sqrt(f)/np.sqrt(n_traces))
                    except:
                        import pdb
                        import traceback
                        print(traceback.format_exc())
                        pdb.set_trace()

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

                    print(f"Simultaneous spikes: {np.mean(readout):.2f} (expected {expected_mean:.2f}) "
                          f"- correlation {correlation}")
                    try:
                        self.assertTrue(expected_mean * 0.9 < np.mean(readout) < expected_mean * 1.1)

                    except:
                        import traceback

                        t_str = traceback.format_exc()
                        print(t_str)
                        import pdb

                        pdb.set_trace()

    def test_arbitrary_function(self):

        func_lambda = lambda t: t*100
        func_str = "t*100"

        si_empty = SnuddaInput()

        for func in [func_str, func_lambda]:

            # We run this twice, for string functions and for lambda functions

            rng = np.random.default_rng(112)
            spikes = si_empty.generate_spikes_function(frequency_function=func, time_range=[1, 10], dt=1e-4, rng=rng)
            isi = np.diff(spikes)

            self.assertTrue((isi > 0).all(), "Resulting spikes should be sorted")

            with self.subTest("Checking no spikes before t=1"):
                t_idx = np.where(spikes < 1)[0]
                self.assertTrue(len(t_idx) == 0, "There should be no spikes before t=1")

            # Check average frequency at around 3,5,7,9 seconds.

            for t_check in [3, 5, 7, 9]:
                t_range = [t_check-0.1, t_check+0.1]
                freq = self.find_freq_in_range(spikes, t_range)

                with self.subTest(f"Checking frequency at {t_check} (expecting around {t_check*100} Hz)"):
                    self.assertTrue(t_check*80 <= freq <= t_check*120,
                                    f"Found frequency {freq} Hz at {t_check}s, expected {t_check*100} Hz")

    def find_spikes_in_range(self, spikes, time_range):
        t_idx = np.where(np.logical_and(time_range[0] <= spikes, spikes <= time_range[1]))[0]
        return spikes[t_idx]

    def find_freq_in_range(self, spikes, time_range):
        return len(self.find_spikes_in_range(spikes, time_range)) / (time_range[1] - time_range[0])

    def test_arbitrary_function_range(self):

        func_lambda = lambda t: t*100
        func_str = "t*100"

        t_range = [[1, 4], [2, 5]]

        si_empty = SnuddaInput()

        for func in [func_str, func_lambda]:

            # We run this twice, for string functions and for lambda functions

            rng = np.random.default_rng(112)
            spikes = si_empty.generate_spikes_function(frequency_function=func, time_range=t_range, dt=1e-4, rng=rng)

            with self.subTest("Freq test"):
                self.assertTrue(100 <= self.find_freq_in_range(spikes, [1, 2]) <= 200,
                                f"Expected frequency 150Hz, found {self.find_freq_in_range(spikes, [1, 2])} Hz")

                self.assertTrue(self.find_freq_in_range(spikes, [2, 4]) == 0,
                                f"Expected frequency 0Hz, found {self.find_freq_in_range(spikes, [2, 4])} Hz")

                self.assertTrue(400 <= self.find_freq_in_range(spikes, [4, 5]) <= 500,
                                f"Expected frequency 500Hz, found {self.find_freq_in_range(spikes, [4, 5])} Hz")


if __name__ == '__main__':
    unittest.main()
