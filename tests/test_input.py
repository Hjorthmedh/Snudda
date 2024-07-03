import unittest
import os
import h5py
import json
import numpy as np

from copy import deepcopy

from snudda.detect.detect import SnuddaDetect
from snudda.input.input import SnuddaInput
from snudda.detect.prune import SnuddaPrune


class InputTestCase(unittest.TestCase):

    def setUp(self):

        print("RUNNING SETUP")
        os.chdir(os.path.dirname(__file__))

        self.network_path = os.path.join("networks", "network_testing_input")
        self.config_file = os.path.join(self.network_path, "network-config.json")
        self.position_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        self.save_file = os.path.join(self.network_path, "voxels", "network-putative-synapses.hdf5")
        self.network_file = os.path.join(self.network_path, "network-synapses.hdf5")

        # Setup network so we can test input generation
        from snudda.init.init import SnuddaInit
        cell_spec = os.path.join(os.path.dirname(__file__), "validation")
        cnc = SnuddaInit(struct_def={}, config_file=self.config_file, random_seed=1234)
        cnc.define_striatum(num_dSPN=5, num_iSPN=0, num_FS=5, num_LTS=0, num_ChIN=0,
                            volume_type="cube", neurons_dir=cell_spec)

        cnc.add_population_unit_random("Striatum", ["dSPN", "FS"], 0.4, unit_id=1)
        cnc.add_population_unit_random("Striatum", ["dSPN", "FS"], 0.4, unit_id=2)

        cnc.write_json(self.config_file)

        from snudda import Snudda
        snd = Snudda(network_path=self.network_path)
        snd.create_network()

    def test_generate(self):

        print("Checking generate_spikes.")

        rng = np.random.default_rng(777)
        si2 = SnuddaInput(verbose=True)
        freq = [100, 200, 300]
        start_times = [2, 5, 7]
        end_times = [3, 6, 9]
        spike_times = si2.generate_poisson_spikes(freq=freq, time_range=(start_times, end_times), rng=rng)

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

        p_keep2 = np.array([0.2, 0.5, 1])
        t_range = (start_times, end_times)
        culled_spike_times2 = si2.cull_spikes(spikes=spike_times, p_keep=p_keep2, rng=rng, time_range=t_range)

        for st, et, pk, f in zip(start_times, end_times, p_keep2, freq):
            n_spikes = len(np.where(np.logical_and(st <= culled_spike_times2,
                                                   culled_spike_times2 <= et))[0])
            n_expected = (et - st) * f * pk

            self.assertTrue(0.8*n_expected <= n_spikes <= 1.2*n_expected, f"Found {n_spikes}, expected {n_expected} (+/- 20%)")

        spike_times2 = si2.generate_poisson_spikes(freq=10, time_range=[0, 10], rng=rng)
        mixed_spike_times = si2.mix_spikes([spike_times, spike_times2])

        self.assertTrue((np.diff(mixed_spike_times) >= 0).all())
        self.assertTrue(len(spike_times) + len(spike_times2) == len(mixed_spike_times))

        jitter_dt = 10e-3
        jittered_spikes = si2.jitter_spikes(spike_trains=[spike_times], dt=jitter_dt, rng=rng)

        self.assertTrue((np.abs(spike_times - jittered_spikes[0]) < 4*jitter_dt).all())

    def test_input_1(self):

        # This tests Poisson inputs

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

        # Loop through all inputs, and verify them

        for neuron_id_str in input_data["input"].keys():
            neuron_id = int(neuron_id_str)
            neuron_name = si.network_data["neurons"][neuron_id]["name"]
            neuron_type = neuron_name.split("_")[0]

            # Check frequency is as advertised...
            for input_type in input_data["input"][neuron_id_str]:
                input_info = input_data["input"][neuron_id_str][input_type]

                start_time = input_info["spikes"].attrs["start"].copy()
                end_time = input_info["spikes"].attrs["end"].copy()
                freq = input_info["spikes"].attrs["freq"].copy()
                spikes = input_info["spikes"][()]
                n_traces = spikes.shape[0]

                if "num_inputs" in config_data[neuron_type][input_type]:
                    if "cluster_size" in config_data[neuron_type][input_type]:
                        cluster_size = config_data[neuron_type][input_type]["cluster_size"]
                    else:
                        cluster_size = 1

                    if isinstance(config_data[neuron_type][input_type]['num_inputs'], dict):
                        config_n_inputs = config_data[neuron_type][input_type]['num_inputs'][neuron_name]
                    else:
                        config_n_inputs = config_data[neuron_type][input_type]['num_inputs']
                    print(f"Checking number of inputs is {config_n_inputs} (cluster size used: {cluster_size})")
                    self.assertEqual(config_n_inputs, n_traces)

                    # TODO: We can no longer assume that section_id is the same for all inputs in a cluster
                    #       the new code also works at branch points, so cluster can be spread over different sections.
                    # if cluster_size > 1:
                    #     # Verify that all the clusters have the right size
                    #     for ctr in range(0, cluster_size-1):
                    #         self.assertTrue(np.all(np.diff(input_info["section_id"])[ctr::cluster_size] == 0))

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

                    if np.max(input_info["spikes"].attrs["correlation"]) == 0:
                        self.assertTrue(f_gen > f - 5*np.sqrt(f)/np.sqrt(n_traces))
                        self.assertTrue(f_gen < f + 5*np.sqrt(f)/np.sqrt(n_traces))
                    else:
                        # For high correlations and short durations we have huge fluctuations, so skip those
                        pass

                if "population_unit_correlation" in config_data[neuron_type][input_type]:
                    correlation = config_data[neuron_type][input_type]["population_unit_correlation"]

                    if "jitter" in config_data[neuron_type][input_type]:
                        jitter = config_data[neuron_type][input_type]["jitter"]
                    else:
                        jitter = 0

                    p_keep = np.sqrt(correlation)
                    if np.size(p_keep) == 1:
                        p_keep = np.full(np.size(start_time), p_keep)

                    # Is correlation what we expect?
                    bin_size = 6*jitter + 1e-3
                    n_bins = int(np.ceil(input_time / bin_size)) + 1
                    binned_data = np.zeros((n_bins,))

                    for t_idx in (spikes.flatten() / bin_size).astype(int):
                        if t_idx >= 0:
                            binned_data[t_idx] += 1

                    readout = np.zeros((spikes.size, ))
                    ctr = 0
                    for t_idx in (spikes.flatten() / bin_size).astype(int):
                        if t_idx > 0:
                            readout[ctr] = binned_data[t_idx]
                            ctr += 1

                    readout = readout[:ctr]

                    if np.size(freq) == 1:
                        mean_freq = freq[0]
                    else:
                        # Note this is the mean freq during period of spiking (since we dont sample silent periods)
                        mean_freq = np.sum(np.multiply(end_time - start_time, freq)) / np.sum(end_time - start_time)

                    # If we look at a spike in a spike train, then with P=p_keep it is a mother spike,
                    # and then there should be (N-1) * p_keep + 1 spikes in that bin.
                    # With P=(1-p_keep) it is just a normal spike, and then there should be 1 + f*dt*(N-1) spikes
                    # in the bin

                    # REMOVE THIS
                    #if np.size(freq) == 1:
                    #    expected_mean = (p_keep * ((n_traces - 1) * p_keep + 1 + freq[0] * bin_size * n_traces)
                    #                     + (1 - p_keep) * (1 + freq[0] * bin_size * (n_traces - 1)))
                    # END REMOVE

                    # When calculating expected mean number of simultaneous spikes for a bin with a spike
                    # we need to take into account that high freq periods are more likely, and they also have
                    # higher freq during that period
                    picked_ctr = 0
                    spike_cnt = 0

                    for st, et, f, p_k in zip(start_time, end_time, freq, p_keep):
                        picked_ctr += f*(et-st)  # Number of readouts in this time interval
                        spike_cnt += f*(et-st) * (p_k * ((n_traces - 1) * p_k + 1 + f * bin_size * n_traces * (1 - p_k))
                                                  + (1 - p_k) * (1 + f * bin_size * (n_traces - 1)))

                    expected_mean = spike_cnt / picked_ctr

                    print(f"Simultaneous spikes: {np.mean(readout):.2f} (expected {expected_mean:.2f}) "
                          f"- correlation {correlation}")
                    if jitter <= 0.001:
                        # Only do check for non-jittered input
                        self.assertTrue(expected_mean * 0.75 < np.mean(readout) < expected_mean * 1.25)

    def test_input_2(self):

        # This tests function based frequency input

        input_time = 0.5
        input_config = os.path.join(self.network_path, "input-test-2.json")
        spike_file = os.path.join(self.network_path, "input-spikes-2.hdf5")

        si = SnuddaInput(input_config_file=input_config,
                         hdf5_network_file=self.network_file,
                         spike_data_filename=spike_file,
                         time=input_time, verbose=True)
        si.generate()

        input_data = h5py.File(spike_file, 'r')
        config_data = json.loads(input_data["config"][()])

        # OBS, population unit 0 does not get any of the extra mother spikes specified
        # So we need to check FS neuron that belongs to population unit 1 or 2.
        some_spikes = input_data["input/4/Cortical/spikes"][()].flatten()
        some_spikes = some_spikes[some_spikes >= 0]
        n_trains = input_data["input/4/Cortical/spikes"][()].shape[0]

        for extra_spike in [0.2, 0.3, 0.45]:

            self.assertTrue(np.sum(np.abs(some_spikes - extra_spike) < 1e-4)
                            >= n_trains)
            self.assertTrue(np.sum(np.abs(some_spikes - extra_spike + 0.05) < 1e-3) < 50)

        some_spikes2 = input_data["input/4/Thalamic/spikes"][()].flatten()
        some_spikes2 = some_spikes2[some_spikes2 >= 0]

        for spike in [0.1, 0.2, 0.3]:
            self.assertTrue(np.sum(np.abs(some_spikes2 - spike) < 1e-4) == 2000)

        self.assertTrue(np.size(some_spikes2) == 6000)

        # Check input generated, this focuses on the frequency function generation
        # and also checks input correlation

        # TODO: New cell numbering, so need to pick other cell numbers
        some_spikes_c0 = input_data["input/0/CorticalSignal/spikes"][()]
        some_spikes_c1 = input_data["input/1/CorticalSignal/spikes"][()]

        pop0 = input_data["input/0/CorticalSignal/population_unit_spikes"][()]
        pop1 = input_data["input/1/CorticalSignal/population_unit_spikes"][()]

        # TODO: Add checks

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
                    self.assertTrue((t_check-1)*80 <= freq <= (t_check-1)*120,
                                    f"Found frequency {freq} Hz at {t_check}s, expected {t_check*100} Hz")


    def test_arbitrary_function_range(self):

        func_lambda = lambda t: t*100
        func_str = "t*100"  # OBS, t=0 at the start of each new stimulus

        t_range = [[1, 4], [2, 7]]

        si_empty = SnuddaInput()

        for func in [func_str, func_lambda]:

            # We run this twice, for string functions and for lambda functions

            rng = np.random.default_rng(112)
            spikes = si_empty.generate_spikes_function(frequency_function=func, time_range=t_range, dt=1e-4, rng=rng)

            with self.subTest("Freq test"):
                self.assertTrue(40 <= self.find_freq_in_range(spikes, [1, 2]) <= 60,
                                f"Expected frequency 50Hz, found {self.find_freq_in_range(spikes, [1, 2])} Hz")

                self.assertTrue(self.find_freq_in_range(spikes, [2, 4]) == 0,
                                f"Expected frequency 0Hz, found {self.find_freq_in_range(spikes, [2, 4])} Hz")

                self.assertTrue(135 <= self.find_freq_in_range(spikes, [4, 7]) <= 165,
                                f"Expected frequency 150Hz, found {self.find_freq_in_range(spikes, [4, 5])} Hz")

    def test_fraction_mixing(self):

        si_empty = SnuddaInput()
        rng = np.random.default_rng(112)

        spikes_a = np.arange(0, 10, 0.1)
        spikes_b = np.arange(0.01, 10, 0.1)
        fraction_a = [0.9, 0.1]
        fraction_b = [0.1, 0.9]
        time_range = [[0, 5], [5, 8]]

        mixed_spikes = SnuddaInput.mix_fraction_of_spikes(spikes_a=spikes_a, spikes_b=spikes_b,
                                                          fraction_a=fraction_a, fraction_b=fraction_b,
                                                          rng=rng, time_range=time_range)

        # No spikes after t=8s, since outside time range
        self.assertEqual(np.sum(mixed_spikes > 8), 0)

        n_a_1 = np.logical_and((mixed_spikes + 1e-7) % 0.1 < 1e-5, mixed_spikes <= 5)
        n_b_1 = np.logical_and((mixed_spikes + 1e-7) % 0.1 > 1e-5, mixed_spikes <= 5)
        n_a_2 = np.logical_and((mixed_spikes + 1e-7) % 0.1 < 1e-5, mixed_spikes >= 5)
        n_b_2 = np.logical_and((mixed_spikes + 1e-7) % 0.1 > 1e-5, mixed_spikes >= 5)

        self.assertTrue(40 < np.sum(n_a_1) < 50)
        self.assertTrue(1 < np.sum(n_b_1) < 10)
        self.assertTrue(1 < np.sum(n_a_2) < 6)
        self.assertTrue(20 < np.sum(n_b_2) < 30)

    def find_spikes_in_range(self, spikes, time_range):
        t_idx = np.where(np.logical_and(time_range[0] <= spikes, spikes <= time_range[1]))[0]
        return spikes[t_idx]

    def find_freq_in_range(self, spikes, time_range):
        return len(self.find_spikes_in_range(spikes, time_range)) / (time_range[1] - time_range[0])


if __name__ == '__main__':
    unittest.main()
