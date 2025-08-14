import unittest
import numpy as np
from snudda.input import SnuddaInput


class MyTestCase(unittest.TestCase):

    def setUp(self) -> None:

        self.input = SnuddaInput()
        self.rng = np.random.default_rng(113)

    def test_lognormal_basic_functionality(self):
        """Test basic spike generation with single time range"""
        mean_freq = np.array([10.0])  # 10 Hz
        std_freq = np.array([2.0])  # 2 Hz std
        time_ranges = (np.array([0.0]), np.array([1.0]))  # 1 second

        spikes = self.input.generate_log_normal_spikes(mean_freq, std_freq, time_ranges, self.rng)

        # Check output properties
        self.assertIsInstance(spikes, np.ndarray)
        self.assertTrue(len(spikes) > 0, "Should generate some spikes")
        self.assertTrue(np.all(spikes >= 0.0), "All spikes should be >= start time")
        self.assertTrue(np.all(spikes < 1.0), "All spikes should be < end time")
        self.assertTrue(np.all(np.diff(spikes) > 0), "Spikes should be sorted")

    def test_lognormal_multiple_time_ranges(self):
        """Test with multiple time ranges"""
        mean_freq = np.array([5.0, 15.0])
        std_freq = np.array([1.0, 3.0])
        time_ranges = (np.array([0.0, 2.0]), np.array([1.0, 3.0]))

        spikes = self.input.generate_log_normal_spikes(mean_freq, std_freq, time_ranges, self.rng)

        # Check that spikes exist in both ranges
        spikes_range1 = spikes[(spikes >= 0.0) & (spikes < 1.0)]
        spikes_range2 = spikes[(spikes >= 2.0) & (spikes < 3.0)]

        no_spikes = spikes[(spikes > 1.0) & (spikes < 2.0)]

        self.assertTrue(len(spikes_range1) > 0, "Should have spikes in first range")
        self.assertTrue(len(spikes_range2) > 0, "Should have spikes in second range")
        self.assertTrue(np.all(np.diff(spikes) > 0), "All spikes should be sorted")
        self.assertTrue(len(no_spikes) == 0, "Middle part should be empty")

    def test_lognormal_empty_time_range(self):
        """Test with zero duration time range"""
        mean_freq = np.array([10.0])
        std_freq = np.array([2.0])
        time_ranges = (np.array([1.0]), np.array([1.0]))  # Zero duration

        spikes = self.input.generate_log_normal_spikes(mean_freq, std_freq, time_ranges, self.rng)

        self.assertEqual(len(spikes), 0, "Should generate no spikes for zero duration")

    def test_negative_duration(self):
        """Test with negative duration (end < start)"""
        mean_freq = np.array([10.0])
        std_freq = np.array([2.0])
        time_ranges = (np.array([2.0]), np.array([1.0]))  # Negative duration

        spikes = self.input.generate_log_normal_spikes(mean_freq, std_freq, time_ranges, self.rng)

        self.assertEqual(len(spikes), 0, "Should generate no spikes for negative duration")

    def test_high_frequency(self):
        """Test with high frequency to ensure many spikes are generated"""
        mean_freq = np.array([100.0])  # 100 Hz
        std_freq = np.array([20.0])
        time_ranges = (np.array([0.0]), np.array([1.0]))

        spikes = self.input.generate_log_normal_spikes(
            mean_freq, std_freq, time_ranges, self.rng
        )

        # Should have roughly 100 spikes (give or take due to randomness)
        self.assertGreater(len(spikes), 50, "Should generate many spikes at high frequency")
        self.assertLess(len(spikes), 200, "Shouldn't generate too many spikes")

    def test_low_frequency(self):
        """Test with very low frequency"""
        mean_freq = np.array([0.5])  # 0.5 Hz (one spike every 2 seconds on average)
        std_freq = np.array([0.1])
        time_ranges = (np.array([0.0]), np.array([1.0]))  # 1 second window

        spikes = self.input.generate_log_normal_spikes(
            mean_freq, std_freq, time_ranges, self.rng
        )

        # Might have 0 or 1 spike in 1 second window
        self.assertLessEqual(len(spikes), 3, "Should generate few spikes at low frequency")

    def test_reproducibility(self):
        """Test that same seed produces same results"""
        mean_freq = np.array([10.0])
        std_freq = np.array([2.0])
        time_ranges = (np.array([0.0]), np.array([1.0]))

        # Generate spikes twice with same seed
        rng1 = np.random.RandomState(123)
        rng2 = np.random.RandomState(123)

        spikes1 = self.input.generate_log_normal_spikes(mean_freq, std_freq, time_ranges, rng1)
        spikes2 = self.input.generate_log_normal_spikes(mean_freq, std_freq, time_ranges, rng2)

        np.testing.assert_array_equal(spikes1, spikes2, "Same seed should produce same spikes")

    def test_frequency_approximation(self):
        """Test that generated spike rate approximately matches input frequency"""
        mean_freq = np.array([20.0])  # 20 Hz
        std_freq = np.array([4.0])
        time_ranges = (np.array([0.0]), np.array([10.0]))  # Long duration for better statistics

        spikes = self.input.generate_log_normal_spikes(
            mean_freq, std_freq, time_ranges, self.rng
        )

        # Calculate actual frequency
        duration = 10.0
        actual_freq = len(spikes) / duration

        # Should be within reasonable range (allow for randomness)
        self.assertGreater(actual_freq, 15.0, "Frequency should be roughly correct (lower bound)")
        self.assertLess(actual_freq, 30.0, "Frequency should be roughly correct (upper bound)")

    def test_spike_time_bounds(self):
        """Test that all spikes fall within specified time bounds"""
        mean_freq = np.array([25.0, 10.0])
        std_freq = np.array([5.0, 2.0])
        time_ranges = (np.array([1.5, 5.0]), np.array([3.5, 8.0]))

        spikes = self.input.generate_log_normal_spikes(
            mean_freq, std_freq, time_ranges, self.rng
        )

        # All spikes should be within the union of time ranges
        valid_range1 = (spikes >= 1.5) & (spikes < 3.5)
        valid_range2 = (spikes >= 5.0) & (spikes < 8.0)
        all_valid = valid_range1 | valid_range2

        self.assertTrue(np.all(all_valid), "All spikes should be within specified time ranges")


if __name__ == '__main__':
    unittest.main()
