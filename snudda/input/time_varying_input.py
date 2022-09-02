import numpy as np


class TimeVaryingInput:

    def __init__(self):

        pass

    def get_stretch_time(self, frequency_function, end_time, dt=0.01, check_positive=True):

        """ We want to stretch the time, so that a Poisson process with frequency 1Hz will result in a time varying
            Poisson process with the instantaneous frequncy_function from 0 to end_time, with time resolution dt.

            Any frequency values < 0 are treated as 0

            Args:
                frequency_function : a numpy compatible function that takes an array and returns an array
                end_time : end time of the valid time range
                dt : time resolution of the sampling of the frequency function

        """

        time = np.arange(0, end_time, dt)
        frequency = frequency_function(time)

        if check_positive:
            frequency[frequency <= 0] = 0

        stretched_time = np.cumsum(frequency*dt) - frequency[0]*dt  # We want stretched time to start at 0
        func = lambda t: np.interp(t, stretched_time, time)
        stretch_end_time = stretched_time[-1]

        return func, stretch_end_time

    def _poisson_helper(self, end_time, rng):

        t_diff = -np.log(1.0 - rng.random(int(np.ceil(max(1, end_time)))))  # Generate 1 Hz spiking
        t_spikes = [np.cumsum(t_diff)]

        # Is last spike after end of duration
        while t_spikes[-1][-1] <= end_time:
            t_diff = -np.log(1.0 - rng.random(int(np.ceil(end_time * 0.1))))
            t_spikes.append(t_spikes[-1][-1] + np.cumsum(t_diff))

        # Prune away any spikes after end
        if len(t_spikes[-1]) > 0:
            t_spikes[-1] = t_spikes[-1][t_spikes[-1] <= end_time]

        # Return spike times
        return np.concatenate(t_spikes)

    def generate_spikes(self, frequency_function, end_time, rng, n_spike_trains=1):

        stretch_func, stretched_end_time = self.get_stretch_time(frequency_function=frequency_function, end_time=end_time)
        spike_trains = []

        for idx in range(0, n_spike_trains):
            t_spike = self._poisson_helper(end_time=stretched_end_time, rng=rng)
            spike_trains.append(stretch_func(t_spike))  # Stretch the spike times

        return spike_trains

    def test_spike_frequency(self):

        func = lambda t: 5*np.cos(10*2*np.pi*t) + 6
        rng = np.random.default_rng()

        spikes = self.generate_spikes(frequency_function=func, end_time=3, n_spike_trains=10000, rng=rng)

        all_spikes = np.concatenate(spikes)

        import matplotlib.pyplot as plt
        plt.figure()

        # We generate 3 seconds of data, 1000 trains, 3000 seconds of total data
        # we place that in 1000 bins, so if we weight result by 1/3 we should show approximate instantaneous frequency
        plt.hist(all_spikes, bins=1000, weights=np.full(all_spikes.shape, 1/30))
        plt.xlabel("Time (s)")
        plt.ylabel("Frequency (Hz)")
        plt.title("Oscillating between 1 and 11 Hz, periodicity 0.1s")
        plt.ion()
        plt.show()
        plt.pause(0.1)

        import pdb
        pdb.set_trace()


if __name__ == "__main__":

    tvi = TimeVaryingInput()

    tvi.test_spike_frequency()


