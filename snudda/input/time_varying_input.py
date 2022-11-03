import numpy as np


class TimeVaryingInput:

    def __init__(self):

        pass

    @staticmethod
    def get_stretch_time(frequency_function, end_time, start_time=0, dt=0.001, check_positive=True, start_at_zero=True):

        """ We want to stretch the time, so that a Poisson process with frequency 1Hz will result in a time varying
            Poisson process with the instantaneous frequency_function from 0 to end_time, with time resolution dt.

            Any frequency values < 0 are treated as 0.

            Args:
                frequency_function : a numpy compatible function that takes an array and returns an array
                end_time : end time of the valid time range
                start_time : start_time that function is defined for
                dt : time resolution of the sampling of the frequency function
                start_at_zero (bool) : If start_at_zero is True, then the frequency function is f(t=0)
                                       at the start of each interval, if False the simulation time is used.

        """

        time = np.arange(start_time, end_time, dt)

        if start_at_zero:
            func_time = np.arange(0, end_time - start_time, dt)
        else:
            func_time = time

        frequency = frequency_function(func_time)

        # If frequency was a scalar, extend it to be that frequency in entire time range
        if frequency.size == 1:
            frequency = np.full(time.shape, frequency)

        if check_positive:
            frequency[frequency <= 0] = 0

        try:
            stretched_time = np.cumsum(frequency*dt) - frequency[0]*dt  # We want stretched time to start at 0
            func = lambda t: np.interp(t, stretched_time, time)
            stretch_end_time = stretched_time[-1]
        except:
            import traceback
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()

        return func, stretch_end_time

    @staticmethod
    def _poisson_helper(end_time, rng):

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

    @staticmethod
    def generate_spikes(frequency_function, end_time, rng, start_time=0, n_spike_trains=1, start_at_zero=True):

        """
            Generates spikes with frequency f(t) where f is specified by frequency_function.

            Args:
                frequency_function : Either a python function, or a string of a function of t, e.g. "sin(2*pi*5*t)"
                end_time (float) : End time of spike train
                start_time (float) : Start time of spike train (default 0)
                n_spike_trains (int) : Number of spike trains, default 1
                start_at_zero (bool) : Is t=0 at the start of the time interval the function is active within?
                                       Default True. If set to False the simulation t is sent directly to the
                                       frequency_function.
        """

        stretch_func, stretched_end_time = TimeVaryingInput.get_stretch_time(frequency_function=frequency_function,
                                                                             start_time=start_time, end_time=end_time,
                                                                             start_at_zero=start_at_zero)
        spike_trains = []

        for idx in range(0, n_spike_trains):
            t_spike = TimeVaryingInput._poisson_helper(end_time=stretched_end_time, rng=rng)
            spike_trains.append(stretch_func(t_spike))  # Stretch the spike times

        return spike_trains

    @staticmethod
    def test_spike_frequency():

        func = lambda t: 5*np.cos(10*2*np.pi*t) + 6 + 10*t
        rng = np.random.default_rng()

        spikes = TimeVaryingInput.generate_spikes(frequency_function=func, start_time=1, end_time=3,
                                                  n_spike_trains=10000, rng=rng)

        all_spikes = np.concatenate(spikes)

        import matplotlib.pyplot as plt
        plt.figure()

        # We generate 3 seconds of data, 1000 trains, 3000 seconds of total data
        # we place that in 1000 bins, so if we weight result by 1/3 we should show approximate instantaneous frequency
        plt.hist(all_spikes, bins=1000, weights=np.full(all_spikes.shape, 1/20))
        plt.xlabel("Time (s)")
        plt.ylabel("Frequency (Hz)")
        plt.title("Oscillating between 1 and 11 Hz, periodicity 0.1s")
        plt.ion()
        plt.show()
        plt.pause(0.1)


if __name__ == "__main__":

    TimeVaryingInput.test_spike_frequency()
    import pdb
    pdb.set_trace()


