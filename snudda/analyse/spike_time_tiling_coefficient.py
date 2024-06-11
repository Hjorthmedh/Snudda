# This code is taken from https://github.com/NeuralEnsemble/elephant/blob/master/elephant/spike_train_correlation.py
# and modified to use numba for faster run time.

import numpy as np
from numba import jit


@jit(nopython=True, fastmath=True, cache=True)
def spike_time_tiling_coefficient(spiketrain_a, spiketrain_b, end_time, start_time=0, dt=0.005):
    """
    Calculates the Spike Time Tiling Coefficient (STTC) as described in
    :cite:`correlation-Cutts2014_14288` following their implementation in C.
    The STTC is a pair wise measure of correlation between spike trains.
    It has been proposed as a replacement for the correlation index as it
    presents several advantages (e.g. it's not confounded by firing rate,
    appropriately distinguishes lack of correlation from anti-correlation,
    periods of silence don't add to the correlation and it's sensitive to
    firing patterns).
    The STTC is calculated as follows:
    .. math::
        STTC = 1/2((PA - TB)/(1 - PA*TB) + (PB - TA)/(1 - PB*TA))
    Where `PA` is the proportion of spikes from train 1 that lie within
    `[-dt, +dt]` of any spike of train 2 divided by the total number of spikes
    in train 1, `PB` is the same proportion for the spikes in train 2;
    `TA` is the proportion of total recording time within `[-dt, +dt]` of any
    spike in train 1, TB is the same proportion for train 2.
    For :math:`TA = PB = 1`and for :math:`TB = PA = 1`
    the resulting :math:`0/0` is replaced with :math:`1`,
    since every spike from the train with :math:`T = 1` is within
    `[-dt, +dt]` of a spike of the other train.
    This is a Python implementation compatible with the elephant library of
    the original code by C. Cutts written in C and available `here
    <https://github.com/CCutts/Detecting_pairwise_correlations_in_spike_trains/
    blob/master/spike_time_tiling_coefficient.c>`_:
    """

    N1 = len(spiketrain_a)
    N2 = len(spiketrain_b)

    if N1 == 0 or N2 == 0:
        index = np.nan
    else:
        TA = run_t(spiketrain_a, end_time=end_time, start_time=start_time, dt=dt)
        TB = run_t(spiketrain_b, end_time=end_time, start_time=start_time, dt=dt)
        PA = run_p(spiketrain_a, spiketrain_b, dt=dt)
        PA = PA / N1
        PB = run_p(spiketrain_b, spiketrain_a, dt=dt)
        PB = PB / N2
        # check if the P and T values are 1 to avoid division by zero
        # This only happens for TA = PB = 1 and/or TB = PA = 1,
        # which leads to 0/0 in the calculation of the index.
        # In those cases, every spike in the train with P = 1
        # is within dt of a spike in the other train,
        # so we set the respective (partial) index to 1.
        if PA * TB == 1:
            if PB * TA == 1:
                index = 1.
            else:
                index = 0.5 + 0.5 * (PB - TA) / (1 - PB * TA)
        elif PB * TA == 1:
            index = 0.5 + 0.5 * (PA - TB) / (1 - PA * TB)
        else:
            index = 0.5 * (PA - TB) / (1 - PA * TB) + 0.5 * (PB - TA) / (1 - PB * TA)

    return index


@jit(nopython=True, fastmath=True, cache=True)
def run_p(spiketrain_a, spiketrain_b, dt=0.005):
    """
    Check every spike in train 1 to see if there's a spike in train 2
    within dt
    """
    N2 = len(spiketrain_b)

    # Search spikes of spiketrain_i in spiketrain_j
    # ind will contain index of
    ind = np.searchsorted(spiketrain_b, spiketrain_a)

    # To prevent IndexErrors
    # If a spike of spiketrain_i is after the last spike of spiketrain_j,
    # the index is N2, however spiketrain_j[N2] raises an IndexError.
    # By shifting this index, the spike of spiketrain_i will be compared
    # to the last 2 spikes of spiketrain_j (negligible overhead).
    # Note: Not necessary for index 0 that will be shifted to -1,
    # because spiketrain_j[-1] is valid (additional negligible comparison)
    ind[ind == N2] = N2 - 1

    # Compare to nearest spike in spiketrain_j BEFORE spike in spiketrain_i
    close_left = np.abs(spiketrain_b[ind - 1] - spiketrain_a) <= dt
    # Compare to nearest spike in spiketrain_j AFTER (or simultaneous)
    # spike in spiketrain_j
    close_right = np.abs(spiketrain_b[ind] - spiketrain_a) <= dt

    # spiketrain_j spikes that are in [-dt, dt] range of spiketrain_i
    # spikes are counted only ONCE (as per original implementation)
    close = close_left + close_right

    # Count how many spikes in spiketrain_i have a "partner" in
    # spiketrain_j
    return np.count_nonzero(close)


@jit(nopython=True, fastmath=True, cache=True)
def run_t(spiketrain, end_time, start_time=0, dt=0.005):
    """
    Calculate the proportion of the total recording time 'tiled' by spikes.
    """
    N = len(spiketrain)
    time_A = 2 * N * dt  # maximum possible time

    if N == 1:  # for only a single spike in the train

        # Check difference between start of recording and single spike
        if spiketrain[0] - start_time < dt:
            time_A += - dt + spiketrain[0] - start_time

        # Check difference between single spike and end of recording
        elif spiketrain[0] + dt > end_time:
            time_A += - dt - spiketrain[0] + end_time

    else:  # if more than a single spike in the train

        # Calculate difference between consecutive spikes
        diff = np.diff(spiketrain)

        # Find spikes whose tiles overlap
        idx = np.where(diff < 2 * dt)[0]
        # Subtract overlapping "2*dt" tiles and add differences instead
        time_A += - 2 * dt * len(idx) + diff[idx].sum()

        # Check if spikes are within +/-dt of the start and/or end
        # if so, subtract overlap of first and/or last spike
        if spiketrain[0] - start_time < dt:
            time_A += spiketrain[0] - start_time - dt
        if (end_time - spiketrain[N - 1]) < dt:
            time_A += - spiketrain[-1] - dt + end_time

    # Calculate the proportion of total recorded time to "tiled" time
    T = time_A / end_time
    return T


if __name__ == "__main__":

    time_a = np.sort(np.random.uniform(low=0, high=10, size=10))
    time_b = np.sort(np.random.uniform(low=0, high=10, size=10))

    run_p(time_a, time_b, dt=0.005)
    run_t(time_a, end_time=10, dt=0.005)

    spike_time_tiling_coefficient(time_a, time_b, end_time=10)

