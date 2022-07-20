import math
import os
import numpy as np
from scipy.ndimage.filters import uniform_filter1d
import pathlib
import quantities as pq
import neo
import elephant as elp

def compare_normal_distribution(volt, mu, sigma, title_name):
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy.stats as stats
    import math

    ##computing the bin properties (same for both distributions)
    num_bin = 100
    bin_lims = np.linspace(-100, 0, num_bin + 1)
    bin_centers = 0.5 * (bin_lims[:-1] + bin_lims[1:])
    bin_widths = bin_lims[1:] - bin_lims[:-1]

    ##computing the histograms
    hist1, _ = np.histogram(volt, bins=bin_lims)

    ##normalizing
    hist1b = hist1 / np.max(hist1)

    fig, ax = plt.subplots(nrows=1, ncols=1)
    dist2 = np.linspace(mu - 3 * sigma, mu + 3 * sigma, 100)
    ax.plot(dist2, stats.norm.pdf(dist2, mu, sigma) / max(stats.norm.pdf(dist2, mu, sigma)), c='red')
    ax.bar(bin_centers, hist1b, width=bin_widths, align='center')
    plt.title(title_name)
    plt.savefig(os.path.join("figures", f"membrane_potential_{title_name}.png"))

    plt.clf()
    plt.close()


def load_voltage_data(file_name):
    data = np.loadtxt(file_name, delimiter=",")
    t = data[0, 1:]
    cell_id = data[1:, 0].astype(int)
    volt = data[1:, 1:]
    folder = pathlib.Path(file_name).parent
    cell_id_np = np.array([cell_id], dtype=int)
    voltage_np = np.array(volt)
    time_np = np.array(t)

    np.save(str(folder / "cell_id_simulation.npy"), cell_id_np)
    np.save(str(folder / "voltage.npy"), voltage_np)
    np.save(str(folder / "time.npy"), time_np)

    return cell_id, t, volt


def load_voltage_npy(file_name):
    f = np.load(file_name, mmap_mode="r")

    return f


def is_spiking(volt):
    if max(volt) > 0:
        return True
    else:
        return False


def statistics_spiking(volt):

    freq = get_frequency(volt)

    return freq.take(0)

def get_frequency(volt, end=10):

    s = get_spikes(volt)
    print(s)
    if len(s)!= 0:
        av = elp.statistics.mean_firing_rate(s, t_start=0 * pq.s, t_stop=end * pq.s)
    else:
        av = [0] * 1/pq.s


    return av

# already exits in root/spiking.py

def get_spikes(volt, sampling_period=25e-5):

        neo_voltage = neo.AnalogSignal(volt, units='mV', sampling_period=sampling_period * pq.s)
        spike = elp.spike_train_generation.threshold_detection(neo_voltage, threshold=-0.030 * pq.mV)

        return spike

def z_statistics_spiking(volt, mean, std):

    m_s = statistics_spiking(volt)
    print(m_s)
    z_statistics = (m_s - mean) / (math.sqrt(std))

    return z_statistics

def check_spiking(volt, mean, std, z_threshold=1):

    z_stat = abs(z_statistics_spiking(volt, mean, std))
    print(z_stat)
    if z_stat < z_threshold:
        return True
    else:
        return False

def z_statistics_membrane_potential(volt, mean, std):

    m_s, std_s = statistics_resting_membrane(volt)

    z_statistics = (m_s - mean) / (math.sqrt(std_s + std))

    return z_statistics


def check_resting(volt, mean, std, z_threshold=1):
    z_stat = abs(z_statistics_membrane_potential(volt, mean, std))

    if z_stat < z_threshold:
        return True
    else:
        return False


def statistics_resting_membrane(volt):
    v_m = np.mean(volt)
    v_std = np.std(volt)

    return v_m, v_std


def filter_volt(volt, size=1000):
    y = uniform_filter1d(volt, size=size)

    return y


def check_depolarization(volt, depolarization_block=-0.030):
    y = filter_volt(volt)
    #print("volt:")
    #print(volt)
    #print("y:")
    #print(y)
    if max(y) > depolarization_block:
        return True
    else:
        return False


def split_voltage_trace(volt, num_pieces):
    # split volt for each trial!

    split_array = np.array_split(volt, num_pieces)

    return split_array


def analyse_tuning(volt, num_trials, mean_resting_membrane, std_resting_membrane):
    trials = split_voltage_trace(volt, num_trials)

    results = list()
    for v in trials:
        #print("check_depolarization(v):")
        #print(check_depolarization(v))
        #print("len(get_spikes(v)) == 0:")
        #print(len(get_spikes(v)) == 0)
        #print(get_spikes(v))
        if not (check_depolarization(v)) and len(get_spikes(v)) == 0:

            res = check_resting(v, mean_resting_membrane, std_resting_membrane)
            results.append(res)
        else:
            results.append(False)

    return results


def analyse_tuning_spiking(volt, num_trials, mean_frq, std_frq):
    trials = split_voltage_trace(volt, num_trials)
    results = list()
    for v in trials:

        if not (check_depolarization(v)):
            res = check_spiking(v, mean=mean_frq, std=std_frq)
            print(res)
            results.append(res)
        else:
            results.append(False)

    return results
