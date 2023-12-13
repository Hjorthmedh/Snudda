import os
import numpy as np


def get_spn_fs_freq():

    data_dir = os.path.dirname(__file__)
    data_file = os.path.join(data_dir, "Kim et al., 2019 - Henry Yin - Extended2a.csv")

    data = np.genfromtxt(data_file, delimiter=",", filling_values=np.nan)
    spn_freq, fs_freq = data[:, 0], data[:, 1]
    fs_freq = fs_freq[~np.isnan(fs_freq)]

    return spn_freq, fs_freq


def resample_spn_freq(size, rng=None):

    spn_freq, _ = get_spn_fs_freq()

    if rng is None:
        rng = np.random.default_rng()

    resampled_freq = rng.choice(a=spn_freq, size=size, replace=True)

    return resampled_freq


def resample_fs_freq(size, rng=None):

    _, fs_freq = get_spn_fs_freq()

    if rng is None:
        rng = np.random.default_rng()

    resampled_freq = rng.choice(a=fs_freq, size=size, replace=True)

    return resampled_freq


def plot_freq_hist(data, title=None, bins=None):

    import matplotlib.pyplot as plt

    plt.figure()
    plt.hist(data, bins=bins)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Frequency of occurence")
    if title:
        plt.title(title)

    plt.ion()
    plt.show()
    plt.pause(0.1)


if __name__ == "__main__":

    spn_freq, _ = get_spn_fs_freq()
    resampled_spn_freq = resample_spn_freq((10000,1))

    plot_freq_hist(spn_freq, "SPN frequency", bins=100)
    plot_freq_hist(resampled_spn_freq, "Resampled SPN frequency", bins=100)

    input("Press a key to continue")

