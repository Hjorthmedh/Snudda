import numpy as np
import h5py
from collections import OrderedDict
import matplotlib.pyplot as plt
from matplotlib import cm


class PlotInput(object):

    def __init__(self, input_file):

        self.input_data = None

        if input_file:
            self.load_input(input_file)

    def load_input(self, input_file):
        self.input_data = h5py.File(input_file, "r")

    def extract_input(self, input_target):

        data = OrderedDict()

        for input_type in self.input_data["input"][input_target]:
            input_info = self.input_data["input"][input_target][input_type]

            data[input_type] = input_info["spikes"][()]

        return data

    def plot_input(self, input_target):

        data = self.extract_input(input_target)
        viridis = cm.get_cmap('viridis', len(data))
        plt_ctr = 0

        plt.figure()
        for input_ctr, input_type in enumerate(data):
            for spike_train in data[input_type]:
                idx = np.where(spike_train > 0)
                plt.scatter(spike_train[idx], plt_ctr*len(idx), color=viridis[input_ctr, :])

        plt.ion()
        plt.show()

