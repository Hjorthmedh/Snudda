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

    # TODO: Add labels to spike rasters
    def plot_input(self, input_target):

        if type(input_target) != list:
            input_target = [input_target]

        # Make sure each target is a str
        input_target = [str(x) for x in input_target]
        colours = cm.get_cmap('tab20', len(input_target) * 2)

        synapse_ctr = 0
        input_ctr = 0
        plt.figure()

        for it in input_target:

            data = self.extract_input(it)

            for input_type in data:
                for spike_train in data[input_type]:
                    idx = np.where(spike_train > 0)[0]
                    plt.scatter(spike_train[idx], synapse_ctr * np.ones((len(idx),)), color=colours(input_ctr))
                    synapse_ctr += 1

                input_ctr += 1
                    
        plt.ion()
        plt.show()

