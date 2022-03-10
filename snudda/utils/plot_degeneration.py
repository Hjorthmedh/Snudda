import h5py
import numpy as np
import matplotlib.pyplot as plt

from snudda import SnuddaLoad


class PlotDegeneration:

    def __init__(self, original_network_file, original_input_file, degenerated_network_file, degenerated_input_file):

        self.original_network_file = original_network_file
        self.original_input_file = original_input_file

        self.degenerated_network_file = degenerated_network_file
        self.degenerated_input_file = degenerated_input_file

        self.original_network_loader = SnuddaLoad(self.original_network_file)
        self.degenerated_network_loader = SnuddaLoad(self.degenerated_network_file)

        self.original_data = self.original_network_loader.data
        self.degenerated_data = self.degenerated_network_loader.data

        self.original_input = h5py.File(self.original_input_file, 'r')
        self.degenerated_input = h5py.File(self.degenerated_input_file, 'r')

