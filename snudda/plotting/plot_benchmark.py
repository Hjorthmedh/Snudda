import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import os
import json


class PlotBenchmark:

    def __init__(self, network_path_list=None, number_of_nodes=None, show_only_latest=True):

        if network_path_list:
            pass

    def load_benchmark(self, network_path):
        file_path = os.path.join(network_path, "benchmark_log.json")
        data = json.load(file_path)


