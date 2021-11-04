import json
import numpy as np
from snudda.utils.load import SnuddaLoad
from snudda.utils.numpy_encoder import NumpyEncoder
from collections import OrderedDict


class FakeLoad(SnuddaLoad):

    def __init__(self, network_file=None, load_synapses=True, verbose=False):

        super(FakeLoad, self).__init__(network_file=network_file, load_synapses=load_synapses, verbose=verbose)

    def import_hdf5(self, network_file):

        self.load_hdf5(network_file)

    def export_json(self, json_file_name):
        if self.data:
            data_copy = self.data.copy()

            # Fix for connectivityDistributions
            data_copy["connectivityDistributions"] = OrderedDict()
            for pre_id, post_id in self.data["connectivityDistributions"].keys():
                data_copy["connectivityDistributions"][f"{pre_id}$${post_id}"] =\
                    self.data["connectivityDistributions"][pre_id,post_id]

            data_copy["connectivityDistributions"]

            with open(json_file_name, "w") as f:
                json.dump(data_copy, f, indent=4, cls=NumpyEncoder)
        else:
            print("You need to load data before exporting.")

    def import_json(self, json_file_name):

        with open(json_file_name, "r") as f:
            self.data = json.load(f, object_pairs_hook=OrderedDict)

        con_dist_copy = self.data["connectivityDistributions"].copy()
        self.data["connectivityDistributions"] = OrderedDict()

        for keys in con_dist_copy:
            (pre_type, post_type) = keys.split("$$")
            self.data["connectivityDistributions"][pre_type, post_type] = con_dist_copy[keys]

        # Also fix format
        format_fix = ["neuronPositions", "synapseCoords", "neuronID", "synapses", "gapJunctions",
                      "simulationOrigo", "populationUnit"]
        for fix_data in format_fix:
            self.data[fix_data] = np.array(self.data[fix_data])

        for idx in range(len(self.data["neurons"])):
            self.data["neurons"][idx]["position"] = np.array(self.data["neurons"][idx]["position"])
            self.data["neurons"][idx]["rotation"] = np.array(self.data["neurons"][idx]["rotation"])

        if "config" in self.data:
            self.config = json.loads(self.data["config"], object_pairs_hook=OrderedDict)

