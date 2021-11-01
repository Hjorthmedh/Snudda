import os
from glob import glob
import json
from collections import OrderedDict


class CreateParameterMorphologyInputMap:

    def __init__(self):


        pass

    # TODO: Parse through all the neuron types, and neuron directories

    def parse_neuron(self, neuron_path):
        param_file = os.path.join(neuron_path, "parameters.json")
        mapping_file = os.path.join(neuron_path, "morphology_input_mapping.json")

        mapping_data = OrderedDict()

        with open(param_file, "r") as f:
            param_data = json.load(f, object_pairs_hook=OrderedDict)

        # First add "key" to each parameter set
        for idx, p in param_data:
            if "key" not in p:
                param_data[idx]["key"] = str(idx)

        for p in param_data:
            key = p["key"]
            if key not in mapping_data:
                mapping_data[key] = OrderedDict()

            for m in p["morphology"]:

                mapping_data[key][m] = OrderedDict()
                mapping_data[key][m]["Cortical"] = 100
                mapping_data[key][m]["Thalamic"] = 100

        # TODO: We should remove the old morphology mapping from the param_file

        with open(mapping_file, "w") as fo:
            json.dump(mapping_data, fo)
