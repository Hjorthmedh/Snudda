# 1. Load master_config into master dictionary
# 2. Identify slave configs (configs referred in master_config), and load them into master dictionary
# 3. Enrich master dictionary, e.g with additional neuron information
# 4. Additional modifications?
# 5. Export network-config.json

import os
import json
from snudda.utils import snudda_parse_path
from snudda.utils import NumpyEncoder


class ConfigParser:

    def __init__(self, path=None, snudda_data=None, parse=True):

        self.path = path
        self.network_info = dict()  # master dictionary
        self.snudda_data = snudda_data
        self.config_data = {}

        self.exclude_parse_keys = ["parameter_file", "projection_file", "rotation_file"]
        self.exclude_parse_values = ["parameters.json", "mechanisms.json", "modulation.json", "meta.json"]

        if path is not None and parse:
            self.load()

    def load(self):
        if os.path.isfile(self.path):
            with open(self.path) as f:
                self.config_data = json.load(f)

        if self.snudda_data is None:
            if "meta" in self.config_data and "snudda_data" in self.config_data["meta"]:
                self.snudda_data = self.config_data["meta"]["snudda_data"]

    def parse_config(self):

        self.config_data = self.parse_subtree(self.config_data)

    def parse_subtree(self, config_dict):

        for key, value in config_dict.items():

            if key in self.exclude_parse_keys:
                continue

            if isinstance(value, str) and value.endswith(".json"):
                putative_path = snudda_parse_path(value, snudda_data=self.snudda_data)
                if os.path.isfile(putative_path):

                    # This allows us to exclude parsing of certain json files
                    if os.path.basename(putative_path) in self.exclude_parse_values:
                        continue

                    with open(putative_path) as f:
                        sub_tree = json.load(f)

                    config_dict[key] = self.parse_subtree(sub_tree)

        return config_dict

    def write_config(self, output_path=None):

        with open(output_path, "wt") as f:
            json.dump(self.config_data, f, indent=4, cls=NumpyEncoder)

        pass


