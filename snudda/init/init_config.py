# 1. Load master_config into master dictionary
# 2. Identify slave configs (configs referred in master_config), and load them into master dictionary
# 3. Enrich master dictionary, e.g with additional neuron information
# 4. Additional modifications?
# 5. Export network-config.json

import os
import json
from copy import deepcopy
from snudda.utils import snudda_parse_path
from snudda.utils import NumpyEncoder


class ConfigParser:

    def __init__(self, path=None, snudda_data=None, parse=True):

        self.path = path
        self.network_info = dict()  # master dictionary
        self.snudda_data = snudda_data
        self.config_data = {}
        self.rng = None

        self.exclude_parse_keys = ["parameter_file", "projection_file", "rotation_file", "rotation_field_file"]
        self.exclude_parse_values = ["parameters.json", "mechanisms.json", "modulation.json", "meta.json"]

        if path is not None and parse:
            self.load()

    def load(self):
        print(f"Loading {self.path}")

        if os.path.isfile(self.path):
            with open(self.path) as f:
                self.config_data = json.load(f)
        else:
            print(f"File not found!")

        if self.snudda_data is None:
            if "meta" in self.config_data and "snudda_data" in self.config_data["meta"]:
                self.snudda_data = self.config_data["meta"]["snudda_data"]

    def parse_config(self):

        self.config_data = self.parse_subtree(self.config_data)

        self.setup_random_seeds()

    def substitute_json(self, putative_file, key=None):

        if not (isinstance(putative_file, str) and putative_file.endswith(".json")):
            return putative_file

        putative_path = snudda_parse_path(putative_file, snudda_data=self.snudda_data)
        if not os.path.isfile(putative_path):
            raise ValueError(f"File not found {putative_path}")

        # This allows us to exclude parsing of certain json files
        if os.path.basename(putative_path) in self.exclude_parse_values:
            return putative_file

        with open(putative_path) as f:
            print(f"Loading {putative_path}")
            sub_tree = json.load(f)

        sub_tree = self.parse_subtree(sub_tree)

        if key in sub_tree:
            return sub_tree[key]

        return sub_tree

    def parse_subtree(self, config_dict):

        updated_config = deepcopy(config_dict)

        for key, value in config_dict.items():

            if key in self.exclude_parse_keys:
                continue

            # Keys beginning with ! are excluded
            elif key[0] == "!":
                del updated_config[key]

            elif isinstance(value, dict):
                updated_config[key] = self.parse_subtree(value)

            elif isinstance(value, list):
                updated_config[key] = [self.substitute_json(x) for x in value]

            else:
                updated_config[key] = self.substitute_json(value, key)

        return updated_config

    def setup_random_seeds(self):

        if "random_seed" in self.config_data:
            if "masterseed" in self.config_data["random_seed"]:
                master_seed = self.config_data["random_seed"]["masterseed"]
            else:
                master_seed = None

        from snudda.init import SnuddaInit

        rand_seed_dict, init_rng = SnuddaInit.setup_random_seeds(random_seed=master_seed)

        self.config_data["random_seed"] = rand_seed_dict
        self.rng = init_rng

    def write_config(self, output_path=None):

        with open(output_path, "wt") as f:
            json.dump(self.config_data, f, indent=4, cls=NumpyEncoder)


if __name__ == "__main__":

    config_path = "/home/hjorth/HBP/BasalGangliaData/data/connectivity/network.json"

    conf = ConfigParser(path=config_path, snudda_data="/home/hjorth/HBP/BasalGangliaData/data")
    conf.parse_config()
    conf.write_config("test.json")

    # import pdb
    # pdb.set_trace()