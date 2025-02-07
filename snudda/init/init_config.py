# 1. Load master_config into master dictionary
# 2. Identify slave configs (configs referred in master_config), and load them into master dictionary
# 3. Enrich master dictionary, e.g with additional neuron information
# 4. Additional modifications?
# 5. Export network-config.json

import os
import json
from copy import deepcopy
from collections import ChainMap
from snudda.utils import snudda_parse_path
from snudda.utils import NumpyEncoder


class ConfigParser:

    def __init__(self, config_file=None, snudda_data=None, parse=True):

        self.config_file = config_file
        self.network_info = dict()  # master dictionary
        self.snudda_data = snudda_data
        self.config_data = {}
        self.rng = None

        # By default, JSON files which are included in the config file are parsed and inserted in place
        # if you instead want to keep the reference to the file, then add it below:
        self.exclude_parse_keys = ["parameter_file", "projection_file", "rotation_file", "rotation_field_file",
                                   "reaction_diffusion", "modulation", "extracellular_config"]
        self.exclude_parse_values = ["parameters.json", "mechanisms.json", "modulation.json", "meta.json"]

        if config_file is not None and parse:
            self.load()

    def load(self):
        print(f"Loading {self.config_file}")

        if os.path.isfile(self.config_file):
            with open(self.config_file) as f:
                self.config_data = json.load(f)
        else:
            raise ValueError(f"File not found: {self.config_file}")

        if self.snudda_data is None:
            if "snudda_data" in self.config_data:
                self.snudda_data = os.path.abspath(self.config_data["snudda_data"])
                self.config_data["snudda_data"] = self.snudda_data  # Make sure abspath is stored
            else:
                snudda_data = ""
        else:
            if "snudda_data" in self.config_data:
                print(f"Warning, overriding snudda_data with {os.path.abspath(self.snudda_data)}")

            self.config_data["snudda_data"] = os.path.abspath(self.snudda_data)

        if not self.snudda_data:
            print(f"Warning, snudda_data is not set!")

    def parse_config(self):

        self.config_data = self.parse_subtree(self.config_data, parent_file=self.config_file)

        self.setup_random_seeds()

    def replace_network_path(self, network_path):
        print(f"Setting {network_path = }")
        self.config_data["network_path"] = network_path

    def skip_item(self, item):

        return isinstance(item, str) and len(item) > 0 and item[0] == '!'

    def get_putative_path(self, putative_file, parent_file=None):

        # This allows us to exclude parsing of certain json files
        if os.path.basename(putative_file) in self.exclude_parse_values:
            return putative_file

        parent_dir = os.path.dirname(parent_file)
        base_dir = os.path.dirname(self.config_file)

        # Version of file name that we check to see if we can find the file
        putative_files = [putative_file,
                          os.path.join(parent_dir, putative_file),
                          os.path.join(base_dir, putative_file),
                          snudda_parse_path(putative_file, snudda_data=self.snudda_data)]
        putative_path = None

        for pf in putative_files:
            if os.path.isfile(pf):
                putative_path = pf
                break

        if putative_path is None:
            raise ValueError(f"File not found {putative_file}")

        return putative_path

    def substitute_json(self, putative_file, key=None, parent_file=None):

        if not (isinstance(putative_file, str) and putative_file.endswith(".json")):
            return putative_file

        # First we look for files in the same directory as the parent directory
        # after that we do SNUDDA_DATA substitution if no match

        putative_path = self.get_putative_path(putative_file=putative_file,
                                               parent_file=parent_file)

        with open(putative_path) as f:
            print(f"Loading {putative_path}")
            sub_tree = json.load(f)

        sub_tree = self.parse_subtree(sub_tree, parent_file=putative_path)

        if key in sub_tree:
            return sub_tree[key]

        return sub_tree

    def parse_subtree(self, config_dict, parent_file=None):

        updated_config = deepcopy(config_dict)

        for key, value in config_dict.items():

            if key in self.exclude_parse_keys:
                continue

            # Keys beginning with ! are excluded
            elif key[0] == "!":
                del updated_config[key]

            elif isinstance(value, dict):
                updated_config[key] = self.parse_subtree(value, parent_file=parent_file)

            elif isinstance(value, list):
                updated_config[key] = [self.substitute_json(x, parent_file=parent_file) for x in value if not self.skip_item(x)]

                # If the first element is a dict, assume all are dict and merge them
                if isinstance(updated_config[key][0], dict):
                    updated_config[key] = dict(ChainMap(*reversed(updated_config[key])))

            else:
                if not self.skip_item(item=value):
                    updated_config[key] = self.substitute_json(value, key, parent_file=parent_file)
                else:
                    del updated_config[key]

        return updated_config

    def setup_random_seeds(self):

        master_seed = None

        if "random_seed" in self.config_data and "master_seed" in self.config_data["random_seed"]:
            master_seed = self.config_data["random_seed"]["master_seed"]

        from snudda.init import SnuddaInit

        rand_seed_dict, init_rng = SnuddaInit.setup_random_seeds(random_seed=master_seed)

        self.config_data["random_seed"] = rand_seed_dict
        self.rng = init_rng

    def write_config(self, output_path=None):

        if output_path is None:
            output_path = os.path.join(self.config_data["network_path"], "network-config.json")

        print(f"Writing to file {output_path}")
        network_dir = os.path.dirname(output_path)
        if network_dir and not os.path.isdir(network_dir):
            os.makedirs(network_dir)

        with open(output_path, "wt") as f:
            json.dump(self.config_data, f, indent=4, cls=NumpyEncoder)


if __name__ == "__main__":

    config_path = "/home/hjorth/HBP/Snudda/new_config/network.json"

    conf = ConfigParser(config_file=config_path, snudda_data="/home/hjorth/HBP/BasalGangliaData/data")
    conf.parse_config()
    conf.write_config()

    # import pdb
    # pdb.set_trace()