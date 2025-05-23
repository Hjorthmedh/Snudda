#!/usr/bin/env python3

import argparse
import os
import glob
import json
import numpy as np

from snudda import SnuddaLoad
from snudda.utils import snudda_parse_path, SnuddaLoadSimulation


class PurgeBadParameters:

    """ This code assumes you have manually moved all the figures
        corresponding to 'bad' parametersets to a 'bad' folder """

    def __init__(self, network_path, bad_figure_path=None, bad_meta_key_list=None,
                 network_simulation_path=None, snudda_data=None):

        if os.path.isfile(network_path):
            self.network_path = os.path.dirname(network_path)
            self.network_file = network_path
        elif os.path.isdir(network_path):
            self.network_path = network_path
            self.network_file = os.path.join(network_path, "network-synapses.hdf5")
        else:
            raise ValueError(f"No such file or directory {network_path}")

        self.snudda_data = snudda_data
        self.bad_keys = dict()

        if bad_meta_key_list is not None:
            self.parse_meta_key_list(bad_meta_key_list)

        if bad_figure_path is not None:
            self.get_bad_keys_in_dir(path=bad_figure_path)

        if network_simulation_path is not None:

            if "*" in network_simulation_path:
                for sim_files in glob.glob(network_simulation_path):
                    print(f"Processing {sim_files}")

                    try:
                        self.remove_depolarisation_blocked_neurons(sim_files)
                    except Exception as e:
                        import traceback
                        print(traceback.format_exc())
                        print(f"Failed to process file: {sim_files}")

            else:
                print(f"Processing {network_simulation_path}")
                self.remove_depolarisation_blocked_neurons(network_simulation_path)

    def process(self, update_files=True):
        neuron_paths = self.load_network_config(network_path=self.network_path)
        self.purge_bad_parameters(neuron_paths=neuron_paths, update_files=update_files)

    def remove_depolarisation_blocked_neurons(self, network_simulation_path):

        print(f"Processing: {network_simulation_path} ({self.network_path})")
        sls = SnuddaLoadSimulation(network_simulation_output_file=network_simulation_path,
                                   network_path=self.network_path,
                                   do_test=True)

        sl = SnuddaLoad(network_file=self.network_path)
        neuron_info = sl.data["neurons"]

        meta_data = sls.network_simulation_file["meta_data"]

        for neuron_id, parameter_key, morphology_key \
            in zip(meta_data["id"][()].copy(),
                   meta_data["parameter_key"][()].copy(),
                   meta_data["morphology_key"][()].copy()):

            if SnuddaLoad.to_str(morphology_key) != sl.data["neurons"][neuron_id]["morphology_key"] \
                    or SnuddaLoad.to_str(parameter_key) != sl.data["neurons"][neuron_id]["parameter_key"]:
                raise ValueError(f"Neurons do not match in network file {self.network_file} "
                                 f"and simulation file {network_simulation_path}.\n"
                                 f"Neuron id: {neuron_id}: "
                                 f"({sl.data['neurons'][neuron_id]['morphology_key']}, "
                                 f"{sl.data['neurons'][neuron_id]['parameter_key']}) "
                                 f"vs ({morphology_key} {parameter_key})")

        bad_cell_id = sorted(list(set([x for x, ts, te in sls.depolarisation_block])))

        for b_id in bad_cell_id:
            param_key = neuron_info[b_id]["parameter_key"]
            morph_key = neuron_info[b_id]["morphology_key"]

            print(f"Purging neuron {b_id} ({param_key}, {morph_key})")

            if param_key not in self.bad_keys:
                self.bad_keys[param_key] = []

            self.bad_keys[param_key].append(morph_key)

    def parse_meta_key_list(self, bad_meta_key_list):
        for morph_key, param_key in bad_meta_key_list:
            if param_key not in self.bad_keys:
                self.bad_keys[param_key] = []

            self.bad_keys[param_key].append(morph_key)

    def get_bad_keys_in_dir(self, path, file_extension=".png"):
        file_list = glob.glob(os.path.join(path, f"*{file_extension}"))

        for f in file_list:
            f_parts = os.path.basename(f).split("-")
            morph_key = f_parts[0]
            param_key = f_parts[1]

            if param_key not in self.bad_keys:
                self.bad_keys[param_key] = []

            self.bad_keys[param_key].append(morph_key)

    def load_network_config(self, network_path):

        network_config_path = os.path.join(network_path, "network-config.json")

        with open(network_config_path, "r") as f:
            self.network_config = json.load(f)

        if self.snudda_data is None:
            self.snudda_data = self.network_config["snudda_data"]

        neuron_paths = []

        for region_name, region_data in self.network_config["regions"].items():
            for neuron_name, neuron_data in region_data["neurons"].items():
                for name, path in neuron_data["neuron_path"].items():
                    neuron_paths.append(snudda_parse_path(path, self.snudda_data))

        return neuron_paths

    def purge_bad_parameters(self, neuron_paths, update_files=True):

        parsed_meta = []

        for neuron_path in neuron_paths:
            meta_file = os.path.join(neuron_path, "meta.json")

            if meta_file in parsed_meta:
                continue
            else:
                parsed_meta.append(meta_file)

            purge_counter = 0
            param_sets_left = 0

            with open(meta_file, "rt") as f:
                meta_data = json.load(f)

            print(f"Parsing {meta_file}")

            for param_key in list(meta_data.keys()):

                if param_key in self.bad_keys:
                    for morph_key in list(meta_data[param_key].keys()):
                        if morph_key in self.bad_keys[param_key]:
                            print(f"Removing {param_key = }, {morph_key = }")
                            del meta_data[param_key][morph_key]
                            purge_counter += 1

                if len(meta_data[param_key]) == 0:
                    print(f"--> No morphologies left for {param_key = }")
                    del meta_data[param_key]
                    meta_changed = True

            for param_key, param_data in meta_data.items():
                param_sets_left += len(param_data.keys())

            if purge_counter > 0 and update_files:
                print(f"Writing {meta_file} (purged {purge_counter} parameter sets, {param_sets_left} left)")
                with open(meta_file, "wt") as f:
                    json.dump(meta_data, f, indent=4)


def cli():

    import argparse
    parser = argparse.ArgumentParser(description="Purge the bad input parameters from meta.json")
    parser.add_argument("network_path", help="Path to network folder")
    parser.add_argument("--bad_figure_path", help="Path to 'bad' figure folder", default=None, type=str)
    parser.add_argument("--bad_key_list_file", help="Path to 'bad' meta key list file", default=None, type=str)
    parser.add_argument("--snudda_data", type=str, default=None)
    parser.add_argument("--network_simulation_path", type=str, default=None)
    parser.add_argument("--mock_run", action="store_true")
    args = parser.parse_args()

    if args.bad_key_list_file is not None:
        bad_key_list =np.genfromtxt(args.bad_key_list_file, delimiter=",", dtype=None, encoding="utf-8", skip_header=0)
    else:
        bad_key_list = None

    pbp = PurgeBadParameters(network_path=args.network_path,
                             bad_figure_path=args.bad_figure_path,
                             bad_meta_key_list=bad_key_list,
                             network_simulation_path=args.network_simulation_path,
                             snudda_data=args.snudda_data)

    pbp.process(update_files=not args.mock_run)


if __name__ == "__main__":
    cli()
