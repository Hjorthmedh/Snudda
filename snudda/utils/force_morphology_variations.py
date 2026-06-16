import os
import h5py
import shutil
import glob
import json
import numpy as np

from snudda.utils.load import SnuddaLoad
from snudda.utils import snudda_parse_path
from snudda.utils.snudda_path import snudda_simplify_path

class ForceMorphologyVariations(object):

    def __init__(self, network_path, original_network_path):

        self.verbose = True

        self.network_path = network_path
        self.original_network_path = original_network_path
        self.original_snudda_data = None

        self.original_network = None
        self.updated_network = None
        self.get_data_from_updated_config()


    def max_len(self, my_list):
        return np.max([len(x) for x in my_list])

    def swap(self):

        # Copy the original network file
        source_file = os.path.join(self.original_network_path, "network-neuron-positions.hdf5")
        updated_file = os.path.join(self.network_path, "network-neuron-positions.hdf5")
        shutil.copyfile(source_file, updated_file)

        print(f"Template: {source_file}\nUpdated place file: {updated_file}")

        self.original_network = h5py.File(source_file, "r")
        self.updated_network = h5py.File(updated_file, "r+")

        self.original_snudda_data = SnuddaLoad.to_str(self.original_network["meta/snudda_data"][()])

        # The new morphology paths might be longer, so we need to do this carefully.

        for idx, neuron_id in enumerate(self.original_network["network/neurons/neuron_id"][()]):
            param_key, morph_key, neuron_path, new_morph_name, param_id, morph_id = self.find_morpology(neuron_id)

            self.updated_network[f"network/neurons/parameter_key"][idx] = param_key
            self.updated_network[f"network/neurons/morphology_key"][idx] = morph_key
            self.updated_network[f"network/neurons/neuron_path"][idx] = neuron_path
            self.updated_network[f"network/neurons/morphology"][idx] = new_morph_name


            if SnuddaLoad.to_str(self.updated_network[f"network/neurons/neuron_path"][idx]) != neuron_path:
                print(f"Warning truncated names for neuron_id {neuron_id}")
                print(SnuddaLoad.to_str(self.updated_network[f"network/neurons/neuron_path"][idx]))
                print(neuron_path)
                import pdb
                pdb.set_trace()


            # TODO: Make sure that axon_density, and friends,... reaction_diffusion_file etc
            #       are not set, because if so they might differ between WT and PD/SZ networks
            #

        # Write the new data to file
        self.updated_network["meta/snudda_data"][()] = self.updated_snudda_data
        self.updated_network["meta/config_file"][()] = self.updated_network_config_file
        self.updated_network["meta/config"][()] = json.dumps(self.updated_config)

        self.original_network.close()
        self.original_network = None

        self.updated_network.close()
        self.updated_network = None



    def get_data_from_updated_config(self):

        config_file = os.path.join(self.network_path, "network-config.json")
        with open(config_file, "r") as f:
            config = json.load(f)

        self.updated_snudda_data = config["snudda_data"]
        self.updated_network_config_file = config_file
        self.updated_config = config

    # 1. Do snudda place on the network_path
    # 2. Replace the morphology key and parameter key, to match the original variations of the morphologies
    # 3. Perform touch detection and pruning (bonus)

    def find_morpology(self, neuron_id):

        """ Given neuron_id it looks up what was old morphology, then looks for a new morphology with the same name
            and finds the corresponding morphology_key. It also tries to pick a parameter_key """

        orig_morph_key = SnuddaLoad.to_str(self.original_network["network/neurons"]["morphology_key"][neuron_id])
        orig_param_key = SnuddaLoad.to_str(self.original_network["network/neurons"]["parameter_key"][neuron_id])

        # We assume the new morphology is in the same relative path, but using a different SNUDDA_DATA
        orig_simple_path = SnuddaLoad.to_str(self.original_network["network/neurons"]["neuron_path"][neuron_id])
        orig_neuron_path = snudda_parse_path(orig_simple_path, os.path.realpath(self.original_snudda_data))
        new_neuron_path = orig_neuron_path.replace(os.path.realpath(self.original_snudda_data),
                                                   os.path.realpath(self.updated_snudda_data))

        if not os.path.isdir(new_neuron_path):
            # We did not find an exact match, see if we find a dir that contains a substring
            last_dir = os.path.basename(os.path.realpath(new_neuron_path))
            parent_dir = os.path.dirname(os.path.realpath(new_neuron_path))

            # Find all the directories in that, and see if any of them contains 'last_dir'
            potential_dir = glob.glob(os.path.join(parent_dir, '*'))
            candidate_dir = [os.path.join(last_dir, x) for x in potential_dir if
                             last_dir in os.path.basename(x)]

            if len(candidate_dir) == 0:
                # Let's do one last hail mary try
                partial_last_dir = last_dir.rsplit("-", 1)[
                    0]  # Here we remove everything after last -, and treat that as the stub
                assert len(partial_last_dir) > 0

                candidate_dir = [os.path.join(last_dir, x) for x in potential_dir if
                                 partial_last_dir in os.path.basename(x)]

            if len(candidate_dir) != 1:
                raise ValueError(
                    f"Unable to find {new_neuron_path}, also tried looking for similar folders in {parent_dir}")

            else:
                print(
                    f"Did not find exact directory match, using {candidate_dir[0]} instead of {new_neuron_path}")
                new_neuron_path = candidate_dir[0]

        if orig_morph_key is None or orig_morph_key == '':
            # Only a single morphology

            original_morphology_id = 0  # self.old_data["neurons"][neuron_id]["morphology_id"]
            original_parameter_id = 0  # self.old_data["neurons"][neuron_id]["parameter_id"]
            new_morph_name = os.path.basename(new_neuron_path)

            return '', '', new_neuron_path, new_morph_name, original_parameter_id, original_morphology_id

        # Here we assume there is a meta.json file
        with open(os.path.join(orig_neuron_path, "meta.json"), "r") as f:
            orig_meta_info = json.load(f)

        with open(os.path.join(new_neuron_path, "meta.json"), "r") as f:
            new_meta_info = json.load(f)

        # Find the parameter_key and morphology_key corresponding to the morphology used. Note that there might be
        # multiple parameter_key:s valid for a morphology, so then we randomly pick one of those

        orig_morph_name = orig_meta_info[orig_param_key][orig_morph_key]["morphology"]

        possible_keys = []

        for param_key, param_data in new_meta_info.items():
            for morph_key, morph_data in param_data.items():
                morph_name = morph_data["morphology"]
                if self.verbose:
                    print(f"-- Comparing {orig_morph_name} {morph_name}")

                if orig_morph_name == morph_name:
                    possible_keys.append((param_key, morph_key, morph_name))
                    if self.verbose:
                        print(f"Matching {orig_morph_name} with {morph_name}")

                # We also need to be able to handle Treem adding tags to the new filename
                elif os.path.splitext(os.path.basename(orig_morph_name))[0] in morph_name:

                    # Hack, so we do not match a non-var morphology to a var morphology
                    # We allow non-var to match with var0 (which is alias for non-var)
                    if "var" not in os.path.splitext(os.path.basename(orig_morph_name))[0] \
                            and "var" in os.path.splitext(os.path.basename(morph_name))[0] \
                            and not "var0" in os.path.splitext(os.path.basename(morph_name))[0]:
                        print(f"Skipping match : {orig_morph_name} -- {morph_name}")
                        continue

                    possible_keys.append((param_key, morph_key, morph_name))
                    if self.verbose:
                        print(f"Matching (close) {orig_morph_name} with {morph_name}")

                elif "var0" in os.path.splitext(os.path.basename(orig_morph_name))[0]:
                    # Also need to check if var0 has a corresponding morphology without var0 in name
                    new_candidate = orig_morph_name.replace("-var0", "")
                    if new_candidate == morph_name:
                        possible_keys.append((param_key, morph_key, morph_name))
                        if self.verbose:
                            print(f"Matching (close2) {orig_morph_name} with {morph_name}")

        if len(possible_keys) == 0:
            raise ValueError(f"No morphology matching for {orig_morph_name}, "
                             f"unable to pick parameter_key and morphology_key")

        # Pick one of the key pairs
        idx = np.random.randint(low=0, high=len(possible_keys))
        new_param_key, new_morph_key, new_morph_name = possible_keys[idx]
        new_morph_name = snudda_simplify_path(os.path.join(new_neuron_path, "morphology", new_morph_name),
                                              self.updated_snudda_data)

        # We also need parameter_id and morphology_id
        parameter_id = np.where([x == new_param_key for x in new_meta_info.keys()])[0][0]
        morphology_id = np.where([x == new_morph_key for x in new_meta_info[new_param_key].keys()])[0][0]

        return new_param_key, new_morph_key, new_neuron_path, new_morph_name, parameter_id, morphology_id


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("network_path", type=str)
    parser.add_argument("original_network_path", type=str)
    args = parser.parse_args()

    fmv = ForceMorphologyVariations(args.network_path, args.original_network_path)
    fmv.swap()

