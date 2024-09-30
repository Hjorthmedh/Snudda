import itertools
import os
import numpy as np

from snudda.utils.load import SnuddaLoad
from snudda.utils.load_network_simulation import SnuddaLoadSimulation
from snudda.utils.export_connection_matrix import SnuddaExportConnectionMatrix
from collections import OrderedDict


class SnuddaAnalyseTopology:

    def __init__(self, network_file):

        self.network_file = network_file
        self.snudda_load = SnuddaLoad(network_file=network_file)
        self.simplex_data = dict()

        secm = SnuddaExportConnectionMatrix(in_file=self.network_file, out_file="dummy_file", save_on_init=False)
        self.connection_matrix = secm.create_connection_matrix()

    def load_simplex_file(self, simplex_file_name):
        data = np.loadtxt(simplex_file_name, dtype=int)
        simplex_dimension = data.shape[1] - 1
        self.simplex_data[simplex_dimension] = data

        print(f"Loaded simplex data of dimension {simplex_dimension} from {simplex_file_name}")

    def verify_source_sink_order(self):

        for dim in self.simplex_data.keys():
            print(f"Verifying dimension {dim}")
            ctr = 0
            for simplex in self.simplex_data[dim]:
                sub_connection_matrix = self.connection_matrix[simplex, :][:, simplex]
                try:
                    assert (sub_connection_matrix[0, 1:] > 0).all(), f"First element not source in simplex {simplex}"
                    assert (sub_connection_matrix[:-1, -1] > 0).all(), f"Last element not sink in simplex {simplex}"
                except:
                    import traceback
                    print(traceback.format_exc())
                    import pdb
                    pdb.set_trace()
                ctr += 1
            print(f"Verified {ctr} simplices")

    def get_multiplicity(self):

        multiplicity = OrderedDict()

        for dim in self.simplex_data.keys():
            mult_dict = dict()
            for simplex in self.simplex_data[dim]:
                idx = tuple(np.sort(simplex))
                if idx not in mult_dict:
                    mult_dict[idx] = 1
                else:
                    mult_dict[idx] += 1

            multiplicity[dim] = mult_dict

        return multiplicity

    def get_fixed_multiplicity(self):

        multiplicity = OrderedDict()

        for dim in self.simplex_data.keys():
            mult_dict = dict()
            for simplex in self.simplex_data[dim]:
                simplex_copy = simplex.copy()
                simplex_copy[1:-1].sort()  # Sorts the part of the array in place, wow!
                idx = tuple(simplex_copy)
                if idx not in mult_dict:
                    mult_dict[idx] = 1
                else:
                    mult_dict[idx] += 1

            multiplicity[dim] = mult_dict

        return multiplicity

    def get_fixed_multiplicity_ALT(self):

        multiplicity = OrderedDict()

        for dim in self.simplex_data.keys():
            mult_dict = dict()
            for simplex in self.simplex_data[dim]:
                simplex_copy = simplex.copy()
                simplex_copy[1:-1].sort()  # Sorts the part of the array in place, wow!
                idx = tuple(simplex_copy)
                if idx not in mult_dict:
                    mult_dict[idx] = [1, [simplex]]
                else:
                    mult_dict[idx][0] += 1
                    mult_dict[idx][1].append(simplex)

            multiplicity[dim] = mult_dict

        return multiplicity

    def filter_multiplicity(self, multiplicity, dimension, neuron_type_list, multiplicity_requirement=None):

        """

        Args:
            multiplicity (OrderedDict) : Multiplicity dictionary
            dimension (int) : Dimension of data to look at
            neuron_type_list (list) : Filtering information, e.g. [('FS', 2), ('dSPN', 1)] means that
                                      only cliques that contain exactly 2 FS and 1 dSPN are kept
            multiplicity_requirement (int) : Multiplicity of the simplex kept, default = None (no filtering)

        Returns:
            filtered_multiplicity (OrderedDict) : Only returns dimension data
                                  """

        filtered_multiplicity = OrderedDict()

        for neurons_key, mult in multiplicity[dimension].items():
            neuron_types = [self.snudda_load.data["neurons"][x]["type"] for x in neurons_key]

            keep_simplex = True if multiplicity_requirement is None else mult == multiplicity_requirement

            if keep_simplex is False:
                continue

            for neuron_type, neuron_type_number in neuron_type_list:
                if np.sum([neuron_type == x for x in neuron_types]) != neuron_type_number:
                    keep_simplex = False

            if keep_simplex:
                filtered_multiplicity[neurons_key] = mult

        return filtered_multiplicity

    def print_multiplicity(self, fixed=True):

        if fixed:
            # Fix source and sink in the ordering
            multiplicity = self.get_fixed_multiplicity()
        else:
            multiplicity = self.get_multiplicity()

        for dim in multiplicity.keys():
            print(f"-- Analysing dimension {dim}")
            #import pdb
            #pdb.set_trace()
            repeats = np.bincount(np.array([x for x in multiplicity[dim].values()]))
            for reps, count in enumerate(repeats):
                if count > 0:
                    print(f"Multiplicity {reps} for {count} simplices")

            print("")
        #import pdb
        #pdb.set_trace()

    def get_clique_composition_combination(self, dimension, neuron_types=None):

        if neuron_types is None:
            neuron_types = self.snudda_load.get_neuron_types()

        clique_combinations = [x for x in itertools.combinations_with_replacement(neuron_types, dimension+1)]

        return clique_combinations
        
    def get_clique_composition_disposition(self, dimension, neuron_types=None):

        if neuron_types is None:
            neuron_types = self.snudda_load.get_neuron_types()

        clique_combinations = [x for x in itertools.product(neuron_types, repeat=dimension+1)]

        return clique_combinations

    def get_clique_neuron_type_composition_statistics(self, multiplicity, dimension):

        count = OrderedDict()
        neuron_type_lookup = dict()

        for neuron_id, neuron_type in zip([x["neuron_id"] for x in self.snudda_load.data["neurons"]],
                                          [x["type"] for x in self.snudda_load.data["neurons"]]):
            neuron_type_lookup[neuron_id] = neuron_type

        for neuron_keys, mult in multiplicity[dimension].items():
            neuron_type_list = tuple(sorted([neuron_type_lookup[nk] for nk in neuron_keys]))
            if neuron_type_list not in count:
                count[neuron_type_list] = 1
            else:
                count[neuron_type_list] += 1

        return count


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Load snudda network file (hdf5)")
    parser.add_argument("network_file", help="Network file (hdf5)", type=str)
    parser.add_argument("simplex_files", nargs="*", action="append", type=str)

    args = parser.parse_args()

    at = SnuddaAnalyseTopology(network_file=args.network_file)
    for simplex in args.simplex_files[0]:
        at.load_simplex_file(simplex)

    at.verify_source_sink_order()
    at.print_multiplicity(fixed=True)

    import pdb
    pdb.set_trace()
