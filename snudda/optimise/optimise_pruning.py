import os
import json
import collections
import numpy as np
import h5py

from snudda.utils import SnuddaLoad
from snudda.detect import SnuddaPrune

from scipy.stats import binomtest
from scipy.spatial import distance_matrix
from scipy.optimize import differential_evolution

import uuid


class OptimisePruning:

    """ Optimises pruning parameters. First it creates a relatively small network and does touch detection on it.
        After that it does multiple pruning attempts (while keeping the files original detection files). """

    def __init__(self, network_path, pop_size=10, epochs=10):

        self.network_path = network_path
        self.pop_size = pop_size
        self.epochs = epochs

        self.prune = SnuddaPrune(network_path=network_path, keep_files=True)

        self.merge_files_syn = None
        self.merge_neuron_range_syn = None
        self.merge_syn_ctr = None
        self.merge_files_gj = None
        self.merge_neuron_range_gj = None
        self.merge_gj_ctr = None

        self.optimisation_info = dict()

    def merge_putative_synapses(self):

        merge_info = self.prune.get_merge_info()

        if merge_info:
            merge_files_syn, merge_neuron_range_syn, merge_syn_ctr, \
            merge_files_gj, merge_neuron_range_gj, merge_gj_ctr = merge_info
        else:
            # From the hyper voxels gather all synapses (and gap junctions) belonging to specific neurons
            merge_files_syn, merge_neuron_range_syn, merge_syn_ctr, \
            merge_files_gj, merge_neuron_range_gj, merge_gj_ctr = self.prune.gather_neuron_synapses()

            self.prune.save_merge_info(merge_files_syn=merge_files_syn,
                                       merge_neuron_range_syn=merge_neuron_range_syn,
                                       merge_syn_ctr=merge_syn_ctr,
                                       merge_files_gj=merge_files_gj,
                                       merge_neuron_range_gj=merge_neuron_range_gj,
                                       merge_gj_ctr=merge_gj_ctr)

        self.merge_files_syn = merge_files_syn
        self.merge_neuron_range_syn = merge_neuron_range_syn
        self.merge_syn_ctr = merge_syn_ctr
        self.merge_files_gj = merge_files_gj
        self.merge_neuron_range_gj = merge_neuron_range_gj
        self.merge_gj_ctr = merge_gj_ctr

    def prune_synapses(self, pre_type, post_type, con_type, pruning_parameters, output_file):

        """ Prunes the network keeping only synapses between pre_type and post_type, using pruning_parameters

        Args:
            pre_type (str) : presynaptic neuron type (e.g. "FS")
            post_type (str) : postsynaptic neuron type (e.g. "dSPN")
            con_type (str) : type of connection (e.g. "GABA")
            pruning_parameters (dict) :
        """

        self.prune.load_pruning_information()

        orig_connectivity_distributions = json.loads(self.prune.hist_file["meta/connectivityDistributions"][()],
                                                     object_pairs_hook=collections.OrderedDict)

        # Next overwrite the pruning parameters
        pre_type_id = self.prune.type_id_lookup[pre_type]
        post_type_id = self.prune.type_id_lookup[post_type]
        orig_key = f"{pre_type}$${post_type}"
        synapse_type_id = orig_connectivity_distributions[orig_key][con_type]["channelModelID"]

        pruning = self.prune.complete_pruning_info(pruning_parameters)
        pruning_other = None

        self.prune.connectivity_distributions = dict([])
        self.prune.connectivity_distributions[pre_type_id, post_type_id, synapse_type_id] = (pruning, pruning_other)

        assert len(self.merge_files_syn) == 1, f"merge_Files_syn should be a list with one file only"

        with h5py.File(self.merge_files_syn[0], "r") as f_syn:
            print(f"Writing file {output_file}")
            # We need to make sure that we keep the pruned data files separate
            self.prune.prune_synapses(synapse_file=f_syn,
                                      output_filename=output_file,
                                      row_range=None,
                                      close_out_file=True,
                                      close_input_file=True,
                                      merge_data_type="synapses")

            if not os.path.exists(output_file):
                print(f"Output file missing {output_file}")
                import pdb
                pdb.set_trace()

    def evaluate_fitness(self, pre_type, post_type, output_file, experimental_data):

        """

            Args:
                pre_type
                post_type
                output_file: path to output file from prune
                experiment_data: [(bin start, bin end, P)]

        """

        print(f"Opening {output_file}")
        snudda_load = SnuddaLoad(network_file=output_file)
        snudda_data = snudda_load.data

        connection_matrix = np.zeros((snudda_data["nNeurons"], snudda_data["nNeurons"]))

        pre_id = snudda_load.get_neuron_id_of_type(neuron_type=pre_type)
        post_id = snudda_load.get_neuron_id_of_type(neuron_type=post_type)

        pre_mask = np.zeros((snudda_data["nNeurons"],), dtype=bool)
        post_mask = np.zeros((snudda_data["nNeurons"],), dtype=bool)

        pre_mask[pre_id] = True
        post_mask[post_id] = True

        for row in snudda_data["synapses"]:
            if pre_mask[row[0]] and post_mask[row[1]]:
                # Only include connections between the right pre and post types
                connection_matrix[row[0], row[1]] += 1

        pos = snudda_data["neuronPositions"]
        dist_matrix = distance_matrix(pos, pos)

        n_connected = np.zeros((len(experimental_data),), dtype=int)
        n_total = np.zeros((len(experimental_data),), dtype=int)

        for dist, con in zip(dist_matrix.flatten(), connection_matrix.flatten()):
            for idx, (bin_start, bin_end, _) in enumerate(experimental_data):
                if bin_start <= dist <= bin_end:
                    n_total[idx] += 1
                    if con > 0:
                        n_connected[idx] += 1

        p_hyp = np.zeros((len(experimental_data), ))

        for idx, (n_con, n_tot, p_exp) in enumerate(zip(n_connected, n_total, [x[2] for x in experimental_data])):
            # test = binomtest(n_con, n_tot, p_exp)
            # p_hyp[idx] = test.pvalue  # gave 0 when p values are too far apart, not informative
            p_hyp[idx] = abs(n_con/n_tot - p_exp)

        return np.mean(p_hyp)

    def helper_func(self, x):

        pruning_parameters = dict()
        pruning_parameters |= self.optimisation_info["extra_pruning_parameters"]

        pruning_parameters["f1"] = x[0]
        pruning_parameters["softMax"] = x[1]
        pruning_parameters["mu2"] = x[2]
        pruning_parameters["a3"] = x[3]

        output_file = os.path.join(self.network_path, f"network-synapses-{uuid.uuid4()}.hdf5")
        print(f"Output file {output_file}")

        self.prune_synapses(pre_type=self.optimisation_info["pre_type"],
                            post_type=self.optimisation_info["post_type"],
                            con_type=self.optimisation_info["con_type"],
                            pruning_parameters=pruning_parameters,
                            output_file=output_file)

        fitness = self.evaluate_fitness(pre_type=self.optimisation_info["pre_type"],
                                        post_type=self.optimisation_info["post_type"],
                                        output_file=output_file,
                                        experimental_data=self.optimisation_info["exp_data"])

        self.optimisation_info["ctr"] += 1

        print(f"({self.optimisation_info['ctr'] }) Evaluating f1 = {x[0]}, SM = {x[1]}, mu2 = {x[2]}, a3 = {x[3]}, fitness: {fitness}")

        return fitness

    def optimize(self, pre_type, post_type, con_type, experimental_data, extra_pruning_parameters):

        self.optimisation_info["pre_type"] = pre_type
        self.optimisation_info["post_type"] = post_type
        self.optimisation_info["con_type"] = con_type
        self.optimisation_info["exp_data"] = experimental_data
        self.optimisation_info["extra_pruning_parameters"] = extra_pruning_parameters
        self.optimisation_info["ctr"] = 0

        bounds = [(0, 1), (0, 20), (0, 5), (0, 1)]

        res = differential_evolution(func=self.helper_func, bounds=bounds)

        return res