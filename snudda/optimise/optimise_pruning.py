import os
import json
import collections
import numpy as np
import h5py

import timeit

from snudda.utils import SnuddaLoad
from snudda.detect import SnuddaPrune

from scipy.stats import binomtest
from scipy.spatial import distance_matrix
from scipy.optimize import differential_evolution
from snudda.utils.numpy_encoder import NumpyEncoder

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

        self.log_file = None

    def merge_putative_synapses(self, force_merge=False):

        if force_merge:
            merge_info = None
        else:
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

        orig_connectivity_distributions = json.loads(self.prune.hist_file["meta/connectivity_distributions"][()],
                                                     object_pairs_hook=collections.OrderedDict)

        # Next overwrite the pruning parameters
        pre_type_id = self.prune.type_id_lookup[pre_type]
        post_type_id = self.prune.type_id_lookup[post_type]
        orig_key = f"{pre_type}$${post_type}"
        synapse_type_id = orig_connectivity_distributions[orig_key][con_type]["channel_model_id"]

        pruning = self.prune.complete_pruning_info(pruning_parameters)
        pruning_other = None

        self.prune.connectivity_distributions = dict([])
        self.prune.connectivity_distributions[pre_type_id, post_type_id, synapse_type_id] = (pruning, pruning_other)

        # Update the config file in self.prune so the data gets written to network_synapses

        assert len(self.merge_files_syn) == 1, f"merge_Files_syn should be a list with one file only"

        # print(f"Writing to {output_file} (*)")

        # Clear the out file
        self.prune.out_file = None

        out_file = h5py.File(output_file, "w")

        if self.merge_files_syn[0] is not None:
            with h5py.File(self.merge_files_syn[0], "r") as f_syn:
                # print(f"Writing file {output_file}")
                # We need to make sure that we keep the pruned data files separate

                num_syn, num_syn_kept = self.prune.prune_synapses(synapse_file=f_syn,
                                                                  output_filename=out_file,
                                                                  row_range=None,
                                                                  close_out_file=False,
                                                                  close_input_file=True,
                                                                  merge_data_type="synapses")
        else:
            num_syn, num_syn_kept = self.prune.prune_synapses(synapse_file=None,
                                                              output_filename=out_file,
                                                              row_range=None,
                                                              close_out_file=False,
                                                              close_input_file=True,
                                                              merge_data_type="synapses")

        if self.merge_files_gj[0] is not None:
            with h5py.File(self.merge_files_gj[0], "r") as f_gj:

                num_gj, num_gj_kept = self.prune.prune_synapses(synapse_file=f_gj,
                                                                output_filename=out_file,
                                                                row_range=None,
                                                                close_out_file=False,
                                                                close_input_file=True,
                                                                merge_data_type="gap_junctions")
        else:
            num_gj, num_gj_kept = self.prune.prune_synapses(synapse_file=None,
                                                            output_filename=out_file,
                                                            row_range=None,
                                                            close_out_file=False,
                                                            close_input_file=True,
                                                            merge_data_type="gap_junctions")

        if not os.path.exists(output_file):
            print(f"Output file missing {output_file}")
            import pdb
            pdb.set_trace()

        try:
            assert output_file == self.prune.out_file.filename, f"Expected name {output_file}, had name {self.prune.out_file.name}"
        except:
            import traceback
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()

        n_synapses = out_file["network/num_synapses"][()]
        n_gj = out_file["network/num_gap_junctions"][()]

        out_file["network/synapses"].resize((n_synapses, out_file["network/synapses"].shape[1]))
        out_file["network/gap_junctions"].resize((n_gj, out_file["network/gap_junctions"].shape[1]))

        out_file.close()

    def evaluate_fitness(self, pre_type, post_type, output_file, experimental_data,
                         avg_num_synapses_per_pair=None, con_type="synapses"):

        """

            Args:
                pre_type
                post_type
                output_file: path to output file from prune
                experiment_data: [(bin start, bin end, P)]
                avg_num_synapses_per_pair: avg_num

        """

        # print(f"Opening {output_file}")
        snudda_load = SnuddaLoad(network_file=output_file)
        snudda_data = snudda_load.data

        connection_matrix = np.zeros((snudda_data["num_neurons"], snudda_data["num_neurons"]))

        pre_id = snudda_load.get_neuron_id_of_type(neuron_type=pre_type)
        post_id = snudda_load.get_neuron_id_of_type(neuron_type=post_type)

        pre_mask = np.zeros((snudda_data["num_neurons"],), dtype=bool)
        post_mask = np.zeros((snudda_data["num_neurons"],), dtype=bool)

        pre_mask[pre_id] = True
        post_mask[post_id] = True

        # print(f"snudda_data.keys = {snudda_data.keys()}")

        for row in snudda_data[con_type]:
            if pre_mask[row[0]] and post_mask[row[1]]:
                # Only include connections between the right pre and post types
                connection_matrix[row[0], row[1]] += 1

        pos = snudda_data["neuron_positions"]
        # We need to extract the parts relevant
        dist_matrix = distance_matrix(pos[pre_mask, :], pos[post_mask, :])

        n_connected = np.zeros((len(experimental_data),), dtype=int)
        n_total = np.zeros((len(experimental_data),), dtype=int)
        n_syn = 0
        n_pairs = 0
        n_syn_list = []

        # We also need to extract the relevant corresponding parts in connection_matrix
        con_mat = connection_matrix[pre_mask, :][:, post_mask]

        assert con_mat.shape == dist_matrix.shape, "Mismatch in shape"

        for dist, con in zip(dist_matrix.flatten(), con_mat.flatten()):

            if con > 0:
                n_syn += con
                n_pairs += 1
                n_syn_list.append(con)

            for idx, (bin_start, bin_end, _) in enumerate(experimental_data):
                if bin_start <= dist <= bin_end:
                    n_total[idx] += 1
                    if con > 0:
                        n_connected[idx] += 1

        n_syn_list = np.array(n_syn_list)

        # print(f"P={n_connected/n_total} .. {n_connected} .. {n_total}")

        p_hyp = np.zeros((len(experimental_data), ))

        for idx, (n_con, n_tot, p_exp) in enumerate(zip(n_connected, n_total, [x[2] for x in experimental_data])):
            # test = binomtest(n_con, n_tot, p_exp)
            # p_hyp[idx] = test.pvalue  # gave 0 when p values are too far apart, not informative
            p_hyp[idx] = abs(n_con/n_tot - p_exp) * 100

        if avg_num_synapses_per_pair is not None:
            if n_pairs > 0:
                #per_pair_error = np.sum(abs(avg_num_synapses_per_pair - n_syn_list))/n_pairs
                per_pair_error = 5*np.sum(np.square(avg_num_synapses_per_pair - n_syn_list))/n_pairs

            else:
                per_pair_error = 5*abs(avg_num_synapses_per_pair)

            error = np.sum(p_hyp) + per_pair_error
        else:
            error = np.sum(p_hyp)

        return error

    @staticmethod
    def opt_helper_func(x, *args):

        pruning_parameters = dict()

        # args must be (optimisation_info)
        if args is None:
            raise ValueError("args must contain optimisation_info")
        else:
            optimisation_info = args[0]
            param_names = optimisation_info["param_names"]

            if "extra_pruning_parameters" in optimisation_info:
                pruning_parameters |= optimisation_info["extra_pruning_parameters"]

        for p_name, p_value in zip(param_names, x):
            pruning_parameters[p_name] = p_value

        if "output_file" in optimisation_info:
            output_file = optimisation_info["output_file"]
        else:
            output_file = os.path.join(optimisation_info["network_path"], "temp", f"network-synapses-{uuid.uuid4()}.hdf5")
        # print(f"Output file {output_file}")

        # This trick allows us to reuse the same OptimisePruning object, will be faster
        op = OptimisePruning.get_op(optimisation_info)

        op.prune_synapses(pre_type=optimisation_info["pre_type"],
                          post_type=optimisation_info["post_type"],
                          con_type=optimisation_info["con_type"],
                          pruning_parameters=pruning_parameters,
                          output_file=output_file)

        if optimisation_info["con_type"].lower() == "gap_junction":
            con_type = "gap_junctions"
        else:
            con_type = "synapses"

        fitness = op.evaluate_fitness(pre_type=optimisation_info["pre_type"],
                                      post_type=optimisation_info["post_type"],
                                      output_file=output_file,
                                      experimental_data=optimisation_info["exp_data"],
                                      avg_num_synapses_per_pair=optimisation_info["avg_num_synapses_per_pair"],
                                      con_type=con_type)

        # print(f"Evaluating f1 = {x[0]}, fitness: {fitness}\n{output_file}\n")
        # print(f"Fitness: {fitness}")

        OptimisePruning.report_fitness(fitness)

        if "output_file" not in optimisation_info:
            os.remove(output_file)

        return fitness

    @staticmethod
    def report_fitness(fitness):

        OptimisePruning.ctr += 1
        OptimisePruning.fitness = min(OptimisePruning.fitness, fitness)

        if OptimisePruning.ctr % 100 == 0:
            print(f"Worker iter: {OptimisePruning.ctr}, fitness {OptimisePruning.fitness}")

    @staticmethod
    def get_op(optimisation_info):

        # I am sorry, this is ugly... but need to get this pickeable

        if "op" in vars(OptimisePruning) and OptimisePruning.op is not None:
            op = OptimisePruning.op
        else:
            op = OptimisePruning(network_path=optimisation_info["network_path"])
            merge_info = op.prune.get_merge_info()
            op.merge_files_syn, op.merge_neuron_range_syn, op.merge_syn_ctr, \
                op.merge_files_gj, op.merge_neuron_range_gj, op.merge_gj_ctr = merge_info

            OptimisePruning.op = op
            OptimisePruning.ctr = 0
            OptimisePruning.fitness = np.inf

        return op

    def optimize(self, pre_type, post_type, con_type,
                 experimental_data,
                 param_names, param_bounds,
                 extra_pruning_parameters, avg_num_synapses_per_pair=None,
                 workers=1, maxiter=50, tol=0.001, pop_size=None):

        start = timeit.default_timer()

        self.log_file = open(os.path.join(self.network_path, f"{pre_type}-{post_type}-{con_type}-optimisation-log.txt"), "wt")

        if pop_size is None:
            pop_size = self.pop_size

        if param_bounds == "default":
            param_bounds = []
            for p_name in param_names:
                if p_name == "f1":
                    param_bounds.append((0, 1))
                elif p_name == "softMax":
                    param_bounds.append((0, 20))
                elif p_name == "mu2":
                    param_bounds.append((0, 5))
                elif p_name == "a3":
                    param_bounds.append((0, 1))
                else:
                    raise ValueError(f"No default parameter bounds for {p_name} (f1, softMax, mu2, a3)")

        self.optimisation_info["pre_type"] = pre_type
        self.optimisation_info["post_type"] = post_type
        self.optimisation_info["con_type"] = con_type
        self.optimisation_info["exp_data"] = experimental_data
        self.optimisation_info["avg_num_synapses_per_pair"] = avg_num_synapses_per_pair
        self.optimisation_info["extra_pruning_parameters"] = extra_pruning_parameters
        self.optimisation_info["ctr"] = 0
        self.optimisation_info["network_path"] = self.network_path
        self.optimisation_info["param_names"] = param_names
        self.optimisation_info["param_bounds"] = param_bounds

        optimisation_info = self.optimisation_info

        res = differential_evolution(func=OptimisePruning.opt_helper_func, args=(optimisation_info,),
                                     bounds=optimisation_info["param_bounds"], workers=workers,
                                     maxiter=maxiter, tol=tol, popsize=pop_size)

        # Rerun the best parameters, and keep data as network-synapses.hdf5
        optimisation_info["output_file"] = os.path.join(self.network_path, "network-synapses.hdf5")
        OptimisePruning.opt_helper_func(res.x, optimisation_info)

        duration = timeit.default_timer() - start
        self.log_file.write(f"Duration: {duration} s\n")
        print(f"Duration: {duration} s")

        self.log_file.close()

        # Clear the cache of OptimisePruning, so people can do other networks in same program
        if "op" in vars(OptimisePruning):
            del OptimisePruning.op

        return res

    @staticmethod
    def get_parameters(res, optimisation_info, truncate_decimals=False):

        param_names = optimisation_info["param_names"]

        pruning_parameters = dict()
        if "extra_pruning_parameters" in optimisation_info:
            pruning_parameters |= optimisation_info["extra_pruning_parameters"]

        for p_name, p_value in zip(param_names, res.x):
            if truncate_decimals:
                p_value = np.around(p_value, decimals=truncate_decimals)
            pruning_parameters[p_name] = p_value

        return pruning_parameters

    def export_json(self, file_name, res, append=False, truncate_decimals=4):

        pre_type = self.optimisation_info["pre_type"]
        post_type = self.optimisation_info["post_type"]
        connection_type = self.optimisation_info["con_type"]

        param_names = self.optimisation_info["param_names"]

        if append and os.path.isfile(file_name):
            with open(file_name, "r") as f:
                config_data = json.load(f)
        else:
            config_data = {"connectivity": {}}

        n_params = len(res.x)

        if f"{pre_type},{post_type}" in config_data["connectivity"]:
            con_data = config_data["connectivity"][f"{pre_type},{post_type}"]
        else:
            con_data = dict()

        if connection_type not in con_data:
            con_data[connection_type] = dict()

        con_data[connection_type]["pruning"] = self.get_parameters(res=res, optimisation_info=self.optimisation_info,
                                                                   truncate_decimals=truncate_decimals)
        if "pruning_other" in con_data[connection_type]:
            del con_data[connection_type]["pruning_other"]

        config_data["connectivity"][f"{pre_type},{post_type}"] = con_data

        with open(file_name, "w") as f:
            json.dump(config_data, f, cls=NumpyEncoder, indent=4)

