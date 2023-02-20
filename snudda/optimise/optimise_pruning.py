import json
import collections
import numpy as np
from scipy.spatial import distance_matrix
from snudda.utils import SnuddaLoad
from snudda.detect import SnuddaPrune


class OptimisePruning:

    """ Optimises pruning parameters. First it creates a relatively small network and does touch detection on it.
        After that it does multiple pruning attempts (while keeping the files original detection files). """

    def __init__(self, network_path, pop_size=10, epochs=10):

        self.pop_size = pop_size
        self.epochs = epochs

        self.prune = SnuddaPrune(network_path=network_path, keep_files=True)

        self.merge_files_syn = None
        self.merge_neuron_range_syn = None
        self.merge_syn_ctr = None
        self.merge_files_gj = None
        self.merge_neuron_range_gj = None
        self.merge_gj_ctr = None

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

        # We need to make sure that we keep the pruned data files separate
        self.prune.prune_synapses(synapse_file=self.merge_files_syn,
                                  output_filename=output_file,
                                  row_range=None,
                                  close_out_file=False,
                                  close_input_file=True,
                                  merge_data_type="synapses")

    def evaluate_fitness(self, pre_type, post_type, output_file, experiment_data):

        """

            Args:
                pre_type
                post_type
                output_file: path to output file from prune
                experiment_data: [(bin start, bin end, n_con_pairs, n_tot_pairs, P)]

        """

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

        # !!! TODO: Calculate the number of connected and non-connected pairs in the experimental bins
        #           and then calculate how well the model data fits the original experimental data



    def
        # For plotting
        n_bins = 10
        bin_size = 20
        connected_pairs = np.zeros((n_bins, ))
        total_pairs = np.zeros((n_bins, ))

        for connected, d in zip(connection_matrix.flatten(), dist_matrix.flatten()):

            bin = int(d / bin_size)

            if bin < n_bins:
                if connected:
                    connected_pairs[bin] += 1

                total_pairs[bin] += 1




        """ Evaluate the fitness of the connection between pre_type and post_type """


        pass

