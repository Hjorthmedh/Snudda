

class ReplaceMorphology:

    def __init__(self, original_network, new_network, new_snudda_data_path):

        """ This function takes a network-neuron-positions.hdf5 file and replaces the morphologies by switching
            SNUDDA_DATA path to new location, and remapping morphologies to point to new morphologies. """


class MergeSynapseFiles:

    def __init__(self, parent_network_A, parent_network_B, new_network, merge_type=None):

        """ Merging the synapses in parent_network_A and parent_network_B, writing the new synapses to new_network.

            There are different merge styles available:

                "addition" - The neurons in the two networks are part of different networks,
                             just add all neurons together. This will renumber all neuron_id:s.

                "structural_changes" - The neurons in the parent networks must be related in a 1-1 fashion.
                                       This can be used for example to degenerate an existing network by replacing
                                       morphologies with degenerated versions. Important that neuron location is
                                       identical in the two networks.

                                       Synapses on the old morphology are retained if the dendrites and axons remain.
                                       New synapses are added if the new morphology dendrites or axon did not exist
                                       in the old morphologies.





        """

        # TODO: We need to start storing SNUDDA_DATA in the save file (meta),
        #  so the user does not need to specify it all the time

        # Workflow
        # 1. Generate parent_network_A as normal (init, place, detect, prune)
        # 2. Use network-neuron-position.hdf5 from parent_network_A and generate a network-neuron-position.hdf5
        #    for parent_network_B, then do (detect, prune) on that network.
        # 3. Merge synapse files together.
        #

        pass