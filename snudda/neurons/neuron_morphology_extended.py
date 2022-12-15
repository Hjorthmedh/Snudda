import numexpr
import numpy as np

from snudda.neurons.morphology_data import MorphologyData, SectionMetaData
from snudda.utils.snudda_path import snudda_parse_path


class NeuronMorphologyExtended:

    def __init__(self, name, position, rotation, swc_filename, snudda_data,
                 parameter_key, morphology_key, modulation_key, logfile=None):

        self.log_file = logfile

        self.name = name
        self.position = position
        self.rotation = rotation
        self.swc_filename = swc_filename
        self.snudda_data = snudda_data

        self.parameter_key = parameter_key
        self.morphology_key = morphology_key
        self.modulation_key = modulation_key

        self.morphology_data = dict()

    def add_morphology(self, swc_file, name="neuron", position=None, rotation=None, parent=None):

        self.morphology_data[name] = MorphologyData(swc_file=swc_file, parent=parent)
        self.morphology_data[name].place(position=position, rotation=rotation)

    def section_iterator(self, section_type=None):
        for subtree in self.morphology_data.values():
            for section in subtree.section_iterator(section_type=section_type):
                yield section

    def place(self, rotation=None, position=None, name="neuron"):
        self.morphology_data[name].place(position=position, rotation=rotation)
        return self

    def get_section_coordinates(self, section_id, section_x):
        raise NotImplementedError("Apologies.")

    def write_log(self, text, flush=True, is_error=False, force_print=False):  # Change flush to False in future, debug

        """
        Writes to log file. Use setup_log first. Text is only written to screen if self.verbose=True,
        or is_error = True, or force_print = True.

        test (str) : Text to write
        flush (bool) : Should all writes be flushed to disk directly?
        is_error (bool) : Is this an error, always written.
        force_print (bool) : Force printing, even if self.verbose=False.
        """

        if self.logfile is not None:
            self.logfile.write(f"{text}\n")
            if flush:
                self.logfile.flush()

        if self.verbose or is_error or force_print:
            print(text, flush=True)

    def plot_neuron(self,
                    axis=None,
                    plot_axon=True,
                    plot_dendrite=True,
                    line_style='-',
                    alpha=1.0,
                    plot_origo=None,  # Only use this when plotting hyper voxels
                    plot_scale=1.0,
                    axon_colour=None,
                    dend_colour=None,
                    soma_colour=None,
                    show_plot=True):

        raise NotImplementedError("This function will move to separate plot class.")

    def get_weighted_synapse_density(self, synapse_density_str):

        """ Given synapse_density it returns expected number of synapses in all dendrite compartments """
        section_data = self.morphology_data["neuron"].section_data
        geometry_data = self.morphology_data["neuron"].geometry
        dend_idx = section_data[:, 2] == 3
        d = geometry_data[dend_idx, 4]
        synapse_density = np.full((geometry_data.shape[0],), np.nan)

        try:
            synapse_density[dend_idx] = numexpr.evaluate(synapse_density_str)
        except:
            self.write_log(f"Bad synapse density string: {synapse_density_str}")
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str)
            raise ValueError(f"Bad synapse density {synapse_density_str}")

        if synapse_density.ndim == 0:
            # If synapse_density is a constant, create a np.array of correct size
            synapse_density = np.full((geometry_data.shape[0],), synapse_density)

        return synapse_density, dend_idx

    def dendrite_input_locations(self, synapse_density_str, rng, num_locations,
                                 cluster_size=None, cluster_spread=20e-6):

        synapse_density, dend_idx = self.get_weighted_synapse_density(synapse_density_str=synapse_density_str)

        # Iterate over all dendrites.
        geometry = self.morphology_data["neuron"].geometry
        section_data = self.morphology_data["neuron"].section_data
        soma_dist = geometry[:, 4]
        parent_idx = section_data[:, 4]
        comp_len = soma_dist - soma_dist[parent_idx]
        comp_synapse_density = (synapse_density - synapse_density[parent_idx]) / 2
        expected_synapses = np.multiply(comp_len, comp_synapse_density)
        expected_sum = np.sum(expected_synapses[dend_idx])

        if (expected_synapses[dend_idx] < 0).any():
            raise ValueError(f"Found negative synapse densities using {synapse_density_str}")

        if expected_sum <= 0:
            raise ValueError(f"All compartments have zero synapse density: {synapse_density_str}")

        syn_idx = dend_idx[rng.choice(dend_idx, size=num_locations, replace=True,
                                      p=expected_synapses[dend_idx] / expected_sum)]

        if cluster_size is None or cluster_size == 1:
            comp_x = rng.random(num_locations)
            xyz = comp_x * geometry[syn_idx, :3] + (1-comp_x) * geometry[parent_idx[syn_idx], :3]
            sec_id = section_data[syn_idx, 0]
            sec_x = comp_x * section_data[syn_idx, 1] + (1-comp_x) * section_data[parent_idx[syn_idx], 1]
            dist_to_soma = comp_x * geometry[syn_idx, 4] + (1-comp_x) * geometry[parent_idx[syn_idx], 4]

        else:
            # If either end point of the section is within cluster_spread/2 then we need to
            # include that section as well (and in case of parent, potentially its children).

            kd_tree = self.morphology_data["neuron"].get_kd_tree(compartment_type=3)
            lust_of_closest_point_idx = kd_tree.query_ball_point(x=geometry[syn_idx, :3], r=cluster_spread)

            list_cluster_syn_idx = []

            for closest_point_idx in lust_of_closest_point_idx:
                list_cluster_syn_idx.append(rng.choice(closest_point_idx, size=cluster_size, replace=True))

            cluster_syn_idx = np.concatenate(list_cluster_syn_idx)

            num_locations = len(cluster_syn_idx)
            comp_x = rng.random(num_locations)
            xyz = comp_x * geometry[cluster_syn_idx, :3] + (1 - comp_x) * geometry[parent_idx[cluster_syn_idx], :3]
            sec_id = section_data[cluster_syn_idx, 0]
            sec_x = comp_x * section_data[cluster_syn_idx, 1] + (1 - comp_x) * section_data[parent_idx[cluster_syn_idx], 1]
            dist_to_soma = comp_x * geometry[cluster_syn_idx, 4] + (1 - comp_x) * geometry[parent_idx[cluster_syn_idx], 4]

            # TODO: Check that this is correct!!

        return xyz, sec_id, sec_x, dist_to_soma
