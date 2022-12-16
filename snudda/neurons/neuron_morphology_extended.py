import numexpr
import numpy as np

from snudda.neurons.morphology_data import MorphologyData, SectionMetaData
from snudda.utils.snudda_path import snudda_parse_path


class NeuronMorphologyExtended:

    def __init__(self,
                 name=None,
                 position=None,
                 rotation=None,
                 swc_filename=None,
                 snudda_data=None,
                 param_data=None,
                 mech_filename=None,
                 neuron_path=None,
                 parameter_key=None,
                 morphology_key=None,
                 modulation_key=None,
                 load_morphology=True,
                 virtual_neuron=False,
                 axon_stump_id_flag=True,
                 colour=None,
                 logfile=None,
                 verbose=False):

        self.log_file = logfile
        self.verbose = verbose

        self.name = name

        if isinstance(rotation, (list, tuple)):
            rotation = np.array(rotation)
        self._rotation = rotation

        if isinstance(position, (list, tuple)):
            position = np.array(position)
        self._position = position

        self.swc_filename = swc_filename
        self.snudda_data = snudda_data

        self.param_data = param_data
        self.mech_filename = mech_filename
        self.neuron_path = neuron_path

        self.parameter_key = parameter_key
        self.morphology_key = morphology_key
        self.modulation_key = modulation_key

        self.load_morphology = load_morphology
        self.virtual_neuron = virtual_neuron
        self.axon_stump_id_flag = axon_stump_id_flag
        self.colour = colour

        self.morphology_data = dict()

        # Can we remove these:
        self.axon_density_type = None
        self.dend_density = None
        self.axon_density = None
        self.axon_density_bounds_xyz = None
        self.voxel_size = 5e6
        self.density_bin_size = 10e-6
        # Or are they required for axon densities?

        if self.load_morphology:
            self.add_morphology(swc_file=swc_filename, position=position, rotation=rotation)

    @property
    def position(self):
        if "neuron" in self.morphology_data and not (self.morphology_data["neuron"].position == self._position).all():
            raise ValueError(f"Internal inconsistency, position {self._position} differs "
                             f"from morphology_data position {self.morphology_data['neuron'].position}")

        return self._position

    @property
    def rotation(self):
        if "neuron" in self.morphology_data and not (self.morphology_data["neuron"].rotation == self._rotation).all():
            raise ValueError(f"Internal inconsistency, rotation {self._rotation} differs "
                             f"from morphology_data rotation {self.morphology_data['neuron'].rotation}")

        return self._rotation

    @position.setter
    def position(self, position):
        self._position = position

    def rotation(self, rotation):
        self._rotation = rotation

    def add_morphology(self, swc_file, name="neuron", position=None, rotation=None, parent_tree_info=None):

        """
            MorphologyData

            This can hold an entire neuron, or a part of a neuron.

            Args:
                swc_file (str): Path to SWC file
                name (str): Label of subtree, default "neuron" = main neuron tree
                position (np.ndarray): x,y,z coordinates
                rotation (np.ndarray): 3x3 rotation matrix
                parent_tree_info (tuple, optional): Specify subtree attachment point
                                                    (MorphologyData, parent_label, parent_point_idx, arc_factor)

        """

        self.morphology_data[name] = MorphologyData(swc_file=swc_file, parent_tree_info=parent_tree_info)
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
        dend_idx = np.where(section_data[:, 2] == 3)[0]
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

        return synapse_density, dend_idx

    def dendrite_input_locations(self, synapse_density_str, rng, num_locations,
                                 cluster_size=None, cluster_spread=20e-6):

        synapse_density, dend_idx = self.get_weighted_synapse_density(synapse_density_str=synapse_density_str)

        # Iterate over all dendrites.
        geometry = self.morphology_data["neuron"].geometry
        section_data = self.morphology_data["neuron"].section_data
        soma_dist = geometry[:, 4]
        parent_idx = section_data[:, 3]

        comp_len = soma_dist
        comp_len[1:] -= soma_dist[parent_idx[1:]]

        import pdb
        pdb.set_trace()

        assert (comp_len[1:] > 0).all(), "Internal error. Zero or negative compartment lengths."

        comp_synapse_density = (synapse_density - synapse_density[parent_idx]) / 2
        expected_synapses = np.multiply(comp_len, comp_synapse_density)
        expected_sum = np.sum(expected_synapses[dend_idx])

        if (expected_synapses[dend_idx] < 0).any():
            raise ValueError(f"Found negative synapse densities using {synapse_density_str}")

        if expected_sum <= 0:
            raise ValueError(f"All compartments have zero synapse density: {synapse_density_str}")

        import pdb
        pdb.set_trace()

        syn_idx = dend_idx[rng.choice(a=dend_idx, size=num_locations, replace=True,
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


if __name__ == "__main__":

    swc_file = "/home/hjorth/HBP/Snudda/snudda/data/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20190508/WT-0728MSN01-cor-rep-ax.swc"
    nme = NeuronMorphologyExtended(swc_filename=swc_file)
    nme.place(position=[0,0,0], rotation=np.array([[1,0,0],[0,1,0],[0,0,1]]))
    nme.dendrite_input_locations("1*d", rng=np.random.default_rng(1), num_locations=20)

    import pdb
    pdb.set_trace()