import os.path
from copy import deepcopy

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
                 colour=None,
                 logfile=None,
                 verbose=False):

        self.logfile = logfile
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
        self.colour = colour

        self.morphology_data = dict()

        # Can we remove these:
        self.axon_density_type = None
        self.dend_density = None
        self.axon_density = None
        self.max_axon_radius = None
        self.axon_density_bounds_xyz = None
        self.voxel_size = 5e6
        self.density_bin_size = 10e-6
        # Or are they required for axon densities?

        if self.load_morphology and swc_filename is not None:
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

    @rotation.setter
    def rotation(self, rotation):
        self._rotation = rotation

    def add_morphology(self, swc_file, name="neuron", position=None, rotation=None, parent_tree_info=None,
                       overwrite=False, morphology_data=None):

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
                morphology_data (optional): MorphologyData object to use, then swc_file is None

        """

        if not overwrite and name in self.morphology_data:
            raise KeyError(f"Error when loading {swc_file}, key {name} already exists in morphology_data")

        if morphology_data is None:
            # No cached morphology data, load it from file (slow)
            self.morphology_data[name] = MorphologyData(swc_file=swc_file, parent_tree_info=parent_tree_info,
                                                        snudda_data=self.snudda_data)
            if position is not None:
                self.morphology_data[name].place(position=position, rotation=rotation)

        else:
            # Use provided morphology data
            self.morphology_data[name] = morphology_data.clone(position=position, rotation=rotation)

    def section_iterator(self, section_type=None, subtree=None):

        if subtree is not None:
            return subtree.section_iterator(section_type=section_type)

        else:
            for subtree in self.morphology_data.values():
                for section in subtree.section_iterator(section_type=section_type):
                    yield section

    def section_iterator_selective(self, section_type, section_id, subtree="neuron"):
        return self.morphology_data[subtree].section_iterator_selective(section_type=section_type,
                                                                        section_id=section_id)

    def place(self, rotation=None, position=None, name="neuron"):

        if name in self.morphology_data:
            self.morphology_data[name].place(position=position, rotation=rotation)

        # If load_morphology is False, we allow "neuron" not to be defined yet.
        if name not in self.morphology_data and (name != "neuron" and not self.load_morphology):
            raise ValueError(f"Neuron morphology data '{name}' not loaded.")

        self.position = position
        self.rotation = rotation

        return self

    def clone(self,
              position=None,
              rotation=None,
              parameter_key=None,
              morphology_key=None,
              modulation_key=None):

        """
        Creates a clone copy of a neuron.

        Args:
            position (float,float,float) : x,y,z coordinate of clone
            rotation (rotation matrix) : Rotation matrix for clone

            parameter_key (str): Parameter Key for clone
            morphology_key (str): Morphology Key for clone
            modulation_key (str): Modulation Key for clone

        """

        # np.set_printoptions(precision=2)
        # print(f"rot {rotation.flatten()}, place pos {position}")

        new_neuron = NeuronMorphologyExtended(name=self.name,
                                              position=None,
                                              rotation=None,
                                              swc_filename=self.swc_filename,
                                              snudda_data=self.snudda_data,
                                              param_data=self.param_data,
                                              mech_filename=self.mech_filename,
                                              neuron_path=self.neuron_path,
                                              parameter_key=self.parameter_key,
                                              morphology_key=self.morphology_key,
                                              modulation_key=self.modulation_key,
                                              load_morphology=False,
                                              virtual_neuron=self.virtual_neuron,
                                              colour=self.colour,
                                              logfile=self.logfile,
                                              verbose=self.verbose)

        new_neuron.position = position.copy() if position is not None else None
        new_neuron.rotation = rotation.copy() if rotation is not None else None

        # Copy over old morphology data
        for md_key, md_value in self.morphology_data.items():
            if md_key == "neuron":
                new_neuron.morphology_data[md_key] = md_value.clone(position=position, rotation=rotation)
            else:
                new_neuron.morphology_data[md_key] = md_value.clone()

        if parameter_key is not None:
            new_neuron.parameter_key = parameter_key

        if modulation_key is not None:
            new_neuron.modulation_key = modulation_key

        if morphology_key != self.morphology_key:
            print("PROBLEM!!")
            import pdb
            pdb.set_trace()
            raise ValueError(f"Not allowed to change morphology_key when cloning: {self.morphology_key} -> {morphology_key}")

        new_neuron.load_morphology = self.load_morphology

        return new_neuron

    def get_section_coordinates(self, section_id, section_x):

        if section_id == -1:
            return self.position

        if "neuron" not in self.morphology_data or 3 not in self.morphology_data["neuron"].sections:
            raise ValueError(f"No dendrites loaded for neuron {self.swc_filename}")

        section = self.morphology_data["neuron"].sections[3][section_id]
        sec_x = section.section_x
        pos = section.position

        closest_idx = np.argmin(np.abs(sec_x - section_x))

        try:
            if sec_x[closest_idx] == section_x:
                coords = pos[closest_idx, :]

            elif sec_x[closest_idx] < section_x and closest_idx < len(sec_x)-1:
                x = (section_x - sec_x[closest_idx]) / (sec_x[closest_idx+1] - sec_x[closest_idx])
                coords = x * pos[closest_idx + 1, :] + (1-x) * pos[closest_idx, :]

            else:
                x = (sec_x[closest_idx] - section_x) / (sec_x[closest_idx] - sec_x[closest_idx - 1])
                coords = x * pos[closest_idx - 1, :] + (1-x) * pos[closest_idx, :]
        except:
            import traceback
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()

        return coords

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

        if plot_origo is None:
            plot_origo = np.array([0, 0, 0])

        self.write_log(f"Plotting neuron {self.swc_filename}")

        if axon_colour is None:
            axon_colour = self.colour if self.colour is not None else 'orange'
        if dend_colour is None:
            dend_colour = self.colour if self.colour is not None else 'black'
        if soma_colour is None:
            soma_colour = self.colour if self.colour is not None else 'black'

        import matplotlib.pyplot as plt

        if axis is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = axis

        section_colour = {2: axon_colour, 3: dend_colour}

        skip_section_types = []
        if not plot_axon:
            skip_section_types.append(2)
        if not plot_dendrite:
            skip_section_types.append(3)

        for morph_key, morphology in self.morphology_data.items():
            for section_type in morphology.sections.keys():

                if section_type in skip_section_types:
                    continue

                for section in morphology.section_iterator(section_type=section_type):
                    if section_type == 1:
                        # Draw soma
                        xyz = section.position[0]
                        radie = section.radie[0]

                        u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
                        x = (radie * np.cos(u) * np.sin(v) + xyz[0] - plot_origo[0]) * plot_scale
                        y = (radie * np.sin(u) * np.sin(v) + xyz[1] - plot_origo[1]) * plot_scale
                        z = (radie * np.cos(v) + xyz[2] - plot_origo[2]) * plot_scale

                        ax.plot_wireframe(x, y, z, color=soma_colour, alpha=alpha)
                    elif section_type == 0:
                        # Do not draw section type 0
                        pass
                    else:
                        xyz = (section.position - plot_origo) * plot_scale
                        ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], linestyle=line_style,
                                alpha=alpha, c=section_colour[section_type])

        if axis is None:
            plt.title(f"Neuron: {os.path.basename(self.swc_filename)}")

            if show_plot:
                plt.ion()
                plt.show()
                plt.draw()
                plt.pause(0.001)

        return ax

        # raise NotImplementedError("This function will move to separate plot class.")

    def get_weighted_synapse_density(self, synapse_density_str):

        """ Given synapse_density it returns expected number of synapses in all dendrite compartments """
        section_data = self.morphology_data["neuron"].section_data
        geometry_data = self.morphology_data["neuron"].geometry

        # We need to evaluate the synapse density at all dendrites and at soma (since we need all dendrite parents)
        keep_mask = section_data[:, 2] == 3
        dend_idx = np.where(keep_mask)[0]
        soma_idx = np.where(section_data[:, 2] <= 1)[0]  # We also evaluate at section_type = 0,
                                                         # ie removed 1-point sections that might be
                                                         # parent points potential secton end points
        keep_mask[soma_idx] = True
        d_idx = np.where(keep_mask)[0]

        d = geometry_data[d_idx, 4]
        synapse_density = np.full((geometry_data.shape[0],), np.nan)

        try:
            synapse_density[d_idx] = numexpr.evaluate(synapse_density_str)
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

        comp_len = soma_dist.copy()
        comp_len[1:] -= soma_dist[parent_idx[1:]]

        assert (comp_len[1:] > 0).all(), "Internal error. Zero or negative compartment lengths."

        comp_synapse_density = (synapse_density + synapse_density[parent_idx]) / 2
        expected_synapses = np.multiply(comp_len, comp_synapse_density)
        expected_sum = np.sum(expected_synapses[dend_idx])

        if (expected_synapses[dend_idx] < 0).any():
            raise ValueError(f"Found negative synapse densities using {synapse_density_str}")

        if expected_sum <= 0:
            raise ValueError(f"All compartments have zero synapse density: {synapse_density_str}")

        if num_locations is not None:
            try:
                syn_idx = rng.choice(a=dend_idx, size=num_locations, replace=True,
                                     p=expected_synapses[dend_idx] / expected_sum)
            except:
                import traceback
                self.write_log(traceback.format_exc(), is_error=True)
                import pdb
                pdb.set_trace()
        else:
            if not (cluster_size is None or cluster_size == 1):
                raise ValueError(f"If cluster_size is set, then num_locations must be set.")

            # Expected synapses is synapses per micrometers, multiply by 1e6 to get per meter (SI)
            syn_idx = dend_idx[np.where(rng.uniform(size=dend_idx.shape) < 1e6*expected_synapses[dend_idx])[0]]
            num_locations = len(syn_idx)

        if cluster_size is None or cluster_size == 1:
            try:
                comp_x = rng.random(num_locations)
                xyz = comp_x[:, None] * geometry[syn_idx, :3] + (1-comp_x[:, None]) * geometry[parent_idx[syn_idx], :3]
                sec_id = section_data[syn_idx, 0]
                sec_x = comp_x * section_data[syn_idx, 1] + (1-comp_x) * section_data[parent_idx[syn_idx], 1]
                dist_to_soma = comp_x * geometry[syn_idx, 4] + (1-comp_x) * geometry[parent_idx[syn_idx], 4]
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

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
            xyz = comp_x[:, None] * geometry[cluster_syn_idx, :3] + (1 - comp_x[:, None]) * geometry[parent_idx[cluster_syn_idx], :3]
            sec_id = section_data[cluster_syn_idx, 0]
            sec_x = comp_x * section_data[cluster_syn_idx, 1] + (1 - comp_x) * section_data[parent_idx[cluster_syn_idx], 1]
            dist_to_soma = comp_x * geometry[cluster_syn_idx, 4] + (1 - comp_x) * geometry[parent_idx[cluster_syn_idx], 4]

            # TODO: Check that this is correct!!

        return xyz, sec_id, sec_x / 1e3, dist_to_soma

    def cluster_synapses(self, sec_id, sec_x, count, distance, rng):

        """
             Randomize sec_x for cluster of 'count' synapses with centre placed at 'sec_id', 'sec_x'
             spread over 'distance' (but constrained to current section extent).

             Args:
                 sec_id : Section id of cluster centre
                 sec_x : Section x of cluster centre
                 count : Number of synapses in cluster
                 distance : Maximal spread of cluster along dendrite
                 rng : Numpy random stream

             Returns:
                 cluster_sec_x : Section x for cluster synapse
                 coords : Coordinates for synapse (in meters)
                 soma_dist : Distance to soma (in meters)
         """

        section = self.morphology_data["neuron"].sections[3][sec_id]  # 3 is dendrite
        geometry = self.morphology_data["neuron"].geometry[section.point_idx, :]
        section_data = self.morphology_data["neuron"].section_data[section.point_idx, :]
        sec_len = section.section_length()

        min_sec_x = np.maximum(1e-3, sec_x - 0.5 * distance / sec_len)
        max_sec_x = np.minimum(1 - 1e-3, sec_x + 0.5 * distance / sec_len)

        cluster_sec_x = rng.uniform(low=min_sec_x, high=max_sec_x, size=count)

        syn_coords = np.zeros((count, 3))
        soma_dist = np.zeros((count,))

        if section.section_type == 1:
            # We are at soma
            syn_coords[:, :] = geometry[0, :3]
            soma_dist[:] = 0
        else:
            # geometry stores section_x * 1e3 (as int), so we need to scale up before comparing
            for i, cx in enumerate(cluster_sec_x*1e3):
                idx = np.sum(cx > section_data[:, 1])
                if idx == 0:
                    syn_coords[i, :] = geometry[idx, :3]
                    soma_dist[i] = geometry[idx, 4]
                else:
                    f = (cx - section_data[idx-1, 1]) / (section_data[idx, 1] - section_data[idx-1, 1])
                    syn_coords[i, :] = (1-f)*geometry[idx-1, :3] + f*geometry[idx, :3]
                    soma_dist[i] = (1-f)*geometry[idx-1, 4] + f*geometry[idx, 4]

        # OBS, syn_coords in meters, and soma dist in meters also
        return cluster_sec_x, syn_coords, soma_dist

    def set_axon_voxel_radial_density(self, density, max_axon_radius):

        """
        Sets axon radial density

        Args:
            density: Axon density f(r), r = radius from soma
            max_axon_radius: Axon density is calculated within a sphere of radius max_axon_radius

        """

        self.write_log("Only saving equation now")

        self.axon_density_type = "r"
        self.axon_density = density
        self.max_axon_radius = max_axon_radius

    def set_axon_voxel_xyz_density(self,
                                   density,
                                   axon_density_bounds_xyz):

        """
        Sets axon density

        Args:
            density: Axon density f(x,y,z), x,y,z = SWC coordinates in relative to soma
            axon_density_bounds_xyz: Bounding box for the axon density in x,y,z
        """

        self.write_log("Only saving equation now")

        self.axon_density_type = "xyz"
        self.axon_density = density
        self.axon_density_bounds_xyz = axon_density_bounds_xyz


if __name__ == "__main__":

    swc_file = "/home/hjorth/HBP/Snudda/snudda/data/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20190508/WT-0728MSN01-cor-rep-ax.swc"
    nme = NeuronMorphologyExtended(swc_filename=swc_file)
    nme.place(position=[0,0,0], rotation=np.array([[1,0,0],[0,1,0],[0,0,1]]))
    nme.dendrite_input_locations("1*d", rng=np.random.default_rng(1), num_locations=20)

    import pdb
    pdb.set_trace()