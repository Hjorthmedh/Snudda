import os
import pickle

import numpy as np
from scipy.spatial import cKDTree

import snudda.utils

# TODO: Move constants like 1000 * sec_x to separate file


class SectionMetaData:

    """ Holds parent_id, children_id, points_id"""

    __slots__ = ["section_id", "parent_section_id", "parent_point_idx", "parent_section_type",
                 "child_section_id", "point_idx", "section_type",
                 "morphology_data", "neuron_id"]

    section_id: int
    parent_section_id: int
    parent_point_idx: int
    parent_section_type: int
    child_section_id: np.ndarray  # First row is child_section_id, second row is child_section_type
    point_idx: np.ndarray
    section_type: int
    morphology_data: object
    neuron_id: int  # Is this set, is it used?

    def __init__(self, section_id, section_type, morphology_data, build_section=True):

        self.morphology_data = morphology_data
        self.section_id = section_id
        self.section_type = section_type

        self.point_idx = None
        self.parent_section_id = None
        self.parent_point_idx = None
        self.parent_section_type = None
        self.child_section_id = None

        if build_section:
            self.build_section()

    def build_section(self):

        idx = np.where((self.morphology_data.section_data[:, 0] == self.section_id)
                       & (self.morphology_data.section_data[:, 2] == self.section_type))[0].astype(np.int32)

        if len(idx) == 0:
            raise ValueError(f"Section id {self.section_id} has no points in morphology_data")

        if not (np.diff(idx) == 1).all():
            raise ValueError(f"Points on section must be consecutive")

        parent_idx = self.morphology_data.section_data[idx[0], 3]
        if parent_idx == -1:
            # Special case, section is soma
            self.point_idx = idx
            self.parent_point_idx = -1
            self.parent_section_type = -1
            self.parent_section_id = -9999
        elif self.morphology_data.section_data[parent_idx, 2] != self.morphology_data.section_data[idx[0], 2]:
            # Special case, root node -- parent section is of different type (e.g. soma -- dend)
            self.point_idx = idx
            self.parent_point_idx = int(self.morphology_data.section_data[self.point_idx[0], 3])
            self.parent_section_type = int(self.morphology_data.section_data[parent_idx, 2])
            self.parent_section_id = int(self.morphology_data.section_data[self.parent_point_idx, 0])
        else:
            self.point_idx = np.concatenate(([self.morphology_data.section_data[idx[0], 3]], idx), dtype=np.int32)
            self.parent_point_idx = int(self.morphology_data.section_data[self.point_idx[0], 3])
            self.parent_section_type = int(self.morphology_data.section_data[parent_idx, 2])
            self.parent_section_id = int(self.morphology_data.section_data[self.parent_point_idx, 0])

        # By definition only the last point in a section can be a parent to other sections
        child_idx = np.where(self.morphology_data.section_data[:, 3] == idx[-1])[0]

        temp_child_section_id = dict()
        for child_section_id, child_type in zip(self.morphology_data.section_data[child_idx, 0],
                                                self.morphology_data.section_data[child_idx, 2]):
            # self.child_section_id is a dict, which holds the child_section_id for different
            # types (e.g. 1=soma, 2=axon, 3=dend, 4=apical)
            if child_type not in temp_child_section_id:
                temp_child_section_id[child_type] = []

            temp_child_section_id[child_type].append(child_section_id)

#        for child_type in self.child_section_id.keys():
#            self.child_section_id[child_type] = np.array(self.child_section_id[child_type])

        section_id = []
        section_type = []
        for child_type in temp_child_section_id.keys():
            section_id += temp_child_section_id[child_type]
            section_type += [child_type]*len(temp_child_section_id[child_type])

        # Convert dictionary to np.array to conserve memory
        self.child_section_id = np.array([section_id, section_type], dtype=np.int16)

        # Also check it above
        bastard_idx = np.where((idx[0] <= self.morphology_data.section_data[:, 3])
                               & (self.morphology_data.section_data[:, 3] < idx[-1]))[0]

        if not (bastard_idx == idx[1:]).all():
            raise ValueError(f"Only last point in section may have children outside section.")

    def section_length(self):

        comp_len = self.morphology_data.geometry[self.point_idx[-1], 4] \
                   - self.morphology_data.geometry[self.point_idx[0], 4]

        # Comment: Soma should not be part of the point_idx, so this is not needed.
        # if self.morphology_data.section_data[self.point_idx[0], 2] == 1:
        #     # If first point (parent point) is soma then we have to be careful and subtract soma radie if inside.
        #     comp_len -= self.morphology_data.geometry[self.point_idx[0], 3]

        if comp_len <= 0:
            raise ValueError(f"Negative section length detected: {self.morphology_data.swc_file}")

        return comp_len

    def clone(self, new_morphology_data, share_memory=True):
        new_smd = SectionMetaData(section_id=self.section_id, section_type=self.section_type,
                                  morphology_data=new_morphology_data, build_section=False)

        if share_memory:
            new_smd.point_idx = self.point_idx
            new_smd.parent_point_idx = self.parent_point_idx
            new_smd.child_section_id = self.child_section_id
            new_smd.parent_section_id = self.parent_section_id

            # Prevent the user from changing these now that the memory is shared
            self.point_idx.setflags(write=False)

            # for cs in self.child_section_id.values():
            #     cs.setflags(write=False)

            self.child_section_id.setflags(write=False)

        else:
            new_smd.point_idx = self.point_idx.copy()
            new_smd.parent_point_idx = self.parent_point_idx
            new_smd.parent_section_id = self.parent_section_id

            new_smd = dict()
            for key, val in self.child_section_id.items():
                new_smd.child_section_id[key] = val.copy()

        return new_smd

    @property
    def position(self):
        return self.morphology_data.geometry[self.point_idx, :3]

    @property
    def radie(self):
        return self.morphology_data.geometry[self.point_idx, 3]

    @property
    def soma_distance(self):
        return self.morphology_data.geometry[self.point_idx, 4]

    @property
    def section_x(self):
        # Double check that this creates a copy of the data before overwriting first element with 0
        sec_x = self.morphology_data.section_data[self.point_idx, 1] / 1e3
        if self.parent_section_type != 1:
            # If parent is not soma, set first section_x to 0
            sec_x[0] = 0
        return sec_x

    def soma_distance_at(self, section_x):
        return np.interp(section_x, self.section_x, self.soma_distance)


class MorphologyData:

    """
        MorphologyData

        This can hold an entire neuron, or a part of a neuron.

        Args:
            swc_file (str): Path to SWC file
                parent_tree_info (tuple, optional): Specify subtree attachment point
                                                (MorphologyData, parent_label, parent_point_idx, arc_factor)

    """

    __slots__ = ["swc_file", "snudda_data", "verbose", "cache_version",
                 "geometry", "section_data", "sections",
                 "point_lookup", "rotation", "position", "parent_tree_info", "is_loaded", "kd_tree_lookup"]

    def __init__(self, swc_file=None, parent_tree_info=None, snudda_data=None,
                 verbose=False, use_cache=True, lazy_loading=False):

        self.swc_file = swc_file
        self.snudda_data = snudda_data
        self.verbose = verbose
        self.cache_version = 1.4

        self.geometry = None      # x, y, z, r, soma_dist (float)
        self.section_data = None  # section_id, section_x (*1000), section_type (int), parent_point_id (int)
        self.sections = None      # dictionary section_id --> SectionMetaData

        self.point_lookup = dict()    # "dend" --> np.array of point_id for dend points

        self.rotation = None
        self.position = None

        self.parent_tree_info = parent_tree_info     # parent tree, if subtree
        self.is_loaded = False
        
        if not lazy_loading and self.swc_file is not None:
            self.load_swc_file(swc_file=self.swc_file, use_cache=use_cache)

        self.kd_tree_lookup = dict()

    def section_iterator_selective(self, section_type, section_id):

        """ Iterates over all sections of a specific type.

        Args:
            section_type: 1 = soma, 2 = axon, 3 = dend
            section_id: ID of sections to iterate over"""

        if section_id is None:
            section_id = self.sections[section_type].keys()

        for sid in section_id:
            yield self.sections[section_type][sid]

    def has_axon(self):
        return len(self.sections[2]) > 0

    def has_dendrite(self):
        return len(self.sections[3]) > 0

    def load_swc_file(self, swc_file=None, remapping_types={4: 3}, use_cache=True):

        """ Loads SWC morphology, not SNUDDA_DATA aware (file must exist).

            Args:
                swc_file (str): Path to swc file
                remapping_types (dict): Remapping of compartment types (default: 4 (apical) -> 3 (normal dendrites))
                use_cache (bool): Save and load neuron morphology to cache
        """
        if swc_file is None:
            swc_file = self.swc_file

        swc_file = snudda.utils.snudda_parse_path(swc_file, self.snudda_data)

        if not os.path.isfile(swc_file):
            raise FileNotFoundError(f"Missing SWC file '{swc_file}'")

        if use_cache:
            cache_file, valid_cache = self.get_cache_file()
            if valid_cache and self.load_cache():
                self.is_loaded = True
                return

        data = np.loadtxt(swc_file)

        if len(data.shape) == 1:
            # This is to handle case when we only have a soma
            data = data.reshape([1, 7])

        if any(np.diff(data[:, 0]) != 1):
            raise IndexError(f"SWC file has gaps in ID numbering ({swc_file})")

        if data[0, 0] != 1:
            raise IndexError(f"ID does not start from 1 ({swc_file})")

        if not (data[0, 2:5] == [0, 0, 0]).all():
            raise ValueError(f"Does not have root centered at origo ({swc_file})")

        if data[0, 6] != -1:
            raise ValueError(f"First element must be root, and have parent ID -1 ({swc_file})")

        if not (data[:, 0] > data[:, 6]).all():
            raise ValueError(f"Parent ID must be lower than row ID ({swc_file})")

        item, count = np.unique(data[:, 0], return_counts=True)
        if (count > 1).any():
            raise ValueError(f"Duplicate index: {item[count > 1]} ({swc_file})")

        if self.parent_tree_info is not None and (data[:, 1] != 2).any():
            raise ValueError(f"Only axonal compartments allowed when subtree of neuron")

        self.geometry = np.zeros((data.shape[0], 5), dtype=np.single)  # 2023-04-24: float -> single, to save memory
        self.geometry[:, :4] = data[:, 2:6] * 1e-6  # x, y, z, r -- converted to meter

        # Store metadata for points
        self.section_data = np.full((data.shape[0], 4), -1, dtype=np.int32)
        self.section_data[:, 2] = data[:, 1]
        self.section_data[0, 3] = -1
        parent_row_id = data[1:, 6].astype(int) - 1
        self.section_data[1:, 3] = parent_row_id

        # This remaps apical dendrites to normal dendrites 4 --> 3 (by default)
        for old_key, new_key in remapping_types.items():
            key_idx = self.section_data[:, 2] == old_key
            self.section_data[key_idx] = new_key

        if (np.abs(self.section_data[:, 2] - data[:, 1]) > 1e-12).any():
            raise ValueError(f"Internal error, non integer ID numbers detected ({swc_file})")

        # We previously removed the points inside the soma, but now keep them
        # Got problems with short dendrites.
        # self.delete_points_inside_soma()

        # OBS, parent_row_id is updated when we delete_points_inside_soma
        parent_row_id = self.section_data[1:, 3]

        # Calculate distance to soma and store in self.geometry
        comp_length = np.linalg.norm(self.geometry[parent_row_id, :3] - self.geometry[1:, :3], axis=1)

        for comp_id, parent_id, c_len in zip(range(1, len(parent_row_id)+1), parent_row_id, comp_length):

            self.geometry[comp_id, 4] = self.geometry[parent_id, 4] + c_len

            # if data[0, 1] == 1 and parent_id == 0:
            #     # We need to subtract soma radius from first compartment connecting to soma
            #     self.geometry[comp_id, 4] = max(0, c_len - self.geometry[0, 3])
            # else:
            #     # distance to soma = parents distance to soma + compartment length
            #     self.geometry[comp_id, 4] = self.geometry[parent_id, 4] + c_len

        if (self.geometry[1:, 4] < 0).any():
            raise ValueError("Found compartments with 0 or negative distance to soma.")

        self.build_tree()

        if use_cache:
            self.save_cache(skip_check=True)  # skip_check since we have not done any rotations

        self.is_loaded = True

    def build_tree(self):

        # New sections are triggered when:
        # -- At branch points
        # -- Change of section type
        parent_id, counts = np.unique(self.section_data[:, 3], return_counts=True)
        branch_id = parent_id[counts > 1]
        type_switch_id = np.where(self.section_data[self.section_data[:, 3], 2] - self.section_data[:, 2] != 0)[0]
        type_switch_id = type_switch_id[type_switch_id != 0]
        edge_id = np.union1d(branch_id, self.section_data[type_switch_id, 3])
        edge_flag = np.zeros((self.section_data.shape[0],), dtype=bool)
        edge_flag[edge_id] = True

        section_counter = dict()

        # Assign section id to all points in section_data
        for idx, row in enumerate(self.section_data):
            section_type = row[2]
            parent_id = row[3]

            if parent_id == -1 or edge_flag[parent_id]:

                # https://github.com/neuronsimulator/nrn/blob/5038de0b79ddf7da9b536639989da4c10dbae7f7/share/lib/hoc/import3d/read_swc.hoc?fbclid=IwAR2kEJOcWkbze8i6G2t9uUVZn5MfmxdSHtm3yzWdP240guJY9KFCalUMvug#L304
                if parent_id == 0 and idx in branch_id:

                    # Special case, parent is soma, and the point itself is a branch point
                    # then mark it as section_type = 0, to not create a one point section
                    self.section_data[idx, 2] = 0
                    section_type = 0

                # Parent point is edge, create new section
                if section_type not in section_counter:
                    section_counter[section_type] = 0
                else:
                    section_counter[section_type] += 1

                section_id = section_counter.get(section_type)
                self.section_data[idx, 0] = section_id
            else:
                # Parent was not an edge, inherit section id
                self.section_data[idx, 0] = self.section_data[parent_id, 0]

        # Calculate section_x for all points in section_data
        for section_type in section_counter:
            for section_id in range(section_counter[section_type]+1):
                idx = np.where((self.section_data[:, 0] == section_id) & (self.section_data[:, 2] == section_type))[0]

                if len(idx) == 1 and self.section_data[idx, 2] == 1:
                    self.section_data[idx, 1] = 1000 * 0.5
                    continue

                if not (np.diff(idx) == 1).all():
                    try:
                        raise ValueError(f"Points on a {section_type} (1=soma, 2=axon, 3=dend) section must be consecutive")
                    except:
                        import traceback
                        print(traceback.format_exc())
                        import pdb
                        pdb.set_trace()

                parent_idx = self.section_data[idx, 3]

                # Use soma distance as a shortcut for calculating length of compartments
                dx = self.geometry[idx, 4] - self.geometry[parent_idx, 4]

                # Special case, if parent is soma, then distance to soma should be set to 0
                if self.section_data[parent_idx[0], 2] == 1:
                    dx[0] = 0

                if len(dx) == 1 and dx[0] == 0:
                    self.section_data[idx, 1] = 0
                else:
                    self.section_data[idx, 1] = 1000 * np.cumsum(dx) / np.sum(dx)

        # Build the actual tree
        self.sections = dict()
        for section_type in section_counter:
            if section_type not in self.sections:
                self.sections[section_type] = dict()

            for section_id in range(0, section_counter[section_type] + 1):
                self.sections[section_type][section_id] = SectionMetaData(section_id=section_id,
                                                                          section_type=section_type,
                                                                          morphology_data=self)
        # Create point lookup
        for section_type in self.sections:
            idx = np.where(self.section_data[:, 2] == section_type)[0]
            self.point_lookup[section_type] = idx

    def save_cache(self, skip_check=False):

        cache_file, _ = self.get_cache_file()

        if not skip_check and (self.rotation is not None or self.position is not None):
            raise ValueError(f"Position and rotation must be None when calling save_cache: {self.swc_file}")

        if self.parent_tree_info is not None:
            raise ValueError(f"Parent tree info must be None when calling save_cache: {self.swc_file}")

        data = dict()
        data["cache_version"] = self.cache_version
        data["swc_file"] = self.swc_file
        data["geometry"] = self.geometry
        data["section_data"] = self.section_data

        data["sections"] = dict()
        for sect_type in self.sections.keys():
            data["sections"][sect_type] = dict()
            for sect_key, sect_value in self.sections[sect_type].items():
                data["sections"][sect_type][sect_key] = dict()
                data["sections"][sect_type][sect_key]["section_id"] = sect_value.section_id
                data["sections"][sect_type][sect_key]["section_type"] = sect_value.section_type

                data["sections"][sect_type][sect_key]["parent_point_idx"] = sect_value.parent_point_idx
                data["sections"][sect_type][sect_key]["parent_section_id"] = sect_value.parent_section_id
                data["sections"][sect_type][sect_key]["parent_section_type"] = sect_value.parent_section_type
                data["sections"][sect_type][sect_key]["point_idx"] = sect_value.point_idx

                # data["sections"][sect_type][sect_key]["child_section_id"] = dict()
                # for ch_key, ch_value in sect_value.child_section_id.items():
                #     data["sections"][sect_type][sect_key]["child_section_id"][ch_key] = ch_value

                data["sections"][sect_type][sect_key]["child_section_id"] = sect_value.child_section_id

        try:
            with open(cache_file, "wb") as f:
                pickle.dump(data, f)
        except:
            if self.verbose:
                print(f"Unable to save cache file {cache_file} -- do you have write permission?")

    def load_cache(self):

        cache_file, valid_cache = self.get_cache_file()
        cache_loaded = False

        if valid_cache:
            try:
                with open(cache_file, "rb") as f:
                    data = pickle.load(f)

                if self.cache_version != data["cache_version"]:
                    raise ValueError(f"Cache version mismatch: {data['cache_version']} (required {self.cache_version})")

                if self.swc_file != data["swc_file"]:
                    raise ValueError(f"Cache mismatch. Different paths:\nRequested: {self.swc_file}\nCached: {data['swc_file']}")

                self.geometry = data["geometry"]
                self.section_data = data["section_data"]

                self.sections = dict()
                for sect_type in data["sections"].keys():
                    assert np.issubdtype(type(sect_type), np.integer), \
                        f"sec_type must be int, found {sect_type} ({type(sect_type)})"
                    self.sections[sect_type] = dict()
                    for sect_id, sect_val in data["sections"][sect_type].items():
                        assert np.issubdtype(type(sect_id), np.integer), \
                            f"sec_id key must be int, found sec_id = {sect_id} ({type(sect_id)}"
                        sec = SectionMetaData(section_id=sect_id, section_type=sect_type,
                                              morphology_data=self, build_section=False)
                        sec.point_idx = sect_val["point_idx"]

                        assert sec.point_idx.dtype == np.int32, f"Old format. New format is 32-bit integer"

                        sec.parent_point_idx = int(sect_val["parent_point_idx"])
                        sec.parent_section_id = int(sect_val["parent_section_id"])

                        sec.parent_section_type = int(sect_val["parent_section_type"])

                        # sec.child_section_id = dict()
                        # for ch_key, ch_val in sect_val["child_section_id"].items():
                        #     sec.child_section_id[ch_key] = ch_val

                        sec.child_section_id = sect_val["child_section_id"]

                        self.sections[sect_type][sect_id] = sec

                cache_loaded = True

            except:
                if self.verbose:
                    import traceback
                    print(traceback.format_exc())
                    print(f"Failed to load cache from {cache_file}")

        return cache_loaded

    def get_cache_file(self):

        swc_file = snudda.utils.snudda_parse_path(self.swc_file, self.snudda_data)
        cache_file = f"{swc_file}-cache.pickle"

        if os.path.isfile(cache_file):
            swc_time = os.path.getmtime(swc_file)
            cache_time = os.path.getmtime(cache_file)

            if cache_time > swc_time:
                # Cache file is newer than swc file
                return cache_file, True

        # No valid cache file found
        return cache_file, False

    def delete_points_inside_soma(self):

        # If no soma, do nothing.
        if self.section_data[0, 2] != 1:
            return

        # We assume soma is in centre
        soma_radius = self.geometry[0, 3]
        dist_to_soma = np.linalg.norm(self.geometry[:, :3], axis=1)
        remove_idx = np.where(dist_to_soma < soma_radius)[0]
        remove_idx = remove_idx[remove_idx > 0]  # Do not remove the soma please

        dont_remove_idx = []
        safety_ctr = 0

        marked_parents = set()  # Checked nodes

        for r_idx in remove_idx[::-1]:
            if r_idx in marked_parents:
                continue

            p_chain = [r_idx]

            p_idx = self.section_data[r_idx, 3]
            while p_idx > 0 and safety_ctr < 1000:
                p_chain.append(p_idx)
                marked_parents.add(p_idx)
                if dist_to_soma[p_idx] > soma_radius:
                    dont_remove_idx += p_chain
                    break

                p_idx = self.section_data[p_idx, 3]
                safety_ctr += 1

        if safety_ctr >= 1000:
            raise ValueError("delete_points_inside_soma: Iteration limit hit, bad swc file?")

        if len(dont_remove_idx) > 0:
            remove_idx = np.array(sorted(list(set(remove_idx) - set(dont_remove_idx))))

            if self.verbose:
                print(f"Found dendrite in {self.swc_file} that goes out, and in and out of soma. "
                      f"Points inside soma {[x+1 for x in remove_idx]}, "
                      f"but these have grandparents outside soma, "
                      f"so keep {[x+1 for x in dont_remove_idx]} (SWC numbering)")

        for r_idx in remove_idx:
            update_parent_idx = np.where(self.section_data[:, 3] == r_idx)[0]
            self.section_data[update_parent_idx, 3] = 0  # Children with removed parents, have soma as stepparent

        # Reindex the parents
        for r_idx in remove_idx[::-1]:
            self.section_data[self.section_data[:, 3] >= r_idx, 3] -= 1

        if len(remove_idx) > 0:

            if self.verbose:
                print(f"Removing {len(remove_idx)} dendrite points inside soma from {self.swc_file}: "
                      f"{[x+1 for x in remove_idx]} (SWC numbering)")

            self.section_data = np.delete(self.section_data, remove_idx, axis=0)
            self.geometry = np.delete(self.geometry, remove_idx, axis=0)

    def shifting_dendrites_onto_soma(self):

        raise DeprecationWarning("This was not how NEURON did it.")

        """ Dendrites should start right at the surface of the soma. Here we delete interior points, and add points
            on surface. Branches which do not touch soma get their innermost point shifted to the soma surface. """

        # If the first point is not a soma, skip this step
        if self.section_data[0, 2] != 1:
            return

        soma_radius = self.geometry[0, 3]

        # We assume soma is at origo, there are checks for that in the load_swc
        dist_to_soma = np.linalg.norm(self.geometry[:, :3], axis=1)
        inside_soma_idx = np.where(dist_to_soma < soma_radius)[0]
        parent_inside_idx = self.section_data[inside_soma_idx, 3]
        remove_idx = np.unique(parent_inside_idx[parent_inside_idx > 0])

        for r_idx in remove_idx[::-1]:
            self.section_data[self.section_data[:, 3] >= r_idx, 3] -= 1

        if len(remove_idx) > 0:
            self.section_data = np.delete(self.section_data, remove_idx, axis=0)
            self.geometry = np.delete(self.geometry, remove_idx, axis=0)

        # New points, with soma as parent, should be shifted to surface of soma
        shift_idx = np.where(self.section_data[:, 3] == 0)[0]

        dx = (self.geometry[shift_idx, :3].T / np.linalg.norm(self.geometry[shift_idx, :3], axis=1)[None, :]).T
        self.geometry[shift_idx, :3] = dx * soma_radius

    def clone(self, position, rotation, parent_tree_info=None, share_memory=True):

        """
            Clone the neuron, and place it.

            Args:
                position (np.array): x,y,z coordinate of cloned neuron
                rotation (np.array): 3x3 rotation matrix
                parent_tree_info (optional): Parent tree
                share_memory (bool): If True some numpy arrays are shared with parent neuron

        """

        if self.position is not None or self.rotation is not None:
            raise ValueError("Not allowed to rotate or position a neuron that has already been rotated or positioned")

        new_md = MorphologyData(verbose=self.verbose)
        new_md.swc_file = self.swc_file

        # We can never share memory for geometry, since repositioning a neuron updates geometry
        new_md.geometry = self.geometry.copy()

        if share_memory:
            # Assuming topology of neuron does not change, these values will be constant
            new_md.section_data = self.section_data
            self.section_data.setflags(write=False)

            for p_key, p_value in self.point_lookup.items():
                new_md.point_lookup[p_key] = p_value
                p_value.setflags(write=False)

        else:
            new_md.section_data = self.section_data.copy()

            for p_key, p_value in self.point_lookup.items():
                new_md.point_lookup[p_key] = p_value.copy()  # !!! Should this be copy.deepcopy?

        new_md.sections = dict()

        for sec_type, sec_type_data in self.sections.items():
            new_md.sections[sec_type] = dict()

            for sec_key, sec_value in sec_type_data.items():
                new_md.sections[sec_type][sec_key] = sec_value.clone(new_morphology_data=new_md,
                                                                     share_memory=share_memory)

        new_md.kd_tree_lookup = dict()
        new_md.parent_tree_info = parent_tree_info

        new_md.place(position=position, rotation=rotation)

        return new_md

    def place(self, position=None, rotation=None, parent_tree_info=None, lazy=False):
        if self.position is not None or self.rotation is not None:
            raise ValueError("Not allowed to rotate or position a neuron that has already been rotated or positioned")

        if isinstance(rotation, (list, tuple)):
            rotation = np.array(rotation)
        else:
            rotation = rotation.copy() if rotation is not None else None

        self.rotation = rotation

        if isinstance(position, (list, tuple)):
            position = np.array(position)
        else:
            position = position.copy() if position is not None else None

        self.position = position

        if lazy:
            return

        if not self.is_loaded:
            self.load_swc_file()

        if parent_tree_info is not None:
            self.parent_tree_info = parent_tree_info  # (MorphologyData, point_idx, arc_factor)

        # Here we assume soma is only a point
        if 1 in self.point_lookup:
            soma_position = self.geometry[self.point_lookup[1], :3]
            if not (soma_position == 0).all():
                raise ValueError("Soma must be centered at origo before placement.")
        elif 2 in self.point_lookup:
            # We have no soma, so it is probably an axonal tree, check that it is centered
            axon_root_position = self.geometry[self.point_lookup[2][0], :3]
            if not (axon_root_position == 0).all():
                raise ValueError("Axon root must be centered at origo before placement.")

        if rotation is not None:
            if not np.abs(np.linalg.det(rotation) - 1) < 1e-10 \
                    or not np.abs(np.matmul(rotation, rotation.T) - np.eye(3) < 1e-10).all():
                raise ValueError(f"Not a valid rotation matrix {rotation}")

            self.geometry[:, :3] = np.matmul(self.rotation, self.geometry[:, :3].T).T

        if self.position is not None:
            self.geometry[:, :3] += self.position

        if self.parent_tree_info is not None:
            # We need to update soma distance for subtree based on distance to parent
            # parent_tree_info = (MorphologyData, point_idx, arc_factor) -- attachment point

            parent_object, parent_point_idx, arc_factor = self.parent_tree_info
            parent_position = parent_object.geometry[parent_point_idx, :3]
            parent_soma_distance = parent_object.geometry[parent_point_idx, 4]

            dist_to_parent = np.linalg.norm(self.position - parent_position)
            self.geometry[:, 4] += parent_soma_distance + dist_to_parent * arc_factor

    def section_iterator(self, section_type=None):

        """ Iterates over all sections of a specific type.

        Args:
            section_type: 1 = soma, 2 = axon, 3 = dend """

        if section_type is None:
            for section_dict in self.sections.values():
                for section in section_dict.values():
                    yield section

        elif section_type in self.sections:
            for section in self.sections[section_type].values():
                yield section

        return

    def get_kd_tree(self, compartment_type):

        if compartment_type not in self.kd_tree_lookup:
            comp_idx = np.where(self.section_data[:, 2] == compartment_type)[0]
            coords = self.geometry[comp_idx, :3]

            self.kd_tree_lookup[compartment_type] = cKDTree(coords)

        return self.kd_tree_lookup[compartment_type]

    def get_closest_point(self, coords, compartment_type):

        """ """

        kd_tree = self.get_kd_tree(compartment_type=compartment_type)
        closest_dist, closest_point_idx = kd_tree.query(coords)
        return closest_dist, closest_point_idx


# Separate method to create a random rotation matrix that distributes rotations evenly (this is non-trivial)
def rand_rotation_matrix(deflection=1.0, rand_nums=None):
    """
    Creates a random rotation matrix.

    deflection: the magnitude of the rotation. For 0, no rotation; for 1, completely random
    rotation. Small deflection => small perturbation.
    rand_nums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
    """

    # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c

    if rand_nums is None:
        rand_nums = np.random.uniform(size=(3,))

    theta, phi, z = rand_nums

    theta = theta * 2.0 * deflection * np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0 * np.pi  # For direction of pole deflection.
    z = z * 2.0 * deflection  # For magnitude of pole deflection.

    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.

    r = np.sqrt(z)
    vv = (np.sin(phi) * r, np.cos(phi) * r, np.sqrt(2.0 - z))

    st = np.sin(theta)
    ct = np.cos(theta)

    rr = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))

    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    mm = (np.outer(vv, vv) - np.eye(3)).dot(rr)
    return mm


if __name__ == "__main__":

    # file_name = "/home/hjorth/HBP/BasalGangliaData/data/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026/morphology/WT-0728MSN01-cor-rep-ax-res3.swc"
    file_name = "/home/hjorth/HBP/BasalGangliaData/Parkinson/20220225/PD0/neurons/striatum/dspn/26/WT-1215MSN03-cor-rep-ax-res3-var8.swc"
    md = MorphologyData(swc_file=file_name)

    import pdb
    pdb.set_trace()
