import os
import numpy as np
from copy import deepcopy
from scipy.spatial import cKDTree

import snudda.utils


# TODO: Move constants like 1000 * sec_x to separate file


class SectionMetaData:

    """ Holds parent_id, children_id, points_id"""

    __slots__ = ["section_id", "parent_section_id",
                 "child_section_id", "point_idx", "section_type", "morphology_data"]

    section_id: int
    parent_section_id: int
    child_section_id: dict
    point_idx: np.ndarray
    section_type: int
    morphology_data: object

    def __init__(self, section_id, section_type, morphology_data):

        self.morphology_data = morphology_data
        self.section_id = section_id
        self.section_type = section_type

        idx = np.where((self.morphology_data.section_data[:, 0] == section_id)
                       & (self.morphology_data.section_data[:, 2] == section_type))[0]

        if len(idx) == 0:
            raise ValueError(f"Section id {section_id} has no points in morphology_data")

        if not (np.diff(idx) == 1).all():
            raise ValueError(f"Points on section must be consecutive")

        if self.morphology_data.section_data[idx[0], 3] < 0:
            # Special case, root node -- do not include soma point in section
            self.point_idx = idx
            self.parent_section_id = -1
        else:
            self.point_idx = np.concatenate(([self.morphology_data.section_data[idx[0], 3]], idx))
            self.parent_section_id = self.morphology_data.section_data[self.point_idx[0], 0]

        # By definition only the last point in a section can be a parent to other sections
        child_idx = np.where(self.morphology_data.section_data[:, 3] == idx[-1])[0]

        self.child_section_id = dict()
        for child_section_id, child_type in zip(self.morphology_data.section_data[child_idx, 0],
                                        self.morphology_data.section_data[child_idx, 2]):
            # self.child_section_id is a dict, which holds the child_section_id for different
            # types (e.g. 1=soma, 2=axon, 3=dend, 4=apical)
            if child_type not in self.child_section_id:
                self.child_section_id[child_type] = []

            self.child_section_id[child_type].append(child_section_id)

        for child_type in self.child_section_id.keys():
            self.child_section_id[child_type] = np.array(self.child_section_id[child_type])

        # Also check it above
        bastard_idx = np.where((idx[0] <= self.morphology_data.section_data[:, 3])
                               & (self.morphology_data.section_data[:, 3] < idx[-1]))[0]

        if not (bastard_idx == idx[1:]).all():
            raise ValueError(f"Only last point in section may have children outside section.")

    def section_length(self):

        comp_len = self.morphology_data.geometry[self.point_idx[-1], 4] \
                   - self.morphology_data.geometry[self.point_idx[0], 4]

        if self.morphology_data.section_data[self.point_idx[0], 2] == 1:
            # If first point (parent point) is soma then we have to be careful and subtract soma radie if inside.
            comp_len -= self.morphology_data.geometry[self.point_idx[0], 3]

        if comp_len <= 0:
            raise ValueError(f"Negative section length detected: {self.morphology_data.swc_file}")

        return comp_len

    def clone(self, new_morphology_data):
        new_smd = SectionMetaData(section_id=self.section_id, section_type=self.section_type,
                                  morphology_data=new_morphology_data)
        return new_smd

    @property
    def position(self):
        return self.morphology_data.geometry[self.point_idx, :3]

    @property
    def radie(self):
        return self.morphology_data.geometry[self.point_idx, 3]


class MorphologyData:

    """
        MorphologyData

        This can hold an entire neuron, or a part of a neuron.

        Args:
            swc_file (str): Path to SWC file
            parent_tree_info (tuple, optional): Specify subtree attachment point
                                                (MorphologyData, parent_label, parent_point_idx, arc_factor)

    """

    def __init__(self, swc_file=None, parent_tree_info=None, snudda_data=None):

        self.swc_file = swc_file
        self.snudda_data = snudda_data

        self.geometry = None      # x, y, z, r, soma_dist (float)
        self.section_data = None  # section_id, section_x (*1000), section_type (int), parent_point_id (int)
        self.sections = None      # dictionary section_id --> SectionMetaData

        self.point_lookup = dict()    # "dend" --> np.array of point_id for dend points

        self.rotation = None
        self.position = None

        self.parent_tree_info = parent_tree_info     # parent tree, if subtree

        if swc_file is not None:
            self.load_swc_file(swc_file=swc_file)

        self.kd_tree_lookup = dict()

    def section_iterator(self, section_type):

        """ Iterates over all sections of a specific type.

        Args:
            section_type: 1 = soma, 2 = axon, 3 = dend """

        for section in self.sections[section_type].values():
            yield section

    def load_swc_file(self, swc_file, remapping_types={4: 3}):

        """ Loads SWC morphology, not SNUDDA_DATA aware (file must exist).

            Args:
                swc_file (str): Path to swc file
                remapping_types (dict): Remapping of compartment types (default: 4 (apical) -> 3 (normal dendrites))
        """

        swc_file = snudda.utils.snudda_parse_path(swc_file, self.snudda_data)

        if not os.path.isfile(swc_file):
            raise FileNotFoundError(f"Missing SWC file '{swc_file}'")

        data = np.loadtxt(swc_file)

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

        self.geometry = np.zeros((data.shape[0], 5), dtype=float)
        self.geometry[:, :4] = data[:, 2:6] * 1e-6  # x, y, z, r -- converted to meter

        # Calculate distance to soma and store in self.geometry
        parent_row_id = data[1:, 6].astype(int) - 1
        comp_length = np.linalg.norm(self.geometry[parent_row_id, :3] - self.geometry[1:, :3], axis=1)

        for comp_id, parent_id, c_len in zip(range(1, len(parent_row_id)+1), parent_row_id, comp_length):
            if data[0, 1] == 1 and parent_id == 0:
                # We need to subtract soma radius from first compartment connecting to soma
                # If the first point is inside the soma, set its compartment length to 0
                if c_len < self.geometry[0, 3]:
                    print(f"Warning: Branch starts inside soma: {swc_file}, line id {comp_id+1}"
                          f" -- will truncate length to 1 micrometer")

                self.geometry[comp_id, 4] = max(1e-6, c_len - self.geometry[0, 3])

            else:
                # distance to soma = parents distance to soma + compartment length
                self.geometry[comp_id, 4] = self.geometry[parent_id, 4] + c_len

        if (self.geometry[1:, 4] <= 0).any():
            import pdb
            pdb.set_trace()
            raise ValueError("Found compartments with 0 or negative length.")

        # Store metadata for points
        self.section_data = np.full((data.shape[0], 4), -1, dtype=int)
        self.section_data[:, 2] = data[:, 1]
        self.section_data[0, 3] = -1
        self.section_data[1:, 3] = parent_row_id

        # This remaps apical dendrites to normal dendrites 4 --> 3 (by default)
        for old_key, new_key in remapping_types.items():
            key_idx = self.section_data[:, 2] == old_key
            self.section_data[key_idx] = new_key

        if (np.abs(self.section_data[:, 2] - data[:, 1]) > 1e-12).any():
            raise ValueError(f"Internal error, non integer ID numbers detected ({swc_file})")

        # TODO: We need to remove all points within the soma. Move first point to soma boundary.
        #       Be careful, there might be multiple points inside the soma

        self.build_tree()

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
            for section_id in range(0, section_counter[section_type]+1):
                idx = np.where((self.section_data[:, 0] == section_id) & (self.section_data[:, 2] == section_type))[0]

                if len(idx) == 1:
                    self.section_data[idx, 1] = 0  # 1000 * 0.5 -- since soma parent to dendrites, set to 0
                    continue

                if not (np.diff(idx) == 1).all():
                    raise ValueError(f"Points on a {section_type} section must be consecutive")

                parent_idx = self.section_data[idx, 3]
                comp_length = np.linalg.norm(self.geometry[idx, :3] - self.geometry[parent_idx, :3], axis=1)

                # If parent compartment is soma, then we need to remove soma radius from the distance
                # because no part of the dendrite is inside the soma.
                if self.section_data[parent_idx[0], 2] == 1:

                    comp_length[0] -= self.geometry[parent_idx[0], 3]

                    if comp_length[0] < 0:
                        comp_length[0] = 1e-6  # set compartments inside soma to length 1
                        # raise ValueError(f"Internal error, compartment length {comp_length[0]} invalid.")

                self.section_data[idx, 1] = 1000 * np.cumsum(comp_length) / np.sum(comp_length)

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

    def clone(self, position, rotation, parent_tree_info=None):

        if self.position is not None or self.rotation is not None:
            raise ValueError("Not allowed to rotate or position a neuron that has already been rotated or positioned")

        new_md = MorphologyData()
        new_md.swc_file = self.swc_file
        new_md.geometry = self.geometry.copy()
        new_md.section_data = self.section_data.copy()
        new_md.sections = dict()

        for sec_type, sec_type_data in self.sections.items():
            new_md.sections[sec_type] = dict()

            for sec_key, sec_value in sec_type_data.items():
                new_md.sections[sec_type][sec_key] = sec_value.clone(new_morphology_data=new_md)

        for p_key, p_value in self.point_lookup.items():
            new_md.point_lookup[p_key] = p_value.copy()

        new_md.kd_tree_lookup = dict()
        new_md.parent_tree_info = parent_tree_info

        new_md.place(position=position, rotation=rotation)

        return new_md

    def place(self, position=None, rotation=None, parent_tree_info=None):

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

        if self.position is not None or self.rotation is not None:
            raise ValueError("Not allowed to rotate or position a neuron that has already been rotated or positioned")

        if isinstance(rotation, (list, tuple)):
            rotation = np.array(rotation)
        else:
            rotation = rotation.copy()

        self.rotation = rotation

        if isinstance(position, (list, tuple)):
            position = np.array(position)
        else:
            position = position.copy()

        self.position = position

        if rotation is not None:
            if not np.abs(np.linalg.det(rotation) - 1) < 1e-10 \
                    or not np.abs(np.matmul(rotation, rotation.T) - np.eye(3) < 1e-10).all():
                raise ValueError(f"Not a valid rotation matrix {rotation}")

            self.geometry[:, :3] = np.matmul(self.rotation, self.geometry[:, :3].T).T

        self.geometry[:, :3] += self.position

        if self.parent_tree_info is not None:
            # We need to update soma distance for subtree based on distance to parent
            # self.parent_tree = (MorphologyData, point_idx, arc_factor) -- attachment point

            parent_object, parent_point_idx, arc_factor = self.parent_tree
            parent_position = self.parent_object.geometry[parent_point_idx, :3]
            parent_soma_distance = self.parent_object.geometry[parent_point_idx, 4]

            dist_to_parent = np.linalg.norm(self.position - parent_position)
            self.geometry[:, 4] += parent_soma_distance + dist_to_parent * arc_factor

    def section_iterator(self, section_type=None):

        if section_type is None:
            for section_dict in self.sections.values():
                for section in section_dict:
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

    file_name = "/home/hjorth/HBP/BasalGangliaData/data/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026/morphology/WT-0728MSN01-cor-rep-ax-res3.swc"

    md = MorphologyData(swc_file=file_name)

    import pdb
    pdb.set_trace()