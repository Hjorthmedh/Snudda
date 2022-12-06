import os
import numpy as np


class SectionMetaData:

    """ Holds parent_id, children_id, points_id"""

    __slots__ = ["parent_id", "children_id", "points_id"]

    section_id: int
    parent_id: int
    children_id: list[int]
    point_id: np.ndarray[int]
    section_type: str

    def __init__(self, section_id, parent_id, point_id, children_id=None):

        self.section_id = section_id
        self.parent_id = parent_id
        self.point_id = point_id

        if children_id is not None:
            self.children_id = children_id
        else:
            self.children_id = []


class MorphologyData:

    """
        MorphologyData

        This can hold an entire neuron, or a part of a neuron.

    """

    def __init__(self, swc_file=None):

        self.swc_file = swc_file

        self.geometry = None      # x, y, z, r, soma_dist (float)
        self.point_data = None    # section_id, section_x (*1000), section_type (int), parent_id (int)
        self.sections = None      # dictionary section_id --> SectionMetaData

        self.point_lookup = dict()    # "dend" --> np.array of point_id for dend points
        self.section_lookup = dict()  # 'dend' --> np.array of section_id for dend sections

        self.rotation = None
        self.position = None

        if swc_file is not None:
            self.load_swc_file(swc_file=swc_file)

    def section_iterator(self, section_type):
        section_id = self.section_lookup[section_type]

        for sid in section_id:
            yield self.sections[sid]

    def load_swc_file(self, swc_file):
        # This function is not SNUDDA_DATA aware, the file must exist

        if not os.path.isfile(swc_file):
            raise FileNotFoundError(f"Missing SWC file '{swc_file}'")

        data = np.loadtxt(swc_file)

        if any(np.diff(data[:, 0]) != 1):
            raise IndexError(f"SWC file has gaps in ID numbering ({swc_file})")

        if data[0, 0] != 1:
            raise IndexError(f"ID does not start from 1 ({swc_file})")

        if not data[0, 2:5] == [0, 0, 0]:
            raise ValueError(f"Does not have root centered at origo ({swc_file})")

        if data[0, 6] != -1:
            raise ValueError(f"First element must be root, and have parent ID -1 ({swc_file})")

        if not (data[:, 0] > data[:, 6]).all():
            raise ValueError(f"Parent ID must be lower than row ID ({swc_file})")

        self.geometry = np.zeros((data.shape[0], 5), dtype=float)
        self.geometry[:, :4] = data[:, 2:6] * 1e-6  # x, y, z, r -- converted to meter

        # Calculate distance to soma and store in self.geometry
        parent_row_id = data[:, 6] - 1
        comp_length = np.linalg.norm(self.geometry[parent_row_id, :3] - self.geometry[1:, :3], axis=1)

        for comp_id, parent_id, c_len in enumerate(zip(parent_row_id, comp_length)):
            self.geometry[comp_id, 5] = self.geometry[parent_id, 5] + comp_length

        # Store metadata for points
        self.point_data = np.full((data.shape[0], 4), -1, dtype=int)
        self.point_data[:, 2] = data[:, 1]
        self.point_data[:, 3] = parent_row_id

        if (np.abs(self.point_data[:, 2] - data[:, 1]) > 1e-12).any():
            raise ValueError(f"Internal error, non integer ID numbers detected ({swc_file})")

        self.build_tree()

    def build_tree(self):

        # Find branch and leaves
        child_count = np.array((self.geometry.shape[0],), dtype=int)

        parent_id, counts = np.unique(self.point_data[:, 3], return_counts=True)

        branch_id = parent_id[counts > 1]
        leaf_id = parent_id[counts == 0]

        edge_mask = np.zeros((self.geometry.shape[0],), dtype=bool)
        edge_mask[branch_id] = True
        edge_mask[leaf_id] = True

        # TO BE CONTINUED HERE
        assert False


        for parent in self.point_data[:, 3]:

        pass


    def clone(self):
        # Implement
        pass

    def place(self, position, rotation=None):

        # Here we assume soma is only a point
        soma_position = self.geometry[self.section_lookup["soma"], :3]
        if not (soma_position == 0).all():
            raise ValueError("Soma must be centered at origo before placement.")

        if self.position is not None and self.rotation is not None:
            raise ValueError("Not allowed to rotate or position a neuron that has already been rotated or positioned")

        self.rotation = rotation
        self.position = position

        if rotation is not None:
            self.geometry[:, :3] = np.matmul(self.rotation, self.geometry[:, :3].T).T + self.position



