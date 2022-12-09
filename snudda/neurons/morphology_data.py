import os
import numpy as np

# TODO: Move constants like 1000 * sec_x to separate file


class SectionMetaData:

    """ Holds parent_id, children_id, points_id"""

    __slots__ = ["section_id", "parent_section_id", "parent_point_id",
                 "child_section_id", "point_range", "section_type", "morphology_data"]

    section_id: int
    parent_section_id: int
    parent_point_id: int
    child_section_id: dict
    point_range: np.ndarray
    section_type: int
    morphology_data: object

    def __init__(self, section_id, section_type, morphology_data):

        print(f"section_id = {section_id}, section_type = {section_type}")

        self.morphology_data = morphology_data
        self.section_id = section_id
        self.section_type = section_type

        idx = np.where((self.morphology_data.point_data[:, 0] == section_id)
                       & (self.morphology_data.point_data[:, 2] == section_type))[0]

        if len(idx) == 0:
            raise ValueError(f"Section id {section_id} has no points in morphology_data")

        if not (np.diff(idx) == 1).all():
            raise ValueError(f"Points on section must be consecutive")

        self.point_range = slice(idx[0], idx[-1] + 1)  # need + 1 at end for python slice to get inclusive
        self.parent_point_id = self.morphology_data.point_data[idx[0], 3]
        self.parent_section_id = self.morphology_data.point_data[self.parent_point_id, 2]

        # By definition only the last point in a section can be a parent to other sections
        child_idx = np.where(self.morphology_data.point_data[:, 3] == idx[-1])[0]

        self.child_section_id = dict()
        for child_id, child_type in zip(self.morphology_data.point_data[child_idx, 0],
                                        self.morphology_data.point_data[child_idx, 2]):
            if child_type not in self.child_section_id:
                self.child_section_id[child_type] = []

            self.child_section_id[child_type].append(child_id)

        for child_type in self.child_section_id.keys():
            self.child_section_id[child_type] = np.array(self.child_section_id[child_type])

        # Also check it above
        bastard_idx = np.where((idx[0] <= self.morphology_data.point_data[:, 3])
                               & (self.morphology_data.point_data[:, 3] < idx[-1]))[0]

        if not (bastard_idx == idx[1:]).all():
            raise ValueError(f"Only last point in section may have children outside section.")


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

        if not (data[0, 2:5] == [0, 0, 0]).all():
            raise ValueError(f"Does not have root centered at origo ({swc_file})")

        if data[0, 6] != -1:
            raise ValueError(f"First element must be root, and have parent ID -1 ({swc_file})")

        if not (data[:, 0] > data[:, 6]).all():
            raise ValueError(f"Parent ID must be lower than row ID ({swc_file})")

        item, count = np.unique(data[:, 0], return_counts=True)
        if (count > 1).any():
            raise ValueError(f"Duplicate index: {item[count > 1]} ({swc_file})")

        self.geometry = np.zeros((data.shape[0], 5), dtype=float)
        self.geometry[:, :4] = data[:, 2:6] * 1e-6  # x, y, z, r -- converted to meter

        # Calculate distance to soma and store in self.geometry
        parent_row_id = data[1:, 6].astype(int) - 1
        comp_length = np.linalg.norm(self.geometry[parent_row_id, :3] - self.geometry[1:, :3], axis=1)

        for comp_id, (parent_id, c_len) in enumerate(zip(parent_row_id, comp_length)):
            self.geometry[comp_id, 4] = self.geometry[parent_id, 4] + c_len

        # Store metadata for points
        self.point_data = np.full((data.shape[0], 4), -1, dtype=int)
        self.point_data[:, 2] = data[:, 1]
        self.point_data[0, 3] = -1
        self.point_data[1:, 3] = parent_row_id

        if (np.abs(self.point_data[:, 2] - data[:, 1]) > 1e-12).any():
            raise ValueError(f"Internal error, non integer ID numbers detected ({swc_file})")

        self.build_tree()

    def build_tree(self):

        # New sections are triggered when:
        # -- At branch points
        # -- Change of section type
        parent_id, counts = np.unique(self.point_data[:, 3], return_counts=True)
        branch_id = parent_id[counts > 1]
        type_switch_id = np.where(self.point_data[self.point_data[:, 3], 2] - self.point_data[:, 2] != 0)[0]
        type_switch_id = type_switch_id[type_switch_id != 0]
        edge_id = np.union1d(branch_id, self.point_data[type_switch_id, 3])
        edge_flag = np.zeros((self.point_data.shape[0],), dtype=bool)
        edge_flag[edge_id] = True

        section_counter = dict()

        # Assign section id to all points in point_data
        for idx, row in enumerate(self.point_data):
            section_type = row[2]
            parent_id = row[3]

            if parent_id == -1 or edge_flag[parent_id]:
                # Parent point is edge, create new section
                if section_type not in section_counter:
                    section_counter[section_type] = 0
                else:
                    section_counter[section_type] += 1

                section_id = section_counter.get(section_type)
                self.point_data[idx, 0] = section_id
            else:
                # Parent was not an edge, inherit section id
                self.point_data[idx, 0] = self.point_data[parent_id, 0]

        # Calculate section_x for all points in point_data
        for section_type in section_counter:
            for section_id in range(0, section_counter[section_type]):
                idx = np.where((self.point_data[:, 0] == section_id) & (self.point_data[:, 2] == section_type))[0]

                if len(idx) == 1:
                    self.point_data[idx, 1] = 1000 * 0.5
                    continue

                if not (np.diff(idx) == 1).all():
                    raise ValueError(f"Points on a {section_type} section must be consecutive")

                parent_idx = self.point_data[idx, 3]
                comp_length = np.linalg.norm(self.geometry[idx, :3] - self.geometry[parent_idx, :3], axis=1)

                # If parent compartment is soma, then we need to remove soma radius from the distance
                # because no part of the dendrite is inside the soma.
                if self.point_data[parent_idx[0], 3] == 1:
                    comp_length[0] -= self.geometry[parent_idx[0], 3]

                self.point_data[idx, 1] = 1000 * np.cumsum(comp_length) / np.sum(comp_length)

        # Build the actual tree
        self.sections = dict()
        for section_type in section_counter:
            if section_type not in self.sections:
                self.sections[section_type] = dict()

            for section_id in range(0, section_counter[section_type] + 1):
                self.sections[section_type][section_id] = SectionMetaData(section_id=section_id,
                                                                          section_type=section_type,
                                                                          morphology_data=self)

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


if __name__ == "__main__":

    file_name = "/home/hjorth/HBP/BasalGangliaData/data/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20211026/morphology/WT-0728MSN01-cor-rep-ax-res3.swc"

    md = MorphologyData(swc_file=file_name)

    import pdb
    pdb.set_trace()