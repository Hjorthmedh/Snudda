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

    def __init__(self, swc_fil):

        self.geometry = None  # x, y, z, r, soma_dist (float)
        self.section_data = None  # section_id, section_x (*1000), section_type (int)
        self.sections = None  # dictionary section_id --> SectionMetaData

        self.point_lookup = dict()  # "dend" --> np.array of point_id for dend points
        self.section_lookup = dict()  # 'dend' --> np.array of section_id for dend sections

        self.rotation = None
        self.position = None

    def section_iterator(self, section_type):
        section_id = self.section_lookup[section_type]

        for sid in section_id:
            yield self.sections[sid]

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



