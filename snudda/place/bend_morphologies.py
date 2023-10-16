import numpy as np
from scipy.spatial.transform import Rotation

from snudda.neurons import NeuronMorphologyExtended
from snudda.neurons.morphology_data import MorphologyData, SectionMetaData
from snudda.place.region_mesh_redux import RegionMeshRedux


class BendMorphologies:

    def __init__(self, region_mesh: RegionMeshRedux, rng):

        self.region_mesh = region_mesh
        self.rng = rng

    def check_if_inside(self, morphology: NeuronMorphologyExtended):

        coords = morphology.morphology_data["neuron"].geometry
        inside_flag = self.region_mesh.check_inside(points=coords[:, :3])

        return inside_flag

    def bend_morphology(self, morphology: NeuronMorphologyExtended, k=50e-6):

        # TODO: Parent point idx is included if parent section is of same type as section!
        #       So we should not rotate the first point!

        # k -- decay constant
        n_random = 10
        candidate_pos = np.zeros((n_random, 3))

        parent_rotation_matrices = dict()  # indexed by parent_point_idx in SectionMetaData

        # Iterate over all parts of neuron
        for section in morphology.section_iterator():

            # We need to track the rotation of each point, in particular save rotations
            # for each branch point, so that all children can start with that rotation

            if section.parent_section_idx in parent_rotation_matrices:
                rotation_matrix = parent_rotation_matrices[section.parent_section_idx]
            else:
                rotation_matrix = np.eye(3)

            parent_coords = section.morphology_data.position[section.parent_section_idx, :]

            # Loop over all points in section
            for point_idx in section.point_idx:
                coords = section.morphology_data.geometry[point_idx, :]
                dist = self.region_mesh.distance_to_border(points=coords)

                P = 1 / (1 + np.exp(-k * dist))

                if self.rng.uniform(1) < P:
                    # We need to randomize new rotation matrix
                    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.html
                    angles = self.rng.uniform(size=(n_random, 3), low=-0.2, high=0.2)  # Angles in radians
                    rots = Rotation.from_euler(seq="XYZ", angles=angles)

                    for idx, rot in enumerate(rots):
                        candidate_pos[idx, :] = parent_coords[idx, :] + rots.apply(vectors=coords-parent_coords)

                    candidate_dist = self.region_mesh.distance_to_border(points=candidate_pos)

                    P_candidate = np.divide(1, 1 + np.exp(-k * candidate_dist))
                    picked_idx = self.rng.choice(n_random, p=P_candidate)

                    new_coords = candidate_pos[picked_idx, :]

                    picked_rotation = rots[picked_idx]
                    rotation_point = parent_coords

                    # TODO: We also need to apply this rotation to ALL the other points in this section,
                    #       and to all child sections of this branch.

            else:
                # Keep the old coords
                new_coords = coords



                # Old point is next points parent
                parent_coords = new_coords

            # Gradient is calculated based on distance to mesh
            # Calculate delta_gradient = gradient_self - gradient_parent
            # Amount of angle to bend f(delta_gradient)
            # Randomize the rotation matrix

            # After section is done, store the last rotation matrix, so its children can get their parent rotation
            # store in dictionary?

            parent_rotation_matrices[point_idx] = rotation_matrix

    def get_full_rotation_representation(self, morphology: MorphologyData):

        rotation_representation = dict()
        parent_direction = dict()

        for section in morphology.section_iterator():

            if (section.section_id, section.section_type) in parent_direction:
                parent_dir = parent_direction[section.section_id, section.section_type]
            else:
                parent_dir = np.array([[1, 0, 0]])

            try:
                rot_and_len, last_direction = self.rotation_representation(section=section, parent_direction=parent_dir)
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

            rotation_representation[section.section_id, section.section_type] = rot_and_len

            for child_id, child_type in section.child_section_id.T:
                parent_direction[child_id, child_type] = last_direction

        return rotation_representation

    def apply_rotation(self, morphology: MorphologyData, rotation_representation):

        parent_direction = dict()

        # new_coords = np.zeros((morphology.geometry.shape[0], 3))
        new_coords = np.full((morphology.geometry.shape[0], 3), np.nan)  # Use nan as fill value, to see problems

        for section in morphology.section_iterator():
            if (section.section_id, section.section_type) in parent_direction:
                parent_dir, parent_pos = parent_direction[section.section_id, section.section_type]
            else:
                if morphology.rotation:
                    parent_dir = np.matmul(morphology.rotation, np.array([[1, 0, 0]]))
                else:
                    parent_dir = np.array([[1, 0, 0]])

                if morphology.position:
                    parent_pos = morphology.position
                else:
                    parent_pos = np.zeros((3, ))

            if section.section_type == section.parent_section_type or section.section_type == 1:
                # The section includes the parent point, if both sections are of the same type
                # or if the section is the soma (since then it has no rotations, so position gets passed as parent point)
                include_parent_point = True
            else:
                include_parent_point = False

            rot_rep = rotation_representation[section.section_id, section.section_type]
            coords, last_dir = self.coordinate_representation(rotation_representation=rot_rep,
                                                              parent_direction=parent_dir,
                                                              parent_point=parent_pos,
                                                              return_last_direction=True,
                                                              include_parent_point=include_parent_point)

            new_coords[section.point_idx, :3] = coords

            for child_id, child_type in section.child_section_id.T:
                parent_direction[child_id, child_type] = (last_dir, coords[-1, :3])

        # TODO: This just returns the coords for now, add option to update coords in morphology?
        #       OBS! Then rotation should also be reset, since it is now included in the coordinates

        return new_coords

    def rotation_representation(self, section: SectionMetaData, parent_direction=None):

        """ Represent each section as a series of length of segment, and rotations relative the parent segment."""

        if parent_direction is None:
            parent_direction = np.array([[1, 0, 0]])

        rotations_and_length = []
        parent_direction = parent_direction / np.linalg.norm(parent_direction)
        segment_direction = parent_direction  # To handle soma...

        # If parent compartment is of a different type (e.g. soma parent, for the dendrite) then we need
        # to make sure that the first point is also included. So we need to artificially add the soma.

        if section.section_type != section.parent_section_type and section.parent_section_type != -1:
            coords = np.zeros((len(section.point_idx) + 1, 3))
            coords[0, :3] = section.morphology_data.geometry[section.parent_point_idx, :3]
            coords[1:, :3] = section.morphology_data.geometry[section.point_idx, :3]
        else:
            coords = section.morphology_data.geometry[section.point_idx, :3]

        delta = np.diff(coords, axis=0)
        delta_length = np.linalg.norm(delta, axis=1)
        delta_direction = delta / delta_length.reshape((delta.shape[0], 1))

        for segment_direction, segment_length in zip(delta_direction, delta_length):
            segment_direction = segment_direction.reshape((1, 3))
            try:
                rotation, _ = Rotation.align_vectors(segment_direction, parent_direction)
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()


            rotations_and_length.append((rotation, segment_length))
            parent_direction = segment_direction

        return rotations_and_length, segment_direction

    def coordinate_representation(self, rotation_representation,
                                  parent_direction=None,
                                  parent_point=None,
                                  return_last_direction=False,
                                  include_parent_point=True):

        if parent_direction is None:
            parent_direction = np.array([1, 0, 0])

        if parent_point is None:
            parent_point = np.zeros(3)

        parent_direction = parent_direction / np.linalg.norm(parent_direction)

        coords = np.zeros((len(rotation_representation)+include_parent_point, 3))

        if include_parent_point:
            coords[0, :] = parent_point

        for idx, (rotation, length) in enumerate(rotation_representation):
            segment_direction = rotation.apply(parent_direction)
            parent_point = coords[idx+include_parent_point, :] = segment_direction * length + parent_point
            parent_direction = segment_direction

        if return_last_direction:
            return coords, parent_direction

        return coords


def test_rotation_representation():

    file_path = "../data/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20190508/WT-0728MSN01-cor-rep-ax.swc"

    md = MorphologyData(swc_file=file_path)

    bm = BendMorphologies(None, rng=np.random.default_rng())
    # sec = md.sections[3][0]
    # rot, _ = bm.rotation_representation(sec)
    # coords = bm.coordinate_representation(rotation_representation=rot, parent_point=sec.position[0, :])

    rot_rep = bm.get_full_rotation_representation(morphology=md)
    coords = bm.apply_rotation(morphology=md, rotation_representation=rot_rep)

    if not (np.abs(md.geometry[:, :3] - coords) < 1e-10).all():
        print(f"Geometry mismatch, problem with representation.")
        import pdb
        pdb.set_trace()
    else:
        print(f"Geometry matches")


if __name__ == "__main__":

    test_rotation_representation()

    import pdb
    pdb.set_trace()