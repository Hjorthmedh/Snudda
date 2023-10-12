import numpy as np
from scipy.spatial.transform import Rotation

from snudda.neurons import NeuronMorphologyExtended
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



    def rotate_subtree(self, section, from_idx, rotation, rotation_point):
        
        pass