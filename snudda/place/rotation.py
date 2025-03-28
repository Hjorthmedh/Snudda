import json

import numpy as np
from scipy.interpolate import griddata


class SnuddaRotate:
    """ Rotation object. """

    def __init__(self, config_file=None):

        """
        Constructor

        Args:
            config_file: Path to network config file
        """

        self.rotation_lookup = dict()
        self.config = None

        if config_file:
            self.parse_config_file(config_file=config_file)

    def parse_config_file(self, config_file):

        with open(config_file, "r") as f:
            self.config = json.load(f)

        # Parse the config
        for region_name, region_data in self.config["regions"].items():

            if "neuron_orientation" in region_data["volume"]:
                for neuron_type in region_data["volume"]["neuron_orientation"]:

                    orientation_info = region_data["volume"]["neuron_orientation"][neuron_type]
                    rotation_mode = orientation_info["rotation_mode"]
                    rotation_field_file = orientation_info["rotation_field_file"] if "rotation_field_file" in orientation_info else None
                    if rotation_field_file:
                        position, rotation = self.load_rotation_field(rotation_field_file, region_name)
                    else:
                        position, rotation = None, None

                    self.rotation_lookup[region_name, neuron_type] = (rotation_mode, position, rotation)

    @staticmethod
    def random_z_rotate(rng):
        """ Helper method, rotate around z-axis (of SWC coordinates)"""
        ang = 2 * np.pi * rng.uniform()

        return np.array([[np.cos(ang), -np.sin(ang), 0],
                         [np.sin(ang), np.cos(ang), 0],
                         [0, 0, 1]])

    def get_rotations(self, volume_name, neuron_type, neuron_positions, rng):
        """ Gets rotations for neuron_type in volumne name at neuron_positions """

        if (volume_name, neuron_type) in self.rotation_lookup:
            rotation_mode, field_position, field_rotation = self.rotation_lookup[volume_name, neuron_type]
        else:
            rotation_mode, field_position, field_rotation = "random", None, None

        if not rotation_mode or rotation_mode.lower() == "none":
            rotation_matrices = [np.eye] * neuron_positions.shape[0]

        elif rotation_mode in ["random", "default"]:
            rotation_matrices = [self.rand_rotation_matrix(rand_nums=rng.random(size=(3,)))
                                 for x in range(0, neuron_positions.shape[0])]

        elif "vector_field" in rotation_mode:
            rotation_vectors = griddata(points=field_position,
                                        values=field_rotation,
                                        xi=neuron_positions, method="linear")

            assert not np.isnan(np.sum(rotation_vectors)), \
                (f"Invalid rotation vector for volume {volume_name}, neuron {neuron_type}, "
                 f"is neuron position outside the field?\nNeuron positions: {neuron_positions}"
                 f" (must be inside convex hull of the field's positions points)")

            if rotation_mode == "vector_field_and_z":
                rotation_matrices = [np.matmul(SnuddaRotate.rotation_matrix_from_vectors(np.array([0, 0, 1]), rv),
                                               self.random_z_rotate(rng))
                                     for rv in rotation_vectors]
            else:
                # We need to rotate z-axis to point to rotation_vector
                rotation_matrices = [SnuddaRotate.rotation_matrix_from_vectors(np.array([0, 0, 1]), rv)
                                     for rv in rotation_vectors]
        else:
            raise TypeError(f"Unknown rotation mode {rotation_mode}")

        return rotation_matrices

    @staticmethod
    def rotation_matrix_from_vectors(vec1, vec2):
        """ Creates a rotation matrix to rotate vec1 into vec2."""

        # Special case, which gives cross product zero
        if (vec1 == vec2).all():
            return np.eye(3)

        # Taken from stackoverflow:
        # https://stackoverflow.com/questions/45142959/calculate-rotation-matrix-to-align-two-vectors-in-3d-space
        """ Find the rotation matrix that aligns vec1 to vec2
        :param vec1: A 3d "source" vector
        :param vec2: A 3d "destination" vector
        :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
        """
        a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
        v = np.cross(a, b)
        c = np.dot(a, b)
        s = np.linalg.norm(v)
        k_mat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + k_mat + k_mat.dot(k_mat) * ((1 - c) / (s ** 2))
        return rotation_matrix

    @staticmethod
    def test_rotation_matrix():

        vec1 = [2, 3, 2.5]
        vec2 = [-3, 1, -3.4]

        mat = SnuddaRotate.rotation_matrix_from_vectors(vec1, vec2)
        vec1_rot = mat.dot(vec1)
        assert np.allclose(vec1_rot / np.linalg.norm(vec1_rot), vec2 / np.linalg.norm(vec2))

    @staticmethod
    def load_rotation_field(rotation_field_file, volume_name):
        """ Loads rotation field for volumne_name from rotation_field_file """

        with open(rotation_field_file, "r") as f:
            rotation_field_data = json.load(f)

        if volume_name in rotation_field_data:
            assert "position" in rotation_field_data[volume_name] \
                   and "rotation" in rotation_field_data[volume_name], \
                f"Missing position and/or rotation tag in volume {volume_name}"

            return np.array(rotation_field_data[volume_name]["position"]) * 1e-6, \
                   np.array(rotation_field_data[volume_name]["rotation"])

        else:
            print(f"No volume name, assuming position and rotation is for {volume_name}")
            return np.array(rotation_field_data["position"]) * 1e-6, np.array(rotation_field_data["rotation"])

    @staticmethod
    def rand_rotation_matrix(deflection=1.0, rand_nums=None, rng=None):
        """
        Creates a random rotation matrix.

        deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
        rotation. Small deflection => small perturbation.
        rand_nums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
        """

        # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c

        if rand_nums is None:
            if rng is not None:
                rand_nums = rng.uniform(size=(3,))
            else:
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
