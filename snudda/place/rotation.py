import numpy as np
import json
from scipy.interpolate import griddata

from snudda.neurons.neuron_morphology import NeuronMorphology


class SnuddaRotate(object):

    def __init__(self, config_file=None, default_rotation_mode=None):

        self.rotation_lookup = dict()
        self.rotation_field_cache = dict()
        self.config = None

        if not default_rotation_mode:
            default_rotation_mode = SnuddaRotate.default_rotation

        self.rotation_lookup["default"] = default_rotation_mode

        if config_file:
            self.parse_config_file(config_file=config_file)

    def parse_config_file(self, config_file):

        with open(config_file, "r") as f:
            self.config = json.load(f)

        # Parse the config
        for volume in self.config["Volume"]:
            if "neuronOrientation" in self.config["Volume"][volume]:
                for neuron_types in self.config["Volume"][volume]["neuronOrientation"]:

                    orientation_info = self.config["Volume"][volume]["neuronOrientation"][neuron_types]
                    rotation_mode = orientation_info["rotationMode"]
                    vector_file = orientation_info["rotationFieldFile"] \
                        if "rotationFieldFile" in orientation_info else None

                    neuron_types = [x.strip().lower() for x in neuron_types.split(",")]
                    self.add_neuron_rotation(volume_name=volume, neuron_types=neuron_types,
                                             rotation_mode=rotation_mode, rotation_field_file=vector_file)

    @staticmethod
    def default_rotation(position, rng):
        return NeuronMorphology.rand_rotation_matrix(rand_nums=rng.random(size=(3,)))

    @staticmethod
    def no_rotation(position, rng):
        return np.eye(3)

    def add_neuron_rotation(self, volume_name: str, neuron_types, rotation_mode, rotation_field_file=None):

        for neuron_type in neuron_types:

            if rotation_mode == "random":
                self.rotation_lookup[volume_name, neuron_type] = SnuddaRotate.default_rotation

            elif "vectorField" in rotation_mode:
                assert rotation_field_file, \
                    f"You need to specify rotation_field_file (rotationFieldFile) for {neuron_type} in {volume_name}"

                position, rotation = self.load_rotation_field(rotation_field_file, volume_name)
                # Since our volume obj files are given in micrometers, we also specify positions in the
                # rotation file in micrometer, so correct for that when we read in so internally SI-units
                self.rotation_field_cache[volume_name, neuron_type] = np.array(position)*1e-6, np.array(rotation)

                if rotation_mode == "vectorFieldAndZ":
                    # Combine rotation around z-axis with rotation based on location in volume
                    self.rotation_lookup[volume_name, neuron_type] = \
                        lambda neuron_pos, rng: np.matmul(self.get_rotation(volume_name=volume_name,
                                                                            neuron_type=neuron_type,
                                                                            neuron_position=neuron_pos),
                                                          self.random_z_rotate(rng))
                else:
                    # Rotation based on position in volume
                    self.rotation_lookup[volume_name, neuron_type] = \
                        lambda neuron_pos, rng: self.get_rotation(volume_name=volume_name,
                                                                  neuron_type=neuron_type,
                                                                  neuron_position=neuron_pos)
            elif rotation_mode is None or rotation_mode == "" or rotation_mode == "none":
                    self.rotation_lookup[volume_name, neuron_type] = SnuddaRotate.no_rotation()
            else:
                assert False, f"Unknown rotation mode {rotation_mode}"

    @staticmethod
    def random_z_rotate(rng):

        ang = 2*np.pi*rng.uniform()

        return np.array([[np.cos(ang), -np.sin(ang), 0],
                         [np.sin(ang), np.cos(ang), 0],
                         [0, 0, 1]])

    def get_rotation(self, volume_name, neuron_type, neuron_position):

        field_position, field_rotation = self.rotation_field_cache[volume_name, neuron_type.lower()]

        rotation_vector = griddata(points=field_position,
                                   values=field_rotation,
                                   xi=neuron_position, method="linear")

        assert not np.isnan(np.sum(rotation_vector)), \
            (f"Invalid rotation vector for volume {volume_name}, neuron {neuron_type}, "
             f"is neuron position outside the field?\nNeuron position: {neuron_position}"
             f" (must be inside convex hull of the field's positions points)")

        # We need to rotate z-axis to point to rotation_vector
        rotation_matrix = SnuddaRotate.rotation_matrix_from_vectors(np.array([0, 0, 1]), rotation_vector)

        return rotation_matrix

    @staticmethod
    def rotation_matrix_from_vectors(vec1, vec2):

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

        with open(rotation_field_file, "r") as f:
            rotation_field_data = json.load(f)

        if volume_name in rotation_field_data:
            assert "position" in rotation_field_data[volume_name] \
                and "rotation" in rotation_field_data[volume_name], \
                f"Missing position and/or rotation tag in volume {volume_name}"

            return rotation_field_data[volume_name]["position"], rotation_field_data[volume_name]["rotation"]

        else:
            print(f"No volume name, assuming position and rotation is for {volume_name}")
            return rotation_field_data["position"], rotation_field_data["rotation"]

    def rotate_neuron(self, volume_name, neuron_type, position, rng):

        if (volume_name, neuron_type) in self.rotation_lookup:
            return self.rotation_lookup[volume_name, neuron_type](position, rng)
        else:
            return self.rotation_lookup["default"](position, rng)