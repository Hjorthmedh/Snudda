import numpy as np
import json
from collections import OrderedDict
from scipy.interpolate import griddata


class ProjectionDetection:

    def __init__(self, snudda_detect):

        self.snudda_detect = snudda_detect

        self.projections = dict()

    def neurons_projections_in_hyper_voxels(self):

        """ For each hyper voxel, list which neurons project to that hypervoxel"""

        # Good way to calculate intersection between elipsoid and hyper voxel cubes?
        pass

    def voxelise_projections(self, pre_neuron_list, rng):

        """ Add the projection of each neuron in pre_neuron_list to the axon space in currently active hyper voxel
         """

        for neuron_id in pre_neuron_list:
            neuron_type = self.snudda_detect.neurons[neuron_id]["type"]

            for proj in self.projections[neuron_type]:
                target_pos, target_rotation, axon_dist = proj["target_info"]
                num_points = proj["num_points"]
                if len(proj["radius"]) == 1:
                    rx = ry = rz = proj["radius"]
                else:
                    rx, ry, rz = proj["radius"]

                theta = 2 * np.pi * rng.random(num_points)
                phi = np.arccos(2 * rng.random(num_points) - 1)
                r_scale = rng.random(num_points) ** (1 / 3)

                x_coord = np.multiply(rx * r_scale, np.multiply(np.sin(phi), np.cos(theta)))
                y_coord = np.multiply(ry * r_scale, np.multiply(np.sin(phi), np.sin(theta)))
                z_coord = np.multiply(rz * r_scale, np.cos(phi))

                pos = np.hstack((x_coord, y_coord, z_coord))

                if "rotation" in proj:
                    pos = np.matmul(proj["rotation"], pos)

                # Convert to coordinates in the hyper voxel
                voxel_coords = np.round((pos - self.snudda_detect.hyper_voxel_origo) / self.snudda_detect.voxel_size)

                inside_idx = np.where(np.sum(0 <= voxel_coords < self.snudda_detect.hyper_voxel_size, axis=1) == 3)[0]

                for x, y, z in voxel_coords[inside_idx, :]:
                    if self.snudda_detect.axon_voxel_ctr[x, y, z] < self.snudda_detect.max_axon:
                        self.snudda_detect.axon_voxels[x, y, z, self.snudda_detect.axon_voxel_ctr[x, y, z]] = neuron_id
                        self.snudda_detect.axon_voxel_ctr[x, y, z] += 1
                        self.snudda_detect.axon_soma_dist[x, y, z] = axon_dist  # This is an underestimation
                    else:
                        self.snudda_detect.voxel_overflow_ctr += 1

    def get_neurons_of_type(self, neuron_type):
        neuron_id = np.array([x["neuronID"] for x in self.snudda_detect.neurons if x["type"] == neuron_type])
        return neuron_id

    def get_neuron_type(self, neuron_id):
        return self.snudda_detect.neurons[neuron_id]["type"]

    def get_neuron_positions(self, neuron_id):
        return self.neuron_positions[neuron_id, :]

    def delete_projection(self, projection_name, pre_neuron_type):

        projection_list = self.projections[pre_neuron_type]
        new_projection_list = []
        for proj in projection_list:
            if proj["name"] != projection_name:
                new_projection_list.append(proj)

        self.projections[pre_neuron_type] = new_projection_list

    def add_projection(self, projection_name, pre_neuron_type, projection_file, projection_radius=30e-6,
                       num_points=50):

        """

            Args:
                projection_name (str) : Name of projection
                pre_neuron_type (str) : Neuron type of presynaptic neuron
                projection_file (str) : Path to projection file, used by griddata to find target region coordinates

            Projection mapping and rotations are specified in the projection_file.

        """

        with open(projection_file, "r") as f:
            projection_data = json.load(f, object_pairs_hook=OrderedDict)

        if projection_name:
            proj_file_info = projection_data[projection_name]
        else:
            proj_file_info = projection_data

        proj_info = dict()
        proj_info["name"] = projection_name
        proj_info["source"] = np.array(proj_file_info["source"]) * 1e-6
        proj_info["dest"] = np.array(proj_file_info["destination"]) * 1e-6

        proj_info["radius"] = proj_file_info["radius"]  # Either single radius, or (rx, ry, rz)
        proj_info["num_points"] = proj_file_info["numPoints"]  # Number of axon points placed

        if "rotation" in projection_data[projection_name]:
            proj_info["rotation"] = np.array(proj_file_info["rotation"])  # n x 9 matrix

        # Pre-compute target centres for each neuron of pre_neuron_type
        neuron_id = self.get_neurons_of_type(pre_neuron_type)
        pre_positions = self.get_neuron_positions(neuron_id)
        target_centres = griddata(points=proj_info["source"], values=proj_info["dest"],
                                  xi=pre_positions, method="linear")

        axon_dist = np.linalg.norm(target_centres - pre_positions, axis=1)

        if "rotation" in proj_info:
            target_rotation = griddata(points=proj_info["dest"], values=["rotation"],
                                       xi=target_centres, method="linear")
        else:
            target_rotation = [None for x in neuron_id]

        proj_info["target_info"] = dict()
        for nid, pos, rot, ad in zip(neuron_id, target_centres, target_rotation, axon_dist):
            if rot is not None:
                proj_info["target_info"][nid] = pos, rot.reshape((3, 3)), ad
            else:
                proj_info["target_info"][nid] = pos, None, ad

        proj_info["radius"] = projection_radius

        if pre_neuron_type not in self.projections:
            self.projections[pre_neuron_type] = [proj_info]
        else:
            self.projections[pre_neuron_type].append(proj_info)


