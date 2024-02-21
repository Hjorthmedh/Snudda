import json
from collections import OrderedDict

import numpy as np
from scipy.interpolate import griddata

from snudda.utils import snudda_parse_path


class ProjectionDetection:

    def __init__(self, snudda_detect, rc=None, role=None):

        """

            Args:
                snudda_detect : SnuddaDetect object
                rc (default None): Ipyparallel remote client

        """

        self.snudda_detect = snudda_detect
        self.rc = rc
        self.d_view = None
        self.workers_initialised = False

        self.hyper_voxel_projections = dict()
        self.projections = dict()

        if not role:
            self.role = "master"
        else:
            self.role = role

        self.parse_config()

        if self.rc is not None:
            self.d_view = self.rc.direct_view(targets='all')
            self.setup_parallel()

    def setup_parallel(self):

        if self.d_view is None:
            return

        """ Prepares workers for parallel execution if d_view is not None. """

        assert self.role == "master", \
            "setup_parallel: Should only be called by master node"

        if self.workers_initialised:
            self.write_log("Workers already initialised.")
            return

        with self.d_view.sync_imports():
            from snudda.detect.projection_detection import ProjectionDetection

        cmd_str = "spd = ProjectionDetection(snudda_detect=sd, role='worker')"
        self.d_view.execute(cmd_str, block=True)
        cmd_str2 = "sd.projection_detection = spd"
        self.d_view.execute(cmd_str2, block=True)

        self.write_log(f"Workers set up.")

    def sync_projection_info(self):

        if self.role == "master" and self.d_view is not None:

            self.write_log("Synchronising projection info to workers.")
            for var_name in ["projections", "hyper_voxel_projections"]:
                self.d_view.push({"tmp_val": eval(f"self.{var_name}")}, block=True)
                cmd_str = f"spd.{var_name} = tmp_val"
                self.d_view.execute(cmd_str, block=True)

    def write_log(self, *args, **kvargs):
        # Reuse the log files for snudda_detect
        self.snudda_detect.write_log(*args, **kvargs)

    def find_neurons_projections_in_hyper_voxels(self):

        """ For each hyper voxel, list which neurons project to that hypervoxel"""

        if self.role != "master":
            # Master node sets this up
            return

        try:
            # In case the user has been naughty and added extra projections in the code.
            self.sync_projection_info()

            ss = np.random.SeedSequence(self.snudda_detect.random_seed + 202222)
            all_seeds = ss.generate_state(len(self.snudda_detect.neurons))
            all_neuron_id = sorted([n["neuron_id"] for n in self.snudda_detect.neurons])

            seed_lookup = dict()
            for s, n in zip(all_seeds, all_neuron_id):
                seed_lookup[n] = s

            proj_list = []

            self.write_log("!! about to loop projections")

            for neuron_type in self.projections:

                neuron_id_list = self.get_neurons_of_type(neuron_type)
                for neuron_id in neuron_id_list:
                    proj_list.append((neuron_id, neuron_type, seed_lookup[neuron_id]))

            self.hyper_voxel_projections = dict()

            if self.d_view is not None:
                self.d_view.scatter("proj_list", proj_list, block=True)

                cmd_str = "neuron_hv_list = spd.find_hyper_voxel_helper_parallel(proj_info=proj_list)"
                self.d_view.execute(cmd_str, block=True)

                neuron_hv_list = self.d_view.gather("neuron_hv_list", block=True)  # List of list of tuples

                for n_hv in neuron_hv_list:
                    neuron_id = n_hv[0]
                    neuron_hv = n_hv[1]

                    for hid in neuron_hv:
                        if hid not in self.hyper_voxel_projections:
                            self.hyper_voxel_projections[hid] = set([neuron_id])
                        else:
                            self.hyper_voxel_projections[hid].update([neuron_id])

            else:
                neuron_hv_list = self.find_hyper_voxel_helper_parallel(proj_list)

                for n_hv in neuron_hv_list:
                    neuron_id = n_hv[0]
                    neuron_hv = n_hv[1]

                    for hid in neuron_hv:
                        if hid not in self.hyper_voxel_projections:
                            self.hyper_voxel_projections[hid] = set([neuron_id])
                        else:
                            self.hyper_voxel_projections[hid].update([neuron_id])

            # Make sure the workers get the updated info
            self.sync_projection_info()

        except:
            # TODO: Remove this logging code
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str, is_error=True)
            print(t_str)
            import pdb
            pdb.set_trace()

    def find_hyper_voxel_helper_parallel(self, proj_info):

        neuron_hv = []

        for neuron_id, neuron_type, random_seed in proj_info:
            n_id, hyper_voxels = self.find_hyper_voxel_helper(neuron_id=neuron_id, neuron_type=neuron_type,
                                                              random_seed=random_seed)

            assert n_id == neuron_id
            neuron_hv.append((n_id, hyper_voxels))

        return neuron_hv

    def find_hyper_voxel_helper(self, neuron_id, neuron_type, random_seed):

        """ Returns hyper voxels that neuron_id (of neuron_type) projects to. """

        try:
            rng = np.random.default_rng(random_seed)

            hyper_voxels = set()

            if neuron_type in self.projections:
                for proj in self.projections[neuron_type]:
                    target_pos, target_rotation, axon_dist = proj["target_info"][neuron_id]
                    num_points = proj["num_points"]

                    if isinstance(proj["radius"], (np.ndarray, list)):
                        rx, ry, rz = proj["radius"]
                    else:
                        rx = ry = rz = proj["radius"]

                    # Find better way to calculate intersection between ellipsoid and hyper voxel cubes?
                    pos = self.ellipsoid_coordinates(target_pos, rx, ry, rz, target_rotation, num_points, rng)

                    hv_idx = ((pos - self.snudda_detect.simulation_origo)
                              / (self.snudda_detect.hyper_voxel_size * self.snudda_detect.voxel_size)).astype(int)

                    hv_id = map(lambda x, y, z: self.snudda_detect.hyper_voxel_id_lookup[x, y, z],
                                hv_idx[:, 0], hv_idx[:, 1], hv_idx[:, 2])

                    hyper_voxels.update(hv_id)

        except:
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str, is_error=True)

        return neuron_id, hyper_voxels

    @staticmethod
    def ellipsoid_coordinates(target_pos, rx, ry, rz, rotation, num_points, rng):

        """
            Return num points centred around target_pos, within an ellipsoid with axis rx, ry, rz rotated by rotation

            Args:
                target_pos : 3 x 1 numpy array
                rx, ry, rz : float values, axis of ellipsoid (before rotation)
                rotation : 3x3 rotation matrix, or None
                num_points : number of points to place
                rng : numpy random stream

            Returns:
                coordinates in meters
        """

        rnd_values = rng.random((num_points, 3))  # Use one call for all random values
        theta = 2 * np.pi * rnd_values[:, 0]
        phi = np.arccos(2 * rnd_values[:, 1] - 1)
        r_scale = rnd_values[:, 2] ** (1 / 3)

        x_coord = np.multiply(rx * r_scale, np.multiply(np.sin(phi), np.cos(theta)))
        y_coord = np.multiply(ry * r_scale, np.multiply(np.sin(phi), np.sin(theta)))
        z_coord = np.multiply(rz * r_scale, np.cos(phi))

        pos_offset = np.vstack((x_coord, y_coord, z_coord))

        if rotation is not None:
            pos_offset = np.matmul(rotation, pos_offset)

        pos = target_pos + pos_offset.T

        return pos

    def voxelise_projections(self):

        """ Add the projection of each neuron in pre_neuron_list to the axon space in currently active hyper voxel
         """

        if self.snudda_detect.hyper_voxel_id not in self.hyper_voxel_projections:
            # No projections in this hyper voxel, skip it
            return

        rng = self.snudda_detect.hyper_voxel_rng

        pre_neuron_list = np.sort([x for x in self.hyper_voxel_projections[self.snudda_detect.hyper_voxel_id]])

        for neuron_id in pre_neuron_list:
            neuron_type = self.snudda_detect.neurons[neuron_id]["type"]

            for proj in self.projections[neuron_type]:
                target_pos, target_rotation, axon_dist = proj["target_info"][neuron_id]
                num_points = proj["num_points"]
                if isinstance(proj["radius"], (np.ndarray, list)):
                    rx, ry, rz = proj["radius"]
                else:
                    rx = ry = rz = proj["radius"]

                pos = self.ellipsoid_coordinates(target_pos, rx, ry, rz, target_rotation, num_points, rng)

                # Convert to coordinates in the hyper voxel
                voxel_coords = np.round((pos - self.snudda_detect.hyper_voxel_origo)
                                        / self.snudda_detect.voxel_size).astype(int)

                inside_idx = np.where(np.sum(np.logical_and(0 <= voxel_coords,
                                                            voxel_coords < self.snudda_detect.hyper_voxel_size),
                                             axis=1) == 3)[0]

                for x, y, z in voxel_coords[inside_idx, :]:
                    if self.snudda_detect.axon_voxel_ctr[x, y, z] < self.snudda_detect.max_axon:
                        self.snudda_detect.axon_voxels[x, y, z, self.snudda_detect.axon_voxel_ctr[x, y, z]] = neuron_id
                        self.snudda_detect.axon_voxel_ctr[x, y, z] += 1
                        self.snudda_detect.axon_soma_dist[x, y, z] = axon_dist  # This is an underestimation
                    else:
                        self.snudda_detect.voxel_overflow_counter += 1

    def get_neurons_of_type(self, neuron_type):
        neuron_id = np.array([x["neuron_id"] for x in self.snudda_detect.neurons if x["type"] == neuron_type])
        return neuron_id

    def get_neuron_type(self, neuron_id):
        return self.snudda_detect.neurons[neuron_id]["type"]

    def get_neuron_positions(self, neuron_id):
        return self.snudda_detect.neuron_positions[neuron_id, :]

    def delete_projection(self, projection_name, pre_neuron_type):

        projection_list = self.projections[pre_neuron_type]
        new_projection_list = []
        for proj in projection_list:
            if proj["name"] != projection_name:
                new_projection_list.append(proj)

        self.projections[pre_neuron_type] = new_projection_list

    def parse_config(self):

        for region_name, region_data in self.snudda_detect.config["regions"].items():
            for con_name, con_info in region_data["connectivity"].items():
                for con_type, con_config in con_info.items():
                    if "projection_config_file" in con_config:
                        pre_neuron_type = con_name.split(",")[0]

                        if "projection_name" in con_config:
                            projection_name = con_config["projection_name"]
                        else:
                            projection_name = f"{con_name},{con_type}"

                        self.add_projection(projection_name=projection_name, pre_neuron_type=pre_neuron_type,
                                            projection_file=snudda_parse_path(con_config["projection_config_file"],
                                                                              snudda_data=self.snudda_detect.snudda_data))

    def add_projection(self, projection_name, pre_neuron_type, projection_file):

        """

            Args:
                projection_name (str) : Name of projection
                pre_neuron_type (str) : Neuron type of presynaptic neuron
                projection_file (str) : Path to projection file, used by griddata to find target region coordinates

            Projection mapping and rotations are specified in the projection_file.

        """

        with open(projection_file, "r") as f:
            projection_data = json.load(f, object_pairs_hook=OrderedDict)

        if projection_name and projection_name in projection_data:
            proj_file_info = projection_data[projection_name]
        else:
            proj_file_info = projection_data

        if "axon_morphology" in proj_file_info:
            self.write_log("Axon morphology projections not handled by projection_detection.py")
            self.write_log("UPDATING THE CODE, MAKE SURE detect.py DOES INCLUDE IT IN NEURON MORPHOLOGIES FOR NORMAL DETECTION")
            self.write_log("This code might become obsolete...")
            return

        assert "source" in proj_file_info, f"'source' must exist in {projection_file}"
        assert "destination" in proj_file_info, f"'destination' must exist in {projection_file}"
        assert "radius" in proj_file_info, f"'radius' must exist in {projection_file}"
        assert "num_points" in proj_file_info, f"'num_points' must exist in {projection_file}"

        proj_info = dict()
        proj_info["name"] = projection_name
        proj_info["source"] = np.array(proj_file_info["source"]) * 1e-6
        proj_info["dest"] = np.array(proj_file_info["destination"]) * 1e-6

        proj_info["radius"] = proj_file_info["radius"]  # Either single radius, or (rx, ry, rz)
        proj_info["num_points"] = proj_file_info["num_points"]  # Number of axon points placed

        if "rotation" in projection_data[projection_name]:
            proj_info["rotation"] = np.array(proj_file_info["rotation"])  # n x 9 matrix

        # Pre-compute target centres for each neuron of pre_neuron_type
        neuron_id = self.get_neurons_of_type(pre_neuron_type)
        pre_positions = self.get_neuron_positions(neuron_id)
        target_centres = griddata(points=proj_info["source"], values=proj_info["dest"],
                                  xi=pre_positions, method="linear")

        axon_dist = np.linalg.norm(target_centres - pre_positions, axis=1)

        if "rotation" in proj_info:
            target_rotation = griddata(points=proj_info["dest"], values=proj_info["rotation"],
                                       xi=target_centres, method="linear")
        else:
            target_rotation = [None for x in neuron_id]

        proj_info["target_info"] = dict()
        for nid, pos, rot, ad in zip(neuron_id, target_centres, target_rotation, axon_dist):
            if rot is not None:
                rot_mat = rot.reshape((3, 3))
                assert np.abs(np.linalg.det(rot_mat) - 1) < 1e-6, \
                    f"Invalid rotation matrix for {projection_name}: {rot_mat}"
                proj_info["target_info"][nid] = pos, rot_mat, ad
            else:
                proj_info["target_info"][nid] = pos, None, ad

        if pre_neuron_type not in self.projections:
            self.projections[pre_neuron_type] = [proj_info]
        else:
            self.projections[pre_neuron_type].append(proj_info)
