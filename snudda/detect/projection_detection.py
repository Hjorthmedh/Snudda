import numpy as np
import json
from collections import OrderedDict
from scipy.interpolate import griddata


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

        if self.rc is not None:
            self.d_view = self.rc.direct_view(targets='all')

        if not role:
            self.role = "master"
        else:
            self.role = role

        self.projections = dict()

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

        self.d_view.push({"random_seed": self.random_seed}, block=True)

        cmd_str = "spd = ProjectionDetection(snudda_detect=sd, role='worker', random_seed=random_seed)"
        self.d_view.execute(cmd_str, block=True)

        self.write_log(f"Workers set up.")

    def sync_projection_info(self):

        if self.role == "master" and self.dview is not None:
            self.write_log("Synchronising projection info to workers.")
            self.dview.push({"proj_info": self.projections}, block=True)
            cmd_str = "spd.projections = proj_info"
            self.dview.execute(cmd_str, block=True)

    def write_log(self, *args, **kvargs):
        # Reuse the log files for snudda_detect
        self.snudda_detect.write_log(*args, **kvargs)

    def neurons_projections_in_hyper_voxels(self):

        if self.role != "master":
            # Master node sets this up
            return

        self.sync_projection_info()

        """ For each hyper voxel, list which neurons project to that hypervoxel"""

        ss = np.random.SeedSequence(self.snudda_detect.random_seed + 202222)
        all_seeds = ss.generate_state(len(self.snudda_detect.neurons))
        all_neuron_id = sorted([n["neuronID"] for n in self.snudda_detect.neurons])

        seed_lookup = dict()
        for s, n in zip(all_seeds, all_neuron_id):
            seed_lookup[n] = s

        proj_list = []

        for neuron_type in self.projections:

            neuron_id = self.get_neurons_of_type(neuron_type)
            proj_list.append((neuron_id, neuron_type, seed_lookup[neuron_id]))

        self.hyper_voxel_projections = set()

        if self.d_view is not None:
            self.d_view.scatter("proj_list", proj_list, block=True)
            self.d_view.execute("neuron_hv_list = self.find_hyper_voxel_helper_parallel(proj-info=proj_list",
                                block=True)
            neuron_hv_list = self.d_view.gather("neuron_hv_list", block=True)  # List of list of tuples

            for nl_worker in neuron_hv_list:
                for n_hv in nl_worker:
                    neuron_id = n_hv[0]
                    neuron_hv = n_hv[1]

                    for hid in neuron_hv:
                        if hid not in self.hyper_voxel_projections:
                            self.hyper_voxel_projections = set(neuron_id)
                        else:
                            self.hyper_voxel_projections[hid].add(neuron_id)

        else:
            neuron_hv_list = map(self.find_hyper_voxel_helper_parallel, proj_list)

            for n_hv in neuron_hv_list:
                neuron_id = n_hv[0]
                neuron_hv = n_hv[1]

                for hid in neuron_hv:
                    if hid not in self.hyper_voxel_projections:
                        self.hyper_voxel_projections = set(neuron_id)
                    else:
                        self.hyper_voxel_projections[hid].add(neuron_id)

        # Better way to calculate intersection between ellipsoid and hyper voxel cubes?

    def find_hyper_voxel_helper_parallel(self, proj_info):
        neuron_hv = []

        for neuron_id, neuron_type, random_seed in proj_info:
            neuron_hv.append(self.find_hyper_voxel_helper(neuron_id=neuron_id, neuron_type=neuron_type,
                                                          random_seed=random_seed))

        return neuron_hv

    def find_hyper_voxel_helper(self, neuron_id, neuron_type, random_seed):

        """ Returns hyper voxels that neuron_id (of neuron_type) projects to. """

        rng = np.random.default_rng(random_seed)

        hyper_voxels = set()

        for proj in self.projections[neuron_type]:
            target_pos, target_rotation, axon_dist = proj["target_info"][neuron_id]
            num_points = proj["num_points"]
            if len(proj["radius"]) == 1:
                rx = ry = rz = proj["radius"]
            else:
                rx, ry, rz = proj["radius"]

            if "rotation" in proj:
                rotation = proj["rotation"]
            else:
                proj = None

            pos = self.ellipsoid_coordinates(target_pos, rx, ry, rz, rotation, num_points, rng)

            hv_idx = ((pos - self.snudda_detect.simulation_origo)
                      / (self.snudda_detect.hyper_voxel_size * self.snudda_detect.voxel_size)).astype(int)

            hv_id = map(self.snudda_detect.hyper_voxel_id_lookup, hv_idx[:, 0], hv_idx[:, 1], hv_idx[:, 2])

            hyper_voxels += set(hv_id)

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

        pos_offset = np.hstack((x_coord, y_coord, z_coord))

        if rotation is not None:
            pos_offset = np.matmul(rotation, pos_offset)

        pos = target_pos + pos_offset

        return pos

    def voxelise_projections(self, pre_neuron_list, rng):

        """ Add the projection of each neuron in pre_neuron_list to the axon space in currently active hyper voxel
         """

        for neuron_id in pre_neuron_list:
            neuron_type = self.snudda_detect.neurons[neuron_id]["type"]

            for proj in self.projections[neuron_type]:
                target_pos, target_rotation, axon_dist = proj["target_info"][neuron_id]
                num_points = proj["num_points"]
                if len(proj["radius"]) == 1:
                    rx = ry = rz = proj["radius"]
                else:
                    rx, ry, rz = proj["radius"]

                if "rotation" in proj:
                    rotation = proj["rotation"]
                else:
                    proj = None

                pos = self.ellipsoid_coordinates(target_pos, rx, ry, rz, rotation, num_points, rng)

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
                num_points (int) : Number of axon voxels placed

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


