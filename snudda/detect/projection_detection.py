import json
import os
from collections import OrderedDict

import numpy as np
from scipy.interpolate import griddata
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
from snudda.utils import snudda_parse_path
from snudda.place.region_mesh_redux import RegionMeshRedux


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
       
        #self.n_terminals_gmm = self.load_gmm_from_json(os.path.join('/Users/wst/Desktop/Karolinska/ReferenceData/GMMs', 'n_terminals_gmm.json'))
        #self.mst_gmm = self.load_gmm_from_json(os.path.join('/Users/wst/Desktop/Karolinska/ReferenceData/GMMs', 'mst_gmm.json'))
        
        #self.sec_gmm = self.load_gmm_from_json(os.path.join('/Users/wst/Desktop/Karolinska/ReferenceData/GMMs', 'sec_gmm.json'))
        #self.branches_gmm = self.load_gmm_from_json( os.path.join('/Users/wst/Desktop/Karolinska/ReferenceData/GMMs', 'branches_gmm.json'))
        
        self.synth_coords = dict()
        #self.target_mesh = RegionMeshRedux(mesh_path= '/Users/wst/Desktop/Karolinska/Meshes/SNr/SNr_total.obj')
        
        
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
                    #pos = self.ellipsoid_coordinates(target_pos, rx, ry, rz, target_rotation, num_points, rng, jitter = True)

                    pos = self.synthesize_coordinates(target_pos)
                    print(len(pos))
                    self.synth_coords[neuron_id] = pos

                    hv_idx = ((pos - self.snudda_detect.simulation_origo)/ (self.snudda_detect.hyper_voxel_size * self.snudda_detect.voxel_size)).astype(int)

                    hv_id = map(lambda x, y, z: self.snudda_detect.hyper_voxel_id_lookup[x, y, z],
                                hv_idx[:, 0], hv_idx[:, 1], hv_idx[:, 2])
                    hyper_voxels.update(hv_id)

        except:
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str, is_error=True)

        return neuron_id, hyper_voxels

    @staticmethod
    def ellipsoid_coordinates(target_pos, rx, ry, rz, rotation, num_points, rng, jitter = True):

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
        
        if jitter:
            rotation = rotation * rng.uniform(0.90, 1.10)
            rx = rx * rng.uniform(0.7, 1.1)
            ry = ry * rng.uniform(0.8, 1.50)
            rz = rz * rng.uniform(0.8, 1.50)
        
        x_coord = np.multiply(rx * r_scale, np.multiply(np.sin(phi), np.cos(theta)))
        y_coord = np.multiply(ry * r_scale, np.multiply(np.sin(phi), np.sin(theta)))
        z_coord = np.multiply(rz * r_scale, np.cos(phi))

        pos_offset = np.vstack((x_coord, y_coord, z_coord))

        if rotation is not None:
            pos_offset = np.matmul(rotation, pos_offset)

        pos = target_pos + pos_offset.T

        return pos

    # def synthesize_coordinates(self, target_pos):
        
    #     # current = np.array([0.0075, target_pos[1], target_pos[2]])
    #     # current = np.array([8200e-6, 5000e-6, 7300e-6])
    #     current = target_pos - np.array([500e-6, 0, 0]) 
    #     synth_terminals = []
    #     n = max(2, int(round(self.n_terminals_gmm.sample(1)[0][0][0])))
    #     print('Number of terminals: ' + str(n))
        
    #     sampled_vectors = self.mst_gmm.sample(n)[0] * 1e-6
    #     reference_vector = np.mean(sampled_vectors, axis = 0)
    #     angles = [self.angle_with_reference(vec, reference_vector) for vec in sampled_vectors]
    #     sorted_indices = np.argsort(angles)
    #     sorted_vectors = sampled_vectors[sorted_indices]
    #     sampled_vectors = self.partial_shuffle(sorted_vectors, fraction=0.2)

    #     v = 0
    #     while len(synth_terminals) < len(sampled_vectors): 
    
    #         new_terminal = current + sampled_vectors[v]
    #         synth_terminals.append(new_terminal)
    #         current = new_terminal
    #         if v > 1: 
    #             # if not self.bounding_box_check(current, min_point, max_point):
    #             #     current = synth_terminals[np.random.choice(len(synth_terminals))]
    #             if self.target_mesh.distance_to_border(current.reshape(1, -1)) > 25e-6:
    #                 current = synth_terminals[np.random.choice(len(synth_terminals))]
    #         v+=1
             
    #     return np.array(synth_terminals)
    
    
    '''
    def synthesize_coordinates(self, target_pos):
        
        n = max(1, int(round(self.branches_gmm.sample(1)[0][0][0])))
        
        p = np.array(target_pos) - np.array([500e-6, 0, 0])   ## p = np.array([7500e-6, target_pos[1], target_pos[2]])
        branch_points = [(p, 0)]
        
        vs = self.sec_gmm.sample(n)[0]*1e-6
        # np.random.shuffle(vs)
        vs = self.partial_shuffle(vs, fraction=0.2)
        lengths = np.linalg.norm(vs, axis =1)
        
        all_coords = []
        i = 0
        depth = 0

        while i < n:
            v = vs[i]
            all_coords.extend(self.generate_tortuous_path(start = p, end = p+v, num_points = int(round(lengths[i]/3e-6)), tortuosity = lengths[i]/20))
            depth +=1
            if lengths[i] > 25e-6:
                p = all_coords[-1]
                branch_points.append((p, depth))

            if depth >= 7:
                bp_idx = np.random.choice(range(len(branch_points)))
                p = branch_points[bp_idx][0]
                depth = branch_points[bp_idx][1]

            if self.target_mesh.distance_to_border(p.reshape(1, -1)) > 25e-6:
                bp_idx = np.random.choice(range(len(branch_points)))
                p = branch_points[bp_idx][0]
                depth = branch_points[bp_idx][1]

            i+=1
        all_coords = self.jitter_coordinates(np.array(all_coords), scale = 1e-6)

        return all_coords
    
    def resample_coordinates(self, coords, resolution):
        """
        Resample a list of coordinates to a specified resolution.

        Args:
            coords (list of tuple): List of (x, y, z) coordinates.
            resolution (float): Desired spacing between resampled points.

        Returns:
            list of tuple: Resampled coordinates.
        """
        resampled_coords = []
        
        for i in range(len(coords) - 1):
            # Get the current and next point
            p1 = np.array(coords[i])
            p2 = np.array(coords[i + 1])
            
            # Calculate the distance and direction vector
            distance = np.linalg.norm(p2 - p1)
            direction = (p2 - p1) / distance
            
            # Add the starting point of the segment
            if not resampled_coords or not np.allclose(resampled_coords[-1], p1):
                resampled_coords.append(tuple(p1))
            
            # Interpolate points along the segment
            num_points = int(np.floor(distance / resolution))
            for j in range(1, num_points + 1):
                new_point = p1 + direction * j * resolution
                resampled_coords.append(tuple(new_point))
        
        # Add the last point
        resampled_coords.append(tuple(coords[-1]))
        
        return resampled_coords

    def jitter_coordinates(self, points, scale= 3e-6):
        """
        Add random noise (jitter) to 3D coordinates.

        Parameters:
        - points (np.ndarray): Array of shape (N, 3), where N is the number of 3D points.
        - scale (float): The standard deviation of the jitter noise.

        Returns:
        - np.ndarray: Jittered coordinates of the same shape as `points`.
        """
        # Ensure the points are in a numpy array
        points = np.asarray(points)
        
        # Generate random noise
        noise = np.random.normal(loc=0.0, scale=scale, size=points.shape)
        
        # Add noise to the original points
        jittered_points = points + noise
        
        return jittered_points
    
    def generate_tortuous_path(self, start, end, num_points, tortuosity=0.1):
        """
        Generate a tortuous path between two points in 3D.
    
        Args:
            start (tuple): Starting coordinate (x, y, z).
            end (tuple): Ending coordinate (x, y, z).
            num_points (int): Number of points in the path.
            tortuosity (float): Factor controlling the amount of deviation from a straight line.
    
        Returns:
            numpy.ndarray: Array of shape (num_points, 3) with the path coordinates.
        """
        # Create a straight line path
        straight_path = np.linspace(start, end, num_points)
        
        # Calculate the direction vector of the straight line
        direction = np.array(end) - np.array(start)
        
        # Generate two orthogonal vectors to the direction
        if np.allclose(direction, [0, 0, 0]):
            raise ValueError("Start and end points cannot be the same.")
        arbitrary_vector = np.array([1, 0, 0]) if direction[0] == 0 else np.array([0, 1, 0])
        normal1 = np.cross(direction, arbitrary_vector)
        normal1 /= np.linalg.norm(normal1)
        normal2 = np.cross(direction, normal1)
        normal2 /= np.linalg.norm(normal2)
        
        # Add sinusoidal or random deviations for tortuosity
        deviations1 = np.sin(np.linspace(0, 2 * np.pi, num_points)) * tortuosity
        deviations2 = np.cos(np.linspace(0, 2 * np.pi, num_points)) * tortuosity
        tortuous_path = straight_path + deviations1[:, None] * normal1 + deviations2[:, None] * normal2
        
        return tortuous_path
    '''

        
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

                #pos = self.ellipsoid_coordinates(target_pos, rx, ry, rz, target_rotation, num_points, rng, jitter = True)
                pos = self.synth_coords[neuron_id]
                
                # Convert to coordinates in the hyper voxel
                voxel_coords = np.round((pos - self.snudda_detect.hyper_voxel_origo)
                                        / self.snudda_detect.voxel_size).astype(int)

                inside_idx = np.where(np.sum(np.logical_and(0 <= voxel_coords, voxel_coords < self.snudda_detect.hyper_voxel_size),axis=1) == 3)[0]
                print(len(inside_idx))

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
                        # self.add_projection(projection_name=projection_name, pre_neuron_type=pre_neuron_type,
                        #                     projection_file=con_config["projection_config_file"])

    def angle_with_reference(self, v, ref):
        cos_theta = np.dot(v, ref) / (np.linalg.norm(v) * np.linalg.norm(ref))
        # Clip to avoid numerical errors
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        return np.arccos(cos_theta)
    
    def load_gmm_from_json(self, json_path):
        
        with open(json_path, 'r') as f:
            gmm_params = json.load(f)
        
        gmm = BayesianGaussianMixture(covariance_type = 'full', n_components = len(gmm_params['weights']))    
        gmm.weights_ = np.array(gmm_params['weights'])
        gmm.means_ = np.array(gmm_params['means'])
        gmm.covariances_ = np.array(gmm_params['covariances'])
        return gmm
    
    
    def partial_shuffle(self, array, fraction=0.1):
        array = array.copy()
        num_swaps = int(len(array) * fraction)
        for _ in range(num_swaps):
            i, j = np.random.randint(0, len(array), size=2)
            array[i], array[j] = array[j], array[i]
        return array

    def bounding_box_check(self, point, min_point, max_point):
        x, y, z = point
        xmin, ymin, zmin = min_point
        xmax, ymax, zmax = max_point
        return xmin <= x <= xmax and ymin <= y <= ymax and zmin <= z <= zmax


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
                                  xi=pre_positions, method="nearest")

        axon_dist = np.linalg.norm(target_centres - pre_positions, axis=1)

        if "rotation" in proj_info:
            target_rotation = griddata(points=proj_info["dest"], values=proj_info["rotation"],
                                       xi=target_centres, method="nearest")
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
