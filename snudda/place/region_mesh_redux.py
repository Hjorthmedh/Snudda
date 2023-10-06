
import numpy as np
from scipy.spatial import cKDTree
import open3d as o3d
from numba import jit


class RegionMeshRedux:

    def __init__(self, mesh_path):

        self.mesh_path = mesh_path
        self.mesh = o3d.io.read_triangle_mesh(mesh_path)

        # Convert from micrometers to meters to get SI units
        scale_factor = 1e-6
        self.mesh.scale(scale_factor, center=self.mesh.get_center()*scale_factor)

        self.min_coord = self.mesh.get_min_bound()
        self.max_coord = self.mesh.get_max_bound()

        self.scene = o3d.t.geometry.RaycastingScene()
        legacy_mesh = o3d.t.geometry.TriangleMesh.from_legacy(self.mesh)
        self.scene.add_triangles(legacy_mesh)

    def check_inside(self, points):

        # http://www.open3d.org/docs/latest/tutorial/geometry/distance_queries.html
        query_point = o3d.core.Tensor(points, dtype=o3d.core.Dtype.Float32)
        is_inside = self.scene.compute_occupancy(query_point)
        return is_inside.numpy().astype(bool)

    def distance_to_border(self, points):

        """ Positive values are distance to mesh (outside), and negative (inside)"""

        # http://www.open3d.org/docs/latest/tutorial/geometry/distance_queries.html
        query_point = o3d.core.Tensor(points, dtype=o3d.core.Dtype.Float32)
        signed_distance = self.scene.compute_signed_distance(query_point).numpy()

        return signed_distance


class NeuronPlacer:

    # Improvement suggestion. Randomize only points within a mesh-voxel, then randomise
    # within which of the mesh voxels the points should be placed.

    def __init__(self, region_mesh: RegionMeshRedux, d_min: float, seed=None):

        self.region_mesh = region_mesh
        self.d_min = d_min

        self.rng = np.random.default_rng(seed)

        # We generate a cube of points, obs that we pad it with d_min on all sides to avoid
        # edge effects (which would give higher densities at borders)
        self.cube_side = self.region_mesh.max_coord-self.region_mesh.min_coord + self.d_min*2
        self.cube_offset = self.region_mesh.min_coord - self.d_min

    def get_point_cloud(self, n):
        points = self.rng.uniform(size=(n, 3), low=0, high=self.cube_side) + self.cube_offset

        return points

    def remove_close_neurons(self, points):

        done = False

        while not done:
            points, done = self._remove_close_neurons_helper(points)

        return points

    def _remove_close_neurons_helper(self, points, remove_fraction=0.05):

        close_pairs = cKDTree(data=points).query_pairs(r=self.d_min)

        if len(close_pairs) == 0:
            return points, True

        close_neurons = np.array(list(close_pairs), dtype=int)
        unique, counts = np.unique(close_neurons, return_counts=True)

        if max(counts) == 1:
            # Only single pairs left, be smart about removing them, just pick one of each pair
            remove_idx = close_neurons[np.arange(0, close_neurons.shape[0]),
                                       self.rng.integers(low=0, high=2, size=close_neurons.shape[0])]

        else:

            # First remove worst offenders, with the most neighbours
            # (but none with count 1, we do those separately
            sort_idx = np.argsort(-counts)
            sorted_offenders = unique[sort_idx]
            sorted_counts = counts[sort_idx]

            first_pair = np.argmax(sorted_counts == 1)
            remove_fraction_idx = int(np.ceil(remove_fraction*len(sorted_offenders)))

            remove_idx = sorted_offenders[:min(first_pair, remove_fraction_idx)]

        points = np.delete(points, remove_idx, axis=0)

        print(f"n_points = {points.shape[0]}, previous close_pairs = {len(close_pairs)}")

        return points, False

    def remove_outside(self, points):

        print(f"Filtering {points.shape[0]} points..")
        # keep_flag = self.region_mesh.point_inside(points=points)
        keep_flag = self.region_mesh.check_inside(points=points)

        print(f"Filtering, keeping inside points: {np.sum(keep_flag)} / {len(keep_flag)}")

        return points[keep_flag, :]

    def set_neuron_density(self, neuron_type, neuron_density):

        pass

    def place_neurons(self, neuron_type, number_of_neurons):

        pass


class NeuronBender:

    def __init__(self):
        pass

# Goal for today:
# - Function that returns N positions that are not within d_min of each other
# - Within a given volume
# - Also not colliding with already previously placed points.

# Idea:
# We want to have N positions.
# Randomize M > N positions
# Build a cKDTree
# Find too close neighbours query_pairs(d_min)
# Remove the worst offenders -- how?



if __name__ == "__main__":

    rmr = RegionMeshRedux(mesh_path="/home/hjorth/HBP/Snudda/snudda/data/mesh/Striatum-dorsal-right-hemisphere.obj")

    nep = NeuronPlacer(region_mesh=rmr, d_min=10e-6)
    points = nep.get_point_cloud(n=1000000)
    new_points = nep.remove_close_neurons(points)
    new_inside_points = nep.remove_outside(new_points)

    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    x,y,z = new_inside_points.T
    ax.scatter(x, y, z, marker='.')
    plt.show()

    # Lets calculate the distance to border
    distance = nep.region_mesh.distance_to_border(new_inside_points)

    import pdb
    pdb.set_trace()