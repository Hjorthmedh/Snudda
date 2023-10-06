import numexpr
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
        self.mesh.scale(scale_factor, center=(0, 0, 0))

        self.mesh.remove_non_manifold_edges()

        self.min_coord = self.mesh.get_min_bound()
        self.max_coord = self.mesh.get_max_bound()

        self.scene = o3d.t.geometry.RaycastingScene()
        legacy_mesh = o3d.t.geometry.TriangleMesh.from_legacy(self.mesh)

        # SHALL WE PICk NUMBER OF PUTATIVE POINTS BASED ON VOLUME???

        self.scene.add_triangles(legacy_mesh)

    def check_inside(self, points):
        """ Check if points are inside, returns bool array."""

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

    def __init__(self, mesh_path: str, d_min: float, seed=None, rng=None, n_putative_points=1000000):

        self.region_mesh = RegionMeshRedux(mesh_path=mesh_path)
        self.d_min = d_min

        if rng:
            if seed:
                raise ValueError("If rng is set, then seed should not be set.")
            self.rng = rng
        else:
            self.rng = np.random.default_rng(seed)

        # We generate a cube of points, obs that we pad it with d_min on all sides to avoid
        # edge effects (which would give higher densities at borders)
        self.cube_side = self.region_mesh.max_coord-self.region_mesh.min_coord + self.d_min*2
        self.cube_offset = self.region_mesh.min_coord - self.d_min

        putative_points = self.get_point_cloud(n=n_putative_points)
        putative_points = self.remove_close_neurons(putative_points)
        putative_points = self.remove_outside(putative_points)

        self.putative_points = putative_points
        self.allocated_points = np.zeros(shape=(putative_points.shape[0],), dtype=bool)

    def plot_putative_points(self):

        self.plot_points(points=self.putative_points)

    def plot_points(self, points, colour=None):
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        x, y, z = points.T
        ax.scatter(x, y, z, marker='.', color=colour)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.ion()
        plt.show()

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
            # First remove the worst offenders, with the most neighbours
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

    def get_neuron_positions(self, n_positions, neuron_density=None):

        """ neuron_density either None (even), or a str representing a function f(x,y,z)
            where x,y,z are the coordinates in meters """

        # We have the putative_points, pick positions from them, based on neuron density
        # then update the allocated points

        # 1. Calculate the distance to the closest free (non-allocated) neighbour for all points
        # This is so that if we allocated neurons, we can correct for that, and get flat gradients

        free_positions = self.putative_points[~self.allocated_points]

        # k=2, since we don't want distance to point itself, but closest neighbour
        # closest_distance, _ = cKDTree(data=free_positions).query(x=free_positions, k=2)
        closest_distance, _ = cKDTree(data=free_positions).query(x=free_positions, k=2)

        # Volume is proportional to distance**3, so scale probabilities to pick position by that
        free_volume = np.power(np.mean(closest_distance[:, 1:2], axis=1), 3)
        x, y, z = free_positions.T

        if neuron_density:
            # TODO: Temp disabled volume... still does not seem to work
            P_neuron = np.multiply(numexpr.evaluate(neuron_density), free_volume)
            # P_neuron = numexpr.evaluate(neuron_density)
        else:
            P_neuron = free_volume

        P_neuron /= np.sum(P_neuron)

        idx = self.rng.choice(len(free_positions), n_positions, p=P_neuron, replace=False)

        neuron_positions = free_positions[idx, :]
        used_idx = np.where(~self.allocated_points)[0][idx]
        self.allocated_points[used_idx] = True

        return neuron_positions



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

    mesh_path="/home/hjorth/HBP/Snudda/snudda/data/mesh/Striatum-dorsal-right-hemisphere.obj"

    nep = NeuronPlacer(mesh_path=mesh_path, d_min=10e-6, n_putative_points=10000000)
    # nep.plot_putative_points()

    points_flat = nep.get_neuron_positions(5000)
    nep.plot_points(points_flat, colour="black")

    points = nep.get_neuron_positions(200000, neuron_density="exp((y-0.0025)*2000)")
    nep.plot_points(points, colour="red")

    points_flat2 = nep.get_neuron_positions(5000)
    nep.plot_points(points_flat, colour="blue")

    import matplotlib.pyplot as plt
    plt.figure()
    plt.hist(points_flat[:,1], color="black")
    plt.hist(points[:,1], color="red", alpha=0.5)
    plt.show()

    plt.figure()
    plt.hist(points_flat[:,1], color="black")
    plt.hist(points_flat2[:,1], color="blue", alpha=0.5)
    plt.show()


    import pdb
    pdb.set_trace()