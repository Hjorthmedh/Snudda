
import numpy as np
from scipy.spatial import cKDTree


class RegionMeshRedux:

    def __init__(self):

        pass

    def load_mesh(self):

        pass

    def point_inside(self, point):

        pass

    def distance_to_border(self, point):

        pass


class NeuronPlacer:

    # Improvement suggestion. Randomize only points within a mesh-voxel, then randomise
    # within which of the mesh voxels the points should be placed.

    def __init__(self, region_mesh: RegionMeshRedux, d_min: float, seed=None):

        self.region_mesh = region_mesh
        self.d_min = d_min

        self.rng = np.random.default_rng(seed)

        # Temp values, should be read from the 3d mesh + padding
        self.cube_side = 1000e-6
        self.cube_offset = np.array([1, 2, 5000])

        # Use kdtree to avoid d_min overlaps

        pass

    def get_point_cloud(self, n):
        points = self.rng.uniform(size=(n, 3), low=0, high=self.cube_side) + self.cube_offset

        return points

    def remove_close_neurons(self, points):

        done = False

        while not done:
            points, done = self._remove_close_neurons_helper(points)

        return points

    def _remove_close_neurons_helper(self, points):

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
            ten_percent = int(np.ceil(0.05*len(sorted_offenders)))

            remove_idx = sorted_offenders[:min(first_pair, ten_percent)]

        points = np.delete(points, remove_idx, axis=0)

        print(f"n_points = {points.shape[0]}, prev_close_pairs = {len(close_pairs)}")

        return points, False


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

    nep = NeuronPlacer(region_mesh=None, d_min=10e-6)
    points = nep.get_point_cloud(n=1000000)
    new_points = nep.remove_close_neurons(points)

    import pdb
    pdb.set_trace()