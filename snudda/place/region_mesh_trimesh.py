import numexpr
import numpy as np
from scipy.spatial import cKDTree
import trimesh
from numba import jit

# Teseting trimesh as alterantive to open3d, but it uses more memory and is slower...

class RegionMeshRedux:

    def __init__(self, mesh_path, verbose=False):

        self.mesh_path = mesh_path
        self.mesh = trimesh.load(mesh_path)
        self.verbose = verbose

        # Convert from micrometers to meters to get SI units
        scale_factor = 1e-6
        self.mesh.apply_scale(scale_factor)

        # Trimesh equivalent of removing non-manifold edges
        # Remove degenerate faces (updated API)
        self.mesh.update_faces(self.mesh.nondegenerate_faces())
        # Remove duplicate faces (updated API)
        self.mesh.update_faces(self.mesh.unique_faces())

        self.min_coord = self.mesh.bounds[0]  # minimum bounds
        self.max_coord = self.mesh.bounds[1]  # maximum bounds

        # Get volume - trimesh can compute volume directly
        self.volume = abs(self.mesh.volume)  # abs in case mesh has flipped normals

    def check_inside(self, points):
        """ Check if points are inside, returns bool array."""

        # Process points in batches to reduce memory usage
        batch_size = 10000  # Adjust based on available memory
        n_points = len(points)
        is_inside = np.zeros(n_points, dtype=bool)

        if self.verbose:
            print(f"Processing {n_points} points in batches of {batch_size}")

        for i in range(0, n_points, batch_size):
            end_idx = min(i + batch_size, n_points)
            batch_points = points[i:end_idx]

            if self.verbose and i % (batch_size * 10) == 0:
                print(f"Processing batch {i // batch_size + 1}/{(n_points - 1) // batch_size + 1}")

            # Use trimesh contains method on batch
            is_inside[i:end_idx] = self.mesh.contains(batch_points)

        return is_inside

    def distance_to_border(self, points):
        """ Positive values are distance to mesh (outside), and negative (inside)"""

        # Process in batches for memory efficiency
        batch_size = 50000
        n_points = len(points)
        signed_distance = np.zeros(n_points)

        for i in range(0, n_points, batch_size):
            end_idx = min(i + batch_size, n_points)
            batch_points = points[i:end_idx]

            # Use trimesh proximity for signed distance on batch
            signed_distance[i:end_idx] = trimesh.proximity.signed_distance(self.mesh, batch_points)

        return signed_distance

    def plot(self, line_set=None, neurons=None, show_axis=False, show_faces=True):

        # Press w to see wireframe...

        if show_faces:
            plot_list = [self.mesh]
        else:
            # Just plot wireframe
            plot_list = [self.mesh.as_open3d.wireframe()]

        if line_set:
            plot_list += line_set

        if neurons:
            plot_list += self.get_line_set(neurons)

        if show_axis:
            # The x, y, z axis will be rendered as red, green, and blue arrows respectively.
            axis_length = 1000e-6
            # Create simple coordinate frame using trimesh
            origin = [0, 0, 0]
            # For now, simplified - you can enhance this later
            pass

        # Use trimesh's built-in show
        scene = trimesh.Scene(plot_list)
        scene.show()

    def get_line_set(self, neurons):

        if type(neurons) != list:
            return self.get_line_set(neurons=[neurons])

        line_sets = []

        for neuron in neurons:

            morph_data = neuron.morphology_data["neuron"]

            lines = []
            for section in morph_data.section_iterator():
                if len(section.point_idx) > 1:
                    for start_point, end_point in zip(section.point_idx[:-1], section.point_idx[1:]):
                        lines.append([start_point, end_point])

            if len(lines) > 0:
                # Create line data for trimesh
                points = morph_data.geometry[:, :3]
                # Convert to trimesh Path3D format
                path_data = {
                    'vertices': points,
                    'entities': [trimesh.path.entities.Line(lines)]
                }
                line_sets.append(path_data)

        return line_sets


class NeuronPlacer:

    def __init__(self, mesh_path: str, d_min: float, random_seed=None, rng=None,
                 n_putative_points=None, putative_density=None, verbose=False):

        """ Args:
            mesh_path (str): Path to wavefront obj file
            d_min (float): Minimum distance between neurons
            random_seed (int): Random seed
            rng: Numpy rng object, either rng or random_seed is given
            n_putative_points (int): Number of putative positions to place within volume (before d_min filtering)"""

        self.verbose = verbose
        self.region_mesh = RegionMeshRedux(mesh_path=mesh_path, verbose=verbose)
        self.d_min = d_min
        self.density_functions = dict()

        if rng:
            self.rng = rng

            if random_seed:
                raise ValueError("If rng is set, then seed should not be set.")
        else:
            self.rng = np.random.default_rng(random_seed)

        # We generate a cube of points, obs that we pad it with d_min on all sides to avoid
        # edge effects (which would give higher densities at borders)
        self.cube_side = self.region_mesh.max_coord - self.region_mesh.min_coord + self.d_min * 2
        self.cube_offset = self.region_mesh.min_coord - self.d_min

        if n_putative_points is None:
            # The volume of the cube multiplied by a density estimated by d_min

            if putative_density:
                n_putative_points = int(np.ceil(np.prod(self.cube_side) * putative_density * 1e9))
            else:

                # n_putative_points = min(int(np.ceil(np.prod(self.cube_side) * (1/self.d_min) ** 3)), 1000000)
                n_putative_points = min(int(np.ceil(np.prod(self.cube_side) * 300e3 * 1e9)), 1000000)

                print(f"No n_putative_points and putative_density, setting {n_putative_points = }"
                      f"\n(this must be larger than the number of neurons you want to place)")
        else:
            # We need to compenate n_putative_points for fact that we sample points outside volume also
            n_putative_points *= np.prod(self.cube_side) / self.region_mesh.volume
            n_putative_points = int(np.ceil(n_putative_points))

        print(f"Generating {n_putative_points} points for {mesh_path}")

        putative_points = self.get_point_cloud(n=n_putative_points)
        putative_points = self.remove_close_neurons(putative_points)
        putative_points = self.remove_outside(putative_points)

        if putative_points.shape[0] < 0.05 * n_putative_points:
            print(f"Managed to create {putative_points.shape[0]} putative points within the volume.\n"
                  f"  WARNING --> is the volume too small? You can create new cube mesh using create_cube_mesh.py\n\n"
                  f"Example how to use create_cube_mesh.py:\n"
                  f"from snudda.place.create_cube_mesh import create_cube_mesh\n"
                  f"create_cube_mesh(file_name='your_cube_mesh_name.obj', \n"
                  f"                 centre_point=(0,0,0), side_len=300e-6,\n"
                  f"                 description='Adjust side_len to get correct neuron density')\n"
                  )

        self.putative_points = putative_points
        self.allocated_points = np.zeros(shape=(putative_points.shape[0],), dtype=bool)

    def define_density(self, neuron_type, density_function):
        if neuron_type in self.density_functions:
            print(f"Warning, overwriting {neuron_type} density with {density_function}")

        self.density_functions[neuron_type] = density_function

    def place_neurons(self, num_neurons, neuron_type=None):

        if neuron_type is None or neuron_type not in self.density_functions:
            density_function = None
        else:
            density_function = self.density_functions[neuron_type]

        return self.get_neuron_positions(n_positions=num_neurons, neuron_density=density_function)

    def plot_putative_points(self):

        self.plot_points(points=self.putative_points)

    def plot_placed_points(self):

        self.plot_points(points=self.putative_points[self.allocated_points, :])

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
            if sorted_counts[first_pair] != 1:
                first_pair = len(sorted_counts) - 1  # Basically use remove_fraction_idx

            remove_fraction_idx = int(np.ceil(remove_fraction * len(sorted_offenders)))
            remove_idx = sorted_offenders[:min(first_pair, remove_fraction_idx)]

        points = np.delete(points, remove_idx, axis=0)

        if self.verbose:
            print(f"n_points = {points.shape[0]}, previous close_pairs = {len(close_pairs)}")

        return points, False

    def remove_outside(self, points):

        if self.verbose:
            print(f"Filtering {points.shape[0]} points..")

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
        if len(free_positions) > 1:
            closest_distance, _ = cKDTree(data=free_positions).query(x=free_positions, k=2)
            # Volume is proportional to distance**3, so scale probabilities to pick position by that
            free_volume = np.power(np.mean(closest_distance[:, 1:2], axis=1), 3)
            P_neuron = free_volume
        else:
            # Special case, when there is only one
            P_neuron = [1]
            assert neuron_density is None, "You can not specify neuron_density if there is only one free position."

        x, y, z = free_positions.T

        if neuron_density:
            # TODO: Temp disabled volume... still does not seem to work
            if type(neuron_density) == str:
                P_neuron = np.multiply(numexpr.evaluate(neuron_density), free_volume)
            else:
                P_neuron = np.multiply(neuron_density(x=x, y=y, z=z), free_volume)
            # P_neuron = numexpr.evaluate(neuron_density)

        P_neuron /= np.sum(P_neuron)

        try:
            idx = self.rng.choice(len(free_positions), n_positions, p=P_neuron, replace=False)
        except ValueError as ve:
            print("Error: Increase n_putative_points or putative_density, too few putative points set.")
            print(f"      Mesh: {self.region_mesh.mesh_path}")

            import pdb
            pdb.set_trace()

            raise ve

        neuron_positions = free_positions[idx, :]
        used_idx = np.where(~self.allocated_points)[0][idx]
        self.allocated_points[used_idx] = True

        return neuron_positions


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
    mesh_path = "../data/mesh/Striatum-d-right.obj"

    # nep = NeuronPlacer(mesh_path=mesh_path, d_min=10e-6, n_putative_points=10000000)
    nep = NeuronPlacer(mesh_path=mesh_path, d_min=10e-6,
                       n_putative_points=None, putative_density=100e3,
                       verbose=True)
    # nep.plot_putative_points()

    points_flat = nep.get_neuron_positions(5000)
    nep.plot_points(points_flat, colour="black")

    points = nep.get_neuron_positions(200000, neuron_density="exp((y-0.0025)*2000)")
    nep.plot_points(points, colour="red")

    points_flat2 = nep.get_neuron_positions(5000)
    nep.plot_points(points_flat2, colour="blue")

    import matplotlib.pyplot as plt

    plt.figure()
    plt.hist(points_flat[:, 1], color="black")
    plt.hist(points[:, 1], color="red", alpha=0.5)
    plt.show()

    plt.figure()
    plt.hist(points_flat[:, 1], color="black")
    plt.hist(points_flat2[:, 1], color="blue", alpha=0.5)
    plt.show()

    import pdb

    pdb.set_trace()