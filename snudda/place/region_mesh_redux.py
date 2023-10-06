
import numpy as np
from scipy.spatial import cKDTree
from numba import jit


class RegionMeshRedux:

    def __init__(self, mesh_path):

        self.mesh_path = mesh_path

        # Set by load_mesh
        self.mesh_vertices = None
        self.mesh_faces = None
        self.mesh_normals = None
        self.min_coord = None
        self.max_coord = None
        self.point_outside = None

        # Set by pre_compute
        self.mesh_u = None
        self.mesh_v = None
        self.mesh_v0 = None
        self.mesh_uv = None
        self.mesh_vv = None
        self.mesh_uu = None
        self.mesh_denom = None
        self.mesh_nrm = None

        if self.mesh_path is not None:
            self.load_mesh(file_path=self.mesh_path)
            self.pre_compute()

    def load_mesh(self, file_path):

        """ Reads wavefront obj files. Updates mesh_vertices, mesh_faces, mesh_normals,
            min_coord, max_coord, point_outside """

        vertices = []
        normals = []
        faces = []

        with open(file_path, 'r') as obj_file:
            for line in obj_file:
                parts = line.strip().split()
                if not parts:
                    continue  # Skip empty lines

                if parts[0][0] == 'v':
                    if len(parts[0]) > 1:
                        if parts[0][1] == 'n':
                            # Parse normal vectors
                            nx, ny, nz = map(float, parts[1:4])
                            normals.append((nx, ny, nz))
                        else:
                            # Skip "vt" lines and others
                            continue
                    else:
                        # Parse vertex coordinates in micrometers and convert to meters
                        x, y, z = [float(coord) * 1e-6 for coord in parts[1:4]]
                        vertices.append((x, y, z))
                elif parts[0] == 'f':
                    # Parse face definitions with vertex and normal indices
                    face_indices = []
                    for vertex in parts[1:]:
                        try:
                            # Only take first value of each triplet 1/?/? 2/?/? 3/?/?
                            indices = int(vertex.split('/')[0]) - 1
                        except:
                            import traceback
                            print(traceback.format_exc())
                            import pdb
                            pdb.set_trace()
                        # indices = [int(idx) - 1 if idx else 0 for idx in vertex.split('/')]
                        face_indices.append(indices)
                    faces.append(face_indices)

        self.mesh_vertices = np.array(vertices)
        self.mesh_faces = np.array(faces)
        self.mesh_normals = np.array(normals)

        self.min_coord = np.min(self.mesh_vertices, axis=0)
        self.max_coord = np.max(self.mesh_vertices, axis=0)

        # Used by ray casting when checking if another point is interior
        self.point_outside = self.max_coord + np.array([1e-1, 1e-2, 1e-3])

    def pre_compute(self):

        """ Helper function, precomputes values for raytracing. """

        i0 = self.mesh_faces[:, 0]
        i1 = self.mesh_faces[:, 1]
        i2 = self.mesh_faces[:, 2]

        self.mesh_u = self.mesh_vertices[i1, :] - self.mesh_vertices[i0, :]
        self.mesh_v = self.mesh_vertices[i2, :] - self.mesh_vertices[i0, :]

        self.mesh_v0 = self.mesh_vertices[i0, :]

        self.mesh_uv = np.sum(np.multiply(self.mesh_u, self.mesh_v), axis=1)
        self.mesh_vv = np.sum(np.multiply(self.mesh_v, self.mesh_v), axis=1)
        self.mesh_uu = np.sum(np.multiply(self.mesh_u, self.mesh_u), axis=1)

        self.mesh_denom = np.multiply(self.mesh_uv, self.mesh_uv) - np.multiply(self.mesh_uu, self.mesh_vv)

        # Normal of triangle
        self.mesh_nrm = np.cross(self.mesh_u, self.mesh_v)

        # We need to normalise it
        nl = np.repeat(np.reshape(self.mesh_uv, [self.mesh_uv.shape[0], 1]), 3, axis=1)
        self.mesh_nrm = np.divide(self.mesh_nrm, nl)

    def ray_casting(self, point):

        """ Ray-casting, to determine if a point is inside or outside of mesh. """

        return RegionMeshRedux._ray_casting_helper(point=point,
                                                   self_mesh_faces=self.mesh_faces,
                                                   self_mesh_nrm=self.mesh_nrm,
                                                   self_mesh_v0=self.mesh_v0,
                                                   self_point_out=self.point_outside,
                                                   self_mesh_denom=self.mesh_denom,
                                                   self_mesh_uv=self.mesh_uv,
                                                   self_mesh_uu=self.mesh_uu,
                                                   self_mesh_vv=self.mesh_vv,
                                                   self_mesh_u=self.mesh_u,
                                                   self_mesh_v=self.mesh_v)

    @staticmethod
    @jit(nopython=True)
    def _ray_casting_helper(point,
                            self_mesh_faces, self_mesh_nrm, self_point_out,
                            self_mesh_v0, self_mesh_denom,
                            self_mesh_uv, self_mesh_vv, self_mesh_uu, self_mesh_u, self_mesh_v):

        """
        Helper function for ray-casting, to determine if a point is inside the 3D-mesh.
        It draws a line from the point given, and a second point defined as outside.
        If that line intersects the surface of the 3D mesh an odd number of times, then the first point is inside.

        Uses values pre-computed by pre_compute function.
        """

        # print(f"Processing {point}")

        n_tri = self_mesh_faces.shape[0]

        p = self_point_out - point
        # rn = nominator, rd = denominator
        rn = np.sum(np.multiply(self_mesh_nrm, self_mesh_v0 - point), axis=1)
        rd = np.dot(self_mesh_nrm, p)

        # If rd == 0 and rn != 0 --> r = -1, parallel to plane, but outside, mark -1 to avoid counting
        # If rd == 0 and rn == 0 --> r = 0, parallel and lies in plane
        idx0 = (rd == 0)
        idx1 = np.logical_and(idx0, rn != 0)

        rn[idx1] = -1
        rd[idx0] = 1
        r = np.divide(rn, rd)

        w = point + r.reshape(len(r), 1) * p.reshape(1, 3) - self_mesh_v0
        n_points = len(r)

        s = np.divide(np.multiply(self_mesh_uv, np.sum(np.multiply(w, self_mesh_v), axis=1).reshape(n_points, ))
                      - np.multiply(self_mesh_vv, np.sum(np.multiply(w, self_mesh_u), axis=1).reshape(n_points, )),
                      self_mesh_denom)

        t = np.divide(np.multiply(self_mesh_uv, np.sum(np.multiply(w, self_mesh_u), axis=1).reshape(n_points, ))
                      - np.multiply(self_mesh_uu, np.sum(np.multiply(w, self_mesh_v), axis=1).reshape(n_points, )),
                      self_mesh_denom)

        intersect_count = np.sum((0 <= r) * (r <= 1) * (0 <= s) * (s <= 1) * (0 <= t) * (s + t <= 1))

        return np.mod(intersect_count, 2) == 1

    def point_inside(self, points):
        # Let's see how fast it is if we do ray casting directly

        inside_flag = np.zeros(shape=(points.shape[0],), dtype=bool)

        for ctr, p in enumerate(points):
            inside_flag[ctr] = self.ray_casting(p)

        return inside_flag

    def distance_to_border(self, point):

        pass


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
        keep_flag = self.region_mesh.point_inside(points=points)

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
    points = nep.get_point_cloud(n=100000)
    new_points = nep.remove_close_neurons(points)
    new_inside_points = nep.remove_outside(new_points)

    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    x,y,z = new_inside_points.T
    ax.scatter(x, y, z, marker='.')
    plt.show()

    import pdb
    pdb.set_trace()