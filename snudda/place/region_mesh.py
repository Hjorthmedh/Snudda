import sys

import numpy as np
import scipy
from scipy import ndimage
from numba import jit
import re
import os
import pickle
import timeit


class RegionMesh(object):

    """ Handles neuron placement within a 3D mesh. """

    def __init__(self, filename, d_view=None, role="master",
                 use_cache=True, pickle_version=-1, raytrace_borders=True,
                 d_min=15e-6, bin_width=1e-4,
                 logfile_name=None, log_file=None,
                 random_seed=112, verbose=False):

        """
        Constructor.

        Args:
            filename (str) : Path to wavefront mesh file
            d_view : ipyparallel direct view object
            role (str) : Role, ie. "master" or "worker"
            use_cache (bool) : The meshes voxel representation is cached, use the cache?
            pickle_version (int) : Which version of pickle to use? (-1 latest)
            raytrace_borders (bool) : Raytrace positions in border regions? Slower, but more accurate.
            d_min (float) : Closest distance between soma
            bin_width (float) : Mesh size
            logfile_name (str) : Path to logfile
            log_file : File pointer to log file
            random_seed (int) : Random seed value
            verbose (bool) : Flag to be verbose?

        """

        self.d_view = d_view

        self.role = role
        self.workers_initialised = False

        self.verbose = verbose
        self.close_log_file = False

        if log_file is not None:
            self.logfile = log_file
            self.logfile_name = log_file.name
        elif logfile_name is not None and len(logfile_name) > 0:
            self.logfile = open(logfile_name, 'wt')
            self.logfile_name = logfile_name
            self.close_log_file = True
        else:
            self.logfile = None
            self.logfile_name = None

        self.bin_width = bin_width
        self.d_min = d_min
        self.padding = max(self.bin_width, d_min)
        self.num_bins = None

        self.density_function = dict()
        self.density_voxel_sum = dict()
        self.density_total_sum = dict()
        self.placed_voxel = dict()
        self.placed_total = dict()
        self.density_voxel_n_sample = dict()
        self.density_total_n_sample = dict()

        self.last_neuron_type_added = None
        self.density_lookup = None

        self.random_seed = random_seed
        self.random_generator = np.random.default_rng(self.random_seed)

        # Set by load_mesh
        self.mesh_vec = None
        self.mesh_faces = None
        self.mesh_norm = None
        self.min_coord = None
        self.max_coord = None
        self.voxel_mask_inner = None
        self.voxel_mask_border = None
        self.point_out = None

        # Set by pre_compute
        self.mesh_u = None
        self.mesh_v = None
        self.mesh_v0 = None
        self.mesh_uv = None
        self.mesh_vv = None
        self.mesh_uu = None
        self.mesh_denom = None
        self.mesh_nrm = None

        # Used by setup_place_neurons
        self.max_rand = 10000
        self.max_neurons = 3000000
        self.max_reject = 100e6

        # Set by setup_place_neurons
        self.rand_ctr = None
        self.random_pool = None
        self.neuron_coords = None
        self.neuron_ctr = None
        self.padding_coords = None
        self.padding_ctr = None
        self.neuron_types = None
        self.neuron_type = None
        self.next_neuron_type = 1
        self.reject_ctr = None

        # Used or set by setup_voxel_list
        self.max_neurons_voxel = int(np.ceil(300000*(self.bin_width/1e-3)**3))  # We assume no more than 300k neurons per mm3
        self.voxel_next_neuron = None
        self.voxel_neurons = None

        # Set by update_padding_mask
        self.voxel_mask_padding = None

        # This determines if we ray trace the border voxels, for finer detail
        # or not (activating this is SLOW)
        self.raytrace_borders = raytrace_borders

        if raytrace_borders:
            self.write_log("Ray tracing points in border voxels, this is slow.")
            rt_str = "-RTB"
        else:
            rt_str = ""

        # binWidth 5e-4 (94s) --> 10.8 % border voxels
        # binWidth 2.5e-4 (577) --> 6.2 % border voxels
        # binWidth 1e-4 (8090s) --> 2.7 % border voxels
        # binWidth 0.5e-4 (??? s) --> ?? % border voxels

        self.filename = filename
        self.cache_file = f"{filename}-{int(1e6 * self.bin_width)}{rt_str}-cache.pickle"
        self.pickle_version = pickle_version

        self.debug_flag = False

        cache_loaded = False

        if use_cache and self.cache_exist():
            try:
                self.load_cache()
                cache_loaded = True
            except:
                self.write_log("Failed to load cache.")
                cache_loaded = False

        if not cache_loaded:
            self.load_mesh(filename=filename)

            t_a = timeit.default_timer()

            self.pre_compute()
            self.setup_voxel_filter()

            t_b = timeit.default_timer()

            self.write_log(f"Calculation time: {t_b - t_a} s")

        self.setup_voxel_list()

        self.setup_place_neurons(d_min=d_min)

        self.write_log(f"Inner voxel bin volume: {np.round(self.inner_voxel_volume() * 1e9, 1)} mm³")
        # !!! RUN THIS
        # self.verifyInside(100)

    ############################################################################

    def mark_borders(self):

        """ Mark border voxels in 3D mesh. """

        self.write_log("Marking borders")

        for row in self.mesh_vec:
            ic = np.floor((row - self.min_coord) / self.bin_width)
            self.voxel_mask_border[int(ic[0]), int(ic[1]), int(ic[2])] = 1

        max_num = 0

        for row in self.mesh_faces:
            coord1 = self.mesh_vec[row[0], :]
            coord2 = self.mesh_vec[row[1], :]
            coord3 = self.mesh_vec[row[2], :]

            vx = coord2 - coord1
            vy = coord3 - coord1

            dx = np.linalg.norm(coord2 - coord1)
            dy = np.linalg.norm(coord3 - coord1)

            nx = int(2 * np.ceil(dx / self.bin_width)) + 1
            ny = int(2 * np.ceil(dy / self.bin_width)) + 1

            max_num = max(max_num, max(nx, ny))

            for x_step in np.linspace(0, 1, nx):
                for y_step in np.linspace(0, 1 - x_step, ny):
                    x_point = coord1 + vx * x_step + vy * y_step
                    ic = np.floor((x_point - self.min_coord) / self.bin_width)
                    self.voxel_mask_border[int(ic[0]), int(ic[1]), int(ic[2])] = 1

        self.write_log(f"maxN = {max_num}")

    def __del__(self):
        if self.close_log_file:
            self.logfile.close()

    ############################################################################

    def setup_parallel(self, d_view=None):

        """ Setup workers for parallel execution. """

        if d_view is None:

            if self.d_view is None:
                self.write_log("Running in serial")
                return

            d_view = self.d_view

        self.write_log("Setting up parallel")

        assert self.role == "master", \
            "setupParallel should only be called by master"

        if self.workers_initialised:
            self.write_log("setupParallel: workers already initialised")
            return

        with d_view.sync_imports():
            from snudda.place.region_mesh import RegionMesh

        try:
            self.write_log("Setting up RegionMesh on workers")
            d_view.push({"filename": self.filename}, block=True)

            if self.logfile_name is not None:
                engine_log_file = [self.logfile_name + "-"
                                   + str(x) for x in range(0, len(d_view))]
            else:
                engine_log_file = [[] for x in range(0, len(d_view))]

            d_view.scatter('logfile_name', engine_log_file, block=True)

            cmd_str = "sm = RegionMesh(filename=filename,role='worker',logfile_name=logfile_name[0])"

            d_view.execute(cmd_str, block=True)
            self.write_log("Worker RegionMesh setup done.")

        except Exception as e:
            import uuid
            import traceback
            tstr = traceback.format_exc()
            tmp = open("save/tmp-striatum-mesh-log-file-" + str(uuid.uuid4()), 'w')
            tmp.write("Exception: " + str(e))
            tmp.write("Trace:" + tstr)
            tmp.close()
            self.write_log(tstr, is_error=True)
            import pdb
            pdb.set_trace()

        self.workers_initialised = True

    ############################################################################

    def cache_exist(self):

        """ Check if cache for 3D mesh exists. Returns True or False. """

        cache_flag = False

        if os.path.isfile(self.cache_file):

            obj_time = os.path.getmtime(self.filename)
            cache_time = os.path.getmtime(self.cache_file)

            if cache_time > obj_time:
                self.write_log(f"Found cache file {self.cache_file}")
                cache_flag = True
            else:
                self.write_log(f"Found old cache file ({self.cache_file}), ignoring.")
        else:
            self.write_log(f"No mesh cache file found ({self.cache_file})")

        return cache_flag

    ############################################################################

    def load_cache(self):

        """ Loading 3D mesh cache. """

        self.write_log(f"Loading cache file: {self.cache_file}")

        with open(self.cache_file, 'rb') as f:
            data = pickle.load(f)

        self.mesh_vec = data["meshVec"]
        self.mesh_faces = data["meshFaces"]
        self.mesh_norm = data["meshNorm"]
        self.min_coord = data["minCoord"]
        self.max_coord = data["maxCoord"]
        self.voxel_mask_inner = data["voxelMaskInner"]
        self.voxel_mask_border = data["voxelMaskBorder"]

        self.point_out = self.max_coord + np.array([1e-1, 1e-2, 1e-3])
        # self.point_out = self.max_coord + 1e-2

        assert self.bin_width == data["binWidth"], \
            f"Mismatch binWidth: {self.bin_width} vs {data['binWidth']}"
        assert self.padding == data["padding"], \
            f"Mismatch padding: {self.padding} vs {data['padding']}"

        assert self.raytrace_borders == data["raytraceBorders"], \
            f"Mismatch raytraceBorders: {self.raytrace_borders} vs {data['raytraceBorders']}"

        self.num_bins = self.voxel_mask_inner.shape
        self.pre_compute()

    ############################################################################

    def save_cache(self):

        """ Save 3D mesh cache. """

        if self.role != "master":
            return

        data = dict([])
        data["meshVec"] = self.mesh_vec
        data["meshFaces"] = self.mesh_faces
        data["meshNorm"] = self.mesh_norm
        data["minCoord"] = self.min_coord
        data["maxCoord"] = self.max_coord
        data["voxelMaskInner"] = self.voxel_mask_inner
        data["voxelMaskBorder"] = self.voxel_mask_border
        data["padding"] = self.padding
        data["binWidth"] = self.bin_width
        data["raytraceBorders"] = self.raytrace_borders

        self.write_log(f"Saving mesh cache file {self.cache_file}")
        with open(self.cache_file, 'wb') as f:
            pickle.dump(data, f, self.pickle_version)

    ############################################################################

    def load_mesh(self, filename):

        """ Load 3D mesh. """

        self.filename = filename

        all_vec = []
        all_faces = []
        all_norm = []

        # https://stackoverflow.com/questions/4703390/how-to-extract-a-floating-number-from-a-string
        numeric_const_pattern = r"""
    [-+]? # optional sign
    (?:
    (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
    |
    (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
    )
    # followed by optional exponent part if desired
    (?: [Ee] [+-]? \d+ ) ?
    """
        rx = re.compile(numeric_const_pattern, re.VERBOSE)

        with open(filename, 'rt') as f:
            for row in f:
                if row[0:2] == 'v ':
                    digits = rx.findall(row)
                    # Convert to SI units
                    all_vec.append([float(d) * 1e-6 for d in digits])

                if row[0:2] == 'f ':
                    # Only take first value of each triplet 1/?/? 2/?/? 3/?/?
                    digits = re.findall(r'f\s+(\d+)/\d*/\d*\s+(\d+)//\d*\s+(\d+)//\d*',
                                        row)
                    # Subtract one, to get python indexing
                    try:
                        all_faces.append([int(d) - 1 for d in digits[0]])
                    except:
                        self.write_log("Problem with reading digits", is_error=True)
                        self.write_log(f"{row}\nread: {digits}", is_error=True)
                        import pdb
                        pdb.set_trace()

                if row[0:2] == 'vn':
                    digits = rx.findall(row)
                    all_norm.append([float(d) for d in digits])

            self.mesh_vec = np.zeros((len(all_vec), 3))
            self.mesh_faces = np.zeros((len(all_faces), 3), dtype=int)
            self.mesh_norm = np.zeros((len(all_norm), 3))

            for ir, row in enumerate(all_vec):
                self.mesh_vec[ir, :] = row

            for ir, row in enumerate(all_faces):
                self.mesh_faces[ir,] = row

            for ir, row in enumerate(all_norm):
                self.mesh_norm[ir, :] = row

            try:
                self.min_coord = np.min(self.mesh_vec, axis=0)
                self.max_coord = np.max(self.mesh_vec, axis=0)
            except:
                self.write_log("Problem calculating min and max coordinates")
                import pdb
                pdb.set_trace()

            # Used by ray casting when checking if another point is interior
            self.point_out = self.max_coord + np.array([1e-1, 1e-2, 1e-3])
            # self.point_out = self.max_coord + 1e-2

    ############################################################################

    def pre_compute(self):

        """ Helper function, precomputes values for raytracing. """

        i0 = self.mesh_faces[:, 0]
        i1 = self.mesh_faces[:, 1]
        i2 = self.mesh_faces[:, 2]

        self.mesh_u = self.mesh_vec[i1, :] - self.mesh_vec[i0, :]
        self.mesh_v = self.mesh_vec[i2, :] - self.mesh_vec[i0, :]

        self.mesh_v0 = self.mesh_vec[i0, :]

        self.mesh_uv = np.sum(np.multiply(self.mesh_u, self.mesh_v), axis=1)
        self.mesh_vv = np.sum(np.multiply(self.mesh_v, self.mesh_v), axis=1)
        self.mesh_uu = np.sum(np.multiply(self.mesh_u, self.mesh_u), axis=1)

        self.mesh_denom = np.multiply(self.mesh_uv, self.mesh_uv) - np.multiply(self.mesh_uu, self.mesh_vv)

        # Normal of triangle
        self.mesh_nrm = np.cross(self.mesh_u, self.mesh_v)

        # We need to normalise it
        nl = np.repeat(np.reshape(self.mesh_uv, [self.mesh_uv.shape[0], 1]), 3, axis=1)
        self.mesh_nrm = np.divide(self.mesh_nrm, nl)

    ############################################################################

    def setup_voxel_filter(self):

        """ Setup voxel filter for 3D mesh. """

        if self.role == "master":
            self.setup_parallel()

        self.min_coord = np.floor((np.min(self.mesh_vec, axis=0)
                                   - self.padding) / self.bin_width) * self.bin_width
        self.max_coord = np.ceil((np.max(self.mesh_vec, axis=0)
                                  + self.padding) / self.bin_width) * self.bin_width

        self.num_bins = np.array(np.ceil((self.max_coord - self.min_coord) / self.bin_width + 1),
                                 dtype=int)

        self.write_log(f"Voxel mask: {self.num_bins[0]} x {self.num_bins[1]} x {self.num_bins[2]}")

        self.voxel_mask_inner = np.zeros(self.num_bins, dtype=bool)
        self.voxel_mask_border = np.zeros(self.num_bins, dtype=bool)

        # All voxels with a mesh point in them are "border voxels"
        # For all remaining points, check if inner or outer
        # If raytraceBorders is false, we dont raytace for border points
        # at run time, this gives a bit jagged edges, but is MUCH faster
        # when placing cells (no ray tracing then)

        if self.raytrace_borders:
            self.mark_borders()

        num_bins_total = self.num_bins[0] * self.num_bins[1] * self.num_bins[2]
        iter_ctr = 0

        # This second part is only run by the master, it calls the workers
        # to perform part of the computation

        if self.role == "master":
            # This should only be done by master

            if self.d_view is None:

                # No workers, do all work ourselves
                # The worker function adds a dimension (so gather works in parallel
                # case), here we just need to reshape results.
                vm_inner = self._voxel_mask_helper(range(0, self.num_bins[0]))
                self.voxel_mask_inner = np.reshape(vm_inner, self.num_bins)

            else:
                # Distribute the work to the workers
                # Randomize order, to spread work load a bit better -- order should not affect computation
                # as computation is deterministic
                all_x = np.random.permutation(np.arange(0, self.num_bins[0]))

                self.d_view.scatter("x_range", all_x, block=True)
                self.write_log("Starting parallel job")
                self.d_view.execute("innerMask = sm._voxel_mask_helper(x_range)", block=True)
                self.write_log("Gathering results")
                inner_mask = self.d_view.gather("innerMask", block=True)

                for m in inner_mask:
                    self.voxel_mask_inner = np.logical_or(self.voxel_mask_inner, m)

        self.write_log(f"Fraction of border voxels: " 
                       f"{np.sum(self.voxel_mask_border)/np.prod(self.voxel_mask_border.shape)}")

        self.save_cache()

        if np.sum(self.voxel_mask_inner) == 0:
            self.write_log(f"!!! Warning no inner voxels in mesh, is your meshBinWidth={self.bin_width} too large?"
                           f"\nThis will prevent neurons from being placed in the volume.",
                           is_error=True)
            self.write_log(f"mesh file: {self.filename}", is_error=True)

    ############################################################################

    def check_inside(self, coords):

        """ Check if coordinates are inside 3D mesh. """

        idx = np.array(np.floor((coords - self.min_coord) / self.bin_width), dtype=int)

        if self.voxel_mask_inner[idx[0], idx[1], idx[2]]:
            # We know it is an inner voxel
            return True
        elif self.voxel_mask_border[idx[0], idx[1], idx[2]]:
            # We are in a border voxel, need to ray cast this
            return self.ray_casting(coords)
        else:
            # We are outside structure
            return False

    ############################################################################

    def _voxel_mask_helper(self, x_range):

        """ Helper function. """

        try:

            # Need the extra dimension at the top for "gather" work
            vm_inner = np.zeros((1, self.num_bins[0], self.num_bins[1], self.num_bins[2]), dtype=bool)

            for ix in x_range:
                self.write_log(f"Processing x = {ix}")

                for iy in range(0, self.num_bins[1]):
                    # print(f"Processing x = {ix}/{self.num_bins[0]}, y = {iy}/{self.num_bins[1]}")

                    for iz in range(0, self.num_bins[2]):

                        if not self.voxel_mask_border[ix, iy, iz]:
                            # Inner or outer point, check centre
                            xyz = np.array([self.min_coord[0] + (ix + 0.5) * self.bin_width,
                                            self.min_coord[1] + (iy + 0.5) * self.bin_width,
                                            self.min_coord[2] + (iz + 0.5) * self.bin_width])

                            vm_inner[0, ix, iy, iz] = self.ray_casting(xyz)

        except Exception as e:
            # Write error to log file to help trace it.
            import traceback
            t_str = traceback.format_exc()
            self.write_log(t_str, is_error=True)

            sys.exit(-1)

        return vm_inner

    ############################################################################

    # Cast a ray, see how many times it intersects the triangles of the structure
    # if an odd number of intersections, then it is an inner point
    #
    # Based on: https://www.erikrotteveel.com/python/three-dimensional-ray-tracing-in-python/
    # which is based on : source: http://geomalgorithms.com/a06-_intersect-2.html

    # TODO: When the line between the interior and exterior point crosses the line between two vertexes this code
    #       might incorrectly say the line is outside by considering it crosses both lines

    # TODO: Vectorise this to speed it up
    def ray_casting_OLD(self, point):

        n_tri = self.mesh_faces.shape[0]

        P = self.point_out - point
        # rn = nominator, rd = denominator
        rn = np.sum(np.multiply(self.mesh_nrm, self.mesh_v0 - point), axis=1)
        rd = np.dot(self.mesh_nrm, P)

        # r = np.divide(rn,rd)

        intersect_count = 0

        for i in range(0, n_tri):

            if rd[i] == 0:
                if rn[i] == 0:
                    # Parallel and lies in the plane
                    rI = 0.0
                else:
                    # Parallel to plane, but outside. Mark by -1 to avoid counting
                    rI = -1
            else:
                rI = rn[i] / rd[i]

            if 0 <= rI <= 1:
                # Crosses the plane, but is it within triangle?

                w = point + rI * P - self.mesh_v0[i, :]

                si = (self.mesh_uv[i] * np.inner(w, self.mesh_v[i, :])
                      - self.mesh_vv[i] * np.inner(w, self.mesh_u[i, :])) / self.mesh_denom[i]

                if si < 0 or si > 1:
                    # outside of triangle
                    continue

                ti = (self.mesh_uv[i] * np.inner(w, self.mesh_u[i, :])
                      - self.mesh_uu[i] * np.inner(w, self.mesh_v[i, :])) / self.mesh_denom[i]

                if ti < 0 or (si + ti) > 1:
                    # outside of triangle
                    continue

                # print("intersects face i = " + str(i))
                # print("si = " +str(si) + ", ti = " + str(ti))

                intersect_count += 1

        print(f"ray_casting_OLD - intersection count {intersect_count}")
        return np.mod(intersect_count, 2) == 1

    ############################################################################

    def ray_casting(self, point):

        """ Ray-casting, to determine if a point is inside or outside of mesh. """

        return RegionMesh.ray_casting_helper(point=point,
                                             self_mesh_faces=self.mesh_faces,
                                             self_mesh_nrm=self.mesh_nrm,
                                             self_mesh_v0=self.mesh_v0,
                                             self_point_out=self.point_out,
                                             self_mesh_denom=self.mesh_denom,
                                             self_mesh_uv=self.mesh_uv,
                                             self_mesh_uu=self.mesh_uu,
                                             self_mesh_vv=self.mesh_vv,
                                             self_mesh_u=self.mesh_u,
                                             self_mesh_v=self.mesh_v)

    @staticmethod
    @jit(nopython=True)
    def ray_casting_helper(point,
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

        intersect_count = 0

        idx = np.where(np.logical_and(0 <= r, r <= 1))[0]

        w = point + r.reshape(len(r), 1) * p.reshape(1, 3) - self_mesh_v0
        n_points = len(r)

        s = np.divide(np.multiply(self_mesh_uv, np.sum(np.multiply(w, self_mesh_v), axis=1).reshape(n_points, ))
                      - np.multiply(self_mesh_vv, np.sum(np.multiply(w, self_mesh_u), axis=1).reshape(n_points, )),
                      self_mesh_denom)

        t = np.divide(np.multiply(self_mesh_uv, np.sum(np.multiply(w, self_mesh_u), axis=1).reshape(n_points, ))
                      - np.multiply(self_mesh_uu, np.sum(np.multiply(w, self_mesh_v), axis=1).reshape(n_points, )),
                      self_mesh_denom)

        intersect_count = np.sum((0 <= r) * (r <= 1) * (0 <= s) * (s <= 1) * (0 <= t) * (s + t <= 1))

        # print(f"{point} intersection count {intersect_count}")

        # if np.random.uniform() < 0.05:
        #     if (np.mod(intersect_count, 2) == 1) != self.ray_casting_OLD(point):
        #         print(f"New function differs from old for point {point}... check why")
        #         import pdb
        #         pdb.set_trace()

        return np.mod(intersect_count, 2) == 1

    ############################################################################

    def verify_inside(self, num_points=1000):

        """ Verify the check-inside method against the ray-casting method to make sure they give same results. """

        if self.role != "master":
            return

        # Use default random number generator for verification, so different each run -- more chances to catch error
        x_test = np.random.uniform(self.min_coord[0], self.max_coord[0], num_points)
        y_test = np.random.uniform(self.min_coord[1], self.max_coord[1], num_points)
        z_test = np.random.uniform(self.min_coord[2], self.max_coord[2], num_points)

        self.write_log("Verifying our method")

        for i in range(0, num_points):

            self.write_log(f"{i}/{num_points}")
            coords = np.array([x_test[i], y_test[i], z_test[i]])

            try:
                assert (self.check_inside(coords) == self.ray_casting(coords))
            except:
                self.write_log(f"Mismatch for coordinates: {coords}", is_error=True)
                self.write_log(f"Cached: {self.check_inside(coords)}", is_error=True)
                self.write_log(f"RC: {self.ray_casting(coords)}", is_error=True)
                import pdb
                pdb.set_trace()

    ############################################################################

    def setup_place_neurons(self, d_min=None):

        """
        Initialises variables for neuron placement.

        Args:
            d_min (float) : Minimal distance between neuron somas, in SI units (meters)
        """

        self.write_log("Setup place neurons")

        self.d_min = d_min

        self.rand_ctr = self.max_rand + 1

        self.random_pool = np.zeros((self.max_rand, 3))
        self.density_lookup = np.zeros((self.max_rand,))

        self.neuron_coords = np.zeros((self.max_neurons, 3))
        self.neuron_ctr = 0

        self.padding_coords = np.zeros((self.max_neurons, 3))
        self.padding_ctr = 0

        self.neuron_types = dict([])
        self.neuron_type = np.zeros((self.max_neurons,))
        self.next_neuron_type = 1

        self.reject_ctr = 0

        self.num_bins = self.voxel_mask_inner.shape

        self.update_padding_mask()
        self.update_random_pool()

        self.write_log("Setup done")

    ############################################################################

    def setup_voxel_list(self):

        """ Setup voxel list, these are voxels that needs to be checked for neurons close by. """

        self.write_log("Setup voxel list")

        self.voxel_next_neuron = np.zeros(self.num_bins, dtype=int)
        self.voxel_neurons = np.zeros((self.num_bins[0], self.num_bins[1],
                                       self.num_bins[2], self.max_neurons_voxel, 3))

    ############################################################################

    def update_padding_mask(self):

        """ Updates padding mask. We add neurons outside our region of interest, to avoid artificially inflating
            neuron density at the edges (which would happen without the padding region). """

        self.write_log("Update padding mask")

        # Only save padding for border voxels, and voxels nearby
        n_dist = int(np.ceil(2 * self.d_min / self.bin_width))
        s = ndimage.generate_binary_structure(3, 3)

        if np.sum(self.voxel_mask_border):
            # If we do raytracing then the border exists
            self.voxel_mask_padding = np.copy(self.voxel_mask_border)

            self.voxel_mask_padding = \
                scipy.ndimage.binary_dilation(self.voxel_mask_padding,
                                              structure=s,
                                              iterations=n_dist)
        else:
            # No ray tracing, we need to create a padding region

            dilated_mask = scipy.ndimage.binary_dilation(self.voxel_mask_inner,
                                                         structure=s,
                                                         iterations=n_dist)
            self.voxel_mask_padding = \
                np.logical_xor(dilated_mask, self.voxel_mask_inner)

    ############################################################################

    def check_padding_zone(self, coords):

        # TODO: Check/remember why min was taken here :)
        idx = np.array(np.floor((coords - self.min_coord) / self.bin_width), dtype=int)

        return self.voxel_mask_padding[idx[0], idx[1], idx[2]]

    ############################################################################

    def update_random_pool(self):

        """ Refills the random pool with new random numbers. """

        if self.rand_ctr >= self.max_rand:

            # self.write_log("Regenerating new random pool")
            for i in range(0, 3):
                self.random_pool[:, i] = self.random_generator.uniform(low=self.min_coord[i],
                                                                       high=self.max_coord[i],
                                                                       size=self.max_rand)

            self.rand_ctr = 0

    ############################################################################

    # density_function is either None, or a function of pos = [x,y,z] (in SI units)

    def define_density(self, neuron_type, density_function):

        """
        Defines density for neuron type.

        Args:
            neuron_type (str): Neuron type
            density_function (str): density_function is either None, or a function of pos = [x,y,z] (in SI units)

        """

        self.density_function[neuron_type] = density_function
        self.density_voxel_sum[neuron_type] = np.zeros(self.num_bins, dtype=float)
        self.density_total_sum[neuron_type] = 0
        self.density_voxel_n_sample[neuron_type] = np.zeros(self.num_bins, dtype=int)
        self.density_total_n_sample[neuron_type] = 0
        self.placed_voxel[neuron_type] = np.zeros(self.num_bins, dtype=int)
        self.placed_total[neuron_type] = 0

    def place_neurons(self, num_cells, neuron_type=None, d_min=None):

        if d_min is None:
            d_min = self.d_min

        d_min2 = d_min ** 2

        # If this is not fulfilled, then we need to update the range values below
        assert 2 * d_min < self.bin_width, \
            "2*dMin (2 * " + str(d_min) + ") must be smaller than binWidth (" + str(self.bin_width) + ")"

        if neuron_type in self.neuron_types:
            neuron_type_id = self.neuron_types[neuron_type]
        else:
            neuron_type_id = self.next_neuron_type
            self.next_neuron_type += 1

            self.neuron_types[neuron_type] = neuron_type_id

        start_ctr = self.neuron_ctr
        end_ctr = self.neuron_ctr + num_cells

        t_a = timeit.default_timer()

        if neuron_type in self.density_function and neuron_type != self.last_neuron_type_added:
            # Precalculate density function for potential putative locations
            xv = self.random_pool[self.rand_ctr:, 0]
            yv = self.random_pool[self.rand_ctr:, 1]
            zv = self.random_pool[self.rand_ctr:, 2]
            self.density_lookup[self.rand_ctr:] = self.density_function[neuron_type](x=xv, y=yv, z=zv)
            self.last_neuron_type_added = neuron_type

            if np.sum(self.density_lookup[self.rand_ctr:]) == 0:
                self.write_log(f"Density zero for {len(xv)} {neuron_type} neurons -- error with density?",
                               is_error=True)

        while self.neuron_ctr < end_ctr and self.reject_ctr < self.max_reject:

            putative_loc = self.random_pool[self.rand_ctr, :]
            df = self.density_lookup[self.rand_ctr]

            self.rand_ctr += 1

            if self.rand_ctr % 100000 == 0:
                self.write_log(f"Neurons: {self.neuron_ctr} Rejected: {self.reject_ctr} Padding: {self.padding_ctr}")

                if self.neuron_ctr == 0:
                    self.write_log("No neurons placed, check why!", is_error=True)
                    import pdb
                    pdb.set_trace()

            if self.rand_ctr >= self.max_rand:
                self.update_random_pool()

                if neuron_type in self.density_function:
                    xv = self.random_pool[:, 0]
                    yv = self.random_pool[:, 1]
                    zv = self.random_pool[:, 2]
                    self.density_lookup[:] = self.density_function[neuron_type](x=xv, y=yv, z=zv)

                # TODO: We should recalculate griddata density for these points

            inside_flag = self.check_inside(coords=putative_loc)

            if not inside_flag and not self.check_padding_zone(putative_loc):
                # We are outside, go to next neuron
                self.reject_ctr += 1
                continue

            # Check that we are not too close to existing points
            # Only check the neighbouring voxels, to speed things up
            voxel_idx = np.array(np.floor((putative_loc - self.min_coord)
                                          / self.bin_width), dtype=int)

            # Density check is fast, do that to get an early rejection if needed
            # TODO: When a function is defined it is relatively fast, but when griddata is used
            #       it is slooow...
            if inside_flag and neuron_type in self.density_function:
                # Update: df, density function value is now precomputed above
                # xp, yp, zp = putative_loc
                # df = self.density_function[neuron_type](x=xp, y=yp, z=zp)
                assert df >= 0, f"Error your density for {neuron_type} is incorrect, value={df} at {putative_loc}"
                vx, vy, vz = voxel_idx

                self.density_voxel_sum[neuron_type][vx, vy, vz] += df
                self.density_total_sum[neuron_type] += df
                self.density_voxel_n_sample[neuron_type][vx, vy, vz] += 1
                self.density_total_n_sample[neuron_type] += 1

                try:
                    n_expected = (self.density_voxel_sum[neuron_type][vx, vy, vz]
                                  / self.density_total_sum[neuron_type]
                                  * (self.placed_total[neuron_type] + 1))
                except:
                    import traceback
                    t_str = traceback.format_exc()
                    self.write_log(t_str)
                    import pdb
                    pdb.set_trace()

                # This assumes all of mesh voxels have same volume
                # n_expected = ((self.density_voxel_sum[neuron_type][vx, vy, vz]
                #               / self.density_voxel_n_sample[neuron_type][vx, vy, vz])
                #               / (self.density_total_sum[neuron_type] / self.density_total_n_sample[neuron_type])
                #               * (self.placed_total[neuron_type] + 1))

                if self.placed_voxel[neuron_type][vx, vy, vz] > np.ceil(n_expected):
                    # We have too many neurons in this part of the volume already, reject
                    self.reject_ctr += 1
                    continue

            voxel_idx_list = [voxel_idx]

            border_voxel = np.zeros((3,), dtype=int)

            for idx in range(0, 3):  # Looping over x,y,z coordinates
                if (putative_loc[idx] - self.min_coord[idx]) % self.bin_width < d_min:
                    border_voxel[idx] = -1
                    new_idx = np.copy(voxel_idx)
                    new_idx[idx] -= 1
                    voxel_idx_list.append(new_idx)

                elif (putative_loc[idx] - self.min_coord[idx]) % self.bin_width > self.bin_width - d_min:
                    border_voxel[idx] = 1
                    new_idx = np.copy(voxel_idx)
                    new_idx[idx] += 1
                    voxel_idx_list.append(new_idx)

            n_border = np.sum(np.abs(border_voxel))

            if n_border == 2:
                # Along one of the lines, need to check diagonal voxel also
                voxel_idx_list.append(voxel_idx + border_voxel)
            elif n_border == 3:
                # Close to corner, need to check 8 voxels in total (ouch!)
                voxel_idx_list.append(voxel_idx + [border_voxel[0], border_voxel[1], 0])
                voxel_idx_list.append(voxel_idx + [border_voxel[0], 0, border_voxel[2]])
                voxel_idx_list.append(voxel_idx + [0, border_voxel[1], border_voxel[2]])
                voxel_idx_list.append(voxel_idx + border_voxel)

            min_dist2 = 1e6
            for vox_idx in voxel_idx_list:
                if (vox_idx < 0).any() or (vox_idx > self.num_bins).any():
                    # Voxel outside bounds, ignore
                    continue

                if min_dist2 < d_min2:
                    # No need to calculate further, we are too close
                    break

                if self.voxel_next_neuron[vox_idx[0], vox_idx[1], vox_idx[2]] > 0:
                    tmp = self.voxel_neurons[vox_idx[0], vox_idx[1], vox_idx[2],
                                             0:self.voxel_next_neuron[vox_idx[0], vox_idx[1], vox_idx[2]],
                                             :] - putative_loc

                    min_dist2 = min(min_dist2, np.min(np.sum(np.square(tmp), axis=1)))

            if d_min2 < min_dist2:
                # Ok neuron is not too close to any neighbours

                if inside_flag:
                    # We are inside, add to inside points
                    self.neuron_coords[self.neuron_ctr, :] = putative_loc
                    self.neuron_type[self.neuron_ctr] = neuron_type_id
                    self.neuron_ctr += 1
                    # self.writeLog("Placed neuron " + str(self.neuronCtr))

                    # Update counts if we have a density function defined
                    if neuron_type in self.density_function:
                        self.placed_voxel[neuron_type][vx, vy, vz] += 1
                        self.placed_total[neuron_type] += 1
                else:
                    self.padding_ctr += 1

                # Also save the point in the specific voxel, this way we can ignore
                # lots of distance comparisons
                try:
                    self.voxel_neurons[voxel_idx[0], voxel_idx[1], voxel_idx[2],
                                       self.voxel_next_neuron[voxel_idx[0], voxel_idx[1], voxel_idx[2]], :] = putative_loc
                    self.voxel_next_neuron[voxel_idx[0], voxel_idx[1], voxel_idx[2]] += 1
                except:
                    self.write_log(f"If you see this error you probably need to increase " 
                                   f"self.max_neurons_voxel={self.max_neurons_voxel}")
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)
                    import pdb
                    pdb.set_trace()

            else:
                self.reject_ctr += 1

        t_b = timeit.default_timer()
        self.write_log(f"Placed {num_cells} in {t_b - t_a} s")

        if neuron_type in self.placed_voxel:
            if np.max(self.placed_voxel[neuron_type]) < 5 \
                    and neuron_type in self.density_function and self.density_function[neuron_type]:
                self.write_log(f"Warning, mesh_bin_width might be too small to setup accurate {neuron_type} density. "
                               f"Max {neuron_type} in any mesh voxel is {np.max(self.placed_voxel[neuron_type])}",
                               is_error=True)

        return self.neuron_coords[start_ctr:end_ctr, :]

        # Store point in global list, but also reference it in specific voxel

    ############################################################################

    # shape = "cube" or "sphere"
    # radius = radius of sphere, or half the length of side of the cube

    def get_subset(self, centre, radius=None, num_neurons=None, shape="cube",
                   return_idx_flag=False):

        """
        Returns subset of positions, either within a radius, or the closest num_neurons neurons.

        Args:
            centre (float,float,float): Centre of space
            radius (float): Radius if sphere, half-side if cube
            num_neurons (int): Number of neurons (if None, all within radius are returned)
            shape (str): "sphere" or "cube"
            return_idx_flag (bool): Return the indexes instead of coordinates (default False)

        """

        assert ((radius is None) ^ (num_neurons is None)), \
            "Specify one of radius or nNeurons."

        coords = self.neuron_coords[:self.neuron_ctr, :]
        nrn_type = self.neuron_type[:self.neuron_ctr, :]

        if shape == "cube":
            dist = np.amax(abs(coords - centre), axis=1)
        elif shape == "sphere":
            dist = np.sqrt(np.sum(np.square(coords - centre), axis=1))
        else:
            assert False, f"Unknown shape {shape} use cube or sphere"

        sorted_idx = np.argsort(dist)

        if radius is not None:
            self.write_log(f"Using radius {radius}")
            idx = sorted_idx[np.where(dist[sorted_idx] <= radius)]
        else:
            self.write_log(f"Selecting {num_neurons} closest neuron")
            idx = sorted_idx[:num_neurons]

        # Next we need to return them in the order they were originally sorted
        keep_mask = np.zeros((self.neuron_ctr,), dtype=bool)
        keep_mask[idx] = True

        if return_idx_flag:
            return np.where(keep_mask)
        else:
            return coords[keep_mask, :], nrn_type[keep_mask, :]

    ############################################################################

    def plot_struct(self, pdf_name=None):

        """ Plot structure.

        Args:
            pdf_name (str) : Save file name (default None)
        """

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(self.mesh_vec[:, 0],
                   self.mesh_vec[:, 1],
                   self.mesh_vec[:, 2],
                   'black')

        plt.ion()
        plt.show()
        plt.pause(0.001)

        ax.view_init(150, 40)
        plt.draw()
        plt.axis("off")

        if pdf_name is not None:
            plt.savefig(pdf_name)

    ############################################################################

    def plot_neurons(self, plot_idx=None, pdf_name=None):

        """ Plot neurons.

        Args:
            plot_idx (list): Neuron ID to plot
            pdf_name (str): Name of file to save figure to (default None)
            """

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        if plot_idx is None:
            plot_idx = range(0, self.neuron_ctr)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(self.neuron_coords[plot_idx, 0],
                   self.neuron_coords[plot_idx, 1],
                   self.neuron_coords[plot_idx, 2],
                   'black')

        plt.ion()
        plt.show()
        plt.pause(0.001)

        ax.view_init(150, 40)
        plt.draw()
        plt.axis("off")

        if pdf_name is not None:
            plt.savefig(pdf_name)

    ############################################################################

    def test_plot(self):

        """ Test plot"""

        # !!! NEXT ADD dMIN TO THIS

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        n_points = 300

        x_test = np.random.uniform(self.min_coord[0], self.max_coord[0], n_points)
        y_test = np.random.uniform(self.min_coord[1], self.max_coord[1], n_points)
        z_test = np.random.uniform(self.min_coord[2], self.max_coord[2], n_points)

        for i in range(0, n_points):

            self.write_log(f"Checking {i + 1}/{n_points}")

            if self.ray_casting(np.array([x_test[i], y_test[i], z_test[i]])):
                self.write_log("Inside!")
                color = 'red'
                ax.scatter(x_test[i], y_test[i], z_test[i], color=color)
            else:
                color = 'black'

            # ax.scatter(xTest[i],yTest[i],zTest[i],color=color)

        plt.show()
        plt.pause(0.001)

    ############################################################################

    def test_plot_cached(self):

        """ Test plot cached. """

        # !!! NEXT ADD dMIN TO THIS

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        n_points = 1000

        x_test = np.random.uniform(self.min_coord[0], self.max_coord[0], n_points)
        y_test = np.random.uniform(self.min_coord[1], self.max_coord[1], n_points)
        z_test = np.random.uniform(self.min_coord[2], self.max_coord[2], n_points)

        for i in range(0, n_points):

            coord = np.array([x_test[i], y_test[i], z_test[i]])
            ray_val = self.ray_casting(coord)
            cached_val = self.check_inside(coord)

            if ray_val == cached_val:
                self.write_log(f"Checking {i + 1}/{n_points}")

            elif ray_val:
                # Inside, but cached was wrong
                color = 'red'
                ax.scatter(x_test[i], y_test[i], z_test[i], color=color)
            else:
                # Outside, but cached wrong
                color = 'blue'
                ax.scatter(x_test[i], y_test[i], z_test[i], color=color)

        plt.show()
        plt.pause(0.001)

    ############################################################################

    def verify_d_min(self):

        """ Verify that d_min constraint is met. """

        self.write_log("Verifying that dMin constraint is met")

        min_dist = np.zeros((self.neuron_ctr,))

        if self.neuron_ctr < 100000:
            neuron_range = range(0, self.neuron_ctr)
        else:
            self.write_log("Too many to check all, picking random neurons to check")
            neuron_range = np.random.randint(0, self.neuron_ctr, size=(100000, 1))
            self.write_log(str(neuron_range))

        ctr = 0

        for i_neuron in range(0, self.neuron_ctr):

            ctr = ctr + 1
            if ctr % 10000 == 0:
                self.write_log(f"{ctr}/{len(neuron_range)}")

            d = np.sqrt(np.sum(np.square(self.neuron_coords[:self.neuron_ctr, :]
                                         - self.neuron_coords[i_neuron, :]),
                               axis=1))
            d[i_neuron] = 1e6  # Dont count self distance
            min_dist[i_neuron] = np.min(d)

        self.write_log(f"Closest neighbour, min = {np.min(min_dist)}")

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        if True:
            fig = plt.figure()
            plt.hist(min_dist)
            plt.title(f"Closest neighbour {np.min(min_dist)}")
            plt.ylabel("Count")
            plt.xlabel("Distance")
            plt.ion()
            plt.show()
            plt.pause(0.001)

        # Also plot where the neurons are relative to the voxel boundaries
        nr = [x for x in neuron_range]
        bad_idx = [nr[i] for i in np.where(min_dist < self.d_min)[0]]

        fig2 = plt.figure()
        ax = fig2.add_subplot(111, projection='3d')

        xc = (self.neuron_coords[bad_idx, :] - self.min_coord) % self.bin_width
        ax.scatter(xc[:, 0], xc[:, 1], xc[:, 2], 'black')
        plt.title("Bad neuron locations: " + str(xc.shape[0]))
        plt.axis('tight')
        plt.xlabel('X')
        plt.ylabel('Y')
        # plt.zlabel('Z')
        plt.ion()
        plt.show()
        plt.pause(0.001)

        try:
            assert self.d_min <= np.min(min_dist), "dMin criteria not fulfilled: " \
                                                  + str(np.min(min_dist)) + " < dMin = " + str(self.d_min)
        except:
            self.write_log("dMin not fulfilled")
            import pdb
            pdb.set_trace()

    ############################################################################

    def simple_test_case(self):

        self.write_log("This redefines the object, please restart after")

        self.filename = "cube.obj"
        self.cache_file = "cube.obj-cached.pickle"
        self.load_mesh(self.filename)
        self.pre_compute()
        self.setup_voxel_filter()

        self.point_out = np.array([0.5, 0.5, -10.5]) * 1e-6

        for i in range(0, 1000):
            test_point = np.random.uniform(0, 1e-6, 3) \
                        + np.array([0, 0, 0]) * 1e-6

            # testPoint = np.array([0.5,0.61,0.5])*1e-6

            # self.debugFlag = True
            inside_flag = self.ray_casting(test_point)

            if not inside_flag:
                self.write_log(f"Test point = {test_point}", is_error=True)
                self.write_log("wrong!", is_error=True)
                import pdb
                pdb.set_trace()

        self.write_log("All correct")
        self.write_log("Debug mode")
        import pdb
        pdb.set_trace()

    ############################################################################

    def write_log(self, text, flush=True, is_error=False):  # Change flush to False in future, debug

        """
        Writes to log file. Use setup_log first. Text is only written to screen if self.verbose=True,
        or is_error = True, or force_print = True.

        test (str) : Text to write
        flush (bool) : Should all writes be flushed to disk directly?
        is_error (bool) : Is this an error, always written.
        force_print (bool) : Force printing, even if self.verbose=False.
        """

        if self.logfile is not None:
            self.logfile.write(text + "\n")
            if flush:
                self.logfile.flush()

        if self.verbose or is_error:
            print(text)

    ############################################################################

    def inner_voxel_volume(self):

        """ Volume of all inner voxels. """

        return np.sum(self.voxel_mask_inner) * (self.bin_width ** 3)

    ############################################################################


if __name__ == "__main__":

    # sm = RegionMesh("cube.obj",useCache=False)
    # sm.simpleTestCase()

    if os.getenv('IPYTHON_PROFILE') is not None:
        from ipyparallel import Client

        rc = Client(profile=os.getenv('IPYTHON_PROFILE'),
                    # sshserver='127.0.0.1',
                    debug=False)
        print('Client IDs: ' + str(rc.ids))

        # http://davidmasad.com/blog/simulation-with-ipyparallel/
        # http://people.duke.edu/~ccc14/sta-663-2016/19C_IPyParallel.html
        print("Client IDs: " + str(rc.ids))
        d_view = rc.direct_view(targets='all')  # rc[:] # Direct view into clients
    else:
        print("No IPYTHON_PROFILE enviroment variable set, running in serial")
        d_view = None
        rc = None

    meshFile = '../data/mesh/Striatum-d.obj'
    # meshFile = "mesh/cortex-mesh-200.obj"
    sm = RegionMesh(meshFile, d_view=d_view, raytrace_borders=False, verbose=True)

    # import cProfile
    # cProfile.run("neuronPos = sm.placeNeurons(1000)")

    # sm.plotStruct()

    nNeurons = 1730000
    neuronPos = sm.place_neurons(nNeurons)
    # sm.verify_d_min()
    sm.plot_neurons(pdf_name="figures/striatum-fig-somas.png")

    sm.plot_struct(pdf_name="figures/striatum-fig-struct.png")

    # sm.testPlot()
    # sm.testPlotCached()

    # tp = (sm.minCoord + sm.maxCoord)/2
    # sm.rayCasting(np.array(tp))

    if d_view and rc:
        rc.shutdown(hub=True)
