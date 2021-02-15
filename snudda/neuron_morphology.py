# This open source software code was developed in part or in whole in
# the Human Brain Project, funded from the European Union’s Horizon
# 2020 Framework Programme for Research and Innovation under Specific
# Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
# and SGA2).
#
# etc...
#

import numpy as np


class NeuronMorphology(object):

    # axonStumpIDFlag should be True if running Network_simulate.py
    # it should be False if we are running Neurodamus simulation.

    def __init__(self,
                 name=None,
                 position=np.zeros((1, 3)),
                 rotation=None,  # np.eye(3),
                 swc_filename=None,
                 param_filename=None,
                 param_data=None,
                 mech_filename=None,
                 modulation=None,
                 parameter_id=None,
                 modulation_id=None,
                 verbose=False,
                 load_morphology=True,
                 hoc=None,
                 colour=None,
                 use_cache=True,
                 pickle_version=-1,
                 logfile=None,
                 virtual_neuron=False,
                 axon_stump_id_flag=False):

        self.cache_version = 0.9

        self.position = np.copy(np.array(position))

        if rotation is not None:
            self.rotation = np.copy(np.array(rotation))
        else:
            self.rotation = None

        self.soma = np.array([])
        self.axon = np.array([])
        self.dend = np.array([])  # 0,1,2: x,y,z  3: radie, 4: dist to soma (all in meters)

        self.axon_density_type = None
        self.dend_density = None
        self.axon_density = None
        self.axon_density_bounds_xyz = None

        self.voxel_size = 5e6
        self.density_bin_size = 10e-6

        self.load_morphology = load_morphology

        # This tells if axon is indexed fully or if only as a stump
        # this affects sectionID for both axon and dendrites
        self.axon_stump_id_flag = axon_stump_id_flag

        # Meta data
        self.name = name
        self.swc_filename = swc_filename
        self.param_filename = param_filename
        self.param_data = param_data
        self.mech_filename = mech_filename
        self.modulation = modulation
        self.parameter_id = parameter_id
        self.modulation_id = modulation_id
        self.verbose = verbose
        self.use_cache = use_cache
        self.pickle_version = pickle_version
        self.logFile = logfile
        self.virtual_neuron = virtual_neuron

        self.rotated_flag = False

        self.cache_filename = swc_filename.replace('.swc', '-cache.pickle')
        assert (self.cache_filename != swc_filename)

        # This is used for Neurodamus, which instantiates through hoc files
        if hoc is None:
            hoc = ""

        self.hoc = hoc

        # This is useful when determining connectivity, to exclude pairs
        self.max_axon_radius = 0
        self.max_dend_radius = 0

        # Telling how the different points link together into lines
        self.axon_links = np.array((2, 0))  # These should never be changed after CLONE
        self.dend_links = np.array((2, 0))

        self.dend_sec_id = np.array((1,))
        self.dend_sec_x = np.array((2, 0))

        if colour is None:
            self.colour = np.random.random((3,))
        else:
            self.colour = colour

        if load_morphology:
            # This loads, rotates and places neuron
            self.load_neuron_morphology()

    ############################################################################

    def load_neuron_morphology(self):

        if self.use_cache:
            if self.cache_exist():
                # Load existing cache
                try:
                    self.load_cache()
                except Exception as e:

                    # import traceback
                    # tstr = traceback.format_exc()
                    # print(tstr)

                    self.write_log("Failed to read cache file, loading: " + self.swc_filename)
                    self.load_swc(self.swc_filename)
                    self.save_cache()

            else:
                self.write_log("No cache found, create it.")
                # Load SWC and save cache file
                self.load_swc(self.swc_filename)
                self.save_cache()
        else:
            # Load SWC file
            self.write_log("Ignoring old cache, rewriting new cache file")
            self.load_swc(self.swc_filename)
            self.save_cache()

        self.place()  # Updates position and rotation

        # Remove axonStumpIDFlag completely later...
        assert not self.axon_stump_id_flag, \
            "axonStumpFlag is depricated, should be off"

    ############################################################################

    def clone(self,
              load_morphology=None,  # True or False, None = same as parent
              position=np.zeros((1, 3)),
              rotation=None,
              parameter_id=None,
              modulation_id=None):

        if load_morphology is None:
            load_morphology = self.load_morphology

        # If these are explicitly set to None, reuse to original coordinates
        # and rotation
        if position is None:
            position = self.position

        if rotation is None:
            rotation = self.rotation

        new_neuron = NeuronMorphology(name=self.name,
                                      position=position,
                                      rotation=rotation,
                                      swc_filename=self.swc_filename,
                                      param_filename=self.param_filename,
                                      param_data=self.param_data,
                                      mech_filename=self.mech_filename,
                                      modulation=self.modulation,
                                      parameter_id=parameter_id,
                                      modulation_id=modulation_id,
                                      verbose=self.verbose,
                                      load_morphology=False,
                                      hoc=self.hoc,
                                      virtual_neuron=self.virtual_neuron)

        if load_morphology:
            # Set the flag
            new_neuron.load_morphology = load_morphology

            # Copy the data
            new_neuron.axon = np.copy(self.axon)
            new_neuron.dend = np.copy(self.dend)
            new_neuron.soma = np.copy(self.soma)

            # Warn the user if the neuron is already rotated
            new_neuron.rotated_flag = self.rotated_flag

            new_neuron.place()

            # These dont change either, so skip np.copy
            new_neuron.axon_links = self.axon_links
            new_neuron.dend_links = self.dend_links
            new_neuron.dend_sec_x = self.dend_sec_x
            new_neuron.dend_sec_id = self.dend_sec_id

            new_neuron.axon_stump_id_flag = self.axon_stump_id_flag

        new_neuron.max_axon_radius = self.max_axon_radius
        new_neuron.max_dend_radius = self.max_dend_radius

        if self.dend_density is not None:
            new_neuron.dend_density = self.dend_density

        if self.axon_density is not None:
            new_neuron.axon_density = self.axon_density

        new_neuron.voxel_size = self.voxel_size

        if self.axon_density_type is not None:
            new_neuron.axon_density_type = self.axon_density_type

        if self.axon_density_bounds_xyz is not None:
            new_neuron.axon_density_bounds_xyz = self.axon_density_bounds_xyz

        return new_neuron

    ############################################################################

    def write_log(self, text):
        if self.logFile is not None:
            self.logFile.write(f"{text}\n")
            print(text)
        else:
            if self.verbose:
                print(text)

    ############################################################################

    # http://blog.lostinmyterminal.com/python/2015/05/12/random-rotation-matrix.html
    @staticmethod
    def rand_rotation_matrix(deflection=1.0, rand_nums=None):
        """
    Creates a random rotation matrix.
    
    deflection: the magnitude of the rotation. For 0, no rotation; for 1, competely random
    rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, they will be auto-generated.
    """

        # from http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c

        if rand_nums is None:
            rand_nums = np.random.uniform(size=(3,))

        theta, phi, z = rand_nums

        theta = theta * 2.0 * deflection * np.pi  # Rotation about the pole (Z).
        phi = phi * 2.0 * np.pi  # For direction of pole deflection.
        z = z * 2.0 * deflection  # For magnitude of pole deflection.

        # Compute a vector V used for distributing points over the sphere
        # via the reflection I - V Transpose(V).  This formulation of V
        # will guarantee that if x[1] and x[2] are uniformly distributed,
        # the reflected points will be uniform on the sphere.  Note that V
        # has length sqrt(2) to eliminate the 2 in the Householder matrix.

        r = np.sqrt(z)
        vv = (np.sin(phi) * r, np.cos(phi) * r, np.sqrt(2.0 - z))

        st = np.sin(theta)
        ct = np.cos(theta)

        rr = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))

        # Construct the rotation matrix  ( V Transpose(V) - I ) R.
        mm = (np.outer(vv, vv) - np.eye(3)).dot(rr)
        return mm

    ############################################################################

    # We can specify a position and rotation
    def place(self, rotation=None, position=None):

        if self.rotated_flag:
            self.write_log("!!! WARNING, rotating a rotated neuron...")

        if rotation is None:
            rotation = self.rotation
        elif type(rotation) is not np.ndarray:
            rotation = np.array(rotation)

        if position is None:
            position = self.position
        elif type(position) is not np.ndarray:
            position = np.array(position)

        # print("Place called! pos = " + str(position) + ", rot = " + str(rotation))

        # rotation = self.randRotationMatrix()

        # We subtract soma before rotating to centre neuron
        if rotation is not None:

            assert np.abs(np.linalg.det(rotation) - 1) < 1e-6, \
                "place: determinant of rotation matrix should be 1 (did you miss matmul when multiplying?)"

            self.rotated_flag = True

            if len(self.axon) > 0:
                self.axon[:, 0:3] = \
                    np.transpose(np.matmul(rotation, np.transpose(self.axon[:, 0:3] - self.soma[0, 0:3])))
            if len(self.dend) > 0:
                self.dend[:, 0:3] = \
                    np.transpose(np.matmul(rotation, np.transpose(self.dend[:, 0:3] - self.soma[0, 0:3])))
            if len(self.soma) > 0:
                self.soma[:, 0:3] = \
                    np.transpose(np.matmul(rotation, np.transpose(self.soma[:, 0:3] - self.soma[0, 0:3])))

        # Place neuron in correct position
        if len(self.axon) > 0:
            self.axon[:, 0:3] = self.axon[:, 0:3] - self.soma[0, 0:3] + position

        if len(self.dend) > 0:
            self.dend[:, 0:3] = self.dend[:, 0:3] - self.soma[0, 0:3] + position

        if len(self.soma) > 0:
            self.soma[:, 0:3] = self.soma[:, 0:3] - self.soma[0, 0:3] + position

        # Track rotation and location
        self.rotation = rotation
        self.position = position

        # Plot neuron post rotation
        # self.plot_neuron()

        return self

    ############################################################################

    def save_cache(self, cache_file=None):

        if cache_file is None:
            cache_file = self.cache_filename

        assert not self.rotated_flag, \
            "saveCache: The neuron should not be rotated when saving cache"

        morph = dict([])

        morph["swc_filename"] = self.swc_filename
        morph["soma"] = self.soma
        morph["axon"] = self.axon
        morph["dend"] = self.dend
        morph["axonLinks"] = self.axon_links
        morph["dendLinks"] = self.dend_links
        morph["dendSecX"] = self.dend_sec_x
        morph["dendSecID"] = self.dend_sec_id
        morph["axonStumpIDFlag"] = self.axon_stump_id_flag
        morph["maxAxonRadius"] = self.max_axon_radius
        morph["maxDendRadius"] = self.max_dend_radius
        morph["dendDensity"] = self.dend_density
        morph["axonDensity"] = self.axon_density
        morph["version"] = self.cache_version

        assert (cache_file != self.swc_filename)
        print("Saving cache file: " + cache_file)

        import pickle
        with open(cache_file, 'wb') as cache_file:
            pickle.dump(morph, cache_file, self.pickle_version)

    ############################################################################

    def cache_exist(self, cache_file=None):

        if cache_file is None:
            cache_file = self.cache_filename

        cache_flag = False

        import os

        if os.path.isfile(cache_file):

            swc_time = os.path.getmtime(self.swc_filename)
            cache_time = os.path.getmtime(cache_file)

            if cache_time > swc_time:
                print("Found cache file: " + cache_file)
                cache_flag = True
            else:
                print("Found old cache file: " + cache_file)

        else:
            print("No cache file found.")

        return cache_flag

    ############################################################################

    def load_cache(self, cache_file=None):

        if cache_file is None:
            cache_file = self.cache_filename

        import pickle
        with open(cache_file, 'rb') as cache_file:
            morph = pickle.load(cache_file)

        assert self.swc_filename == morph["swc_filename"], \
            "Cached file had different path. Saving new version of cache."

        assert self.axon_stump_id_flag == morph["axonStumpIDFlag"], \
            "axonStumpIDFlag must match cached version"

        # axonStumpIDFlag affects the section ID for the dendrites (and axon)
        # True when running Network_simulate.py and False if running Neurodamus.

        # self.axonStumpIDFlag = morph["axonStumpIDFlag"] # True or False

        self.soma = np.copy(morph["soma"])
        self.axon = np.copy(morph["axon"])
        self.dend = np.copy(morph["dend"])

        self.axon_links = morph["axonLinks"]
        self.dend_links = morph["dendLinks"]
        self.dend_sec_x = morph["dendSecX"]
        self.dend_sec_id = morph["dendSecID"]

        assert morph["version"] == self.cache_version, \
            "Cache version mismatch, regenerating cache"

        self.max_axon_radius = morph["maxAxonRadius"]
        self.max_dend_radius = morph["maxDendRadius"]

        if morph["dendDensity"] is not None:
            self.dend_density = morph["dendDensity"]
        else:
            self.dend_density = None

        if morph["axonDensity"] is not None:
            self.axon_density = morph["axonDensity"]
        else:
            self.axon_density = None

        # Place neuron -- Do not place neuron, loadNeuronMorphology does that
        # self.place()

    ############################################################################

    # self.actionStumpIDFlag only affects the section ID.
    # If it is set to False, all sectionID are computed normally
    # if it is set to True, each axon will have the same sectionID throughout
    # if there are multiple axons they will have separate sectionIDs

    def load_swc(self, swc_file):

        with open(swc_file, 'r') as f:
            lines = f.readlines()

        comp_type = {1: "soma", 2: "axon", 3: "dend", 4: "apic"}

        swc_vals = np.zeros(shape=(len(lines), 7))

        num_comps = 0
        for ss in lines:
            if ss[0] != '#':
                swc_vals[num_comps, :] = [float(s) for s in ss.split()]
                num_comps = num_comps + 1

        # swcVals -- 0: compID, 1: type, 2,3,4: xyz coords, 5: radius, 6: parentID
        assert (1 <= swc_vals[:num_comps, 1]).all() and (swc_vals[:num_comps, 1] <= 4).all(), \
            f"loadMorphology: Only types 1,2,3,4 are supported: {swc_file}"

        # Subtract 1 from ID and parentID, so we get easier indexing
        swc_vals[:, 0] -= 1
        swc_vals[:, 6] -= 1

        swc_vals[:, 2:6] *= 1e-6  # Convert to meter x,y,z, radie

        # Columns:
        # 0: ID, 1,2,3: x,y,z 4: radie, 5: type, 6: parent, 7: somaDist,
        # 8: nodeParent, 9: childCount, 10: sectionID, 11: sectionLen,
        # 12: segmentLen

        # -- careful with sectionID and sectionLen at branch points,
        #    they belong to the parent section
        # -- also dont confuse sectionLen and segmentLen (the latter is
        #    for the segment, which is a part of the larger section)
        points = np.zeros((num_comps, 13))
        points[:num_comps, 0] = swc_vals[:num_comps, 0]  # ID
        points[:num_comps, 1:5] = swc_vals[:num_comps, 2:6]  # x,y,z,r
        points[:num_comps, 5] = swc_vals[:num_comps, 1]  # type
        points[:num_comps, 6] = swc_vals[:num_comps, 6]  # parent

        assert points[0, 5] == 1, \
            "First compartment must be a soma: " + str(swc_file)

        # Create list of the links,
        # exclude soma -> first comp link (should be within soma radius)
        # Columns: 0: ID1, 1: ID2, 2: sectionID, 3: sectionX0, 4: sectionX1
        # 5: nodeParent, 6:type

        links = np.zeros((num_comps, 7))

        link_idx = 0
        for idx in range(0, num_comps):
            id0 = int(points[idx, 6])  # parent
            id1 = int(points[idx, 0])  # point

            if id0 <= 0:
                # No parent or soma is parent, skip link
                continue

            links[link_idx, 0:2] = [id0, id1]
            links[link_idx, 5] = points[idx, 5]

            link_idx += 1

        # Trim link list
        links = links[:link_idx, :]

        # Count children each node has
        for idx in range(1, num_comps):
            try:
                # Increment parents child counter
                points[int(points[idx, 6]), 9] += 1
            except:
                print("Are there gaps in the numbering of the compartments in the SWC file: " + str(swc_file))
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr)
                import pdb
                pdb.set_trace()

        # Make sure soma has a child count > 1 --- no we dont want some as node
        # if(points[0,9] == 0):
        #   points[0,9] = 100

        # Also make sure all points with soma as parent get child count > 1
        # (Child Count > 1 ==> start or end of segment)
        soma_child_idx = np.where(points[:, 6] == 0)[0]
        points[soma_child_idx, 9] += 50

        # Mark node parent, and assign sectionID
        # -- this is used to set sectionID for links, the link
        # will use the end points sectionID
        # !!! Make sure sectionID is correct, and match what Neuron uses internally

        # Nodes are branch points (> 1 child), or end points (0 children)
        # and not soma
        node_idx = np.where((points[:, 9] != 1) & (points[:, 5] != 1))[0]

        # soma is section 0, but we dont include connection soma to first node
        # so let the first dend node be 0, since the section ID is taken from
        # the child ID
        section_id = 1

        # Sonata specifies first axon, then basal, then apical sections
        axon_idx = node_idx[np.where(points[node_idx, 5] == 2)[0]]
        basal_idx = node_idx[np.where(points[node_idx, 5] == 3)[0]]
        apical_idx = node_idx[np.where(points[node_idx, 5] == 4)[0]]

        # If simulation will use an axon stump, where each axon branch is shortened
        # to a stump with the same section ID, then we need to make sure the
        # numbering is correct for the dendrites.

        # Update, set axonID to -1
        for nIdx in axon_idx:
            points[nIdx, 10] = -1

        # Set soma ID to 0
        points[0, 10] = 0

        # Calculate sectionID for dendrites
        section_id = 1

        # Axon dealt with, only loop over dendrites next
        node_loop_list = [basal_idx, apical_idx]

        for idxList in node_loop_list:
            for nIdx in idxList:
                if points[nIdx, 6] > 0:
                    # Set section ID, exclude soma, and compartments bordering to soma
                    points[nIdx, 10] = section_id
                    section_id += 1

        # Assign node parents
        for nIdx in node_idx:

            # Find node parent
            parent_idx = int(points[nIdx, 6])
            # While one child (= no node), keep following parent
            # But stop if parent is soma, or if grandparent is soma
            # !!! Here last link node to soma is not included in neurite morphology
            #     since we assume it is inside the soma
            while points[parent_idx, 9] == 1 and parent_idx > 0 and points[parent_idx, 6] > 0:
                parent_idx = int(points[parent_idx, 6])

            node_parent_idx = parent_idx
            points[nIdx, 8] = node_parent_idx

            section_id = points[nIdx, 10]
            parent_idx = int(points[nIdx, 6])
            while points[parent_idx, 9] == 1 and parent_idx > 0:
                points[parent_idx, 8] = node_parent_idx
                assert points[parent_idx, 10] == 0, "SectionID should be unset prior"
                points[parent_idx, 10] = section_id
                parent_idx = int(points[parent_idx, 6])

        for idx in range(1, num_comps):
            parent_idx = int(points[idx, 6])

            # Calculate soma dist (and also save segLen)
            seg_len = np.sqrt(np.sum((points[idx, 1:4] - points[parent_idx, 1:4]) ** 2))
            points[idx, 7] = points[parent_idx, 7] + seg_len
            points[idx, 12] = seg_len

        # Calculate section length (length between nodes)
        for idx in node_idx:
            node_parent_idx = int(points[idx, 8])

            # Difference in soma distance is section length
            section_len = points[idx, 7] - points[node_parent_idx, 7]
            points[idx, 11] = section_len

            if section_len == 0:
                self.write_log("Section length is zero --- !!! ")
                import pdb
                pdb.set_trace()

            prev_idx = int(points[idx, 6])
            while prev_idx > node_parent_idx:
                points[prev_idx, 11] = section_len
                prev_idx = int(points[prev_idx, 6])

        # Calculate sectionX
        for idx in range(0, links.shape[0]):
            id0 = int(links[idx, 0])
            id1 = int(links[idx, 1])
            links[idx, 2] = points[id1, 10]  # Section ID from point (not parent)

            node_parent = int(points[id1, 8])
            node_parent_soma_dist = points[node_parent, 7]
            section_len = points[id1, 11]

            # segX0 and segX1
            links[idx, 3] = (points[id0, 7] - node_parent_soma_dist) / section_len
            links[idx, 4] = (points[id1, 7] - node_parent_soma_dist) / section_len

            links[idx, 5] = node_parent
            links[idx, 6] = points[id0, 5]  # type (use parent,
            # to avoid soma to dend link)

        # Store the soma, axon, dend and links in the object

        self.soma = np.zeros((1, 4))
        self.soma[0, :] = swc_vals[0, 2:6]  # save x,y,z,r

        dend_idx = np.where((points[:, 5] == 3) | (points[:, 5] == 4))[0]
        axon_idx = np.where(points[:, 5] == 2)[0]

        dend_link_idx = np.where((links[:, 6] == 3) | (links[:, 6] == 4))[0]
        axon_link_idx = np.where(links[:, 6] == 2)[0]

        # 0,1,2: x,y,z  3: radie, 4: dist to soma
        self.dend = np.zeros((len(dend_idx), 5))
        self.axon = np.zeros((len(axon_idx), 5))

        self.dend_links = np.zeros((len(dend_link_idx), 2), dtype=int)  # ID0,ID1
        self.axon_links = np.zeros((len(axon_link_idx), 2), dtype=int)  # ID0,ID1

        self.dend_sec_id = np.zeros((len(dend_link_idx),), dtype=int)  # SectionID
        self.dend_sec_x = np.zeros((len(dend_link_idx), 2))  # SecX0, SecX1

        dend_lookup = dict([])
        axon_lookup = dict([])

        for idx in range(0, len(dend_idx)):
            dend_lookup[dend_idx[idx]] = idx

        for idx in range(0, len(axon_idx)):
            axon_lookup[axon_idx[idx]] = idx

        for idx, dIdx in enumerate(dend_idx):
            self.dend[idx, 0:4] = points[dIdx, 1:5]  # x,y,z,r
            self.dend[idx, 4] = points[dIdx, 7]  # dist to soma

        for idx, aIdx in enumerate(axon_idx):
            self.axon[idx, 0:4] = points[aIdx, 1:5]  # x,y,z,r
            self.axon[idx, 4] = points[aIdx, 7]  # dist to soma

        for idx, dIdx in enumerate(dend_link_idx):
            self.dend_links[idx, 0] = dend_lookup[int(links[dIdx, 0])]  # ID0 - parent
            self.dend_links[idx, 1] = dend_lookup[int(links[dIdx, 1])]  # ID1

            self.dend_sec_id[idx] = links[dIdx, 2]
            self.dend_sec_x[idx, :] = links[dIdx, 3:5]

        for idx, aIdx in enumerate(axon_link_idx):
            self.axon_links[idx, 0] = axon_lookup[links[aIdx, 0]]
            self.axon_links[idx, 1] = axon_lookup[links[aIdx, 1]]
            # We also have sectionID, secX0 and secX1 saved in links[:,2:5]
            # if needed in the future

        if self.virtual_neuron:
            # For virtual neurons, skip the dendrites (save space)
            self.dend = np.zeros((0, self.dend.shape[1]))
            self.dend_links = np.zeros((0, 2))
            self.dend_sec_id = np.zeros((0,))
            self.dend_sec_x = np.zeros((0, 2))

        # self.dendriteDensity() # -- depricated
        self.find_radius()
        self.place()

        # self.debug_plot()

    ############################################################################

    def find_radius(self):

        if len(self.axon) > 0:
            self.max_axon_radius = \
                np.max(np.linalg.norm(self.axon[:, 0:3] - self.soma[0, 0:3], axis=1))

        if len(self.dend) > 0:
            self.max_dend_radius = \
                np.max(np.linalg.norm(self.dend[:, 0:3] - self.soma[0, 0:3], axis=1))
        else:
            self.max_dend_radius = 0

        if self.verbose:
            print("Max axon radius = " + str(self.max_axon_radius))
            print("Max dend radius = " + str(self.max_dend_radius))

    ############################################################################

    def plot_neuron(self,
                    axis=None,
                    plot_axon=True,
                    plot_dendrite=True,
                    line_style='-',
                    alpha=1.0,
                    plot_origo=np.array([0, 0, 0]),  # Only use this when plotting hyper voxels
                    plot_scale=1.0,
                    axon_colour=None,
                    dend_colour=None,
                    soma_colour=None):

        if self.verbose:
            print("Plotting neuron " + self.swc_filename)

        if axon_colour is None:
            axon_colour = self.colour
        if dend_colour is None:
            dend_colour = self.colour
        if soma_colour is None:
            soma_colour = self.colour

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        if axis is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = axis

        if len(self.axon) > 0 and plot_axon:
            ax_links = []

            for row in self.axon_links[:, :2].astype(int):
                if len(ax_links) == 0:
                    ax_links = list(row)
                elif row[0] == ax_links[-1]:
                    ax_links.append(row[1])
                elif row[1] == ax_links[-1]:
                    ax_links.append(row[0])
                else:
                    ax.plot((self.axon[ax_links, 0] - plot_origo[0]) * plot_scale,
                            (self.axon[ax_links, 1] - plot_origo[1]) * plot_scale,
                            (self.axon[ax_links, 2] - plot_origo[2]) * plot_scale,
                            linestyle=line_style,
                            marker=',',
                            alpha=alpha,
                            c=axon_colour)

                    ax_links = list(row)

            if len(ax_links) > 0:
                ax.plot((self.axon[ax_links, 0] - plot_origo[0]) * plot_scale,
                        (self.axon[ax_links, 1] - plot_origo[1]) * plot_scale,
                        (self.axon[ax_links, 2] - plot_origo[2]) * plot_scale,
                        linestyle=line_style,
                        marker=',',
                        alpha=alpha,
                        c=axon_colour)

                # TODO: Also connect first axon compartment to soma

        if plot_dendrite:
            dend_links = []
            for row in self.dend_links[:, :2].astype(int):
                if len(dend_links) == 0:
                    dend_links = list(row)
                elif row[0] == dend_links[-1]:
                    dend_links.append(row[1])
                elif row[1] == dend_links[-1]:
                    dend_links.append(row[0])
                else:
                    ax.plot((self.dend[dend_links, 0] - plot_origo[0]) * plot_scale,
                            (self.dend[dend_links, 1] - plot_origo[1]) * plot_scale,
                            (self.dend[dend_links, 2] - plot_origo[2]) * plot_scale,
                            linestyle=line_style,
                            marker=',',
                            alpha=alpha,
                            c=dend_colour)

                    dend_links = list(row)

            if len(dend_links) > 0:
                ax.plot((self.dend[dend_links, 0] - plot_origo[0]) * plot_scale,
                        (self.dend[dend_links, 1] - plot_origo[1]) * plot_scale,
                        (self.dend[dend_links, 2] - plot_origo[2]) * plot_scale,
                        linestyle=line_style,
                        marker=',',
                        alpha=alpha,
                        c=dend_colour)

                # TODO: Also connect dendrites to soma for plot

        if len(self.soma) > 0:

            if self.soma.shape[0] == 1:
                u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
                x = (self.soma[0, 3] * np.cos(u) * np.sin(v) + self.soma[0, 0] - plot_origo[0])*plot_scale
                y = (self.soma[0, 3] * np.sin(u) * np.sin(v) + self.soma[0, 1] - plot_origo[1])*plot_scale
                z = (self.soma[0, 3] * np.cos(v) + self.soma[0, 2] - plot_origo[2])*plot_scale

                ax.plot_wireframe(x, y, z, color=soma_colour, alpha=alpha)
            else:
                ax.scatter((self.soma[:, 0] - plot_origo[0]) * plot_scale,
                           (self.soma[:, 1] - plot_origo[1]) * plot_scale,
                           (self.soma[:, 2] - plot_origo[2]) * plot_scale,
                           c=soma_colour, alpha=alpha)

        if axis is None:
            plt.title("Neuron: " + self.swc_filename.split("/")[-3] + "_" + self.swc_filename.split('/').pop())
            plt.ion()
            plt.show()
            plt.draw()
            plt.pause(0.001)

        return ax

    ############################################################################

    def set_axon_voxel_radial_density(self, density, max_axon_radius):

        print("Only saving equation now")

        self.axon_density_type = "r"
        self.axon_density = density
        self.max_axon_radius = max_axon_radius

    ############################################################################

    def set_axon_voxel_xyz_density(self,
                                   density,
                                   axon_density_bounds_xyz):

        print("Only saving equation now")

        self.axon_density_type = "xyz"
        self.axon_density = density
        self.axon_density_bounds_xyz = axon_density_bounds_xyz

    ############################################################################

    def compartment_length(self, comp_type="dend"):
        if comp_type == "dend":
            links = self.dend_links
            coords = self.dend
        elif comp_type == "axon":
            links = self.axon_links
            coords = self.axon
        else:
            assert False, "Unknown compartment type: " + str(comp_type) \
                          + ", valid types are 'axon' and 'dend'"

        comp_len = np.linalg.norm(coords[links[:, 0], :][:, :3] - coords[links[:, 1], :][:, :3], axis=1)

        return comp_len

    ############################################################################

    def dendrite_input_locations(self, synapse_density, rng, num_locations=None, return_density=False):

        # Calculate the input density at each point in dendrite morphology
        d = self.dend[:, 4]
        try:
            # d is now distance from some, so synapseDensity is a func of d
            i_density = eval(synapse_density)
        except:
            self.write_log("Bad synapse density string: " + str(synapse_density))
            import traceback
            tstr = traceback.format_exc()
            self.write_log(tstr)
            assert False, f"Problem with synapse density {synapse_density}"

        if type(i_density) in (int, float):
            # If iDensity is a constant, we need to set it for all points
            i_density *= np.ones(d.shape)

        comp_density = (i_density[self.dend_links[:, 0]] + i_density[self.dend_links[:, 1]]) / 2
        comp_len = self.compartment_length(comp_type="dend")

        # compDensity is in synapses per micrometer, multiply by 1e6
        expected_synapses = comp_density * comp_len * 1e6

        if num_locations is not None:
            print("Trying to set nLocations = " + str(num_locations) + " (approx)")
            expected_synapses *= num_locations / np.sum(expected_synapses)

        # Number of input synapses on each compartment
        number_of_synapses = (expected_synapses + ((expected_synapses % 1)
                                                   > rng.random(len(expected_synapses)))).astype(int)

        n_syn_tot = np.sum(number_of_synapses)
        dist_syn_soma = []

        # x,y,z, secID, secX
        input_loc = np.zeros((n_syn_tot, 5))

        # Iterate over each compartment
        syn_ctr = 0
        for i_comp, nSyn in enumerate(number_of_synapses):

            # Add synapses to that compartment
            for j in range(0, nSyn):
                # print('Compartment containing a synapse',iComp)
                # print('Distance from soma',self.dend[iComp][4]*1e6,'$mum$')
                input_loc[syn_ctr, 3] = self.dend_sec_id[i_comp]

                # Cant have at endpoints 0 or 1
                comp_x = rng.random()
                dist_syn_soma = np.append(dist_syn_soma,
                                          d[self.dend_links[i_comp, 0]] * (1 - comp_x)
                                          + d[self.dend_links[i_comp, 1]] * comp_x)

                coords = (self.dend[self.dend_links[i_comp, 0], :3] * (1 - comp_x)
                          + self.dend[self.dend_links[i_comp, 1], :3] * comp_x)

                input_loc[syn_ctr, :3] = coords

                # Use compX (where between comp endpoints) to calculate sectionX
                # (where between section end points)
                input_loc[syn_ctr, 4] = self.dend_sec_x[i_comp, 0] * (1 - comp_x) + comp_x * self.dend_sec_x[i_comp, 1]

                syn_ctr += 1

        if return_density:
            # Return xyz,secID,secX,iDensity,distSynSoma
            return input_loc[:, :3], input_loc[:, 3], input_loc[:, 4], i_density, dist_syn_soma

        # Return xyz,secID,secX, dist_to_soma (update: now also added distance synapse to soma)
        return input_loc[:, :3], input_loc[:, 3], input_loc[:, 4], dist_syn_soma

    ############################################################################

    def debug_plot(self, wait_flag=True, plot_step=1, plot_axon_flag=False):

        ax = self.plot_neuron(plot_axon=plot_axon_flag)

        if plot_axon_flag:
            for a in self.axon_links:
                x0 = self.axon[int(a[0]), 0:3]
                x1 = self.axon[int(a[1]), 0:3]
                x = (x0 + x1) / 2

                # ax.text(x=x0[0],y=x0[1],z=x0[2],s=str(np.around(a[3],2)),color='blue')
                # ax.text(x=x1[0],y=x1[1],z=x1[2],s=str(np.around(a[4],2)),color='red')
                ax.text(x=x[0], y=x[1], z=x[2], s=str(a[2]), color='black')

                print("ID: " + str(a[2]))
                input(" ")

        ctr = 0
        for (d, dID, dX) in zip(self.dend_links, self.dend_sec_id, self.dend_sec_x):
            x0 = self.dend[int(d[0]), 0:3]
            x1 = self.dend[int(d[1]), 0:3]
            x = (x0 + x1) / 2

            # ax.text(x=x0[0],y=x0[1],z=x0[2],s=str(np.around(dX[0],2)),color='blue')
            # ax.text(x=x1[0],y=x1[1],z=x1[2],s=str(np.around(dX[1],2)),color='red')

            if ctr % plot_step == 0:
                ax.text(x=x[0], y=x[1], z=x[2], s=str(dID), color='black')
            ctr += 1

            print("ID: " + str(dID) + " X = " + str(np.around(dX[0], 2)) + " - "
                  + str(np.around(dX[1], 2)))

            if wait_flag:
                input(" ")

        return ax

    ############################################################################


if __name__ == "__main__":
    # The lines below are just for testing purposes

    fName = "data/cellspecs/dspn/str-dspn-e150917_c9_d1-mWT-1215MSN03-v20190521/WT-1215MSN03-cor-rep-ax2.swc"
    # fName = "data/cellspecs/lts/LTS_Experiment-9862_20181211/lts_morp_2019-11-07_centered_noAxon.swc"

    nm = NeuronMorphology(swc_filename=fName, verbose=True, use_cache=False)

    nm.place(rotation=nm.rand_rotation_matrix(), position=np.array([0, 0, 0]))

    nm.debug_plot()

    # nm.setAxonDensity("3e9*np.exp(-d/100e-6)",300e-6)

    # nm.plotDensity()

    ax1 = nm.plot_neuron()

    print("In main function")
    import pdb

    pdb.set_trace()

    nm2 = nm.clone(rotation=nm.rand_rotation_matrix(), position=np.array([0.001, 0.001, 0.001]))
    nm2.plot_neuron(ax1)

    # raw_input("Test")

    import pdb

    pdb.set_trace()
