# This script is custom written to handle very large datasets. It does so by
# not keeping all the information in memory, instead parsing the HDF5
# piece by piece

import numpy as np
import scipy.sparse as sps
import h5py
import timeit
import time
from glob import glob
import os
import sys
import json
import pickle

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D

from .load import SnuddaLoad


# !!! We need to parallelise the analysis script also!

class SnuddaAnalyse(object):

    # saveCache = should we save a pickled file with connection matrix?
    # loadCache = should we load cache file if available
    # lowMemory = if false, uses dense matrix which is faster (assuming lots of memory)

    def __init__(self,
                 hdf5_file=None,
                 load_cache=False,
                 save_cache=True,
                 low_memory=False,
                 side_len=250e-6,
                 volume_type="cube",
                 n_max_analyse=None,
                 show_plots=False,
                 close_plots=True):  # "cube" or "full"

        self.debug = False
        self.show_plots = show_plots

        print("Assuming volume type: " + str(volume_type) + "[cube or full]")

        self.volume_type = volume_type
        self.close_plots = close_plots

        if n_max_analyse is None:
            if volume_type == "cube":
                n_max_analyse = 20000
            elif volume_type == "full":
                n_max_analyse = 20000

        print("Only using " + str(n_max_analyse) + "neurons of the connection data")

        self.num_max_analyse = n_max_analyse

        self.populations = None
        self.dend_position_bin = None
        self.allTypes = None
        self.neuronTypeID = None

        if hdf5_file is None or hdf5_file == "last":
            hdf5_file = self.find_latest_file()

        self.fig_dir = os.path.dirname(hdf5_file) + "/figures"
        if not os.path.exists(self.fig_dir):
            os.makedirs(self.fig_dir)

        # First load all data but synapses
        self.network_load = SnuddaLoad(hdf5_file, load_synapses=False)
        self.network = self.network_load.data

        if "config" in self.network:
            self.config = json.loads(self.network["config"])
        self.side_len = side_len

        self.low_memory = low_memory

        self.neuron_name_remap = {"FSN": "FS"}

        cache_loaded = False
        if load_cache:
            try:
                cache_loaded = self.load_cache_data(hdf5_file)
            except:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)
                assert not cache_loaded, "Load failed, cacheLoaded flag should be False"

        self.data = h5py.File(hdf5_file, 'r')

        if not cache_loaded:
            self.num_neurons = self.network["nNeurons"]
            print("Number of neurons: " + str(self.num_neurons))

            # GABA connection matrix (synType = 1) (ignore AMPA/NMDA = 2, GJ = 3)
            # self.connectionMatrix = self.createConnectionMatrix(synType=1,
            #                                                    lowMemory=lowMemory)

            self.connection_matrix = self.create_connection_matrix(low_memory=low_memory)
            self.connection_matrix_gj = self.create_connection_matrix_gj()

            # self.connectionMatrix = self.createConnectionMatrixSLOW(synType=1)
            self.make_pop_dict()
            self.positions = self.network["neuronPositions"]

            self.synapse_dist()

            if save_cache:
                self.save_cache_data(hdf5_file)

        self.worker_data = []

        self.neuron_colors = {"dSPN": (77. / 255, 151. / 255, 1.0),
                              "iSPN": (67. / 255, 55. / 255, 181. / 255),
                              "FSN": (6. / 255, 31. / 255, 85. / 255),
                              "ChIN": (252. / 266, 102. / 255, 0.0),
                              "LTS": (150. / 255, 63. / 255, 212. / 255),
                              "default": [0.4, 0.4, 0.4]}

    ############################################################################

    def neuron_name(self, neuron_type):

        if neuron_type in self.neuron_name_remap:
            return self.neuron_name_remap[neuron_type]
        else:
            return neuron_type

    ############################################################################

    def get_neuron_color(self, neuron_type):

        if neuron_type in self.neuron_colors:
            return self.neuron_colors[neuron_type]
        else:
            return self.neuron_colors["default"]

    ############################################################################

    # Reading the HDF5 files takes a lot of time, this stores a cached copy
    # of the connection matrix, positions and populations

    def save_cache_data(self, hdf5_file):

        import h5py

        cache_file = hdf5_file + "-cache"

        print("Saving cache to " + cache_file)
        out_file = h5py.File(cache_file, 'w', libver='latest')

        # Connection Matrix
        out_file.create_dataset("conMat_data", data=self.connection_matrix.data,
                                compression='gzip')
        out_file.create_dataset("conMat_indices",
                                data=self.connection_matrix.indices,
                                compression='gzip')
        out_file.create_dataset("conMat_indptr",
                                data=self.connection_matrix.indptr,
                                compression='gzip')
        out_file.create_dataset("conMat_shape",
                                data=self.connection_matrix.shape)

        # GJ connection matrix
        out_file.create_dataset("conMatGJ_data", data=self.connection_matrix_gj.data,
                                compression='gzip')
        out_file.create_dataset("conMatGJ_indices",
                                data=self.connection_matrix_gj.indices,
                                compression='gzip')
        out_file.create_dataset("conMatGJ_indptr",
                                data=self.connection_matrix_gj.indptr,
                                compression='gzip')
        out_file.create_dataset("conMatGJ_shape",
                                data=self.connection_matrix_gj.shape)

        pop_group = out_file.create_group("populations")
        for k in self.populations:
            v = self.populations[k]
            pop_group.create_dataset(k, data=v)

        out_file["nNeurons"] = self.num_neurons
        out_file.create_dataset("positions", data=self.positions)

        try:

            dend_pos_bin = dict([])
            for prePost in self.dend_position_bin:
                pp = self.allTypes[prePost[0]] + "_" + self.allTypes[prePost[1]]
                dend_pos_bin[pp] = list(self.dend_position_bin[prePost])

            out_file.create_dataset("dendPositionBin",
                                    data=json.dumps(dend_pos_bin))
            out_file.create_dataset("dendPositionEdges",
                                    data=self.dend_position_edges)

            all_types = [x.encode("ascii", "ignore") for x in self.allTypes]
            out_file.create_dataset("allTypes",
                                   data=all_types)
            out_file.create_dataset("neuronTypeID", data=self.neuronTypeID)

            out_file.create_dataset("nMaxAnalyse", data=self.num_max_analyse)

        except Exception as e:

            import traceback
            tstr = traceback.format_exc()
            print(tstr)

            print("Problem with writing hdf5")
            import pdb
            pdb.set_trace()

        out_file.close()

    ############################################################################

    def load_cache_data(self, hdf5_file):

        import os
        import h5py

        cache_file = hdf5_file + "-cache"
        data_loaded = False

        if os.path.exists(cache_file):
            t_orig = os.path.getmtime(hdf5_file)
            t_cache = os.path.getmtime(cache_file)

            # Make sure cache file is newer than data file
            if t_cache > t_orig:
                print("Loading from " + cache_file)

                try:
                    with h5py.File(cache_file, 'r') as data:

                        assert self.num_max_analyse == data["nMaxAnalyse"][()], \
                            "nMaxAnalyse has changed, have to reload connection matrix"

                        self.connection_matrix = sps.csr_matrix((data["conMat_data"],
                                                                 data["conMat_indices"],
                                                                 data["conMat_indptr"]),
                                                                data["conMat_shape"])

                        self.connection_matrix_gj = sps.csr_matrix((data["conMatGJ_data"],
                                                                    data["conMatGJ_indices"],
                                                                    data["conMatGJ_indptr"]),
                                                                   data["conMatGJ_shape"])

                        self.populations = dict([])

                        for k in data["populations"].keys():
                            self.populations[k] = data["populations"][k][:]

                            self.num_neurons = data["nNeurons"][()]
                            self.positions = data["positions"][:]

                        dend_pos_bin = json.loads(data["dendPositionBin"][()])
                        self.dend_position_bin = dict([])

                        all_types = list(data["allTypes"][()])
                        self.allTypes = [x.decode() for x in all_types]
                        self.neuronTypeID = data["neuronTypeID"][()]

                        for pp in dend_pos_bin:
                            p_str = pp.split("_")
                            pre_type = self.allTypes.index(p_str[0])
                            post_type = self.allTypes.index(p_str[1])

                            self.dend_position_bin[(pre_type, post_type)] \
                                = np.array(dend_pos_bin[pp])

                        self.dend_position_edges = data["dendPositionEdges"][()]

                        # import pdb
                        # pdb.set_trace()

                    data_loaded = True
                    print("Loading done.")

                except Exception as e:
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)
                    print("Failed to load cache file.")
                    data_loaded = False
                    # import pdb
                    # pdb.set_trace()

        return data_loaded

    ############################################################################

    def create_connection_matrix(self, chunk_size=1000000,
                                 syn_type=None,
                                 min_dend_dist=None,
                                 max_dend_dist=None,
                                 low_memory=False):

        t0 = timeit.default_timer()

        if low_memory:
            print("Trying to conserve memory, this is slower.")
            connection_matrix = sps.lil_matrix((self.num_neurons, self.num_neurons),
                                               dtype=np.int16)
        else:
            try:
                connection_matrix = np.zeros((self.num_neurons, self.num_neurons),
                                             dtype=np.int16)
            except:
                print("Unable to allocate full matrix, using sparse matrix instead")
                connection_matrix = sps.lil_matrix((self.num_neurons, self.num_neurons),
                                                   dtype=np.int16)
        last_src_id = 0
        last_dest_id = 0
        last_count = 0

        row_ctr = 0
        num_syn_total = self.network["nSynapses"]

        for synapses in self.network_load.synapse_iterator(chunk_size=chunk_size):

            print("Synapse row " + str(row_ctr)
                  + " - " + str(100 * row_ctr / float(num_syn_total)) + " %"
                  + " time: " + str(timeit.default_timer() - t0) + " seconds")

            for synRow in synapses:

                row_ctr += 1

                src_id = synRow[0]
                dest_id = synRow[1]  # New format!
                dend_dist = synRow[8]

                if ((min_dend_dist is not None and dend_dist < min_dend_dist)
                        or (max_dend_dist is not None and dend_dist > max_dend_dist)):
                    # Not correct distance to soma on dendrite
                    continue

                if syn_type is None or syn_type == synRow[6]:
                    # Only include specific synapse type

                    if last_src_id == src_id and last_dest_id == dest_id:
                        last_count += 1
                    else:
                        # Write the previous set of synapses to matrix
                        # For first iteration of loop, lastCount is zero
                        connection_matrix[last_src_id, last_dest_id] += last_count

                        last_src_id = src_id
                        last_dest_id = dest_id
                        last_count = 1

                    # connectionMatrix[srcID,destID] += 1

        # Update the last row also
        connection_matrix[last_src_id, last_dest_id] += last_count
        last_count = 0

        t1 = timeit.default_timer()

        print("Created connection matrix " + str(t1 - t0) + " seconds")

        return sps.csr_matrix(connection_matrix, dtype=np.int16)

    ############################################################################

    def create_connection_matrix_gj(self):

        t0 = timeit.default_timer()

        connection_matrix_gj = sps.lil_matrix((self.num_neurons, self.num_neurons),
                                              dtype=np.int16)

        last_src_id = 0
        last_dest_id = 0
        last_count = 0

        row_ctr = 0
        num_gj_total = self.network["nGapJunctions"]

        for gjList in self.network_load.gap_junction_iterator(chunk_size=100000):

            if num_gj_total > 0:
                print(f"GJ row : {row_ctr} - {100 * row_ctr / float(num_gj_total)} % "
                      f" time : {timeit.default_timer() - t0} seconds")

            for gjRow in gjList:
                row_ctr += 1

                src_id = gjRow[0]
                dest_id = gjRow[1]

                if last_src_id == src_id and last_dest_id == dest_id:
                    last_count += 1
                else:
                    connection_matrix_gj[last_src_id, last_dest_id] += last_count

                    last_src_id = src_id
                    last_dest_id = dest_id
                    last_count = 1

            connection_matrix_gj[last_src_id, last_dest_id] += last_count
            last_count = 0

        t1 = timeit.default_timer()
        print("Created gap junction connection matrix " + str(t1 - t0) + " seconds")

        return sps.csr_matrix(connection_matrix_gj, dtype=np.int16)

    ############################################################################

    def make_pop_dict(self):

        print("Creating population dictionary")

        self.populations = dict([])

        for nid, neuron in enumerate(self.network["neurons"]):

            assert (nid == neuron["neuronID"])
            name = neuron["name"].split("_")[0]

            if name not in self.populations:
                self.populations[name] = []

            self.populations[name].append(neuron["neuronID"])

        print("Done.")

    ############################################################################

    def get_sub_pop(self, volume_type="cube", volume_part="centre", side_len=None,
                    neuron_id=None, volume_id="Striatum", num_max_analyse=None):

        # print("volumeType=" + volumeType + ",volumePart=" + volumePart + ",sideLen=" +str(sideLen))

        if volume_type == "full":
            # return all neurons

            if volume_id is not None:
                idx = np.where([x["volumeID"] == volume_id
                                for x in self.network["neurons"]])[0]

                if neuron_id is None:
                    neuron_id = idx

        elif volume_type == "cube":

            if volume_part == "centre":
                try:
                    neuron_id = self.centre_neurons(side_len=side_len,
                                                    neuron_id=neuron_id,
                                                    volume_id=volume_id)
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)
                    import pdb
                    pdb.set_trace()

            elif volume_part == "corner":
                neuron_id = self.corner_neurons(side_len=side_len,
                                                neuron_id=neuron_id,
                                                volume_id=volume_id)
        else:

            print("Unknown volume type: " + str(volume_type))
            import pdb
            pdb.set_trace()

        if num_max_analyse is None:
            num_max_analyse = self.num_max_analyse

        if num_max_analyse is not None:

            if num_max_analyse < len(neuron_id):

                try:
                    keep_idx = np.linspace(0, len(neuron_id), num_max_analyse,
                                           endpoint=False, dtype=int)
                    print("Returning subset of neurons to analyse:"
                          + str(len(keep_idx)) + "/" + str(len(neuron_id)))
                    neuron_id = np.array([neuron_id[x] for x in keep_idx])
                except:
                    import traceback
                    tstr = traceback.format_exc()
                    print(tstr)

                    print("no no no wrong")
                    import pdb
                    pdb.set_trace()

        return neuron_id

    ############################################################################

    def centre_neurons(self, side_len=None, neuron_id=None, volume_id="Striatum"):

        if side_len is None:
            side_len = self.side_len

        if volume_id is None:
            idx = np.arange(0, self.network["nNeurons"])
        else:
            idx = np.where([x["volumeID"] == volume_id
                            for x in self.network["neurons"]])[0]

        min_coord = np.min(self.network["neuronPositions"][idx, :], axis=0)
        max_coord = np.max(self.network["neuronPositions"][idx, :], axis=0)

        x_min = min_coord[0]
        y_min = min_coord[1]
        z_min = min_coord[2]

        x_max = max_coord[0]
        y_max = max_coord[1]
        z_max = max_coord[2]

        x_centre = (x_max + x_min) / 2.0
        y_centre = (y_max + y_min) / 2.0
        z_centre = (z_max + z_min) / 2.0

        if neuron_id is None:
            neuron_id = idx

        if side_len < 0:
            print("We want to have all but the outher layer, assume a cube")
            buf_len = -side_len
            side_len = np.min([x_max - x_centre - buf_len,
                               y_max - y_centre - buf_len,
                               z_max - z_centre - buf_len])

            assert side_len > 0, "Unable to autodetect a good side len"

        if side_len is None:
            return neuron_id

        c_id = []

        for nid in neuron_id:
            # pos = self.network["neurons"][nid]["position"]
            pos = self.positions[nid, :]

            assert volume_id is None \
                   or self.network["neurons"][nid]["volumeID"] == volume_id, \
                   "Neuron " + str(nid) + " does not belong to volumeID " + str(volume_id)

            if (abs(pos[0] - x_centre) <= side_len
                    and abs(pos[1] - y_centre) <= side_len
                    and abs(pos[2] - z_centre) <= side_len):
                c_id.append(nid)

        print("Centering in " + str(volume_id) + " : Keeping " + str(len(c_id)) + "/" + str(len(neuron_id)))

        return c_id

    ############################################################################

    # If we use a corner, and account for missing 1/8th of the surrounding
    # synapses then we can get larger distances.
    #
    # <--->

    def corner_neurons(self, side_len=None, neuron_id=None, volume_id="Striatum"):

        if side_len is None:
            side_len = self.side_len

        if volume_id is None:
            idx = np.arange(0, self.network["nNeurons"])
        else:
            idx = np.where([x["volumeID"] == volume_id
                            for x in self.network["neurons"]])[0]

        if len(idx) == 0:
            print("No neurons found in volume " + str(volume_id))

            import pdb
            pdb.set_trace()

        min_coord = np.min(self.network["neuronPositions"][idx, :], axis=0)
        max_coord = np.max(self.network["neuronPositions"][idx, :], axis=0)

        x_min = min_coord[0]
        y_min = min_coord[1]
        z_min = min_coord[2]

        x_max = max_coord[0]
        y_max = max_coord[1]
        z_max = max_coord[2]

        if (side_len > x_max - x_min
                or side_len > y_max - y_min
                or side_len > z_max - z_min):
            print("Warning: the analysis cube specified by sideLen is too large.")
            print("!!! Setting sideLen to None")

            side_len = None

        if neuron_id is None:
            neuron_id = idx

        if side_len is None:
            return neuron_id

        c_id = []

        for nid in neuron_id:
            # pos = self.network["neurons"][nid]["position"]
            pos = self.positions[nid, :]

            # We assume centre is at zero
            if ((pos[0] <= x_min + side_len or pos[0] >= x_max - side_len)
                    and (pos[1] <= y_min + side_len or pos[1] >= y_max - side_len)
                    and (pos[2] <= z_min + side_len or pos[2] >= z_max - side_len)):
                c_id.append(nid)

        print("Taking corner neurons: Keeping " + str(len(c_id)) + "/" + str(len(neuron_id)))

        if False:

            # Debug plot
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            pos = self.data["network"]["neurons"]["position"][()]

            ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], 'black')
            ax.scatter(pos[c_id, 0], pos[c_id, 1], pos[c_id, 2], 'red')

            plt.ion()
            # plt.draw()

            if self.show_plots:
                plt.show()

            plt.pause(0.001)
            import pdb
            pdb.set_trace()

        return c_id

    ############################################################################

    def save_figure(self, plt, fig_name, fig_type="pdf"):

        if not os.path.isdir(self.fig_dir):
            os.mkdir(self.fig_dir)

        full_fig_name = self.fig_dir + "/" + fig_name + "." + fig_type

        # Remove part of the frame
        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["top"].set_visible(False)

        plt.tight_layout()
        plt.pause(0.001)
        plt.savefig(full_fig_name)
        # plt.savefig(full_fig_name.replace('.pdf', '.eps'))

        print("Wrote " + full_fig_name)

        if self.close_plots:
            time.sleep(1)
            plt.close()

    ############################################################################

    def plot_num_synapses_per_pair(self, pre_type, post_type, side_len=None,
                                   name_str="", volume_id="Striatum",
                                   connection_type="synapses"):

        if side_len is None:
            side_len = self.side_len

        print("Plotting number of connections")

        if pre_type not in self.populations:
            print(f"plotNumSynapsesPerPair: {pre_type} is not in the simulation")
            return

        if post_type not in self.populations:
            print(f"plotNumSynapsesPerPair: {post_type} is not in the simulation")
            return

        pre_pop = self.populations[pre_type]
        post_pop = self.populations[post_type]

        if connection_type == "synapses":
            con_mat = self.connection_matrix
        elif connection_type == "gapjunctions":
            con_mat = self.connection_matrix_gj
        else:
            con_mat = None
            print(f"Unknown connection_type: {connection_type}")
            print("Please use 'synapses' or 'gapjunctions'")
            sys.exit(-1)

        if side_len is not None:
            # We are only looking at post synaptic neurons at the centre,
            # to avoid edge effects
            print(f"Only analysing centre post synaptic neurons, sideLen = {side_len}")
            # postPop = self.centreNeurons(neuronID=postPop,sideLen=sideLen)
            post_pop = self.get_sub_pop(volume_type=self.volume_type,
                                        volume_part="centre",
                                        side_len=side_len,
                                        neuron_id=post_pop,
                                        volume_id=volume_id)

        print("Calculating max synapses")
        max_synapses = con_mat[pre_pop, :][:, post_pop].max()

        print("Calculating mean synapses")

        # The prune tuning func might set data to 0, we want to exclude those
        mean_synapses = float(np.sum(con_mat[pre_pop, :][:, post_pop].data)) \
                              / np.sum(con_mat[pre_pop, :][:, post_pop].data != 0)

        con = con_mat[pre_pop, :][:, post_pop]

        existing_con = con[np.nonzero(con)].transpose()

        # Any connections? Otherwise skip plot
        if ((type(existing_con) == np.matrixlib.defmatrix.matrix
             and len(existing_con) == 0)
                or (type(existing_con) != np.matrixlib.defmatrix.matrix
                    and existing_con.getnnz() == 0)):
            return

        print(f"Plotting {existing_con.shape[0]} connections")

        plt.figure()
        matplotlib.rcParams.update({'font.size': 22})

        plt.hist(existing_con,
                 range(0, 1 + max_synapses),
                 density=True,
                 align="left",
                 color=self.get_neuron_color(pre_type))

        plt.xlabel("Number of " + connection_type)
        plt.ylabel('Probability density')
        plt.title(self.neuron_name(pre_type) + " to " + self.neuron_name(post_type))

        plt.tight_layout()
        plt.draw()

        fig_name = f"Network-number-of-{connection_type}-from-{pre_type}-to-{post_type}-per-cell"

        self.save_figure(plt, fig_name)

        if self.show_plots:
            plt.ion()
            plt.show()
            plt.pause(0.001)

    ############################################################################

    def plot_connection_probability_parallel(self,
                                             pre_type=None,
                                             post_type=None,
                                             num_bins=86,
                                             dist_3d=True):

        assert pre_type is not None
        assert post_type is not None

        print(f"Plotting connection probability {pre_type} to {post_type}")

        if pre_type not in self.populations:
            print(f"plotConnectionProbability: {pre_type} is not in the simulation")
            return

        if post_type not in self.populations:
            print(f"plotConnectionProbability: {post_type} is not in the simulation")
            return

        pre_id = self.populations[pre_type]
        post_id = self.populations[post_type]

        # We need to split the work between multiple workers
        import threading
        num_threads = 4
        threads = []

        self.worker_data = []

        for iThread in range(0, num_threads):
            worker_pre_id = pre_id[range(iThread, len(pre_id), num_threads)]
            print(f"Worker {iThread} PreID: {worker_pre_id}")
            t = threading.Thread(target=self.connection_probability_wrapper,
                                 args=(worker_pre_id, post_id, num_bins, 1000000.0, dist_3d))
            threads.append(t)
            print(f"Starting {t.getName()}")
            t.start()

        for t in threads:
            print(f"Joining {t.getName()}")
            t.join()

        # Gather all the data
        dist = self.worker_data[0][0]
        count_con = np.zeros((num_bins, 1))
        count_all = np.zeros((num_bins, 1))

        for data in self.worker_data:
            count_con += data[2]
            count_all += data[3]

        p_con = np.divide(count_con, count_all)

        # Now let's plot it
        matplotlib.rcParams.update({'font.size': 22})
        plt.figure()
        plt.plot(dist * 1e6, p_con)
        plt.xlabel("Distance ($\mu$m)")
        plt.ylabel("Connection probability")

        plt.title(f"{self.neuron_name(pre_type)} to {self.neuron_name(post_type)} connections")
        plt.tight_layout()

        plt.xlim([0, 250])
        plt.ylim([0, 1])

        plt.ion()
        plt.draw()

        plt.pause(0.001)
        fig_name = f"Network-distance-dependent-connection-probability-{pre_type}-to-{post_type}"

        self.save_figure(plt, fig_name)

        if self.show_plots:
            plt.show()

    ############################################################################

    def plot_connection_probability(self,
                                    pre_type=None,
                                    post_type=None,
                                    num_bins=86,
                                    name_str="",
                                    side_len=None,
                                    exp_max_dist=None,
                                    exp_data=None,
                                    exp_data_detailed=None,
                                    dist_3d=True,
                                    volume_id=None,
                                    x_max=250,
                                    y_max=None,
                                    connection_type="synapses",
                                    draw_step=False):

        assert pre_type is not None
        assert post_type is not None

        if not exp_max_dist:
            exp_max_dist = []

        if not exp_data:
            if exp_data_detailed:
                exp_data = [x[0]/x[1] for x in exp_data_detailed]
            else:
                exp_data = []

        if not exp_data_detailed:
            exp_data_detailed = None
        else:
            assert (np.array(exp_data) == np.array([x[0]/x[1] for x in exp_data_detailed])).all(), \
                f"exp_data = {exp_data }and exp_data_detailed {exp_data_detailed} do not match"

        if side_len is None:
            side_len = self.side_len

        if pre_type not in self.populations or post_type not in self.populations:
            print(f"Missing {pre_type} or {post_type} in network, skipping plot with their connectivity")
            return

        print(f"Plotting connection probability {pre_type} to {post_type} ({connection_type})")

        pre_id = self.populations[pre_type]
        post_id = self.populations[post_type]

        # We can in principle use all pairs, but here we restrict to just the
        # pairs who's post partner are in the centre
        # postID = self.centreNeurons(neuronID=postID,sideLen=sideLen)
        post_id = self.get_sub_pop(volume_type=self.volume_type,
                                   volume_part="centre",
                                   side_len=side_len,
                                   neuron_id=post_id,
                                   volume_id=volume_id)

        if pre_type not in self.populations:
            print(f"plotConnectionProbabilityChannels: {pre_type} is not in the simulation")
            return

        if post_type not in self.populations:
            print(f"plotConnectionProbabilityChannels: {post_type} is not in the simulation")
            return

        if len(exp_data) == 0 and len(exp_data_detailed) > 0:
            exp_data = []
            for x in exp_data_detailed:
                exp_data.append(x[0] / float(x[1]))

        if len(exp_data_detailed) == 0 and len(exp_data) > 0:
            exp_data_detailed = [None for x in exp_data]

        (dist, p_con, count_con, count_all) = \
            self.connection_probability(pre_id, post_id, num_bins, dist_3d=dist_3d,
                                        connection_type=connection_type)

        # Now let's plot it

        # fig = plt.figure()
        fig, ax = plt.subplots(1)

        matplotlib.rcParams.update({'font.size': 24})

        plt_ctr = 0

        # Add lines for experimental data and matching data for model
        for (d_limit, p_exp, exp_num) in zip(exp_max_dist, exp_data, exp_data_detailed):
            cnt = 0
            cnt_all = 0

            for (d, c, ca) in zip(dist, count_con, count_all):
                if d <= d_limit:
                    cnt += c
                    cnt_all += ca

            # Hack to avoid divide by zero
            cnt_all[cnt_all == 0] = 1

            p_model = float(cnt) / float(cnt_all)

            print("P(d<" + str(d_limit) + ")=" + str(p_model))
            # ax = fig.get_axes()

            # Also add errorbars
            if exp_num is not None:
                P = exp_num[0] / float(exp_num[1])

                # Exp data specified should match
                assert p_exp is None or p_exp == P

                # stdExp = np.sqrt(P*(1-P)/expNum[1])

                # https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval
                # https://www.tandfonline.com/doi/abs/10.1080/01621459.1927.10502953
                # Wilson score
                # Wilson, Edwin B. "Probable inference, the law of succession, and statistical inference." Journal of the American Statistical Association 22.158 (1927): 209-212.
                ns = exp_num[0]
                n = exp_num[1]
                z = 1.96  # This gives us 95% confidence intervall
                bar_centre = (ns + (z ** 2) / 2) / (n + z * 2)
                bar_height = z / (n + z ** 2) * np.sqrt((ns * (n - ns) / n + (z ** 2) / 4))

                # plt.errorbar(dLimit*1e6/2,P,stdExp,color="gray",
                #             elinewidth=1,capsize=5)
                plt.errorbar(d_limit * 1e6 / 2, bar_centre, bar_height, color="gray",
                             elinewidth=1, capsize=5)

                # import pdb
                # pdb.set_trace()

                # rectExpStd = patches.Rectangle((0,P-stdExp),
                #                               width=dLimit*1e6,height=2*stdExp,
                #                               alpha=0.2,linewidth=0,
                #                               color="red",fill=True)
                # ax.add_patch(rectExpStd)
            else:
                std_exp = 0

            if p_exp is not None:
                plt.plot([0, d_limit * 1e6], [p_exp, p_exp],
                         color=(0.8, 0.3 * plt_ctr, 0.3 * plt_ctr), linewidth=2)

                # Add a star also
                plt.plot(d_limit * 1e6 / 2, p_exp,
                         color=(0.8, 0.3 * plt_ctr, 0.3 * plt_ctr),
                         marker="D",
                         markersize=10)

                plt_ctr += 1
                plt.ion()
                plt.draw()

                if self.show_plots:
                    plt.show()

                # rectExp = patches.Rectangle((0,0), dLimit*1e6, Pexp, \
                #                            linewidth=2,color='red',fill=False)
                # ax.add_patch(rectExp)

            if False:  #
                # Show the binning from the model data
                plt.plot([0, d_limit * 1e6], [p_model, p_model], color="blue", linewidth=2)

            # rect = patches.Rectangle((0,0), dLimit*1e6, P, \
            #                         linewidth=1,color='blue',fill=False)
            # ax.add_patch(rect)

            if False:
                plt.text(d_limit * 1e6, p_model,
                         "P(d<" + str(d_limit * 1e6) + ")=%.3f" % p_model, size=9)

        # Draw the curve itself
        if draw_step:
            plt.step(dist * 1e6, p_con, color='black', linewidth=2, where="post")
        else:
            d_half_step = (dist[1] - dist[0]) / 2
            plt.plot((dist + d_half_step) * 1e6, p_con, color='black', linewidth=2)

        plt.xticks(fontsize=14, rotation=0)
        plt.yticks(fontsize=14, rotation=0)

        label_size = 22

        # if(dist_3d):
        #  plt.xlabel("Distance ($\mu$m)",fontsize=labelSize)
        # else:
        #  plt.xlabel("2D Distance ($\mu$m)",fontsize=labelSize)

        # Hack to avoid divide by zero
        count_all_b = count_all.copy()
        count_all_b[count_all_b == 0] = 1.0

        # This gives us 95% confidence intervall
        z = 1.96

        p_centre = np.array([(ns + (z ** 2) / 2) / (n + z * 2)
                            for (ns, n) in zip(count_con, count_all_b)]).flatten()
        p_height = np.array([z / (n + z ** 2)
                            * np.sqrt((ns * (n - ns) / n + (z ** 2) / 4))
                            for (ns, n) in zip(count_con, count_all_b)]).flatten()

        # Use the last bin larger than xMax as the end
        d_idx = np.where(dist * 1e6 > x_max)[0][0]

        p_min = p_centre - p_height
        p_max = p_centre + p_height

        if draw_step:
          plt.fill_between(dist[:d_idx] * 1e6, p_min[:d_idx], p_max[:d_idx],
                           color='grey', step="post",
                           alpha=0.4)
        else:
            plt.fill_between((dist[:d_idx] + d_half_step) * 1e6, p_min[:d_idx], p_max[:d_idx],
                             color='grey', step=None,
                             alpha=0.4)

        plt.xlabel("Distance ($\mu$m)", fontsize=label_size)
        plt.ylabel("Con Prob (%)", fontsize=label_size)

        # Adjust axis if needed
        # import pdb
        # pdb.set_trace()

        # Set plot limits
        # if(any(x is not None for x in expData)):
        #  if(max(plt.ylim()) < max(expData+stdExp)):
        #    plt.ylim([0, np.ceil(max(expData+stdExp)*10)/10])

        if x_max is not None:
            plt.xlim([0, x_max])

        if y_max is not None:
            plt.ylim([0, y_max])

        locs, labels = plt.yticks()
        new_labels = ["{0:g}".format(yy * 100) for yy in locs]

        if locs[0] < 0:
            locs = locs[1:]
            new_labels = new_labels[1:]

        if locs[-1] > plt.ylim()[1]:
            locs = locs[:-1]
            new_labels = new_labels[:-1]

        plt.yticks(locs, new_labels)

        plt.title(f"{self.neuron_name(pre_type)} to {self.neuron_name(post_type)}")
        plt.tight_layout()
        plt.ion()
        plt.draw()

        if dist_3d:
            proj_text = '-3D-dist'
        else:
            proj_text = '-2D-dist'

        fig_name = (f"Network-distance-dependent-connection-probability-{pre_type}"
                    f"-to-{post_type}-{connection_type}{proj_text}")

        self.save_figure(plt, fig_name)

        if self.show_plots:
            plt.show()

        plt.pause(0.001)

    ############################################################################

    # expMaxDist=[50e-6,100e-6]
    # expData=[None, None] (none if no exp data, or values)
    # expDataDetailed=[(conA,totA),(conB,totB)]

    def plot_connection_probability_channels(self,
                                             pre_type=None,
                                             post_type=None,
                                             num_bins=86,
                                             name_str="",
                                             side_len=None,  # buffer around edges if neg
                                             exp_max_dist=None,
                                             exp_data=None,
                                             exp_data_detailed=None,
                                             dist_3d=True,
                                             volume_id="Striatum"):

        if pre_type not in self.populations:
            print(f"plotConnectionProbabilityChannels: {pre_type} is not in the simulation")
            return

        if post_type not in self.populations:
            print(f"plotConnectionProbabilityChannels: {post_type} is not in the simulation")
            return

        if not exp_max_dist:
            exp_max_dist = []

        if not exp_data:
            exp_data = []

        if not exp_data_detailed:
            exp_data_detailed = None

        if side_len is None:
            side_len = self.side_len

        if len(exp_data) == 0 and len(exp_data_detailed) > 0:
            exp_data = []
            for x in exp_data_detailed:
                exp_data.append(x[0] / float(x[1]))

        if len(exp_data_detailed) == 0:
            exp_data_detailed = []
            for x in exp_data:
              exp_data_detailed.append(None)

        assert pre_type is not None
        assert post_type is not None

        print(f"Plotting connection probability {pre_type} to {post_type}")

        pre_id = self.populations[pre_type]
        post_id = self.populations[post_type]

        # We can in principle use all pairs, but here we restrict to just the
        # pairs who's post partner are in the centre
        # postID = self.centreNeurons(neuronID=postID,sideLen=sideLen)
        post_id = self.get_sub_pop(volume_type=self.volume_type,
                                   volume_part="centre",
                                   side_len=side_len,
                                   neuron_id=post_id,
                                   volume_id=volume_id)

        (dist, p_con_within, p_con_between,
         count_con_within, count_con_between,
         count_all_within, count_all_between) = \
            self.connection_probability_population_units(pre_id, post_id, num_bins, dist_3d=dist_3d)

        # Now let's plot it

        # fig = plt.figure()
        fig, ax = plt.subplots(1)
        matplotlib.rcParams.update({'font.size': 22})

        # Draw the curve itself
        plt.plot(dist * 1e6, p_con_within, color='black', linewidth=2)
        plt.plot(dist * 1e6, p_con_between, color='grey', linewidth=2)

        label_size = 14
        if dist_3d:
            plt.xlabel("Distance ($\mu$m)", fontsize=label_size)
        else:
            plt.xlabel("2D Distance ($\mu$m)", fontsize=label_size)
        plt.ylabel("Connection probability", fontsize=label_size)

        plt.xticks(fontsize=12, rotation=0)
        plt.yticks(fontsize=12, rotation=0)

        # plt.tick_params(axis='both', which='major', labelsize=10)
        # plt.tick_params(axis='both', which='minor', labelsize=8)

        # import pdb
        # pdb.set_trace()

        # Add lines for experimental data and matching data for model
        for (dLimit, Pexp, expNum) in zip(exp_max_dist, exp_data, exp_data_detailed):
            cnt_within = 0
            cnt_all_within = 0
            cnt_between = 0
            cnt_all_between = 0

            for (d, c, ca) in zip(dist, count_con_within, count_all_within):
                if d <= dLimit:
                    cnt_within += c
                    cnt_all_within += ca

            for (d, c, ca) in zip(dist, count_con_between, count_all_between):
                if d <= dLimit:
                    cnt_between += c
                    cnt_all_between += ca

            # Hack to avoid divide by zero, since the corresponding cntWithin
            # bin will also be zero, the resulting division is zero
            cnt_all_within[cnt_all_within == 0] = 1
            cnt_all_between[cnt_all_between == 0] = 1

            p_within = float(cnt_within) / float(cnt_all_within)
            p_between = float(cnt_between) / float(cnt_all_between)

            p_total = float(cnt_within + cnt_between) / float(cnt_all_within + cnt_all_between)

            print(f"Pwithin(d<{dLimit}) = {p_within}")
            print(f"Pbetween(d<{dLimit}) = {p_between}")
            print(f"Ptotal(d<{dLimit}) = {p_total}")
            # ax = fig.get_axes()

            if Pexp is not None:
                plt.plot([0, dLimit * 1e6], [Pexp, Pexp], color="red", linewidth=2)

                old_rect_flag = False
                if old_rect_flag:
                    # Old rectangle
                    rect_exp = patches.Rectangle((0, 0), dLimit * 1e6, Pexp,
                                                 linewidth=2, color='red', fill=False,
                                                 linestyle="--")
                # Also add errorbars
                if expNum is not None:
                    assert False, "Should change the error bars to wilson score!"
                    p = expNum[0] / float(expNum[1])
                    std_exp = np.sqrt(p * (1 - p) / expNum[1])
                    rect_exp_std = patches.Rectangle((0, Pexp - std_exp),
                                                     width=dLimit * 1e6, height=2 * std_exp,
                                                     alpha=0.2, linewidth=0,
                                                     color="red", fill=True)

                    ax.add_patch(rect_exp_std)

                if old_rect_flag:
                    # Add P-line on top
                    ax.add_patch(rect_exp)

            if False:
                rect_within = patches.Rectangle((0, 0), dLimit * 1e6, p_within,
                                                linewidth=1, color='blue', fill=False)
                rect_between = patches.Rectangle((0, 0), dLimit * 1e6, p_between,
                                                 linewidth=1, color='lightblue',
                                                 fill=False)
                rect_total = patches.Rectangle((0, 0), dLimit * 1e6, p_total,
                                               linewidth=2, color='blue', fill=False)

                ax.add_patch(rect_within)
                ax.add_patch(rect_between)
                ax.add_patch(rect_total)

            if True:
                # Draw binned data as lines instead
                plt.plot([0, dLimit * 1e6], [p_within, p_within],
                         color="blue", linewidth=1)
                plt.plot([0, dLimit * 1e6], [p_between, p_between],
                         color="lightblue", linewidth=1)
                plt.plot([0, dLimit * 1e6], [p_total, p_total],
                         color="blue", linewidth=2)

            if False:
                Plt(dLimit * 1e6, p_within, "Pw(d<" + str(dLimit * 1e6) + ")=%.3f" % p_within, size=9)
                plt.text(dLimit * 1e6, p_between, "Pb(d<" + str(dLimit * 1e6) + ")=%.3f" % p_between, size=9)

        # Adjust axis if needed
        # import pdb
        # pdb.set_trace()

        if any(x is not None for x in exp_data):
            if max(plt.ylim()) < max(exp_data):
                plt.ylim([0, np.ceil(max(exp_data) * 10) / 10])

        plt.xlim([0, 250])
        # plt.xlim([0, 1000])

        plt.title(self.neuron_name(pre_type) + " to "
                  + self.neuron_name(post_type) + " connections")

        plt.tight_layout()
        plt.ion()
        plt.draw()

        if dist_3d:
            proj_text = '-3D-dist'
        else:
            proj_text = '-2D-dist'

        fig_name = 'Network-distance-dependent-connection-probability-channels-' \
                   + str(pre_type) + "-to-" + str(post_type) \
                   + proj_text

        self.save_figure(plt, fig_name)

        if self.show_plots:
            plt.show()

        plt.pause(0.001)

    ############################################################################

    def connection_probability_wrapper(self, pre_id, post_id, num_bins=86, dist_3d=True):

        print("Worker started, preID: " + str(pre_id))

        (dist, p_con, count_con, count_all) = self.connection_probability(pre_id,
                                                                          post_id,
                                                                          num_bins,
                                                                          dist_3d=dist_3d)

        self.worker_data.append((dist, p_con, count_con, count_all))

    ############################################################################

    # connection_type: "synapses" or "gapjunctions"

    def connection_probability(self,
                               pre_id,
                               post_id,
                               num_bins=86,
                               num_points=10000000.0,
                               dist_3d=True,
                               connection_type="synapses"):

        # Count the connected neurons
        print("Counting connections")
        dist = np.linspace(0.0, 1700.0e-6, num=num_bins)
        delta_dist = dist[1] - dist[0]

        count_con = np.zeros((num_bins, 1))
        count_all = np.zeros((num_bins, 1))
        count_rejected = 0

        num_per_pair = num_points / len(pre_id)

        if connection_type == "synapses":
            con_mat = self.connection_matrix
        elif connection_type == "gapjunctions":
            con_mat = self.connection_matrix_gj
        else:
            assert False,"Unknown connection_type: " + str(connection_type)

        # Make this loop use threads, to speed it up

        for xi, x in enumerate(pre_id):
            t_a = timeit.default_timer()

            if num_per_pair - np.floor(num_per_pair) > np.random.rand():
                n_pts = int(np.ceil(num_per_pair))
            else:
                n_pts = int(np.floor(num_per_pair))

            po_id = np.random.permutation(post_id)[0:min(n_pts, len(post_id))]

            for y in po_id:  # postID:

                if x == y:
                    # Do not count self connections in statistics!!
                    # This can lead to what looks like an artificial drop in
                    # connectivity proximally
                    continue

                if dist_3d:
                    d = np.sqrt(np.sum((self.positions[x, :] - self.positions[y, :]) ** 2))
                else:
                    d = np.sqrt(np.sum((self.positions[x, 0:2] - self.positions[y, 0:2]) ** 2))
                    # We also need to check that z-distance is not too large

                    # Gilad email 2017-11-21:
                    # The 100 um is the lateral (XY) distance between somata of
                    # recorded cells. In terms of the Z axis, the slice
                    # thickness is 250 um but the recordings are all done in the
                    # upper half of the slice due to visibility of the cells. So
                    # lets say in a depth of roughly 40-110 um from the upper
                    # surface of the slice.
                    dz = np.abs(self.positions[x, 2] - self.positions[y, 2])

                    # Using dzMax = 70, see comment above from Gilad.
                    if dz > 70e-6:
                        # Skip this pair, too far away in z-depth
                        count_rejected += 1
                        continue

                idx = int(np.floor(d / delta_dist))
                assert (idx < num_bins), "Idx too large " + str(idx)

                if con_mat[x, y] > 0:
                    count_con[idx] += 1

                count_all[idx] += 1

            t_b = timeit.default_timer()

            if self.debug:
                print(str(xi + 1) + "/" + str(len(pre_id)) + " " + str(t_b - t_a) + "s")

        p_con = np.divide(count_con, count_all)

        print("Requested: " + str(num_points) + " calculated " + str(sum(count_all)))

        if not dist_3d:
            print("Rejected (too large z-depth): " + str(count_rejected))

        # import pdb
        # pdb.set_trace()

        return dist, p_con, count_con, count_all

    ############################################################################

    def connection_probability_population_units(self,
                                                pre_id,
                                                post_id,
                                                num_bins=86,
                                                num_points=5000000.0,
                                                dist_3d=True):

        population_unit = self.network["populationUnit"]

        # Count the connected neurons
        print("Counting connections")
        dist = np.linspace(0.0, 1700.0e-6, num=num_bins)
        delta_dist = dist[1] - dist[0]

        count_con_within_channel = np.zeros((num_bins, 1))
        count_all_within_pop_unit = np.zeros((num_bins, 1))

        count_con_between_pop_units = np.zeros((num_bins, 1))
        count_all_between_pop_units = np.zeros((num_bins, 1))

        count_rejected = 0

        num_per_pair = num_points / len(pre_id)

        # Make this loop use threads, to speed it up

        for xi, x in enumerate(pre_id):
            t_a = timeit.default_timer()

            if num_per_pair - np.floor(num_per_pair) > np.random.rand():
                n_pts = int(np.ceil(num_per_pair))
            else:
                n_pts = int(np.floor(num_per_pair))

            po_id = np.random.permutation(post_id)[0:min(n_pts, len(post_id))]

            for y in po_id:  # postID:

                if x == y:
                    # Dont include self-self
                    continue

                if dist_3d:
                    d = np.sqrt(np.sum((self.positions[x, :] - self.positions[y, :]) ** 2))
                else:
                    d = np.sqrt(np.sum((self.positions[x, 0:2] - self.positions[y, 0:2]) ** 2))
                    # We also need to check that z-distance is not too large

                    # Gilad email 2017-11-21:
                    # The 100 um is the lateral (XY) distance between somata of
                    # recorded cells. In terms of the Z axis, the slice
                    # thickness is 250 um but the recordings are all done in the
                    # upper half of the slice due to visibility of the cells. So
                    # lets say in a depth of roughly 40-110 um from the upper
                    # surface of the slice.
                    dz = np.abs(self.positions[x, 2] - self.positions[y, 2])

                    # Using dzMax = 70, see comment above from Gilad.
                    if dz > 70e-6:
                        # Skip this pair, too far away in z-depth
                        count_rejected += 1
                        continue

                idx = int(np.floor(d / delta_dist))

                if idx >= num_bins:
                    count_rejected += 1
                    continue

                # assert (idx < nBins), "Idx too large " + str(idx)

                same_population_unit_flag = (population_unit[x] == population_unit[y])

                if self.connection_matrix[x, y] > 0:
                    if same_population_unit_flag:
                        count_con_within_channel[idx] += 1
                    else:
                        count_con_between_pop_units[idx] += 1

                if same_population_unit_flag:
                    count_all_within_pop_unit[idx] += 1
                else:
                    count_all_between_pop_units[idx] += 1

            t_b = timeit.default_timer()

            if self.debug:
                print(str(xi + 1) + "/" + str(len(pre_id)) + " " + str(t_b - t_a) + "s")

        pcon_within_pop_unit = np.divide(count_con_within_channel, count_all_within_pop_unit)
        pcon_between_pop_unit = np.divide(count_con_between_pop_units,
                                          count_all_between_pop_units)

        print("Requested: " + str(num_points) + " calculated " + str(
            sum(count_all_within_pop_unit) + sum(count_all_between_pop_units)))

        if not dist_3d:
            print("Rejected (too large z-depth): " + str(count_rejected))

        return (dist, pcon_within_pop_unit, pcon_between_pop_unit,
                count_con_within_channel, count_con_between_pop_units,
                count_all_within_pop_unit, count_all_between_pop_units)

    ############################################################################

    def num_incoming_connections(self, neuron_type, pre_type, side_len=None,
                                 volume_id="Striatum",
                                 connection_type="synapses"):

        if side_len is None:
            side_len = 100e-6

        print("Calculating number of incoming connections " + pre_type
              + " -> " + neuron_type)

        # Only use post synaptic cell in central part of structure,
        # to minimize edge effect
        # neuronID = self.centreNeurons(neuronID=self.populations[neuronType],
        #                              sideLen=sideLen)

        neuron_id = self.get_sub_pop(volume_type=self.volume_type,
                                     volume_part="centre",
                                     side_len=side_len,
                                     neuron_id=self.populations[neuron_type],
                                     volume_id=volume_id)

        pre_id = self.populations[pre_type]

        print("#pre = " + str(len(pre_id)) + ", #post = " + str(len(neuron_id)))

        if connection_type == "synapses":
            con_mat = self.connection_matrix
        elif connection_type == "gapjunctions":
            con_mat = self.connection_matrix_gj
        else:
            assert "Unknown connection_type: " + str(connection_type)
            con_mat = None  # To get pycharm to shut up ;)

        n_con = np.zeros((len(neuron_id), 1))
        n_syn = np.zeros((len(neuron_id), 1))

        for i, nID in enumerate(neuron_id):
            cons = con_mat[:, nID][pre_id, :]
            n_con[i] = np.sum(cons > 0)
            n_syn[i] = np.sum(cons)

        return n_con, n_syn

    ############################################################################

    def plot_incoming_connections(self, pre_type, neuron_type, side_len=None,
                                  name_str="",
                                  connection_type="synapses",
                                  fig=None, colour=None, hist_range=None,
                                  bin_size=None,
                                  num_bins=None):

        if pre_type not in self.populations:
            print("plotIncomingConnections: " + str(pre_type)
                  + " is not in the simulation")
            return

        if neuron_type not in self.populations:
            print("plotIncomingConnections: " + str(neuron_type)
                  + " is not in the simulation")
            return

        (n_con, nSyn) = self.num_incoming_connections(neuron_type=neuron_type,
                                                      pre_type=pre_type,
                                                      side_len=side_len,
                                                      connection_type=connection_type)

        # Plotting number of connected neighbours
        if fig is None:
            fig = plt.figure()
            hist_type = 'bar'
        else:
            plt.sca(fig.axes[0])
            hist_type = 'step'

        matplotlib.rcParams.update({'font.size': 22})

        if sum(n_con) > 0:

            if bin_size is None:
                if max(n_con) > 100:
                    bin_size = 10
                else:
                    bin_size = 1

            if colour is None:
                colour = self.get_neuron_color(pre_type)

            if hist_range is None:
                hist_range = range(0, int(np.max(n_con)), bin_size)

            if num_bins is None:
                plt.hist(n_con, hist_range,
                         align="left", density=True,
                         color=colour, histtype=hist_type)
            else:
                plt.hist(n_con, bins=num_bins,
                         align="left", density=True,
                         color=colour, histtype=hist_type)

        plt.xlabel("Number of connected neighbours")
        plt.ylabel("Probability density")
        plt.title(self.neuron_name(pre_type) + " connecting to "
                  + self.neuron_name(neuron_type))
        plt.tight_layout()
        xleft, xright = plt.xlim()
        plt.xlim(0, xright)

        plt.ion()
        plt.draw()

        fig_name = "Network-" + connection_type + "-input-to-" \
                   + neuron_type + "-from-" + pre_type

        self.save_figure(plt, fig_name)

        if self.show_plots:
            plt.show()

        plt.pause(0.001)

        # Plotting number of input synapses
        plt.figure()
        matplotlib.rcParams.update({'font.size': 22})

        if sum(nSyn) > 0:
            if max(nSyn) > 700:
                bin_size = 100
            elif max(nSyn) < 10:
                bin_size = 1
            elif max(nSyn) < 40:
                bin_size = 5
            else:
                bin_size = 10
            plt.hist(nSyn, range(0, int(np.max(nSyn)), bin_size),
                     align="left", density=True,
                     color=self.get_neuron_color(pre_type))

        plt.xlabel("Number of incoming " + connection_type)
        plt.ylabel("Probability density")
        plt.title(self.neuron_name(pre_type) + " " + connection_type
                  + " on " + self.neuron_name(neuron_type))
        plt.tight_layout()

        xleft, xright = plt.xlim()
        plt.xlim(0, xright)

        plt.ion()
        plt.draw()

        fig_name = "Network-" + connection_type + "-to-" + neuron_type + "-from-" + pre_type

        self.save_figure(plt, fig_name)

        if self.show_plots:
            plt.show()

        plt.pause(0.001)

        return fig

    ############################################################################

    def find_synapses(self, neuron_id, max_synapses=10000, chunk_size=1000000):

        print("Finding synapses between: " + str(neuron_id))

        synapses = np.zeros((max_synapses, 8))
        synapse_ctr = 0

        num_rows = self.data["network/synapses"].shape[0]
        chunk_size = min(num_rows, chunk_size)

        num_steps = int(np.ceil(num_rows / chunk_size))
        row_start = 0

        for rowEnd in np.linspace(chunk_size, num_rows, num_steps, dtype=int):

            print("Rows: " + str(row_start) + ":" + str(rowEnd)
                  + " total: " + str(num_rows))

            syn = self.data["network/synapses"][row_start:rowEnd, :]

            for synRow in syn:
                if synRow[0] in neuron_id and synRow[1] in neuron_id:
                    synapses[synapse_ctr, :] = synRow
                    synapse_ctr += 1

        print("Found " + str(synapse_ctr) + " synapses")

        print(str(synapses[0:synapse_ctr, :]))

        return synapses[0:synapse_ctr, :]

    ############################################################################

    # Plots neuronID neuron, and all presynaptic partners

    def plot_neurons(self, neuron_id, show_synapses=True, plot_pre_neurons=True):

        axis = None

        # Finding synapses, this might take time
        (synapses, synapse_coords) = self.network_load.find_synapses(postID=neuron_id)

        assert (synapses[:, 1] == neuron_id).all(), "!!!! Wrong synapses extracted"

        neurons = dict([])

        post_neuron = self.network_load.load_neuron(neuron_id)
        axis = post_neuron.plot_neuron(axis=axis, plot_axon=False,
                                       plot_dendrite=True)

        plotted_neurons = [neuron_id]

        if plot_pre_neurons:
            for synRow in synapses:

                src_id = int(synRow[0])

                if src_id not in plotted_neurons:
                    neurons[src_id] = self.network_load.load_neuron(src_id)
                    axis = neurons[src_id].plot_neuron(axis=axis, plot_axon=True,
                                                       plot_dendrite=False)
                    plotted_neurons.append(src_id)

        if show_synapses:
            axis.scatter(synapse_coords[:, 0],
                         synapse_coords[:, 1],
                         synapse_coords[:, 2], c='red')

        plt.ion()
        plt.draw()

        if self.show_plots:
            plt.show()

        plt.pause(0.001)

        return axis

    ############################################################################

    def find_latest_file(self):

        files = glob('save/network-connect-voxel-pruned-synapse-file-*.hdf5')
        #    files = glob('save/network-connect-synapse-file-*.hdf5')

        mod_time = [os.path.getmtime(f) for f in files]
        idx = np.argsort(mod_time)

        print("Using the newest file: " + files[idx[-1]])

        return files[idx[-1]]

    ############################################################################

    # Loop through all synapses, create distance histogram

    def synapse_dist(self, side_len=None, volume_id="Striatum"):

        if side_len is None:
            side_len = self.side_len

        t_a = timeit.default_timer()

        # cornerID = self.cornerNeurons(sideLen=sideLen)
        corner_id = self.get_sub_pop(volume_type=self.volume_type,
                                     volume_part="corner",
                                     side_len=side_len,
                                     neuron_id=None,
                                     volume_id=volume_id)

        is_corner = np.zeros((self.num_neurons,), dtype=bool)
        is_corner[corner_id] = 1

        print("Calculating synapse distance histogram")

        n_bins = 200
        max_dist = 1000e-6
        bin_width = max_dist / n_bins
        bin_scaling = 1e-6 / bin_width  # To avoid a division

        neuron_type = [n["name"].split("_")[0] for n in self.network["neurons"]]
        self.allTypes = np.unique(neuron_type).tolist()

        # To be quicker, we temporary define neuronTypeIDs
        self.neuronTypeID = [self.allTypes.index(x) for x in neuron_type]

        self.dend_position_bin = dict([])
        self.dend_position_edges = np.linspace(0, max_dist, n_bins + 1)

        for pre_type_id in range(0, len(self.allTypes)):
            for post_type_id in range(0, len(self.allTypes)):
                self.dend_position_bin[(pre_type_id, post_type_id)] = np.zeros((n_bins,))

        print("Creating dist histogram")

        if self.volume_type == "full":
            syn_increment = 1
        elif self.volume_type == "cube":
            syn_increment = 8  # Becaues we use 1/8th of volume,
            # count all synapses 8 times
        else:
            print("Unknown volume type: " + str(self.volume_type))
            import pdb
            pdb.set_trace()

        n_synapses = self.data["network"]["synapses"].shape[0]
        chunk_size = 1000000
        start_idx = 0
        end_idx = min(chunk_size, n_synapses)

        # Outer while loop is to split the data processing into smaller chunks
        # otherwise we use too much memory.
        while start_idx < n_synapses:

            assert start_idx <= end_idx, "startIdx = " + str(start_idx) + ", endIdx = " + str(end_idx)

            all_pre_id = self.data["network"]["synapses"][start_idx:end_idx, 0]
            all_post_id = self.data["network"]["synapses"][start_idx:end_idx, 1]

            # Dist to soma on dendrite
            all_dist_idx = np.array(np.floor(self.data["network"]["synapses"][start_idx:end_idx, 8]
                                             * bin_scaling),
                                    dtype=np.int)

            last_pre = None
            last_post = None

            print("nSynapses = " + str(n_synapses) + ", at " + str(start_idx))

            for i in range(0, all_pre_id.shape[0]):

                if not is_corner[all_post_id[i]]:
                    # We only want to include the corner neurons
                    continue

                pre_type_id = self.neuronTypeID[all_pre_id[i]]
                post_type_id = self.neuronTypeID[all_post_id[i]]
                idx = all_dist_idx[i]

                if pre_type_id != last_pre or post_type_id != last_post:
                    last_pre = pre_type_id
                    last_post = post_type_id
                    bin_count = self.dend_position_bin[(pre_type_id, post_type_id)]

                # +8 for cube corners, and +1 for full
                bin_count[idx] += syn_increment

            start_idx += chunk_size
            end_idx += chunk_size
            end_idx = min(end_idx, n_synapses)

        t_b = timeit.default_timer()

        print("Created distance histogram (optimised) in " + str(t_b - t_a) + " seconds")

    ############################################################################

    def plot_synapse_cum_dist_summary(self, pair_list):

        matplotlib.rcParams.update({'font.size': 22})
        plt.figure()

        fig_name = "SynapseCumDistSummary"
        legend_text = []

        for pair in pair_list:

            try:
                pair_id = (self.allTypes.index(pair[0]),
                           self.allTypes.index(pair[1]))
            except:
                print("Missing pair: " + str(pair))
                continue

            if pair_id not in self.dend_position_bin:
                print("Missing cum dist information for " + str(pair))
                continue

            if sum(self.dend_position_bin[pair_id]) == 0:
                print("Empty cum dist data for " + str(pair))
                continue

            cum_dist = np.cumsum(self.dend_position_bin[pair_id]) / np.sum(self.dend_position_bin[pair_id])

            # Select range to plot
            end_idx = np.where(self.dend_position_edges <= 400e-6)[0][-1]

            try:
                pre_type = pair[0]
                post_type = pair[1]

                if pre_type in self.neuron_colors:
                    plot_col = self.neuron_colors[pre_type]
                else:
                    plot_col = self.neuron_colors["default"]

                # Hack: dSPN and iSPN are overlapping, need to make dSPN visible
                if pre_type == "dSPN":
                    line_width = 7
                else:
                    line_width = 3

                plt.plot(self.dend_position_edges[:end_idx] * 1e6,
                         cum_dist[:end_idx],
                         linewidth=line_width, color=plot_col)

                fig_name += "_" + pre_type + "-" + post_type

                legend_text.append(self.neuron_name(pre_type) + "-" + self.neuron_name(post_type))

            except:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)
                import pdb
                pdb.set_trace()

        font_p = matplotlib.font_manager.FontProperties()
        font_p.set_size('small')
        plt.legend(legend_text, prop=font_p)
        plt.xlabel('Distance from soma ($\mu$m)')
        plt.ylabel('Cumulative distrib.')

        fig_name += ".pdf"

        self.save_figure(plt, fig_name)

    ############################################################################

    def plot_synapse_cum_dist(self):

        # import pdb
        # pdb.set_trace()

        for pair in self.dend_position_bin:

            if sum(self.dend_position_bin[pair]) == 0:
                # Nothing in the bins, skip this plot
                continue

            # import pdb
            # pdb.set_trace()

            cum_dist = np.cumsum(self.dend_position_bin[pair]) / np.sum(self.dend_position_bin[pair])

            idx = np.where(cum_dist < 0.5)[0][-1]
            print(self.allTypes[pair[0]] + " to " + self.allTypes[pair[1]]
                  + " 50% of synapses are within "
                  + str(self.dend_position_edges[idx] * 1e6) + " micrometer")

            # Dont plot the full range
            end_idx = np.where(self.dend_position_edges <= 300e-6)[0][-1]

            try:
                pre_type = pair[0]
                post_type = pair[1]

                matplotlib.rcParams.update({'font.size': 22})

                plt.figure()
                plt.plot(self.dend_position_edges[:end_idx] * 1e6,
                         cum_dist[:end_idx],
                         linewidth=3)
                plt.xlabel('Distance from soma ($\mu$m)')
                plt.ylabel('Cumulative distrib.')
                plt.title('Synapses ' + self.neuron_name(self.allTypes[pre_type])
                          + " to " + self.neuron_name(self.allTypes[post_type]))
                plt.tight_layout()

                plt.ion()

                if self.show_plots:
                    plt.show()

                plt.draw()
                plt.pause(0.0001)
                fig_name = "SynapseCumulativeDistribution-" \
                           + self.neuron_name(self.allTypes[pre_type]) + "-to-" \
                           + self.neuron_name(self.allTypes[post_type])

                self.save_figure(plt, fig_name)

            except Exception:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)

                print("fucked up")
                import pdb
                pdb.set_trace()

    ############################################################################

    def plot_synapse_dist(self, density_flag=False, side_len=None):

        if side_len is None:
            side_len = self.side_len

        # import pdb
        # pdb.set_trace()

        dend_hist_tot = self.dendrite_density(num_bins=len(self.dend_position_edges),
                                              bin_width=self.dend_position_edges[1],
                                              side_len=side_len)

        # self.dendPositionBin calculated by synapseDist(), done at time of init

        for pair in self.dend_position_bin:

            if sum(self.dend_position_bin[pair]) == 0:
                # Nothing in the bins, skip this plot
                continue

            end_idx = np.where(self.dend_position_edges <= 400e-6)[0][-1]

            try:
                pre_type = self.allTypes[pair[0]]
                post_type = self.allTypes[pair[1]]

                plt.rcParams.update({'font.size': 16})

                plt.figure()
                if density_flag:

                    # Skip first bin if we plot densities, since synapses from soma
                    # are in first bin, but contribute no dendritic length

                    plt.plot(self.dend_position_edges[1:end_idx] * 1e6,
                             np.divide(self.dend_position_bin[pair][1:end_idx],
                                       dend_hist_tot[post_type][1:end_idx] * 1e6))
                    plt.ylabel('Synapse/micrometer')
                    plt.xlabel('Distance from soma ($\mu$m)')

                    plt.title('Synapse density ' + self.neuron_name(pre_type)
                              + " to " + self.neuron_name(post_type))

                else:
                    plt.plot(self.dend_position_edges[:end_idx] * 1e6,
                             self.dend_position_bin[pair][:end_idx])
                    plt.ylabel('Synapse count')
                    plt.xlabel('Distance from soma ($\mu$m)')
                    plt.ylim([0, np.ceil(np.max(self.dend_position_bin[pair][:end_idx]))])

                    plt.title('Synapses ' + self.neuron_name(pre_type)
                              + " to " + self.neuron_name(post_type))

                plt.tight_layout()

                plt.ion()

                plt.draw()
                plt.pause(0.0001)

                if density_flag:
                    fig_name = "SynapseDistribution-density-dend-" \
                              + pre_type + "-to-" + post_type
                else:
                    fig_name = "SynapseDistribution-dend-" \
                              + pre_type + "-to-" + post_type

                self.save_figure(plt, fig_name)

                if self.show_plots:
                    plt.show()

            except Exception:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)

                print("fucked up")
                import pdb
                pdb.set_trace()

    ############################################################################

    ############################################################################

    # We want to get the number of synapses per unit length
    # We have the number of synapses at different distances from some (arc length)
    # So we want to calculate dendritic length, and then divide the synapses by
    # that

    def dendrite_density(self, num_bins, bin_width, side_len=None, volume_id="Striatum"):

        if side_len is None:
            side_len = self.side_len

        print("Using nBins = " + str(num_bins) + ", binWidth = " + str(bin_width))

        # 1. Loop through all morphologies
        morph_files = set([self.data["morphologies"][name]["location"][()]
                          for name in self.data["morphologies"]])

        morph_hist = dict([])
        for m_file in morph_files:
            morph_hist[m_file] = self._dend_density(m_file, num_bins, bin_width)

        neuron_types = set([m.split("_")[0]
                           for m in self.data["morphologies"]])

        # 2. Sum all the histograms together, with the right multiplicity

        dend_hist_tot = dict([])
        for nt in neuron_types:
            dend_hist_tot[nt] = np.zeros((num_bins,))

        # cornerID = self.cornerNeurons(sideLen=sideLen)
        corner_id = self.get_sub_pop(volume_type=self.volume_type,
                                     volume_part="corner",
                                     side_len=side_len,
                                     neuron_id=None,
                                     volume_id=volume_id)

        is_corner = np.zeros((self.num_neurons,), dtype=bool)
        is_corner[corner_id] = 1

        # OBS, neuronID unique, but many neurons can have same
        # morphology and properties and thus same neuronName
        # So the same name can appear twice in this list
        for nID, name in zip(self.data["network"]["neurons"]["neuronID"],
                             self.data["network"]["neurons"]["name"]):

            if not is_corner[nID]:
                # Only post synaptic corner neurons included
                continue

            m_file = self.data["morphologies"][name]["location"][()]
            n_type = name.decode().split("_")[0]
            dend_hist_tot[n_type] += morph_hist[m_file]

        return dend_hist_tot

    ############################################################################

    def _dend_density(self, swc_file, num_bins, bin_width):

        from .neuron_morphology import NeuronMorphology

        dend_hist = np.zeros((num_bins,))

        n_morph = NeuronMorphology(swc_filename=swc_file)

        print("Parsing dendrite histogram : " + swc_file)

        # 0,1,2: x,y,z  3: radie, 4: dist to soma
        dist_to_soma = n_morph.dend[:, 4]
        compartment_length = np.zeros((n_morph.dend.shape[0],))

        for idx, link in enumerate(n_morph.dend_links):
            compartment_length[idx] = np.linalg.norm(n_morph.dend[link[0], :3]
                                                     - n_morph.dend[link[1], :3])
            dist_to_soma[idx] = (n_morph.dend[link[0], 4] + n_morph.dend[link[1], 4]) / 2

        for d, cl in zip(dist_to_soma, compartment_length):
            idx = int(d / bin_width)
            dend_hist[idx] += cl

        return dend_hist

    ############################################################################

    # This

    def virtual_axon_synapses(self, post_neuron_type):

        virt_idx = np.where([n["virtualNeuron"] for n in self.network["neurons"]])[0]

        post_idx = self.populations[post_neuron_type]

        virt_syn = dict([])

        for v_idx in virt_idx:
            n_type = self.network["neurons"][v_idx]["name"].split("_")[0]

            syn_mat = self.connection_matrix[v_idx, :][:, post_idx]

            if syn_mat.count_nonzero() == 0:
                # No synapses here, skip plot
                continue

            syn_mat = syn_mat.todense()
            syn_mat = syn_mat[np.where(syn_mat > 0)]

            if n_type not in virt_syn:
                virt_syn[n_type] = [syn_mat]
            else:
                virt_syn[n_type].append(syn_mat)

        for axonType in virt_syn:
            plt.figure()

            v_syn = np.concatenate(virt_syn[axonType]).flatten().transpose()

            # We dont want to show the count for zeros
            n_bins = np.max(v_syn) + 1
            virt_syn_bins = np.array(range(1, n_bins))

            plt.hist(v_syn, bins=virt_syn_bins, align="left")
            plt.xlabel("Synapses (only connected pairs)")
            plt.ylabel("Count")
            plt.title("Synapses from " + self.neuron_name(axonType)
                      + " to " + self.neuron_name(post_neuron_type)
                      + " (nSyn=" + str(np.sum(v_syn)) + ")")

            plt.tight_layout()

            plt.ion()
            plt.draw()

            fig_name = "VirtuaAxon-synapses-" + axonType + "-to-" + post_neuron_type

            self.save_figure(plt, fig_name)

            if self.show_plots:
                plt.show()

            # import pdb
            # pdb.set_trace()

    ############################################################################

    def count_motifs(self, type_a, type_b, type_c, n_repeats=1000000):

        print("Counting motivs between " + str(type_a) + ", " + str(type_b) \
              + " " + str(type_c) + ". " + str(n_repeats) + " repeats.")

        ida = self.populations[type_a]
        idb = self.populations[type_b]
        idc = self.populations[type_c]

        # Init counter
        motif_ctr = np.zeros((64,), dtype=int)

        # bit 1 : A->B (1)
        # bit 2 : A<-B (2)
        # bit 3 : A->C (4)
        # bit 4 : A<-C (8)
        # bit 5 : B->C (16)
        # bit 6 : B<-C (32)

        # Faster to generate all int at once
        i_a_all = np.random.randint(len(ida), size=(n_repeats,))
        i_b_all = np.random.randint(len(idb), size=(n_repeats,))
        i_c_all = np.random.randint(len(idc), size=(n_repeats,))

        for i_rep in range(0, n_repeats):

            if i_rep % 100000 == 0 and i_rep > 0:
                print("rep: " + str(i_rep))

            i_a = i_a_all[i_rep]
            i_b = i_b_all[i_rep]
            i_c = i_c_all[i_rep]

            # In case the same neuron was picked twice, redo sampling
            while i_a == i_b or i_b == i_c or i_c == i_a:
                i_a = np.random.randint(len(ida))
                i_b = np.random.randint(len(idb))
                i_c = np.random.randint(len(idc))

            idx = (int(self.connection_matrix[i_a, i_b] > 0) * 1
                   + int(self.connection_matrix[i_b, i_a] > 0) * 2
                   + int(self.connection_matrix[i_a, i_c] > 0) * 4
                   + int(self.connection_matrix[i_c, i_a] > 0) * 8
                   + int(self.connection_matrix[i_b, i_c] > 0) * 16
                   + int(self.connection_matrix[i_c, i_b] > 0) * 32)

            motif_ctr[idx] += 1

        return motif_ctr, type_a, type_b, type_c

    ############################################################################

    def analyse_single_motifs(self, neuron_type, num_repeats=10000000):

        (motif_ctr, tA, tB, tC) = self.count_motifs(type_a=neuron_type,
                                                    type_b=neuron_type,
                                                    type_c=neuron_type,
                                                    n_repeats=num_repeats)

        # !!! For debug
        # motifCtr = [bin(x).count('1') for x in np.arange(0,64)]

        # No connections
        no_con = motif_ctr[0]

        # One connection
        one_con = (motif_ctr[1] + motif_ctr[2] + motif_ctr[4] + motif_ctr[8]
                   + motif_ctr[16] + motif_ctr[32])

        # Two connections
        two_con_diverge = motif_ctr[1 + 4] + motif_ctr[2 + 16] + motif_ctr[8 + 32]
        two_con_converge = motif_ctr[2 + 8] + motif_ctr[1 + 32] + motif_ctr[4 + 16]
        two_con_line = (motif_ctr[1 + 16] + motif_ctr[2 + 4] + motif_ctr[4 + 32]
                        + motif_ctr[8 + 1] + motif_ctr[16 + 8] + motif_ctr[32 + 2])
        two_con_bi = motif_ctr[1 + 2] + motif_ctr[4 + 8] + motif_ctr[16 + 32]

        # Three connections
        three_circle = motif_ctr[1 + 16 + 8] + motif_ctr[2 + 4 + 32]
        three_circle_flip = (motif_ctr[2 + 16 + 8] + motif_ctr[1 + 32 + 8] + motif_ctr[1 + 16 + 4]
                             + motif_ctr[1 + 4 + 32] + motif_ctr[2 + 8 + 32] + motif_ctr[2 + 4 + 16])
        three_bi_div = (motif_ctr[1 + 2 + 4] + motif_ctr[1 + 2 + 16]
                        + motif_ctr[4 + 8 + 1] + motif_ctr[4 + 8 + 32]
                        + motif_ctr[16 + 32 + 2] + motif_ctr[16 + 32 + 8])
        three_bi_conv = (motif_ctr[1 + 2 + 8] + motif_ctr[1 + 2 + 32]
                         + motif_ctr[4 + 8 + 2] + motif_ctr[4 + 8 + 16]
                         + motif_ctr[16 + 32 + 1] + motif_ctr[16 + 32 + 4])

        # Four connections
        four_double_bi = motif_ctr[1 + 2 + 4 + 8] + motif_ctr[1 + 2 + 16 + 32] + motif_ctr[4 + 8 + 16 + 32]
        four_bi_diverge = motif_ctr[1 + 2 + 4 + 16] + motif_ctr[4 + 8 + 1 + 32] + motif_ctr[16 + 32 + 2 + 8]
        four_bi_converge = motif_ctr[1 + 2 + 8 + 32] + motif_ctr[4 + 8 + 2 + 16] + motif_ctr[16 + 32 + 1 + 4]
        four_bi_cycle = (motif_ctr[1 + 2 + 4 + 32] + motif_ctr[1 + 2 + 16 + 8]
                         + motif_ctr[4 + 8 + 1 + 16] + motif_ctr[4 + 8 + 32 + 2]
                         + motif_ctr[16 + 32 + 2 + 4] + motif_ctr[16 + 32 + 8 + 1])

        # Five connections
        five_con = (motif_ctr[2 + 4 + 8 + 16 + 32]
                    + motif_ctr[1 + 4 + 8 + 16 + 32]
                    + motif_ctr[1 + 2 + 8 + 16 + 32]
                    + motif_ctr[1 + 2 + 4 + 16 + 32]
                    + motif_ctr[1 + 2 + 4 + 8 + 32]
                    + motif_ctr[1 + 2 + 4 + 8 + 16])

        # Six connections
        six_con = motif_ctr[1 + 2 + 4 + 8 + 16 + 32]

        con_data = [("No connection", no_con),
                    ("One connection", one_con),
                    ("Two connections, diverge", two_con_diverge),
                    ("Two connections, converge", two_con_converge),
                    ("Two connections, line", two_con_line),
                    ("Two connections, bidirectional", two_con_bi),
                    ("Three connections, circle", three_circle),
                    ("Three connections, circular, one flipped", three_circle_flip),
                    ("Three connections, one bidirectional, one away", three_bi_div),
                    ("Three connections, one bidirectional, one towards",
                    three_bi_conv),
                    ("Four connections, two bidirectional", four_double_bi),
                    ("Four connections, one bi, two away", four_bi_diverge),
                    ("Four connections, one bi, two towards", four_bi_converge),
                    ("Four connections, cycle with one bidirectional", four_bi_cycle),
                    ("Five connections", five_con),
                    ("Six connections", six_con)]

        print("Motif analysis:")
        for (name, data) in con_data:
            print(name + " " + str(100 * data / num_repeats) + "% (" + str(data) + ")")

    ############################################################################

    # If A is connected to B and C, what is probability that
    # B and C are connected?

    def simple_motif(self, type_a, type_b, type_c, n_rep=1000):

        ida = self.populations[type_a]
        idb = self.populations[type_b]
        idc = self.populations[type_c]

        # We need to reduce the IDA population, otherwise it will be too slow
        if n_rep < len(ida):
            try:
                r_perm_idx = np.random.permutation(len(ida))[:n_rep]
                ida = [ida[x] for x in r_perm_idx]
            except:
                import traceback
                tstr = traceback.format_exc()
                print(tstr)
                import pdb
                pdb.set_trace()

        con_bc = 0
        con_cb = 0
        con_bi = 0
        all_ctr = 0

        n_a = len(ida)
        ctr_a = 0

        for iA in ida:

            ctr_a += 1

            if ctr_a % 100 == 0:
                print(str(ctr_a) + "/" + str(n_a))

            mat_row = self.connection_matrix[iA, :]

            i_b_all = idb[mat_row[0, idb].nonzero()[1]]
            i_c_all = idc[mat_row[0, idc].nonzero()[1]]

            for iB in i_b_all:
                for iC in i_c_all:
                    bc = self.connection_matrix[iB, iC]
                    cb = self.connection_matrix[iC, iB]

                    if bc > 0:
                        con_bc += 1
                        if cb > 0:
                            con_cb += 1
                            con_bi += 1
                    elif cb > 0:
                        con_cb += 1

                    all_ctr += 1

        print("If " + str(type_a) + " connected to " + str(type_b) + " " \
              + str(type_c) + ":")
        print("P(" + type_b + "->" + type_c + ") = " + str((100.0 * con_bc) / all_ctr) + "%")
        print("P(" + type_c + "->" + type_b + ") = " + str((100.0 * con_cb) / all_ctr) + "%")
        print("P(" + type_b + "<->" + type_c + ") = " + str((100.0 * con_bi) / all_ctr) + "%")

        return con_bc, con_cb, con_bi, all_ctr

    ############################################################################

    # Pick a post synaptic neuron, find out the distance to its closest
    # presynaptic neighbour

    def nearest_pre_neighbour_distance(self, pre_type, post_type, rabies_rate=1.0,
                                       name_str=""):

        if pre_type not in self.populations:
            print("nearestPreNeighbourDistance: " + str(pre_type)
                  + " is not in the simulation")
            return

        if post_type not in self.populations:
            print("nearestPreNeighbourDistance: " + str(post_type)
                  + " is not in the simulation")
            return

        # postPop = self.getSubPop(neuronID=self.populations[postType],
        #                         volumePart="centre")

        post_pop = self.populations[post_type]
        pre_pop = self.populations[pre_type]

        # Assume 300 micrometer slice
        max_z = np.max(self.positions[post_pop, 2])
        min_z = np.max(self.positions[post_pop, 2])

        slice_min_z = (max_z + min_z) / 2 - 150e-6
        slice_max_z = (max_z + min_z) / 2 + 150e-6

        slice_pop = np.where(np.bitwise_and(slice_min_z <= self.positions[:, 2],
                                           self.positions[:, 2] <= slice_max_z))[0]

        # Only look at neurons in slice
        pre_pop = np.intersect1d(pre_pop, slice_pop)
        post_pop = np.intersect1d(post_pop, slice_pop)

        # conIdx = np.sum(self.connectionMatrix[prePop,postPop],axis=1)

        nearest_dist = np.nan * np.zeros((len(post_pop),))

        for c, postID in enumerate(post_pop):

            # Calculate distance to all connected neighbours
            idx = np.where(self.connection_matrix[pre_pop, postID].todense() > 0)[0]
            con_idx = pre_pop[idx]

            if rabies_rate is not None and rabies_rate < 1:
                keep_idx = np.where(np.random.rand(len(con_idx)) < rabies_rate)[0]
                con_idx = con_idx[keep_idx]

            d = np.sqrt(np.sum((self.positions[con_idx, :3]
                                - self.positions[postID, :3]) ** 2,
                               axis=1))

            if len(d) > 0:
                nearest_dist[c] = np.min(d) * 1e6  # micrometers

        max_dist = np.nanmax(nearest_dist)

        plt.figure()
        matplotlib.rcParams.update({"font.size": 22})

        plt.hist(nearest_dist, np.arange(0, max_dist + 25, 25))

        plt.xlabel("Distance")
        if rabies_rate is not None and rabies_rate < 1:
            plt.ylabel("Count (Rabies rate: " + str(rabies_rate * 100) + "%)")
        else:
            plt.ylabel("Count")
        plt.title("Nearest presynaptic neighbour "
                  + self.neuron_name(pre_type) + " to " + self.neuron_name(post_type))

        # Data from Sabatini 2016
        if pre_type == "LTS" and (post_type == "dSPN" or post_type == "iSPN"):
            sabatini_lts = np.genfromtxt("DATA/LTS-nearest-neighbour-points-Sabatini2016.csv", delimiter=",")
            lts_points = sabatini_lts[:, 1] * 1e3  # Get in micrometer
            plt.hist(lts_points, color='r', histtype="step")

        plt.ion()
        plt.draw()

        fig_name = "Nearest-presynaptic-slice-neighbour-to-" \
                  + str(post_type) + "-from-" + str(pre_type) + "-ID-" \
                  + str(self.network["SlurmID"]) + name_str + ".pdf"

        self.save_figure(plt, fig_name)

        if self.show_plots:
            plt.show()

        plt.pause(0.001)

        if self.close_plots:
            time.sleep(1)
            plt.close()

    ############################################################################

    # !!! Doh, just realised using the connection matrix might have been smarter

    def nearest_pre_neighbour_distance_old(self, pre_type, post_type, num_max=1000,
                                           name_str=""):

        post_pop = self.get_sub_pop(neuron_id=self.populations[post_type],
                                    volume_part="centre")

        if len(post_pop) > num_max:
            post_pop = np.random.permutation(post_pop)[1:num_max]

        nearest_neighbour_dist = dict([])
        for idx in post_pop:
            nearest_neighbour_dist[idx] = np.inf

        potential_pre_pop = self.populations[pre_type]

        last_syn_pair = (np.nan, np.nan)

        for synapses in self.network_load.synapse_iterator():
            for synRow in synapses:
                pre_idx = synRow[0]
                post_idx = synRow[1]
                if (pre_idx, post_idx) == last_syn_pair:
                    # Already did this pair
                    continue

                if post_idx in post_pop and pre_idx in potential_pre_pop:
                    d = np.sqrt(np.sum((self.positions[pre_idx, :]
                                        - self.positions[post_idx, :]) ** 2))
                    if d < nearest_neighbour_dist[post_idx]:
                        nearest_neighbour_dist[post_idx] = d

                last_syn_pair = (pre_idx, post_idx)

        nn_dist = np.array(list(nearest_neighbour_dist.values()))
        nn_dist = nn_dist[nn_dist < np.inf] * 1e6
        max_dist = max(nn_dist)

        plt.figure()
        matplotlib.rcParams.update({"font.size": 22})

        plt.hist(nn_dist, np.arange(0, max_dist + 25, 25))

        plt.xlabel("Distance")
        plt.ylabel("Count")
        plt.title("Nearest presynaptic neighbour "
                  + self.neuron_name(pre_type) + " to " + self.neuron_name(post_type))

        # Data from Sabatini 2016
        if pre_type == "LTS" and (post_type == "dSPN" or post_type == "iSPN"):
            sabatini_lts = np.genfromtxt("../DATA/LTS-nearest-neighbour-points-Sabatini2016.csv", delimiter=",")
            lts_points = sabatini_lts[:, 1] * 1e3  # Get in micrometer
            plt.hist(lts_points, color='r', histtype="step")

        plt.ion()
        plt.draw()

        fig_name = "figures/Nearest-presynaptic-neighbour-to-" \
                   + str(post_type) + "-from-" + str(pre_type)

        self.save_figure(plt, fig_name)

        if self.show_plots:
            plt.show()

        plt.pause(0.001)

    ############################################################################

    # Inspired by:
    # Nao Chuhma, Kenji F. Tanaka, Rene Hen and Stephen Rayport 2011
    #
    # 10% of a neuron type are marked, fraction of presynaptic neurons
    # out of total population

    def chuhma_virtual_experiment(self, tagged_type=["dSPN", "iSPN"], tag_fraction=0.1):

        print("Doing Chuma experiments: " + str(tagged_type))

        idx = np.concatenate([self.populations[x] for x in tagged_type])
        num_tagged = np.round(len(idx) * tag_fraction).astype(int)

        tagged_neurons = np.sort(np.random.permutation(idx)[:num_tagged])

        # Find all presynaptic neurons
        pre_idx = np.where(np.sum(self.connection_matrix[:, tagged_neurons],
                                 axis=1) > 0)[0]

        # Ooops, no fun... pretty much all neurons got tagged.

        import pdb
        pdb.set_trace()

    ############################################################################

    # Number of ChINs connected to each MS

    ############################################################################


if __name__ == "__main__":

    assert False, "Do you want to run Network_analyse_striatum.py instead?"

    if len(sys.argv) > 1:
        hdf5File = sys.argv[1]
        print("Using user supplied HDF5 file: " + hdf5File)

    else:
        hdf5File = None  # Auto detect, use latest file
        # hdf5File = "save/network-connect-synapse-file-0.hdf5"

    # import cProfile
    # cProfile.run('na = NetworkAnalyse(hdf5File,loadCache=False,lowMemory=False)')
    # import pdb
    # pdb.set_trace()

    na = NetworkAnalyse(hdf5File, loadCache=True, lowMemory=False, sideLen=250e-6,
                        volumeType="full")  # "cube"

    # na = NetworkAnalyse(hdf5File,loadCache=False)

    # na.plotNeurons(neuronID=5,plotPreNeurons=False)
    # na.plotNeurons(neuronID=5)

    enableAnalysis = True  # True #False

    # No exp data for this
    #  dist_3d = False
    #  na.plotConnectionProbabilityChannels("FSN","FSN", \
    #                                       dist_3d=dist_3d, \
    #                                       expMaxDist=[],\
    #                                       expData=[])
    #
    #  import pdb
    #  pdb.set_trace()

    # na.cornerNeurons(sideLen = 100e-6)

    # na.plotSynapseDist(densityFlag=True)
    # na.plotSynapseCumDist()

    # na.virtualAxonSynapses("dSPN")
    # na.virtualAxonSynapses("iSPN")
    # na.virtualAxonSynapses("FSN")

    # na.simpleMotif("dSPN","dSPN","dSPN")
    # na.simpleMotif("iSPN","iSPN","iSPN")
    # na.simpleMotif("dSPN","dSPN","iSPN")
    # na.analyseSingleMotifs("dSPN")

    na.nearest_pre_neighbour_distance("LTS", "dSPN")
    na.nearest_pre_neighbour_distance("LTS", "iSPN")

    # na.ChuhmaVirtualExperiment(taggedType=["dSPN","iSPN"])

    dist_3d = False

    # 3/21 LTS->MS, Basal Ganglia book --- distance??
    # Ibanez-Sandoval, ..., Tepper  2011 3/21 -- if patching around visual axon
    # but 2/60 when patching blind
    # !!! Use the focused 3/21 statistics for validation!! --- please :)
    na.plot_connection_probability("LTS", "dSPN",
                                   dist_3d=dist_3d,
                                   exp_max_dist=[250e-6],
                                   exp_data=[2 / 60.0],
                                   exp_data_detailed=[(2, 60)])

    na.plot_connection_probability("LTS", "iSPN",
                                   dist_3d=dist_3d,
                                   exp_max_dist=[250e-6],
                                   exp_data=[2 / 60.0],
                                   exp_data_detailed=[(2, 60)])

    # Silberberg et al 2013, 2/12 FS-> LTS connected --- distance??
    na.plot_connection_probability("FSN", "LTS",
                                   dist_3d=dist_3d,
                                   exp_max_dist=[250e-6],
                                   exp_data=[2.0 / 12],
                                   exp_data_detailed=[(2, 12)])

    # import pdb
    # pdb.set_trace()

    na.plot_connection_probability("dSPN", "ChIN",
                                   dist_3d=dist_3d)
    na.plot_connection_probability("iSPN", "ChIN",
                                   dist_3d=dist_3d)

    na.plot_connection_probability("ChIN", "LTS",
                                   dist_3d=dist_3d)
    na.plot_connection_probability("ChIN", "iSPN",
                                   dist_3d=dist_3d)
    na.plot_connection_probability("ChIN", "dSPN",
                                   dist_3d=dist_3d)

    print("Check the ChIN stuff")
    # import pdb
    # pdb.set_trace()

    # import pdb
    # pdb.set_trace()

    na.nearest_pre_neighbour_distance("FSN", "dSPN")
    na.nearest_pre_neighbour_distance("FSN", "iSPN")

    na.plot_num_synapses_per_pair("ChIN", "dSPN")
    na.plot_num_synapses_per_pair("ChIN", "iSPN")

    na.plot_num_synapses_per_pair("dSPN", "ChIN")
    na.plot_num_synapses_per_pair("iSPN", "ChIN")

    # 2-5 ChIN should connect to each MS (approx)
    na.plot_incoming_connections(neuron_type="dSPN", pre_type="ChIN")
    na.plot_incoming_connections(neuron_type="iSPN", pre_type="ChIN")

    na.plot_incoming_connections(neuron_type="ChIN", pre_type="dSPN")
    na.plot_incoming_connections(neuron_type="ChIN", pre_type="iSPN")

    # LTS plots
    na.plot_num_synapses_per_pair("LTS", "dSPN")
    na.plot_num_synapses_per_pair("LTS", "iSPN")
    na.plot_num_synapses_per_pair("LTS", "ChIN")

    na.plot_num_synapses_per_pair("ChIN", "LTS")
    na.plot_num_synapses_per_pair("FSN", "LTS")

    na.plot_incoming_connections(neuron_type="dSPN", pre_type="LTS")
    na.plot_incoming_connections(neuron_type="iSPN", pre_type="LTS")
    na.plot_incoming_connections(neuron_type="ChIN", pre_type="LTS")

    na.plot_incoming_connections(neuron_type="LTS", pre_type="ChIN")
    na.plot_incoming_connections(neuron_type="LTS", pre_type="FSN")

    if True or enableAnalysis:
        dist_3d = False

        # 100e-6 from Planert 2010, and 250e-6 data from Gittis 2010
        # 150e-6 from Gittis 2011 (actually 100 +/- 50 micrometers)
        na.plot_connection_probability("FSN", "iSPN",
                                       dist_3d=dist_3d,
                                       exp_max_dist=[100e-6, 150e-6, 250e-6],
                                       exp_data=[6 / 9.0, 21 / 54.0, 27 / 77.0],
                                       exp_data_detailed=[(6, 9), (21, 54), (27, 77)])
        na.plot_connection_probability("FSN", "dSPN",
                                       dist_3d=dist_3d,
                                       exp_max_dist=[100e-6, 150e-6, 250e-6],
                                       exp_data=[8 / 9.0, 29 / 48.0, 48 / 90.0],
                                       exp_data_detailed=[(8, 9), (29, 48), (48, 90)])
        na.plot_connection_probability("dSPN", "iSPN",
                                       dist_3d=dist_3d,
                                       exp_max_dist=[50e-6, 100e-6],
                                       exp_data=[3 / 47.0, 3 / 66.0],
                                       exp_data_detailed=[(3, 47), (3, 66)])
        na.plot_connection_probability("dSPN", "dSPN",
                                       dist_3d=dist_3d,
                                       exp_max_dist=[50e-6, 100e-6],
                                       exp_data=[5 / 19.0, 3 / 43.0],
                                       exp_data_detailed=[(5, 19), (3, 43)])
        na.plot_connection_probability("iSPN", "dSPN",
                                       dist_3d=dist_3d,
                                       exp_max_dist=[50e-6, 100e-6],
                                       exp_data=[13 / 47.0, 10 / 80.0],
                                       exp_data_detailed=[(13, 47), (10, 80)])
        na.plot_connection_probability("iSPN", "iSPN",
                                       dist_3d=dist_3d,
                                       exp_max_dist=[50e-6, 100e-6],
                                       exp_data=[14 / 39.0, 7 / 31.0],
                                       exp_data_detailed=[(14, 39), (7, 31)])

        # No exp data for this -- Gittis,...,Kreitzer 2010 (p2228) -- 7/12 (and 3/4 reciprocal) -- distance?
        # FS->FS synapses weaker, 1.1 +/- 1.5nS
        na.plot_connection_probability("FSN", "FSN",
                                       dist_3d=dist_3d,
                                       exp_max_dist=[250e-6],
                                       exp_data=[7 / 12.0],
                                       exp_data_detailed=[(7, 12)])

        # Do we have ChINs?
        # ChIN data, Johanna had ref. ????
        # Janickova

        if True:
            # REF???!?!?!?!
            # na.plotConnectionProbability("ChIN","iSPN", \
            #                                     dist_3d=dist_3d,
            #                                     expMaxDist=[200e-6],
            #                                     expData=[62/89.0],
            #                                     expDataDetailed=[(62,89)])
            # na.plotConnectionProbability("ChIN","dSPN", \
            #                                     dist_3d=dist_3d,
            #                                     expMaxDist=[200e-6],
            #                                     expData=[62/89.0],
            #                                     expDataDetailed=[(62,89)])

            # Derived from Janickova H, ..., Bernard V 2017
            na.plot_connection_probability("ChIN", "iSPN",
                                           dist_3d=dist_3d,
                                           exp_max_dist=[200e-6],
                                           exp_data=[0.05])
            na.plot_connection_probability("ChIN", "dSPN",
                                           dist_3d=dist_3d,
                                           exp_max_dist=[200e-6],
                                           exp_data=[0.05])

            na.plot_connection_probability("ChIN", "FSN",
                                           dist_3d=dist_3d)

            na.plot_incoming_connections(neuron_type="dSPN", pre_type="ChIN")
            na.plot_incoming_connections(neuron_type="iSPN", pre_type="ChIN")

        if True:
            na.plot_connection_probability("LTS", "ChIN",
                                           dist_3d=dist_3d)

            na.plot_connection_probability("ChIN", "LTS",
                                           dist_3d=dist_3d)

    if True:
        print("The synapse dist function needs a density func, which currently not working since we no longer "
              + "include compartment length in the dend data, so need to calculate it")
        na.plot_synapse_dist()
        na.plot_synapse_cum_dist()
        na.plot_synapse_dist(density_flag=True)

        # import pdb
        # pdb.set_trace()

    if True and enableAnalysis:
        na.plot_num_synapses_per_pair("FSN", "dSPN")
        na.plot_num_synapses_per_pair("FSN", "iSPN")
        na.plot_num_synapses_per_pair("dSPN", "dSPN")
        na.plot_num_synapses_per_pair("dSPN", "iSPN")
        na.plot_num_synapses_per_pair("iSPN", "dSPN")
        na.plot_num_synapses_per_pair("iSPN", "iSPN")

    if True and enableAnalysis:
        na.plot_incoming_connections(neuron_type="dSPN", pre_type="iSPN")
        na.plot_incoming_connections(neuron_type="dSPN", pre_type="dSPN")
        na.plot_incoming_connections(neuron_type="dSPN", pre_type="FSN")

    if True and enableAnalysis:
        na.plot_incoming_connections(neuron_type="iSPN", pre_type="iSPN")
        na.plot_incoming_connections(neuron_type="iSPN", pre_type="dSPN")
        na.plot_incoming_connections(neuron_type="iSPN", pre_type="FSN")

    print("All done, exiting")
    # import pdb
    # pdb.set_trace()

