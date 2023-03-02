import os
import numpy as np

from snudda import SnuddaPlace
from snudda import SnuddaDetect
from snudda.optimise.optimise_pruning import OptimisePruning
from snudda.analyse.analyse import SnuddaAnalyse

con_types = [('FS', 'FS', 4),
             ('FS', 'dSPN', 4),
             ('FS', 'iSPN', 4),
             ('dSPN', 'dSPN', 3),
             ('dSPN', 'iSPN', 3),
             ('iSPN', 'iSPN', 3)]


n_neurons = 150
con_type = "GABA"

SPN2SPNdistDepPruning = "1-exp(-(0.4*d/60e-6)**2)"
FS_dist_dep_pruning = "exp(-(0.5*d/60e-6)**2)"

all_experimental_data = dict()
all_experimental_data["dSPN", "iSPN"] = [(0, 50e-6, 3 / 47.0), (0, 100e-6, 3 / 66.0)]
all_experimental_data["dSPN", "dSPN"] = [(0, 50e-6, 5 / 19.0), (0, 100e-6, 3 / 43.0)]
all_experimental_data["iSPN", "iSPN"] = [(0, 50e-6, 14 / 39.0), (0, 100e-6, 7 / 31.0)]
all_experimental_data["iSPN", "dSPN"] = [(0, 50e-6, 13 / 47.0), (0, 100e-6, 10 / 80.0)]
all_experimental_data["FS", "FS"] = [(0, 250e-6, 7 / 12.0)]
all_experimental_data["FS", "iSPN"] = [(0, 100e-6, 6 / 9.0), (0, 150e-6, 21 / 54.0), (0, 250e-6, 27 / 77.0)]
all_experimental_data["FS", "dSPN"] = [(0, 100e-6, 8 / 9.0), (0, 150e-6, 29 / 48.0), (0, 250e-6, 48 / 90.0)]

for ct in con_types:
    for dist_dep in [True, False]:
        for num_params in [1, 2, 3, 4]:

            pre_type = ct[0]
            post_type = ct[1]
            avg_num_synapses_per_pair = ct[2]

            distance_dependent_pruning = dist_dep

            extra_pruning_parameters = {}

            if distance_dependent_pruning:
                if pre_type == "FS" and (post_type == "dSPN" or post_type == "iSPN"):
                    extra_pruning_parameters = {'distPruning': FS_dist_dep_pruning}
                elif "SPN" in pre_type and "SPN" in post_type:
                    extra_pruning_parameters = {'distPruning': SPN2SPNdistDepPruning}
                else:
                    # No distance dependent pruning available for this neuron type, set it to False
                    distance_dependent_pruning = False

            dd_str = "_dd" if distance_dependent_pruning else ""
            network_path = os.path.join("networks", f"{pre_type}_to_{post_type}_np{num_params}{dd_str}")

            experimental_data = all_experimental_data[pre_type, post_type]

            cube_side = 1.5 * np.max([x[1] for x in experimental_data])

            # Place cube mesh and setup network config
            from snudda.place.create_cube_mesh import create_cube_mesh

            mesh_file = os.path.join(network_path, "cube-mesh.obj")
            create_cube_mesh(mesh_file, [0, 0, 0], cube_side)

            from snudda.init import SnuddaInit

            si = SnuddaInit(network_path=network_path, random_seed=123,
                            snudda_data="../../../../BasalGangliaData/data/")

            si.define_structure(struct_name="Cube", struct_mesh=mesh_file, d_min=12e-6, mesh_bin_width=25e-6)

            if pre_type == post_type:
                si.add_neurons(name=pre_type, num_neurons=n_neurons, volume_id="Cube",
                               neuron_dir=os.path.join("$DATA", "neurons", "striatum", "dspn"))
            else:
                si.add_neurons(name=pre_type, num_neurons=int(n_neurons / 2), volume_id="Cube",
                               neuron_dir=os.path.join("$DATA", "neurons", "striatum", "dspn"))
                si.add_neurons(name=post_type, num_neurons=int(n_neurons / 2), volume_id="Cube",
                               neuron_dir=os.path.join("$DATA", "neurons", "striatum", "dspn"))

            # The parameters here does not matter, they will be set during optimisation
            si.add_neuron_target(neuron_name=pre_type,
                                 target_name=post_type,
                                 connection_type=con_type,
                                 dist_pruning=SPN2SPNdistDepPruning,
                                 f1=None, soft_max=None, mu2=None,
                                 a3=None,
                                 conductance=[0.24e-9, 0.1e-9],
                                 mod_file="tmGabaA")

            si.write_json()

            # Place neurons
            sp = SnuddaPlace(network_path=network_path, verbose=False)
            sp.place()

            sd = SnuddaDetect(network_path=network_path, hyper_voxel_size=100)
            sd.detect()

            op = OptimisePruning(network_path=network_path)
            op.merge_putative_synapses(force_merge=True)

            print(op.prune.connectivity_distributions)
            print(op.prune.hist_file["meta/connectivityDistributions"][()])

            res = op.optimize(pre_type=pre_type, post_type=post_type, con_type=con_type,
                              experimental_data=experimental_data,
                              avg_num_synapses_per_pair=avg_num_synapses_per_pair,
                              extra_pruning_parameters=extra_pruning_parameters,
                              workers=8, maxiter=1000, tol=0.00001, num_params=num_params)

            if num_params == 1:
                param_str = f"f1 = %f" % (res.x[0])
            elif num_params == 2:
                param_str = f"f1 = %f, mu2 = %f" % (res.x[0], res.x[1])
            elif num_params == 3:
                param_str = f"f1 = %f, mu2 = %f, a3 = %f" % (res.x[0], res.x[1], res.x[2])
            elif num_params == 4:
                param_str = f"f1 = %f, softMax = %f, mu2 = %f, a3 = %f" % (res.x[0], res.x[1], res.x[2], res.x[3])
            else:
                param_str = res.x

            if "distPruning" in extra_pruning_parameters:
                param_str += f" ({extra_pruning_parameters['distPruning']})"

            print(param_str)
            print(res)

            # Get the last file
            # list_of_files = glob.glob(os.path.join(network_path, "temp", "network-synapses-*hdf5"))
            # network_file = max(list_of_files, key=os.path.getctime)

            network_file = os.path.join(network_path, "network-synapses.hdf5")

            dist3D = False
            y_max_H = None

            sa = SnuddaAnalyse(network_file)

            if pre_type == "dSPN" and post_type == "iSPN":
                sa.plot_connection_probability("dSPN", "iSPN", dist_3d=True, exp_max_dist=[50e-6, 100e-6],
                                               exp_data_detailed=[(3, 47), (3, 66)], sub_title=param_str)
                sa.plot_num_synapses_per_pair("dSPN", "iSPN", sub_title=param_str)
            elif pre_type == "dSPN" and post_type == "dSPN":
                sa.plot_connection_probability("dSPN", "dSPN", dist_3d=True, exp_max_dist=[50e-6, 100e-6],
                                               exp_data_detailed=[(5, 19), (3, 43)], sub_title=param_str)
                sa.plot_num_synapses_per_pair("dSPN", "dSPN", sub_title=param_str)
            elif pre_type == "iSPN" and post_type == "iSPN":
                sa.plot_connection_probability("iSPN", "iSPN", dist_3d=True, exp_max_dist=[50e-6, 100e-6],
                                               exp_data_detailed=[(14, 39), (7, 31)], sub_title=param_str)
                sa.plot_num_synapses_per_pair("iSPN", "iSPN", sub_title=param_str)
            elif pre_type == "iSPN" and post_type == "dSPN":
                sa.plot_connection_probability("iSPN", "dSPN", dist_3d=True, exp_max_dist=[50e-6, 100e-6],
                                               exp_data_detailed=[(13, 47), (10, 80)], sub_title=param_str)
                sa.plot_num_synapses_per_pair("iSPN", "dSPN", sub_title=param_str)
            elif pre_type == "FS" and post_type == "FS":
                sa.plot_connection_probability("FS", "FS", dist_3d=True, exp_max_dist=[250e-6],
                                               exp_data_detailed=[(7, 12)], sub_title=param_str)
                sa.plot_num_synapses_per_pair("FS", "FS", sub_title=param_str)
            elif pre_type == "FS" and post_type == "iSPN":
                sa.plot_connection_probability("FS", "iSPN", dist_3d=True, exp_max_dist=[100e-6, 150e-6, 250e-6],
                                               exp_data_detailed=[(6, 9), (21, 54), (27, 77)], sub_title=param_str)
                sa.plot_num_synapses_per_pair("FS", "iSPN", sub_title=param_str)
            elif pre_type == "FS" and post_type == "dSPN":
                sa.plot_connection_probability("FS", "dSPN", dist_3d=True, exp_max_dist=[100e-6, 150e-6, 250e-6],
                                               exp_data_detailed=[(8, 9), (29, 48), (48, 90)], sub_title=param_str)
                sa.plot_num_synapses_per_pair("FS", "dSPN", sub_title=param_str)


