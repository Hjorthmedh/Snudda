# This writes the JSON config for the network
import sys
import os

from snudda.place import create_cube_mesh
from snudda import SnuddaInit

os.environ["SNUDDA_DATA"] = "../../../../BasalGangliaData/data/"
assert os.path.isdir(os.getenv("SNUDDA_DATA")), f"You need to have BasalGangliaData installed for this example. You can run this example without it, but then do not execute this cell."

if len(sys.argv) > 1:
    network_path = sys.argv[1]
else:
    network_path = os.path.join("FS_network_1")

print(f"Using network_path = {network_path}")

mesh_file = os.path.join(network_path, "mesh", "volume.obj")
create_cube_mesh(mesh_file, [0, 0, 0], 1e-3, "FS network volume")

si = SnuddaInit(network_path=network_path, random_seed=123)

si.define_structure(struct_name="StriatalVolume", struct_mesh=mesh_file, d_min=15e-6, mesh_bin_width=50e-6)

# Should be 1050 neurons, temp reducing it to 50 neurons for runtime of simulation while developing
si.add_neurons(name="FS", num_neurons=1050, volume_id="StriatalVolume",
               neuron_dir=os.path.join("$DATA", "neurons", "striatum", "fs"))

cluster_FS_synapses = True
FS_gGABA = [1.1e-9, 1.5e-9]
FS_gGapJunction = [0.5e-9, 0.1e-9]

si.add_neuron_target(neuron_name="FS",
                     target_name="FS",
                     connection_type="GABA",
                     dist_pruning=None,
                     f1=0.15, soft_max=5, mu2=2, a3=1,
                     conductance=FS_gGABA,
                     cluster_synapses=cluster_FS_synapses,
                     mod_file="tmGabaA",
                     channel_param_dictionary={"tau1": (1.33e-3, 1e3),
                                               "tau2": (5.7e-3, 1e3)})
    
si.add_neuron_target(neuron_name="FS",
                     target_name="FS",
                     connection_type="GapJunction",
                     dist_pruning=None,
                     f1=0.7, soft_max=8, mu2=2, a3=1.0,
                     conductance=FS_gGapJunction,
                     cluster_synapses=False,
                     channel_param_dictionary=None)

si.write_json()


