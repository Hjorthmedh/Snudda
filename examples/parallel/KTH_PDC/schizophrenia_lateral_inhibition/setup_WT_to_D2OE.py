import os
from snudda.place.create_cube_mesh import create_cube_mesh
from snudda import Snudda

ipython_profile = os.getenv("IPYTHON_PROFILE", "default")
ipython_dir = os.getenv("IPYTHONDIR")

parallel_flag = ipython_dir is not None

# See HaloperidolEffectSimulation.ipynb in Alex's szmod repository

network_path = os.getenv("NETWORK_PATH", default=os.path.join("networks", "WT_to_D2OE"))
snudda_data_wt = os.getenv("SNUDDA_DATA_WT", default="/home/hjorth/HBP/szmod/snudda/snudda_data/wt")
snudda_data_d2oe = os.getenv("SNUDDA_DATA_D2OE", default="/home/hjorth/HBP/szmod/snudda/snudda_data/d2oe")
connectivity_info = os.path.join("..", "snudda_data", "wt", "connectivity", "striatum", "striatum-connectivity.json")

original_input = os.path.join(network_path, "input-wt-spikes.hdf5")
new_input = os.path.join(network_path, "input-d2oe-spikes.hdf5")

wt_output = os.path.join(network_path, "simulations", "output-wt.hdf5")
d2oe_output = os.path.join(network_path, "simulations", "output-d2oe.hdf5")
halo_output = os.path.join(network_path, "simulations", "output-d2oe-haloperidol.hdf5")

print(f"Using {network_path = }")
print(f"Using {snudda_data_wt = }")
print(f"Using {snudda_data_d2oe = }")

if parallel_flag:
    print(f"Running in parallel")
else:
    print(f"Running in serial")


neuron_paths = [os.path.join(snudda_data_wt, "neurons", "striatum", "dspn"),
                os.path.join(snudda_data_wt, "neurons", "striatum", "ispn"),
                os.path.join(snudda_data_wt, "neurons", "striatum", "fs")]

snd = Snudda(network_path=network_path)
snd_init = snd.init_tiny(neuron_paths=neuron_paths,
                         neuron_names=["dSPN", "iSPN", "FS"], number_of_neurons=[50, 50, 2],
                         connection_config=connectivity_info,
                         density=80500, d_min=15e-6,
                         snudda_data=snudda_data_wt, random_seed=123)


network_config_file = os.path.join("config.SZ", "network.json")

mesh_file = os.path.join("mesh", "cube_mesh.obj")
create_cube_mesh(file_name=mesh_file,
                 centre_point=(0, 0, 0),
                 side_len=(10000.0/80500)**(1/3)*0.001)

snd = Snudda(network_path=network_path, parallel=parallel_flag, ipython_profile=ipython_profile)
snd.import_config(network_config_file=network_config_file, snudda_data=snudda_data, overwrite=True)
snd.create_network()

snd.setup_input(input_config="input.json")



