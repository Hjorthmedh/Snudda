#!/usr/bin/env python3
import json

import os
import uuid
import numpy as np

from snudda.place.create_cube_mesh import create_cube_mesh
from snudda import Snudda
from snudda.utils import SnuddaLoad


ipython_profile = os.getenv("IPYTHON_PROFILE", "default")
ipython_dir = os.getenv("IPYTHONDIR")

parallel_flag = ipython_dir is not None

snudda_data = "../../../../../../BasalGangliaData/data/"

network_path = os.path.join("networks", "small_network")

random_seed = 12345

n_dspn = 500  # 1000
n_ispn = 0  # 1000
n_fs = 20  # 40

n_total = n_dspn + n_ispn + n_fs
connectivity_info = os.path.join(snudda_data, "connectivity", "striatum", "striatum-connectivity.json")

neuron_paths = [os.path.join(snudda_data, "neurons", "striatum", "dspn"),
                os.path.join(snudda_data, "neurons", "striatum", "ispn"),
                os.path.join(snudda_data, "neurons", "striatum", "fs")]


# Important that the mesh has the correct size, because distance dependent connectivity is dependent on the density
mesh_file = os.path.join("mesh", "cube_mesh.obj")
create_cube_mesh(file_name=mesh_file,
                 centre_point=(0, 0, 0),
                 side_len=(n_total/80500.0)**(1/3)*0.001)  # 0.001 convert to meters




snd = Snudda(network_path=network_path, parallel=parallel_flag)
snd_init = snd.init_tiny(neuron_paths=neuron_paths,
                         neuron_names=["dSPN", "iSPN", "FS"],
                         number_of_neurons=[n_dspn, n_ispn, n_fs],
                         connection_config=connectivity_info,
                         density=80500,
                         d_min=15e-6,
                         snudda_data=snudda_data,
                         random_seed=random_seed)

snd_init.network_data["regions"]["Cube"]["neurons"]["dSPN"]["reaction_diffusion"] = os.path.join(snudda_data, "neurons/striatum/modulation","reaction_diffusion_D1.json")
snd_init.network_data["regions"]["Cube"]["neurons"]["dSPN"]["modulation"] = os.path.join(snudda_data, "neurons/striatum/modulation", "modulation_parameters.json")
snd_init.network_data["regions"]["Cube"]["neurons"]["dSPN"]["modulation_key"] = "abc"

if "iSPN" in snd_init.network_data["regions"]["Cube"]["neurons"]:
    snd_init.network_data["regions"]["Cube"]["neurons"]["iSPN"]["reaction_diffusion"] = os.path.join(snudda_data, "neurons/striatum/modulation", "reaction_diffusion_D2.json")
    snd_init.network_data["regions"]["Cube"]["neurons"]["iSPN"]["modulation"] =  os.path.join(snudda_data, "neurons/striatum/modulation", "modulation_parameters.json")
    snd_init.network_data["regions"]["Cube"]["neurons"]["iSPN"]["modulation_key"] = "abc"
else:
    print(f"Skipping iSPN")
    
snd_init.write_json()
    
snd.create_network()

snd.setup_input(input_config="input.json")

