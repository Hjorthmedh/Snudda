

import os

network_path = os.getenv("NETWORK_PATH", default=os.path.join("networks", "SZ_scan"))
snudda_data = os.getenv("SNUDDA_DATA", default="/home/hjorth/HBP/szmod/snudda/snudda_data/wt")

print(f"Using {network_path = }")
print(f"Using {snudda_data = }")

network_config_file = os.path.join("config.SZ", "network.json")

from snudda.place.create_cube_mesh import create_cube_mesh

mesh_file = os.path.join("mesh", "cube_mesh.obj")
create_cube_mesh(file_name=mesh_file,
                 centre_point=(0, 0, 0),
                 side_len=(10000.0/80500)**(1/3)*0.001)

from snudda import Snudda
snd = Snudda(network_path=network_path)
snd.import_config(network_config_file=network_config_file, snudda_data=snudda_data, overwrite=True)
snd.create_network()



