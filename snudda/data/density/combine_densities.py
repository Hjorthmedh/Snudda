import numpy as np
from collections import OrderedDict
import json
from snudda.utils.numpy_encoder import NumpyEncoder

volume_name = "Striatum"

file_lookup = { "LTS" : "1001_STRd_Sst_dens_smooth.csv",
                "ChIN" : "252_STRd_Chat_dens_smooth.csv",
                "dSPN" : "352_STRd_Drd1_dens_smooth.csv",
                "iSPN" : "72109410_STRd_Adora2a_dens_smooth.csv",
                "FS" : "79556738_STRd_Pvalb_dens_smooth.csv"}

output_file = "dorsal_striatum_density.json"

density_data = OrderedDict()

for neuron_type, density_file in file_lookup.items():
    print(f"Loading {neuron_type} data from {density_file}")
    file_data = np.genfromtxt(density_file, delimiter=' ')
    
    coord = file_data[:,:3]
    density = file_data[:,3]

    density_data[neuron_type] = OrderedDict()
    density_data[neuron_type]["Coordinates"] = coord
    density_data[neuron_type]["Density"] = density

volume_density_data = OrderedDict()
volume_density_data[volume_name] = density_data
    
with open(output_file, "w") as f:
    json.dump(volume_density_data, f, indent=4, cls=NumpyEncoder)

print(f"Wrote {output_file}")
