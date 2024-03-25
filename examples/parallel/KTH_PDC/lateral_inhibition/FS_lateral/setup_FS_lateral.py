import os
import sys
import numpy as np

if len(sys.argv) > 1:
    network_path = sys.argv[1]
else:
    sys.exit("No network path specified!")
    network_path="../networks/FS_lateral"

# network_path = "networks/lateral_1"
# snudda_data = "$HOME/BasalGangliaData/data"
snudda_data = "../../../../../../BasalGangliaData/data/"

print(f"Network_path = {network_path}, snudda data = {snudda_data}")


duration=18

import snudda.init

n_FS = 1000   # 1.3% of total population, obs need to set density correctly below!

print("Starting SnuddaInit")
si = snudda.init.SnuddaInit(network_path=network_path, snudda_data=snudda_data, random_seed=12345, honor_stay_inside=False)
si.define_striatum(num_dSPN=0, num_iSPN=0, num_FS=n_FS, num_LTS=0, num_ChIN=0,
                   volume_type="cube", neuron_density=85000*0.013)

print("Adding population units")

# The centre of the cube is [0.00475, 0.004, 0.00775]. num_neurons is optional
si.add_population_unit_density(structure_name="Striatum", neuron_types=["FS"], 
                               unit_centre=np.array([0.00475, 0.004, 0.00775]) -np.array([0, 0, 0e-6]),
                               probability_function="(d < 300e-6) * 1", num_neurons=50)
si.add_population_unit_density(structure_name="Striatum", neuron_types=["FS"], 
                               unit_centre=np.array([0.00475, 0.004, 0.00775]) -np.array([0, 0, -0e-6]),
                               probability_function="(d > 300e-6) * (d < 550e-6) * 1 ", num_neurons=500)

print("Writing json")

si.write_json()
