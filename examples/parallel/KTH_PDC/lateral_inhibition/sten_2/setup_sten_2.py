import os
import sys
import numpy as np

if len(sys.argv) > 1:
    network_path = sys.argv[1]
else:
    sys.exit("No network path specified!")
    network_path="../networks/sten_2"

# network_path = "networks/lateral_1"
# snudda_data = "$HOME/BasalGangliaData/data"
snudda_data = "../../../../../../BasalGangliaData/data/"

print(f"Network_path = {network_path}, snudda data = {snudda_data}")


duration=18

import snudda.init

n_DSPN = 20000
n_ISPN = 20000
n_FS = 0
n_LTS = 0
n_ChIN = 0

print("Starting SnuddaInit")
si = snudda.init.SnuddaInit(network_path=network_path, snudda_data=snudda_data, random_seed=12345, honor_stay_inside=False)
si.define_striatum(num_dSPN=n_DSPN, num_iSPN=n_ISPN, num_FS=n_FS, num_LTS=n_LTS, num_ChIN=n_ChIN,
                   volume_type="cube")

print("Adding population units")

# si.add_population_unit_random(structure_name="Striatum", neuron_types=["dSPN", "iSPN"],
#                               fraction_of_neurons=0.5, unit_id=1)
# si.add_population_unit_random(structure_name="Striatum", neuron_types=["dSPN", "iSPN"],
#                               fraction_of_neurons=0.5, unit_id=2)

# The centre of the cube is [0.00475, 0.004, 0.00775]. num_neurons is optional
si.add_population_unit_density(structure_name="Striatum", neuron_types=["dSPN", "iSPN"], 
                               unit_centre=np.array([0.00475, 0.004, 0.00775]) -np.array([0, 0, 100e-6]),
                               probability_function="(d < 300e-6) * 1", num_neurons=4000)
si.add_population_unit_density(structure_name="Striatum", neuron_types=["dSPN", "iSPN"], 
                               unit_centre=np.array([0.00475, 0.004, 0.00775]) -np.array([0, 0, -100e-6]),
                               probability_function="(d < 300e-6) * 1", num_neurons=4000)

print("Writing json")

si.write_json()
