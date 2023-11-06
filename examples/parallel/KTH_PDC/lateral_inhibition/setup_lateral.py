import os
import numpy as np
network_path = "networks/lateral_1"
# snudda_data = "$HOME/BasalGangliaData/data"
snudda_data = "../../../../../BasalGangliaData/data/"

duration=5

# import pdb
# pdb.set_trace()

# from snudda import SnuddaInit -- this works on local machine, but not on dardel
import snudda.init

n_DSPN = 30000
n_ISPN = 30000
n_FS = 0
n_LTS = 0
n_ChIN = 0

print("Starting SnuddaInit")
si = snudda.init.SnuddaInit(network_path=network_path, snudda_data=snudda_data, random_seed=12345, honor_stay_inside=False)
si.define_striatum(num_dSPN=n_DSPN, num_iSPN=n_ISPN, num_FS=n_FS, num_LTS=n_LTS, num_ChIN=n_ChIN,
                   volume_type="cube")

print("Adding population units")

# The centre of the cube is [0.00475, 0.004, 0.00775]. num_neurons is optional
si.add_population_unit_density(structure_name="Striatum", neuron_types=["dSPN", "iSPN"], 
                               unit_centre=np.array([0.00475, 0.004, 0.00775]) -np.array([0,50e-6,150e-6]),
                               probability_function="(d < 200e-6)*1", num_neurons=200)
si.add_population_unit_density(structure_name="Striatum", neuron_types=["dSPN", "iSPN"], 
                               unit_centre=np.array([0.00475, 0.004, 0.00775]) -np.array([150e-6,0,0]),
                               probability_function="(d < 200e-6) * 1", num_neurons=200)

print("Writing json")

si.write_json()
