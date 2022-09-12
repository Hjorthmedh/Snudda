original_network = "/home/hjorth/HBP/DELME/networks/pd0_1k_swap"
updated_network = "/home/hjorth/HBP/DELME/networks/pd2_1k_swap"
output_network =  "/home/hjorth/HBP/DELME/networks/pd2_output_network"

import os

original_network_file = os.path.join(original_network, "network-synapses.hdf5")
updated_network_file = os.path.join(updated_network, "network-synapses.hdf5")
output_network_file =  os.path.join(output_network, "network-synapses.hdf5")

original_snudda_data_dir = "/home/hjorth/HBP/BasalGangliaData/Parkinson/20220225/PD0"
updated_snudda_data_dir = "/home/hjorth/HBP/BasalGangliaData/Parkinson/20220225/PD2"

# original_snudda_data_dir = "/home/hjorth/HBP/BasalGangliaData/Parkinson/20211105/PD0"
# updated_snudda_data_dir = "/home/hjorth/HBP/BasalGangliaData/Parkinson/20211105/PD2"


from snudda.plotting.plot_degeneration_and_growth import PlotDegenerationAndGrowth
pdg = PlotDegenerationAndGrowth(original_network_path=original_network, degenerated_network_path=output_network, neuron_id=500)
fig, ax = pdg.plot_synapses()


import pdb
pdb.set_trace()
