import os
import sys

if len(sys.argv) > 1:
    network_path = sys.argv[1]
else:
    sys.exit("No network path specified!")
    network_path="networks/lateral_1"
    
modified_network_file=os.path.join(network_path, "network-synapses-virtual.hdf5")

print(f"Network_path = {network_path}, modified file = {modified_network_file}")

from snudda.utils.ablate_network import SnuddaAblateNetwork

ab = SnuddaAblateNetwork(network_file=network_path)
pop_unit_0 = ab.snudda_load.get_population_unit_members(population_unit=0)  # Here, surrounding neurons
pop_unit_1 = ab.snudda_load.get_population_unit_members(population_unit=1)
pop_unit_2 = ab.snudda_load.get_population_unit_members(population_unit=2)

#ab.only_keep_neuron_id(neuron_id=set(pop_unit_1).union(set(pop_unit_2)))
ab.make_virtual(pop_unit_0)

ab.write_network(out_file_name=modified_network_file)



# Lets also do the virtual input while we are at it...
# since we have access to pop_unit_0 index here

from snudda.input.virtual_input import VirtualInput

duration = 5.5

spike_file = os.path.join(network_path, "virtual_input_spikes.txt")
mapping_file = os.path.join(network_path, "virtual_input_mapping.txt")

vi = VirtualInput(spike_file=spike_file, mapping_file=mapping_file)

from snudda.data.input_config.Kim2019.get_experimental_firing_freq import resample_spn_freq

spn_freq = resample_spn_freq(vidx.shape, rng=None)

for vidx, freq in zip(pop_unit_0, spn_freq):
    # vi.add_input(neuron_id=vidx, spike_times = vi.poisson_spikes(frequency=5, max_time=duration))
    vi.add_input(neuron_id=vidx, spike_times = vi.poisson_spikes(frequency=freq, max_time=duration))
    
vi.write_data()
