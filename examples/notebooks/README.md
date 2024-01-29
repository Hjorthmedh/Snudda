# Notebook examples for Snudda

Here is a collection of Jupyter Notebooks, some of the workflows are split over multiple notebooks, and you have to run them in a specific order. Snudda uses HDF5 files and they can sometimes be locked by one Notebook, so make sure to shutdown the kernel in the old Notebook before proceeding to the next one.

## Network creation
* [simple_network_creation](simple_network_creation.ipynb) first example of network creation
* [analyse_network_connectivity](analyse_network_connectivity.ipynb) analyses network connectivity created by [simple_network_creation](simple_network_creation.ipynb).
* [simple_network_parallel](simple_network_parallel.ipynb) how to use ipyparallel when running Snudda
* [custom_slice_example](custom_slice_example.ipynb) shows how to create custom slice and define your own connectivity rules for neuron types.
* [population_unit_network](population_unit_network.ipynb) how to define population units.
* [example_of_density_function](example_of_density_function.ipynb) how to specify density variations using a function of (x,y,z) in a volume.
* [example_of_neuron_rotations](example_of_neuron_rotations.ipynb) shows how to rotate neurons based on position.
* [bend_morphologies](bend_morphologies.ipynb) shows how to make the neurons bend the axons and dendrites at the edge of the mesh, to keep them constrained to the volume.
* [connect_structures_example](connect_structures_example.ipynb) shows how to create neuron projections between volumes when no-axon data is available ([parallel version](connect_structures_example_parallel.ipynb)). There is also an [alternative version](connect_structures_example_projection_detection.ipynb) that places axon voxels randomly within the projection zone and then applies touch detection.
* [virtual_neurons](VirtualNeurons/VirtualNeurons.ipynb) shows how to only simulate the core of a volume of neurons, by turning the outer neurons to virtual neurons that spike at predetermined times. This is useful to avoid edge effects in your simulations.

## Input creation
* [input_generation_example_1](input_generation_example_1.ipynb) generate constant Poisson input (uses [simple_network_creation](simple_network_creation.ipynb))
* [input_generation_example_2_frequency_vectors](input_generation_example_2_frequency_vectors.ipynb) define Poisson input with multiple start/stop times (uses [simple_network_creation](simple_network_creation.ipynb)).
* [input_generation_example_3_correlation](input_generation_example_3_correlation.ipynb) finer control in input targeting (uses [population_unit_network](population_unit_network.ipynb))
* [input_tuning_example](input_tuning_example.ipynb) explore what input number and frequency are good neurons, e.g to avoid depolarisation block.
* [Advanced input tuning](../parallel/KTH_PDC/input_tuning/) using the Dardel super computer. Here both background and signal level are optimised, and the result is written to meta.json in BasalGangliaData.
* [input_frequency](InputFrequency/InputFrequency.ipynb) time varying input frequency example
* [basal_ganglia_data_input](Input-generation-using-BasalGangliaData.ipynb) builds upon the input tuned in BasalGangliaData.

## Striatum example
* [striatum_example](striatum_example.ipynb) creates a small striatal network, increase number of neurons for the full version.
* [striatum_example_simulate](striatum_example_simulate.ipynb) sets up input and simulates the [striatum_example](striatum_example.ipynb) network.

## Complete Snudda example for EBRAINS
* [Create custom network](EBRAINS/CreateMicrocircuit.ipynb) with more control, and [simulate network](EBRAINS/SimulateMicrocircuitry.ipynb)


## Visualisation
* [blender_example](blender_example.ipynb) - Use Blender to visualise a network, or part of a network
* [Hjorth, Hellgren Kotaleski, Kozlov](../Neuroinformatics2021/), Neuroinformatics 2021 figures


## Additional notebooks
* [Paired current injection](validation/synapses/network_pair_pulse_simulation.ipynb) in a simulated network
* [Neuromodulation examples](https://github.com/jofrony/Neuromodulation-software/tree/main/examples) in a network