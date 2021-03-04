# Notebook examples for Snudda

Here is a collection of Jupyter Notebooks, some of the workflows are split over multiple notebooks, and you have to run them in a specific order. Snudda uses HDF5 files and they can sometimes be locked by one Notebook, so make sure to shutdown the kernel in the old Notebook before proceeding to the next one.

## Network creation
* [simple_network_creation](simple_network_creation.ipynb) first example of network creation
* [analyse_network_connectivity](analyse_network_connectivity.ipynb) analyses network connectivity created by [simple_network_creation](simple_network_creation.ipynb).

* [simple_network_parallel](simple_network_parallel.ipynb) how to use ipyparallel when running Snudda

* [population_unit_network](population_unit_network.ipynb) how to define population units.
* [custom_slice_example](custom_slice_example.ipynb) shows how to create custom slice and define your own connectivity rules for neuron types.

## Input creation
* [input_generation_example_1](input_generation_example_1.ipynb) generate constant Poisson input (uses [simple_network_creation](simple_network_creation.ipynb))
* [input_generation_example_2_frequency_vectors](input_generation_example_2_frequency_vectors.ipynb) define Poisson input with multiple start/stop times (uses [simple_network_creation](simple_network_creation.ipynb)).
* [input_generation_example_3_correlation](input_generation_example_3_correlation.ipynb) finer control in input targeting (uses [population_unit_network](population_unit_network.ipynb))

* [input_tuning_example](input_tuning_example.ipynb) explore what input number and frequency are good neurons, e.g to avoid depolarisation block.


## Striatum example
* [striatum_example](striatum_example.ipynb) creates a small striatal network, increase number of neurons for the full version.
* [striatum_example_simulate](striatum_example_simulate.ipynb) sets up input and simulates the [striatum_example](striatum_example.ipynb) network.
* [striatum_example_plot](striatum_example_plot.ipynb) plots spike raster from [striatum_example_simulate](striatum_example_simulate.ipynb). Figure 4 in Methods paper has a larger version.

## Touch detection and pruning figures
* [touch_detection_hypervoxel_illustration](../illustrations/touch_detection_hypervoxel_illustration.ipynb) - Panel 2A,B
* [touch_detection_illustration](../illustrations/touch_detection_illustration.ipynb) - Panel 2C
* [pruning_illustration](../illustrations/pruning_illustration.ipynb) - Panel 2D-G
