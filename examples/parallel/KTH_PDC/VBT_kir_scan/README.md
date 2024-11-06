# Running detailed simulations of Striatum for VBT

This code creates a striatal network with 10000 neurons based on touch detection using Snudda. The network is simulated, first in the WT configuration, then with decreased KIR current in dSPN and iSPN neurons.

## Snudda installation instructions

[Installation instructions](https://github.com/Hjorthmedh/Snudda/wiki/1.-User-installation)

This code also requires [BasalGangliaData](https://github.com/hjorthmedh/BasalGangliaData)


## Quick start:

To build the network and run it (note this SLURM file only uses one node, so will take 11 hours per simulation on Dardel).

```
Dardel_run_VBT_kir_scan.job
```

and to run only the simulations with 10-20 nodes for shorter run time:

```
Dardel_run_VBT_kir_scan.job-sim-only
```

## Code structure

The program ```setup_kir_scan.py``` takes ```network_path``` as an argument, and optionally ```--python_profile default``` and creates the network that is defined in the ```config/network.json``` file. It is a striatal network with 10000 neurons containing dSPN, iSPN, FS, LTS and ChIN. The input is specified in ```input.json```, here a steady 2Hz input to each synapse on the neurons. 

The program ```simulate_kir_scan.py``` takes ```network_path```, ```--kir_factor``` (between 0.0 and 1.0, ie WT), ```--output``` specifying the simulation output file and ```--time``` the duration (in seconds) of the simulation.

Both these programs are called from ```Dardel_run_VBT_kir_scan.job```. In the default configuration the kir_factor used is 1.0, 0.8, 0.6, 0.4, 0.2 and 0.

The network structure can be inspected with ```snudda_load your_simulation_file.hdf5```, to get instructions use ```snudda_load --help```. This uses ```snudda/utils/load.py``` which you can also call directly in python.

The simulation data can be loaded using ```snudda/utils/load_network_simulation.py```, if you want to plot individual traces you can use ```snudda_plot_traces```, use ```--help``` as an option to get instructions.


## Load example

More information about how to extract the simulated data can be found [here](LoadDataExample.ipynb).
