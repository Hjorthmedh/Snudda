# Setting up Snudda on CSCS

Create a Snudda environment:

```
CSCS_setup.sh
```

Snudda uses NEURON to run simulations. To install NEURON run:

```
CSCS_setup_neuron.sh
```

## Test your Snudda installation

Generate a striatal network using Snudda:

```
sbatch CSCS_daint_runSnudda.job
```

Simulate the striatal network:

```
sbatch CSCS_daint_simulate.job
```

The network files are in ```/scratch/snx3000/$USER/networks/CSCS_Network```, and the simulation output in the ```simulations``` subdirectory under that folder.


