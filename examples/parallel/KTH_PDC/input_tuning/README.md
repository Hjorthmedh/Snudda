# Optimising synaptic input to neurons

A population of multi-compartmental neuron models usually have a
variation in excitability. This code helps you generate two types of
input background ( ```cortical_background``` and
```thalamic_background```) and signal input ( ```cortical``` and
```thalamic```).

## Prepare SNUDDA_DATA

You need to point ```SNUDDA_DATA``` to the ```data``` directory where
your neurons are located. The ```meta.json``` files belonging to all
neurons in that directory will be modified. In this example we use the
data from the ```BasalGangliaData``` repository, but you can also use
other locations.

This script reads information about the input parameters from
```input_config/input_info.json``` in the SNUDDA_DATA folder. If you
want to optimize for a neuron or an input_type not previously defined,
then you need to update this file. This is also where you can change
the input to be clustered etc.

## Run background simulations

The first step is to tune the background input. On Dardel we are able to run all combinations of parametersets and morphologies for a particular neuron type in one go.

```slurm Dardel_run_input_tuning_background_dspn.job```

This will call ```setup_input_tuning_background.py``` that generates
the network, and the input configuration file and input data. Here
each morphologykey-parameterkey combination is repeated multiple
times, to vary the number of input synapses. The simulation is then
repeated a number of times, to get different locations of the
synapses, as this also affects excitability.

This particular code optimises the input for the ```dspn``` neurons by
setting the ```SNUDDA_TUNE_NEURON``` to ```dspn```. The
```SEED_LIST``` indicates how many times the simulation is repeated,
and with what seed (this must match the ```SEED_LIST``` in the
analysis script later. It is a string that has a list of numbers, e.g. "[10,20,55]"

## Update background input in ```meta.json```

The notebook ```Analyse_input_tuning_background_dspn.ipynb``` contains
the code for doing the update of the ```meta.json``` files. Note that
you should only run the first part of the notebook that updates the
background input.

Next you need to copy over the now updated ```BasalGangliaData``` to Dardel.

## Run signal simulations

After that we want to run the cortical and thalamic signal simulations while still having the background input active.

To run the cortical and thalamic signal runs use:

```sbatch Dardel_run_input_tuning_cortical_signal_dspn.job ```

```sbatch Dardel_run_input_tuning_thalamic_signal_dspn.job ```

This will vary the number of synaptic inputs, to try and get a 10Hz
input frequency to result in a 10Hz output frequency. If needed you
can change the input frequency and the requested output frequency.

## Update ```meta.json```

For the final step, we want to write the cortical and thalamic
(signal) input to the ```meta.json``` files. Here we will set the
```input_frequency``` to 0Hz for the signal, so the user can modify it
with the standard ```input-config.json``` file when needed.

