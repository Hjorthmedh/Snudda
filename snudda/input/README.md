# Input JSON file format specification

The top level key of the input JSON file is the
```target_identifier```, which can be ```neuron_type``` (e.g. "dSPN",
"iSPN"), ```neuron_sub_type``` or ```neuron_id``` (e.g. "437", obs has
to be a string due to JSON format). For each neuron Snudda first
checks if the ```neuron_id``` exists in the config file, failing that
it checks ```neuron_name``` and finally ```neuron_type``` to see if
there are any inputs. Only the most specific block it finds in the
input JSON file will be used.

Note: Internally ```neuron_sub_type``` is called ```neuron_name```,
eg. "dSPN_0", "dSPN_1", and has its own subdirectory in the neuron
folder.

The hierarchy below specifies the input name for each of the inputs,
eg. ("cortical", "thalamic").

Example format:

```
{
  "target1_identifier": {
    "input1_name": {
       ...
       },
    "input2_name": {
       ...
       }
       
    },
  "target2_identifier": {
     ...
  }
}
```

There are several options that can be specified for a specific input.

Mandatory fields:

| ```generator``` | ```poisson```, ```lognormal```, ```csv``` |
| ```type``` | "AMPA_NMDA", "GABA", ... |
| ```num_inputs``` | integer > 0 |
| ```frequency``` | float >= 0 (Hz) |
| ```conductance``` | float >= 0 (S)|
| ```mod_file``` | string, eg. "tmGlut" |

Other fields fields:

| ```synapse_density``` | string, eg. "1.15*0.05/(1+exp(-(d-30e-6)/5e-6))" |
| ```correlation``` | float >= 0 |
| ```jitter``` | float >= 0 (s)|
| ```parameter_file``` | str, path to JSON file with parameters for MOD file |
| ```random_seed``` | int |
| ```spikes``` | list of floats |
| ```spike_file``` | CSV file that specify spikes
| ```num_soma_synapses``` | int, extra synapses on the soma |


If ```input_type``` is ```virtual_neuron``` then this specifies the spike times of the virtual neuron. 