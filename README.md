## Summary of Snudda
Snudda creates the connectivity for realistic networks of simulated neurons in silico in a bottom up fashion that can then be simulated using the NEURON software. Neurons are placed within 3D meshes representing the structures of interest, with neural densities as seen in experiments. Based on reconstructed morphologies and neuron placement we can infer locations of putative synapses based on proximity between axon and dendrites. Projections between different structures can be added either using axon reconstructions, or by defining a connectivity map between regions. Putative synapses are pruned to match experimental pair-wise data on connectivity. Networks can be simulated either on desktop machines, or on super computers.

## Contact details
Johannes Hjorth, Royal Institute of Technology (KTH)
Human Brain Project
hjorth@kth.se

## Funding
Simulations were also performed on resources provided by the National Academic Infrastructure for Supercomputing in Sweden (NAISS) at PDC KTH partially funded by the Swedish Research Council through grant agreement no. 2022-06725.

The study was supported by the Swedish Research Council (VR-M-2020-01652), Swedish e-Science Research Centre (SeRC), Science for Life Lab, EU/Horizon 2020 no. 945539 (HBP SGA3) and No. 101147319 (EBRAINS 2.0 Project), European Union's Research and Innovation Program Horizon Europe under grant agreement No 101137289 (the Virtual Brain Twin Project), and KTH Digital Futures.

Horizon 2020 Framework Programme (785907, HBP SGA2); Horizon 2020 Framework Programme (945539, HBP SGA3); Vetenskapsrådet (VR-M-2017-02806, VR-M-2020-01652); Swedish e-science Research Center (SeRC); KTH Digital Futures. The computations are enabled by resources provided by the Swedish National Infrastructure for Computing (SNIC) at PDC KTH partially funded by the Swedish Research Council through grant agreement no. 2018-05973. We acknowledge the use of Fenix Infrastructure resources, which are partially funded from the European Union's Horizon 2020 research and innovation programme through the ICEI project under the grant agreement No. 800858. Snudda is supported and featured on EBRAINS.

## Citation
Please cite the first paper for the general Snudda network creation and simulation methods, and the second paper for the Striatal microcircutiry model.

* *Predicting Synaptic Connectivity for Large-Scale Microcircuit Simulations Using Snudda.* J. J. Johannes Hjorth, Jeanette Hellgren Kotaleski, Alexander Kozlov. Neuroinform (2021). https://doi.org/10.1007/s12021-021-09531-w

* *The microcircuits of striatum in silico.* J. J. Johannes Hjorth, Alexander Kozlov, Ilaria Carannante, Johanna Frost Nylén, Robert Lindroos, Yvonne Johansson, Anna Tokarska, Matthijs C. Dorst, Shreyas M. Suryanarayana, Gilad Silberberg, Jeanette Hellgren Kotaleski, Sten Grillner. Proceedings of the National Academy of Sciences (2020). https://doi.org/10.1073/pnas.2000671117

## Installation

To install Snudda:

```
pip3 install snudda
```

For more information, see Github:

https://github.com/Hjorthmedh/Snudda/wiki/1.-User-installation

## Jupyter Notebook examples

There are a number of examples for how to create and run networks on github which illustrates the functionality of Snudda. Several of these are created as short notebooks to showcase a particular feature or function.

https://github.com/Hjorthmedh/Snudda/tree/master/examples/notebooks

## Command line example

Once installed Snudda can also be run from the command line, using the snudda command. Below is a small list of the relevant commands that can be used.

Creates an a json config file:
```
snudda init <networkPath> --size XXX
```

Cell placement within volumes specified:
```
snudda place <networkPath>
```

Touch detection of putative synapses:
```
snudda detect <networkPath> [--hvsize hyperVoxelSize]
```

Prune the synapses
```
snudda prune <networkPath> [--mergeonly]
```

Setup the input, obs you need to manually pick a input config file
```
snudda input <networkPath> [--input yourInputConfig]
```

Run the network simulation using neuron
```
snudda simulate <networkPath>
```

Plot figurs with some network analysis:
```
snudda analyse <networkPath>
```

Show this help text
```
snudda help me
```


## Additional information:

https://snudda.readthedocs.io/
https://snudda.readthedocs.io/en/dev
