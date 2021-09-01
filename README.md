# Summary of Snudda
Software to create realistic networks of simulated neurons in silico in a bottom up fashion. Neurons are placed within 3D meshes representing the structures of interest, with neural densities as seen in experiments. Based on reconstructed morphologies and neuron placement we can infer locations of putative synapses based on proximity between axon and dendrites. Projections between different structures can be added either using axon reconstructions, or by defining a connectivity map between regions. Putative synapses are pruned to match experimental pair-wise data on connectivity. Networks can be simulated either on desktop machines, or on super computers.

# Contact details
Johannes Hjorth, Royal Institute of Technology (KTH)
Human Brain Project
hjorth@kth.se

# Funding
Horizon 2020 Framework Programme (785907, HBP SGA2); Horizon 2020 Framework Programme (945539, HBP SGA3); Vetenskapsrådet (VR-M-2017-02806); Swedish e-science Research Center (SeRC); KTH Digital Futures.

# Citation
Please cite the first paper for the general Snudda network creation and simulation methods, and the second paper for the Striatal microcircutiry model.

* J. J. Johannes Hjorth, Jeanette Hellgren Kotaleski, Alexander Kozlov. Predicting Synaptic Connectivity for Large-Scale Microcircuit Simulations Using Snudda. Neuroinform (2021). https://doi.org/10.1007/s12021-021-09531-w

* J. J. Johannes Hjorth, Alexander Kozlov, Ilaria Carannante, Johanna Frost Nylén, Robert Lindroos, Yvonne Johansson, Anna Tokarska, Matthijs C. Dorst, Shreyas M. Suryanarayana, Gilad Silberberg, Jeanette Hellgren Kotaleski, Sten Grillner
The microcircuits of striatum in silico. Proceedings of the National Academy of Sciences (2020). https://doi.org/10.1073/pnas.2000671117

# Installation

https://github.com/Hjorthmedh/Snudda/wiki/1.-User-installation

# Jupyter Notebook examples

https://github.com/Hjorthmedh/Snudda/tree/master/examples/notebooks

# Usage

  Usage:

  snudda init <networkPath> --size XXX
  -- Creates an a json config file

  snudda place <networkPath>
  -- Cell placement within volumes specified

  snudda detect <networkPath> [--hvsize hyperVoxelSize]
  -- Touch detection of putative synapses

  snudda prune <networkPath> [--mergeonly]
  -- Prune the synapses

  snudda input <networkPath> [--input yourInputConfig]
  -- Setup the input, obs you need to manually pick a input config file

  snudda simulate <networkPath>
  -- Run the network simulation using neuron

  snudda analyse <networkPath>

  snudda help me
  -- Show this help text

