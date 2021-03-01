# Snudda
Create realistic networks of neurons, synapses placed using touch detection between axons and dendrites

# Contact details
Johannes Hjorth, Royal Institute of Technology (KTH)
Human Brain Project 2019
hjorth@kth.se

# Funding
Horizon 2020 Framework Programme (785907, HBP SGA2); Horizon 2020 Framework Programme (945539, HBP SGA3); Vetenskapsrådet (VR-M-2017-02806); Swedish e-science Research Center (SeRC); KTH Digital Futures.

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

