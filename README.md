# Snudda
Create realistic networks of neurons, synapses placed using touch detection between axons and dendrites

# Contact details
Johannes Hjorth, Royal Institute of Technology (KTH)
Human Brain Project 2019
hjorth@kth.se

# Funding
This open source software code was developed in part or in whole in
the Human Brain Project, funded from the European Union's Horizon
2020 Framework Programme for Research and Innovation under Specific
Grant Agreements No. 720270 and No. 785907 (Human Brain Project SGA1
and SGA2).

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

  snudda export <networkPath>
  -- Export to SONATA format (optional)

  snudda simulate <networkPath>
  -- Run the network simulation using neuron

  snudda analyse <networkPath>

  snudda help me
  -- Show this help text

Example of a small simulation of the striatal circuitry can be found in Jupyter notebook [StriatumScaffoldExample-tiny.ipynb](./snudda/examples/StriatumScaffoldExample-tiny.ipynb).
