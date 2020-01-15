def snudda_help_text():
  return """Snudda - Create a network of neurons, with connectivity constrained by morphology

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
"""
