def snudda_help_text():
  return """Snudda - Create a network of neurons, with connectivity constrained by morphology

  Usage:

  snudda init <network_path> --size XXX
  -- Creates an a json config file

  snudda create <network_path>
  -- Create a network defined by json config file
     (this performs place, detect and prune)

  snudda input <networkPath> [--input yourInputConfig]
  -- Setup the input, obs you need to manually pick a input config file

  snudda simulate <networkPath>
  -- Run the network simulation using NEURON

  Legacy commands:
     
  snudda place <networkPath>
  -- Cell placement within volumes specified

  snudda detect <networkPath> [--hvsize hyperVoxelSize]
  -- Touch detection of putative synapses

  snudda prune <networkPath> [--mergeonly]
  -- Prune the synapses


  snudda export <networkPath>
  -- Export to SONATA format (optional)

  snudda analyse <networkPath>

  snudda help me
  -- Show this help text
"""
