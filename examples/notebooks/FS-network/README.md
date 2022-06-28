# Jupyter Notebooks

* [FS_network_1](FS_network_1.ipynb) creates and simulates a small network of gap junction coupled striatal fast spiking interneurons, driven by synaptic input.
* [FS_network_2](FS_network_2.ipynb) - ORPHANED NOTEBOOK
* [FS_network_3](FS_network_3.ipynb) creates a network of 100 striatal fast spiking interneurons, driven by current injection. First the network is run without gap junctions, then it is rerun with gap junctions, leading to increased synchronisation.
* [FS_network_4](FS_network_4.ipynb) creates a small network of striatal fast spiking interneurons, the network is driven by synaptic input. First the simulation is run without gap junctions, then the gap junctions are added.


# Simulations on Dardel Super Computer

# FS network driven by current injections

Setup and run the network:
```
Dardel_create_FS_network-cur-inj.job
Dardel_simulate_FS_network_2-cur-inj.job
```

To analyse the network run:
```
python3 ../../snudda/simulate/network_pair_pulse_simulation.py analyse Planert2010 FS_network_2-cur-inj --pre FS
```

# FS network simulation with synaptic input

Setup and run FS network connected by gap junctions

```
Dardel_create_FS_network.job
Dardel_simulate_FS_network.job
```

