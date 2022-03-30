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

