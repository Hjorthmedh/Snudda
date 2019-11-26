mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run networks/SynTest-v5/network-cut-slice.hdf5 dSPN iSPN

mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run networks/SynTest-v5/network-cut-slice.hdf5 iSPN dSPN

mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run networks/SynTest-v5/network-cut-slice.hdf5 FSN ALL
