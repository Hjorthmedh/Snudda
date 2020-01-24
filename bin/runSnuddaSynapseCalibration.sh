# Run in parallel
python3 snudda_init_custom.py

python3 snudda.py place networks/SynTest-v1/
python3 snudda.py detect networks/SynTest-v1/
python3 snudda.py prune networks/SynTest-v1/

mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run networks/SynTest-v1/ dSPN iSPN

python3 snudda_calibrate_synapses.py analyse networks/SynTest-v1/ dSPN iSPN 

python3 Network_plot_traces.py networks/SynTest-v1/synapse-calibration-volt-dSPN-iSPN.txt networks/SynTest-v1/network-pruned-synapses.hdf5 
