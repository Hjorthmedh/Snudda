# Calibrate the synapses in a network, edit snudda_init_custom.py

export simName="networks/SynTest-v28-LTS"

python3 snudda_init_custom.py $simName

./snudda.py place $simName 

export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

ipcluster start -n 12 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&

sleep 10

./snudda.py detect $simName --volumeID Striatum
./snudda.py prune $simName

export IPYTHONDIR=""
export IPYTHON_PROFILE=""

ipcluster stop

echo "Type Ctrl+D after inspecting the cut"

python3 snudda_cut.py $simName/network-pruned-synapses.hdf5 "abs(z)<100e-6"
#python3 snudda_cut.py $simName/network-pruned-synapses.hdf5 "abs(z-0.00511)<100e-6"

mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run $simName/network-cut-slice.hdf5 LTS dSPN

mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run $simName/network-cut-slice.hdf5 LTS iSPN

mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run $simName/network-cut-slice.hdf5 LTS ChIN

mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run $simName/network-cut-slice.hdf5 FS LTS

mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run $simName/network-cut-slice.hdf5 ChIN LTS
