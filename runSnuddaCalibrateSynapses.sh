# Calibrate the synapses in a network, edit snudda_init_custom.py

# Make sure the defineStriatum line looks like this:
#   cnc.defineStriatum(nMSD1=120,nMSD2=120,nFS=20,nLTS=0,nChIN=0,volumeType="slice")


export simName="networks/Planert2010-v1"

python3 snudda_calibrate_synapses.py setup Planert2010 $simName

unset IPYTHONDIR
unset IPYTHON_PROFILE

./snudda.py place $simName 

export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

ipcluster stop
ipcluster start -n 12 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&

sleep 10

./snudda.py detect $simName --volumeID Striatum
./snudda.py prune $simName

unset IPYTHONDIR
unset IPYTHON_PROFILE

ipcluster stop

echo "Type Ctrl+D after inspecting the cut"

python3 snudda_cut.py $simName/network-pruned-synapses.hdf5 "abs(z)<100e-6"
#python3 snudda_cut.py $simName/network-pruned-synapses.hdf5 "abs(z-0.00511)<100e-6"

mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run Planert2010 $simName/network-cut-slice.hdf5 --pre dSPN --post iSPN

mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run Planert2010 $simName/network-cut-slice.hdf5 --pre iSPN --post dSPN

mpiexec -n 12 -map-by socket:OVERSUBSCRIBE python3 snudda_calibrate_synapses.py run Planert2010 $simName/network-cut-slice.hdf5 --pre FSN --post ALL

python3 snudda_calibrate_synapses.py analyse Planert2010 $simName --pre dSPN --post iSPN
python3 snudda_calibrate_synapses.py analyse Planert2010 $simName --pre dSPN --post dSPN
python3 snudda_calibrate_synapses.py analyse Planert2010 $simName --pre iSPN --post iSPN
python3 snudda_calibrate_synapses.py analyse Planert2010 $simName --pre iSPN --post dSPN
python3 snudda_calibrate_synapses.py analyse Planert2010 $simName --pre FSN --post iSPN
python3 snudda_calibrate_synapses.py analyse Planert2010 $simName --pre FSN --post dSPN
