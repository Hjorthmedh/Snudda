export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=default


# OBS, currently init is commented out, so numneurons not used
export numNeurons=3000
export simNamePart=3k_swap
export PD_DATA_DIR="$HOME/HBP/BasalGangliaData/Parkinson/20220225"

ipcluster start -n 3 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

######################


simName=networks/pd0_$simNamePart


# Note that we need to start ipcluster after the SNUDDA_DATA environment
# variable is set, otherwise the workers will not know of it.

snudda init $simName --size $numNeurons --overwrite --seed 1234 --snudda_data "$PD_DATA_DIR/PD0"
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

# python3 ../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5



####################

simName=networks/pd2_$simNamePart

snudda init $simName --size $numNeurons --overwrite --connectionFile $PD_DATA_DIR/connectivity/network-config-PD.json --seed 1234 --snudda_data "$PD_DATA_DIR/PD2"
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

# python3 ../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5

ipcluster stop

exit

# !!! We are skipping PD1 and PD3 here

####################

simName=networks/pd1_$simNamePart

snudda init $simName --size $numNeurons --overwrite --connectionFile $PD/connectivity/network-config-PD.json --seed 1234 "$PD_DATA_DIR/PD1"
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

# python3 ../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5



#####################

simName=networks/pd3_$simNamePart

snudda init $simName --size $numNeurons --overwrite --connectionFile $PD_DATA_DIR/connectivity/network-config-PD.json --seed 1234 "$PD_DATA_DIR/PD3"
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

# python3 ../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5


ipcluster stop

