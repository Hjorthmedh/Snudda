export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=default


export numNeurons=1000
export simNamePart=1k_test_v0
export PD_DATA_DIR="../../../../BasalGangliaData/Parkinson/20220225"

ipcluster start -n 3 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

######################


simName=networks/pd0_$simNamePart

# Note that we need to start ipcluster after the SNUDDA_DATA environment
# variable is set, otherwise the workers will not know of it.

snudda init $simName --size $numNeurons --overwrite --connectionFile $PD_DATA_DIR/connectivity/network-config.json --seed 1234 --snudda_data "$PD_DATA_DIR/PD0"
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

python3 ../../../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5

####################

simName=networks/pd1_$simNamePart

snudda init $simName --size $numNeurons --overwrite --connectionFile $PD_DATA_DIR/connectivity/network-config.json --seed 1234 --snudda_data "$PD_DATA_DIR/PD1"
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

python3 ../../../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5

####################

simName=networks/pd2_$simNamePart

snudda init $simName --size $numNeurons --overwrite --connectionFile $PD_DATA_DIR/connectivity/network-config.json --seed 1234 --snudda_data "$PD_DATA_DIR/PD2"
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

python3 ../../../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5

#####################

simName=networks/pd3_$simNamePart

snudda init $simName --size $numNeurons --overwrite --connectionFile $PD_DATA_DIR/connectivity/network-config.json --seed 1234 --snudda_data "$PD_DATA_DIR/PD3"
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

python3 ../../../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5


ipcluster stop

