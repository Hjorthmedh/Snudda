export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=default


# OBS, currently init is commented out, so numneurons not used
numNeurons=10000
simNamePart=10k_20211105-tune-FS-v1
#cellspecDir=data/parkinson-2020-12-17
#cellspecDir=data/parkinson-2021-01-07
#cellspecDir=../snudda/data/parkinson-2021-03-19
#cellspecDir=../../BasalGangliaData/Parkinson.OLD
export PD_DATA_DIR="../../../BasalGangliaData/Parkinson/20221213"

ipcluster start --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

######################


simName=networks/pd0_$simNamePart

export SNUDDA_DATA="$PD_DATA_DIR/PD0"
echo "Setting SNUDDA_DATA to $SNUDDA_DATA"

# Note that we need to start ipcluster after the SNUDDA_DATA environment
# variable is set, otherwise the workers will not know of it.

snudda init $simName --size $numNeurons --overwrite
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

python3 ../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5


####################

if false; then

  simName=networks/pd1_$simNamePart

  export SNUDDA_DATA="$PD_DATA_DIR/PD1"
  echo "Setting SNUDDA_DATA to $SNUDDA_DATA"

  snudda init $simName --size $numNeurons --overwrite --connectionFile $SNUDDA_DATA/connectivity/network-config-PD.json
  snudda place $simName --parallel
  snudda detect $simName --parallel
  snudda prune $simName --parallel

  python3 ../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5
fi

####################

simName=networks/pd2_$simNamePart

export SNUDDA_DATA="$PD_DATA_DIR/PD2"
echo "Setting SNUDDA_DATA to $SNUDDA_DATA"

snudda init $simName --size $numNeurons --overwrite --connectionFile $SNUDDA_DATA/connectivity/network-config-PD.json
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

python3 ../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5


#####################

if false; then

  simName=networks/pd3_$simNamePart

  export SNUDDA_DATA="$PD_DATA_DIR/PD3"
  echo "Setting SNUDDA_DATA to $SNUDDA_DATA"

  snudda init $simName --size $numNeurons --overwrite --connectionFile $SNUDDA_DATA/connectivity/network-config-PD.json
  snudda place $simName --parallel
  snudda detect $simName --parallel
  snudda prune $simName --parallel

  python3 ../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5

fi

