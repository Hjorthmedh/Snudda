
# OBS, currently init is commented out, so numneurons not used
export numNeurons=1000
export simNamePart=1k_cluster_close
export PD_DATA_DIR="$HOME/HBP/BasalGangliaData/Parkinson/20220225"

######################


simName=networks/pd0_$simNamePart


# Note that we need to start ipcluster after the SNUDDA_DATA environment
# variable is set, otherwise the workers will not know of it.

snudda init $simName --size $numNeurons --overwrite --seed 1234 --snudda_data "$PD_DATA_DIR/PD0"
snudda place $simName
snudda detect $simName
snudda prune $simName

# python3 ../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5



####################

simName=networks/pd2_$simNamePart

snudda init $simName --size $numNeurons --overwrite --connectionFile $PD_DATA_DIR/connectivity/network-config-PD.json --seed 1234 --snudda_data "$PD_DATA_DIR/PD2"
snudda place $simName 
snudda detect $simName
snudda prune $simName 

# python3 ../snudda/analyse/analyse_striatum.py $simName/network-synapses.hdf5

