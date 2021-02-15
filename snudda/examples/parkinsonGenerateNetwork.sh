export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

# OBS, currently init is commented out, so numneurons not used
numNeurons=100000
simNamePart=20210107
#cellspecDir=data/parkinson-2020-12-17
cellspecDir=data/parkinson-2021-01-07

ipcluster start -n 4 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

# Before you start, make sure that num_dSPN=num_iSPN=4750, num_FS=130,
# num_LTS=num_ChIN = 0 in init_custom.py


######################

simName=networks/pd0_100k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 init_custom.py $simName --cellspec $cellspecDir/pd0/
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

simName=networks/pd1_100k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 init_custom.py $simName --cellspec $cellspecDir/pd1/
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

simName=networks/pd2_100k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 init_custom.py $simName --cellspec $cellspecDir/pd2/
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel


simName=networks/pd3_100k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 init_custom.py $simName --cellspec $cellspecDir/pd3/
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel


ipcluster stop


