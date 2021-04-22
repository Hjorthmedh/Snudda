export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

# OBS, currently init is commented out, so numneurons not used
numNeurons=50000
simNamePart=20210420
#cellspecDir=data/parkinson-2020-12-17
#cellspecDir=data/parkinson-2021-01-07
#cellspecDir=../snudda/data/parkinson-2021-03-19
cellspecDir=../snudda/data/parkinson-2021-04-20

ipcluster start -n 4 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

# Before you start, make sure that num_dSPN=num_iSPN=4750, num_FS=130,
# num_LTS=num_ChIN = 0 in init_custom.py


######################

simName=networks/pd0_50k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 ../snudda/init/init_custom.py $simName --cellspec $cellspecDir/pd0/
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

simName=networks/pd1_50k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 ../snudda/init/init_custom.py $simName --cellspec $cellspecDir/pd1/
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel

simName=networks/pd2_50k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 ../snudda/init/init_custom.py $simName --cellspec $cellspecDir/pd2/
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel


simName=networks/pd3_50k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 ../snudda/init/init_custom.py $simName --cellspec $cellspecDir/pd3/
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel


ipcluster stop


