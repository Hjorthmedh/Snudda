export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=default

# OBS, currently init is commented out, so numneurons not used
numNeurons=50000
simNamePart=20211022_Henri_1
#cellspecDir=data/parkinson-2020-12-17
#cellspecDir=data/parkinson-2021-01-07
#cellspecDir=../snudda/data/parkinson-2021-03-19
cellspecDir=../../BasalGangliaData/Parkinson



ipcluster start -n 4 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

# Before you start, make sure that num_dSPN=num_iSPN=4750, num_FS=130,
# num_LTS=num_ChIN = 0 in init_custom.py


######################

simName=networks/pd0_50k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 ../snudda/init/init_custom.py $simName --cellspec $cellspecDir/PD0/neurons
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel --keepfiles --savePutative

simName=networks/pd1_50k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 ../snudda/init/init_custom.py $simName --cellspec $cellspecDir/PD1/neurons
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel --keepfiles --savePutative

simName=networks/pd2_50k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 ../snudda/init/init_custom.py $simName --cellspec $cellspecDir/PD2/neurons/
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel --keepfiles --savePutative


simName=networks/pd3_50k_$simNamePart

# snudda init $simName --size $numNeurons --overwrite
python3 ../snudda/init/init_custom.py $simName --cellspec $cellspecDir/PD3/neurons
snudda place $simName --parallel
snudda detect $simName --parallel
snudda prune $simName --parallel --keepfiles --savePutative


ipcluster stop


