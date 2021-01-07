export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

numNeurons=3000
simNamePart=20210107

ipcluster start -n 4 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

# Before you start, make sure that num_dSPN=num_iSPN=4750, num_FS=130,
# num_LTS=num_ChIN = 0 in init_custom.py


######################

simName=networks/pd0_$simNamePart

snudda init $simName --size $numNeurons --overwrite
python3 init_custom.py $simName --cellspec data/parkinson-2020-12-17/pd0/
snudda place $simName 
snudda detect $simName --volumeID Striatum
snudda prune $simName

simName=networks/pd1_$simNamePart

snudda init $simName --size $numNeurons --overwrite
python3 init_custom.py $simName --cellspec data/parkinson-2020-12-17/pd1/
snudda place $simName 
snudda detect $simName --volumeID Striatum
snudda prune $simName

simName=networks/pd2_$simNamePart

snudda init $simName --size $numNeurons --overwrite
python3 init_custom.py $simName --cellspec data/parkinson-2020-12-17/pd2/
snudda place $simName 
snudda detect $simName --volumeID Striatum
snudda prune $simName


simName=networks/pd3_$simNamePart

snudda init $simName --size $numNeurons --overwrite
python3 init_custom.py $simName --cellspec data/parkinson-2020-12-17/pd3/
snudda place $simName 
snudda detect $simName --volumeID Striatum
snudda prune $simName


ipcluster stop


