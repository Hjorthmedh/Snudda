export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

ipcluster start -n 4 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

# Before you start, make sure that num_dSPN=num_iSPN=4750, num_FS=130,
# num_LTS=num_ChIN = 0 in init_custom.py


simName=networks/pd0_20201218

snudda init $simName --size 10000 --overwrite
python3 init_custom.py $simName --cellspec data/parkinson-2020-12-17/pd0/
snudda place $simName 
snudda detect $simName --volumeID Striatum
snudda prune $simName

simName=networks/pd1_20201218

snudda init $simName --size 10000 --overwrite
python3 init_custom.py $simName --cellspec data/parkinson-2020-12-17/pd1/
snudda place $simName 
snudda detect $simName --volumeID Striatum
snudda prune $simName

simName=networks/pd2_20201218

snudda init $simName --size 10000 --overwrite
python3 init_custom.py $simName --cellspec data/parkinson-2020-12-17/pd2/
snudda place $simName 
snudda detect $simName --volumeID Striatum
snudda prune $simName


simName=networks/pd3_20201218

snudda init $simName --size 10000 --overwrite
python3 init_custom.py $simName --cellspec data/parkinson-2020-12-17/pd3/
snudda place $simName 
snudda detect $simName --volumeID Striatum
snudda prune $simName


ipcluster stop


