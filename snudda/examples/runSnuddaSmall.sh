export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

ipcluster start -n 12 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

# simName=networks/article100000-v1
#simName=networks/res5test
# simName=networks/article2019-v6
#simName=networks/article2019-v13
#simName=networks/SPN2SPN-distDepPruning-v5
#simName=networks/LTS-con-v7
#simName=networks/Robert-bug-v1
#simName=networks/NewConfig-v4
#simName=networks/FullStriatum-v2
#simName=networks/Article10000-v1
simName=networks/Article100000-v3

if [ -d "$simName" ]; then
  echo "Directory $simName already exists!!"
  exit 666    
fi

#snudda init $simName --size 1760000
snudda init $simName --size 100000

snudda place $simName 
#snudda detect $simName --hvsize 50 --volumeID Striatum
snudda detect $simName --volumeID Striatum
snudda prune $simName

# Copy over template input
cp -a config/input-tinytest-v5.json $simName/input.json

# Uncomment this to generate input
#snudda input $simName --input $simName/input.json

ipcluster stop

# Uncomment this to run simulation
#snudda simulate $simName

