export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

ipcluster start -n 4 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
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
#simName=networks/Article-jitter-v6
#simName=networks/Article-var-v2-100000
simName=networks/tinySim10
#simName=networks/Article-nojitter

#snudda init $simName --size 1760000
#snudda init $simName --size 100000
snudda init $simName --size 1000 --overwrite

snudda place $simName 
#snudda detect $simName --hvsize 50 --volumeID Striatum
snudda detect $simName --volumeID Striatum
snudda prune $simName

# Copy over template input
cp -a data/config/input-tinytest-v8.json $simName/input.json
echo "Make sure the input config file was found, otherwise provide your own"

# TODO, maybe use to get snudda base install dir:
# import inspect
# import snudda
# inspect.getfile(snudda) <--- can use this for path


# Uncomment this to generate input
snudda input $simName --input $simName/input.json

ipcluster stop

# Uncomment this to run simulation
# Remember you need to run "nrnivmodl data/cellspecs/mechanisms"
# first to create the mechanisms
mpiexec -n 6 snudda simulate $simName

