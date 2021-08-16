export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

ipcluster start -n 6 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

simName=networks/test-10k
# simName=networks/test-20k
# simName=networks/test-50k

#snudda init $simName --size 1760000
#snudda init $simName --size 100000
snudda init $simName --size 10000 --overwrite
# snudda init $simName --size 20000 --overwrite
#snudda init $simName --size 50000 --overwrite

snudda place $simName --parallel
snudda detect $simName --volumeID Striatum --parallel
snudda prune $simName --parallel

# Copy over template input, you might need to update the path here if not
# run from the examples directory
cp -a ../snudda/data/input_config/input-tinytest-v9-freq-vectors.json $simName/input.json
echo "Make sure the input config file was found, otherwise provide your own"

# TODO, maybe use to get snudda base install dir:
# import inspect
# import snudda
# inspect.getfile(snudda) <--- can use this for path


# !!! Uncomment the line below to generate input

# snudda input $simName --input $simName/input.json --parallel

ipcluster stop


# Remember you need to run "nrnivmodl data/cellspecs/mechanisms"
# first to create the mechanisms

# !!! Uncomment the line below to run simulation
# mpiexec -n 6 snudda simulate $simName

