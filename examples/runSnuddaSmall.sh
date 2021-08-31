export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

# If the BasalGangliaData directory exists, then use that for our data
if [[ -d "../../BasalGangliaData/data" ]]; then
    export SNUDDA_DATA="../../BasalGangliaData/data"
    echo "Setting SNUDDA_DATA to $SNUDDA_DATA"
else
    echo "SNUDDA_DATA environment variable not changed (may be empty): $SNUDDA_DATA"
fi
    
ipcluster start --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

simName=networks/test-100

#snudda init $simName --size 1760000
#snudda init $simName --size 100000
snudda init $simName --size 100 --overwrite

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


# Uncomment this to generate input
snudda input $simName --input $simName/input.json --parallel

ipcluster stop

# Uncomment this to run simulation
# Remember you need to run "nrnivmodl data/cellspecs/mechanisms"
# first to create the mechanisms
mpiexec -n 6 -x SNUDDA_DATA snudda simulate $simName

