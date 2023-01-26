export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=default

# If the BasalGangliaData directory exists, then use that for our data
    
ipcluster start --n=6 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

simName=networks/test-10k

#snudda init $simName --size 1760000
#snudda init $simName --size 100000
if [[ -d "../../BasalGangliaData/data" ]]; then
    snudda init $simName --size 10000 --overwrite --seed 1234 --snudda_data ../../BasalGangliaData/data
else
    snudda init $simName --size 10000 --overwrite --seed 1234
fi

snudda place $simName --parallel
snudda detect $simName --volumeID Striatum --parallel
snudda prune $simName --parallel

# Copy over template input, you might need to update the path here if not
# run from the examples directory
cp -a ../snudda/data/input_config/external-input-dSTR-scaled-v4.json $simName/input.json
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
mpiexec -x SNUDDA_DATA snudda simulate $simName --verbose

