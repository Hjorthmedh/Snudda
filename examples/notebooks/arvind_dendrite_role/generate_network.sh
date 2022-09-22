export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=default

# If the BasalGangliaData directory exists, then use that for our data
if [[ -d "../../BasalGangliaData/data" ]]; then
    export SNUDDA_DATA="../../BasalGangliaData/data"
    echo "Setting SNUDDA_DATA to $SNUDDA_DATA"
else
    echo "SNUDDA_DATA environment variable not changed (may be empty): $SNUDDA_DATA"
fi
    
ipcluster start -n 2 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20


for sim_size in 1000 5000 10000 50000; do

    sim_name=networks/arvind-$sim_size

    snudda init $sim_name --size $sim_size --overwrite
    snudda place $sim_name --parallel
    snudda detect $sim_name --volumeID Striatum --parallel
    snudda prune $sim_name --parallel

done


ipcluster stop


