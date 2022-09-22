export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=default

echo "This code assumes you have BasalGangliaData installed"

ipcluster start -n 2 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20


for sim_size in 1000 5000 10000 50000; do

    sim_name=networks/arvind-$sim_size

    snudda init $sim_name --size $sim_size --overwrite --snudda_data ../../../../BasalGangliaData/data --seed 123456
    snudda place $sim_name --parallel
    snudda detect $sim_name --volumeID Striatum --parallel
    snudda prune $sim_name --parallel

done


ipcluster stop


