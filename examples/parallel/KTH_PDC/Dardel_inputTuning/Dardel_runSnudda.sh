#!/bin/bash

SNUDDA_DIR=/cfs/klemming/scratch/${USER:0:1}/$USER/Projects/SnuddaProj09/Snudda/snudda

JOBDIR=networks/$1_1k_spn04

echo $JOBDIR

SIMSIZE=1000
#SIMSIZE=120000
#SIMSIZE=1760000

# If the BasalGangliaData directory exists, then use that for our data
#/cfs/klemming/scratch/${USER:0:1}/$USER/BasalGangliaData/data
#BasalGangliaData/Parkinson/PD0

#Parkinson/20211105/PD0     /neurons/striatum
#BasalGangliaData/data    /neurons/striatum
if [[ -d "bgd01/Parkinson/20211105-PD0/$1" ]]; then
    export SNUDDA_DATA="bgd01/Parkinson/20211105-PD0/$1"
    echo "Setting SNUDDA_DATA to $SNUDDA_DATA"
else
    echo "SNUDDA_DATA environment variable not changed (may be empty): $SNUDDA_DATA"
fi



mkdir -p $JOBDIR

if [ $SLURM_PROCID -gt 0 ]; then
	mock_string="Not main process"

else

    # For debug purposes:                                                         
    echo "PATH: "$PATH
    echo "IPYTHONDIR: "$IPYTHONDIR
    echo "PYTHONPATH: "$PYTHONPATH
    echo "LD_LIBRARY_PATH: "$LD_LIBRARY_PATH

    echo ">>>>>> Main process starting ipcluster"
    echo

    echo "Start time: " > start_time_network_connect.txt
    date >> start_time_network_connect.txt

    # Change to: snudda init ${JOBDIR} --size ${SIMSIZE}
    # need to figure out how to get it ti find snudda so we can call it directly
    # instead of calling core.py
    echo ">>> Init: "`date`
    #snudda init ${JOBDIR} --size ${SIMSIZE} --overwrite --randomseed  3385905604
    #python3 $SNUDDA_DIR/init/init_bosse.py ${JOBDIR} --cellspec Parkinson/PD0/
    #python3 $SNUDDA_DIR/init/init_custom.py ${JOBDIR} --cellspec /cfs/klemming/scratch/b/bobek/Projects/SnuddaProj03/Snudda/examples/parallel/BasalGangliaData/Parkinson/PD0/neurons
    python3 $SNUDDA_DIR/init/init_custom.py ${JOBDIR} --neurons $SNUDDA_DATA/neurons
    
    if [ $? != 0 ]; then
	echo "Something went wrong during init, aborting!"
	ipcluster stop	
	exit -1
    fi
    #################################################################################################################
    python removeConnections.py ${JOBDIR} #Removes all synapses between neurons, but does not touch cortical/thalamic input
    #################################################################################################################
    echo "SLURM_NODELIST = $SLURM_NODELIST"
    
    #.. Obtain infiniband ip - this is faster and also internal
    ifconfig > ifconfig-info.txt
    ifconfig ipogif0 | head -n 2 | tail -n 1 | awk '{print $2}' > controller_ip.txt
    CONTROLLERIP=$(<controller_ip.txt)

    let NWORKERS="$SLURM_NTASKS - 1"

    echo ">>> NWORKERS " $NWORKERS
    echo ">>> Starting ipcluster `date`"
    
    ipcluster start -n ${NWORKERS} \
	      --ip='*' \
	      --HeartMonitor.max_heartmonitor_misses=1000 \
	      --HubFactory.registration_timeout=600 \
	      --HeartMonitor.period=10000 & 

    
    #.. Sleep to allow engines to start
    echo ">>> Wait 120s to allow engines to start"
    sleep 60 #60, 120

    echo ">>> Place: "`date`
    snudda place ${JOBDIR} --verbose --parallel --randomseed 1132564546

    if [ $? != 0 ]; then
	echo "Something went wrong during placement, aborting!"
	# ipcluster stop	
	exit -1
    fi

    echo ">>> Detect: "`date`
    snudda detect ${JOBDIR} --hvsize 50 --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during detection, aborting!"
	# ipcluster stop	
	exit -1
    fi
echo "Sleeping 120 sec"
sleep 10 #60, 120
echo "Sleeping done"
    echo ">>> Prune: "`date`
    snudda prune ${JOBDIR} --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during pruning, aborting!"
	# ipcluster stop	
	exit -1
    fi

    # Disable input generation at the moment

    #echo ">>> Input: "`date`
    #cp -a $SNUDDA_DIR/data/input_config/input-v10-scaled.json ${JOBDIR}/input.json
    # cp -a $SNUDDA_DIR/data/input_config/external-input-dSTR-scaled.json ${JOBDIR}/input.json

    #snudda input ${JOBDIR} --parallel --time 5

    
    #.. Shut down cluster
    # ipcluster stop	

    date
    echo "JOB END "`date` start_time_network_connect.txt

fi