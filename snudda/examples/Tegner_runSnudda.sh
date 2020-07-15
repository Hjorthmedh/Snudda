#!/bin/bash -l

SNUDDA_DIR=/cfs/klemming/nobackup/${USER:0:1}/$USER/Snudda/snudda
$DATA=$SNUDDA_DIR/data


JOBDIR=$SNUDDA_DIR/../networks/TegnerRun.${SLURM_JOBID}
#JOBDIR=$HOME/networks/TegnerRun.${SLURM_JOBID}

#SIMSIZE=120000
SIMSIZE=10000
#SIMSIZE=200000
#SIMSIZE=1760000

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

    #.. Obtain infiniband ip - this is faster and also internal
    /sbin/ifconfig ib0 | head -n 2 | tail -n 1 | awk '{print $2}' > controller_ip.txt
    CONTROLLERIP=$(<controller_ip.txt)

    let NWORKERS="$SLURM_NTASKS - 1"

    echo ">>> NWORKERS " $NWORKERS
    echo ">>> Starting ipcluster"
    
    #.. Start the ipcluster
    ipcluster start -n ${NWORKERS} \
	      --ip=${CONTROLLERIP} \
	      --location=${CONTROLLERIP} \
	      --profile=${IPYTHON_PROFILE} \
	      --engines=MPIEngineSetLauncher --debug \
	      --HeartMonitor.max_heartmonitor_misses=1000 \
	      --HubFactory.registration_timeout=600 \
	      --HeartMonitor.period=10000 &

    #.. Sleep to allow engines to start
    echo ">>> Wait 120s to allow engines to start"
    sleep 120 #60

    #.. Run the self-installed version of python3
    #${ANACONDA_HOME}/bin/python3 badExample.py 
    # Change to: snudda init ${JOBDIR} --size ${SIMSIZE}
    # need to figure out how to get it ti find snudda so we can call it directly
    # instead of calling core.py
    echo ">>> Init: "`date`
    snudda init ${JOBDIR} --size ${SIMSIZE}

    if [ $? != 0 ]; then
	echo "Something went wrong during init, aborting!"
	ipcluster stop	
	exit -1
    fi

    echo ">>> Place: "`date`
    snudda place ${JOBDIR}

    if [ $? != 0 ]; then
	echo "Something went wrong during placement, aborting!"
	ipcluster stop	
	exit -1
    fi

    echo ">>> Detect: "`date`
    snudda detect ${JOBDIR} --hvsize 50 

    if [ $? != 0 ]; then
	echo "Something went wrong during detection, aborting!"
	ipcluster stop	
	exit -1
    fi

    echo ">>> Prune: "`date`
    snudda prune ${JOBDIR}

    if [ $? != 0 ]; then
	echo "Something went wrong during pruning, aborting!"
	ipcluster stop	
	exit -1
    fi

    # echo ">>> Input: "`date`
    # snudda input ${JOBDIR} --input ../data/config/input-tinytest-v6.json
    echo "Use Tegner_setupInput_SERIAL.job to generate the input spikes"

    #.. Shut down cluster
    ipcluster stop	

    date
    echo "JOB END "`date` start_time_network_connect.txt

fi


