#!/bin/bash -l


# How to setup if new intallation:
#
# ipython profile create --parallel --profile=mpi
#
# Add the following lines to .ipython/profile_mpi/ipcluster_config.py
#
# c.IPClusterEngines.engine_launcher_class = 'MPIEngineSetLauncher'
# c.IPClusterStart.controller_launcher_class = 'MPIControllerLauncher'
#
# Update the following lines in .ipython/profile_mpi/ipcontroller_config.py
#
# c.DictDB.record_limit = 10240
# c.DictDB.size_limit = 10737418240
#
# If you run out of memory on the ipcontroller node, use this instead:
# c.HubFactory.db_class = 'SQLiteDB'
#

#.. The idea here is to only launch ipcluster and
#   subsequent python scripts from one mpi task.
#   We still have access to them all, so the launched
#   processes will be able to access the whole allocated
#   network.
#
#   By using infiniband connection, we can connect
#   compute nodes that to not have an external ip.
#   As a bonus, it is much faster than regular ethernet connections. 
#..

SNUDDA_DIR=/cfs/klemming/nobackup/${USER:0:1}/$USER/Snudda/snudda

JOBDIR=$SNUDDA_DIR/../networks/TegnerRun.${SLURM_JOBID}
#JOBDIR=$HOME/networks/TegnerRun.${SLURM_JOBID}

#SIMSIZE=120000
SIMSIZE=50000
#SIMSIZE=200000
#SIMSIZE=1760000

if [ $SLURM_PROCID -gt 0 ]; then
	mock_string="Not main process"

else

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

    echo ">>> Input: "`date`
    snudda input ${JOBDIR} --input ../data/config/input-tinytest.json

    #.. Shut down cluster
    ipcluster stop	

    date
    echo "JOB END "`date` start_time_network_connect.txt

fi


