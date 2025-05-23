#!/bin/bash



SNUDDA_DIR=$HOME/Snudda/snudda
JOBDIR=../networks/lateral_onepop

SIMSIZE=50000

# If the BasalGangliaData directory exists, then use that for our data
#/cfs/klemming/scratch/${USER:0:1}/$USER/BasalGangliaData/data
#BasalGangliaData/Parkinson/PD0
if [[ -d "$HOME/BasalGangliaData/data" ]]; then
    export SNUDDA_DATA="$HOME/BasalGangliaData/data"
    echo "Setting SNUDDA_DATA to $SNUDDA_DATA"
else
    echo "SNUDDA_DATA environment variable not changed (may be empty): $SNUDDA_DATA"
fi

mkdir -p $JOBDIR

echo "Dardel_runSnudda.sh should be started with srun -n 1, to only get one process"

echo "SLURM_PROCID = $SLURM_PROCID"

if [ "$SLURM_PROCID" -gt 0 ]; then
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

    echo ">>> Init: "`date`
    # snudda init ${JOBDIR} --size ${SIMSIZE} --overwrite --randomseed 1234
    python setup_network_one_population.py ${JOBDIR}
    
    if [ $? != 0 ]; then
	echo "Something went wrong during init, aborting!"
	ipcluster stop	
	exit -1
    fi

# WE NOW START IPCLUSTER USING ipcontroller.sh INSTEAD...
#
#    echo "SLURM_NODELIST = $SLURM_NODELIST"
#    let NWORKERS="$SLURM_NTASKS - 1"
#
#    echo ">>> NWORKERS " $NWORKERS
#    echo ">>> Starting ipcluster `date`"
#    
#    #.. Start the ipcluster
#    ipcluster start -n ${NWORKERS} \
#	      --ip='*' \
#	      --HeartMonitor.max_heartmonitor_misses=1000 \
#	      --HubFactory.registration_timeout=600 \
#	      --HeartMonitor.period=10000 & 
#
#    
#    #.. Sleep to allow engines to start
#    echo ">>> Wait 120s to allow engines to start"
#    sleep 120 #60

    echo ">>> Place: "`date`
    snudda place ${JOBDIR} --verbose

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
    
    echo ">>> Prune: "`date`
    snudda prune ${JOBDIR} --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during pruning, aborting!"
	# ipcluster stop	
	exit -1
    fi

    # Disable input generation at the moment

    # echo ">>> Ablate: "`date`
    # python ablate_network.py ${JOBDIR}
    
    echo ">>> Input: "`date`
    snudda input ${JOBDIR} --parallel --time 5 --input input.json
    # snudda input ${JOBDIR} --parallel --time 5 --input input.json --networkFile ${JOBDIR}/network-synapses-minimal.hdf5

    
    #.. Shut down cluster
    # ipcluster stop	
    #.. Shutdown ipcontroller
    echo "Shutting down ipcontroller"

    python ../../ipcontroller_shutdown.py


    date
    #echo "JOB END "`date` start_time_network_connect.txt

    echo "EXITING Dardel_runjob_lateral.sh"

fi
