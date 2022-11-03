#!/bin/bash

L=/cfs/klemming/projects/snic/snic2021-5-492/$USER

SNUDDA_DIR=$HOME/Snudda/snudda
JOBDIR=$L/PD_networks

SIMNAMEPART=50k_v1
SIMSIZE=50000
SIMNAME=$L/pd0_$SIMNAMEPART

#BasalGangliaData directory must exist

#BasalGangliaData/Parkinson/PD0
export SNUDDA_PD_DATA="$HOME/BasalGangliaData/Parkinson/20220225"
echo "Setting SNUDDA_DATA to $SNUDDA_DATA"

mkdir -p $JOBDIR

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

    echo ">>> Init: "`date`
    snudda init ${SIMNAME} --size ${SIMSIZE} --overwrite  --connectionFile $SNUDDA_PD_DATA/connectivity/network-config.json --randomseed 1234 --snudda_data "$SNUDDA_PD_DATA/PD0"

    if [ $? != 0 ]; then
	echo "Something went wrong during init, aborting!"
	ipcluster stop	
	exit -1
    fi

    echo "SLURM_NODELIST = $SLURM_NODELIST"
    let NWORKERS="$SLURM_NTASKS - 1"

    echo ">>> NWORKERS " $NWORKERS
    echo ">>> Starting ipcluster `date`"
    
    #.. Start the ipcluster
    ipcluster start -n ${NWORKERS} \
	      --ip='*' \
	      --HeartMonitor.max_heartmonitor_misses=1000 \
	      --HubFactory.registration_timeout=600 \
	      --HeartMonitor.period=10000 & 

    
    #.. Sleep to allow engines to start
    echo ">>> Wait 120s to allow engines to start"
    sleep 120 #60

    echo ">>> Place: "`date`
    snudda place ${SIMNAME} --verbose --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during placement, aborting!"
	# ipcluster stop	
	exit -1
    fi

    echo ">>> Detect: "`date`
    snudda detect ${SIMNAME} --hvsize 50 --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during detection, aborting!"
	# ipcluster stop	
	exit -1
    fi
    
    echo ">>> Prune: "`date`
    snudda prune ${SIMNAME} --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during pruning, aborting!"
	# ipcluster stop	
	exit -1
    fi

    SIMNAME=$L/pd1_$SIMNAMEPART

    echo ">>> Init: "`date`
    snudda init ${SIMNAME} --size ${SIMSIZE} --overwrite  --connectionFile $SNUDDA_PD_DATA/connectivity/network-config-PD-synapse-recovery.json --randomseed 1234 --snudda_data "$SNUDDA_PD_DATA/PD2"

    if [ $? != 0 ]; then
	echo "Something went wrong during init, aborting!"
	# ipcluster stop	
	exit -1
    fi

    echo ">>> Place: "`date`
    snudda place ${SIMNAME} --verbose --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during placement, aborting!"
	# ipcluster stop	
	exit -1
    fi

    echo ">>> Detect: "`date`
    snudda detect ${SIMNAME} --hvsize 50 --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during detection, aborting!"
	# ipcluster stop	
	exit -1
    fi
    
    echo ">>> Prune: "`date`
    snudda prune ${SIMNAME} --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during pruning, aborting!"
	# ipcluster stop	
	exit -1
    fi
   
    SIMNAME=$L/pd2_$SIMNAMEPART

    echo ">>> Init: "`date`
    snudda init ${SIMNAME} --size ${SIMSIZE} --overwrite  --connectionFile $SNUDDA_PD_DATA/connectivity/network-config-PD-synapse-recovery.json --randomseed 1234 --snudda_data "$SNUDDA_PD_DATA/PD2"

    if [ $? != 0 ]; then
	echo "Something went wrong during init, aborting!"
	# ipcluster stop	
	exit -1
    fi

    echo ">>> Place: "`date`
    snudda place ${SIMNAME} --verbose --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during placement, aborting!"
	# ipcluster stop	
	exit -1
    fi

    echo ">>> Detect: "`date`
    snudda detect ${SIMNAME} --hvsize 50 --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during detection, aborting!"
	# ipcluster stop	
	exit -1
    fi
    
    echo ">>> Prune: "`date`
    snudda prune ${SIMNAME} --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during pruning, aborting!"
	# ipcluster stop	
	exit -1
    fi

    SIMNAME=$L/pd3_$SIMNAMEPART

    echo ">>> Init: "`date`
    snudda init ${SIMNAME} --size ${SIMSIZE} --overwrite  --connectionFile $SNUDDA_PD_DATA/connectivity/network-config-PD-synapse-recovery.json --randomseed 1234 --snudda_data "$SNUDDA_PD_DATA/PD2"

    if [ $? != 0 ]; then
	echo "Something went wrong during init, aborting!"
	# ipcluster stop	
	exit -1
    fi

    echo ">>> Place: "`date`
    snudda place ${SIMNAME} --verbose --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during placement, aborting!"
	# ipcluster stop	
	exit -1
    fi

    echo ">>> Detect: "`date`
    snudda detect ${SIMNAME} --hvsize 50 --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during detection, aborting!"
	# ipcluster stop	
	exit -1
    fi
    
    echo ">>> Prune: "`date`
    snudda prune ${SIMNAME} --parallel

    if [ $? != 0 ]; then
	echo "Something went wrong during pruning, aborting!"
	# ipcluster stop	
	exit -1
    fi

fi
