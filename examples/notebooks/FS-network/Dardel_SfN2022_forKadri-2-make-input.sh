#!/bin/bash

SNUDDA_DIR=$HOME/Snudda/snudda
JOBDIR=$1
INPUT_NAME=$2
DURATION=$3

echo "Using JOBDIR = $JOBDIR"

echo "We do not call init here, the config file has already been prepared"
# SIMSIZE=10000

# If the BasalGangliaData directory exists, then use that for our data
#/cfs/klemming/scratch/${USER:0:1}/$USER/BasalGangliaData/data
#BasalGangliaData/Parkinson/PD0
if [[ -d "$HOME/BasalGangliaData/data" ]]; then
    export SNUDDA_DATA="$HOME/BasalGangliaData/data"
    echo "Setting SNUDDA_DATA to $SNUDDA_DATA"
else
    echo "SNUDDA_DATA environment variable not changed (may be empty): $SNUDDA_DATA"
fi

# Dir should already exist
# mkdir -p $JOBDIR

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


    #echo ">>> Input: "`date`
    # cp -a $SNUDDA_DIR/data/input_config/input-v10-scaled.json ${JOBDIR}/input.json
    cp -a SfN2022-forKadri-background.json ${JOBDIR}/input.json

    snudda input ${JOBDIR} --parallel --time $DURATION --input {$INPUT_NAME}.json --inputFile ${JOBDIR}/{$INPUT_NAME}.hdf5
    
    #.. Shut down cluster
    ipcluster stop	

fi
