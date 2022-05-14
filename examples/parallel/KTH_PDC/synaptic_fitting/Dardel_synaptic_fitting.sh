#!/bin/bash

SNUDDA_DIR=$HOME/Snudda/snudda
JOBDIR=networks/synaptic_fitting

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

    echo "Start time: " > start_time_synaptic_fitting.txt
    date >> start_time_synaptic_fitting.txt

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

    python3 ../../../../snudda/synaptic_fitting/optimise_synapses_full.py ../../../../snudda/data/synapses/example_data/10_MSN12_GBZ_CC_H20.json --synapseParameters ../../../../snudda/data/synapses/example_data/M1LH-contra_dSPN.json --compile

    python3 ../../../../snudda/synaptic_fitting/optimise_synapses_full.py ../../../../snudda/data/synapses/example_data/10_MSN12_GBZ_CC_H20.json --synapseParameters ../../../../snudda/data/synapses/example_data/M1LH-contra_dSPN.json --plot


fi


    
