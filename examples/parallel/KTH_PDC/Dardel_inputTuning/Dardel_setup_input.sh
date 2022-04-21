#!/bin/bash

JOBDIR=networks/$1_1k_spn04

INPUTDIR=input

if [[ -d "bgd01/Parkinson/20211105-PD0/$1" ]]; then
    export SNUDDA_DATA="bgd01/Parkinson/20211105-PD0/$1"
    echo "Setting SNUDDA_DATA to $SNUDDA_DATA"
else
    echo "SNUDDA_DATA environment variable not changed (may be empty): $SNUDDA_DATA"
fi

########## Needed for NEURON

L=/cfs/klemming/scratch/${USER:0:1}/$USER/Projects/SnuddaProj09/snudda_env
LN=$L/neuron

export PATH=$L/bin:$LN/bin:$PATH

export LD_LIBRARY_PATH=$MPICH_DIR/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$L/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=$LN/lib/python:$PYTHONPATH

########### Needed for NEURON end.
echo "PATH: "$PATH
echo "PYTHONPATH: "$PYTHONPATH
echo "LD_LIBRARY_PATH: "$LD_LIBRARY_PATH
echo "SLURM_NODELIST = $SLURM_NODELIST"

let NWORKERS="$SLURM_NTASKS - 1"

echo ">>> Running sundda input"
# snudda input ${JOBDIR}
#snudda input $JOBDIR --input input_config/striatum-test-input.json 
snudda input $JOBDIR --input $INPUTDIR/external-input-dSTR-scaled-v4-bobek11.json --randomseed 2651382891 --useMeta 0
#snudda input $JOBDIR --input $INPUTDIR/noimp.json --randomseed 2651382891

echo "JOB END "`date`