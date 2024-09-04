#!/bin/bash

# This is to prevent NEURON from trying to open display 
unset DISPLAY

export SNUDDA_DATA="$HOME/BasalGangliaData/data"

JOBDIR=../networks/dspn_DA_bath


echo "Running Dardel_setup_neuromodulation_network: $JOBDIR"

echo "SLURM_PROCID = $SLURM_PROCID"

if [ "$SLURM_PROCID" -gt 0 ]; then
    mock_string="Not main process"
else
    echo "Running"
    python Dardel_setup_neuromodulation_network.py
fi
