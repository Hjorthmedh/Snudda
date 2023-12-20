#!/bin/bash

export SNUDDA_TUNE_NEURON="dspn"

NETWORK_DIR=networks/input_tuning_$SNUDDA_TUNE_NEURON

module load snic-env
source $HOME/Snudda/snudda_env/bin/activate
export SNUDDA_DIR=/cfs/klemming/home/"${USER:0:1}"/$USER/Snudda
export SNUDDA_DATA="/cfs/klemming/home/${USER:0:1}/$USER/BasalGangliaData/data"

python setup_input_tuning.py
