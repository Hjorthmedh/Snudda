#!/bin/bash


export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Default

ipcluster start --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

python3 optimise_synapses_full.py ../data/synapses/example_data/10_MSN12_GBZ_CC_H20.json --synapseParameters ../data/synapses/example_data/M1LH-contra_dSPN.json

ipcluster stop
