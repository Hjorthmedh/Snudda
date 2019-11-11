#!/bin/bash

export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

ipcluster start -n 12 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

# Optimize S1

python3 OptimiseSynapses.py --st glut DATA/YvonneJohansson2019/S1_Analysis_191001.h5 

python3 OptimiseSynapses.py --st glut DATA/YvonneJohansson2019/M1LH_Analysis_191001.h5

python3 OptimiseSynapses.py --st glut DATA/YvonneJohansson2019/M1RH_Analysis_190925.h5

python3 OptimiseSynapses.py --st glut DATA/YvonneJohansson2019/TH_Analysis_191001.h5


