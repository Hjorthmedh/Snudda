#!/bin/bash

export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

ipcluster start -n 6 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

# Optimize S1

#python3 OptimiseSynapsesFull.py --st glut DATA/YvonneJohansson2019/S1_Analysis_191001.h5 

# FS
python3 OptimiseSynapsesFull.py --st glut DATA/YvonneJohansson2019/M1LH_Analysis_191001.h5 --id 2,6,63,64,114,116

# SPN
#python3 OptimiseSynapsesFull.py --st glut DATA/YvonneJohansson2019/M1LH_Analysis_191001.h5 --id 3,5,58,15,16,62


#python3 OptimiseSynapsesFull.py --st glut DATA/YvonneJohansson2019/M1RH_Analysis_190925.h5

#python3 OptimiseSynapsesFull.py --st glut DATA/YvonneJohansson2019/TH_Analysis_191001.h5


