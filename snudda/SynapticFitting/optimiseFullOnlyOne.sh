#!/bin/bash


#python3 OptimiseSynapsesFull.py --st glut DATA/YvonneJohansson2019/M1RH_Analysis_190925.h5  --id 318

#exit

export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=Snudda_LOCAL

ipcluster start -n 6 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

python3 OptimiseSynapsesFull.py --st glut DATA/YvonneJohansson2019/M1RH_Analysis_190925.h5  --id 318

ipcluster stop
