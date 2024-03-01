#!/bin/bash

export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=default

if [ -z "$1" ]
then ipcluster start --profile=$IPYTHON_PROFILE --ip=127.0.0.1 &> ipcluster-log.txt &
else ipcluster start -n $1 --profile=$IPYTHON_PROFILE --ip=127.0.0.1 &> ipcluster-log.txt &
fi

echo "Sleeping 20 seconds to wait for workers to start"
sleep 20

echo "To stop ipcluster use:   ipcluster stop"
