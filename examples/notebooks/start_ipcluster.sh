export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=default

ipcluster start --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20


echo "To stop ipcluster use:   ipcluster stop"
