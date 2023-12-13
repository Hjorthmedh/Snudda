#!/bin/bash



SNUDDA_DIR=$HOME/Snudda/snudda
JOBDIR=networks/lateral_1

SIMSIZE=50000

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

echo "Dardel_runSnudda.sh should be started with srun -n 1, to only get one process"

echo "SLURM_PROCID = $SLURM_PROCID"

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

    echo "Start time: " > start_time_network_connect.txt
    date >> start_time_network_connect.txt

    echo ">>> Init: "`date`
    # snudda init ${JOBDIR} --size ${SIMSIZE} --overwrite --randomseed 1234
    python setup_lateral.py
    
    if [ $? != 0 ]; then
	echo "Something went wrong during init, aborting!"
	ipcluster stop	
	exit -1
    fi


    echo ">>> Input: "`date`

    # snudda input ${JOBDIR} --parallel --time 5 --input input.json
    # Generate input for the reduced network...
    
    snudda input ${JOBDIR} --parallel --time 5 --input input.json --networkFile ${JOBDIR}/network-synapses-minimal.hdf5

    
    #.. Shut down cluster
    echo "Shutting down ipcontroller"

    python ../ipcontroller_shutdown.py


    date
    #echo "JOB END "`date` start_time_network_connect.txt

    echo "EXITING Dardel_runjob_lateral_input.sh"

fi
