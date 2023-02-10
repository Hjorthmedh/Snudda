#!/bin/bash -e


	
START=$(date +%s)

echo "Start time: " > start_time_network_connect.txt
date >> start_time_network_connect.txt

#.. Obtain infiniband ip - this is faster and also internal
hostname -i > controller_ip.txt
CONTROLLERIP=$(<controller_ip.txt)


#.. Some old NEURON stuff used in previous tickets
# rm -rf x86_64
# nrnivmodl ./mechanisms >/dev/null

echo ">>>>>>> Starting ipcontroller"

#.. Create the ipython directory

ipython profile create --ipython-dir=${IPYTHONDIR}

ipcontroller --ip=${CONTROLLERIP} --location=${CONTROLLERIP} --profile=${IPYTHON_PROFILE}  \
--HeartMonitor.period=10000 --HeartMonitor.max_heartmonitor_misses=1000 --HubFactory.registration_timeout=60 \
--debug --sqlitedb c.TaskScheduler.scheme = 'pure'  c.TaskScheduler.hwm=0 1>ipc_${SLURM_JOBID}.out 2> ipc_${SLURM_JOBID}.err &


echo ">>> waiting 60s for controller to start"
sleep 60 

#.. Start the engines
echo ">>> starting ${IPNWORKERS} engines "
srun -n ${IPNWORKERS} ipengine --location=${CONTROLLERIP} --profile=${IPYTHON_PROFILE} --mpi \
--ipython-dir=${IPYTHONDIR}  --timeout=30.0 --log-level=DEBUG \
--BaseParallelApplication.verbose_crash=True --IPEngine.verbose_crash=True \
--Kernel.stop_on_error_timeout=1.0 --IPythonKernel.stop_on_error_timeout=1.0 \
Session.buffer_threshold=4096 Session.copy_threshold=250000 \
Session.digest_history_size=250000 c.EngineFactory.max_heartbeat_misses=10  c.MPI.use='mpi4py' \
1> ipe_${SLURM_JOBID}.out 2> ipe_${SLURM_JOBID}.err &


echo ">>> waiting 60s for engines to start"
sleep 60


#.. This is where you run your python script that utilizes the ipyparallel network
echo ">>> Python script"
#python optimize.py 1> python_${SLURM_JOBID}.out 2> python_${SLURM_JOBID}.err

echo ">>> Ending python script, sleep 100s"
#sleep 100

echo "Ending ipcontroller.sh"



