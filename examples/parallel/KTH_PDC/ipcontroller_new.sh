#!/bin/bash -e

ml valgrind4hpc
#ml craype-hugepages1G

	
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

#valgrind --leak-check=full ipcontroller --ip=${CONTROLLERIP} --location=${CONTROLLERIP} --profile=${IPYTHON_PROFILE}  \
#--HeartMonitor.period=10000 --HeartMonitor.max_heartmonitor_misses=10 --HubFactory.registration_timeout=60 \
#--debug --sqlitedb c.TaskScheduler.scheme = 'pure'  c.TaskScheduler.hwm=0 1>ipc_${SLURM_JOBID}.out 2> ipc_${SLURM_JOBID}.err 

ipcontroller --ip=${CONTROLLERIP} --location=${CONTROLLERIP} --profile=${IPYTHON_PROFILE}  \
--HeartMonitor.period=10000 --HeartMonitor.max_heartmonitor_misses=10 --HubFactory.registration_timeout=60 \
--sqlitedb c.TaskScheduler.scheme = 'pure'  c.TaskScheduler.hwm=0 1>ipc_${SLURM_JOBID}.out 2> ipc_${SLURM_JOBID}.err 

#ipcontroller --ip=${CONTROLLERIP} --location=${CONTROLLERIP} --profile=${IPYTHON_PROFILE}  \
#--HeartMonitor.period=10000 --HeartMonitor.max_heartmonitor_misses=10 --HubFactory.registration_timeout=60 \
#--nodb c.TaskScheduler.scheme = 'pure'  c.TaskScheduler.hwm=0 1>ipc_${SLURM_JOBID}.out 2> ipc_${SLURM_JOBID}.err 


echo "Ending ipcontroller.sh"



