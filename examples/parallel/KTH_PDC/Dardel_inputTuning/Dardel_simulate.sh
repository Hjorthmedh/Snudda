#!/bin/bash
#!/bin/bash -l
#SBATCH --partition=main
#SBATCH -t 2:00:00
#SBATCH -J snuddaSimulate
#SBATCH -A snic2021-5-492
#SBATCH --nodes=4
#SBATCH --tasks-per-node=128

NETWORK_DIR=networks/$1_1k_spn04

let N_WORKERS="$SLURM_NNODES * 128"
echo "Number of workers is = $SLURM_NNODES"


if [[ -d "bgd01/Parkinson/20211105-PD0/$1" ]]; then
    export SNUDDA_DATA="bgd01/Parkinson/20211105-PD0/$1"
    echo "Setting SNUDDA_DATA to $SNUDDA_DATA"
else
    echo "SNUDDA_DATA environment variable not changed (may be empty): $SNUDDA_DATA"
fi

# Synapse file
SNUDDA_DIR=/cfs/klemming/scratch/${USER:0:1}/$USER/Projects/SnuddaProj10/Snudda/snudda

NETWORK_INFO_FILE=$NETWORK_DIR/network-synapses.hdf5
NETWORK_INPUT_FILE=$NETWORK_DIR/input-spikes.hdf5
NETWORK_VOLTAGE_FILE=$NETWORK_DIR/simulation/voltage-trace-${SLURM_JOBID}.txt

echo "Network dir: "$NETWORK_DIR


##############
echo "Loading environment"
source ../../../../../snudda_env/bin/activate
module load cray-python
module load cray-mpich-abi

module load snic-env
# --- I have tried with the gnu compiler, and also with the cray compiler
module swap PrgEnv-cray PrgEnv-gnu
module unload cray-libsci atp
export CRAYPE_LINK_TYPE=dynamic
export CRAY_ROOTFS=DSL

########## Needed for NEURON
echo "Setting paths needed for Neuron"
L=/cfs/klemming/scratch/${USER:0:1}/$USER/Projects/SnuddaProj10/snudda_env
LN=$L/neuron
export PATH=$L/bin:$LN/bin:$PATH
export LD_LIBRARY_PATH=$MPICH_DIR/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$L/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=$LN/lib/python:$PYTHONPATH

########### Needed for NEURON end.


# Remove "special" old directory
# This is now done by core.py automatically.

#echo "Compiling mechanisms"
#rm mechanisms
##ln -s ../../snudda/data/neurons/mechanisms
#ln -s BasalGangliaData/data/neurons/mechanisms
#rm -r x86_64
#nrnivmodl mechanisms


# GJ disabled
# srun -n $N_WORKERS $SNUDDA_DIR/../examples/parallel/x86_64/special -mpi -python $SNUDDA_DIR/simulate/simulate.py $NETWORK_INFO_FILE $NETWORK_INPUT_FILE --disableGJ --time 3.5 --voltOut $NETWORK_VOLTAGE_FILE

# GJ active
srun -n $N_WORKERS $SNUDDA_DIR/../examples/parallel/KTH_PDC/Dardel_inputTuning/x86_64/special -mpi -python $SNUDDA_DIR/simulate/simulate.py $NETWORK_INFO_FILE $NETWORK_INPUT_FILE --time 10.0
#srun -n $N_WORKERS $SNUDDA_DIR/../examples/parallel/x86_64/special -mpi -python $SNUDDA_DIR/simulate/simulate.py $NETWORK_INFO_FILE $NETWORK_INPUT_FILE --time 3.5
#srun -n 1 $SNUDDA_DIR/../examples/parallel/x86_64/special -mpi -python $SNUDDA_DIR/simulate/simulate.py $NETWORK_INFO_FILE $NETWORK_INPUT_FILE --time 3.5 --voltOut $NETWORK_VOLTAGE_FILE

# Original version, start neuron from python, does not work on beskow
#aprun -n  $N_WORKERS /cfs/klemming/nobackup/h/hjorth/ChINopt/model/x86_64/special -mpi -python snudda_simulate.py /cfs/klemming/nobackup/h/hjorth/ChINopt/model/save/save/network-connect-synapse-file-1749867.pickle

echo "JOB END "`date` 
