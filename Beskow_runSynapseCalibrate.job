#!/bin/bash -l
#SBATCH -t 0:59:00
#SBATCH --time-min=0:59:00
#SBATCH -J snudda_simulate
#SBATCH -A 2019-3-644
#SBATCH -o save/output-snudda_simulate.o%j
#SBATCH -e save/error-snudda_simulate.e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

HOME=/cfs/klemming/nobackup/"${USER:0:1}"/$USER/Snudda

module swap PrgEnv-cray PrgEnv-intel/6.0.5
module load craype-haswell
module unload cray-libsci atp
module load neuron/7.5-py37

export CRAYPE_LINK_TYPE=dynamic
export CRAY_ROOTFS=DSL

# Remove "special" old directory
rm -r x86_64
nrnivmodl cellspecs/mechanisms

SIMNAME=networks/BeskowSynapseCalibration.${SLURM_JOBID}

python3 snudda_init_custom.py $simName
./snudda.py place $simName 
./snudda.py detect $simName
./snudda.py prune $simName

python3 snudda_cut.py $simName/network-pruned-synapses.hdf5 "abs(z)<100e-6"


srun -n $N_WORKERS /cfs/klemming/nobackup/"${USER:0:1}"/$USER/Snudda/x86_64/special -mpi -python snudda/utils/network_pair_recording_simulation.py run $simName/network-cut-slice.hdf5 dSPN iSPN

srun -n $N_WORKERS /cfs/klemming/nobackup/"${USER:0:1}"/$USER/Snudda/x86_64/special -mpi -python snudda/utils/network_pair_recording_simulation.py run $simName/network-cut-slice.hdf5 iSPN dSPN

srun -n $N_WORKERS /cfs/klemming/nobackup/"${USER:0:1}"/$USER/Snudda/x86_64/special -mpi -python snudda/utils/network_pair_recording_simulation.py run $simName/network-cut-slice.hdf5 FS dSPN

srun -n $N_WORKERS /cfs/klemming/nobackup/"${USER:0:1}"/$USER/Snudda/x86_64/special -mpi -python snudda/utils/network_pair_recording_simulation.py run $simName/network-cut-slice.hdf5 FS iSPN


