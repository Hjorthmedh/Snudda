
module load cray-python
module load daint-mc
module swap PrgEnv-cray PrgEnv-gnu

pushd ~
python3 -m venv snudda_env

module rm cray-python

source snudda_env/bin/activate

MPICC=cc pip install mpi4py
pip install snudda
pip install numba
# salloc -C mc -A ich030 -n 1 -t 1:00:00
# srun -C mc -A ich030 -n 1 --ntasks-per-node=1 ./daint-snudda-venv-build.sh


