# This sets up the environment for Dardel, including compiling NEURON mod files

# IF $NO_ENGINES is set to 1, then ipyparallel setup is skipped
# IF $NO_SIM is set to 1, then NEURON modules are not compiled

N_ENGINES=${N_ENGINES:-${SLURM_NTASKS:-10}}
# let N_ENGINES="100"


# In case we run this interactively, then no SLURM_JOB_NUM_NODES
# Set SLURM_JOB_NUM_NODES to 1 if it is unset, empty, or <= 0
if [ "${SLURM_JOB_NUM_NODES:-0}" -le 0 ]; then
    SLURM_JOB_NUM_NODES=1
fi

# Make sure this is right python version
export SNUDDA_DIR=$HOME/Snudda
# export PYTHONPATH=$SNUDDA_DIR/snudda_env/lib/python3.12/:$PYTHONPATH
unset PYTHONPATH

# This is needed for NEURON
unset DISPLAY

ulimit -s unlimited

# module load snic-env
# module load cray-python
# module swap PrgEnv-cray PrgEnv-gnu
# module load cray-mpich-abi
# module unload cray-libsci

module load PDC
module load neuron

source $HOME/Snudda/snudda_env/bin/activate

# We no longer need to start the engines, that is done by snudda itself

# Only compile NEURON MODules if user sets COMPILE_MOD=1

# We set these, since RxD might ask to compile...
export CC=cc
export CXX=CC
export FC=ftn
export MPICC=cc
export MPICXX=CC

if [  -n "$COMPILE_MOD"  ]; then

    # Only recompile if the MOD files are missing
    if [ ! -d "x86_64" ]; then

      SPECIAL_PATH=x86_64/special

      # export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAY_LD_LIBRARY_PATH

      # Explicitly adding MPI paths (suggested by ChatGTP, why are they not loaded automatically?)
      # export LD_LIBRARY_PATH=/opt/cray/pe/mpich/8.1.31/ofi/gnu/12.3/lib:$LD_LIBRARY_PATH
      # export MPI_LIB_NRN_PATH=/opt/cray/pe/mpich/8.1.31/ofi/gnu/12.3/lib

      # Clear old compilation
      rm -rf x86_64 2>/dev/null || true

      export CXX=CC
      export CC=cc
      export FC=ftn
      export MPICC=cc
      export MPICXX=CC

      CC --version

      export FI_CXI_DEFAULT_VNI=$(od -vAn -N4 -tu < /dev/urandom)
      # srun -n 1 nrnivmodl -incflags "-lltdl=/usr/lib64/libltdl.so.7 -lreadline=/lib64/libreadline.so.7 -lncurses=/lib64/libncurses.so.6.1" -loadflags "-DLTDL_LIBRARY=/usr/lib64/libltdl.so.7 -DREADLINE_LIBRARY=/lib64/libreadline.so.7 -DNCURSES_LIBRARY=/lib64/libncurses.so.6.1" mechanisms/

      srun -n 1 nrnivmodl mechanisms/

    fi

fi
