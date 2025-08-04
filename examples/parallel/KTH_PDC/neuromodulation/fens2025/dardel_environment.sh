# This sets up the environment for Dardel, including compiling NEURON mod files


# Make sure this is right python version
export SNUDDA_DIR=$HOME/Snudda
export PYTHONPATH=$SNUDDA_DIR/snudda_env/lib/python3.11/

# This is needed for NEURON
unset DISPLAY

module load snic-env

# Here we are just running in serial for Snudda network generation 
# export IPYTHONDIR="$network_path/.ipython"

module load cray-python
module swap PrgEnv-cray PrgEnv-gnu
module load cray-mpich-abi
module unload cray-libsci

source $HOME/Snudda/snudda_env/bin/activate

SPECIAL_PATH=x86_64/special

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAY_LD_LIBRARY_PATH

# Explicitly adding MPI paths (suggested by ChatGTP, why are they not loaded automatically?)
export LD_LIBRARY_PATH=/opt/cray/pe/mpich/8.1.31/ofi/gnu/12.3/lib:$LD_LIBRARY_PATH
export MPI_LIB_NRN_PATH=/opt/cray/pe/mpich/8.1.31/ofi/gnu/12.3/lib

# Clear old compilation
rm -r x86_64

export CXX=CC
export CC=cc
export FC=ftn
export MPICC=cc
export MPICXX=CC

CC --version

export FI_CXI_DEFAULT_VNI=$(od -vAn -N4 -tu < /dev/urandom)
srun -n 1 nrnivmodl -incflags "-lltdl=/usr/lib64/libltdl.so.7 -lreadline=/lib64/libreadline.so.7 -lncurses=/lib64/libncurses.so.6.1" -loadflags "-DLTDL_LIBRARY=/usr/lib64/libltdl.so.7 -DREADLINE_LIBRARY=/lib64/libreadline.so.7 -DNCURSES_LIBRARY=/lib64/libncurses.so.6.1" mechanisms/
