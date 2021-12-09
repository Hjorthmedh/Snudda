#!/bin/bash

# Download and install miniconda3
# -- On Dardel compute nodes does not have wget,
# so you have to do this manually then
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh

module load snic-env
L=/cfs/klemming/home/${USER:0:1}/$USER/local/$SNIC_RESOURCE

./Miniconda3-latest-Linux-x86_64.sh -b -p $L/miniconda3

source activate_miniconda.sh
conda activate

conda update -n base conda -y
conda install wget -y
conda install git -y
conda install cmake -y
conda install bison -y
conda install pandoc -y
conda install flex -y
conda install ncurses -y
# conda install openmpi -y
conda update --all -y

# There is a bug in 3.3.2 which does not handle non-numeric host names correctly
# conda install mpich=3.2.1 -y
# Update, we use openmpi instead!

# This is needed to compile mpi4py -- is it really?
if [ $SNIC_RESOURCE == "tegner" ]; then
    module load gcc/9.2.0
    module load openmpi/4.1-gcc-9.2

elif [ $SNIC_RESOURCE == "beskow" ]; then
    echo "On Beskow"

    # Recompile mpi4py using MPICH
    export MPICC=cc
    export LD_LIBRARY_PATH=$MPICH_DIR/lib:$LD_LIBRARY_PATH

    # module load gcc/10.2.0
    # module load ??? # What is openmpi module on Beskow?
   #do something
elif [ $SNIC_RESOURCE == "dardel" ]; then
    echo "On Beskow"

    # Recompile mpi4py using MPICH
    export MPICC=cc
    export LD_LIBRARY_PATH=$MPICH_DIR/lib:$LD_LIBRARY_PATH


else
    echo "Unknown system $SNIC_RESOURCE"
fi

pip install mpi4py --ignore-installed --no-cache-dir




pushd ../../../

# Install Snudda -- only if you do not already have Snudda installed
# cd $L
# git clone git@github.com:Hjorthmedh/Snudda.git
# cd Snudda

pip install -r requirements.txt --no-cache-dir

# Dev installation using local copy
pip install -e .[dev]

pip --version

popd

