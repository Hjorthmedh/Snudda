#!/bin/bash

# Download and install miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh

L=$VIRTUAL_ENV


./Miniconda3-latest-Linux-x86_64.sh -b -p $L/miniconda3

source activate_miniconda_CSCS.sh
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

module load gcc/9.3.0

export MPICC=cc
export LD_LIBRARY_PATH=$MPICH_DIR/lib:$LD_LIBRARY_PATH
pip install mpi4py --ignore-installed




pushd ../../

# Install Snudda -- only if you do not already have Snudda installed
# cd $L
# git clone git@github.com:Hjorthmedh/Snudda.git
# cd Snudda

pip install -r requirements.txt

# Dev installation using local copy
pip install -e .[dev]

pip --version

popd

