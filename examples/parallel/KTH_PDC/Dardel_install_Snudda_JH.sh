#!/bin/bash

if [ -d "$HOME/Snudda" ]; then
    echo "Directory $HOME/Snudda already exists"
else
    pushd $HOME
    git clone -b master https://github.com/Hjorthmedh/Snudda.git
    # git clone -b master git@github.com:Hjorthmedh/Snudda.git    
fi

# If you are low on disk space on $HOME, then you could try
# using $SNIC_TMP/Snudda as location instead

pushd $HOME/Snudda
git pull

module load snic-env
module load cray-python
module swap PrgEnv-cray PrgEnv-gnu
module load cray-mpich-abi
module unload cray-libsci

# Setup virtual environment
python -m venv snudda_env

# Move env files to secondary storage location
mv snudda_env /cfs/klemming/projects/supr/snic2021-5-492/
ln -s /cfs/klemming/projects/supr/snic2021-5-492/snudda_env

source snudda_env/bin/activate

pip install --upgrade pip
MPICC=cc pip install mpi4py

pip install wheel
pip install -r requirements.txt
pip install neuron --upgrade
pip install -e .[dev]






