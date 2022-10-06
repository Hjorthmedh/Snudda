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

module load cray-python
module load snic-env
module swap PrgEnv-cray PrgEnv-gnu
module unload cray-libsci atp

# Setup virtual environment
python -m venv snudda_env
source snudda_env/bin/activate

pip install --upgrade pip
MPICC=cc pip install mpi4py

pip install -r requirements.txt
pip install -e .[dev]






