# Start an interactive node:
#
# salloc --nodes=1 -t 1:10:00 -A snic2022-5-245 --partition=main
#
# Then log on to the node:
#
# ssh $SLURM_NODELIST
#


./Miniconda_install.sh

# Clone the Flagser
git clone --recursive https://github.com/JasonPSmith/flagser-count.git

source activate_miniconda.sh
conda activate


# Create virtual environment
# python3 -m venv tutorial-env
# Activate flagser environment
# source flagser-env/bin/activate

# Setup PDC environemnt
module load snic-env
module load PDC

module swap PrgEnv-cray PrgEnv-gnu
module unload cray-libsci atp

export CRAYPE_LINK_TYPE=dynamic
export CRAY_ROOTFS=DSL

export PYTHONPATH=$L/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=$LN/lib/python:$LM/lib/python3.8/

export CXX=CC
export CC=cc
export FC=ftn
export MPICC=cc
export MPICXX=CC

module load CMake/3.21.2

cd flagser-count
make clean
make

pip install numpy

# TEST
(cd test && python run_test.py && cd ..)

# Finally install...
pip install .
