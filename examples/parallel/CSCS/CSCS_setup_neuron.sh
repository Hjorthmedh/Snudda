module load cray-python
module load daint-mc
module swap PrgEnv-cray PrgEnv-gnu


module rm cray-python

# This is created in CSCS_setup.sh
source ~/snudda_env/bin/activate

L=~/snudda_env
LN=$L/neuron

export CXX=CC
export CC=cc
export FC=ftn
export MPICC=cc
export MPICXX=CC

pip uninstall neuron -y
export PATH=$L/bin:$LN/bin:$PATH

export LD_LIBRARY_PATH=$MPICH_DIR/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$L/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=$LN/lib/python:$PYTHONPATH

export NRN_INSTALL_LOC=$LN


pushd $L
git clone https://github.com/neuronsimulator/nrn -b 8.0.0
  
cd nrn
rm -r build
mkdir build
cd build

echo "About to run"
echo `which cmake`

cmake .. \
      -DNRN_ENABLE_INTERVIEWS=OFF \
      -DNRN_ENABLE_PYTHON=ON   \
      -DNRN_ENABLE_MPI=ON   \
      -DNRN_ENABLE_RX3D=OFF  \
      -DCMAKE_INSTALL_PREFIX=$NRN_INSTALL_LOC \
      -DNRN_ENABLE_BINARY_SPECIAL=ON \
      -DNRN_ENABLE_CORENEURON=OFF \
      -DCMAKE_BUILD_TYPE:STRING=Release \
#      -DCURSES_CURSES_LIBRARY:FILEPATH=$MINIC/lib/libncurses.so \
#      -DCURSES_INCLUDE_PATH:PATH=$MINIC/include/ncurses.h

#cmake --build . \
#	--parallel 1 \
#	--target install 1>$L/build_log_Release.txt 2>$L/build_error_Release.txt

cmake --build . \
      --target install 1>$L/build_log_Release.txt 2>$L/build_error_Release.txt


echo "About to run"
echo `which make`

rm -r $L/share/nrn/{demo,examples}
popd

pip install pypandoc --no-cache-dir

# Already installed:
# pip install ipyparallel --no-cache-dir
# pip install bluepyopt --no-cache-dir
# pip install mpi4py --no-cache-dir

deactivate



