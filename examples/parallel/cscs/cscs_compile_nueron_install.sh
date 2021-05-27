#!/bin/bash

module swap PrgEnv-cray PrgEnv-gnu
module unload cray-libsci atp
export CRAYPE_LINK_TYPE=dynamic

L=$VIRTUAL_ENV/local/
LM=$L/miniconda3
LN=$L/neuron

mkdir -pv $L

export CXX=CC
export CC=cc
export FC=ftn
export MPICC=cc
export MPICXX=CC

conda activate

# neuron is also installed from requirements.txt, remove non-compatible version
pip uninstall neuron

# Is MINIC used?
export MINIC=$LM

export PATH=$LM/bin:$LN/bin:$PATH
# export LD_LIBRARY_PATH=$LN/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$MPICH_DIR/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$L/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=$LN/lib/python:$LM/lib/python3.8/

export NRN_INSTALL_LOC=$LN


# We need to recompile mpi4py to use mpich libraries of beskow
# UPDATE: This is now done in Miniconda_install.sh
# pip install mpi4py --ignore-installed

# Från MDJs gamla buildscript
# export CFLAGS="-dynamic -O3 -funroll-loops -march=corei7-avx -mavx  -ffast-math -DCACHEVEC=1" 

# Maybe git clone does not work on compute nodes?

# rm -rf $L/build
# mkdir -pv $L/build
pushd $L/build

  # build parallel neuron with python interface
  # mkdir neuron
  # pushd neuron
  # OLD: git clone -q https://github.com/nrnhines/nrn

  # You have to do git clone manually!!
  # git clone https://github.com/neuronsimulator/nrn -b 8.0.0
  
  cd nrn
  rm -r build
  mkdir build
  cd build

  echo "About to run"
  echo `which cmake`
#   cmake .. \
#   	  -DCMAKE_BUILD_TYPE:STRING=Debug \
#	  -DNRN_ENABLE_INTERVIEWS=OFF \
#	  -DNRN_ENABLE_PYTHON=ON \
#	  -DNRN_ENABLE_MPI=ON \
#	  -DNRN_ENABLE_RX3D=OFF \
#	  -DCMAKE_INSTALL_PREFIX=$LN \
#	  -DNRN_ENABLE_BINARY_SPECIAL=ON \
#	  -DPYTHON_EXECUTABLE=`which python3` \
#	  -DCMAKE_C_COMPILER:FILEPATH=cc \
#	  -DCMAKE_CXX_COMPILER:FILEPATH=CC \

#          -DCMAKE_C_FLAGS="-DDEBUG -g" \
#          -DCMAKE_CXX_FLAGS="-DDEBUG -g" \
#	  -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
	  # -DCMAKE_C_FLAGS="-mavx2" \
          # -DCMAKE_CXX_FLAGS="-mavx2" \
	  #	  -DNRN_ENABLE_CORENEURON=ON \

	  #-DCMAKE_BUILD_TYPE=Debug \
	  
  cmake .. \
      -DNRN_ENABLE_INTERVIEWS=OFF \
      -DNRN_ENABLE_PYTHON=OFF   \
      -DNRN_ENABLE_MPI=ON   \
      -DNRN_ENABLE_RX3D=OFF  \
      -DCMAKE_INSTALL_PREFIX=$NRN_INSTALL_LOC \
      -DNRN_ENABLE_BINARY_SPECIAL=ON \
      -DNRN_ENABLE_CORENEURON=OFF \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DCURSES_CURSES_LIBRARY:FILEPATH=$MINIC/lib/libncurses.so \
      -DCURSES_INCLUDE_PATH:PATH=$MINIC/include/ncurses.h

  cmake --build . \
	--parallel 1 \
	--target install 1>$L/build_log_Release.txt 2>$L/build_error_Release.txt

    
#    ./build.sh
#    autoreconf --force --install
#    ./configure --prefix=$LM --exec-prefix=$LM \
#      --with-mpi --with-nrnpython --with-paranrn \
#      --without-x \
#      --without-memacs --without-readline --disable-shared --without-nmodl \
#      linux_nrnmech=no always_call_mpi_init=yes java_dlopen=no \
#      --host=x86_64-unknown-linux-gnu --disable-pysetup --without-iv \
#      CC=cc CXX=CC MPICC=cc MPICXX=CC CFLAGS="$CFLAGS" CXXFLAGS="$CFLAGS" \
#      # -LLIBDIR

    echo "About to run"
    echo `which make`
	  
#    # make -j   # -j parallel compilation
#    make -j
#    make install

#    pushd src/nrnpython
#      python setup.py install
#    popd
    rm -r $L/share/nrn/{demo,examples}
  popd

  # install pypandoc
  # pip install pypandoc --install-option="--prefix=$L"
  # pip install pypandoc --prefix=$L
  pip install pypandoc

  # install ipyparallel
  # pip install ipyparallel --install-option="--prefix=$L"
  # pip install ipyparallel --prefix=$L
  pip install ipyparallel

  # install bluepyopt
  # pip install bluepyopt --install-option="--prefix=$L"
  # pip install bluepyopt --prefix=$L
  pip install bluepyopt

  pip install mpi4py

# popd
# rm -rf $L/build


conda deactivate
