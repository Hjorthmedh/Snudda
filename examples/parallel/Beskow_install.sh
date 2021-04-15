#!/bin/bash

source activate_miniconda.txt

module load snic-env

# --- I have tried with the gnu compiler, and also with the cray compiler
module swap PrgEnv-cray PrgEnv-gnu
module unload cray-libsci atp
export CRAYPE_LINK_TYPE=dynamic

# # module load craype-ivybridge

echo "This assumes you have already installed miniconda3 on Tegner"
echo ""
echo "https://github.com/Hjorthmedh/Snudda/wiki/D.-Setting-up-Tegner-@%C2%A0KTH"
echo ""

L=/cfs/klemming/nobackup/${USER:0:1}/$USER/local/$SNIC_RESOURCE
LM=$L/miniconda3
LN=$L/neuron

mkdir -pv $L

export CXX=CC
export CC=cc
export FC=ftn
export MPICC=cc
export MPICXX=CC

conda activate

export PATH=$LM/bin:$LN/bin:$PATH
export LD_LIBRARY_PATH=$LN/lib:$LD_LIBRARY_PATH
# export PYTHONPATH=$L/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=$LN/lib/python:$LM/lib/python3.8/



# Från MDJs gamla buildscript
# export CFLAGS="-dynamic -O3 -funroll-loops -march=corei7-avx -mavx  -ffast-math -DCACHEVEC=1" 

rm -rf $L/build
mkdir -pv $L/build
pushd $L/build

  # build parallel neuron with python interface
  mkdir neuron
  pushd neuron
  # git clone -q https://github.com/nrnhines/nrn
  git clone https://github.com/neuronsimulator/nrn -b 7.8.2
  
  cd nrn
  mkdir build
  cd build

  echo "About to run"
  echo `which cmake`
  cmake .. \
	  -DCMAKE_BUILD_TYPE=Debug \
	  -DNRN_ENABLE_INTERVIEWS=OFF \
	  -DNRN_ENABLE_PYTHON=ON \
	  -DNRN_ENABLE_MPI=ON \
	  -DNRN_ENABLE_RX3D=OFF \
	  -DCMAKE_INSTALL_PREFIX=$LN \
	  -DNRN_ENABLE_BINARY_SPECIAL=ON \
	  -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
	  -DPYTHON_EXECUTABLE=`which python3` \
	  -DCMAKE_C_COMPILER:FILEPATH=cc \
	  -DCMAKE_CXX_COMPILER:FILEPATH=CC \
#	  -DNRN_ENABLE_CORENEURON=ON \




    
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
	  
    make -j
    make install
    pushd src/nrnpython
      # python setup.py install --prefix=$L
      python setup.py install
    popd
    rm -r $L/share/nrn/{demo,examples}
  popd

  # install pandoc
  version=2.1.1
  rm -rf pandoc* $LM/{bin,lib}/pandoc*
  wget https://github.com/jgm/pandoc/releases/download/$version/pandoc-$version-linux.tar.gz
  tar -C $LM/lib -xvf pandoc-$version-linux.tar.gz pandoc-$version/bin/pandoc{,-citeproc}
  ln -sv ../lib/pandoc-$version/bin/pandoc{,-citeproc} $LM/bin/
  rm -fv pandoc*

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

popd
# rm -rf $L/build


conda deactivate
