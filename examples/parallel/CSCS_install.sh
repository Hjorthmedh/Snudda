#!/bin/bash

L=$VIRTUAL_ENV

module swap PrgEnv-cray PrgEnv-gnu
module load cray-python
module load daint-mc
module load gcc/9.3.0

export PATH=$LM/bin:$PATH
export LD_LIBRARY_PATH=$LM/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$L/lib/python3.8/site-packages:$PYTHONPATH

export CXX=CC
export CC=cc
export FC=ftn
export MPICC=cc
export MPICXX=CC

rm -rf $L/build
mkdir -pv $L/build
pushd $L/build

  # build parallel neuron with python interface
  mkdir neuron
  pushd neuron
    git clone -q https://github.com/neuronsimulator/nrn -b 8.0.0
    cd nrn
    ./build.sh
    autoreconf --force --install
    ./configure --prefix=$LM --exec-prefix=$LM \
      --with-mpi --with-nrnpython --with-paranrn \
      --without-x 
    make
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

  MPICC=cc pip install mpi4py --ignore-installed
  
popd
rm -rf $L/build


#conda deactivate
