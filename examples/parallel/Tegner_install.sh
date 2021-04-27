#!/bin/bash

source ~/.profile

module load snic-env

echo "This assumes you have already installed miniconda3 on Tegner"
echo ""
echo "https://github.com/Hjorthmedh/Snudda/wiki/D.-Setting-up-Tegner-@%C2%A0KTH"
echo ""

L=/cfs/klemming/nobackup/${USER:0:1}/$USER/local/$SNIC_RESOURCE
LM=$L/miniconda3

mkdir -pv $L

module load gcc/9.2.0
module load openmpi/4.1-gcc-9.2

conda activate

export PATH=$LM/bin:$PATH
export LD_LIBRARY_PATH=$LM/lib:$LD_LIBRARY_PATH
# export PYTHONPATH=$L/lib/python3.8/site-packages:$PYTHONPATH

rm -rf $L/build
mkdir -pv $L/build
pushd $L/build

  # build parallel neuron with python interface
  mkdir neuron
  pushd neuron
    git clone -q https://github.com/neuronsimulator/nrn
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

popd
rm -rf $L/build


conda deactivate
