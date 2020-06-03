#!/bin/bash -e

# load anaconda with dependencies (modules gcc, openmpi, etc)
module load anaconda/py37/5.0.1

# make sure there are no conflicts
module list

HOST=$(hostname -s | cut -d'-' -f1)

# setup directory for local software
L=/cfs/klemming/nobackup/${USER:0:1}/$USER/local/$HOST
mkdir -pv $L

export PATH=$L/bin:$PATH
export LD_LIBRARY_PATH=$L/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$L/lib/python3.7/site-packages:$PYTHONPATH

rm -rf $L/build
mkdir -pv $L/build
pushd $L/build

  # build parallel neuron with python interface
  mkdir neuron
  pushd neuron
    git clone -q https://github.com/nrnhines/nrn
    cd nrn
    ./build.sh
    autoreconf --force --install
    ./configure --prefix=$L --exec-prefix=$L \
      --with-mpi --with-nrnpython --with-paranrn \
      --without-x 
    make
    make install
    pushd src/nrnpython
      python setup.py install --prefix=$L
    popd
    rm -r $L/share/nrn/{demo,examples}
  popd

  # install pandoc
  version=2.1.1
  rm -rf pandoc* $L/{bin,lib}/pandoc*
  wget https://github.com/jgm/pandoc/releases/download/$version/pandoc-$version-linux.tar.gz
  tar -C $L/lib -xvf pandoc-$version-linux.tar.gz pandoc-$version/bin/pandoc{,-citeproc}
  ln -sv ../lib/pandoc-$version/bin/pandoc{,-citeproc} $L/bin/
  rm -fv pandoc*

  # install pypandoc
  pip install pypandoc --install-option="--prefix=$L"

  # install ipyparallel
  pip install ipyparallel --install-option="--prefix=$L"

  # install bluepyopt
  pip install bluepyopt --install-option="--prefix=$L"

popd
rm -rf $L/build

echo "For dev installation, in Snudda repository's root directory, run:"
echo ""
echo "module load anaconda/py37/5.0.1"
echo ""
echo "Make sure PYTHONPATH points correctly, then run"
echo "PYTHONPATH=/cfs/klemming/nobackup/"${USER:0:1}"/$USER/local/"$HOST"/lib/python3.7/site-packages"
echo ""
echo "pip install -e .[dev] --prefix /cfs/klemming/nobackup/"${USER:0:1}"/$USER/local/"$HOST"/"
