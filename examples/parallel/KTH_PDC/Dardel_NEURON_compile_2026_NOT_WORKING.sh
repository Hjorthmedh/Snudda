#!/bin/bash
set -euo pipefail

# =====================================================================
# Build NEURON 8.2.2 on Dardel using an EXISTING venv (no Miniconda).
#
# Usage:
#   ./Dardel_neuron_compile_2026.sh /path/to/your/venv /path/to/build/dir
#
# If you skip the args, sensible defaults are used below.
# =====================================================================

VENV="${1:?Usage: $0 /path/to/venv /path/to/build/dir}"
L="${2:-/cfs/klemming/projects/snic/<your-project>/$USER/neuron_build}"

export VENV=/cfs/klemming/projects/snic/snic2021-5-492/snudda_env
export L=/cfs/klemming/projects/snic/snic2021-5-492/NEURON/

echo "Using venv:       $VENV"
echo "Build/install dir: $L"
mkdir -pv "$L"

# ---------------------------------------------------------------------
# Activate your existing venv (created beforehand with: python3 -m venv $VENV)
# ---------------------------------------------------------------------
source "$VENV/bin/activate"

# ---------------------------------------------------------------------
# Cray / PDC programming environment
# (module names/versions can change between PDC software stack updates -
#  run `module avail PDC` / `module avail cpeGNU` if this fails)
# ---------------------------------------------------------------------
module load PDC
# module load PDC/24.11        # alternative: pin an explicit CPE release
# module load cpeGNU/24.11     # -> auto-swaps to PrgEnv-gnu for you

module swap PrgEnv-cray PrgEnv-gnu
module unload cray-libsci atp
module load cray-mpich
# module load cray-mpich-abi
module load cmake/3.22.3       # bump if a newer cmake module is available

export CRAYPE_LINK_TYPE=dynamic
export CRAY_ROOTFS=DSL

export CC=cc
export CXX=CC
export FC=ftn
export MPICC=cc
export MPICXX=CC

# ---------------------------------------------------------------------
# A venv doesn't ship its own libpython.so - it points back at the
# Python that created it. Find that real lib and stage the same
# "cray_libs" workaround the old script used for readline/ncurses.
# ---------------------------------------------------------------------
PYLIBDIR=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
PYSOFILE=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('LDLIBRARY'))")

export TMP0_DIR="$L/cray_libs"
mkdir -pv "$TMP0_DIR"
pushd "$TMP0_DIR"
  ln -sf /lib64/libncurses.so.6 libncurses.so
  ln -sf /lib64/libtinfo.so.6 libtinfo.so
  ln -sf /lib64/libreadline.so.7 libreadline.so
  ln -sf "$PYLIBDIR/$PYSOFILE" .
popd

export LD_LIBRARY_PATH="$TMP0_DIR:${MPICH_DIR:-}/lib:${LD_LIBRARY_PATH:-}"

export LIBRARY_PATH=/opt/cray/pe/python/3.12.12/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/cray/pe/python/3.12.12/lib:$LD_LIBRARY_PATH


# Remove any pip-installed neuron package so it doesn't shadow the build
pip uninstall neuron -y || true
pip install mpi4py --no-cache-dir --force-reinstall

# ---------------------------------------------------------------------
# Build NEURON
# ---------------------------------------------------------------------
export NRN_INSTALL_LOC="$L/neuron"

pushd "$L"
  # Clone manually first if this hangs on a compute node (known Beskow/Dardel quirk)
  if [ ! -d nrn ]; then
    git clone https://github.com/neuronsimulator/nrn
  fi

  cd nrn
  rm -rf build
  mkdir build
  cd build

  cmake .. \
    -DCMAKE_SKIP_RPATH:BOOL=YES \
    -DNRN_ENABLE_INTERVIEWS=OFF \
    -DNRN_ENABLE_PYTHON=ON \
    -DNRN_ENABLE_MPI=ON \
    -DNRN_ENABLE_RX3D=OFF \
    -DCMAKE_INSTALL_PREFIX="$NRN_INSTALL_LOC" \
    -DNRN_ENABLE_BINARY_SPECIAL=ON \
    -DNRN_ENABLE_CORENEURON=OFF \
    -DCMAKE_BUILD_TYPE:STRING=Release \
    -DPYTHON_EXECUTABLE="$(which python3)" \
    -DREADLINE_LIBRARY:FILEPATH="$TMP0_DIR/libreadline.so" \
    -DCURSES_EXTRA_LIBRARY:FILEPATH="$TMP0_DIR/libtinfo.so" \
    -DCURSES_CURSES_LIBRARY:FILEPATH="$TMP0_DIR/libncurses.so" \
    -DCURSES_NCURSES_LIBRARY:FILEPATH="$TMP0_DIR/libncurses.so" \
    -DCURSES_INCLUDE_PATH:PATH=/usr/include/curses/

  cmake --build . --parallel 1 --target install \
    1>"$L/build_log_Release.txt" 2>"$L/build_error_Release.txt"

  CC --version
popd

export PATH="$NRN_INSTALL_LOC/bin:$PATH"
export PYTHONPATH="$NRN_INSTALL_LOC/lib/python:${PYTHONPATH:-}"

# ---------------------------------------------------------------------
# Extra Python packages (installed into the venv, not system-wide)
# ---------------------------------------------------------------------
pip install pypandoc --no-cache-dir
pip install ipyparallel --no-cache-dir
pip install bluepyopt --no-cache-dir

# ---------------------------------------------------------------------
# Dirty ldflags fix so nrnivmodl works (needs purge_ldflags.py from the
# original repo, copy it next to this script before running)
# ---------------------------------------------------------------------
if [ -f "$NRN_INSTALL_LOC/bin/nrnmech_makefile" ] && [ -f purge_ldflags.py ]; then
  mv "$NRN_INSTALL_LOC/bin/nrnmech_makefile" "$NRN_INSTALL_LOC/bin/nrnmech_makefile.original"
  python3 purge_ldflags.py "$NRN_INSTALL_LOC/bin/nrnmech_makefile.original" > "$NRN_INSTALL_LOC/bin/nrnmech_makefile"
  chmod u+x "$NRN_INSTALL_LOC/bin/nrnmech_makefile"
fi

echo ""
echo "Build finished. To use NEURON later (e.g. in job scripts), run:"
echo "  source $VENV/bin/activate"
echo "  export PATH=$NRN_INSTALL_LOC/bin:\$PATH"
echo "  export PYTHONPATH=$NRN_INSTALL_LOC/lib/python:\$PYTHONPATH"
echo "  export LD_LIBRARY_PATH=$TMP0_DIR:\$LD_LIBRARY_PATH"
