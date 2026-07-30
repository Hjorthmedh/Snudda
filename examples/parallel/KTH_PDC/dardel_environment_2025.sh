# This sets up the environment for Dardel, including compiling NEURON mod files

# IF $NO_ENGINES is set to 1, then ipyparallel setup is skipped
# IF $NO_SIM is set to 1, then NEURON modules are not compiled



# Will use N_ENGINES if set by user, otherwise default to SLURM_NTASKS - 2
N_ENGINES=${N_ENGINES:-$((SLURM_NTASKS-2))}
# let N_ENGINES="100"


# Make sure this is right python version
export SNUDDA_DIR=$HOME/Snudda
export PYTHONPATH=$SNUDDA_DIR/snudda_env/lib/python3.11/

# This is needed for NEURON
unset DISPLAY

ulimit -s unlimited

module load snic-env

module load cray-python
module swap PrgEnv-cray PrgEnv-gnu
module load cray-mpich-abi
module unload cray-libsci

source $HOME/Snudda/snudda_env/bin/activate

if [ -z "$NO_ENGINES" ]; then

    # IPYTHONDIR defaults to current folder, if network_path is not set
    export IPYTHONDIR="/cfs/klemming/scratch/${USER:0:1}/$USER/.ipython-${SLURM_JOB_ID}"
    export IPYTHON_PROFILE=default
    
    # Start the ipcontroller
    export FI_CXI_DEFAULT_VNI=$(od -vAn -N4 -tu < /dev/urandom)
    srun -n 1 -N 1 -c 2 --exact --overlap --mem=0 ./ipcontroller_new.sh &

    echo ">>> waiting 60s for controller to start"
    sleep 60 

    #.. Read in CONTROLLERIP
    CONTROLLERIP=$(<controller_ip.txt)

    echo ">>> starting ${N_ENGINES} engines "

    export FI_CXI_DEFAULT_VNI=$(od -vAn -N4 -tu < /dev/urandom)
    srun -n ${N_ENGINES} -c 2 -N ${SLURM_JOB_NUM_NODES} --exact --overlap --mem=0 ipengine \
	 --location=${CONTROLLERIP} --profile=${IPYTHON_PROFILE} --mpi \
	 --ipython-dir=${IPYTHONDIR}  --timeout=30.0 c.EngineFactory.max_heartbeat_misses=10  c.MPI.use='mpi4py' \
	 1> ipe_${SLURM_JOBID}.out 2> ipe_${SLURM_JOBID}.err &


    echo ">>> waiting 60s for engines to start"
    sleep 30

fi
    

if [ -z "$NO_SIM" ]; then

    SPECIAL_PATH=x86_64/special

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAY_LD_LIBRARY_PATH

    # Explicitly adding MPI paths (suggested by ChatGTP, why are they not loaded automatically?)
    export LD_LIBRARY_PATH=/opt/cray/pe/mpich/8.1.31/ofi/gnu/12.3/lib:$LD_LIBRARY_PATH
    export MPI_LIB_NRN_PATH=/opt/cray/pe/mpich/8.1.31/ofi/gnu/12.3/lib

    # Clear old compilation
    rm -rf x86_64 2>/dev/null || true

    export CXX=CC
    export CC=cc
    export FC=ftn
    export MPICC=cc
    export MPICXX=CC

    CC --version

    export FI_CXI_DEFAULT_VNI=$(od -vAn -N4 -tu < /dev/urandom)
    srun -n 1 nrnivmodl -incflags "-lltdl=/usr/lib64/libltdl.so.7 -lreadline=/lib64/libreadline.so.7 -lncurses=/lib64/libncurses.so.6.1" -loadflags "-DLTDL_LIBRARY=/usr/lib64/libltdl.so.7 -DREADLINE_LIBRARY=/lib64/libreadline.so.7 -DNCURSES_LIBRARY=/lib64/libncurses.so.6.1" mechanisms/

fi
