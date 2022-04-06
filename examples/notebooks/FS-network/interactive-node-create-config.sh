module load snic-env

export IPYTHONDIR="/cfs/klemming/scratch/${USER:0:1}/$USER/.ipython"
rm -r $IPYTHONDIR
export IPYTHON_PROFILE=default
source $HOME/Snudda/snudda_env/bin/activate
       
if [[ -d "$HOME/BasalGangliaData/data" ]]; then
    export SNUDDA_DATA="$HOME/BasalGangliaData/data"
    echo "Setting SNUDDA_DATA to $SNUDDA_DATA"
else
    echo "SNUDDA_DATA environment variable not changed (may be empty): $SNUDDA_DATA"
fi

mkdir FS_network_2-cur-inj
srun python3 FS_network_make_config_2-cur-inj.py
