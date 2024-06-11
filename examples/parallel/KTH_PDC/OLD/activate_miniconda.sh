##############################################################################
#
# The lines below should be added to .profile in the home directory on Tegner
#
##############################################################################

module load snic-env

# This path must match what is in Dardel_NEURON_isntall.sh
# L=/cfs/klemming/home/${USER:0:1}/$USER/local/$SNIC_RESOURCE
if [ $# -eq 0 ]
  then
    L=/cfs/klemming/projects/snic/snic2021-5-492/$USER/local/$SNIC_RESOURCE
    echo "No argument. Using path: $L"
else
    L=$1
    echo "Using path: $L"
fi


echo "Activating for path $L"

__conda_setup="$('$L/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"

if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "$L/miniconda3/etc/profile.d/conda.sh" ]; then
        . "$L/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="$L/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

##############################################################################
