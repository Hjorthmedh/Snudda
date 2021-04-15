##############################################################################
#
# The lines below should be added to .profile in the home directory on Tegner
#
##############################################################################

module load snic-env
L=/cfs/klemming/nobackup/${USER:0:1}/$USER/local/$SNIC_RESOURCE

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
