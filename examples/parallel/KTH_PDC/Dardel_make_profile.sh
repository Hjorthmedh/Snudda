#!/bin/bash

# You must pass the name of the profile to be created as first argument!

module load snic-env

export IPYTHONDIR="/cfs/klemming/scratch/${USER:0:1}/$USER/.ipython"
export profile_name=$1
source $HOME/Snudda/snudda_env/bin/activate

# Delete old profile
rm -r $IPYTHONDIR/profile_$profile_name

# Create new profile
ipython profile create --parallel --profile=$profile_name

#.. Obtain infiniband ip - this is faster and also internal
hostname -i > controller_ip.txt
CONTROLLERIP=$(<controller_ip.txt)

# Update content of new profile
export replace_str1="s/# c.IPController.registration_timeout = 0/c.IPController.registration_timeout = 300/"
export replace_str2="s/# c.HeartMonitor.max_heartmonitor_misses = 10/c.HeartMonitor.max_heartmonitor_misses = 10000/"
export replace_str3="s/# c.IPController.ip = '127.0.0.1'/c.IPController.ip = '$CONTROLLERIP'/"

sed -i "$replace_str1; $replace_str2; $replace_str3" $IPYTHONDIR/profile_$profile_name/ipcontroller_config.py

echo "Profile written and updated: $IPYTHONDIR/profile_$profile_name/ipcontroller_config.py"
