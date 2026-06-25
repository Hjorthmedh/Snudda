export IPYTHONDIR="`pwd`/.ipython"
export IPYTHON_PROFILE=default

numNeurons=1000
simNamePart='1k_testing'
export PD_DATA_DIR="../../../BasalGangliaData/Parkinson/PD2026"

ipcluster start -n 5 --profile=$IPYTHON_PROFILE --ip=127.0.0.1&
sleep 20

######################################################################
# 1. Wild-type (WT) network: built from scratch
######################################################################
network_path_wt=networks/WT_$simNamePart

# Note that we need to start ipcluster after the SNUDDA_DATA environment
# variable is set, otherwise the workers will not know of it.
snudda init $network_path_wt --size $numNeurons --overwrite \
    --connectionFile $PD_DATA_DIR/WT/connectivity/striatum/striatum-connectivity.json \
    --seed 1234 --snudda_data "$PD_DATA_DIR/WT"
snudda place  $network_path_wt --parallel
snudda detect $network_path_wt --parallel
snudda prune  $network_path_wt --parallel

#python3 ../../snudda/analyse/analyse_striatum.py $network_path_wt/network-synapses.hdf5
######################################################################
# 2. PD network: copy the WT placement and swap morphologies/models
######################################################################
network_path_pd_ref=networks/PD_ref_$simNamePart

snudda init $network_path_pd_ref --size $numNeurons --overwrite \
    --connectionFile $PD_DATA_DIR/WT/connectivity/striatum/striatum-connectivity.json \
    --seed 1234 --snudda_data "$PD_DATA_DIR/PD"
    
#The following is necessary to make sure that the vars match between WT and PD
python3 -m snudda.utils.force_morphology_variations \
    "$network_path_pd_ref" "$network_path_wt"

snudda detect $network_path_pd_ref --parallel
snudda prune  $network_path_pd_ref --parallel

#python3 ../../snudda/analyse/analyse_striatum.py $network_path_pd_ref/network-synapses.hdf5
######################################################################
# 3. Swapped PD network to retain the position
#    of the "surviving" synapses
######################################################################
network_path_pd=networks/PD_$simNamePart
mkdir -p "$network_path_pd"

remap_removed_input=false
remapped_fraction=0.0

NETWORK_PATH_WT="$network_path_wt" \
NETWORK_PATH_PD="$network_path_pd" \
SNUDDA_DATA_WT="$PD_DATA_DIR/WT" \
SNUDDA_DATA_PD="$PD_DATA_DIR/PD" \
REMAP_REMOVED_INPUT="$remap_removed_input" \
REMAPPED_FRACTION="$remapped_fraction" \
python3 - <<'PY'
import os
from snudda.utils.swap_to_degenerated_morphologies import SwapToDegeneratedMorphologies

network_path_wt = os.environ["NETWORK_PATH_WT"]
network_path_pd = os.environ["NETWORK_PATH_PD"]
snudda_data_wt  = os.environ["SNUDDA_DATA_WT"]
snudda_data_pd  = os.environ["SNUDDA_DATA_PD"]
remap_removed_input = os.environ["REMAP_REMOVED_INPUT"].strip().lower() in ("1", "true", "yes")
remapped_fraction = float(os.environ["REMAPPED_FRACTION"])

network_file_wt = os.path.join(network_path_wt, "network-synapses.hdf5")
network_file_pd = os.path.join(network_path_pd, "network-synapses.hdf5")
input_wt = os.path.join(network_path_wt, "input-spikes.hdf5")
input_pd = os.path.join(network_path_pd, "input-spikes.hdf5")

swap = SwapToDegeneratedMorphologies(original_network_file=network_file_wt,
                                     new_network_file=network_file_pd,
                                     original_snudda_data_dir=snudda_data_wt,
                                     new_snudda_data_dir=snudda_data_pd,
                                     original_input_file=input_wt,
                                     new_input_file=input_pd)
swap.write_new_network_file()

# write_new_input_file needs a WT input-spikes.hdf5 (run `snudda input` on WT first).
if os.path.isfile(input_wt):
    swap.write_new_input_file(remap_removed_input=remap_removed_input,
                              remapped_fraction=remapped_fraction)
else:
    print(f"WARNING: {input_wt} not found - skipping input remap. "
          f"Run `snudda input` on the WT network first if you need remapped input.")

swap.close()
PY
#python3 ../../snudda/analyse/analyse_striatum.py $network_path_pd/network-synapses.hdf5

ipcluster stop

