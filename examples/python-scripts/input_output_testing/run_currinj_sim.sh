curinj_network="curinj_network04"
exp_conf_file="../../../snudda/data/experiment_config/pair_recording_output_"$curinj_network"_experiment-config-7-bobek.json"
exp_conf_file_norm="../../../snudda/data/experiment_config/pair_recording_output_"$curinj_network"_pd2_norm_experiment-config-7-bobek.json"
echo $exp_conf_file
echo "networks/"$curinj_network"_pd0"
echo "Starting pd0 sim"
mpiexec -n 3 python ../../../snudda/simulate/pair_recording.py "networks/"$curinj_network"_pd0" $exp_conf_file
echo "---------------------done with pd0"
mpiexec -n 3 python ../../../snudda/simulate/pair_recording.py "networks/"$curinj_network"_pd2" $exp_conf_file
echo "---------------------done with pd2"
mpiexec -n 3 python ../../../snudda/simulate/pair_recording.py "networks/"$curinj_network"_pd2" $exp_conf_file_norm
echo "---------------------done with pd2_norm"
read -p "Press any key to resume ..."
#mpiexec -n 3 python ../../../snudda/simulate/pair_recording.py networks/curinj_network03_pd0 ../../../snudda/data/experiment_config/pair_recording_output_curinj_network03_experiment-config-7-bobek.json