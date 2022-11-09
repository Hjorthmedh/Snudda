#This text describes an example process which can be used to tune your neuron models
#using the Snudda, BasalGangliaData (private) and bgmod (private) repositories
#If you think you should have access to these repositories but don't, 
#then contact Johannes Hjorth (Snudda, BasalGanligaData) and Alexander Kozlov (bgmod)
#The following example uses parkinsonian neurons, so it will say PD0,PD1,PD2,PD3 here and there.
#These are just different progression states of parkinson and can be ignored/adapted to your purpose.

#################################################################################################
#Single cell, current injection-based parameter optimization/modification (bgmod)
#For now, we assume this part has been done already, resulting in a Hall_of_fame.json (among other files)

#################################################################################################
#Before we start. Many of the jobs will be short and not really require Dardels parallel abilities.
#Sometimes the queue takes long time. A way around this is to allocate compute time an get access to a node directly
#This tends to give you access to a node faster than standing in the queue with job.
#So you can write the following to get 1 hour compute time through bash:
#salloc --nodes=1 -t 1:00:00 -A snic2021-5-492 -p main
#Lets say you wanted to run Dardel_plot.job. So first you might have to run: 
chmod +x Dardel_plot.job
#Then just run: 
./Dardel_plot.job
#I've had some problems with the jobs that use parallel computing when running like this in bash.
#So I would only do it for the jobs that use 1 node and 1 task per node (srun -n 1 ....).

#################################################################################################
#A final thing before we start is that. When you are sure everything works after running all the batch scripts individualy 
#and just want to repeat, then you can run the following batch script that includes the entire submission of scripts in one go.
#Dardel_inputTuningSetupSimulatePlot.job
#
#################################################################################################
#Open and read through setup_env.sh, perhaps copy, and adapt it to your purpose

#Transfer Hall_of_fame.json to Parameters, Meta.json etc conversion, Meta base creation (only hash-codes)
sbatch Dardel_createMetaAndParameters.job

#Tuning of background input to neurons (number of input synapses VS neuron membrane pot)
#Creates input.hdf5
#First time. Run it with cortical input (inputType=cortical in setup_env.sh)
#runs input_tuning.py setup
#make sure useMeta=0 in setup_env.sh
sbatch Dardel_inputTuningSetup.job

#Now simulate with the generated input
#runs input_tuning.py simulate
#Check that nodes/tasks match number of neurons (numInputSteps) used
#(takes around 6 min/250neurons if you run on shared partition with 1 node on Dardel)
sbatch Dardel_inputTuningSimulate.job

#Plot e.g. ninputs vs spiking frequency
#(takes around 2 min/250neurons if you run on shared partition with 1 core on Dardel)
sbatch Dardel_inputTuningAnalyse.job

#Plots traces and input spikes
sbatch Dardel_plot.job

#Check which membr pot are in accepted range and inserts the passed voltages in to meta.json
#runs BasalGangliaData/tools/analyse_data_and_create_meta.py
#(takes around 4 min/250neurons if you run on shared partition with 1 core on Dardel)
sbatch Dardel_transfer2Meta.job

#Now run input tuning setup with thalamic input (change inputType to thalamic in setup_env.sh)
#make sure useMeta=0 in setup_env.sh
sbatch Dardel_inputTuningSetup.job
sbatch Dardel_inputTuningSimulate.job 
sbatch Dardel_inputTuningAnalyse.job
sbatch Dardel_plot.job
sbatch Dardel_transfer2Meta.job

#Now make the input tuning setup again , and it will show only the input from meta.json (no step changes) 
#set useMeta=1 in setup_env.sh
sbatch Dardel_inputTuningSetup.job
sbatch Dardel_inputTuningSimulate.job
#Now also set inputType=corticalthalamic in setup_env.sh, in order to name figures differently
sbatch Dardel_plot.job
sbatch Dardel_inputTuningAnalyse.job

#################################################################################################
#Tuning/experimentation with signal input to neurons
#Now if you want you can do "normal" snudda simulations with removed network synapses
#In Dardel_runSnudda.sh there is a line that runs removeConnections.py that does this.
#This part is a work in progress and there are plotting functions available if you contact me (Bo Bekkouche)
#################################################################################################