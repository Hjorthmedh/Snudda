# Install snudda on HPC (Dardel example)

First, here are install instructions adapted for these files. You may use your own method ofcourse. 
Just be aware that you may have to adapt paths in the batch files (and setup_env.sh).

Use Python 3.11 or above. For example on the Dardel cluster write:
```
module load cray-python
module swap PrgEnv-cray PrgEnv-gnu
module load cray-python
module load cray-mpich-abi

``` 
Change directory in terminal to somewhere you want the clone and the run:
```
git clone -b dev https://github.com/Hjorthmedh/Snudda.git
```
You probably want to install it in a virtual environment. If so run the following:
```
python -m venv snudda_env
source snudda_env/bin/activate
```
Before we install Snudda, just cd in to the Snudda clone and install some requirements manually:
```
cd Snudda
MPICC=cc pip install mpi4py
```
You may want to make sure you have latest pip
```
pip install --upgrade pip
```
The follwoing packages needs to be installed due to input_tuning.py
```
pip install quantities
pip install neo
pip install elephant
```
Then just install requirements and Snudda
```
pip install -r requirements.txt
pip install -e .[dev]
```

Now just go in to your script folder, for example /examples/parallel/KTH_PDC/Dardel_inputTuning/, 
and download your neuron data, for example BasalGangliaData and bgmod. In the example files I've used a reorganized version of bgmod to have only the essential folders,
so you will need to adapt that.
Finally, create a link to and compile your neuron mechanisms with the following commands.
```
ln -s BasalGangliaData/data/neurons/mechanisms
rm -r x86_64
nrnivmodl mechanisms
```
Now the installation is done.

# Input tuning introduction
This text describes an example process which can be used to tune your neuron models
using the Snudda, BasalGangliaData (private) and bgmod (private) repositories.
If you think you should have access to these repositories but don't, 
then contact Johannes Hjorth (Snudda, BasalGanligaData) and Alexander Kozlov (bgmod).
The following example uses parkinsonian neurons, so it will say PD0,PD1,PD2,PD3 here and there.
These are just different progression states of parkinson and can be ignored/adapted to your purpose.

#################################################################################################

Single cell, current injection-based parameter optimization/modification (bgmod).
For now, we assume this part has been done already, resulting in a Hall_of_fame.json (among other files).

#################################################################################################

Before we start. Many of the jobs will be short and not really require Dardels parallel abilities.
Sometimes the queue takes long time. A way around this is to allocate compute time an get access to a node directly.
This tends to give you access to a node faster than standing in the queue with job.
So you can write the following to get 1 hour compute time through bash:
```
salloc --nodes=1 -t 1:00:00 -A snic2021-5-492 -p main
```
Lets say you wanted to run Dardel_plot_input_tuning.job. So first you might have to run: 
```
chmod +x Dardel_plot_input_tuning.job
```
Then just run: 
```
./Dardel_plot_input_tuning.job
```
I've had some problems with the jobs that use parallel computing when running like this in bash.
So I would only do it for the jobs that use 1 node and 1 task per node (srun -n 1 ....).

#################################################################################################

A final thing before we start is that, when you are sure everything works after running all the batch scripts individualy 
and just want to repeat, then you can run the following batch script that includes the entire submission of scripts in one go.
```
Dardel_inputTuningSetupSimulatePlot.job
```

# Input tuning

Open and read through setup_env.sh, perhaps copy, and adapt it to your purpose.

Transfer Hall_of_fame.json to Parameters, Meta.json etc conversion, Meta base creation (only hash-codes).
```
sbatch Dardel_createMetaAndParameters.job
```
Tuning of background input to neurons (number of input synapses VS neuron membrane pot).
Creates input.hdf5.
First time, run it with cortical input (inputType=cortical in setup_env.sh).
It runs input_tuning.py setup.
Make sure no_meta_input=--no_meta_input in setup_env.sh.
```
sbatch Dardel_inputTuningSetup.job
```
Now simulate with the generated input.
Runs input_tuning.py simulate.
Check that nodes/tasks match number of neurons (numInputSteps) used
(takes around 6 min/250neurons if you run on shared partition with 1 node on Dardel).
```
sbatch Dardel_inputTuningSimulate.job
```
Plot e.g. ninputs vs spiking frequency
(takes around 2 min/250neurons if you run on shared partition with 1 core on Dardel).
```
sbatch Dardel_inputTuningAnalyse.job
```
Plots traces and input spikes.
```
sbatch Dardel_plot_input_tuning.job
```
Check which membr pot are in accepted range and inserts the passed voltages in to meta.json.
Runs BasalGangliaData/tools/analyse_data_and_create_meta.py
(takes around 4 min/250neurons if you run on shared partition with 1 core on Dardel).
```
sbatch Dardel_transfer2Meta.job
```
Now run input tuning setup with thalamic input (change inputType to thalamic in setup_env.sh).
Make sure no_meta_input=--no_meta_input in setup_env.sh
```
sbatch Dardel_inputTuningSetup.job
sbatch Dardel_inputTuningSimulate.job 
sbatch Dardel_inputTuningAnalyse.job
sbatch Dardel_plot_input_tuning.job
sbatch Dardel_transfer2Meta.job
```
Now make the input tuning setup again , and it will show only the input from meta.json (no step changes).
Set no_meta_input= in setup_env.sh
```
sbatch Dardel_inputTuningSetup.job
sbatch Dardel_inputTuningSimulate.job
```
Now also set inputType=corticalthalamic in setup_env.sh, in order to name figures differently.
```
sbatch Dardel_plot_input_tuning.job
sbatch Dardel_inputTuningAnalyse.job
```
#################################################################################################

Tuning/experimentation with signal input to neurons.
Now if you want you can do "normal" snudda simulations with removed network synapses.
In Dardel_runSnudda.sh there is a line that runs removeConnections.py that does this.
This part is a work in progress and there are plotting functions available if you contact me (Bo Bekkouche).

#################################################################################################
