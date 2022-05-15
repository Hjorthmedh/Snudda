import pathlib
import os
from snudda.input.input_tuning import InputTuning
import json
from mpi4py import MPI

def distributor(metas):

    self.comm = MPI.COMM_WORLD
    self.rank = self.comm.Get_rank()
    self.size = self.comm.Get_size()

    pc = h.ParallelContext()
    self.pc = pc
    self.gidlist = []

    for i in range(int(pc.id()), self.modulation_setup["population"], int(pc.nhost())):
        self.gidlist.append(i)

    if int(pc.id()) == 0:
        "do root stuff"

    else:
        "get instructions on which pair to run"






def running_input_tuning(parameter_id, morphology_id, num_replicas,
                         neuron_types, model_name, input_type, num_input_min,
                         num_input_max, input_duration, input_frequency_range):
    network_path = "_".join(neuron_types, parameter_id, morphology_id)
    pathlib.Path(network_path).mkdir(parents=True, exist_ok=True)

    input_tuning = InputTuning(network_path)
    ### SNUDDA_DATA set in job file...or .sh file

    input_tuning.setup_network(neurons_path=neurons_path,
                               model_name=model_name,
                               num_replicas=num_replicas,
                               neuron_types=neuron_types,
                               parameter_id=parameter_id,
                               morphology_id=morphology_id)
    input_tuning.setup_input(input_type=input_type,  # eg. "cortical" or "thalamic"
                             num_input_min=num_input_min,
                             num_input_max=num_input_max,
                             input_duration=input_duration,
                             input_frequency_range=input_frequency_range)

    return network_path


if __name__ == "__main__":
    print()

    '''
    
    For job file
    
    
    
    srun -n $N_WORKERS x86_64/special -mpi -python $SNUDDA_DIR/input/input_tuning.py simulate 
    
    '''
