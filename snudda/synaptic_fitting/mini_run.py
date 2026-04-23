import os
import json
from snudda.simulate.nrn_simulator_parallel import NrnSimulatorParallel

from snudda.utils import snudda_parse_path
from run_synapse_run import RunSynapseRun

synapse_type = "glut"
#synapse_type = "glut2"

snudda_data = "/home/hjorth/HBP/BasalGangliaData/data/"
neuron_path = os.path.join(snudda_data, "neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20220620")
parameter_key = "p06ae8d19"
morphology_key = "mbb8e5b24"
holding_voltage = -80e-3
holding_current = 2.8045121780451154e-10


synapse_parameter_file = os.path.join(snudda_data, "synapses/striatum/tmGlut_double_config/M1-ipsi_dSPN.json")

with open(synapse_parameter_file, 'r') as f:
    print(f"Reading synapse parameters from {synapse_parameter_file}")
    synapse_parameters = json.load(f)["data"]

if synapse_type == "glut":
    synapse_parameters = {}
    

model_parameters = [ 0.0012717270422199227,
                     0.585003725942023,
                     1.570181110144922,
                     0.4307181070038531,
                     8.050679692977703e-10*0.01]


m_params = { "U": model_parameters[0],
             "tauR": model_parameters[1],
             "tauF": model_parameters[2],
             "tauRatio": model_parameters[3],
             "cond": model_parameters[4] }

c_prop = { "neuron_path": "$DATA/neurons/striatum/dspn/str-dspn-e150602_c1_D1-mWT-0728MSN01-v20220620",
           "neuron_parameter_key": "p06ae8d19",
           "neuron_morphology_key": "mbb8e5b24",
           "synapse_density": "0.05/(1+exp(-(d-150e-6)/5e-6))",
           "num_synapses": 20,
           "baseline_voltage": -0.08,
           "holding_current": 2.8045121780451154e-10 }

sim = NrnSimulatorParallel(cvode_active=False)

rsr = \
            RunSynapseRun(neuron_path=snudda_parse_path(c_prop["neuron_path"],
                                                        snudda_data),
                          neuron_morphology_key=morphology_key,
                          neuron_parameter_key=parameter_key,
                          stim_times=[0.5, 1.0, 1.5],
                          num_synapses=10,
                          synapse_density="1.0",
                          holding_voltage=c_prop["baseline_voltage"],
                          holding_current=c_prop["holding_current"],
                          synapse_type=synapse_type,
                          params=synapse_parameters,
                          time=2,
                          log_file=None,
                          synapse_section_id=None,
                          synapse_section_x=None,
                          sim=sim,
                          random_seed=1234,
                          verbose=True,
                          pc=None)

t, v, i  = rsr.run2(pars=m_params)


import matplotlib.pyplot as plt

# import pdb
# pdb.set_trace()

plt.figure()
plt.ion()
plt.plot(t, v)
plt.show()


plt.figure()
plt.plot(t, i.T)
plt.show()

import pdb
pdb.set_trace()

input()
