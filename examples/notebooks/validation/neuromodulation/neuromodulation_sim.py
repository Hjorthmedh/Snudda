import os
import h5py
from snudda.neuromodulation.modulation_network import Neuromodulation
from snudda.utils.load import SnuddaLoad
import numpy as np

def analyse(network_path):
    
    
    validation = dict(dSPN=dict(mean=6,std=2.8),iSPN=dict(mean=-6,std=3))
    
    sl = SnuddaLoad(os.path.join(network_path, "network-synapses.hdf5"))
    neuron_types = [n["type"] for n in sl.data["neurons"]]
    f = h5py.File(os.path.join(network_path, "simulation", "test.hdf5"))
    fc = h5py.File(os.path.join(network_path, "simulation", "test_control.hdf5"))
    
    tmp = dict()
    
    for i, n in enumerate(f["neurons"].keys()):
        control_spikes = fc["neurons"][n]["spikes"]["data"][()]
        experiment_spikes = f["neurons"][n]["spikes"]["data"][()]
        neurontype = neuron_types[i]
        if experiment_spikes.size > 0 or control_spikes.size>0:
            diff = experiment_spikes.size - control_spikes.size
            z_score = (diff - validation[neurontype]["mean"]) / validation[neurontype]["std"]
            tmp.update({i : z_score})
            
    return tmp

def validation(zscores):
    
    assert len(zscores) > 5
    assert all(abs(value) < 2 for value in zscores.values())
            
def simulate(network_path):
    
    import sys
    import os
    import json
    import timeit
    import numpy as np
    import re
    from snudda.neuromodulation.neuromodulation import SnuddaSimulateNeuromodulation
    import argparse
    from collections import OrderedDict
    
    with open(os.path.join(network_path, "dopamine_modulation.json"), "r") as neuromod_f:
        neuromodulation_dict = json.load(neuromod_f, object_pairs_hook=OrderedDict)
        
    with open(os.path.join(network_path, "current_injection.json"), "r") as f:
        cell_ids_current_injection = json.load(f)
    
    network_file = os.path.join(network_path, "network-synapses.hdf5")
    output_file = os.path.join(network_path, "simulation", "test.hdf5")
    log_file = os.path.join(network_path, "log", "network-simulation-log.txt")
    simulation_config = os.path.join(network_path, "simulation_config.json")
    tSim = 2000
              
    sim = SnuddaSimulateNeuromodulation(network_file=network_file,
                                                    input_file=None,
                                                    output_file=output_file,
                                                    log_file=log_file,
                                                    verbose=True,
                                                    simulation_config=simulation_config)

    sim.setup()
    sim.apply_neuromodulation(neuromodulation_dict)
    sim.neuromodulation_network_wide()
    sim.check_memory_status()
    sim.add_volt_recording_soma()
    
    for neuron_id, data in cell_ids_current_injection.items():

        print(f"Adding to {neuron_id} the amplitude {data}")
        print(f"Within the function to {neuron_id} the amplitude {data*1e9}")
        sim.add_current_injection(neuron_id=int(neuron_id), start_time=0.5, end_time=1.5, amplitude=data)
    
    sim.run(tSim)
    sim.write_output()

def simulate_control(network_path):
    
    import sys
    import os
    import json
    import timeit
    import numpy as np
    import re
    from snudda.neuromodulation.neuromodulation import SnuddaSimulateNeuromodulation
    import argparse
    from collections import OrderedDict
    
    with open(os.path.join(network_path, "dopamine_modulation_control.json"), "r") as neuromod_f:
        neuromodulation_dict = json.load(neuromod_f, object_pairs_hook=OrderedDict)
        
    with open(os.path.join(network_path, "current_injection.json"), "r") as f:
        cell_ids_current_injection = json.load(f)
    
    network_file = os.path.join(network_path, "network-synapses.hdf5")
    output_file = os.path.join(network_path, "simulation", "test_control.hdf5")
    log_file = os.path.join(network_path, "log", "network-simulation-log.txt")
    simulation_config = os.path.join(network_path, "simulation_config.json")
    tSim = 2000
              
    sim = SnuddaSimulateNeuromodulation(network_file=network_file,
                                                    input_file=None,
                                                    output_file=output_file,
                                                    log_file=log_file,
                                                    verbose=True,
                                                    simulation_config=simulation_config)

    sim.check_memory_status()
    sim.distribute_neurons()
    sim.pc.barrier()
    sim.setup_neurons()
    sim.apply_neuromodulation(neuromodulation_dict)
    sim.neuromodulation_network_wide()
    sim.check_memory_status()
    sim.add_volt_recording_soma()
    
    for neuron_id, data in cell_ids_current_injection.items():

        print(f"Adding to {neuron_id} the amplitude {data}")
        print(f"Within the function to {neuron_id} the amplitude {data*1e9}")
        sim.add_current_injection(neuron_id=int(neuron_id), start_time=0.5, end_time=1.5, amplitude=data)
    
    sim.run(tSim)
    sim.write_output()
    
def generate_current_injection(network_path):
    
    sl = SnuddaLoad(os.path.join(network_path, "network-synapses.hdf5"))

    tmp = dict()

    for n in sl.data["neurons"]:
        p = os.path.join(n["neuronPath"], "if_info.json")
        import json
        with open(p, "r") as f:
            pdata = json.load(f)
        p = n['parameterKey']
        m = n['morphologyKey']
        nid = n["neuronID"]
        
        idx = np.argmin(np.array(pdata[p][m]['frequency']) - 10)
        current = pdata[p][m]['current'][idx]
        tmp.update({str(nid): current})
    with open(os.path.join(network_path, "current_injection.json"), "w") as f:
        json.dump(tmp,f)
        
    sim_conf = {"sampleDt":0.5e-3}
    with open(os.path.join(network_path, "simulation_config.json"), "w") as f:
        json.dump(sim_conf,f)
    
    
    
    
def create_modulation():
    
    nl = Neuromodulation()
    nl.set_timestep(dt=0.025)
    nl.set_modulation(neurotransmitter="dopamine", neurotransmitter_key="DA")
    nl.transient(neurotransmitter="dopamine", method="bath_application", duration=3000,
                 parameters={"gmax": 0})

    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="dSPN",
                              section="soma",
                              ion_channels=["kas_ms", "kaf_ms", "can_ms"])
    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="dSPN",
                              section="dendrite",
                              ion_channels=["kas_ms", "kaf_ms"])

    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="iSPN",
                              section="soma",
                              ion_channels=["kir_ms", "kas_ms", "kaf_ms", "naf_ms", "cal12_ms", "cal13_ms",
                                            "can_ms", "car_ms"])
    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="iSPN",
                              section="dendrite",
                              ion_channels=["kir_ms", "kas_ms", "kaf_ms", "naf_ms", "cal12_ms", "cal13_ms",
                                            "can_ms", "car_ms"])
    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="iSPN",
                              section="axon",
                              ion_channels=["kir_ms", "kas_ms", "kaf_ms", "naf_ms", "cal12_ms", "cal13_ms",
                                            "can_ms", "car_ms"])

    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="FSN",
                              section="soma",
                              ion_channels=["kir_fs", "kas_fs", "kaf_fs", "naf_fs"])
    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="FSN",
                              section="dendrite",
                              ion_channels=["kir_fs"])

    nl.save(dir_path=os.path.join(os.path.dirname(__file__), "networks", "dopamine_validation"), name="dopamine_modulation_control.json")
    
    nl = Neuromodulation()
    nl.set_timestep(dt=0.025)
    nl.set_modulation(neurotransmitter="dopamine", neurotransmitter_key="DA")
    nl.transient(neurotransmitter="dopamine", method="bath_application", duration=3000,
                 parameters={"gmax": 1})

    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="dSPN",
                              section="soma",
                              ion_channels=["kas_ms", "kaf_ms", "can_ms"])
    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="dSPN",
                              section="dendrite",
                              ion_channels=["kas_ms", "kaf_ms"])

    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="iSPN",
                              section="soma",
                              ion_channels=["kir_ms", "kas_ms", "kaf_ms", "naf_ms", "cal12_ms", "cal13_ms",
                                            "can_ms", "car_ms"])
    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="iSPN",
                              section="dendrite",
                              ion_channels=["kir_ms", "kas_ms", "kaf_ms", "naf_ms", "cal12_ms", "cal13_ms",
                                            "can_ms", "car_ms"])
    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="iSPN",
                              section="axon",
                              ion_channels=["kir_ms", "kas_ms", "kaf_ms", "naf_ms", "cal12_ms", "cal13_ms",
                                            "can_ms", "car_ms"])

    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="FSN",
                              section="soma",
                              ion_channels=["kir_fs", "kas_fs", "kaf_fs", "naf_fs"])
    nl.ion_channel_modulation(neurotransmitter="dopamine",
                              cell_type="FSN",
                              section="dendrite",
                              ion_channels=["kir_fs"])

    nl.save(dir_path=os.path.join(os.path.dirname(__file__), "networks", "dopamine_validation"), name="dopamine_modulation.json")