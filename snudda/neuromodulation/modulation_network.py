import json
import os
import numpy as np
import copy

class Neuromodulation:

    def __init__(self):

        self.network_wide = dict()
        self.name_to_key = dict()
        self.dt = None
        self.type = 'replay'

    def set_modulation(self, neurotransmitter, neurotransmitter_key):

        """
        neurotransmitter_key is equivalent to the parameter used in the mod files, which marks the
        level, mod and maxMod parameters eg. neurotransmitter_key, ACh, would have parameters levelACh,
        modACh and maxModACh in modulated mod files.
        """

        if neurotransmitter_key in self.network_wide.keys():
            raise KeyError('neurotransmitter already defined')

        else:
            self.name_to_key.update({neurotransmitter: neurotransmitter_key})

            self.network_wide.update({neurotransmitter_key: {'name': neurotransmitter}})

            self.network_wide[neurotransmitter_key].update({'ion_channels': dict(),
                                                            'receptors': dict(),
                                                            'presynaptic': dict()})

    def transient(self, neurotransmitter, method, duration, parameters):

        self.network_wide[self.name_to_key[neurotransmitter]].update({'method': method,
                                                                      'duration': duration,
                                                                      'parameters': parameters})

    def set_timestep(self, dt):

        self.dt = dt

    def ion_channel_modulation(self, neurotransmitter, cell_type, section, ion_channels):

        if cell_type not in self.network_wide[self.name_to_key[neurotransmitter]]['ion_channels'].keys():

            self.network_wide[self.name_to_key[neurotransmitter]]['ion_channels'].update({cell_type : dict()})

            
        self.network_wide[self.name_to_key[neurotransmitter]]['ion_channels'][cell_type].update({section : ion_channels})

    def receptor_modulation(self, neurotransmitter, cell_type, receptor, modulation):

        if cell_type not in self.network_wide[self.name_to_key[neurotransmitter]]['receptors'].keys():
            self.network_wide[self.name_to_key[neurotransmitter]]['receptors'].update({cell_type : dict()})
                
        self.network_wide[self.name_to_key[neurotransmitter]]['receptors'][cell_type].update({receptor : modulation})

    def presynaptic_receptor_modulation(self, neurotransmitter, cell_type, receptor, modulation):

        if cell_type not in self.network_wide[self.name_to_key[neurotransmitter]]['presynaptic'].keys():
            self.network_wide[self.name_to_key[neurotransmitter]]['presynaptic'].update({cell_type: dict()})
                
        self.network_wide[self.name_to_key[neurotransmitter]]['presynaptic'][cell_type].update({receptor: modulation})

    def plot_transient(self, neurotransmitter):

        import snudda.neuromodulation.modulation as modulation
        import matplotlib.pyplot as plt

        temp = copy.deepcopy(self.network_wide[self.name_to_key[neurotransmitter]])
        duration = np.arange(0, temp['duration'], self.dt)
        temp['parameters'].update({"time_step_array": duration})
        method = getattr(modulation, temp['method'])
        modulation_vector = method(temp['parameters'])



        plt.figure()
        plt.title(f" Transient for the modulation using {neurotransmitter} ")
        plt.plot(duration, modulation_vector)
        plt.ylabel("Modulation")
        plt.xlabel( "Time (ms)")
        plt.show()


    def save(self, dir_path, name):

        if not self.dt:
            raise ValueError(' Set time step for simulation')
        else:
            for neurotransmitter in self.network_wide.keys():
                self.network_wide[neurotransmitter].update({'dt': self.dt})


        for neurotransmitter, modulation in self.network_wide.items():

            for cell, sections in modulation["ion_channels"].items():

                for section in ["axon", "dendrite", "soma"]:
                    if section in sections:
                        pass
                    else:
                        sections.update({section : list()})

        temp = dict()
        temp.update({'type': self.type})
        temp.update({'description': self.network_wide})
        with open(os.path.join(dir_path, name), 'w') as out_file:
            json.dump(temp, out_file)

if __name__ == "__main__":
    neurotransmitter = "dopamine"
    neurotransmitter_key = "DA"
    tstart = 300
    tonic = 0.2
    gmax_increase = 0.8
    tau = 500
    nl = Neuromodulation()
    nl.set_timestep(dt=0.025)
    name = "dopamine_modulation.json"
    dir_path = ""
    nl.set_modulation(neurotransmitter=neurotransmitter, neurotransmitter_key=neurotransmitter_key)
    nl.transient(neurotransmitter=neurotransmitter, \
                 method='alpha_background', \
                 duration=3000, \
                 parameters={"tstart": tstart,
                             "tonic": tonic,
                             "gmax_increase": gmax_increase,
                             "tau": tau})

    nl.ion_channel_modulation(neurotransmitter=neurotransmitter, \
                              cell_type="dSPN", \
                              section="soma", \
                              ion_channels=["kas_ms", "kaf_ms", "can_ms"])
    nl.ion_channel_modulation(neurotransmitter=neurotransmitter, \
                              cell_type="dSPN", \
                              section="dendrite", \
                              ion_channels=["kas_ms", "kaf_ms"])


    nl.save(dir_path=dir_path, name=name)
