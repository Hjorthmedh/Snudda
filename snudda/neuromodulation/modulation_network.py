import json
import os

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

    def save(self, dir_path, name):

        if not self.dt:
            raise ValueError(' Set time step for simulation')
        else:
            for neurotransmitter in self.network_wide.keys():
                self.network_wide[neurotransmitter].update({'dt': self.dt})

        temp = dict()
        temp.update({'type': self.type})
        temp.update({'description': self.network_wide})
        with open(os.path.join(dir_path, name), 'w') as out_file:
            json.dump(temp, out_file)
