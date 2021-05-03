import json
import os

class NetworkWideNeuromodulation:

    def __init__(self):

        self.network_wide = dict()
        self.dt = None

    def set_modulation(self, neurotransmitter, neurotransmitter_key):

        """
        neurotransmitter_key is equivalent to the parameter used in the mod files, which marks the
        level, mod and maxMod parameters eg. neurotransmitter_key, ACh, would have parameters levelACh,
        modACh and maxModACh in modulated mod files.
        """

        if neurotransmitter in self.network_wide.keys():
            raise KeyError('neurotransmitter already defined')

        self.network_wide.update({neurotransmitter: {'key': neurotransmitter_key}})

        self.network_wide[neurotransmitter].update({'ion_channels': dict(),
                                                    'receptors': dict(),
                                                    'presynaptic': dict()})

    def transient(self, neurotransmitter, method, duration, parameters):

        self.network_wide[neurotransmitter].update({'method': method,
                                                    'duration': duration,
                                                    'parameters': parameters})

    def set_timestep(self, dt):

        self.dt = dt

    def ion_channel_modulation(self, neurotransmitter, cell_type, section, ion_channels):

        if cell_type not in self.network_wide[neurotransmitter]['ion_channels'].keys():

            self.network_wide[neurotransmitter]['ion_channels'].update({cell_type : dict()})

            
        self.network_wide[neurotransmitter]['ion_channels'][cell_type].update({section : ion_channels})

    def receptor_modulation(self, neurotransmitter, cell_type, receptor, modulation):

        if cell_type not in self.network_wide[neurotransmitter]['receptors'].keys():
            self.network_wide[neurotransmitter]['receptors'].update({cell_type : dict()})
                
        self.network_wide[neurotransmitter]['receptors'][cell_type].update({receptor : modulation})

    def presynaptic_receptor_modulation(self, neurotransmitter, cell_type, receptor, modulation):

        if cell_type not in self.network_wide[neurotransmitter]['presynaptic'].keys():
            self.network_wide[neurotransmitter]['presynaptic'].update({cell_type: dict()})
                
        self.network_wide[neurotransmitter]['presynaptic'][cell_type].update({receptor: modulation})

    def save(self, dir_path, name):

        if not self.dt:
            raise ValueError(' Set time step for simulation')
        else:
            for neurotransmitter in self.network_wide.keys():
                self.network_wide[neurotransmitter].update({'dt': self.dt})

        with open(os.path.join(dir_path, name), 'w') as out_file:
            json.dump(self.network_wide, out_file)
