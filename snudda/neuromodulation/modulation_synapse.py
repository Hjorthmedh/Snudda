import json
import os

class NeuromodulationSynapse:

    def __init__(self):

        self.synapse_modulation = dict()
        self.type = 'adaptive'

    def set_weight(self, weight):

        self.weight = weight

    def set_connection_type(self, neuromodulation_key, connector):

        """
                neurotransmitter_key is equivalent to the parameter used in the mod files, which marks the
                level, mod and maxMod parameters eg. neurotransmitter_key, ACh, would have parameters levelACh,
                modACh and maxModACh in modulated mod files.
        """

        if connector in self.synapse_modulation.keys():
            raise KeyError('connector is already defined')

        self.synapse_modulation.update({neuromodulation_key: {'connector': connector,
                                                    'cells': dict()}})

    def add_cell_modulation(self, neuromodulation_key, cell, ion_channels, receptors, extrinsic, type_connection):

        self.synapse_modulation[neuromodulation_key]['cells'].update({cell: {
                                                             'ion_channels': ion_channels,
                                                             'receptors': receptors,
                                                             'extrinsic': extrinsic,
                                                             'type': type_connection}})
        
    def save(self, dir_path, name):

        temp = dict()
        temp.update({'type': self.type})
        temp.update({'description': self.synapse_modulation})
        temp.update({'weight': self.weight})
        with open(os.path.join(dir_path, name), 'w') as out_file:
            json.dump(temp, out_file, indent=4)

