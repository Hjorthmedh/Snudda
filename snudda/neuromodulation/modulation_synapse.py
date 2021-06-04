import json
import os

class NeuromodulationSynapse:

    def __init__(self):

        self.synapse_modulation = dict()

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

        with open(os.path.join(dir_path, name), 'w') as out_file:
            json.dump(self.synapse_modulation, out_file, indent=4)

