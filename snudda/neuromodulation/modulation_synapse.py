import json
import pathlib

class SynapseWiseNeuromodulation:

    def __init__(self):

        self.synapse_modulation = dict()

    def set_connection_type(self, connector, nt_key):

        self.synapse_modulation.update({connector: {'key' : nt_key,\
                                                     'cells': dict()}})

    def add_cell_modulation(self, connector, cell, ion_channels, receptors, extrinsic,type_connection):

        self.synapse_modulation[connector]['cells'].update({cell: {\
                                                             'ion_channels': ion_channels,\
                                                             'receptors': receptors,\
                                                             'extrinsic': extrinsic,\
                                                             'type': type_connection}})
        
    def save(self, dir_path, name):

        with open(os.path.join(dir_path, name), 'w') as out_file:
            json.dump(self.synapse_modulation, out_file)

