import json
import pathlib

class networkWide:

    def __init__(self):

        self.network_wide = dict()

    def set_modulation(self, neurotransmitter, nt_key):

        self.network_wide.update({neurotransmitter : { 'key' : nt_key}})

        self.network_wide[neurotransmitter].update({'ion_channels' : dict(),\
                                                    'receptors' : dict(),\
                                                    'presynaptic' : dict()})

    def transient(self,neurotransmitter,method,duration,parameters):

        self.network_wide[neurotransmitter].update({'method' : method,\
                                                    'duration' : duration,\
                                                    'parameters' : parameters})

    def ion_channel_modulation(self,neurotransmitter,cell_type,section,ion_channels):

        if cell_type not in self.network_wide[neurotransmitter]['ion_channels'].keys():

            self.network_wide[neurotransmitter]['ion_channels'].update({cell_type : dict()})

            
        self.network_wide[neurotransmitter]['ion_channels'][cell_type].update({section : ion_channels})

    def receptor_modulation(self,neurotransmitter,cell_type,receptor,modulation):

        if cell_type not in self.network_wide[neurotransmitter]['receptors'].keys():
            self.network_wide[neurotransmitter]['receptors'].update({cell_type : dict()})
                
        self.network_wide[neurotransmitter]['receptors'][cell_type].update({receptor : modulation})

    def presynaptic_receptor_modulation(self,neurotransmitter,cell_type,receptor,modulation):

        if cell_type not in self.network_wide[neurotransmitter]['presynaptic'].keys():
            self.network_wide[neurotransmitter]['presynaptic'].update({cell_type : dict()})
                
        self.network_wide[neurotransmitter]['presynaptic'][cell_type].update({receptor : modulation})

    def save(self,dir_path,name):

        out_file = open(pathlib.Path(dir_path) / name,'w')

        json.dump(self.network_wide,out_file)
