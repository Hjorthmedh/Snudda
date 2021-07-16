import numpy as np
import json


class MeasureNeuromodulation:

    def __init__(self, snudda_simulate_obj):

        self.snudda_simulate_obj = snudda_simulate_obj
        self.recording_synapse_gpcr = None

    def recording_gpcr(self):

        self.recording_synapse_gpcr = list()

        self.number_synapses = dict()
        
        for i, syn in enumerate(self.snudda_simulate_obj.syn_gpcrs):
            v = self.snudda_simulate_obj.sim.neuron.h.Vector()

            v.record(syn._ref_concentration)
            self.recording_synapse_gpcr.append(v)
            self.number_synapses.update({ i : {'cell_seg' : [str(syn.get_segment())]}})

    def write_gpcr_synapses(self, filename):

        np.savetxt(filename, self.recording_synapse_gpcr)

        with open(filename.split('.')[0] + 'pos.json', 'w') as df:
            json.dump( self.number_synapses, df)
        

    def check_mod(self, filename):

        cells = dict((k, self.snudda_simulate_obj.neurons[k]) for k in self.snudda_simulate_obj.neuron_id if not self.snudda_simulate_obj.is_virtual_neuron[k])

        data = dict()

        for index, cell in cells.items():

            for tpart in ['axon', 'basal', 'soma']:

                for comp in getattr(cell.icell, tpart):
                    for seg in comp:
                        for mech in seg:
                            if 'ptr' in mech.name():

                                for d in dir(mech):

                                    for k in ['level', 'mod', 'maxMod']:

                                        if k in d:
                                            data.update(
                                                {'_'.join([str(comp), str(seg), str(mech), str(d)]): getattr(mech, d)})

        with open(filename, 'w') as df:
            json.dump(data, df)
