import numpy as np
import json

class MeasureNeuromodulation:

    def __init__(self,obj):

        self.obj = obj


    def recording_gpcr(self):

        self.recording_synapse_gpcr = list()
        for syn in self.obj.syn_gpcrs:
            v = self.obj.sim.neuron.h.Vector()

            v.record(syn._ref_concentration)
            self.recording_synapse_gpcr.append(v)

    def write_gpcr_synapses(self, filename):

        np.savetxt(filename, self.recording_synapse_gpcr)

    def check_mod(self, filename):

        cells = dict((k, self.obj.neurons[k]) for k in self.obj.neuron_id if not self.obj.is_virtual_neuron[k])

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

        

