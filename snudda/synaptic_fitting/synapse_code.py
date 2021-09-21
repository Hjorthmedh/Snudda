from snudda.synaptic_fitting.synapse_class import SynapseClass
from snudda.simulate import SnuddaSimulate
import timeit
import neuron
import numpy as np
from .make_dendrogram import make_dendrogram

class SingleCellNetwork(SynapseClass):

    def __init__(self, neuron_morphology,
                 neuron_mechanisms,
                 neuron_parameters,
                 neuron_modulation,
                 verbose=True):

        super(SingleCellNetwork, self).__init__(
                 neuron_morphology,
                 neuron_mechanisms,
                 neuron_parameters,
                 neuron_modulation,
                 verbose=verbose)

        self.dend_recordings = dict()
        self.syn_v_recordings = dict()
        self.iclamps = dict()
        self.num_iclamps = 0
        self.num_vclamps = 0
        self.vclamps = dict()
        self.current_record = dict()

    def add_parameters(self,parameters):

        for p in parameters:

            if p == "cond":
                cond = self.si_to_natural_units(p, parameters[p])

            else:
                v = self.si_to_natural_units(p, parameters[p])
                for s in self.synapses:
                    setattr(s, p, v)

    def run(self, time=None, cond=1e-8):

        start = timeit.default_timer()

        if time is None:
            time = self.time
        else:
            self.time = time

        self.t_save = neuron.h.Vector()
        self.t_save.record(neuron.h._ref_t)

        neuron.h.finitialize(self.holding_voltage * 1e3)
        for ncs in self.nc_syn:
            ncs.weight[0] = cond

        self.set_resting_voltage(self.holding_voltage * 1e3)

        neuron.h.v_init = self.holding_voltage * 1e3
        neuron.h.tstop = time * 1e3
        self.write_log("About to start NEURON... stay safe")
        neuron.h.run()
        self.write_log("NEURON actually completed?!")

        print(f"{timeit.default_timer()-start} s")
        # Convert results back to SI units
        return (np.array(self.t_save) * 1e-3,
                np.array(self.v_save) * 1e-3)

    def dend_record(self, num, x):

        v_dend_save = neuron.h.Vector()
        v_dend_save.record(getattr(self.neuron.icell.dend[num](x), '_ref_v'))

        if num not in self.dend_recordings.keys():
            self.dend_recordings.update({num: dict()})
            self.dend_recordings[num].update({str(x): v_dend_save})
        else:
            self.dend_recordings[num].update({str(x) : v_dend_save})

    def dend_list(self):

        return [d for d in self.neuron.icell.dend]

    def num_sections(self):

        dendrites = self.dend_list()

        return len(dendrites)

    def synapse_record(self, num, x):

        v_syn_save = neuron.h.Vector()
        v_syn_save.record(self.neuron.icell.dend[num](x)._ref_v)

        if num not in self.dend_recordings.keys():
            self.syn_v_recordings.update({num: dict()})
            self.syn_v_recordings[num].update({str(x): v_syn_save})
        else:
            self.syn_v_recordings[num].update({str(x): v_syn_save})

    def save(self):

        pass

    def define_iclamp(self,amp,delay,dur,num=None,x=None,position=None):

        if position=='soma':
            iclamp = self.sim.neuron.h.IClamp(self.neuron.icell.soma[0](0.5))
            iclamp.amp = amp
            iclamp.delay = delay
            iclamp.dur = dur
        else:
            iclamp = self.sim.neuron.h.IClamp(self.neuron.icell.dend[num](x))
            iclamp.amp = amp
            iclamp.delay = delay
            iclamp.dur = dur
        self.num_iclamps = self.num_iclamps + 1
        self.iclamps.update({self.num_iclamps: iclamp})

    def define_seclamp(self,holding_voltage,dur,num=None,x=None,position=None):

        if position=='soma':

            v_clamp = self.sim.neuron.h.SEClamp(self.neuron.icell.soma[0](0.5))
            v_clamp.rs = 1e-9
            v_clamp.amp1 = holding_voltage
            v_clamp.dur1 = dur
        else:
            v_clamp = self.sim.neuron.h.SEClamp(self.neuron.icell.dend[num](x))
            v_clamp.rs = 1e-9
            v_clamp.amp1 = holding_voltage
            v_clamp.dur1 = dur
        self.num_vclamps = self.num_vclamps + 1
        self.vclamps.update({self.num_vclamps: v_clamp})

    def record_vclamp_current(self):

        for k, vcs in self.vclamps.items():
            vector_v = self.sim.neuron.h.Vector()
            self.current_record.update({k : vector_v.record(vcs._ref_i)})


    def remove_synapse(self,num):

        num = num - 1

        for i, s in enumerate(self.synapses):
            n = int(s.get_segment().sec.name().split('[')[-1].split(']')[0])
            print(n)
            if num == n:
                print('Remove synapse')
                self.synapses.remove(s)
                self.synapse_locations = np.delete(self.synapse_locations, i, axis=0)

    def reset(self):

        self.iclamps = dict()

    def make_dendrogram(self):

        swc = self.nm.get_morphology(self.neuron_parameter_id, 0)

        make_dendrogram(swc)

    def return_variable(self):

            pass


                
                
                
                
                
                
                
                
                
                
                
