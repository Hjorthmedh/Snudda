from snudda.simulate.simulate import SnuddaSimulate
import json
import numpy as np
import snudda.neuromodulation.modulation as modulation
import snudda.neuromodulation.translator as translator


class SnuddaNeuromodulation(SnuddaSimulate):

    def __init__(self,
                 network_path=None,
                 network_file=None,
                 input_file=None,
                 verbose=False,
                 log_file=None,
                 disable_gap_junctions=True,
                 simulation_config=None):

        self.verbose = verbose
        self.neuromodulation = dict()

        super(SnuddaNeuromodulation, self).__init__(network_path=network_path,
                                                    network_file=network_file,
                                                    input_file=input_file,
                                                    verbose=False,
                                                    log_file=log_file,
                                                    disable_gap_junctions=disable_gap_junctions,
                                                    simulation_config=simulation_config)

    def neuron_vector(self, vector):

        return self.sim.neuron.h.Vector(vector)

    def apply_neuromodulation(self, neuromodulation_file):

        from pathlib import Path

        define_neuro_modulation = json.load(open(Path(neuromodulation_file), 'r'))

        for type_modulation, description_neuromodulation in define_neuro_modulation.items():
            duration = np.arange(0, description_neuromodulation['duration'], 0.025)
            method = getattr(modulation, description_neuromodulation['method'])

            description_neuromodulation['parameters'].update({"time_step_array": duration})

            modulation_vector = method(description_neuromodulation['parameters'])

            self.neuromodulation.update({
                description_neuromodulation['key']:
                    {
                        'name': type_modulation,
                        'modulation_vector': self.neuron_vector(modulation_vector),
                        'ion_channels': description_neuromodulation['ion_channels'],
                        'receptors': description_neuromodulation['receptors'],
                        'presynaptic': description_neuromodulation['presynaptic']
                    }
            })

    def neuromodulation_network_wide(self):

        for type_modulation, modulation_items in self.neuromodulation.items():

            if 'ion_channels' in modulation_items.keys():
                self.modulate_ion_channels(modulation=type_modulation, ion_channels=modulation_items['ion_channels'])

            if 'receptors' in modulation_items.keys():
                self.modulate_synapses(modulation=type_modulation, synapses=modulation_items['receptors'],
                                       intrinsic=True)
            if 'presynaptic' in modulation_items.keys():
                self.modulate_synapses(modulation=type_modulation, synapses=modulation_items['presynaptic'],
                                       extrinsic=True)

    def modulate_ion_channels(self, modulation, ion_channels):

        cells = dict((k, self.neurons[k]) for k in self.neuron_id if not self.is_virtual_neuron[k])

        modulation_key = modulation.replace('ulation', '')

        for index, cell in cells.items():

            cell_name = cell.name.split("_")[0]
            cell_modulation = ion_channels[cell_name]

            for part, modulate_section in cell_modulation.items():

                tpart = translator.translate(part)

                for comp in getattr(cell.icell, tpart):
                    for seg in comp:
                        for mech in seg:
                            if mech.name() in modulate_section:
                                setattr(mech, modulation_key, 1)
                                self.neuromodulation[modulation]['modulation_vector'].play(
                                    getattr(mech, "_ref_level" + modulation_key.replace("mod", "")),
                                    self.sim.neuron.h.dt)

    def modulate_synapses(self, modulation, synapses, intrinsic=None, extrinsic=None):

        modulation_key = modulation.replace('ulation', '')

        if extrinsic:
            for neuronID, synlist in self.external_stim.items():
                for syntuple in synlist:

                    cell_name = self.neurons[neuronID].name.split("_")[0]
                    syn = syntuple[3]
                    syn_name = str(syn).split("[")[0]

                    if cell_name in synapses.keys() and syn_name in synapses[cell_name].keys():
                        self.modulate_receptor(syn=syn, modulation=modulation,
                                               modulation_key=modulation_key)

        if intrinsic:
            for syn in self.synapse_list:

                cell_name = str(syn.get_segment()).split("_")[0]
                syn_name = str(syn).split("[")[0]

                if cell_name in synapses.keys() and syn_name in synapses[cell_name].keys():
                    self.modulate_receptor(syn=syn, modulation=modulation,
                                           modulation_key=modulation_key)

    def modulate_receptor(self, syn, modulation, modulation_key):

        setattr(syn, modulation_key, 1)

        self.neuromodulation[modulation]['modulation_vector'].play(
            getattr(syn, "_ref_level" + modulation_key.replace("mod", "")), self.sim.neuron.h.dt)

        ############################################################################
