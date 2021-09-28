from typing import List, Any

from snudda.simulate.simulate import SnuddaSimulate
import snudda.neuromodulation.modulation as modulation
import snudda.neuromodulation.translator as translator
from snudda.simulate.nrn_simulator_parallel import NrnSimulatorParallel
from snudda.utils.snudda_path import snudda_parse_path
from snudda.neurons.neuron_model_extended import NeuronModel
from snudda.utils.load import SnuddaLoad
from neuron import h
import json
import numpy as np
import h5py


class SnuddaSimulateNeuromodulationSynapse(SnuddaSimulate):

    def __init__(self,
                 network_path=None,
                 network_file=None,
                 input_file=None,
                 verbose=False,
                 log_file=None,
                 disable_gap_junctions=True,
                 simulation_config=None,
                 neuromodulator_description=None):
        """

        @type neuromodulation_weight: float

        """
        self.neuromodulator_description = neuromodulator_description['description']
        self.neuromodulation = dict()
        self.current_cell = None
        self.syn_gpcrs = list()
        self.cell_modulator = dict()
        self.neuromodulation_weight = neuromodulator_description['weight']
        self.connector = [info['connector'] for info in self.neuromodulator_description.values()]
        self.module_connector = [k+'()'for k in self.connector]
        self.mod_str = dict(zip(self.module_connector, self.connector))


        super(SnuddaSimulateNeuromodulationSynapse, self).__init__(network_path=network_path,
                                                                   network_file=network_file,
                                                                   input_file=input_file,
                                                                   verbose=verbose,
                                                                   log_file=log_file,
                                                                   disable_gap_junctions=disable_gap_junctions,
                                                                   simulation_config=simulation_config)

        # Change the self.custom_setup from None, and execute the custom setup code within this file
        self.cell_type_ion_channels_per_section = self.ion_channels_per_section()
        self.key_list = self.modulation_keys()

    def setup(self):

        self.check_memory_status()
        self.distribute_neurons()
        self.setup_neurons()
        self.check_memory_status()
        self.pc.barrier()

        # Neuromodulation requires this to be run, before connect_network

        #self.synapse_parameters is loaded in self.setup_neurons, hence it is None before

        self.neuromodulation_synapse_ids = [sid for sid, synapse in self.synapse_parameters.items() if
                                            str(synapse[0]) in self.module_connector]
        self.neuromodulation_setup()

        self.connect_network()
        self.check_memory_status()
        self.pc.barrier()

    def neuromodulation_setup(self):

        # This loops through all the synapses, and connects the relevant ones
        # nextRowSet = [ fromRow, toRow ) -- ie range(fromRow,toRow)
        next_row_set = self.find_next_synapse_group(next_row=0)

        while next_row_set is not None:
            # Add the synapses to the neuron
            self.connect_neuron_synapses_gpcr(start_row=next_row_set[0], end_row=next_row_set[1])

            # Find the next group of synapses
            next_row_set = self.find_next_synapse_group(next_row_set[1])  # 2nd number was not included in range

    def get_neuromodulator_synapses(self, synapse_type_ids):

        return [i for i, synapse_type_id in enumerate(synapse_type_ids) if
                synapse_type_id in self.neuromodulation_synapse_ids]

    def modulation_keys(self):

        key_list = list()
        for neuromodulation_key in self.neuromodulator_description.keys():
            key_list.append('level' + neuromodulation_key)
            key_list.append('mod' + neuromodulation_key)
        return key_list

    def ion_channels_per_section(self):

        cell_type = dict()

        for neuromodulator_key, info_for_neuromodulator in self.neuromodulator_description.items():

            for cell_type_name in info_for_neuromodulator['cells'].keys():

                cell_type.update({cell_type_name: dict()})

                for tpart in info_for_neuromodulator['cells'][cell_type_name]['ion_channels'].keys():

                    if tpart in cell_type[cell_type_name].keys():
                        cell_type[cell_type_name][tpart] = cell_type[cell_type_name][tpart] + \
                                                          info_for_neuromodulator['cells'][cell_type_name][
                                                              'ion_channels'][tpart]
                    else:
                        cell_type[cell_type_name].update({tpart: info_for_neuromodulator['cells'][
                            cell_type_name]['ion_channels'][tpart]})

        return cell_type

    def get_ion_channels_per_section(self, cell_type):
        return self.cell_type_ion_channels_per_section[cell_type]

    def connect_neuron_synapses_gpcr(self, start_row, end_row):

        source_id_list, dest_id, dend_sections, sec_id, sec_x, synapse_type_id, axon_distance, \
            conductance, parameter_id = self.get_synapse_info(start_row=start_row, end_row=end_row)

        gpcr_synapse_index = self.get_neuromodulator_synapses(synapse_type_id)

        # Filter away synapse_type_id not conc* connector and check if the start and end row defines the cell,
        # hence you can send all information to add_gpcr_synapse

        channel_modules = [str(self.synapse_parameters[synapse_type_id[i]][0]) for i in gpcr_synapse_index]

        gpcr_synapse_info = np.take([source_id_list, dend_sections, sec_x, sec_id], gpcr_synapse_index, axis=1).transpose()

        #rewrite code as it is sorted on sec_id, jump in step of sec_id and add to dict
        sort_idx = gpcr_synapse_info[:, -1].argsort()

        gpcr_synapse_info = gpcr_synapse_info[sort_idx]

        self.add_gpcr_synapse(channel_modules, gpcr_synapse_info)

    def add_gpcr_synapse(self, channel_modules, gpcr_info):

        if gpcr_info.shape[0] == 0:
            self.write_log(f'add_gpcr_synapse : Empty gpcr_info for {channel_modules}. Check if the gpcr synapses are correctly defined for all cell types')

        cell = gpcr_info[0][1].cell()
        postcell_type = str(cell).split('_')[0]
        cell_information = dict()
        current_section = gpcr_info[0][-1]
        start_index = 0
        dend_section = gpcr_info[start_index][1]

        for i, dend_info in enumerate(gpcr_info):

            if current_section != dend_info[-1]:

                cell_information.update({dend_section: {'precell': gpcr_info[start_index:i, 0]}})
                cell_information[dend_section].update({'section_dist': gpcr_info[start_index:i, 2]})
                cell_information[dend_section].update({'mod': channel_modules[start_index:i]})

                start_index = i
                dend_section = gpcr_info[start_index][1]
                current_section = gpcr_info[start_index][-1]

        cell_information.update({dend_section: {'precell': gpcr_info[start_index:i, 0]}})
        cell_information[dend_section].update({'section_dist': gpcr_info[start_index:i, 2]})
        cell_information[dend_section].update({'mod': channel_modules[start_index:i]})

        self.connect_gpcr_synapse_to_ion_channels(cell_information, postcell_type)

    def connect_gpcr_synapse_to_ion_channels(self, cell_information, cell_type):

        cell_added_synapses = self.add_gpcrs_in_cell_segments(cell_information=cell_information)

        ion_channels_per_section = self.get_ion_channels_per_section(cell_type)

        for sec, sec_info in cell_information.items():

            tpart = translator.re_translation[sec.name().split('.')[-1].split('[')[0]]

            ion_channels = ion_channels_per_section[tpart]

            for mechanism_name, values in sec.psection()['density_mechs'].items():

                mechanism_name_ptr = mechanism_name + "_ptr"

                if mechanism_name in ion_channels:

                    level_list = [type_level for type_level in [*sec.psection()['density_mechs'][mechanism_name].keys()] if 'level' in type_level]
                    mod_key_list = [f"mod{n.replace('level','')}_{mechanism_name_ptr}" for n in level_list]
                    sec.insert(mechanism_name_ptr)

                    for syn in self.connector:

                        for segment in sec:
                            seg_x_str = str(segment.x)
                            pointer = cell_added_synapses[sec][seg_x_str][syn]._ref_concentration
                            # talk to NEURON maybe they can help change this, so you don't have to replace the mechanisms, crashes with uninitialised pointers
                            for neurotransmitter_level, mod_key in zip(level_list, mod_key_list):

                                setattr(segment, mod_key, 1)
                                self.sim.neuron.h.setpointer(pointer, neurotransmitter_level,
                                                             getattr(segment, mechanism_name_ptr))

                    # Parameterize the pointer version of density_mech, skip level and mod, as that would turn off modulation
                    for param, val in values.items():
                        for i, segmentet in enumerate(sec):
                            if param not in self.key_list:
                                setattr(segmentet, '_'.join([param, mechanism_name_ptr]), val[i])

            for mech_name in ion_channels:
                sec.uninsert(mech_name)

    def add_gpcrs_in_cell_segments(self, cell_information):

        added_synapses = dict()
        # Add all the gpcrs to the segments which have been marked in this cell

        for str_mod_name, syn_in_section in self.mod_str.items():

            for sec, sec_info in cell_information.items():
                step = 1/(sec.nseg*2)
                for seg in sec:

                    synapse_gpcr = getattr(self.sim.neuron.h, syn_in_section)(seg)

                    self.syn_gpcrs.append(synapse_gpcr)

                    str_modules = [str(k) for k in sec_info['mod']]

                    if len(np.where(((seg.x - step) < sec_info['section_dist']) & ((seg.x - step) < sec_info['section_dist']))) > 0:

                        for cell_id_source in np.take(sec_info['precell'], np.where((np.array(str_modules) == str_mod_name) & ((seg.x - step) < sec_info['section_dist']) & ((seg.x - step) < sec_info['section_dist'])))[0]:
                            nc = self.pc.gid_connect(cell_id_source, synapse_gpcr)
                            nc.weight[0] = self.neuromodulation_weight
                            nc.delay = self.synapse_delay
                            nc.threshold = self.spike_threshold

                            self.net_con_list.append(nc)
                    self.synapse_list.append(synapse_gpcr)

                    if sec in added_synapses.keys() and str(seg.x) in added_synapses[sec].keys():

                        added_synapses[sec][str(seg.x)].update({syn_in_section: synapse_gpcr})

                    elif sec in added_synapses.keys() and str(seg.x) not in added_synapses[sec].keys():

                        added_synapses[sec].update({str(seg.x): {syn_in_section: synapse_gpcr}})
                    else:
                        added_synapses.update({sec: {str(seg.x): {syn_in_section: synapse_gpcr}}})

        return added_synapses

    def get_synapse(self, channel_module, dend_compartment, section_dist):

        syn = None
        # add lookup on channel_module to skip split
        if str(channel_module).split('()')[0] not in self.connector:

            channel_module_p = eval(f"self.sim.neuron.h.{str(channel_module).split('()')[0]}_ptr")

            for point_process_in_section in dend_compartment(section_dist).point_processes():

                if str(point_process_in_section).split('[')[0] in self.connector:

                    syn = channel_module_p(dend_compartment(section_dist))

                    level = [x for x in dir(syn) if 'level' in x]

                    pointer = point_process_in_section._ref_concentration

                    for neurotransmitter_key in level:
                        # remove this parameter by setting default 1
                        setattr(syn, f"mod{neurotransmitter_key.replace('level', '')}", 1)
                        self.sim.neuron.h.setpointer(pointer, neurotransmitter_key, syn)

        if syn is None:
            syn = channel_module(dend_compartment(section_dist))

        return syn

    def get_external_input_synapse(self, channel_module, section, section_x):

        syn = None

        for point_process_in_section in section(section_x).point_processes():

            if str(point_process_in_section).split('[')[0] in self.connector and str(channel_module).split('()')[0] not in self.connector:

                channel_module_p = eval('self.sim.neuron.h.' + str(channel_module).split('()')[0] + "_ptr")

                syn = channel_module_p(section(section_x))

                level = [x for x in dir(syn) if 'level' in x]

                pointer = point_process_in_section._ref_concentration

                for neurotransmitter_key in level:
                    setattr(syn, 'mod' + neurotransmitter_key.replace('level', ''), 1)
                    self.sim.neuron.h.setpointer(pointer, neurotransmitter_key, syn)

        if syn is None:
            syn = channel_module(section(section_x))

        return syn

