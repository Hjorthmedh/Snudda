from snudda.simulate.simulate import SnuddaSimulate
import json
import numpy as np
import snudda.neuromodulation.modulation as modulation
import snudda.neuromodulation.translator as translator
from snudda.simulate.nrn_simulator_parallel import NrnSimulatorParallel
from snudda.utils.snudda_path import snudda_parse_path
import h5py
from snudda.neurons.neuron_model_extended import NeuronModel
from snudda.utils.load import SnuddaLoad
from neuron import h


class SnuddaSimulateNeuromodulationSynapse(SnuddaSimulate):

    def __init__(self,
                 network_path=None,
                 network_file=None,
                 input_file=None,
                 verbose=False,
                 log_file=None,
                 disable_gap_junctions=True,
                 simulation_config=None, neuromodulators=None, neuromodulator_description=None,neuromodulation_conductance=None):

        self.verbose = verbose
        self.neuromodulators = neuromodulators
        self.neuro_desc = neuromodulator_description
        self.neuromodulation = dict()
        self.current_cell = None
        self.syn_gpcrs = list()
        self.cell_modulator = dict()
        self.neuromodulation_conductance = neuromodulation_conductance

        super(SnuddaSimulateNeuromodulationSynapse, self).__init__(network_path=network_path,
                                                                   network_file=network_file,
                                                                   input_file=input_file,
                                                                   verbose=False,
                                                                   log_file=log_file,
                                                                   disable_gap_junctions=disable_gap_junctions,
                                                                   simulation_config=simulation_config)
        self.custom_setup_bool = True

    def reorder_cell_info(self, cell_modulator):

        reorder = dict()

        for cell, instructions in cell_modulator.items():

            for i in range(len(instructions['dend_compartment'])):

                dendcomp = instructions['dend_compartment'][i]

                if dendcomp in reorder.keys() and instructions['section_dist'][i] in reorder[dendcomp]['section_dist']:

                    reorder[dendcomp]['section_dist'].append(instructions['section_dist'][i])
                    reorder[dendcomp]['connection'][str(instructions['section_dist'][i])].append(
                        instructions['synapse'][i])
                    reorder[dendcomp]['precell'][str(instructions['section_dist'][i])].append(
                        instructions['precell'][i])
                elif dendcomp in reorder.keys() and instructions['section_dist'][i] not in reorder[dendcomp]['section_dist']:

                    reorder[dendcomp]['section_dist'].append(instructions['section_dist'][i])
                    reorder[dendcomp]['connection'].update(
                        {str(instructions['section_dist'][i]): [instructions['synapse'][i]]})
                    reorder[dendcomp]['precell'].update(
                        {str(instructions['section_dist'][i]): [instructions['precell'][i]]})

                else:
                    reorder.update({dendcomp: {'section_dist': [instructions['section_dist'][i]],
                                               'connection':
                                                   {str(instructions['section_dist'][i]): [instructions['synapse'][i]]},
                                               'precell': {str(instructions['section_dist'][i]): [
                                                   instructions['precell'][i]]}}})

        self.chose_implementation(reorder)

    def chose_implementation(self, reorder):

        added_synapses = dict()

        for syn_in_section in self.neuromodulators:

            for sec, sec_info in reorder.items():

                for seg in sec:

                    sgr = getattr(self.sim.neuron.h, syn_in_section)(seg)

                    self.syn_gpcrs.append(sgr)

                    if str(seg.x) in sec_info['precell'].keys():

                        for cell_id_source in sec_info['precell'][str(seg.x)]:

                            nc = self.pc.gid_connect(cell_id_source, sgr)
                            print(self.neuromodulation_conductance)
                            nc.weight[0] = self.neuromodulation_conductance
                            nc.delay = 1
                            nc.threshold = -20

                            self.net_con_list.append(nc)
                    self.synapse_list.append(sgr)

                    if sec in added_synapses.keys() and str(seg.x) in added_synapses[sec].keys():

                        added_synapses[sec][str(seg.x)].update({syn_in_section: sgr})

                    elif sec in added_synapses.keys() and str(seg.x) not in added_synapses[sec].keys():

                        added_synapses[sec].update({str(seg.x): {syn_in_section: sgr}})
                    else:
                        added_synapses.update({sec: {str(seg.x): {syn_in_section: sgr}}})

        for sec, sec_info in reorder.items():

            cell_name = str(sec).split("_")[0]

            sec_name = sec.name().split('.')[-1].split('[')[0]

            tpart = translator.re_translation[sec_name]

            ion_channels = list()

            for syn in self.neuromodulators:

                if cell_name in self.neuro_desc[syn]['cells'].keys():
                    ion_channels = ion_channels + self.neuro_desc[syn]['cells'][cell_name]['ion_channels'][tpart]

            for seg in sec:

                for mechanism_name, values in sec.psection()['density_mechs'].items():

                    if mechanism_name in ion_channels and mechanism_name + "_ptr" not in sec.psection()['density_mechs'].keys():

                        leve_list = list()

                        for type_level in sec.psection()['density_mechs'][mechanism_name].keys():

                            if 'level' in type_level:
                                leve_list.append(type_level)

                        key_list = list()

                        sec.insert(mechanism_name + "_ptr")

                        for syn in self.neuromodulators:

                            pointer = added_synapses[sec][str(seg.x)][syn]._ref_concentration

                            key_list.append('mod'+self.neuro_desc[syn]["key"])

                            for seg in sec:

                                for type_r in leve_list:
                                    tr = type_r.replace('level', '')
                                    print(type_r)
                                    print(mechanism_name)
                                    print('mod' + tr + "_" + mechanism_name + "_ptr")
                                    setattr(seg, 'mod' + tr + "_" + mechanism_name + "_ptr", 1)
                                    print(getattr(seg, 'mod' + tr + "_" + mechanism_name + "_ptr"))

                                    self.sim.neuron.h.setpointer(pointer, 'level' + tr,
                                                                 getattr(seg, mechanism_name + "_ptr"))

                        for param, val in values.items():

                            for i, seg in enumerate(sec):
                                print(key_list)
                                if 'level' not in param and param not in key_list:
                                    print(param)
                                    setattr(seg, param + "_" + mechanism_name + "_ptr", val[i])

            for mech_name in ion_channels:
                sec.uninsert(mech_name)

    def add_gpcr_synapse(self, channel_module, par_data, cell_id_source, dend_compartment, section_dist, conductance,
                         parameter_id, synapse_type_id, axon_dist):

        cell = dend_compartment.cell()

        if self.current_cell is None:

            self.current_cell = cell
            self.cell_modulator = dict()
            self.cell_modulator.update({cell: dict()})
            self.cell_modulator[cell].update({'precell': list(),
                                              'postcell': list(),
                                              'dend_compartment': list(),
                                              'section_dist': list(),
                                              'method': list(),
                                              'key': list(),
                                              'index': list(),
                                              'synapse': list(),
                                              'cell_name': dend_compartment.cell()})

        elif self.current_cell != cell:

            self.current_cell = None

            self.reorder_cell_info(self.cell_modulator)

        else:

            syn_name = str(channel_module).split('()')[0]
            postcell_name = str(dend_compartment.cell()).split('_')[0]

            if postcell_name in self.neuro_desc[syn_name]["cells"].keys():
                self.cell_modulator[cell]['precell'].append(cell_id_source)
                self.cell_modulator[cell]['postcell'].append(postcell_name)
                self.cell_modulator[cell]['dend_compartment'].append(dend_compartment)
                self.cell_modulator[cell]['section_dist'].append(section_dist)
                self.cell_modulator[cell]['method'].append(self.neuro_desc[syn_name]["cells"][postcell_name])
                self.cell_modulator[cell]['key'].append(self.neuro_desc[syn_name])
                self.cell_modulator[cell]['synapse'].append(channel_module)

    def add_mark_gpcr(self, cell_id_source, dend_compartment, section_dist, conductance,
                      parameter_id, synapse_type_id, axon_dist=None):

        if section_dist == 0.0:
            section_dist = 0.01
        if section_dist == 1.0:
            section_dist = 0.99

        (channel_module, par_data) = self.synapse_parameters[synapse_type_id]

        syn_name = str(channel_module).split('()')[0]

        print(syn_name)

        if syn_name in self.neuromodulators:
            self.add_gpcr_synapse(channel_module, par_data, cell_id_source, dend_compartment, section_dist, conductance, parameter_id, synapse_type_id, axon_dist)

    def custom_setup(self):

        self.write_log("custom setup")

        # This loops through all the synapses, and connects the relevant ones
        next_row = 0
        # nextRowSet = [ fromRow, toRow ) -- ie range(fromRow,toRow)
        next_row_set = self.find_next_synapse_group(next_row)

        while next_row_set is not None:
            # Add the synapses to the neuron
            self.connect_neuron_synapses_gpcr(start_row=next_row_set[0], end_row=next_row_set[1])

            # Find the next group of synapses
            next_row = next_row_set[1]  # 2nd number was not included in range
            next_row_set = self.find_next_synapse_group(next_row)

    def connect_neuron_synapses_gpcr(self, start_row, end_row):

        source_id_list = self.synapses[start_row:end_row, 0]
        dest_id = self.synapses[start_row, 1]
        assert (self.synapses[start_row:end_row, 1] == dest_id).all()

        # Double check mapping
        assert self.pc.gid2cell(dest_id) == self.neurons[dest_id].icell, \
            "GID mismatch: " + str(self.pc.gid2cell(dest_id)) \
            + " != " + str(self.neurons[dest_id].icell)

        synapse_type_id = self.synapses[start_row:end_row, 6]
        axon_distance = self.synapses[start_row:end_row, 7]  # Obs in micrometers

        sec_id = self.synapses[start_row:end_row, 9]
        dend_sections = self.neurons[dest_id].map_id_to_compartment(sec_id)
        sec_x = self.synapses[start_row:end_row, 10] / 1000.0  # Convert to number 0-1

        # conductances are stored in pS (because we use INTs),
        # Neuron wants it in microsiemens??!
        conductance = self.synapses[start_row:end_row, 11] * 1e-6
        parameter_id = self.synapses[start_row:end_row, 12]

        voxel_coords = self.synapses[start_row:end_row, 2:5]
        self.verify_synapse_placement(dend_sections, sec_x, dest_id, voxel_coords)

        for (src_id, section, section_x, s_type_id, axon_dist, cond, p_id) \
                in zip(source_id_list, dend_sections, sec_x, synapse_type_id,
                       axon_distance, conductance, parameter_id):

            try:
                # !!!
                self.add_mark_gpcr(cell_id_source=src_id,
                                   dend_compartment=section,
                                   section_dist=section_x,
                                   synapse_type_id=s_type_id,
                                   axon_dist=axon_dist,
                                   conductance=cond,
                                   parameter_id=p_id)
            except:
                import traceback
                tstr = traceback.format_exc()
                self.write_log(tstr, is_error=True)
                import pdb
                pdb.set_trace()

    def get_synapse(self, channel_module, dend_compartment, section_dist):

        self.write_log('using new function')

        syn_name = str(channel_module).split('()')[0]

        syn = None

        if syn_name not in self.neuromodulators:

            self.write_log('inside the important function')

            for synapse_gpcr in self.neuromodulators:

                for x in dend_compartment(section_dist).point_processes():

                    if synapse_gpcr in str(x).split('[')[0]:

                        channel_module_p = eval('self.sim.neuron.h.' + str(channel_module).split('()')[0] + "_ptr")

                        syn = channel_module_p(dend_compartment(section_dist))

                        level = list()

                        for type_param in dir(syn):

                            if 'level' in type_param:
                                level.append(type_param)

                        pointer = x._ref_concentration

                        # for n in [x for x in dir(np) if "sin" in x]: print(n)

                        for tr in level: # [x for x in dir(syn) if 'level' in x]
                            setattr(syn, 'mod' + tr.replace('level', ''), 1)

                            self.sim.neuron.h.setpointer(pointer, tr, syn)

        if syn is None:
            syn = channel_module(dend_compartment(section_dist))
            self.write_log('doing it normal way')

        return syn

    def get_external_input_synapse(self, eval_str, section, section_x, channel_module):

        syn = None

        self.write_log('ADDING EXTERNAL inside function')

        for synapse_gpcr in self.neuromodulators:

            for x in section(section_x).point_processes():

                if synapse_gpcr in str(x).split('[')[0]:

                    channel_module_p = eval(eval_str + "_ptr")

                    syn = channel_module_p(section(section_x))

                    level = list()

                    for type_param in dir(syn):

                        if 'level' in type_param:
                            level.append(type_param)

                    pointer = x._ref_concentration

                    for tr in level:
                        setattr(syn, 'mod' + tr.replace('level', ''), 1)

                        self.sim.neuron.h.setpointer(pointer, tr, syn)
        if syn is None:
            syn = channel_module(section(section_x))

        return syn

