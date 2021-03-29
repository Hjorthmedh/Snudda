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


class SnuddaNeuromodulation(SnuddaSimulate):

    def __init__(self,
                 network_path=None,
                 network_file=None,
                 input_file=None,
                 verbose=False,
                 log_file=None,
                 disable_gap_junctions=True,
                 simulation_config=None, neuromodulators=None, neuromodulator_description=None,neuromodulation_conductance=None):

        self.input_data = h5py.File(snudda_parse_path(input_file), 'r')
        self.verbose = verbose
        self.neuromodulators = neuromodulators
        self.neuro_desc = neuromodulator_description
        self.neuromodulation = dict()
        self.current_cell = None
        self.syn_gpcrs = list()
        self.cell_modulator = dict()
        self.neuromodulation_conductance = neuromodulation_conductance

        super(SnuddaNeuromodulation, self).__init__(network_path=network_path,
                                                    network_file=network_file,
                                                    input_file=input_file,
                                                    verbose=False,
                                                    log_file=log_file,
                                                    disable_gap_junctions=disable_gap_junctions,
                                                    simulation_config=simulation_config)

    # noinspection PyAttributeOutsideInit
    def setup_neurons(self):

        self.sim = NrnSimulatorParallel(cvode_active=False)

        self.write_log("Setup neurons")

        self.load_synapse_parameters()

        # The neurons this node is responsible for is in self.neuronID
        for ID in self.neuron_id:

            name = self.network_info["neurons"][ID]["name"]

            config = self.config["Neurons"][name]

            morph = snudda_parse_path(config["morphology"])
            param = snudda_parse_path(config["parameters"])
            mech = snudda_parse_path(config["mechanisms"])

            if "modulation" in config:
                modulation = snudda_parse_path(config["modulation"])
            else:
                modulation = None

            # Obs, neurons is a dictionary
            if self.network_info["neurons"][ID]["virtualNeuron"]:

                if self.input_data is None:
                    self.write_log(f"Using {self.input_file} for virtual neurons")

                name = self.network_info["neurons"][ID]["name"]

                if "activity" in self.input_data["input"][str(ID)].keys():
                    spikes = self.input_data["input"][str(ID)]["activity"]["spikes"][:, 0]

                    # Creating NEURON VecStim and vector
                    # https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3125
                    vs = h.VecStim()
                    v = h.Vector(spikes.size)
                    v.from_python(spikes)
                    vs.play(v)

                    self.virtual_neurons[ID] = dict([])
                    self.virtual_neurons[ID]["spikes"] = (v, vs, spikes)
                    self.virtual_neurons[ID]["name"] = name

                    self.pc.set_gid2node(ID, int(self.pc.id()))

                    nc = h.NetCon(vs, None)
                    self.pc.cell(ID, nc, 1)  # The 1 means broadcast spikes to other machines

            else:
                # A real neuron (not a virtual neuron that just provides input)
                parameter_id = self.network_info["neurons"][ID]["parameterID"]
                modulation_id = self.network_info["neurons"][ID]["modulationID"]

                self.neurons[ID] = NeuronModel(param_file=param,
                                               morph_file=morph,
                                               mech_file=mech,
                                               cell_name=name,
                                               modulation_file=modulation,
                                               parameter_id=parameter_id,
                                               modulation_id=modulation_id)

                # Register ID as belonging to this worker node
                self.pc.set_gid2node(ID, int(self.pc.id()))

                self.write_log(f"Node {int(self.pc.id())} - cell {ID} {name}")

                # We need to instantiate the cell
                self.neurons[ID].instantiate(sim=self.sim)
                self.set_resting_voltage(ID)

                # !!! DIRTY FIX for
                # https://github.com/BlueBrain/BluePyOpt/blob/master/bluepyopt/ephys/morphologies.py
                # This is likely the offending line, that pushes a segment to the stack
                # --> sim.neuron.h.execute('create axon[2]', icell)

                self.write_log("!!! Popping extra segment from neuron -- temp fix!")
                h.execute("pop_section()")

                # !!! END OF DIRTY FIX

                # !!! Connect a netcon and register it, taken from ballandstick's
                #     connect2target function
                nc = h.NetCon(self.neurons[ID].icell.axon[0](0.5)._ref_v, None, sec=self.neurons[ID].icell.axon[0])
                nc.threshold = 10

                self.pc.cell(ID, nc, 1)  # The 1 means broadcast spikes to other machines
                # self.pc.outputcell(ID) # -- not needed, cell was called with a 1
                # self.netConList.append(nc) -- Not needed according to Lytton et al 2016

                # Record all spikes
                self.pc.spike_record(ID, self.t_spikes, self.id_spikes)

    def add_external_input(self, input_file=None):

        if input_file is None:
            input_file = self.input_file

        self.write_log(f"Adding external (cortical, thalamic) input from {input_file}")

        self.input_data = h5py.File(input_file, 'r')

        for neuron_id, neuron in self.neurons.items():

            self.external_stim[neuron_id] = []
            name = neuron.name

            if str(neuron_id) not in self.input_data["input"]:
                self.write_log(f"Warning - No input specified for {name}", is_error=True)
                continue

            for inputType in self.input_data["input"][str(neuron_id)]:

                neuron_input = self.input_data["input"][str(neuron_id)][inputType]

                1 * np.ones((neuron_input["sectionID"].shape[0],))

                sections = self.neurons[neuron_id].map_id_to_compartment(neuron_input["sectionID"])

                # Setting individual parameters for synapses
                mod_file = SnuddaLoad.to_str(neuron_input["modFile"][()])
                param_list = json.loads(neuron_input["parameterList"][()])

                # TODO: Sanity check mod_file string
                eval_str = f"self.sim.neuron.h.{mod_file}"
                channel_module = eval(eval_str)

                for inputID, (section, section_x, paramID, nSpikes) \
                        in enumerate(zip(sections,
                                         neuron_input["sectionX"],
                                         neuron_input["parameterID"],
                                         neuron_input["nSpikes"])):
                    # We need to find cellID (int) from neuronID (string, eg. MSD1_3)

                    spikes = neuron_input["spikes"][inputID, :nSpikes] * 1e3  # Neuron uses ms
                    assert (spikes >= 0).all(), \
                        "Negative spike times for neuron " + str(neuron_id) + " " + inputType

                    # Creating NEURON VecStim and vector
                    # https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3125
                    # import pdb
                    # pdb.set_trace()
                    # noinspection PyBroadException
                    try:
                        vs = h.VecStim()
                        v = h.Vector(spikes.size)
                        v.from_python(spikes)
                        vs.play(v)
                    except:
                        import traceback
                        tstr = traceback.format_exc()
                        print(tstr)

                        assert False, "!!! If you see this, make sure that vecevent.mod is included in nrnivmodl compilation"

                    # NEURON: You can not locate a point process at position 0 or 1
                    # if it needs an ion
                    if section_x == 0.0:
                        section_x = 0.01
                    elif section_x == 1.0:
                        section_x = 0.99

                    # !!! Parameters for the tmGlut should be possible to set in the
                    # input specification !!!
                    # syn = self.sim.neuron.h.tmGlut(section(sectionX))

                    syn = None

                    for synapse_gpcr in self.neuromodulators:

                        for x in section(section_x).point_processes():

                            if synapse_gpcr in str(x).split('[')[0]:

                                channel_module_p = eval(eval_str + "_p")

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

                    print(syn)

                    nc = h.NetCon(vs, syn)

                    nc.delay = 0.0
                    # Should weight be between 0 and 1, or in microsiemens?
                    nc.weight[0] = neuron_input["conductance"][()] * 1e6  # !! what is unit? microsiemens?
                    nc.threshold = 0.1

                    # Get the modifications of synapse parameters, specific to
                    # this synapse
                    if param_list is not None and len(param_list) > 0:
                        syn_params = param_list[paramID % len(param_list)]["synapse"]

                        for par in syn_params:
                            if par == "expdata":
                                # Not a parameter
                                continue

                            if par == "cond":
                                # Ignoring cond value specified for synapse, using the
                                # one specified in the input information instead
                                continue

                            setattr(syn, par, syn_params[par])

                    self.external_stim[neuron_id].append((v, vs, nc, syn, spikes))

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

    # noinspection PyAssignmentToLoopOrWithParameter
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

                    if mechanism_name in ion_channels and mechanism_name + "_p" not in sec.psection()['density_mechs'].keys():

                        leve_list = list()

                        for type_level in sec.psection()['density_mechs'][mechanism_name].keys():

                            if 'level' in type_level:
                                leve_list.append(type_level)

                        key_list = list()

                        sec.insert(mechanism_name + "_p")

                        for syn in self.neuromodulators:

                            pointer = added_synapses[sec][str(seg.x)][syn]._ref_concentration

                            key_list.append(self.neuro_desc[syn]["key"])

                            for seg in sec:

                                for type_r in leve_list:
                                    tr = type_r.replace('level', '')

                                    setattr(seg, 'mod' + tr + "_" + mechanism_name + "_p", 1)

                                    self.sim.neuron.h.setpointer(pointer, 'level' + tr,
                                                                 getattr(seg, mechanism_name + "_p"))

                        for param, val in values.items():

                            for i, seg in enumerate(sec):

                                if 'level' not in param and param not in key_list:
                                    setattr(seg, param + "_" + mechanism_name + "_p", val[i])

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
                                              'howto': list(),
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
                self.cell_modulator[cell]['howto'].append(self.neuro_desc[syn_name]["cells"][postcell_name])
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

        if syn_name in self.neuromodulators:
            self.add_gpcr_synapse(channel_module, par_data, cell_id_source, dend_compartment, section_dist, conductance,
                                  parameter_id, synapse_type_id, axon_dist)

    def connect_network_synapses(self):

        self.write_log("connectNetworkSynapses")

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

        # This loops through all the synapses, and connects the relevant ones
        next_row = 0
        # nextRowSet = [ fromRow, toRow ) -- ie range(fromRow,toRow)
        next_row_set = self.find_next_synapse_group(next_row)

        while next_row_set is not None:
            # Add the synapses to the neuron
            self.connect_neuron_synapses(start_row=next_row_set[0], end_row=next_row_set[1])

            # Find the next group of synapses
            next_row = next_row_set[1]  # 2nd number was not included in range
            next_row_set = self.find_next_synapse_group(next_row)

    # noinspection PyBroadException
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

    # noinspection PyBroadException
    def connect_neuron_synapses(self, start_row, end_row):

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
                self.add_synapse(cell_id_source=src_id,
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

    # noinspection PyBroadException
    def add_synapse(self, cell_id_source, dend_compartment, section_dist, conductance,
                    parameter_id, synapse_type_id, axon_dist=None):

        # You can not locate a point process at
        # position 0 or 1 if it needs an ion
        if section_dist == 0.0:
            section_dist = 0.01
        if section_dist == 1.0:
            section_dist = 0.99

        (channel_module, par_data) = self.synapse_parameters[synapse_type_id]

        syn_name = str(channel_module).split('()')[0]

        syn = None

        if syn_name not in self.neuromodulators:

            for synapse_gpcr in self.neuromodulators:

                for x in dend_compartment(section_dist).point_processes():

                    if synapse_gpcr in str(x).split('[')[0]:

                        channel_module_p = eval('self.sim.neuron.h.' + str(channel_module).split('()')[0] + '_p')

                        syn = channel_module_p(dend_compartment(section_dist))

                        level = list()

                        for type_param in dir(syn):

                            if 'level' in type_param:
                                level.append(type_param)

                        pointer = x._ref_concentration

                        for tr in level:
                            setattr(syn, 'mod' + tr.replace('level', ''), 1)

                            self.sim.neuron.h.setpointer(pointer, tr, syn)

            if syn is None:
                syn = channel_module(dend_compartment(section_dist))

            if par_data is not None:
                # Picking one of the parameter sets stored in parData
                par_id = parameter_id % len(par_data)

                par_set = par_data[par_id]
                for par in par_set:
                    if par == "expdata" or par == "cond":
                        # expdata is not a parameter, and cond we take from synapse matrix
                        continue

                    try:
                        # Can be value, or a tuple/list, if so second value is scale factor
                        # for SI -> natural units conversion
                        val = par_set[par]

                        # Do we need to convert from SI to natural units?
                        if type(val) == tuple or type(val) == list:
                            val = val[0] * val[1]

                        setattr(syn, par, val)

                    except:
                        import traceback
                        tstr = traceback.format_exc()
                        print(tstr)
                        import pdb
                        pdb.set_trace()

            if axon_dist is not None:
                # axon dist is in micrometer, want delay in ms
                synapse_delay = (1e3 * 1e-6 * axon_dist) / self.axon_speed + self.synapse_delay
            else:
                synapse_delay = self.synapse_delay

            #    self.write_log(f"Synapse delay: {synapse_delay} ms")

            # What do we do if the GID does not exist?
            # print("GID exists:" + str(self.pc.gid_exists(cellIDsource)))

            if self.is_virtual_neuron[cell_id_source]:
                # Source is a virtual neuron, need to read and connect input

                nc = self.pc.gid_connect(cell_id_source, syn)
                nc.weight[0] = conductance
                nc.delay = synapse_delay
                nc.threshold = self.spike_threshold

                # Prevent garbage collection in python
                self.net_con_list.append(nc)
                self.synapse_list.append(syn)

            else:

                nc = self.pc.gid_connect(cell_id_source, syn)
                nc.weight[0] = conductance
                nc.delay = synapse_delay
                nc.threshold = self.spike_threshold

                self.net_con_list.append(nc)
                self.synapse_list.append(syn)

            return syn
