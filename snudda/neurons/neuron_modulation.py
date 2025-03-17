# TODO: How to handle conductance for DA channel, we have that DA vesicles contain 33000 DA molecules
#       (watch out for overflow of integers in snudda synapse matrix)

# TODO: How to handle co-release?

# TODO: Move DA species to external compartment. IMPORTANT.

import numpy as np
from mpi4py import MPI  # This must be imported before neuron, to run parallel
import neuron.rxd as rxd
import json
from itertools import chain


class NeuronModulation:
    should_update_rxd_nodes = True

    def __init__(self, neuron, config_file=None):

        self.neuron = neuron
        self.compartments = dict()
        self.species = dict()
        self.rates = dict()
        self.reactions = dict()
        self.config_data = {}
        self.config_file = config_file

        self.extracellular_region = None  # Extracellular comparment

        self.node_cache = None

        self.build = {"soma_internal": lambda neuron_dummy: self.set_default_compartments("soma", nrn_region="i"),
                      "dend_internal": lambda neuron_dummy: self.set_default_compartments("dend", nrn_region="i"),
                      "axon_internal": lambda neuron_dummy: self.set_default_compartments("axon", nrn_region="i"),
                      "soma_external": lambda neuron_dummy: self.set_default_compartments("soma", nrn_region="o"),
                      "dend_external": lambda neuron_dummy: self.set_default_compartments("dend", nrn_region="o"),
                      "axon_external": lambda neuron_dummy: self.set_default_compartments("axon", nrn_region="o") }

    def __del__(self):

        # Clear old rxd objects -- this avoids segfault for unittests
        rxd.rxd.byeworld()

    def set_default_compartments(self, compartment, nrn_region="i"):
        if compartment == "axon":
            section_list = self.neuron.icell.axon
        elif compartment == "dend":
            section_list = self.neuron.icell.dend
        elif compartment == "soma":
            section_list = self.neuron.icell.soma

        return rxd.Region(section_list, nrn_region=nrn_region)

    def add_species(self, species_name, diffusion_constant, initial_conc,
                    charge=0,
                    compartment=("soma_internal", "dend_internal"),
                    boundary_condition=False):

        # boundary_condition is used when either a species concentration is fixed, or externally driven (specified)

        # print(f"Adding species: {species_name} to {compartment}")

        if species_name not in self.species:
            self.species[species_name] = dict()

        for comp in compartment:
            if comp not in self.compartments:
                self.compartments[comp] = self.build[comp](self.neuron)

            if comp in self.species[species_name]:
                raise ValueError(f"{species_name = } already defined for {comp = }")

            if boundary_condition:
                # print(f"Fixing {species_name} concentration to constant {initial_conc}")
                # The concentration is fixed
                self.species[species_name][comp] = rxd.Parameter(self.compartments[comp],
                                                                 name=species_name,
                                                                 value=initial_conc,
                                                                 charge=charge)

            else:

                # TODO: Add atol_scale etc...
                self.species[species_name][comp] = rxd.Species(self.compartments[comp],
                                                               d=diffusion_constant,
                                                               initial=initial_conc,
                                                               charge=charge,
                                                               name=species_name)

    def get_species_OLD(self, *species, region_name):
        """ Example usage:
         a, b, c = self.neurons[123].modulation.get_species('A', 'B', 'C')
         """
        return [self.species[x][region_name] for x in species]

    def get_species(self, *species, region_name):
        species_list = []

        for s in species:
            if "__" in s:
                species_name, species_region = s.split("__")
                if species_region != "ecs":
                    raise NotImplementedError("Only 'ecs' region can currently be coupled to")
                species_list.append(self.extracellular_region.species[species_name][species_region])

            else:
                species_list.append(self.species[s][region_name])

        return species_list

    def get_species_with_regions(self, *species, region_name):
        species_list = []

        for s in species:
            _species = None
            _region = None

            if "__" in s:
                species_name, species_region = s.split("__")
                if species_region != "ecs":
                    raise NotImplementedError("Only 'ecs' region can currently be coupled to")

                _species = self.extracellular_region.species[species_name][species_region]
                species_region = self.extracellular_region.compartments[species_region]

            else:
                _species = self.species[s][region_name]
                species_region = self.compartments[region_name]
                # species_region = region_name

            try:
                species_list.append(_species[species_region])
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()


        return species_list


    def add_decay(self, species_name, decay_rate):

        species = self.species[species_name]
        self.rates[f"decay_{species_name}"] = rxd.rate(species, -decay_rate*species)

    def add_rate(self, species_name, left_side, right_side, region_name, overwrite=False):

        # print(f"Add rate {species_name = }, {left_side = }, {right_side = }, {region_name = }")

        if species_name not in self.rates:
            self.rates[species_name] = dict()

        if not overwrite and region_name in self.rates[species_name]:
            raise KeyError(f"Reaction {species_name} is already defined in neuron {self.neuron.name}")

        self.rates[species_name][region_name] = rxd.Rate(left_side, right_side,
                                                         regions=self.compartments[region_name])

    def get_neuron_regions(self, region_list):

        regions = []
        for region in region_list:
            regions.append(self.compartments[region])

        return chain(*regions)

    def add_reaction(self, reaction_name, left_side, right_side, forward_rate, backward_rate, region_name, overwrite=False):

        # print(f"add_reaction {reaction_name = }, {left_side = }, {right_side = }, {forward_rate = }, {backward_rate = }, {region_name = }")

        if reaction_name not in self.reactions:
            self.reactions[reaction_name] = dict()

        if not overwrite and region_name in self.reactions[reaction_name]:
            raise KeyError(f"Reaction {reaction_name} is already defined in neuron {self.neuron.name}")

        if backward_rate is None:
            backward_rate = 0

        try:
            self.reactions[reaction_name][region_name] = rxd.Reaction(left_side,
                                                                      right_side,
                                                                      forward_rate,
                                                                      backward_rate,
                                                                      regions=self.compartments[region_name])
        except:
            import traceback
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()

    def add_multi_compartment_reaction(self, reaction_name,
                                       left_side: rxd.rxdmath._Arithmeticed,
                                       right_side: rxd.rxdmath._Arithmeticed,
                                       forward_rate, backward_rate,
                                       region_name,
                                       overwrite=False):

        if reaction_name not in self.reactions:
            self.reactions[reaction_name] = dict()

        if not overwrite and region_name in self.reactions[reaction_name]:
            raise KeyError(f"Reaction {reaction_name} is already defined in neuron {self.neuron.name}")

        if backward_rate is None:
            backward_rate = 0

        # TODO: 2025-03-04: Replace with include flux instead
        # https://www.neuron.yale.edu/phpBB/viewtopic.php?t=4631
        # https://colab.research.google.com/drive/1d-4bi427QZksinLq46LaYMFpx7Fb1t-N?usp=sharing <--
        #
        # - Identifiera vilken ECS voxel som ett compartment ligger i
        #   Vi har koordinater för NEURON compartment, och vi kan beräkna coordinater för ECS voxlarna
        #
        # - Hur hanterar vi om en section är i två ECS voxlar, tar vi den ena bara?
        #
        # - Vi behöver kolla att MOD filerna inte skriver över eller skrivs över
        #   Dvs att de kan ändra koncentrationerna genom att addera eller ta bort, inte att
        #   de helt ersätter koncentrationen som kan vara modifierad med ett eget värde.
        #
        # Exempel:
        #
        # def flux(a,b, A=surface_area, D=diffusion_coeff,  dist=dist):
        #     return A * D * (a.value - b.value)*1e-18/dist # in mol/ms
        #
        #
        # for prend, postnd in zip(NT.nodes(cytpre), NT.nodes(cytpost)):
        #     postnd.include_flux(lambda:flux(prend,postnd), units='mol/ms')
        #     prend.include_flux(lambda:flux(postnd,prend), units='mol/ms')
        #
        # ... we have problem med rxd.MultiCompartmentReaction:
        # Running Neuron simulator 1000 ms, with dt=5.0
        # Illegal instruction

        # TODO: rxd.MultiCompartmentReaction, 2025-02-21 -- UNDERSTAND THIS!!
        # We are trying to use "neuron_modulation_of_channels.ipynb", to make reactions
        # between intracellular space, and extracellular space... we then need to use this instead
        # of the normal rxd.Reaction:

        # leak = rxd.MultiCompartmentReaction(ca[er]<>ca[cyt], gleak, gleak, membrane=cyt_er_membrane)
        # cyt_er_membrane = rxd.Region(h.allsec(), geometry = rxd.ScalableBorder(1, on_cell_surface=False))

        membrane = rxd.Region([x for x in self.neuron.icell.all], name='membrane', geometry=rxd.membrane)
        try:
            # import pdb
            # pdb.set_trace()
            self.reactions[reaction_name][region_name] = rxd.MultiCompartmentReaction(left_side != right_side,
                                                                                      forward_rate, backward_rate,
                                                                                      membrane=membrane)
        except:
            import traceback
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()

    def _get_nodes(self, species, force_update=None):
        if force_update is None:
            force_update = NeuronModulation.should_update_rxd_nodes

        import neuron.crxd as rxd
        import itertools

        if force_update:
            print("Forcing rxd update...", flush=True)
            NeuronModulation.should_update_rxd_nodes = False

            rxd.species._nrn_shape_update()
            if rxd.initializer.is_initialized():
                print("Updating node data... (takes ≈ 1 microcentury)")
                rxd.rxd._update_node_data()
            else:
                rxd.initializer._do_init()
            print("RxD update completed.")
        else:
            rxd.initializer._do_init()

        species._all_intracellular_nodes = []
        if species._intracellular_nodes:
            for r in species._regions:
                if r in species._intracellular_nodes:
                    species._all_intracellular_nodes += species._intracellular_nodes[r]
        # The first part here is for the 1D -- which doesn't keep live node objects -- the second part is for 3D
        species._all_intracellular_nodes = [
            nd for nd in species._all_intracellular_nodes[:] if nd.sec
        ]
        return rxd.nodelist.NodeList(
            list(itertools.chain.from_iterable([s.nodes for s in species._secs]))
            + species._all_intracellular_nodes
            + species._extracellular_nodes
        )

    def build_node_cache(self):
        print(f"Build node cache {self.neuron.name} ({self.neuron.icell})", flush=True)
        self.node_cache = {}

        for species_name, species_data in self.species.items():
            if species_name not in self.node_cache:
                self.node_cache[species_name] = {}

            for region_name, species in species_data.items():

                if region_name not in self.node_cache[species_name]:
                    self.node_cache[species_name][region_name] = {}

                # This row below might need to be sped up
                #all_nodes = species.nodes

                all_nodes = self._get_nodes(species, force_update=None)

                for node in all_nodes:
                    if node._sec._sec not in self.node_cache[species_name][region_name]:
                         self.node_cache[species_name][region_name][node._sec._sec] = ([], [])

                    self.node_cache[species_name][region_name][node._sec._sec][0].append(node)
                    self.node_cache[species_name][region_name][node._sec._sec][1].append(node._location)

                # Now we want to also extract the sec_x
                for node_key in self.node_cache[species_name][region_name]:
                    node_list, node_x = self.node_cache[species_name][region_name][node_key]
                    self.node_cache[species_name][region_name][node_key] = (node_list, np.array(node_x))

        print(f"Node cache built.", flush=True)

    def get_node_from_cache(self, species_name, seg, region_name):
        try:
            node_list, node_x = self.node_cache[species_name][region_name][seg.sec]
            idx = np.argmin(np.abs(node_x - seg.x))
        except:
            import traceback
            print(f"get_node_from_cache failed: {species_name = }, {seg = }, {region_name = }")
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()
        return node_list[idx]

    def clear_cache(self):
        self.node_cache = None

    def link_synapse(self, species_name, region: str, synapse, flux_variable: str):
        # region: "soma_internal", "soma_external", "dend_internal", "dend_external"

        # print(f"link_synapses {species_name}, {region}")

        if self.node_cache is None:
            self.build_node_cache()

        node = self.get_node_from_cache(species_name=species_name, seg=synapse.get_segment(), region_name=region)
        node.include_flux(synapse, flux_variable)

        # self.species[species_name][region].nodes(synapse.get_segment())[0].include_flux(synapse, flux_variable)

    def load_json(self, config_path=None, extracellular_regions=None, neuron_region=None):

        # print(f"Parsing neuromodulation json: {config_path}")

        if extracellular_regions is not None and neuron_region in extracellular_regions:
            self.extracellular_region = extracellular_regions[neuron_region]
        else:
            self.extracellular_region = None

        if config_path is None:
            config_path = self.config_file

        with open(config_path, "r") as f:
            self.config_data = json.load(f)

        # print(f"Parsing species")

        for species_name, species_data in self.config_data.get("species", {}).items():
            initial_concentration = species_data.get("initial_concentration", 0) * 1e3  # Convert to millimolar for RxD
            diffusion_constant = species_data.get("diffusion_constant", 0)
            charge = species_data.get("charge", 0)
            regions = species_data.get("regions", ("soma_internal", "dendrites_internal"))
            boundary_condition = species_data.get("boundary_condition", False)

            # TODO: Read atol_scale, boundary_boundary_conditions, represents parameters

            self.add_species(species_name=species_name,
                             diffusion_constant=diffusion_constant,
                             initial_conc=initial_concentration,
                             compartment=regions, charge=charge,
                             boundary_condition=boundary_condition)

        # Black magic, set up the species variables, don't remove trailing ","
        species_name_vars = ",".join(self.species.keys()) + ","
        species_name_str = "','".join(self.species.keys())

        if self.extracellular_region is not None:
            ecs_species_name_vars = "__ecs,".join(self.extracellular_region.species.keys()) + "__ecs"
            ecs_species_name_str = "__ecs','".join(self.extracellular_region.species.keys()) + "__ecs"

            species_name_vars = ecs_species_name_vars + "," + species_name_vars
            species_name_str = ecs_species_name_str + "','" + species_name_str

        # Also add extracellular species, if they do exist

        # print(f"Parsing rates.")

        for rate_name, rate_data in self.config_data.get("rates", {}).items():
            if rate_name not in self.species:
                raise ValueError(f"Species {rate_name} is not defined. Available: {self.species.keys()}")

            rates = rate_data["rates"]
            if isinstance(rates, str) and len(rate_data["regions"]) > 1:
                rates = [rates for x in rate_data["regions"]]

            for region, rate in zip(rate_data["regions"], rates):
                # exec(f"{species_name_vars} = self.get_species('{species_name_str}', region_name=region)")
                exec(f"{species_name_vars} = self.get_species_with_regions('{species_name_str}', region_name=region)")

                try:
                    right_side = eval(rate)
                except:
                    print(f"Problem evaluating rate {rate}")
                    import traceback
                    print(traceback.format_exc())
                    import pdb
                    pdb.set_trace()

                self.add_rate(species_name=rate_name,
                              left_side=self.get_species(rate_name, region_name=region)[0],
                              right_side=right_side,
                              region_name=region)

        # print(f"Parsing reactions")

        for reaction_name, reaction_data in self.config_data.get("reactions", {}).items():

            # BE CAREFUL, BECAUSE VALUES ARE UNSCALED AND FED DIRECTLY INTO RxD
            # which uses mmolar, ms, litres as units. So usually forward rate needs
            # to account for the concentration if rescaled.
            # TODO: Add proper unit handling.
            forward_rates = reaction_data["forward_rate"]
            backward_rates = reaction_data["backward_rate"]

            if not isinstance(forward_rates, (tuple, list)):
                forward_rates = [forward_rates] * len(reaction_data["regions"])

            if not isinstance(backward_rates, (tuple, list)):
                backward_rates = [backward_rates] * len(reaction_data["regions"])

            if not (len(reaction_data["regions"]) == len(forward_rates) == len(backward_rates)):
                raise ValueError(f"{reaction_data} incompatible lengths for regions, forward and backward rates")

            for region, forward_rate, backward_rate in zip(reaction_data["regions"], forward_rates, backward_rates):
                # TODO: Sanitise species_name_str before exec call
                # exec(f"{species_name_vars} = self.get_species('{species_name_str}', region_name=region)")
                try:
                    exec(f"{species_name_vars} = self.get_species_with_regions('{species_name_str}', region_name=region)")
                except:
                    import traceback
                    print(traceback.format_exc())
                    import pdb
                    pdb.set_trace()

                left_side = eval(reaction_data["reactants"])
                right_side = eval(reaction_data["products"])

                try:
                    # Ta inte bort ettan (Do not remove the 1) -- or else...
                    n_left = sum((left_side*1)._items.values())
                    n_right = sum((right_side*1)._items.values())
                except:
                    import traceback
                    print(traceback.format_exc())
                    import pdb
                    pdb.set_trace()

                assert n_left >= 1 and n_right >= 1

                # First 1e3 is for 1/second to 1/ms, second term is to correct for concentration
                left_scale_from_SI = 1e-3 * 1e-3 ** (n_left - 1)
                right_scale_from_SI = 1e-3 * 1e-3 ** (n_right - 1)

                scaled_forward_rate = forward_rate * left_scale_from_SI if forward_rate is not None else forward_rate
                scaled_backward_rate = backward_rate * right_scale_from_SI if backward_rate is not None else backward_rate

                print(f"Reaction name: {reaction_name}, left: {left_side}, right: {right_side}")
                print(f"k_forward: {forward_rate} (scaled: {scaled_forward_rate})")
                print(f"k_backward: {backward_rate} (scaled: {scaled_backward_rate})")

                self.add_reaction(reaction_name=reaction_name,
                                  left_side=left_side,
                                  right_side=right_side,
                                  forward_rate=scaled_forward_rate,
                                  backward_rate=scaled_backward_rate,
                                  region_name=region)
                """

                self.add_multi_compartment_reaction(reaction_name=reaction_name,
                                  left_side=left_side,
                                  right_side=right_side,
                                  forward_rate=scaled_forward_rate,
                                  backward_rate=scaled_backward_rate,
                                  region_name=region)
                """

    def concentration_from_vector(self, species_name, concentration_vector, time_vector, interpolate=True):

        # Loops over all nodes in node_cache, and sets a vector to play
        print(f"Playing concentration vector for {species_name} in all neurons.")

        if self.node_cache is None:
            raise ValueError("node_cache not build (build_node_cache)")

        if species_name not in self.node_cache:
            raise ValueError(f"{species_name} not present in node_cache, does {self.neuron.name} have RxD species?")

        for region_name, node_dictionary in self.node_cache[species_name].items():
            for node_name, node_data in node_dictionary.items():
                for nd in node_data[0]:
                    concentration_vector.play(nd._ref_concentration, time_vector, interpolate)