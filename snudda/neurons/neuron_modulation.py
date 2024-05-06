import neuron.crxd as rxd


class NeuronModulation:

    def __init__(self, neuron):

        self.neuron = neuron
        self.internal_compartments = None
        self.species = dict()
        self.rates = dict()
        self.synapse_links = []

    def add_concentration_internal(self, species_name, diffusion_constant, initial_conc):

        if self.internal_compartments is None:
            self.internal_compartments = rxd.Region(self.neuron.icell.all, nrn_region='i')

        self.species[species_name] = rxd.Species(self.neuron_internal_compartments,
                                                 d=diffusion_constant,
                                                 initial=initial_conc)

    def get_species(self, *species):
        """ Example usage:
         a, b, c = self.neurons[123].modulation.get_species('A', 'B', 'C')
         """
        return [self.species[x] for x in species]

    def add_decay(self, species_name, decay_rate):

        species = self.species[species_name]
        self.rates[f"decay_{species_name}"] = rxd.rate(species, -decay_rate*species)

    def add_rate(self, reaction_name, left_side, right_side, overwrite=False):

        if not overwrite and reaction_name in self.rates:
            raise KeyError(f"Reaction {reaction_name} is already defined in neuron {self.neuron.name}")

        self.rates[reaction_name] = rxd.rate(left_side, right_side)

    def link_synapse(self, sim, species_name, synapse):

        segment = synapse.get_segment()
        ref_concentration = self.species[species_name].nodes(segment)._ref_concentration
        syn_link = self.neuron.simulate.setpointer(ref_concentration, f"{species_name}_conc", synapse)
        self.synapse_links.append(syn_link)