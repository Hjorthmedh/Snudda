# This file uses RXD to implement neuro_modulation

import neuron.crxd as rxd


# Användarscenario:
#
# -- Extracellulär diffusion av neurotransmittorer
# -- Koncentration precis utanför/innanför cellmembranet
# -- Aktivering av interna kaskader
# -- Uppdatering av tidskonstanter mm för jonkanaler och synapser
#
#  -- Extracellulär dopamin (med diffusion?), som påverkar synapser och jonkanaler
#  -- Kalciumdynamik, intern concentration, buffering, "decay"
#  --
#
#
#


# TODO: Håll koll på koncentrationen av X
#       Decay av koncentrationen på X
#       Öka koncentrationen av X i samband med synaptisk aktivering (eller motsvarande)
#       Läsa av koncentrationen i en MOD fil, för att modifiera dynamiken på tex jonkanalen
#

class NeuroModulation:

    def __init__(self, simulation):

        self.simulation = simulation
        self.neuron_internal_compartments = dict()
        self.species = dict()

        # rxd.options.enable.extracellular = True

        pass

    def add_neuromodulation(self, neuron_id, neuronmodulation_info):

        pass

    def add_concentration_external(self, neuron_id, species):



        pass

    def add_concentration_internal(self, neuron_id, species_name, diffusion_constant, initial_conc):

        neuron = self.simulation.neurons[neuron_id]

        if neuron_id not in self.neuron_internal_compartments:
            self.neuron_internal_compartments[neuron_id] = rxd.Region(neuron.icell.all, nrn_region='i')

        if neuron_id not in self.species:
            self.species[neuron_id] = dict()

        self.species[neuron_id][species_name] = rxd.Species(self.neuron_internal_compartments[neuron_id],
                                                            d=diffusion_constant,
                                                            initial=initial_conc)




    def add_external_volume_concentration(self, volume, species):

        pass







if __name__ == "__main__":

    from snudda.simulate import SnuddaSimulate

    ss = SnuddaSimulate()