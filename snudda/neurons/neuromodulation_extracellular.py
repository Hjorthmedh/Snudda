from mpi4py import MPI
import neuron.rxd as rxd
import numpy as np


class ExtracellularNeuromodulation:

    def __init__(self, sim, padding=None):
        self.sim = sim
        self.volume_fraction = 0.2
        self.tortuosity = 1.6
        self.dx = np.array([30e-6, 30e-6, 30e-6])
        self.padding = padding

        # We set RxD option for extracellular
        rxd.options.enable.extracellular = True

        self.compartments = dict()

        self.species = dict()

    def get_min_max_coords(self, padding=None):

        """ Returns x_min, y_min, z_min, x_max, y_max, z_max """

        x_min, y_min, z_min = self.sim.network_info["simulation_origo"]
        n_x, n_y, n_z = self.sim.network_info["hyper_voxel_size"]
        hv_width = self.sim.network_info["hyper_voxel_width"]

        x_max = x_min + n_x * hv_width
        y_max = y_min + n_y * hv_width
        z_max = z_min + n_z * hv_width

        if padding:
            x_min -= padding
            y_min -= padding
            z_min -= padding

            x_max += padding
            y_max += padding
            z_max += padding

        return x_min, y_min, z_min, x_max, y_max, z_max

    def add_species(self, species_name, diffusion_constant, initial_conc,
                    compartment= ("ecs"),
                    charge=0, boundary_condition=False, ecs_boundary_condition=None):

        if species_name not in self.species:
            self.species[species_name] = dict()

        for comp in compartment:

            if comp not in self.compartments:
                x_min, y_min, z_min, x_max, y_max, z_max = self.get_min_max_coords(padding=self.padding)

                self.compartments[comp] = rxd.Extracellular(xlo=x_min * 1e6, ylo=y_min * 1e6, zlo=z_min * 1e6,
                                                           xhi=x_max * 1e6, yhi=y_max * 1e6, zhi=z_max * 1e6,
                                                           dx=self.dx*1e6,
                                                           volume_fraction=self.volume_fraction,
                                                           tortuosity=self.tortuosity)

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
                                                               name=species_name,
                                                               ecs_boundary_conditions=ecs_boundary_condition)

    def get_species(self, *species, region_name):
        """ Example usage:
         a, b, c = self.neurons[123].modulation.get_species('A', 'B', 'C')
         """
        return [self.species[x][region_name] for x in species]

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

    def add_reaction(self, reaction_name, left_side, right_side, forward_rate, backward_rate, region_name, overwrite=False):

        # print(f"add_reaction {reaction_name = }, {left_side = }, {right_side = }, {forward_rate = }, {backward_rate = }, {region_name = }")

        if reaction_name not in self.reactions:
            self.reactions[reaction_name] = dict()

        if not overwrite and region_name in self.reactions[reaction_name]:
            raise KeyError(f"Reaction {reaction_name} is already defined in neuron {self.neuron.name}")

        if backward_rate is None:
            backward_rate = 0

        self.reactions[reaction_name][region_name] = rxd.Reaction(left_side,
                                                                  right_side,
                                                                  forward_rate,
                                                                  backward_rate,
                                                                  regions=self.compartments[region_name])


    # TODO: Add functionality to extracellular compartment, similar to what is in neuro_modulation.py
    #       1. Reactions within ECS
    #       2. Reactions between ECS and ICS
    #       How to adapt configuration file, JSON (separate from intracellular, how to deal with ECS <-> ICS"