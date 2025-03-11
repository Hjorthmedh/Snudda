from mpi4py import MPI
import neuron.rxd as rxd
import numpy as np
import json


class ExtracellularNeuromodulation:

    def __init__(self, sim, volume_id=None, padding=None, dx=None):
        self.sim = sim
        self.volume_id = volume_id
        self.volume_fraction = 0.2
        self.tortuosity = 1.6

        if dx is None:
            self.dx = np.array([30e-6, 30e-6, 30e-6])
        else:
            self.dx = np.array(dx)

        self.padding = padding

        self.config_data = None

        # We set RxD option for extracellular
        rxd.options.enable.extracellular = True

        self.compartments = dict()   # "ecs" --> RxD Extracellular object

        self.species = dict()
        self.reactions = dict()

    def get_min_max_coords(self, padding=None, volume_id=None):

        """ Returns x_min, y_min, z_min, x_max, y_max, z_max """

        # TODO: Only neurons belonging to the Volume should be included in the min, max calculation
        #       use volume_id info

        x_min, y_min, z_min = self.sim.network_info["simulation_origo"]

        n_x, n_y, n_z = self.sim.network_info["num_hyper_voxels"]
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
                    compartment=("ecs", ),
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
                self.species[species_name][comp] = rxd.Species([self.compartments[comp]],
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
                                                                  backward_rate) #,
                                                                  #regions=self.compartments[region_name])

    def load_json(self, config_path=None, config_data=None):

        if config_path is not None and config_data is None:
            print(f"Loading extracellular configuration from {config_path =}")

            with open(config_path, "r") as f:
                self.config_data = json.load(f)

        elif config_data is not None and config_path is None:
            print(f"Using extracelluar configuration {config_data =}")
            self.config_data = config_data

        else:
            raise ValueError(f"Set exactly one of these variables: {config_path =}, {config_data =}")

        for species_name, species_data in self.config_data.get("species", {}).items():
            initial_concentration = species_data.get("initial_concentration", 0) * 1e3  # Convert to millimolar for RxD
            diffusion_constant = species_data.get("diffusion_constant", 0)
            charge = species_data.get("charge", 0)
            regions = species_data.get("regions", ("ecs",))
            boundary_condition = species_data.get("boundary_condition", False)

            # TODO: Read atol_scale, boundary_boundary_conditions, represents parameters

            self.add_species(species_name=species_name,
                             diffusion_constant=diffusion_constant,
                             initial_conc=initial_concentration,
                             compartment=regions, charge=charge,
                             boundary_condition=boundary_condition)

        # Black magic, set up the species variables
        species_name_vars = ",".join(self.species.keys()) + ","
        species_name_str = "','".join(self.species.keys())

        for rate_name, rate_data in self.config_data.get("rates", {}).items():
            if rate_name not in self.species:
                raise ValueError(f"Species {rate_name} is not defined. Available: {self.species.keys()}")

            rates = rate_data["rates"]
            if isinstance(rates, str) and len(rate_data["regions"]) > 1:
                rates = [rates for x in rate_data["regions"]]

            for region, rate in zip(rate_data["regions"], rates):
                exec(f"{species_name_vars} = self.get_species('{species_name_str}', region_name=region)")

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
                exec(f"{species_name_vars} = self.get_species('{species_name_str}', region_name=region)")

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


    # TODO: Add functionality to extracellular compartment, similar to what is in neuro_modulation.py
    #       1. Reactions within ECS
    #       2. Reactions between ECS and ICS
    #       How to adapt configuration file, JSON (separate from intracellular, how to deal with ECS <-> ICS"


    # TODO: 2025-02-17
    # Check scaling of ECS, based on how many voxel we have. The way to do that is to use the same number of neurons
    # but just increase the volume they are placed inside. That way only the ECS computation should increase.
    # Alternativly we just reduce the dx, thus increasing the number of voxels.
    #
    # Find out if RxD uses any JIT, are we missing some way to speed up our code?!
