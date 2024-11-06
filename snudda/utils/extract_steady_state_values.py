import json

from snudda.utils import SnuddaLoadSimulation

# TODO: Make sure unit conversion is correct.


class ExtractSteadyStateValues:

    def __init__(self, reaction_diffusion_cfg: str, network_simulation,
                 neuron_id=0, output_cfg=None):

        self.simulation_data = SnuddaLoadSimulation(network_simulation_output_file=network_simulation)

        self.reaction_diffusion_cfg = reaction_diffusion_cfg
        self.neuron_id = neuron_id

        if output_cfg is None or output_cfg == "":
            output_cfg = f"{reaction_diffusion_cfg}-updated"

        self.output_cfg = output_cfg

        with open(self.reaction_diffusion_cfg, "rt") as f:
            self.reaction_diffusion_data = json.load(f)

        self.updated_simulation_config_data = self.reaction_diffusion_data.copy()

    def extract(self, neuron_id=None):

        if neuron_id is None:
            neuron_id = self.neuron_id

        if "species" not in self.updated_simulation_config_data:
            raise ValueError(f"No 'species' defined in JSON file")

        list_species = self.simulation_data.list_data_types(neuron_id=neuron_id)

        for species_name, species_data in self.updated_simulation_config_data["species"].items():

            print(f"Parsing species {species_name}")

            for item_name in ["initial_concentration", "concentration"]:
                if item_name not in species_data:
                    raise KeyError(f"Missing {item_name} for species {species_name}")

            if species_data["concentration"] != species_data["initial_concentration"]:
                raise ValueError(f"{species_name} has concentration = {species_data['concentration']}, "
                                 f"initial_concentration = {species_data['initial_concentration']}")

            if species_name not in list_species:
                raise KeyError(f"Species {species_name} missing in RxD data for neuron {neuron_id}")

            sim_data = self.simulation_data.get_data(data_type=species_name, neuron_id=neuron_id)

            # Check that we read from soma
            if sim_data[1][0][0][0] != -1:
                raise ValueError(f"First data column for {species_name} is not soma!")

            steady_state_conc = sim_data[0][0][-1, 0]

            species_data["concentration"] = steady_state_conc
            species_data["initial_concentration"] = steady_state_conc

        print("Parsing done.")

    def write(self):

        print(f"Writing updated config file to {self.output_cfg}")

        with open(self.output_cfg, "wt") as f:
            json.dump(self.updated_simulation_config_data, f, indent=4)


def cli():

    from argparse import ArgumentParser
    parser = ArgumentParser(description="Extract steady state from simulation, and update concentrations "
                                        "in RxD configuration file.")
    parser.add_argument("reaction_diffusion_cfg", type=str)
    parser.add_argument("network_simulation_file", type=str)
    parser.add_argument("--output_cfg", type=str, default="")
    parser.add_argument("--neuron_id", type=int, default=0)
    args = parser.parse_args()

    print("Args parsed")
    essv = ExtractSteadyStateValues(reaction_diffusion_cfg=args.reaction_diffusion_cfg,
                                    network_simulation=args.network_simulation_file,
                                    neuron_id=args.neuron_id,
                                    output_cfg=args.output_cfg)

    essv.extract()
    essv.write()


if __name__ == "__main__":

    cli()