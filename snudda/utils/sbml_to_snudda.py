# This script reads an SBML file, and creates a snudda neuromodulation file
import os
import json
import libsbml


class ReadSBML:

    def __init__(self, filename=None, out_file=None):

        self.filename = filename
        self.out_file = out_file
        self.reader = None
        self.data = None

        if filename is not None:

            if not os.path.isfile(filename):
                raise ValueError(f"File not found {filename}")
            self.parse()

            if out_file is not None:
                self.write(out_file=out_file)

    def parse(self):

        print(f"Reading {self.filename}")

        self.reader = libsbml.SBMLReader()
        document = self.reader.readSBML(self.filename)

        model = document.getModel()

        species_data = {}
        id_to_name = {}

        # Extract global parameters
        global_parameters = {}
        for param in model.getListOfParameters():
            global_parameters[param.getId()] = param.getValue()

        for species in model.getListOfSpecies():
            species_name = species.getName().replace('*', '_')
            species_id = species.getId()

            if species_id in id_to_name:
                raise KeyError(f"{species_id =} already defined.")
            id_to_name[species_id] = species_name

            initial_concentration = species.getInitialConcentration()
            charge = species.getCharge()
            initial_concentration = species.getInitialConcentration()
            diffusion_constant = 0  # Assuming a default value as it's not in SBML
            regions = ["soma_internal", "dend_internal"]  # Assuming these regions
            atol_scale = None  # Assuming a default value
            ecs_boundary_conditions = None  # Assuming a default value
            # represents = species.getId()  # Assuming species ID represents itself

            species_data[species_name] = {
                "initial_concentration": initial_concentration,
                "diffusion_constant": diffusion_constant,
                "charge": charge,
                "regions": regions,
                "concentration": initial_concentration
                # "represents": represents
            }

            if atol_scale is not None:
                species_data[species_name]["atol_scale"] = atol_scale

            if ecs_boundary_conditions is not None:
                species_data[species_name]["ecs_boundary_conditions"] = ecs_boundary_conditions

        # Using this viewer: https://sv.insysbio.com/online/
        # We have identified that the names with "*" in them are forward and backward rates
        # we need to extract those, and also extract the reaction

        # Extract reactions information
        reactions_data = {}
        for reaction in model.getListOfReactions():
            reaction_id = reaction.getName()
            reactants = " + ".join([f"{reactant.getStoichiometry()} * {id_to_name[reactant.getSpecies()]}"
                                    if reactant.getStoichiometry() != 1 else id_to_name[reactant.getSpecies()]
                                   for reactant in reaction.getListOfReactants()])
            products = " + ".join([f"{product.getStoichiometry()} *{id_to_name[product.getSpecies()]}"
                                   if product.getStoichiometry() != 1 else id_to_name[product.getSpecies()]
                                   for product in reaction.getListOfProducts()])

            import pdb
            pdb.set_trace()

            forward_rate, backward_rate = self.extract_rates(reaction, model, global_parameters)

            reactions_data[reaction_id] = {
                "reactants": reactants,
                "products": products,
                "forward_rate": forward_rate,
                "backward_rate": backward_rate,
                "regions": regions
            }

        # Combine data into the desired format
        self.data = {
            "species": species_data,
            "reactions": reactions_data
        }

    def get_stoichiometry(self, species_iterator):
        reaction_str = " + ".join([])

    def extract_rates(self, reaction, model, global_parameters):

        # TODO: Check that we can handle cAMP ** 2

        compartment_list = [x.getId() for x in model.getListOfCompartments()]
        species_list = [x.getId() for x in model.getListOfSpecies()]

        print(f"{reaction.getKineticLaw().getFormula()}")

        # This function takes a reaction and attempts to extract the forward and backward rates
        current_expression = reaction.getKineticLaw().getMath()
        operator = current_expression.getOperatorName()

        if operator == "times" and current_expression.getLeftChild().getName() in compartment_list:
            print(f"Removing volume {current_expression.getLeftChild().getName()}")
            current_expression = current_expression.getRightChild()

        operator = current_expression.getOperatorName()

        if operator == "minus":
            # We have both forward and backward rate to extract
            forward_rate = self._get_rate_helper(current_expression.getLeftChild(), global_parameters, compartment_list)
            backward_rate = self._get_rate_helper(current_expression.getRightChild(), global_parameters, compartment_list)

        else:
            forward_rate = self._get_rate_helper(current_expression, global_parameters, compartment_list)
            backward_rate = None

        print(f"{forward_rate =}, {backward_rate =}")

        return forward_rate, backward_rate

    def _get_rate_helper(self, expression, global_parameters, compartment_list):

        while expression.getOperatorName() == "times":
            left_child = expression.getLeftChild()

            # If the left child is a volume, then take the right child instead, and end the loop
            if left_child.getName() in compartment_list:
                expression = expression.getRightChild()
                break

            print(f"Removing: {expression.getRightChild().getName()}")
            expression = left_child

        rate_name = expression.getName()
        rate = global_parameters[rate_name]

        return rate

    def write(self, out_file=None):

        if out_file is None:
            out_file = self.out_file

        if out_file is None:
            raise ValueError(f"No outfile specified.")

        print(f"Writing JSON to {out_file}")

        with open(out_file, "wt") as f:
            json.dump(self.data, f, indent=4)


if __name__ == "__main__":

    filename = os.path.join(os.path.dirname(__file__), "..", "..",
                            "examples", "notebooks", "neuromodulation",
                            "data", "MODEL_speedy_reduced2.xml")

    out_file = "converted-neuromodulation-model.json"

    rs = ReadSBML(filename=filename, out_file=out_file)

