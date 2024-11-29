# This script reads an SBML file, and creates a snudda neuromodulation file
import os
import json
import libsbml   # If missing, use:  pip install python-libsbml
import numpy as np


class ReadSBML:

    def __init__(self, filename=None, out_file=None, concentration_scale_factor=None):

        self.filename = filename
        self.out_file = out_file
        self.reader = None
        self.data = None

        if concentration_scale_factor is None:
            self.concentration_scale_factor = 1
        else:
            self.concentration_scale_factor = concentration_scale_factor
            print(f"Using concentration scale factor: {self.concentration_scale_factor}")

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
            boundary_condition = species.getBoundaryCondition()  # Assuming a default value
            # represents = species.getId()  # Assuming species ID represents itself

            species_data[species_name] = {
                "initial_concentration": initial_concentration * self.concentration_scale_factor,
                "diffusion_constant": diffusion_constant,
                "charge": charge,
                "regions": regions,
                "concentration": initial_concentration * self.concentration_scale_factor
                # "represents": represents
            }

            if atol_scale is not None:
                species_data[species_name]["atol_scale"] = atol_scale

            if boundary_condition is not None:
                species_data[species_name]["boundary_condition"] = boundary_condition

        # Using this viewer: https://sv.insysbio.com/online/
        # We have identified that the names with "*" in them are forward and backward rates
        # we need to extract those, and also extract the reaction

        # Extract reactions information
        reactions_data = {}
        for reaction in model.getListOfReactions():
            reaction_id = reaction.getName()
            reactants = " + ".join([f"{self.get_stoichiometry(reactant)} * {id_to_name[reactant.getSpecies()]}"
                                    if reactant.getStoichiometry() != 1 else id_to_name[reactant.getSpecies()]
                                   for reactant in reaction.getListOfReactants()])
            products = " + ".join([f"{self.get_stoichiometry(product)} * {id_to_name[product.getSpecies()]}"
                                   if product.getStoichiometry() != 1 else id_to_name[product.getSpecies()]
                                   for product in reaction.getListOfProducts()])

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

    def get_stoichiometry(self, reactant):
        factor = reactant.getStoichiometry()

        factor_int = int(np.round(factor))

        if np.abs(factor_int - factor) > 1e-8:
            raise ValueError(f"Error, rounding reactant incorrectly ({reactant}): {factor} -> {factor_int}")

        return factor_int

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

        if operator == "times" and current_expression.getRightChild().getName() in compartment_list:
            print(f"Removing volume {current_expression.getRightChild().getName()}")
            current_expression = current_expression.getLeftChild()

        operator = current_expression.getOperatorName()

        if operator == "minus":
            # We have both forward and backward rate to extract
            forward_rate = self._get_rate_helper(current_expression.getLeftChild(),
                                                 global_parameters, compartment_list, model)
            backward_rate = self._get_rate_helper(current_expression.getRightChild(),
                                                  global_parameters, compartment_list, model)

        else:
            forward_rate = self._get_rate_helper(current_expression,
                                                 global_parameters, compartment_list, model)
            backward_rate = None

        print(f"{forward_rate =}, {backward_rate =}")

        return forward_rate, backward_rate

    def _count_reactants(self, expression, model):

        num_species = 0
        species_list = [x.getId() for x in model.getListOfSpecies()]

        if expression.getId() in species_list or expression.getName() in species_list:
            # We found a species, number of reactants is 1
            return 1

        if expression.getOperatorName() == "times":
            left_child = expression.getLeftChild()
            right_child = expression.getRightChild()

            num_species += self._count_reactants(left_child, model)
            num_species += self._count_reactants(right_child, model)
        elif expression.getName() == "power":
            # libsbml.formulaToString(orig_expression)
            left_child = expression.getLeftChild()
            right_child = expression.getRightChild()

            base = self._count_reactants(left_child, model)
            exponent = right_child.getValue()

            num_species = base * exponent

        elif expression.getOperatorName() is None:
            return 0
        else:
            import pdb
            pdb.set_trace()
            raise ValueError(f"Not a multiplication: {expression.getOperatorName()}")

        return num_species

    def _get_rate_helper(self, expression, global_parameters, compartment_list, model):

        orig_expression = expression

        # TODO: We need to rescale the rates with self.concentration_scale_factor for two-factors
        #       and self.concentration_scale_factor**2 for three-factors

        # TODO: How to count reactants?

        while expression.getOperatorName() == "times":
            left_child = expression.getLeftChild()

            # If the left child is a volume (ie has the name of a volume in compartment_list),
            # then take the right child instead, and end the loop
            if left_child.getName() in compartment_list:
                expression = expression.getRightChild()
                break

            print(f"Removing: {expression.getRightChild().getName()}")
            expression = left_child

        rate_name = expression.getName()

        print(f"rate_name = {rate_name}")

        try:
            rate = global_parameters[rate_name]
        except:
            import traceback
            print(traceback.format_exc())
            import pdb
            pdb.set_trace()

        if self.concentration_scale_factor != 1:
            num_species = self._count_reactants(expression=orig_expression, model=model)

            rate /= self.concentration_scale_factor ** (num_species - 1)

        return rate

    def write(self, out_file=None):

        if out_file is None:
            out_file = self.out_file

        if out_file is None:
            raise ValueError(f"No outfile specified.")

        print(f"Writing JSON to {out_file}")

        with open(out_file, "wt") as f:
            json.dump(self.data, f, indent=4)


def cli():

    import argparse

    parser = argparse.ArgumentParser(description="Convert SBML to Snudda reaction diffusion with RxD",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("sbml_file", help="Input SBML file to convert")
    parser.add_argument("json_file", help="Snudda RxD JSON output file")
    parser.add_argument("--conc_scale_factor", default=1, help="Rescale concentrations by factor", type=float)

    args = parser.parse_args()

    rs = ReadSBML(filename=args.sbml_file, out_file=args.json_file, concentration_scale_factor=args.conc_scale_factor)


if __name__ == "__main__":

    cli()

