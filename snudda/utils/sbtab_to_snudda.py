import os
import numpy as np
import pandas as pd
import json
import quantities as pq


# See convert_sbtab_to_json.sh -- example


class ReadSBtab:

    def __init__(self, sbtab_path, out_filename,
                 compounds_filename=None,
                 reactions_filename=None,
                 parameters_filename=None,
                 constants_filename=None):

        self.sbtab_path = sbtab_path

        if compounds_filename is not None:
            self.compounds_filename = compounds_filename
        else:
            self.compounds_filename = "Compound.tsv"

        if reactions_filename is not None:
            self.reactions_filename = reactions_filename
        else:
            self.reactions_filename = "Reaction.tsv"

        if parameters_filename is not None:
            self.parameters_filename = parameters_filename
        else:
            self.parameters_filename = "Parameter.tsv"

        if constants_filename is not None:
            self.constants_filename = constants_filename
        else:
            self.constants_filename = "Constant.tsv"

        self.out_filename = out_filename

        self.compounds_data = None
        self.reactions_data = None
        self.parameters_data = None
        self.constants_data = None

        self.default_concentration_unit = "molar"
        self.unit_dict = dict()

        self.unit_dict["nmol_unit"] = pq.UnitQuantity('nanomole', pq.nano * pq.mole, symbol='nmol')
        self.unit_dict["nM_unit"] = pq.UnitQuantity('nanomolar', pq.nano * pq.molar, symbol='nM')

        self.parameters = dict()

        self.data = dict()

        self.load()

    def load(self):
        compound_path = os.path.join(self.sbtab_path, self.compounds_filename)
        reaction_path = os.path.join(self.sbtab_path, self.reactions_filename)
        parameter_path = os.path.join(self.sbtab_path, self.parameters_filename)
        constant_path = os.path.join(self.sbtab_path, self.constants_filename)

        for dest, source_file in [("compounds_data", compound_path),
                                  ("reactions_data", reaction_path),
                                  ("parameters_data", parameter_path),
                                  ("constants_data", constant_path)]:
            print(f"Loading {source_file}")
            data = pd.read_csv(source_file, sep="\t", skiprows=1)
            setattr(self, dest, data)

    def _get_optional_value(self, row, value_name, default_value=0):

        value = row.get(value_name, default_value)

        if value is None or value == np.nan:
            value = default_value

        return value

    def parse(self):

        self.data["species"] = dict()
        self.data["rates"] = dict()
        self.data["reactions"] = dict()

        original_sbtab_units = dict()

        nmol_unit = self.unit_dict["nmol_unit"]
        nM_unit = self.unit_dict["nM_unit"]

        # conc_unit = nM_unit
        conc_unit = pq.M

        # TODO: Should UNITS be nano Molar or Molar?

        for row_idx, row in self.compounds_data.iterrows():

            species_name = row["!Name"]
            species_unit = pq.CompoundUnit(row["!Unit"])  # str ==> enhet
            rescale_factor = float(species_unit.rescale(conc_unit).base)  # To get in nano molar

            species_data = { "initial_concentration": row["!InitialValue"] * rescale_factor,
                             "diffusion_constant": self._get_optional_value(row, "DiffusionConstant", 0),
                             "charge": self._get_optional_value(row, "Charge", 0),
                             "regions": self._get_optional_value(row, "Regions", ["soma_internal", "dend_internal"]),
                             "concentration": row["!InitialValue"] * rescale_factor,
                             "boundary_condition": row["!IsConstant"] == 1,
                             "unit": conc_unit.symbol
                             }

            self.data["species"][species_name] = species_data
            original_sbtab_units[species_name] = species_unit.symbol 

        self.get_parameters()

        for row_idx, row in self.reactions_data.iterrows():
            reaction_name = row["!Name"]
            reactants, products = row["!ReactionFormula"].split("<=>")
            reactants = self._format_component_str(reactants)
            products = self._format_component_str(products)
            kinetic_law = row["!KineticLaw"].split("-")

            param_name_forward = kinetic_law[0].split("*")[0].strip()
            forward_rate = self.parameters[param_name_forward]["value"]
            forward_unit = self.parameters[param_name_forward]["unit"]

            if len(kinetic_law) == 2:
                param_name_backward = kinetic_law[1].split("*")[0].strip()
                backward_rate = self.parameters[param_name_backward]["value"]
                backward_unit = self.parameters[param_name_backward]["unit"]
            else:
                backward_rate = None
                backward_unit = None

            # TODO: Allow user to specify regions
            reaction_info = {"reactants": reactants,
                             "products": products,
                             "forward_rate": forward_rate,
                             "backward_rate": backward_rate,
                             "forward_rate_unit": forward_unit,
                             "backward_rate_unit": backward_unit,
                             "regions": ["soma_internal", "dend_internal"]}

            self.data["reactions"][reaction_name] = reaction_info

    def _format_component_str(self, component_str):

        new_str = ""
        for comp in component_str.strip().split(" "):
            new_str += f" {comp}"
            if comp.isdigit():
                new_str += " *"

        return new_str.strip()

    def get_parameters(self):

        for row_idx, row in self.parameters_data.iterrows():
            parameter_name = row["!Name"]
            parameter_value = row["!Value:linspace"]


            original_unit_str = row["!Unit"]
            target_unit_str = original_unit_str.replace("milli", "").replace("nano", "")\
                .replace("pico", "").replace("femto", "")

            nmol_unit = self.unit_dict["nmol_unit"]
            nM_unit = self.unit_dict["nM_unit"]

            original_unit = pq.CompoundUnit(original_unit_str)
            target_unit = pq.CompoundUnit(target_unit_str)

            scale = float(original_unit.rescale(target_unit).base)

            print(f"{original_unit_str = }, {target_unit_str = }, {scale = }")

            try:
                self.parameters[parameter_name] = {"value": parameter_value * scale, "unit": target_unit_str}
            except:
                import traceback
                print(traceback.format_exc())
                import pdb
                pdb.set_trace()

    def _get_rates(self, rate_str):
        split_str = rate_str.split("-")
        if len(split_str) == 1:
            forward_rate_name = split_str.split("*")[0]

    def write(self, out_file=None):

        if out_file is None:
            out_file = self.out_filename

        if out_file is None:
            raise ValueError(f"No outfile specified.")

        print(f"Writing JSON to {out_file}")

        with open(out_file, "wt") as f:
            json.dump(self.data, f, indent=4)

def cli():

    import argparse

    parser = argparse.ArgumentParser(description="Convert SBML to Snudda reaction diffusion with RxD",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("sbtab_path", help="SBtab path", type=str)

    parser.add_argument("--compounds_file", help="Input SBtab compound file to convert", default=None, type=str)
    parser.add_argument("--reactions_file", help="Input SBtab reaction file to convert", default=None, type=str)
    parser.add_argument("--parameters_file", help="Input SBtab parameters file to convert", default=None, type=str)
    parser.add_argument("--constants_file", help="Input SBtab constants file to convert", default=None, type=str)

    parser.add_argument("json_file", help="Snudda RxD JSON output file")

    args = parser.parse_args()

    rs = ReadSBtab(sbtab_path=args.sbtab_path,
                   compounds_filename=args.compounds_file,
                   reactions_filename=args.reactions_file,
                   parameters_filename=args.parameters_file,
                   constants_filename=args.constants_file,
                   out_filename=args.json_file)

    rs.parse()
    rs.write(args.json_file)


if __name__ == "__main__":
    cli()
