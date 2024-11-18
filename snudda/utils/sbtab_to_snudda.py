import os
import numpy as np
import pandas as pd


class ReadSBtab:

    def __init__(self, sbtab_path, out_filename,
                 compound_filename=None,
                 reaction_filename=None,
                 parameters_filename=None,
                 constants_filename=None):

        self.sbtab_path = sbtab_path

        if compound_filename is not None:
            self.compound_filename = compound_filename
        else:
            self.compound_filename = "Compound.tsv"

        if reaction_filename is not None:
            self.reaction_filename = reaction_filename
        else:
            self.reaction_filename = "Reaction.tsv"

        if parameters_filename is not None:
            self.parameters_filename = parameters_filename
        else:
            self.parameters_filename = "Parameter.tsv"

        if constants_filename is not None:
            self.constants_filename = constants_filename
        else:
            self.constants_filename = "Constant.tsv"

        self.out_filename = out_filename

        self.compound_data = None
        self.reaction_data = None
        self.parameters_data = None
        self.constants_data = None

        self.load()

    def load(self):
        compound_path = os.path.join(self.sbtab_path, self.compound_filename)
        reaction_path = os.path.join(self.sbtab_path, self.reaction_filename)
        parameter_path = os.path.join(self.sbtab_path, self.parameters_filename)
        constant_path = os.path.join(self.sbtab_path, self.constants_filename)

        for dest, source_file in [(self.compound_data, compound_path),
                                  (self.reaction_data, reaction_path),
                                  (self.parameters_data, parameter_path),
                                  (self.constants_data, constant_path)]:
            print(f"Loading {source_file}")
            dest = pd.read_csv(dest, header=None, sep="\t", skiprows=1)

        import pdb
        pdb.set_trace()

    def parse(self):
        pass


def cli():

    import argparse

    parser = argparse.ArgumentParser(description="Convert SBML to Snudda reaction diffusion with RxD",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("sbtab_path", help="SBtab path", type=str)

    parser.add_argument("--compound_file", help="Input SBtab compound file to convert", default=None, type=str)
    parser.add_argument("--reaction_file", help="Input SBtab reaction file to convert", default=None, type=str)
    parser.add_argument("--parameters_file", help="Input SBtab parameters file to convert", default=None, type=str)
    parser.add_argument("--constants_file", help="Input SBtab constants file to convert", default=None, type=str)

    parser.add_argument("json_file", help="Snudda RxD JSON output file")

    args = parser.parse_args()

    rs = ReadSBtab(sbtab_path=args.sbtab_path,
                   compound_filename=args.compound_file,
                   reaction_filename=args.reaction_file,
                   parameters_filename=args.parameters_file,
                   constants_filename=args.constants_file,
                   out_filename=args.json_file)


if __name__ == "__main__":
    cli()
