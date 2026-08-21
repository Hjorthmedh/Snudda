import os
import json
from snudda.synaptic_fitting.parameter_bookkeeper import ParameterBookkeeper

def export_synapse_parameters(parameter_data_file_name,
                              output_filename,
                              extra_synapse_parameters=None,
                              extra_synapse_parameters_key=None,
                              parameter_key="1",
                              overwrite=True):

    if extra_synapse_parameters is not None:
        raise NotImplementedError(f"We need to read the extra synapse parameters and store them to the output file also")

    parameter_list = ["U", "tauR", "tauF", "tauRatio", "nmda_ratio"]
    parameter_bookkeper = ParameterBookkeeper(old_book_file=parameter_data_file_name, n_max=100)

    best_param = parameter_bookkeper.get_best_parameterset()

    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    if not overwrite and os.path.isfile(output_filename):
        with open(output_filename, "rt") as f:
            data = json.load(f)
    else:
        data = dict()

    syn_param = dict(zip(parameter_list, best_param))

    # Because tauR > tau we use tauRatio during optimisation, convert back to tau in parameter file
    # since MOD file does not know about tauRatio
    if "tauRatio" in syn_param:
        if "tau" in syn_param:
            raise ValueError(f"tau = tauR * tauRatio, however tau already set in {syn_param = }")

        syn_param["tau"] = syn_param["tauR"] * syn_param["tauRatio"]
        del syn_param["tauRatio"]

    data[parameter_key] = {
        "meta": { "parameter_data_file": os.path.basename(parameter_data_file_name)},
        "synapse": syn_param
    }

    with open(output_filename, "wt") as f:
        json.dump(data, f, indent=4)

    print(f"Writing {output_filename}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser("Export synapse parameters")
    parser.add_argument("synapse_parameter_file_name")
    parser.add_argument("output_filename")
    parser.add_argument("--parameter_key", default="1")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    export_synapse_parameters(parameter_data_file_name=args.synapse_parameter_file_name,
                              output_filename=args.output_filename,
                              parameter_key=args.parameter_key,
                              overwrite=args.overwrite)