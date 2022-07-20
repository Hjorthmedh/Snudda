import os
import shutil
import numpy as np
import json
import copy


def add_feature(model, input_type, network_path, ratio):
    meta = os.path.join(network_path, "meta.json")

    input_map = f"input_map_{model}.json"

    whole = f"inputs_population_{model}.json"

    input_definition = "input_definition.json"

    with open(whole, "r") as f:
        r = json.load(f)

    with open(input_map, "r") as f:
        d = json.load(f)

    with open(meta, "r") as f:
        metas = json.load(f)

    with open(input_definition, "r") as f:
        i_d = json.load(f)

    shutil.copy(meta, os.path.join(network_path, "meta_old.json"))

    for p, morph in metas.items():
        for m, morph_name in morph.items():
            name = "_".join([p, m])

            if "input" not in metas[p][m].keys():
                metas[p][m].update({"input": dict()})

            definition = copy.deepcopy(i_d)

            if name in d.keys():
                print("Model in passed, take average value")

                definition.update({"nInputs": int(np.mean(d[name]) * ratio), "frequency": 1.0})
                metas[p][m]["input"].update({input_type: definition})

            else:
                print("Model not in passed, take average of whole population")
                definition.update({"nInputs": int(np.mean(r['1']) * ratio), "frequency": 1.0})
                metas[p][m]["input"].update({input_type: definition})

    with open(meta, "w") as f:
        json.dump(metas, f, indent=4, sort_keys=True)
