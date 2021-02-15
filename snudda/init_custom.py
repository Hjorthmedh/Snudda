from snudda.init import SnuddaInit
from collections import OrderedDict
import json
import os

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Init custom network")
    parser.add_argument("network", type=str, help="Network path")
    parser.add_argument("--cellspec", type=str, help="Cell spec directory", default=None)
    # TODO: Add cell number parameters
    args = parser.parse_args()

    sim_name = args.network
    cell_spec = args.cellspec

    connect_neurons = False

    # simName = "networks/FSmorphTest2orig"
    # simName = "networks/FSmorphTest1b"
    # simName = "LTStest"
    # simName = "networks/twoFS"
    # simName = "networks/FSmorphTest4"
    # simName = "networks/3types-v2"
    # simName = "networks/SynTest-v6" # MSMS tuning
    # simName = "networks/SynTest-v15"

    config_name = os.path.join(sim_name, "network-config.json")
    cnc = SnuddaInit(struct_def={}, config_file=config_name, num_population_units=1)
    # cnc.defineStriatum(nMSD1=500,nMSD2=500,nFS=0,nLTS=0,nChIN=30,volumeType="cube")
    # cnc.defineStriatum(nMSD1=120,nMSD2=120,nFS=20,nLTS=0,nChIN=0,volumeType="slice")
    # cnc.defineStriatum(nMSD1=0,nMSD2=0,nFS=10000,nLTS=0,nChIN=0,volumeType="slice")
    # cnc.defineStriatum(nMSD1=0,nMSD2=0,nFS=10000,nLTS=0,nChIN=0,volumeType="cube")
    # cnc.define_striatum(num_dSPN=0, num_iSPN=0, num_FS=100, num_LTS=0, num_ChIN=0, volume_type="cube")
    # cnc.defineStriatum(nMSD1=10,nMSD2=10,nFS=10,nLTS=10,nChIN=10,volumeType="slice")

    # cnc.defineStriatum(nMSD1=500,nMSD2=500,nFS=0,nLTS=0,nChIN=500,volumeType="cube")
#    cnc.define_striatum(num_dSPN=1500, num_iSPN=1500, num_FS=0, num_LTS=0, num_ChIN=0,
#                        volume_type="cube", cell_spec_dir=cell_spec)
    cnc.define_striatum(num_dSPN=47500, num_iSPN=47500, num_FS=1300, num_LTS=0, num_ChIN=0,
                        volume_type="cube", cell_spec_dir=cell_spec)
    cnc.write_json(config_name)

    if not connect_neurons:
        print("Removing all target information, and rewriting config file")
        # Reopen the config file, and remove all connectivity settings, to
        # get an unconnected network

        with open(config_name, "r") as f:
            con_data = json.load(f, object_pairs_hook=OrderedDict)

        for k in con_data:
            if k.lower() in ["volume", "channels"]:
                continue

            # Remove targets
            if False:
                x = con_data[k]
                del x["targets"]

        with open(config_name, "w") as f:
            print("Writing to file: " + str(config_name))
            json.dump(con_data, f, indent=2)

    print("Now run:\nsnudda place " + sim_name)
    print("snudda detect " + sim_name)
    print("snudda prune " + sim_name)
    print("snudda input " + sim_name + " --input config/input-tinytest-v4.json")
    print("snudda simulate " + sim_name \
          + " --voltOut " + sim_name + "/volt-out.csv --time 10.0")
