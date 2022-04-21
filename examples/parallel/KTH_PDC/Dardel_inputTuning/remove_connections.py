import json,os,sys
from collections import OrderedDict

network_path = sys.argv[1]#os.path.join("networks","testnetwork19")
config_name=os.path.join(network_path, "network-config.json")

with open(config_name, "r") as f:
    con_data = json.load(f, object_pairs_hook=OrderedDict)
con_data["Connectivity"] = {}
with open(config_name, "w") as f:
    print("Writing to file: " + str(config_name))
    json.dump(con_data, f, indent=2)