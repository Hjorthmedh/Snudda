import json

def define_input_config(input_path,input_type):
    
    synapse_data = json.load(open(input_path,"r"))["0"][input_type]
    del synapse_data["start"]
    del synapse_data["end"]
    del synapse_data["frequency"]
    del synapse_data["nInputs"]
    
    with open("input_definition.json", "w") as f:
        json.dump(synapse_data,f)