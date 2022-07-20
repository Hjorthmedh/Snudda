#!/usr/bin/env python
# coding: utf-8

from check_features_traces import load_voltage_data
from input_config import *
import os
import json
import numpy as np
from snudda.utils import SnuddaLoad
from input_tuning_analysis import post_analysis
from number_input import *
from add_feature import *
from count_passing import *
from argparse import ArgumentParser, RawTextHelpFormatter

#########################################################################
class NumpyEncoder(json.JSONEncoder):
   def default(self, obj):
       if isinstance(obj, np.integer):
           return int(obj)
       elif isinstance(obj, np.floating):
           return float(obj)
       elif isinstance(obj, np.ndarray):
           return obj.tolist()
       else:
           return json.JSONEncoder.default(self, obj)
#########################################################################  
parser = ArgumentParser("Analyse data and add to meta.json", formatter_class=RawTextHelpFormatter)   
parser.add_argument("networkPath", help="Network path") 
parser.add_argument("--networkName", help="Network name")
parser.add_argument("--neuronsPath", help="Neurons path")     
parser.add_argument("--input_type", help="Type of external input",choices=["thalamic", "cortical"], default="thalamic")
parser.add_argument("--neuronType", type=str, help="E.g. dspn")
parser.add_argument("--singleNeuronType", default=None, type=str,help="We want to simulate one neuron subtype, eg. FS_1")
parser.add_argument("--snudda_data_tools", default=os.path.join('BasalGangliaData','tools'), type=str,help="")

args = parser.parse_args()
BG_data = os.environ["SNUDDA_DATA"]

#########################################################################
#prep_choose_input
#########################################################################
data = SnuddaLoad(os.path.join(args.networkPath, "network-synapses.hdf5"))
neurons = data.data["neurons"]
neuron_dict ={str(k["neuronID"]): {'parameterKey': k['parameterKey'],
  'morphologyKey': k['morphologyKey']} for k in neurons}

with open(f"neuron_config_{args.networkName}.json","w") as f:
    json.dump(neuron_dict,f,cls=NumpyEncoder)

input_path = os.path.join(args.networkPath,"input_config.json")
define_input_config(input_path,args.input_type)

#########################################################################
#Analysis
#########################################################################
#
#get_ipython().run_cell_magic('time', '', 'result, voltages = post_analysis(args.networkPath)')

result, voltages=post_analysis(args.networkPath,args.snudda_data_tools)

with open(os.path.join(args.networkPath,"input_results.json"),"w") as f:
    json.dump(result, f, indent=4, sort_keys=True,cls=NumpyEncoder)

with open(os.path.join(args.networkPath,"voltages_passed.json"),"w") as f:
    json.dump(voltages, f, indent=4, sort_keys=True,cls=NumpyEncoder)
    
#########################################################################
#Passing models
#########################################################################
#Create input_map_{model}json using neuron_config_{model}.json
count_passing(args.networkName, args.networkPath, args.input_type)

#########################################################################
#Number inputs
#########################################################################
number_inputs(args.networkName, args.networkPath,args.input_type)
voltage_values(args.networkName, args.networkPath,args.input_type)

#########################################################################
#Add feature
#########################################################################
neuron_path = os.path.join(BG_data,"neurons","striatum",args.neuronType,args.singleNeuronType)
#Creates meta.json that is based on  input_map_{model}.json,  inputs_population_{model}.json, and input_definition.json
add_feature(args.networkName,args.input_type,neuron_path, ratio=0.5)




