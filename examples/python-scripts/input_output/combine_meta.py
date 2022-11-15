import json
import numpy as np
import os

from collections import OrderedDict

### ispn ###
#neuron_folder= 'str-ispn-e150908_c4_D2-m51-5-DE'
#neuron_folder= 'str-ispn-e150917_c11_D2-mWT-MSN1'
#neuron_folder= 'str-ispn-e151123_c1_D2-mWT-P270-09'
#neuron_folder= 'str-ispn-e160118_c10_D2-m46-3-DE'
###
### dspn ###
#neuron_folder= 'str-dspn-e150602_c1_D1-mWT-0728MSN01'
#neuron_folder= 'str-dspn-e150917_c6_D1-m21-6-DE'
#neuron_folder= 'str-dspn-e150917_c9_D1-mWT-1215MSN03'
#neuron_folder= 'str-dspn-e150917_c10_D1-mWT-P270-20'
###
meta_pd0 = os.path.join('..', '..', '20220930_3', 'PD0', 'neurons', 'striatum', neuron_folder[4:8], neuron_folder,'meta.json')
meta_pd2 = os.path.join('..', '..', '20220930_3', 'PD2', 'neurons', 'striatum', neuron_folder[4:8], neuron_folder,'meta.json')
meta_path=os.path.join('..', '..', '20220930_3','meta_json_files')

with open(meta_pd0, "r") as f:
    data_pd0 = json.load(f, object_pairs_hook=OrderedDict)

with open(meta_pd2, "r") as f:
    data_pd2 = json.load(f, object_pairs_hook=OrderedDict)
    
merge_data=OrderedDict()

# check that the parameters keys are the same
assert data_pd0.keys()==data_pd2.keys()
print('data.keys() checked')

for par in data_pd0.keys():
    merge_data[par]={}
    for i in range(0,8):
        merge_data[par]['var'+str(i)]=[]
    for mor_pd0 in data_pd0[par].keys():
        var=data_pd0[par][mor_pd0]['morphology'][-5:-4]
        merge_data[par]['var'+str(var)]=[mor_pd0]
    for mor_pd2 in data_pd2[par].keys():

        var=data_pd2[par][mor_pd2]['morphology'][-5:-4]
        merge_data[par]['var'+str(var)].append(mor_pd2) 
        
print('merged_data created, writing to file')

with open(os.path.join(meta_path, neuron_folder+'_merged_meta.json'), 'w') as f:
    json.dump(merge_data, f, indent=4)

