import json
import os
import numpy as np
from datetime import datetime
from collections import OrderedDict
import numpy as np
from neuron import hoc
from snudda.neurons.neuron_model_extended import NeuronModel
from snudda.utils.numpy_encoder import NumpyEncoder

"""
    
    Functions for dumping the contents of a SnuddaSimulate Instance

"""

def extract_vec_stim(vec_stim):

    temp = dict()
    for n in dir(vec_stim):

        if "__" not in n:

            if n in ["loc", "same", "get_loc", "ptr"]:
                pass
                # causes hoc error
            elif isinstance(getattr(vec_stim, n), hoc.HocObject):
                temp.update({n: getattr(vec_stim, n)()})
            else:
                temp.update({n: getattr(vec_stim, n)})

    return temp

def extract_net_con(net_con):
    temp = dict()
    for n in dir(net_con):

        if "__" not in n:

            if isinstance(getattr(net_con, n), hoc.HocObject):

                if n in ["event", "same", "preloc", "postloc"]:
                    pass
                    # causes hoc error
                elif n in ["get_recordvec"]:
                    if getattr(net_con, n)():
                        temp.update({n: getattr(net_con, n)().as_numpy()})
                    else:
                        temp.update({n: getattr(net_con, n)()})
                elif n in ["weight"]:
                    temp.update({n: getattr(net_con, n)[0]})
                else:
                    temp.update({n: getattr(net_con, n)()})
            else:
                temp.update({n: getattr(net_con, n)})
    return temp


def extract_synapse(synapse):
    temp = dict()
    for n in dir(synapse):

        if "__" not in n:

            if n in ["loc", "same", "get_loc"]:
                pass
                # causes hoc error
            elif isinstance(getattr(synapse, n), hoc.HocObject):
                temp.update({n: getattr(synapse, n)()})
            else:
                temp.update({n: getattr(synapse, n)})

    return temp


def extract_external_stim(external_stim):

    temp = dict()
    for num, data in external_stim.items():
        temp.update({num: dict()})
        for i, d in enumerate(data):
            temp[num].update({i: dict()})
            temp[num][i].update(Vector=list(d[0]))
            temp[num][i].update(VecStim=extract_vec_stim(vec_stim=d[1]))
            temp[num][i].update(NetCon=extract_net_con(net_con=d[2]))
            temp[num][i].update(Synapse=extract_synapse(synapse=d[3]))
            temp[num][i].update(SpikeTime=d[4])

    return temp


def extract_synapse_dict(synapse_dict):
    temp = dict()
    k = 0
    for name, data in synapse_dict.items():

        temp.update({k: dict()})
        temp[k].update(tuple=name)
        temp[k].update(data=dict())
        for i, sdata in enumerate(data):
            temp[k]["data"].update({i: dict()})
            temp[k]["data"][i].update({"synapse": extract_synapse(synapse=sdata[0])})
            temp[k]["data"][i].update({"synapse": extract_net_con(net_con=sdata[1])})
            temp[k]["data"][i].update({"synapse": sdata[2]})
            temp[k]["data"][i].update({"synapse": sdata[3]})

    return temp

def extract_synapse_list(synapse_list):
    temp = dict()

    for i, synapse in enumerate(synapse_list):
        temp.update({i: extract_synapse(synapse=synapse)})

    return temp


def extract_net_con_list(net_con_list):
    temp = dict()
    for i, item in enumerate(net_con_list):

        temp.update({i: extract_net_con(net_con=item)})

    return temp


def extract_neuron_data(neurons):
    temp = dict()
    for neuron_id, data in neurons.items():

        temp.update({neuron_id: dict()})
        if isinstance(data, NeuronModel):

            for name, sdata in data.__dict__.items():

                if name == "morphology":
                    temp[neuron_id].update({name: sdata.to_dict()})
                elif name == "mechanisms":
                    to_dict_data = {i: s.to_dict() for i, s in enumerate(sdata)}
                    temp[neuron_id].update({name: to_dict_data})
                elif name == "params":
                    to_dict_data = {pname: nrnobject.to_dict() for pname, nrnobject in sdata.items()}
                    temp[neuron_id].update({name: to_dict_data})
                else:
                    temp[neuron_id].update({name: data})
        else:
            raise TypeError(" Not NeuronModel type")

    return temp


def dump(snudda_simulate_obj):
    date_time = datetime.now()
    d = date_time.strftime("%Y-%m-%d-%H-%M-%S")
    os.makedirs(os.path.join(snudda_simulate_obj.network_path, "dump"), exist_ok=True)

    temp = dict()

    temp.update(types=dict())
    for obj in dir(snudda_simulate_obj):
        temp["types"].update({obj: str(type(getattr(snudda_simulate_obj, obj)))})

    temp.update(__dict__=dict())
    # __dict__ is A dictionary or other mapping object used to store an object’s (writable) attributes
    for name, data in snudda_simulate_obj.__dict__.items():

        # Here we have to cover all the cases represented in the data

        if name == "neurons":
            packaged = extract_neuron_data(neurons=data)
            temp["__dict__"].update({name: packaged})
        elif name == "net_con_list":
            packaged = extract_net_con_list(net_con_list=data)
            temp["__dict__"].update({name: packaged})
        elif name == "synapse_list":
            packaged = extract_synapse_list(synapse_list=data)
            temp["__dict__"].update({name: packaged})
        elif name == "synapse_dict":
            packaged = extract_synapse_dict(synapse_dict=data)
            temp["__dict__"].update({name: packaged})
        elif name == "external_stim":
            packaged = extract_external_stim(external_stim=data)
            temp["__dict__"].update({name: packaged})
        elif name == "i_stim":
            pass

        elif name == "v_clamp_list":
            pass

        elif name == "gap_junction_list":
            pass

        elif name == "check_id_recordings":
            packaged = {v[0]: list(v[1]) for v in data}
            temp["__dict__"].update({name: packaged})

        elif name == "config":
            pass
            print("Fix")
        elif name == "network_info":
            pass
            print("Fix")
        elif name == "connectivityDistributions":
            packaged = dict()
            for connection, information in data.items():
                packaged.update(connection=connection)
                packaged.update(data=information)
            temp["__dict__"].update({name: packaged})
        elif isinstance(data, dict):
            temp["__dict__"].update({name: data})
        elif isinstance(data, bool):
            temp["__dict__"].update({name: data})
        elif isinstance(data, str):
            temp["__dict__"].update({name: data})
        elif isinstance(data, np.ndarray):
            temp["__dict__"].update({name: data})
        elif isinstance(data, int):
            temp["__dict__"].update({name: data})
        elif isinstance(data, OrderedDict):
            temp["__dict__"].update({name: data})
        elif isinstance(data, list):
            temp["__dict__"].update({name: data})
        elif isinstance(data, range):
            temp["__dict__"].update({name: list(data)})
        elif isinstance(data, float):
            temp["__dict__"].update({name: data})

        else:

            print(f" No dump function for {name}")

    # TODO: replace json with hdf5 -- size will be an issue for large networks
    with open(os.path.join(snudda_simulate_obj.network_path, "dump", f"content_{d}.json"), "w") as f:
        json.dump(temp, f, cls=NumpyEncoder)
