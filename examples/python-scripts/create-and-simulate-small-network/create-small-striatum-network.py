import os
network_path = os.path.join("networks","testnetwork1")
input_config_path=os.path.join("..","..","..","snudda","data","input_config","input-v10-scaled.json")

doSnuddaInit=1
doSnuddaPlace=1
doSnuddaDetect=1
doSnuddaPrune=1
doSnuddaInput=1

if doSnuddaInit:
    print('\nSTARTING SnuddaInit')
    from snudda import SnuddaInit
    struct_def = {"Striatum": 10}
    si = SnuddaInit(network_path=network_path, struct_def=struct_def, random_seed=123)

if doSnuddaPlace:
    print('\nSTARTING SnuddaPlace')
    from snudda import SnuddaPlace
    spl = SnuddaPlace(network_path=network_path)
    spl.place()

if doSnuddaDetect:
    print('\nSTARTING SnuddaDetect')
    from snudda import SnuddaDetect
    sd = SnuddaDetect(network_path=network_path)
    sd.detect()

if doSnuddaPrune:
    print('\nSTARTING SnuddaPrune')
    from snudda import SnuddaPrune
    sp = SnuddaPrune(network_path=network_path)
    sp.prune()

if doSnuddaInput:
    print('\nSTARTING SnuddaInput')
    from snudda.input import SnuddaInput
    sii = SnuddaInput(network_path=network_path,input_config_file=input_config_path,verbose=False)
    sii.generate()