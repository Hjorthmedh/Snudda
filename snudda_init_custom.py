from snudda_init import SnuddaInit
from collections import OrderedDict
import json
import os

if __name__ == "__main__":

  connectNeurons = False

  #simName = "networks/FSmorphTest2orig"
  #simName = "networks/FSmorphTest1b"
  #simName = "LTStest"
  #simName = "networks/twoFS"
  #simName = "networks/FSmorphTest4"
  #simName = "networks/3types-v2"
  simName = "networks/SynTest-v5"

  
  cnc = SnuddaInit(structDef={},configName=None,nChannels=1)
  # cnc.defineStriatum(nMSD1=20,nMSD2=20,nFS=10,nLTS=0,nChIN=0,volumeType="cube")  
  cnc.defineStriatum(nMSD1=100,nMSD2=100,nFS=50,nLTS=0,nChIN=0,volumeType="slice")  

  configName= simName + "/network-config.json"

  dirName = os.path.dirname(configName)
  
  if not os.path.exists(dirName):
    os.makedirs(dirName)

  cnc.writeJSON(configName)

  if(not connectNeurons):
    print("Removing all target information, and rewriting config file")
    # Reopen the config file, and remove all connectivity settings, to
    # get an unconnected network

    with open(configName,"r") as f:
      conData = json.load(f,object_pairs_hook=OrderedDict)

    for k in conData:
      if(k.lower() in ["volume", "channels"]):
        continue

      # Remove targets
      if(False):
        x = conData[k]
        del x["targets"]
      
    with open(configName,"w") as f:
      print("Writing to file: " + str(configName))
      json.dump(conData,f,indent=2)
      

  print("Now run:\n./snudda.py place " + simName)
  print("./snudda.py detect " + simName)
  print("./snudda.py prune " + simName)
  print("./snudda.py input " + simName + " --input config/input-tinytest-v4.json")
  print("./snudda.py simulate " + simName \
        + " --voltOut " + simName + "/volt-out.csv --time 10.0")
