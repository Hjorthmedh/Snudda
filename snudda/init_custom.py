from snudda.init import SnuddaInit
from collections import OrderedDict
import json
import os

if __name__ == "__main__":
  
  import argparse
  parser = argparse.ArgumentParser(description="Init custom network")
  parser.add_argument("network",type=str,help="Network path")
  args = parser.parse_args()

  simName = args.network
  
  connectNeurons = False

  #simName = "networks/FSmorphTest2orig"
  #simName = "networks/FSmorphTest1b"
  #simName = "LTStest"
  #simName = "networks/twoFS"
  #simName = "networks/FSmorphTest4"
  #simName = "networks/3types-v2"
  # simName = "networks/SynTest-v6" # MSMS tuning
  #simName = "networks/SynTest-v15"  

  configName= simName + "/network-config.json"
  cnc = SnuddaInit(structDef={},configName=configName,nChannels=1)
  #cnc.defineStriatum(nMSD1=500,nMSD2=500,nFS=0,nLTS=0,nChIN=30,volumeType="cube")  
  #cnc.defineStriatum(nMSD1=120,nMSD2=120,nFS=20,nLTS=0,nChIN=0,volumeType="slice")
  #cnc.defineStriatum(nMSD1=0,nMSD2=0,nFS=10000,nLTS=0,nChIN=0,volumeType="slice")
  #cnc.defineStriatum(nMSD1=0,nMSD2=0,nFS=10000,nLTS=0,nChIN=0,volumeType="cube")
  cnc.defineStriatum(nMSD1=0,nMSD2=0,nFS=100,nLTS=0,nChIN=0,volumeType="cube")
  #cnc.defineStriatum(nMSD1=10,nMSD2=10,nFS=10,nLTS=10,nChIN=10,volumeType="slice")
  
  # cnc.defineStriatum(nMSD1=500,nMSD2=500,nFS=0,nLTS=0,nChIN=500,volumeType="cube")  



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
      

  print("Now run:\nsnudda place " + simName)
  print("snudda detect " + simName)
  print("snudda prune " + simName)
  print("snudda input " + simName + " --input config/input-tinytest-v4.json")
  print("snudda simulate " + simName \
        + " --voltOut " + simName + "/volt-out.csv --time 10.0")
