from .init import SnuddaInit
from collections import OrderedDict
import json
import os

if __name__ == "__main__":
  
  import argparse
  parser = argparse.ArgumentParser(description="Init custom network")
  parser.add_argument("network",type=str,help="Network path")
  parser.add_argument("level",choices=["0","1","2","3","4"],type=str,
                      help="Level of Parkinsson, 0=none, 4=bad")
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

  cellSpecDir = "cellspecs.parkinson/" +str(args.level) + "/"
  
  configName= simName + "/network-config.json"
  cnc = SnuddaInit(structDef={},configName=configName,nChannels=1)
  cnc.defineStriatum(nMSD1=1500,nMSD2=1500,nFS=0,nLTS=0,nChIN=0,
                     cellSpecDir=cellSpecDir,
                     volumeType="cube")  

  # cnc.defineStriatum(nMSD1=0,nMSD2=0,nFS=100,nLTS=100,nChIN=0,volumeType="slice")



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
