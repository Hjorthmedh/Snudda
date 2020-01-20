import json
import os


dataPath = "DATA/Planert2010/d1d2conns/traces/"

fileInfo = [ ("trace_table.txt-parameters.json",
              "trace_table.txt-parameters-FI.json",
              [4,17,18,21]),
             ("trace_table.txt-parameters.json",
              "trace_table.txt-parameters-ID.json",
              [0,1,10,11,12,13,15,16]),
             ("trace_table.txt-parameters.json",
              "trace_table.txt-parameters-FD.json",
              [2,5,6]),
             ("trace_table.txt-parameters.json",
              "trace_table.txt-parameters-II.json",
              [3,9,19,20,22]),
             ("trace_table.txt-parameters.json",
              "trace_table.txt-parameters-LI.json",
              [7]),
             ("trace_table.txt-parameters.json",
              "trace_table.txt-parameters-DD.json",
              [8]),             
             ("trace_table.txt-parameters.json",
              "trace_table.txt-parameters-DI.json",
              [14])             
]

for fi in fileInfo:
  print("Reading from " + str(fi[0]))
  with open(dataPath + fi[0],"r") as fIn:
    dataIn = json.load(fIn)

    dataOut = dict()

    print("Keeping: ", end="")
    for key in dataIn:
      if int(key) in fi[2]:
        dataOut[key] = dataIn[key]
        print(str(key),end=" ")
    print("")
    
    assert not os.path.exists(dataPath + fi[1]), "File already exists " + dataPath + fi[1]

    print("Writing to " + str(dataPath + fi[1]))
    with open(dataPath + fi[1],"w") as fOut:
      json.dump(dataOut,fOut,indent=2)

