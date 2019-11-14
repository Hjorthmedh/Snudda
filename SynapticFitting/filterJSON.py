import json
import os

# M1LH --- NO CHAT?
# TH -- NO LTS (no input from Thalamus to LTS)

dataPath = "DATA/YvonneJohansson2019/"

fileInfo = [ ("M1LH_Analysis_191001.h5-parameters.json",
              "M1LH_Analysis_191001.h5-parameters-FS.json",
              [8,9,97,12,111,2,114,61,116,138,6,63,64]),
             ("M1LH_Analysis_191001.h5-parameters.json",
              "M1LH_Analysis_191001.h5-parameters-MS.json",
              [122,131,20,53,40,132,59,96,41,54,10,73,42,55,11,99,100,13,57,101,78,31,43,117,45,3,58,34,15,47,16,5,107,62,92,37,120,129,72,38,130,51,140]),
             ("M1LH_Analysis_191001.h5-parameters.json",
              "M1LH_Analysis_191001.h5-parameters-LTS.json",
              [21,4,48]),
             
             ("M1RH_Analysis_190925.h5-parameters.json",
              "M1RH_Analysis_190925.h5-parameters-FS.json",
              [161,54,212,88,176,151,179,181,153,91,168,156,157,171,303,173]),
             ("M1RH_Analysis_190925.h5-parameters.json",
              "M1RH_Analysis_190925.h5-parameters-MS.json",
              [318,184,183,120,301,37,47,313,45,70,188,10,175,27,134,56,190,99,316,192,193,137,195,138,1,101,113,141,302,334,90,80,114,154,209,328,39,66,3,335,155,329,40,13,4,104,169,21,5,105,131,170,310,331,50,14,22,6,122,158,311,51,185,23,33,84,119,159,172,43,95,16,85,24,106,124,147,44,53,25,9,187,86]),
             ("M1RH_Analysis_190925.h5-parameters.json",
              "M1RH_Analysis_190925.h5-parameters-LTS.json",
              [308,42,69,55,57,194,49,41,52,64]),
             ("M1RH_Analysis_190925.h5-parameters.json",
              "M1RH_Analysis_190925.h5-parameters-CHAT.json",
              [125,26,189,126,111,121,62]),

             ("S1_Analysis_191001.h5-parameters.json",
              "S1_Analysis_191001.h5-parameters-FS.json",
              [36,345,18,232,43,139,53,225,100,227,167,13,133,21,58,160,379,38,4,439,27,390,51,52]),
             ("S1_Analysis_191001.h5-parameters.json",
              "S1_Analysis_191001.h5-parameters-MS.json",
              [144,162,226,434,47,64,72,223,240,218,328,339,393,155,140,191,340,219,129,261,192,330,341,262,143,319,182,251,375,80,108,16,173,19,233,40,61,68,73,120,252,110,115,131,295,320,287,5,69,101,132,333,377,76,82,106,234,112,159,254,102,177,203,306,378,2,22,268,352,372,150,178,185,307,103,113,117,255,335,84,26,269,63,67,124,179,205,114,135,161,297,317,324,230,354,118,257,298,180,325,277,258,299,119,137,391,222,239,392,138,154,217,96,127,181,259,327,338,87,271]),
             ("S1_Analysis_191001.h5-parameters.json",
              "S1_Analysis_191001.h5-parameters-LTS.json",
              [171,56]),
             ("S1_Analysis_191001.h5-parameters.json",
              "S1_Analysis_191001.h5-parameters-CHAT.json",
              [208,353,188]),
             
             ("TH_Analysis_191001.h5-parameters.json",
              "TH_Analysis_191001.h5-parameters-FS.json",
              [100,74,62,110,181,65,21,103,23,47,27] ),
             ("TH_Analysis_191001.h5-parameters.json",
              "TH_Analysis_191001.h5-parameters-MS.json",
              [122,186,195,20,226,81,225,15,172,60,205,90,173,61,107,206,50,17,207,175,208,63,179,180,66,182,183,67,184,151,164,101,113,154,196,217,44,185,165,53,8,102,155,150,153,144,147,166,54,9,115,156,198,22,228,10,145,167,157,189,220,46,84,168,56,117,200,30,12,57,118,191,201,222,231,25,86,13,169,71,127,223,48,5,58,106,120,224,6,88,14,171,59] ),
             ("TH_Analysis_191001.h5-parameters.json",
              "TH_Analysis_191001.h5-parameters-CHAT.json",
              [192,146,7,89,49,174,91,51,93,95,45,237,230,85,31,240,87] ) ]

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

