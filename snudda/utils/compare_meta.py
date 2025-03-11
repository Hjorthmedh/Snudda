import json
import argparse

def check_missing(dict_a, dict_b, file_b):
  for k in dict_a.keys():
    if k not in dict_b:
      print(f"Parameter key {k} missing in {file_b}")
    else:
      for k2 in dict_a[k].keys():
        if k2 not in dict_b[k]:
          print(f"Morphology key {k2} missing for parameter set {k} in {file_b}")
      
def cli():
  parser = argparse.ArgumentParser()
  parser.add_argument("file_a")
  parser.add_argument("file_b")

  args = parser.parse_args()

  with open(args.file_a, "r") as f:
    dict_a = json.load(f)

  with open(args.file_b, "r") as f:
    dict_b = json.load(f)

  print(f"Comparing {args.file_a} and {args.file_b}")
    
  check_missing(dict_a, dict_b, args.file_b)
  check_missing(dict_b, dict_a, args.file_a)

  
cli()
