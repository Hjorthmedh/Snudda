#!/usr/bin/env python
import sys

if len(sys.argv) <= 1:
    print("Provide name of file to purge last entry in LDFLAGS from")
    sys.exit(-1)
    
with open( sys.argv[1]  ) as file:
    for line in file:
        if "LDFLAGS = $(LINKFLAGS) $(UserLDFLAGS)" in line:
            print(",".join(line.split(",")[:-1]))
        else:
            print(line, end="")
