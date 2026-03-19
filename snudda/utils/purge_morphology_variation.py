#!/usr/bin/env python3

# This is useful if the user wants to purge a particular morphology variation
# from the meta.json file.

import json
import argparse
from pathlib import Path

def filter_morphologies(input_file, output_file, filter_string):
    """Remove morphology keys where morphology contains the filter string."""
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    # Filter out morphology keys with matching morphologies
    for param_key in list(data.keys()):
        for morph_key in list(data[param_key].keys()):
            morphology = data[param_key][morph_key].get('morphology', '')
            if filter_string in morphology:
                del data[param_key][morph_key]
        
        # Remove parameter key if empty after filtering
        if not data[param_key]:
            del data[param_key]
    
    # Write filtered data
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=4)
    
    print(f"Filtered data written to {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Remove morphology variations from meta.json file'
    )
    parser.add_argument('input', help='Input meta.json file')
    parser.add_argument('output', help='Output meta.json file')
    parser.add_argument('filter', help='String to filter (e.g., "var1")')
    
    args = parser.parse_args()
    filter_morphologies(args.input, args.output, args.filter)
