#!/usr/bin/python

import sys
import subprocess


if len(sys.argv) != 2:
    print("You did not specify an input file. Running for the default Cloudy.\n")
    input_file_name = "all_species.txt"    
else:
    input_file_name = str(sys.argv[1])
    

print(input_file_name)

input_file = open(input_file_name,"r")
for current_line in input_file:
        tokens = current_line.split()
        
        if len(tokens) >= 2:
            species = tokens[0]
            num_level_limit = tokens[1]
        elif len(tokens) == 1:
            species = tokens[0]
            num_level_limit = 1000
        else:
            print(f"File {input_file_name} is not properly formatted")
            sys.exit(99)
        
        print(f"\nRunning NIST Extractor on {species} for {num_level_limit} levels")
        
        if sys.platform.startswith('win'):
            rcode = subprocess.call(
                ["python", "NistExtractor.py", species, str(num_level_limit)],
                shell=True
            )
        else:
            rcode = subprocess.call(
                ["python3", "NistExtractor.py", species, str(num_level_limit)]
            )
            



print("\nNEWrapper finished!")
    
