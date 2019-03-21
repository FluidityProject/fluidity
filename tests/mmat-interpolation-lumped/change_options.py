#!/usr/bin/env python

import sys

file_name = ["mmat-interpolation_2_checkpoint.flml"]

for f in range(len(file_name)):
  try:
    flml_file = open(file_name[f], "r")
    flml_options = flml_file.readlines()
    flml_file.close()
  except:
    sys.stderr.write("Error: change_options failed to read options file\n")
    sys.exit(1)
  
  for i in range(len(flml_options))[:len(flml_options) - 1]:
    line = flml_options[i]
    if "<interpolate_at_first_timestep" in line:
      flml_options[i] = ""
  
  for i in range(len(flml_options)):
    if "<checkpointing>" in flml_options[i]:
      while not "</checkpointing>" in flml_options[i] and i < len(flml_options):
        flml_options[i] = ""
        i += 1
      if i < len(flml_options):
        flml_options[i] = ""
      break
  
  try:
    flml_file = open(file_name[f], "w")
    flml_file.writelines(flml_options)
    flml_file.close()
  except:
    sys.stderr.write("Error: change_options failed to write options file\n")
    sys.exit(1)
  


