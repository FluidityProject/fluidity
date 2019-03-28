#!/usr/bin/env python

import sys

file_name = ["mmat-interpolation_2_checkpoint.flml"]

append_options = []
append_options.append("            <adapt_at_first_timestep>\n")
append_options.append("              <number_of_adapts>\n")
append_options.append("                <integer_value rank=\"0\">2</integer_value>\n")
append_options.append("              </number_of_adapts>\n")
append_options.append("            </adapt_at_first_timestep>\n")


for f in range(len(file_name)):
  try:
    flml_file = open(file_name[f], "r")
    flml_options = flml_file.readlines()
    flml_file.close()
  except:
    sys.stderr.write("Error: change_options failed to read options file\n")
    sys.exit(1)
  
  for i in range(len(flml_options)):
    if "<checkpointing>" in flml_options[i]:
      while not "</checkpointing>" in flml_options[i] and i < len(flml_options):
        flml_options[i] = ""
        i += 1
      if i < len(flml_options):
        flml_options[i] = ""
      break

  for i in range(len(append_options)):
    flml_options.append("\n")

  for i in range(len(flml_options)-1,-1,-1):
    if "</hr_adaptivity>" in flml_options[i]:
      flml_options[i+len(append_options):len(flml_options)-len(append_options)] = flml_options[i:-len(append_options)]
      for j in range(len(append_options)):
        flml_options[i+j] = append_options[j]
      break
  
  try:
    flml_file = open(file_name[f], "w")
    flml_file.writelines(flml_options)
    flml_file.close()
  except:
    sys.stderr.write("Error: change_options failed to write options file\n")
    sys.exit(1)
  


