#!/usr/bin/env python3

import sys

file_name = ["explicit-hyperc-shear-conservative_16_checkpoint.flml"]

for f in range(len(file_name)):
  try:
    with open(file_name[f]) as flml_file:
      flml_options = flml_file.readlines()
  except:
    sys.stderr.write("Error: change_options failed to read options file\n")
    sys.exit(1)
  
  for i in range(len(flml_options))[:len(flml_options) - 1]:
    line = flml_options[i]
    if "<current_time" in line:
      flml_options[i+1] = "            <real_value rank=\"0\">0.0</real_value>\n"
    elif "<finish_time" in line:
      flml_options[i+1] = "            <real_value rank=\"0\">7.7176909321518821</real_value>\n"
    elif "Shear rotation about origin." in line:
      flml_options[i] = "                        <string_value lines=\"20\" type=\"python\">def val(X,t):&#x0A;   from math import sin, cos&#x0A;   # Shear rotation about origin.&#x0A;   return (-1.0*sin(X[0])*cos(X[1]), cos(X[0])*sin(X[1]))</string_value>\n"
  
  for i in range(len(flml_options)):
    if "<checkpointing>" in flml_options[i]:
      while not "</checkpointing>" in flml_options[i] and i < len(flml_options):
        flml_options[i] = ""
        i += 1
      if i < len(flml_options):
        flml_options[i] = ""
      break
  
  try:
    with open(file_name[f], "w") as flml_file:
      flml_file.writelines(flml_options)
  except:
    sys.stderr.write("Error: change_options failed to write options file\n")
    sys.exit(1)
  


