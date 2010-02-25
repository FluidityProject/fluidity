#!/usr/bin/env python

import sys
import glob

file_name=glob.glob1('.','*_checkpoint.flml')

for f in range(len(file_name)):
  try:
    flml_file = file(file_name[f], "r")
    flml_options = flml_file.readlines()
    flml_file.close()
  except:
    sys.stderr.write("Error: change_options failed to read options file\n")
    sys.exit(1)
  
  for i in range(len(flml_options))[:len(flml_options) - 1]:
    line = flml_options[i]
    if "<current_time" in line:
      save_time = flml_options[i+1]
      flml_options[i+1] = "            <real_value rank=\"0\">0.0</real_value>\n"
  
  for i in range(len(flml_options))[:len(flml_options) - 1]:
    line = flml_options[i]
    if "<simulation_name" in line:
      flml_options[i+1] = "        <string_value lines=\"1\">"+file_name[f][:-5]+"</string_value>\n"
    elif "<finish_time" in line:
      flml_options[i+1] = save_time
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
    flml_file = file(file_name[f], "w")
    flml_file.writelines(flml_options)
    flml_file.close()
  except:
    sys.stderr.write("Error: change_options failed to write options file\n")
    sys.exit(1)
