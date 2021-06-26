#!/usr/bin/env python3

import sys

file_name = ["explicit-hyperc-superman-conservative_121_checkpoint.flml"]

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
      flml_options[i+1] = "            <real_value rank=\"0\">3.0</real_value>\n"
    elif "Do not remove this comment" in line:
      flml_options[i] = """                        <string_value lines="20" type="python">def val(X,t):&#x0A;   u =  0.25 * ((4*X[0] - 2) + (4*X[1] - 2)**3)&#x0A;   v = -0.25 * ((4*X[1] - 2) + (4*X[0] - 2)**3)&#x0A;   return [-u, -v]</string_value>"""
  
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
  


