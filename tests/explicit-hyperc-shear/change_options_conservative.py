#!/usr/bin/env python3

file_name = ["explicit-hyperc-shear-conservative_16_checkpoint.flml"]

for f in range(len(file_name)):
    with open(file_name[f]) as flml_file:
        flml_options = flml_file.readlines()

    for i in range(len(flml_options))[: len(flml_options) - 1]:
        line = flml_options[i]
        if "<current_time" in line:
            flml_options[i + 1] = '            <real_value rank="0">0.0</real_value>\n'
        elif "<finish_time" in line:
            flml_options[
                i + 1
            ] = '            <real_value rank="0">7.7176909321518821</real_value>\n'
        elif "Shear rotation about origin." in line:
            flml_options[
                i
            ] = """
            <string_value type="code" language="python" lines="20">def val(X,t):
   from math import sin, cos
   # Shear rotation about origin.
   return (-1.0*sin(X[0])*cos(X[1]), cos(X[0])*sin(X[1]))</string_value>
"""

    for i in range(len(flml_options)):
        if "<checkpointing>" in flml_options[i]:
            while "</checkpointing>" not in flml_options[i] and i < len(flml_options):
                flml_options[i] = ""
                i += 1
            if i < len(flml_options):
                flml_options[i] = ""
            break

    with open(file_name[f], "w") as flml_file:
        flml_file.writelines(flml_options)
