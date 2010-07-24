#!/usr/bin/env python

from optparse import OptionParser
import glob
import os 
import sys

optparser=OptionParser(usage='usage: %prog <filename>',
                       add_help_option=True,
                       description="""This takes a flml file and updates the mixing_bin_bounds option""" +
                                   """in line with the schema option to specify with either a constant or a pytion function""")

(options, argv) = optparser.parse_args()

if len(argv)<1:
  optparser.print_help()
  sys.exit(1)

filename = argv[0]

try:
  flml_file = file(filename, 'r')
  flml_options = flml_file.readlines()
  flml_file.close()
except:
  sys.stderr.write('Error: failed to read options file\n')
  sys.exit(1)
  
for i in range(len(flml_options)-1):
  if ('<mixing_bin_bounds>' in flml_options[i]) and ('<constant>' in flml_options[i+1]): 
    sys.stderr.write('Warning: file already fully or partially updated\n no changes made to flml')
    sys.exit(1) 

for i in range(len(flml_options)):
  if ('<mixing_bin_bounds>' in flml_options[i]):
    white_space = flml_options[i].split('<mixing_bin_bounds>')[0] 
    flml_options.insert(i+1,white_space+'  <constant>\n')
    j = i
    while ('</mixing_bin_bounds>' not in flml_options[j]): j = j+1
    white_space = flml_options[j].split('</mixing_bin_bounds>')[0] 
    flml_options.insert(j,white_space+'  </constant>\n')

try:
  flml_file = file(filename, 'w')
  flml_file.writelines(flml_options)
  flml_file.close()
except:
  sys.stderr.write('Error: failed to write options file\n')
  sys.exit(1)
