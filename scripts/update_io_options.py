#!/usr/bin/env python

from optparse import OptionParser
import glob
import os 
import sys

optparser=OptionParser(usage='usage: %prog <filename>',
                       add_help_option=True,
                       description="""This takes a flml file and updates the io options """ +
                                   """for dump period to the new schema""")

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

for i in range(len(flml_options)):
  if ('<dump_period' in flml_options[i]):
    white_space = flml_options[i].split('<dump_period')[0]
    if ('integer' in flml_options[i+1]):
      flml_options[i+1]='  '+flml_options[i+1]
      flml_options.insert(i+1,white_space+'  <constant>\n')
      flml_options.insert(i+3,white_space+'  </constant>\n')
    elif ('real' in flml_options[i+1]):
      flml_options[i] = flml_options[i].split(' replaces')[0]+'>\n'
      flml_options[i+1]='  '+flml_options[i+1]
      flml_options.insert(i+1,white_space+'  <constant replaces="TIMDUM">\n')
      flml_options.insert(i+3,white_space+'  </constant>\n')
      #dump_period = float(flml_options[i+1].split('>')[-2].split('<')[0])
    else:
      sys.stderr.write('Error: dump_period or dump_period_in_timesteps not real or integer, failed to update\n')
      sys.exit(1)

try:
  flml_file = file(filename, 'w')
  flml_file.writelines(flml_options)
  flml_file.close()
except:
  sys.stderr.write('Error: change_options failed to write options file\n')
  sys.exit(1)
