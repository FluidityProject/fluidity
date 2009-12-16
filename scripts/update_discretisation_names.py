#!/usr/bin/env python

from optparse import OptionParser
import glob
import os 
import sys

optparser=OptionParser(usage='usage: %prog <filename>',
                       add_help_option=True,
                       description="""This takes a flml file and updates the discretisation """ +
                                   """names to the new system.""")

(options, argv) = optparser.parse_args()

if len(argv)<1:
  optparser.print_help()
  sys.exit(1)

filename = argv[0]

try:
  flml_file = file(filename, "r")
  flml_options = flml_file.readlines()
  flml_file.close()
except:
  sys.stderr.write("Error: failed to read options file\n")
  sys.exit(1)

for i in range(len(flml_options)):
  if (("<vector_field" in flml_options[i]) and ('name="Velocity"' in flml_options[i])):
      while not "</vector_field>" in flml_options[i] and i < len(flml_options):
        if('<spatial_discretisation>' in flml_options[i]):
          while not "</spatial_discretisation>" in flml_options[i] and i < len(flml_options):
            if('<continuous_galerkin replaces="DISOPT">' in flml_options[i]):
              while not "</continuous_galerkin>" in flml_options[i] and i < len(flml_options):
                if('<continuous_galerkin replaces="DISOPT">' in flml_options[i]):
                  flml_options[i]=flml_options[i].replace("continuous_galerkin", "legacy_continuous_galerkin")
                i = i + 1
              if("</continuous_galerkin>" in flml_options[i]):
                flml_options[i]=flml_options[i].replace("continuous_galerkin", "legacy_continuous_galerkin")
            if('<continuous_galerkin_test>' in flml_options[i]):
              while not "</continuous_galerkin_test>" in flml_options[i] and i < len(flml_options):
                if('<continuous_galerkin_test>' in flml_options[i]):
                  flml_options[i]=flml_options[i].replace("continuous_galerkin_test", "continuous_galerkin")
                i = i + 1
              if("</continuous_galerkin_test>" in flml_options[i]):
                flml_options[i]=flml_options[i].replace("continuous_galerkin_test", "continuous_galerkin")
            if("<continuous_galerkin_test/>" in flml_options[i]):
              flml_options[i]=flml_options[i].replace("continuous_galerkin_test", "continuous_galerkin")
            i = i + 1
        if('<enforce_discrete_properties>' in flml_options[i]):
          while not "</enforce_discrete_properties>" in flml_options[i] and i < len(flml_options):
            i = i + 1 # do nothing
        i = i + 1

for i in range(len(flml_options)):
  if (("<scalar_field" in flml_options[i]) and (not 'name="Pressure"' in flml_options[i])):
      while not "</scalar_field>" in flml_options[i] and i < len(flml_options):
        if('<spatial_discretisation>' in flml_options[i]):
          while not "</spatial_discretisation>" in flml_options[i] and i < len(flml_options):
            if('<continuous_galerkin replaces="DISOTT">' in flml_options[i]):
              while not "</continuous_galerkin>" in flml_options[i] and i < len(flml_options):
                if('<continuous_galerkin replaces="DISOTT">' in flml_options[i]):
                  flml_options[i]=flml_options[i].replace("continuous_galerkin", "legacy_continuous_galerkin")
                i = i + 1
              if("</continuous_galerkin>" in flml_options[i]):
                flml_options[i]=flml_options[i].replace("continuous_galerkin", "legacy_continuous_galerkin")
            if('<continuous_galerkin_test>' in flml_options[i]):
              while not "</continuous_galerkin_test>" in flml_options[i] and i < len(flml_options):
                if('<continuous_galerkin_test>' in flml_options[i]):
                  flml_options[i]=flml_options[i].replace("continuous_galerkin_test", "continuous_galerkin")
                i = i + 1
              if("</continuous_galerkin_test>" in flml_options[i]):
                flml_options[i]=flml_options[i].replace("continuous_galerkin_test", "continuous_galerkin")
            if("<continuous_galerkin_test/>" in flml_options[i]):
              flml_options[i]=flml_options[i].replace("continuous_galerkin_test", "continuous_galerkin")
            if('<mixed_cv_cg' in flml_options[i]):
              while not "</mixed_cv_cg>" in flml_options[i] and i < len(flml_options):
                if('<mixed_cv_cg' in flml_options[i]):
                  flml_options[i]=flml_options[i].replace("mixed_cv_cg", "legacy_mixed_cv_cg")
                i = i + 1
              if("</mixed_cv_cg>" in flml_options[i]):
                flml_options[i]=flml_options[i].replace("mixed_cv_cg", "legacy_mixed_cv_cg")
            i = i + 1
        if('<include_mixing_stats' in flml_options[i]):
          while not "</include_mixing_stats>" in flml_options[i] and i < len(flml_options):
            i = i + 1 # do nothing
        i = i + 1

try:
  flml_file = file(filename, "w")
  flml_file.writelines(flml_options)
  flml_file.close()
except:
  sys.stderr.write("Error: change_options failed to write options file\n")
  sys.exit(1)


