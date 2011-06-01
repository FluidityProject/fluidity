#!/usr/bin/python
import os.path
import numpy
from optparse import OptionParser

# Initial control
if not os.path.isfile("control.npy"):
  m = numpy.array([1.0])
  numpy.save("control.npy", m)


################# Main program ###################
def main():
  parser = OptionParser()
  parser.add_option("-f", "--file", dest="filename", help="the .oml file", metavar="FILE")
  (options, args) = parser.parse_args()
  if not os.path.isfile(options.filename):
    print "File", options.filename, "not found."
    exit()


if '__main__'==__name__:
      main()   
