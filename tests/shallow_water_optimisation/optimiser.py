#!/usr/bin/python
import os.path
import numpy
import argparse
import libspud
import shlex 
from subprocess import Popen, PIPE
import scipy.optimize
import string
from fluidity_tools import stat_parser as stat

# Initial control
#if not os.path.isfile("control.npy"):
#  m = numpy.array([1.0])
#  numpy.save("control.npy", m)

def run_model():
  command_line = libspud.get_option('/model_template/command_line')
  option_file = libspud.get_option('/model_template/option_file')
  args = shlex.split(command_line)
  args.append(option_file)
  p = Popen(args, stdout=PIPE,stderr=PIPE)
  out = string.join(p.stdout.readlines() )
  outerr = string.join(p.stderr.readlines() )
  if p.wait() != 0:
    print "Model execution failed: "
    print outerr
    print "See .. for more information."
    exit()
  print out


def optimisation_loop():
  def J(m):
    # TODO: Update m
    run_model()
    J = stat("wave_A_adjoint_integral_eta_t1.stat")["integral_eta_t1"]["value"][-1]
    print "J = ", J
    return J

  def dJdm(m):
    # TODO: Update m
    run_model()
    djdm = numpy.load("djdm.npy")
    print "dJdm = ", djdm
    return djdm

  # Start the optimisation loop
  algo = libspud.get_option('optimisation_options/optimisation_algorithm[0]/name')
  tol = libspud.get_option('/optimisation_options/tolerance')
  m = numpy.ones(1)
  print "Using ", algo, " as optimisation algorithm."
  if algo == 'BFGS':
      res = scipy.optimize.fmin_bfgs(J, m, dJdm, gtol=tol, full_output=1)
  else:
    print "Unknown optimisation algorithm in option path."
    exit()
  
  print "Functional value J(m): ", res[1]
  print "Control state m: ", res[0]

################# Main program ###################
def main():
  parser = argparse.ArgumentParser(description='Optimisation program for fluidity.')
  parser.add_argument('filename', metavar='FILE', help="the .oml file")
  args = parser.parse_args()
  if not os.path.isfile(args.filename):
    print "File", args.filename, "not found."
    exit()
  libspud.load_options(args.filename)
  if not libspud.have_option('/optimisation_name/'):
    print "File", args.filename, "is not a valid .oml file."
    exit()
  if not os.path.isfile(libspud.get_option('/model_template/option_file')):
      print "Could not find ", libspud.get_option('/model_template/option_file') ," as specified in /model_template/option_file"
      exit()
  

  optimisation_loop() 

if '__main__'==__name__:
      main()   
