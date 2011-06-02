#!/usr/bin/python
import os.path
import numpy
import argparse
import libspud
import shlex, subprocess
import scipy.optimize 

# Initial control
#if not os.path.isfile("control.npy"):
#  m = numpy.array([1.0])
#  numpy.save("control.npy", m)

def optimisation_loop():
  def J(m):
    cmd = libspud.get_option('/model_template/option_file')
    option_template = libspud.get_option('/model_template/command_line')
    args = shlex.split(cmd)
    assert(False)
    subprocess.call()
    u = solve_poisson(n, h, m, q)
    d = numpy.dot(Q, u)
    return get_func(h, m, m_opt, d, dobs, alpha)

  def Jprime(m):
    u = solve_poisson(n, h, m, q)
    return get_Jfuncm(n, h, m, m_opt, q, u, dobs, alpha)


  # Start the optimisation loop
  algo = libspud.get_option('optimisation_options/optimisation_algorithm[0]/name')
  tol = libspud.get_option('/optimisation_options/tolerance')
  print "Using ", algo, " as optimisation algorithm."
  if algo == 'BFGS':
      res = scipy.optimize.fmin_bfgs(J, numpy.ones(n), Jprime, gtol=tol, full_output=1)
  else:
    print "Unknown optimisation algorithm in option path."
    exit()
  
  print "Functional value J(m): ", res[1]
  print "Control state m: ", res[0]
  print "Control state error m - m_opt =  ", res[0] - m_opt



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
  
  optimisation_loop() 

if '__main__'==__name__:
      main()   
