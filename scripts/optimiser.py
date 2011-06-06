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
import time
import pickle

# Hack for libspud to be able to read an option from a second option file.
# A better solution would be to fix libspud or use an alternative implementation like
# https://github.com/gmarkall/manycore_form_compiler/blob/master/mcfc/optionfile.py
def spud_get_option(filename, option_path):
  d = {}
  exec "import libspud" in d
  exec "libspud.load_options('"+filename+"')" in d
  exec "v = libspud.get_option('"+option_path+"')" in d
  return d['v']

# Executes the model specified in the optimiser option tree
# The model stdout is printed to stdout.
def run_model():
  command_line = libspud.get_option('/model/command_line')
  option_file = libspud.get_option('/model/option_file')
  args = shlex.split(command_line)
  args.append(option_file)
  p = Popen(args, stdout=PIPE,stderr=PIPE)
  out = string.join(p.stdout.readlines() )
  outerr = string.join(p.stderr.readlines() )
  if p.wait() != 0:
    print "Model execution failed: "
    print outerr
    exit()
  print out

def initialise_model_controls():
  update_type = libspud.get_option('/control_io/type[0]/name')
  # With the custom type, the user specifies python function to initialise the controls. 
  if update_type == 'custom':
    get_initial_code = libspud.get_option('/control_io/type::custom/get_initial_controls')
    d = {}
    exec get_initial_code in d
    m = d['get_initial_controls']()
  else:
    print "Unknown type ", libspud.get_option('/control_io/type[0]/name'), " in /control_io/type"
    exit()
  return m

def update_model_controls(m):
  update_type = libspud.get_option('/control_io/type[0]/name')
  # With the custom type, the user specifies a python function to upadate the controls. 
  if update_type == 'custom':
    update_code = libspud.get_option('/control_io/type::custom/update_controls')
    d = {}
    exec update_code in d
    d['update_controls'](m)
  else:
    print "Unknown type ", libspud.get_option('/control_io/type[0]/name'), " in /control_io/type"
    exit()

################# Optimisation loop ###################
def optimisation_loop():
  # This function takes in a dictionary m with numpy.array as entries. 
  # From that it creates one serialised numpy.array with all the data.
  # In addition it creates m_shape, a dictionary which is used in unserialise.
  def serialise(m):
    m_serial = numpy.array([])
    m_shape = {}
    for k, v in m.iteritems():
      m_serial = numpy.append(m_serial, v)
      m_shape[k] = len(m_serial)
    return [m_serial, m_shape]

  # Reconstructs the original dictionary of numpy.array's from the serialised version and the shape.
  def unserialise(m_serial, m_shape):
    m = {}
    current_index = 0
    for k, v in m_shape.iteritems():
      m[k] = m_serial[current_index:v]
      current_index = v
    return m

  # Retuns the functional value with the current controls
  def J(m_serial, args):
    m_shape = args
    m = unserialise(m_serial, m_shape)
    update_model_controls(m)
    run_model()
    functional = libspud.get_option('/functional/name')
    option_file = libspud.get_option('/model/option_file')
    simulation_name = spud_get_option(option_file, "/simulation_name")
    stat_file = simulation_name+'_adjoint_'+functional+".stat"
    J = stat(stat_file)[functional]["value"][-1]
    print "J = ", J
    return J

  # Retuns the functional derivative with respect to the controls.
  def dJdm(m_serial, args):
    m_shape = args
    m = unserialise(m_serial, m_shape)
    update_model_controls(m)
    run_model()
    pkl_file = open('func_derivs.pkl', 'rb')
    djdm = pickle.load(pkl_file)
    # Check that the the controls in dJdm are consistent with the ones specified in m
    djdm_keys = djdm.keys()
    djdm_keys.sort()
    m_keys = djdm.keys()
    m_keys.sort()
    if m_keys != djdm_keys:
      print "The specified controls are not consistent with the controls in the derivative in the objective function."
      print "The specified controls are:", m_keys
      print "The controls in dJdm are:", djdm_keys
      print "Check the consistency of the control definition in the model and the optimiser configuration."
      exit()
    for k, v in m.iteritems():
      if len(m[k]) != len(djdm[k]):
        print "The control ", k, " has dimension ", len(m[k]), " but dJd(", k, ") has dimension ", len(djdm[k])
        exit()
    # Serialise djdm in the same order than m_serial
    djdm_serial = [] 
    for k, v in m_shape.iteritems():
      djdm_serial = numpy.append(djdm_serial, djdm[k])
    return djdm_serial

  # Get the optimisation settings
  algo = libspud.get_option('optimisation_options/optimisation_algorithm[0]/name')
  tol = libspud.get_option('/optimisation_options/tolerance')
  # Initialise the controls
  m = initialise_model_controls()
  [m_serial, m_shape] = serialise(m)
  print "Using ", algo, " as optimisation algorithm."
  if algo == 'BFGS':
      res = scipy.optimize.fmin_bfgs(J, m_serial, dJdm, gtol=tol, full_output=1, args=(m_shape, ))
  if algo == 'NCG':
      res = scipy.optimize.fmin_ncg(J, m_serial, dJdm, avextol=tol, full_output=1, args=(m_shape, ))
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
  if not libspud.have_option('/optimisation_options'):
    print "File", args.filename, "is not a valid .oml file."
    exit()
  if not os.path.isfile(libspud.get_option('/model/option_file')):
      print "Could not find ", libspud.get_option('/model/option_file') ," as specified in /model/option_file"
      exit()
  # Start the optimisation loop
  optimisation_loop() 

if '__main__'==__name__:
      start_time = time.time()
      main()   
      print "Optimisation finished in ", time.time() - start_time, "seconds"
