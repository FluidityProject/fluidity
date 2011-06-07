#!/usr/bin/python
import os.path
import numpy
import argparse
import libspud
import shlex 
from subprocess import Popen, PIPE
import scipy.optimize
import string
from fluidity_tools import stat_parser 
from fluidity_tools import stat_creator
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
def run_model(m):
  update_model_controls(m)
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
  # Implement a memoization function to avoid duplicated functional (derivative) evaluations
  class MemoizeMutable:
    def __init__(self, fn):
      self.fn = fn
      self.memo = {}
    def __call__(self, *args, **kwds):
      import cPickle
      str = cPickle.dumps(args, 1)+cPickle.dumps(kwds, 1)
      if not self.memo.has_key(str): 
        self.memo[str] = self.fn(*args, **kwds)
      return self.memo[str]
    # Insert a function value into the cache manually.
    def __add__(self, value, *args, **kwds):
      import cPickle
      str = cPickle.dumps(args, 1)+cPickle.dumps(kwds, 1)
      self.memo[str] = value
   
 
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
  
  # Returns the functional value with the current controls
  def J(m_serial, m_shape, write_stat=True):
    J = mem_pure_J(m_serial, m_shape)
    print "J = ", J 
    if write_stat:
      # Update the functional value in the optimisation stat file
      functional = libspud.get_option('/functional/name')
      stat_writer[(functional, 'value')] = J
    return J

  # A pure version of the computation of J 
  def pure_J(m_serial, m_shape):
    m = unserialise(m_serial, m_shape)
    run_model(m)
    functional = libspud.get_option('/functional/name')
    option_file = libspud.get_option('/model/option_file')
    simulation_name = spud_get_option(option_file, "/simulation_name")
    stat_file = simulation_name+".stat"
    J = stat_parser(stat_file)[functional]["value"][-1]
    return J

  # Returns the functional derivative with respect to the controls.
  def dJdm(m_serial, m_shape, write_stat=True):
    return mem_pure_dJdm(m_serial, m_shape)

  # A pure version of the computation of J 
  def pure_dJdm(m_serial, m_shape):
    m = unserialise(m_serial, m_shape)
    run_model(m)
    # The dJdm we also run the forward model, and in particular we computed the 
    # functional values. In order not to compute the functional values again when calling 
    # J with, we manually add fill the memoize cache.
    functional = libspud.get_option('/functional/name')
    option_file = libspud.get_option('/model/option_file')
    simulation_name = spud_get_option(option_file, "/simulation_name")
    stat_file = simulation_name+".stat"
    J = stat_parser(stat_file)[functional]["value"][-1]
    # Add the functional value the memJ's cache
    mem_pure_J.__add__(J, m_serial, m_shape)
    # Now get the functional derivative information
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

  # This function gets called after each optimisation iteration. 
  # It is currently used to do write the statistics to the stat file.
  def callback(m_serial, m_shape):
    if libspud.have_option("/debugging/check_gradient"):
      grad_err = scipy.optimize.check_grad(lambda x: J(x, m_shape, write_stat = False), lambda x: dJdm(x, m_shape, write_stat = False), m_serial)
      functional = libspud.get_option('/functional/name')
      stat_writer[(functional + "_gradient_error", "l2norm")] = grad_err
    stat_writer.write()

  # Initialise stat file
  stat_writer=stat_creator(libspud.get_option('/name').strip() + '.stat')
  # Get the optimisation settings
  algo = libspud.get_option('optimisation_options/optimisation_algorithm[0]/name')
  tol = libspud.get_option('/optimisation_options/tolerance')
  # Initialise the controls
  m = initialise_model_controls()
  [m_serial, m_shape] = serialise(m)
  print "Using ", algo, " as optimisation algorithm."
  # Create the memoized version of the functional (derivative) evaluation functions
  mem_pure_dJdm = MemoizeMutable(pure_dJdm)
  mem_pure_J = MemoizeMutable(pure_J)
  if algo == 'BFGS':
    res = scipy.optimize.fmin_bfgs(J, m_serial, dJdm, gtol=tol, full_output=1, args=(m_shape, ), callback = lambda m: callback(m, m_shape))
  if algo == 'NCG':
    res = scipy.optimize.fmin_ncg(J, m_serial, dJdm, avextol=tol, full_output=1, args=(m_shape, ), callback = lambda m: callback(m, m_shape))
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
