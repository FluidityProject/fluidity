#!/usr/bin/python
import os.path
import numpy
import argparse
import shlex 
from subprocess import Popen, PIPE
import scipy.optimize
import string
import libspud
from fluidity_tools import stat_parser 
from fluidity_tools import stat_creator
import time
import pickle
import glob

# Hack for libspud to be able to read an option from a different files. 
# A better solution would be to fix libspud or use an alternative implementation like
# https://github.com/gmarkall/manycore_form_compiler/blob/master/mcfc/optionfile.py
def superspud(filename, cmd):
  libspud.load_options(filename)
  r = None
  if hasattr(cmd, '__iter__'):
    for c in cmd:
      exec "try: r = " + c + "\nexcept libspud.SpudNewKeyWarning: pass"
  else:
    exec "try: r = " + cmd + "\nexcept libspud.SpudNewKeyWarning: pass"
  libspud.clear_options()
  return r

# Executes the model specified in the optimiser option tree
# The model stdout is printed to stdout.
def run_model(m, opt_options, model_options):
  update_custom_controls(m, opt_options)
  if (superspud(model_options, "libspud.have_option('/adjoint/controls/load_controls')")):
    # If the model is loading the default controls, we need to make suer the control files are up to date:
    update_default_controls(m, opt_options, model_options)
  command_line = superspud(opt_options, "libspud.get_option('/model/command_line')")
  option_file = superspud(opt_options, "libspud.get_option('/model/option_file')")
  args = shlex.split(command_line)
  args.append(option_file)
  p = Popen(args, stdout=PIPE,stderr=PIPE)
  out = string.join(p.stdout.readlines() )
  outerr = string.join(p.stderr.readlines() )
  if p.wait() != 0:
    print "Model execution failed."
    print "The error was:"
    print outerr
    exit()
  if verbose:
    print "Model output: "
    print out

# Intialises the custom controls using the supplied python code.
def get_custom_controls(opt_options):
  nb_controls = superspud(opt_options, "libspud.option_count('/control_io/control')")
  m = {}
  for i in range(nb_controls):
    cname = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/name')")
    ctype = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/type/name')")
    # With the custom type, the user specifies python function to initialise the controls. 
    if ctype == 'custom':
      initial_control_code = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/type::custom/initial_control')")
      d = {}
      exec initial_control_code in d
      m[cname] = d['initial_control']()
  return m

# Initialse the default controls by reading in the control files/
# This assumes that the model has been run in advance and produced the initial control files.
def read_default_controls(opt_options, model_options):
  simulation_name = superspud(model_options, "libspud.get_option('/simulation_name')")
  nb_controls = superspud(opt_options, "libspud.option_count('/control_io/control')")
  m = {}
  for i in range(nb_controls):
    cname = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/name')")
    ctype = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/type/name')")
    if ctype == 'default':
      act_flag = False # Check that at least one control file exists
      for ctrl_file in glob.iglob('control_'+simulation_name+'_'+cname+ '_[0-9]*.pkl'):
        try:
          timestep = int(ctrl_file.strip()[len('control_'+simulation_name+'_'+ cname+ '_'):len(ctrl_file)-4])
        except:
          print "Error while reading the control files."
          print "The control file ", ctrl_file, " does not conform the standard naming conventions for control files."
          exit()
        f = open(ctrl_file, 'rb')
        m[(cname, timestep)] = pickle.load(f) 
        f.close()
        act_flag = True
      if act_flag == False:
        print "Warning: Found no control derivative file for control ", cname, "."
  return m

# Returns the control derivatives for both the custom and the default controls. 
def read_control_derivatives(opt_options, model_options):
  simulation_name = superspud(model_options, "libspud.get_option('/simulation_name')")
  functional_name = superspud(opt_options, "libspud.get_option('/functional/name')")
  nb_controls = superspud(opt_options, "libspud.option_count('/control_io/control')")
  derivs = {}
  for i in range(nb_controls):
    cname = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/name')")
    ctype = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/type/name')")
    if ctype == 'default':
      act_flag = False # Check that at least one control file exists
      for ctrl_file in glob.iglob('control_'+simulation_name+'_adjoint_'+functional_name+'_'+ cname+ '_TotalDerivative_[0-9]*.pkl'):
        try:
          # The naming convenction is control+simulation_name+control_name+TotalDerivative, but do not forget that
          # the derivatives where produced during the adjoint run in which the simulation name is simulation_name+functional_name
          timestep = int(ctrl_file.strip()[len('control_'+simulation_name+'_adjoint_'+functional_name+'_'+ cname+ '_TotalDerivative_'):len(ctrl_file)-4])
        except:
          print "Error while reading the control derivative files."
          print "The control file ", ctrl_file, " does not conform the standard naming conventions for control files."
          exit()
        f = open(ctrl_file, 'rb')
        derivs[(cname, timestep)] = pickle.load(f) 
        f.close()
        act_flag = True
      if act_flag == False:
        print "Warning: Found no control derivative file for control ", cname, "."
    elif ctype == 'custom':
      control_derivative_code = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/type::custom/control_derivative')")
      d = {}
      exec control_derivative_code in d
      derivs[cname] = d['control_derivative']()
    else:
      print "Unknown control type " + ctype + "."
      exit()
  return derivs

# Writes the custom controls onto disk
def update_custom_controls(m, opt_options):
  nb_controls = superspud(opt_options, "libspud.option_count('/control_io/control')")
  for i in range(nb_controls):
    cname = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/name')")
    ctype = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/type/name')")
    # With the custom type, the user specifies a python function to update the controls. 
    if ctype == 'custom':
      update_control_code = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/type::custom/update_control')")
      d = {}
      exec update_control_code in d
      d['update_control'](m[cname])

# Writes the default controls onto disk
def update_default_controls(m, opt_options, model_options):
  global debug
  simulation_name = superspud(model_options, "libspud.get_option('/simulation_name')")
  nb_controls = superspud(opt_options, "libspud.option_count('/control_io/control')")
  # Loop over default controls
  for i in range(nb_controls):
    cname = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/name')")
    ctype = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/type/name')")
    if ctype == 'default':
      # Loop over controls
      for k in m.keys():
        # Check if that is a control we are looking for
        if k[0] == cname:
          timestep = k[1]
          file_name = 'control_'+simulation_name + '_' + cname + '_' + str(timestep) + '.pkl'
          if not os.path.isfile(file_name):
            print "Error: writing control file ", file_name, " which did not exist before."
            exit()
          if debug:
            # Check that the file we are writing has the same shape than the one we are writing
            f = open(file_name, 'rb')
            m_old = pickle.load(f)
            if m[k].shape != m_old.shape:
              print "Error: The shape of the control in ", file_name, " changed."
              exit()
            f.close()
          f = open(file_name, 'wb')  
          pickle.dump(m[k], f)
          f.close()

# Check the consistency of model and option file
def check_option_consistency(opt_options, model_options):
  nb_controls = superspud(opt_options, "libspud.option_count('/control_io/control')")
  for i in range(nb_controls):
    cname = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/name')")
    ctype = superspud(opt_options, "libspud.get_option('/control_io/control["+str(i)+"]/type/name')")
    # Check that the default controls exist in the model
    # and that custom controls not.
    if ctype == 'custom':
      if superspud(model_options, "libspud.have_option('/adjoint/controls/control::" + cname + "')"):
        print "The custom control " + cname + " is a default control in the model option tree."
        exit()
    elif ctype== 'default':
      if not superspud(model_options, "libspud.have_option('/adjoint/controls/control::" + cname + "')"):
        print "The default control " + cname + " was not found in the model option tree."
        exit()
    else:
      print "Unknown control type " + ctype + "."
      exit()

# Check that the the controls in dJdm are consistent with the ones in m
def check_control_consistency(m, djdm):
    djdm_keys = djdm.keys()
    m_keys = m.keys()
    djdm_keys.sort()
    m_keys.sort()
    if m_keys != djdm_keys:
      print "Error: The controls are not consistent with the controls derivatives."
      print "The controls are:", m_keys
      print "The control derivatives are:", djdm_keys
      print "Check the consistency of the control definition in the model and the optimiser configuration."
      exit()
    for k, v in m.iteritems():
      if m[k].shape != djdm[k].shape:
        assert(False)
        print "The control ", k, " has shape ", m[k].shape, " but dJd(", k, ") has shape ", djdm[k].shape
        exit()

def delete_temporary_files(model_options):
  # remove any control files
  pkl_files = glob.glob('control_*.pkl')
  for f in pkl_files:
    os.remove(f)   
  # remove any stat files from the model
  simulation_name = superspud(model_options, "libspud.get_option('/simulation_name')")
  stat_files = glob.glob(simulation_name+'*.stat')
  for f in stat_files:
    os.remove(f)   

################# Optimisation loop ###################
def optimisation_loop(opt_options, model_options):
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
      else:
        if verbose:
          print "Cache hit for functional evaluation (", self.fn, ")."
      return self.memo[str]
    # Insert a function value into the cache manually.
    def __add__(self, value, *args, **kwds):
      import cPickle
      str = cPickle.dumps(args, 1)+cPickle.dumps(kwds, 1)
      self.memo[str] = value

  # Small test code for the un/serialiser
  def test_serialise():
    x = {'a': numpy.random.rand(3,2), 'b':  numpy.random.rand(3,2,4,5), 'c': numpy.random.rand(1)}
    [m_serial, m_shape] = serialise(x)
    x_re = unserialise(m_serial, m_shape)
    return (x['a'] == x_re['a']).all() and (x['b'] == x_re['b']).all() and (x['c'] == x_re['c']).all()

  # This function takes in a dictionary m with numpy.array as entries. 
  # From that it creates one serialised numpy.array with all the data.
  # In addition it creates m_shape, a dictionary which is used in unserialise.
  def serialise(m):
    m_serial = numpy.array([])
    m_shape = {}
    for k, v in m.iteritems():
      m_serial = numpy.append(m_serial, v.flatten())
      m_shape[k] = v.shape
    return [m_serial, m_shape]

  # Reconstructs the original dictionary of numpy.array's from the serialised version and the shape.
  def unserialise(m_serial, m_shape):
    m = {}
    start_index = 0
    for k, s in m_shape.iteritems():
      offset = 1
      for d in s:
        offset = offset * d
      end_index = start_index + offset
      m[k] = numpy.reshape(m_serial[start_index:end_index], s)
      start_index = end_index
    return m
  
  # Returns the functional value with the current controls
  def J(m_serial, m_shape, write_stat=True):
    J = mem_pure_J(m_serial, m_shape)
    print "J = ", J 
    if write_stat:
      # Update the functional value in the optimisation stat file
      stat_writer[(functional_name, 'value')] = J
    return J

  # A pure version of the computation of J 
  def pure_J(m_serial, m_shape):
    if verbose:
      print "Running forward model for functional evaluation (<function pure_J>)"
    m = unserialise(m_serial, m_shape)
    run_model(m, opt_options, model_options)
    simulation_name = superspud(model_options, "libspud.get_option('/simulation_name')")
    stat_file = simulation_name+".stat"
    J = stat_parser(stat_file)[functional_name]["value"][-1]
    return J

  # Returns the functional derivative with respect to the controls.
  def dJdm(m_serial, m_shape, write_stat=True):
    return mem_pure_dJdm(m_serial, m_shape)

  # A pure version of the computation of J 
  def pure_dJdm(m_serial, m_shape):
    if verbose:
      print "Running forward/adjoint model for functional derivative evaluation (<function pure_dJdm>)"
    m = unserialise(m_serial, m_shape)
    run_model(m, opt_options, model_options)
    # While computing dJdm we run the forward/adjoint model and in particular we compute the 
    # functional values. In order to not compute the functional values again when calling 
    # J, we manually add write it into the memoize cache.
    simulation_name = superspud(model_options, "libspud.get_option('/simulation_name')")
    stat_file = simulation_name+".stat"
    J = stat_parser(stat_file)[functional_name]["value"][-1]
    # Add the functional value the memJ's cache
    mem_pure_J.__add__(J, m_serial, m_shape)
    # Now get the functional derivative information
    djdm = read_control_derivatives(opt_options, model_options)
    check_control_consistency(m, djdm)
    # Serialise djdm in the same order than m_serial
    djdm_serial = [] 
    for k, v in m_shape.iteritems():
      djdm_serial = numpy.append(djdm_serial, djdm[k])
    return djdm_serial

  # This function gets called after each optimisation iteration. 
  # It is currently used to do write the statistics to the stat file.
  def callback(m_serial, m_shape):
    if superspud(opt_options, "libspud.have_option('/debugging/check_gradient')"):
      grad_err = scipy.optimize.check_grad(lambda x: J(x, m_shape, write_stat = False), lambda x: dJdm(x, m_shape, write_stat = False), m_serial)
      stat_writer[(functional_name + "_gradient_error", "l2norm")] = grad_err
    stat_writer.write()
  
  # Initialise stat file
  if verbose:
    print "Initialise stat file"
  stat_writer=stat_creator(superspud(opt_options, "libspud.get_option('/name')").strip() + '.stat')
  # Get the optimisation settings
  if verbose:
    print "Read oml settings"
  algo = superspud(opt_options, "libspud.get_option('optimisation_options/optimisation_algorithm[0]/name')")
  tol = superspud(opt_options, "libspud.get_option('/optimisation_options/tolerance')")
  # Create the memoized version of the functional (derivative) evaluation functions
  mem_pure_dJdm = MemoizeMutable(pure_dJdm)
  mem_pure_J = MemoizeMutable(pure_J)
  # Initialise the controls
  # First we initialise the custom controls
  # This has to be done first since the next step
  # involves running the model and therefore 
  # will need the custom controls to be set.
  if verbose:
    print "Get initial custom controls"
  custom_m = get_custom_controls(opt_options)
  # To get the initial default controls we run the model without the option
  # /adjoint/controls/load_controls
  if verbose:
    print "Get initial default controls"
  model_file = superspud(opt_options, "libspud.get_option('/model/option_file')")
  if (superspud(model_options, "libspud.have_option('/adjoint/controls/load_controls')")):
    superspud(model_options, ["libspud.delete_option('/adjoint/controls/load_controls')", "libspud.write_options('"+ model_file +"')"])
  # Run the forward model including adjoint.
  functional_name = superspud(opt_options, "libspud.get_option('/functional/name')")
  if superspud(opt_options, "libspud.have_option('/adjoint/functional::"+functional_name+"/disable_adjoint_run')"):
    superspud(opt_options, "libspud.delete_option('/adjoint/functional::"+functional_name+"/disable_adjoint_run')")
  [custom_m_serial, custom_m_shape] = serialise(custom_m)
  mem_pure_J(custom_m_serial, custom_m_shape)
  # This should have created all the default initial controls and we can now activate the load_controls flag. 
  superspud(model_options, ["libspud.add_option('/adjoint/controls/load_controls')", "libspud.write_options('"+ model_file +"')"])
  # Load the default controls
  m = read_default_controls(opt_options, model_options)
  nb_controls = len(m) + len(custom_m)
  # And merge them
  m.update(custom_m)
  if (nb_controls != len(m)):
    print "Error: Two controls with the same name defined."
    print "The controls must have all unique names."
    print "Your controls are: ", m.keys()
    exit()
  # Since now all the controls and derivatives are defined, we can check the consistency of the control variables
  djdm = read_control_derivatives(opt_options, model_options)
  check_control_consistency(m, djdm)

  [m_serial, m_shape] = serialise(m)
  if verbose:
    print "Start optimisation loop"
    print "Using ", algo, " as optimisation algorithm."
  if algo == 'BFGS':
    res = scipy.optimize.fmin_bfgs(J, m_serial, dJdm, gtol=tol, full_output=1, args=(m_shape, ), callback = lambda m: callback(m, m_shape))
  elif algo == 'NCG':
    res = scipy.optimize.fmin_ncg(J, m_serial, dJdm, avextol=tol, full_output=1, args=(m_shape, ), callback = lambda m: callback(m, m_shape))
  else:
    print "Unknown optimisation algorithm in option path."
    exit()
  
  if verbose:
    print "End of optimisation loop"
  print "Functional value J(m): ", res[1]
  print "Control state m: ", res[0]

################# Main program ###################
def main():
  global verbose
  global debug
  parser = argparse.ArgumentParser(description='Optimisation program for fluidity.')
  parser.add_argument('filename', metavar='FILE', help="the .oml file")
  parser.add_argument('-v', '--verbose', action='store_true', help='verbose mode')
  parser.add_argument('-d', '--debug', action='store_true', help='the debug mode runs additional internal tests.')
  args = parser.parse_args()
  verbose = args.verbose
  debug = args.debug
  if not os.path.isfile(args.filename):
    print "File", args.filename, "not found."
    exit()
  # Initial spud environments for the optimiser and model options.
  opt_options = args.filename
  if not superspud(opt_options, "libspud.have_option('/optimisation_options')"):
    print "File", args.filename, "is not a valid .oml file."
    exit()
  model_file = superspud(opt_options, "libspud.get_option('/model/option_file')")
  if not os.path.isfile(model_file):
      print "Could not find ", model_file ," as specified in /model/option_file"
      exit()
  # Check consistency of the option files
  check_option_consistency(opt_options, model_file)
  # Start the optimisation loop
  optimisation_loop(opt_options, model_file) 


if '__main__'==__name__:
      start_time = time.time()
      main()   
      print "Optimisation finished in ", time.time() - start_time, "seconds"
