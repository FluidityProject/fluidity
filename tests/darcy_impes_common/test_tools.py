#!/usr/bin/env python

"""Extension of batch_tools.py for test simulations.

More documentation to come."""

import os
import numpy
import sys
from subprocess import call
from string import Template
from re import sub
from collections import OrderedDict
from importlib import import_module
from getpass import getuser
from batch_tools import default_fluidity_path, Tracker, \
    SelfHandlingCommand, DoNothing
from copy import deepcopy
from warnings import warn

default_fluidity_path()
from fluidity_tools import stat_parser


## HELPERS

# since it is likely that the client will want to give all the classes
# in this module the same verbosity, verbosity may as well be a module
# variable.  The verbosity() function below can be used to set/get
# the variable.  It is an integer to allow different levels of verbosity.
# TODO: consider replacing with a singleton
verbosity_level = 0

# output filenames
error_rates_filename = "error_rates.log"
error_norms_filename = "error_norms.log"

xml_skeleton_begin =  """<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE testproblem SYSTEM "regressiontest.dtd">

<testproblem>
  <name>{0}</name>
  <owner userid="{1}"/>
  <tags>{2}</tags>
  <problem_definition length="short" nprocs="1">
    <command_line>
{3}
    </command_line>
  </problem_definition>
  <variables/>
  <pass_tests>
    <test name="Solvers converged" language="python">
import os
files = os.listdir("./")
assert(not "matrixdump" in files and not "matrixdump.info" in files)
        </test>"""

xml_skeleton_end = """
  </pass_tests>
  <warn_tests>
  </warn_tests>
</testproblem>"""

xml_skeleton_assert = """
    <test name="{1}: expect {0} {2} {3:g}" language="python">
from test_tools import find_{0}
assert(find_{0}("{1}") &{2}; {3:g})
    </test>"""

            

## FOR PUBLIC USE

def verbosity(new_verbosity_level=None):
    """Setter/getter of global verbosity_level."""
    global verbosity_level
    if new_verbosity_level is not None:
        verbosity_level = new_verbosity_level
    else:
        return verbosity_level
   
def find_in_open(report, key):
    """Searches a report and returns the value associated with key. """
    report.seek(0)
    v = numpy.nan
    for l in report:
        if not l.strip(): continue
        if l.split()[0] == key: 
            v = float(l.split()[1])
            break
    return v
    
def find(report_filename, key):
    """Opens and searches report_filename and returns the value associated
    with key (ID).
    """
    with open(report_filename, "r") as f:
        v = find_in_open(f, key)
        return v

def find_norm(key): 
    return find(error_norms_filename, key)

def find_rate(key):
    return find(error_rates_filename, key)

            
class DefaultParameteriser:
   """A class for determining filenames and making dictionaries for
   template file expansion.  Clients can extend the class to customise
   any of the methods.  The parameteriser generally needs to be kept up
   to date, via the update() method, before methods are called."""
   # TODO perhaps update() isn't needed, since the underlying component
   # self.tracker will be a reference and not a copy.  

   def __init__(self, options_filename_extension, domain_length_1D=1.):
       self.options_filename_extension_str = options_filename_extension
       self.domain_length_1D_value = domain_length_1D
       # see update_dict() for an explanation of how the following
       # variables are used
       self.master_dict = {}
       self.num_expansions = {}

   def update(self, tracker):
       self.tracker = tracker

   def get(self, what):
        """If what=='level', returns current level name.  If what=='node',
        returns current node name.  Otherwise assumes the argument is a
        level name from which a node name is to be returned.  Return
        value defaults to None."""
        # delegate
        return self.tracker.get(what)

   def at(self, what):
        """Returns True if where=='root' and we are at the base of the
        handler tree, if where=='leaf' and we are at the end of a
        branch, or if where is the current level_name.
        """
        # delegate
        return self.tracker.at(where)

   def make_key(self, level_names=None, exclude=[], substitute={}):
        """Forms a string made up of node names corresponding to the
        supplied list of level names.
        """ 
        # delegate
        return self.tracker.make_key(level_names, exclude, substitute)

   def update_dict(self, dict_name, user_defined_dict={}, nesting_level=1):
       """Registers or updates a named dictionary for use by lookup()
       and expand_template().  Note that this may get merged with a
       corresponding built-in dict (e.g. 'geometry', 'options') and,
       even if the client does not specify a dictionary argument, the
       built-in will still get registered."""
       # in special cases, start with a built-in dict
       if dict_name=='geometry':
           result_dict = self.geometry_dict()
       elif dict_name=='options':
           result_dict = self.options_dict()
       else:
           result_dict = {}
       # update/merge with the optional dictionary in the argument list
       result_dict.update(user_defined_dict)
       # store the result as a sub-dictionary in master, and record the
       # total number of expansions to be performed
       if dict_name not in self.master_dict.keys():
           # here we create a fresh dictionary
           self.master_dict[dict_name] = result_dict
       else:
           # here we update an existing sub-dictionary
           self.master_dict[dict_name].update(result_dict)
           
       
   def lookup(self, dict_name, key):
       """Looks up a variable in the dictionary corresponding to
       dict_name."""
       # update the embedded dict, for safety
       self.update_dict(dict_name)
       sub_dict = self.master_dict[dict_name]
       return sub_dict[key]

   def expand_template(self, buffer, dict_name):
       """In a given buffer, replaces all instances of strings prefixed
       '$' with the corresponding values in the dictionary corresponding
       to dict_name."""
       max_num_loops = 5
       # update the embedded dict, for safety
       self.update_dict(dict_name)
       # substitute until there are no more placeholders
       n = 0
       while '$' in buffer and n < max_num_loops:
           buffer = Template(buffer)
           buffer = buffer.safe_substitute(self.master_dict[dict_name])
           n += 1
       if n == max_num_loops:
           warn("Maximum number of substitutions ({0}) reached".format(n))
       return buffer

   def geometry_dict(self):
      """Built-in dictionary for expanding a geometry template.
      Currently minimal, with meshing support for a 1D domain only.  The
      client may want to override this method to define new entries"""
      L = float(self.domain_length_x)
      n = int(self.get('mesh_res'))
      return {
          'MESH_NAME': self.mesh_name(),
          'DOMAIN_LENGTH_X': str(L),
          'EL_NUM_X': str(n),
          'EL_SIZE_X': str(L/n) }
   
   def options_dict(self):
      """Built-in dictionary for expanding an options template.
      Currently minimal.  The client may want to override this method to
      define new entries.
      """
      dim = self.get('dim')
      if dim=="1":
         mesh_format = 'triangle'
      else:
         mesh_format = 'gmsh'
      return {
          'MESH_NAME': self.mesh_name(),
          'MESH_FORMAT': mesh_format,
          'DOMAIN_DIM': dim }

   def options_filename_extension(self):
       return self.options_filename_extension_str
   
   def options_template_filename(self):
      # assume just one case name
      return "template_" + self.get('case') + self.options_filename_extension()
   
   def options_name(self):
      return self.make_key()
   
   def options_filename(self):
      return self.options_name() + self.options_filename_extension()

   def geometry_filename_extension(self):
      if self.get('dim')=='1':
         return '.sh'
      else:
         return '.geo'

   def mesh_filename_extension_list(self):
      if self.get('dim')=='1':
         mesh_ext_list ['.bound', '.ele', '.node']
      else:
         mesh_ext_list ['.msh']

   def geometry_template_filename(self):
      return 'template_' + self.make_key(['dim', 'mesh_type']) + \
          self.geometry_filename_extension()
   
   def geometry_name(self):
      return self.make_key(['dim', 'mesh_type', 'mesh_res'])
   
   def geometry_filename(self):
      return self.geometry_name() + self.geometry_filename_extension()
   
   def mesh_name(self):
      return self.geometry_name()

   def mesh_filenames(self):
      result = []
      geometry_name = self.geometry_name()
      if self.get('dim')=='1':
         mesh_ext_list = ['.bound', '.ele', '.node']
      else:
         mesh_ext_list = ['.msh']
      for ext in mesh_ext_list:
         result.append(self.geometry_name() + ext)
      return result

      
class TestCommand(SelfHandlingCommand):
    """This abstract class, to be extended by concrete Commands, has
    some test-related implementation and also takes on the role of
    storing level details from which to make labels and filenames."""

    def __init__(self, requisite_level_names,
                 handler, message, parameteriser):
        self.verbosity = verbosity()
        SelfHandlingCommand.__init__(self, requisite_level_names,
                                     handler, message, self.verbosity)
        self.tracker = Tracker({'dim': '{0}d',
                                'norm': 'l{0}',
                                'mesh_res': 'm{0}'})
        self.parameteriser = parameteriser
            
    def inform(self, level_name, node_name, at_leaf, indent, verbose):
        """Records level details, including subtle testing-related details."""

        if level_name == 'mesh_res':
            # mesh_res is a special level; we need to update mesh_res_prev
            res = self.tracker.get('mesh_res')
            if res is not None:
                self.tracker.update('mesh_res_prev', res)

        # in any case record level-node details 
        self.tracker.inform(level_name, node_name, at_leaf, indent)
        self.indent = indent
        self.verbosity = self.verbosity * verbose

        # if this is the top level, infer a name for the test case
        if level_name == self.level_names[0]:
            self.case_name = self.tracker.get(level_name)

        # field is also special; it needs to be reformatted/split into
        # useful bits
        if level_name == 'field':
            pattern = 'Phase(.)::(.*)'
            ph = sub(pattern, '\\1', node_name)
            var = sub(pattern, '\\2', node_name)
            # reformat the field name
            self.tracker.update('field', str.lower(var) + ph)
            # have a separate dictionary so as not to pollute the main one
            self.field_dict = {'full':  node_name,
                               'phase': 'Phase' + ph,
                               'var':   var}
           
    def get(self, what):
        """If what=='level', returns current level name.  If what=='node',
        returns current node name.  Otherwise assumes the argument is a
        level name from which a node name is to be returned.  Return
        value defaults to None."""
        try:
            # try the special field_dict first
            return self.field_dict[what]
        except:
            # delegate to the tracker
            return self.tracker.get(what)

    def at(self, where):
        """Returns True if where=='root' and we are at the base of the
        handler tree, if where=='leaf' and we are at the end of a
        branch, or if where is the current level_name.
        """
        # delegate
        return self.tracker.at(where)

    def make_key(self, level_names=None, exclude=[], substitute={}):
        """Forms a string made up of node names corresponding to the
        supplied list of level names.
        """ 
        # delegate
        return self.tracker.make_key(level_names, exclude, substitute)

    def update_dict(self, dict_name, dict_value):
        """Registers or updates a named dictionary for use by lookup()
        and expand_template().  Note that this may get merged with a
        corresponding built-in dict (e.g. 'geometry', 'options') and,
        even if the client does not specify a dictionary argument, the
        built-in will still get registered."""
        # delegate
        self.parameteriser.update(self.tracker)
        self.parameteriser.update_dict(dict_name, dict_value)

    def lookup(self, dict_name, key):
        """Looks up a key in the dictionary stored under dict_name."""
        # delegate
        self.parameteriser.update(self.tracker)
        return self.parameteriser.lookup(dict_name, key)
   

class WriteXMLFile(TestCommand):
    """Writes out an XML file which will be picked up by Fluidity's test
    harness."""
    
    def __init__(self, handler, parameteriser, command_line,
                 threshold_calculator=None):
        TestCommand.__init__(self, ['field', 'norm', 'mesh_res'], handler, 
                             "Generating XML file", parameteriser)
        self.command_line = command_line
        if threshold_calculator is None:
            self.threshold_calculator = DefaultThresholdCalculator()
        else:
            self.threshold_calculator = threshold_calculator
            
    def __del__(self):
        """Finalises any outstanding open file before expiring."""
        self.write_xml_end()
            
    def execute(self):
        # initialise file if it hasn't been already
        if self.at('root'):
            self.write_xml_begin()
        
        # the last level (at the leaf node) does the important stuff
        if not self.at('leaf'):
            self.mesh_res_prev = None
            return
            
        # write norm assert
        norm_threshold = self.threshold_calculator.\
                         calculate_norm_threshold(self)
        if norm_threshold is not None:
            key = self.tracker.make_key(exclude=['case'])
            self.write_xml_assert('norm', key, 'lt', norm_threshold)
            self.write_message('  norm < {0}'.format(norm_threshold))
                
        # write rate assert
        rate_threshold = self.threshold_calculator.\
                         calculate_rate_threshold(self)
        if rate_threshold is not None:
            mesh_res = self.get('mesh_res')
            if self.mesh_res_prev is None:
                self.mesh_res_prev = mesh_res
                return
            subs_dict = {'mesh_res': '{0}_{1}'.format(
                self.mesh_res_prev, mesh_res)}
            key = self.tracker.make_key(exclude=['case'],
                                        substitute=subs_dict)
            self.write_xml_assert('rate', key, 'gt', rate_threshold)
            self.write_message('  rate > {0}'.format(rate_threshold))
            self.mesh_res_prev = mesh_res

    def write_xml_begin(self):
        self.xml_file = open(self.case_name+'.xml', 'w')
        self.xml_file.write(xml_skeleton_begin.format(
                self.case_name, getuser(),
                self.parameteriser.options_filename_extension(),
                self.command_line))
        
    def write_xml_end(self):
        try:
            self.xml_file.write(xml_skeleton_end)
            self.xml_file.close()
        except AttributeError:
            pass

    def write_xml_assert(self, metric_type, key, rel_op, threshold):
        self.xml_file.write(xml_skeleton_assert.format(
            metric_type, key, rel_op, threshold))

    # def get_norm_threshold(self):
    #     """The client may want to override this if a blanket value isn't
    #     desired.  A rescaling of the norm will be attempted, but the
    #     requirements are very strict: the client needs to have
    #     registered (updated) a dictionary called 'solution', and the
    #     dictionary needs to have a parameter called <field>_scale where
    #     <field> is in the format <variable_name><phase_number>.
    #     """
    #     # TODO: it's too strict.  Refactor this method as strategy
    #     try:
    #         var_name = str.lower(self.get('field')) + '_scale'
    #         # the above variable defines a scale factor
    #         return self.norm_threshold * self.lookup('solution', var_name)
    #     except:
    #         self.norm_threshold

    # def get_rate_threshold(self):
    #     """The client may want to override this if a blanket value isn't
    #     desired"""
    #     # TODO: also refactor this method as strategy
    #     return self.rate_threshold
        
        
class ThresholdCalculator:
    """Strategy to help clients define thresholds in assert statements.
    """
    def calculate_norm_threshold(self, test_command):
        pass
    def calculate_rate_threshold(self, test_command):
        pass
        
class DefaultThresholdCalculator(ThresholdCalculator):
    def __init__(self, norm_threshold=None, rate_threshold=0.8):
        self.norm_threshold = norm_threshold
        self.rate_threshold = rate_threshold
    
    def calculate_norm_threshold(self, test_command):
        """A rescaling of the norm will be attempted, but the
        requirements are very strict: the client needs to have
        registered (updated) a dictionary called 'solution', and the
        dictionary needs to have a parameter called <field>_scale where
        <field> is in the format <variable_name><phase_number>.
        """
        try:
            var_name = str.lower(test_command.get('field')) + '_scale'
            # the above variable defines a scale factor
            return self.norm_threshold * test_command.lookup('solution', var_name)
        except:
            self.norm_threshold

    def calculate_rate_threshold(self, test_command):
        """The client may want to override this if a blanket value isn't
        desired"""
        return self.rate_threshold

        
    
class ExpandOptionsTemplate(TestCommand):
    def __init__(self, handler, parameteriser):
        TestCommand.__init__(self, ['mesh_res'], handler,
                             "Expanding options template", parameteriser)
      
    def execute(self):
        # the last level (at the leaf node) does the important stuff
        if not self.at('leaf'):
            self.mesh_res_prev = None
            return
        self.write_message('...')

        with open(self.parameteriser.options_template_filename()) as src:
            buf = str(src.read())
        buf = self.parameteriser.expand_template(buf, 'options')
        with open(self.parameteriser.options_filename(), 'w') as tgt:
            tgt.write(buf)
        self.write_message('done')
      
            
class ProcessMesh(TestCommand):
    def __init__(self, handler, parameteriser):
        TestCommand.__init__(self, ['mesh_res'], handler, "Meshing",
                             parameteriser)
        # locate the 1D meshing binary
        default_fluidity_path()
        self.interval_path = os.environ["FLUIDITYPATH"] + "bin/interval"
        if not os.path.isfile(self.interval_path): 
            raise RuntimeError("Cannot find " + self.interval_path)

    def execute(self):
        # the last level (at the leaf node) does the important stuff
        if not self.at('leaf'):
            self.mesh_res_prev = None
            return
        self.write_message('...')
            
        # expand mesh template
        self.parameteriser.update(self.tracker)
        with open(self.parameteriser.geometry_template_filename()) as src:
            buf = str(src.read())
        buf = self.parameteriser.expand_template(buf, 'geometry')
        with open(self.parameteriser.geometry_filename(), 'w') as tgt:
            tgt.write(buf)

        # call interval or gmsh
        # TODO: tell user to 'make geo' if needed
        if self.get('dim') == '1':
            lx_str = self.parameteriser.lookup('geometry',
                                               'DOMAIN_LENGTH_X')
            dx_str = self.parameteriser.lookup('geometry',
                                               'EL_SIZE_X')
            call([self.interval_path, '0.0',
                  lx_str, '--dx='+dx_str,
                  self.parameteriser.geometry_name()])
        else:
            call(['gmsh', '-'+self.get('dim'),
                  self.parameteriser.geometry_filename()],
                 stdout=open(os.devnull, 'wb'))
        self.write_message('done')


        
        
class RunSimulation(TestCommand):
    """Runs the simulator."""

    def __init__(self, handler, parameteriser, simulator_name,
                 simulator_verbosity=0):
        TestCommand.__init__(self, ['mesh_res'], handler, \
                                 "Running simulations", parameteriser)
        # locate the simulator binary
        default_fluidity_path()
        self.simulator_path = os.environ["FLUIDITYPATH"] + "bin/" + \
            simulator_name
        if not os.path.isfile(self.simulator_path): 
            raise RuntimeError("Cannot find {0}".format(self.simulator_path))
        self.simulator_verbosity = simulator_verbosity

    def execute(self):
        # exit if we are not at the last level
        if not self.at('leaf'):
            return
        self.write_message('...')

        # start simulation (TODO: guard against absent mesh)
        self.parameteriser.update(self.tracker)
        if self.simulator_verbosity > 0:
            call([self.simulator_path,
                  '-v{0}'.format(simulator_verbosity),
                  '-l {0}.log'.format(self.parameteriser.options_name()),
                  self.parameteriser.options_filename()],
                 stdout=open(os.devnull, 'wb'))
        else:
            call([self.simulator_path,
                  self.parameteriser.options_filename()],
                 stdout=open(os.devnull, 'wb'))
        self.write_message('done')
        
                

class WriteToReport(TestCommand):
    """Writes error norms and/or convergence rates to files."""

    def __init__(self, handler, parameteriser, norm_calculator=None):
        TestCommand.__init__(self, ['field', 'norm', 'mesh_res'], handler, 
                             "Writing report(s)", parameteriser)
        if norm_calculator is None:
            self.norm_calculator = DefaultNormCalculator()
        else:
            self.norm_calculator = norm_calculator

    def __del__(self):
        """Finalises any outstanding open files before expiring."""
        try:
            self.norms_file.close()
            self.rates_file.close()
        except AttributeError:
            pass
        
    def execute(self):
        if self.at('root'):
            # if this is the very first level, create the file objects
            self.norms_file = open(error_norms_filename, 'w')
            self.rates_file = open(error_rates_filename, 'w')

        if not self.at('mesh_res'):
            # for all levels except the resolution level, write out
            # level details
            level_str = '{0}: {1}'.format(self.get('level'),
                                          self.get('node'))
            self.write_message(level_str, self.norms_file, with_newline=True)
            self.write_message(level_str, self.rates_file, with_newline=True)

            # reset these variables which will be computed at the second
            # mesh resolution
            self.mesh_res_prev = None
            self.err_prev = None

        else:
            # for the resolution level, write out keys and values for
            # the norms/rates
            err = self.norm_calculator.calculate_norm(self)
            key = self.tracker.make_key(exclude=['case'])
            # write to file
            self.write_message('{0}{1:12.3e}'.format(key, err),
                               self.norms_file, with_newline=True)
            # write to screen
            self.write_message('   err: {0:.3e}'.format(err))

            # can only start computing rates from the 2nd mesh
            mesh_res = self.get('mesh_res')
            if self.mesh_res_prev is None:
                self.mesh_res_prev = mesh_res
                self.err_prev = err
                return

            subs_dict = {'mesh_res': '{0}_{1}'.format(
                self.mesh_res_prev, mesh_res)}
            key = self.tracker.make_key(exclude=['case'],
                                        substitute=subs_dict)
                                        
            rate = numpy.log2(self.err_prev/err)
            # write to file
            self.write_message('{0}{1:12.6f}'.format(key, rate),
                               self.rates_file, with_newline=True)
            # write to screen
            self.write_message('   rate: {0:.6f}'.format(rate))
                
            self.mesh_res_prev = mesh_res
            self.err_prev = err
            

class NormCalculator:
    """Strategy to help clients define norms their own way."""
    def calculate_norm(self, test_command):
        pass

class DefaultNormCalculator:
    def calculate_norm(self, test_command):
        """if the client intends to examine errors of variable X, he or
        she needs to represent the errors in a diagnostic field called
        XAbsError."""
        filename = test_command.make_key(exclude=['field','norm']) + '.stat'
        if test_command.get('norm')=='1':
            calc_type = "integral"
        elif test_command.get('norm')=='2':
            calc_type = "l2norm"
        return stat_parser(filename)\
            [test_command.get('phase')]\
            [test_command.get('var')+'AbsError']\
            [calc_type][-1]
    
        
            
class CleanUp(TestCommand):
    """Removes generated geometry, mesh, and options files, but leaves
    results intact.  This ensures templates do not get deleted
    (compared with rm -f *.geo, etc)
    """

    def __init__(self, handler, parameteriser):
        TestCommand.__init__(self, [], handler, "Cleaning up",
                             parameteriser)
        
    def execute(self):
        # exit if we are not at the last level
        if not self.at('leaf'):
            return

        # try removing input files, mesh files, results files
        self.parameteriser.update(self.tracker)
        for f in [self.parameteriser.options_filename(),
                  self.parameteriser.geometry_filename()] + \
                  self.parameteriser.mesh_filenames():
                try:
                    os.remove(f)
                    self.write_message('   removed {0}'.format(f),
                                       with_newline=True)
                except OSError:
                    pass


class TestCommandList(TestCommand):
    """Chains multiple Commands together as a single TestCommand.  The
    TestCommandList's handler overrides any handlers associated with the
    contained Commands."""
    
    def __init__(self, handler, message, parameteriser, commands):
        TestCommand.__init__(self, [], handler, message, parameteriser)
        self.commands = commands
        
    def inform(self, level_name, node_name, at_leaf, indent, verbose):
        """Delegates to wrapped commands in addition to calling the
        original superclass method."""
        TestCommand.inform(self, level_name, node_name, at_leaf, indent,
                           verbose)
        for cmd in self.commands:
            cmd.inform(level_name, node_name, at_leaf, indent, False)
            
    def execute(self):
        """Delegates to wrapped commands."""
        # when being verbose, assume all the action takes place on the
        # leaf nodes
        if self.tracker.at('leaf'):
            self.write_message('...')
        for cmd in self.commands:
            cmd.execute()
        if self.tracker.at('leaf'):
            self.write_message('done')

    def update_dict(self, dict_name, dict_value):
        """Delegates to wrapped commands in addition to calling the
        original superclass method."""
        TestCommand.update_dict(self, dict_name, dict_value)
        for cmd in self.commands:
            cmd.update_dict(dict_name, dict_value)


# class TestSuite:
#     """Template for running various TestCommands."""

#     def do(self, what):
#         """General front end method"""
#         cmd_str = str.lower(what[0:3])
#         # decide which commands to run
#         if cmd_str == 'all':
#             cmd_list = [
#                 self.create_generate_symbols_command(),
#                 self.create_write_xml_file_command(),
#                 self.create_preprocess_command(),
#                 self.create_process_command(),
#                 self.create_postprocess_command()]
#         elif cmd_str == 'xml':
#             cmd_list = [self.create_write_xml_file_command()]
#         elif cmd_str == 'gen':
#             cmd_list = [self.generate_symbols()]
#         elif cmd_str == 'pre':
#             cmd_list = [self.create_preprocess_command()]
#         elif cmd_str == 'pro':
#             cmd_list = [self.create_process_command()]
#         elif cmd_str == 'pos':
#             cmd_list = [self.create_postprocess_command()]
#         elif cmd_str == 'cle':
#             cmd_list = [self.create_clean_command()]
#         # execute in turn
#         for cmd in cmd_list:
#             cmd.handle()
    
#     def create_generate_symbols_command(self):
#         return DoNothing("No symbols to be generated")

#     def create_write_xml_file_command(self):
#         return DoNothing("No XML file to be written")
        
#     def create_preprocess_command(self):
#         return DoNothing("No preprocessing to be done")
        
#     def create_process_command(self):
#         return DoNothing("No processing to be done")

#     def create_postprocess_command(self):
#         return DoNothing("No postprocessing to be done")
        
#     def create_clean_command(self):
#         return DoNothing("No cleaning to be done")
    
    
class TestSuite:
    """A front end for running various TestCommands."""
    def __init__(self, simulation_handler, results_handler, 
                 parameteriser=None, 
                 command_line_in_xml=None,
                 threshold_calculator=None,
                 simulator_name='fluidity',
                 options_filename_extension='.flml',
                 simulator_verbosity=0,
                 norm_calculator=None,
                 python_layer_verbosity=0):

        # set this module's verbosity
        verbosity(python_layer_verbosity)

        # set handlers 
        self.simulation_handler = simulation_handler
        self.results_handler = results_handler

        # set parameteriser
        if parameteriser is None:
            self.parameteriser = DefaultParameteriser(
                options_filename_extension,
                domain_length_1D=default_domain_length_1D)
        else:
            self.parameteriser = parameteriser

        # set other stuff
        self.command_line_in_xml = command_line_in_xml
        self.threshold_calculator = threshold_calculator
        self.simulator_name = simulator_name
        self.simulator_verbosity = simulator_verbosity
        self.norm_calculator = norm_calculator
        self.python_layer_verbosity = python_layer_verbosity
            

    def do(self, what):
        """General front end method"""
        cmd_str = str.lower(what[0:3])
        # decide which commands to run
        if cmd_str == 'all':
            cmd_list = [
                self.create_generate_symbols_command(),
                self.create_write_xml_file_command(),
                self.create_preprocess_command(),
                self.create_process_command(),
                self.create_postprocess_command()]
        elif cmd_str == 'xml':
            cmd_list = [self.create_write_xml_file_command()]
        elif cmd_str == 'gen':
            cmd_list = [self.generate_symbols()]
        elif cmd_str == 'pre':
            cmd_list = [self.create_preprocess_command()]
        elif cmd_str == 'pro':
            cmd_list = [self.create_process_command()]
        elif cmd_str == 'pos':
            cmd_list = [self.create_postprocess_command()]
        elif cmd_str == 'cle':
            cmd_list = [self.create_clean_command()]
        else:
            cmd_list = []
        # execute in turn
        for cmd in cmd_list:
            cmd.handle()

            
    def create_generate_symbols_command(self):
        # blank method.  Expect this to be overridden
        return DoNothing("No symbols to be generated")
        
    def create_write_xml_file_command(self):
        if self.command_line_in_xml is None:
            raise RuntimeError("""Tried to execute write_xml_file, but 
TestSuite was initialised without command_line_in_xml 
(a string that tells the test harness what to do). """)
        return WriteXMLFile(self.results_handler,
                            self.parameteriser,
                            self.command_line_in_xml,
                            self.threshold_calculator)
        
    def create_preprocess_command(self):
        cmd1 = ExpandOptionsTemplate(
             self.simulation_handler, self.parameteriser)
        cmd2 = ProcessMesh(self.simulation_handler, self.parameteriser)
        return TestCommandList(self.simulation_handler, 'Preprocessing',
                               self.parameteriser, [cmd1, cmd2])
        
    def create_process_command(self):
        return RunSimulation(self.simulation_handler,
                             self.parameteriser,
                             self.simulator_name,
                             self.simulator_verbosity)
        
    def create_postprocess_command(self):
        return WriteToReport(self.results_handler,
                             self.parameteriser,
                             self.norm_calculator)
        
    def create_clean_command(self):
        return CleanUp(self.simulation_handler, self.parameteriser)
