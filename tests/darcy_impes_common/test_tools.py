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
from batch_tools import SelfHandlingCommand, Validate, DoNothing, Tracker, \
    make_string, default_fluidity_path, set_global_verbosity
from copy import deepcopy
from glob import glob
from warnings import warn

default_fluidity_path()
from fluidity_tools import stat_parser

# domain dimensions may be invariant as far as the client is concerned.
# So make it (optionally) a global variable.
_global_test_dimensions = None
_global_key_formatting_dict = {'dim': '{0}d',
                               'norm': 'l{0}',
                               'mesh_res': 'm{0}'}

# a list of dictionaries to be maintained by parameterisers.  Names not
# appearing here must be explicitly registered by the client.
_built_in_dict_names = [
    'filenames', 'filename_stems', 'filename_extensions',
    'geometry', 'simulation_options']


# output filenames
_error_rates_filename = "error_rates.log"
_error_norms_filename = "error_norms.log"

_xml_skeleton_begin =  """<?xml version="1.0" encoding="UTF-8" ?>
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

_xml_skeleton_end = """
  </pass_tests>
  <warn_tests>
  </warn_tests>
</testproblem>"""

_xml_skeleton_assert = """
    <test name="{1}: expect {0} {2} {3:g}" language="python">
from test_tools import find_{0}
assert(find_{0}("{1}") &{2}; {3:g})
    </test>"""

def get_global_num_dimensions():
    return _global_num_dimensions

def get_global_python_verbosity():
    return _global_python_verbosity

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
    return find(_error_norms_filename, key)

def find_rate(key):
    return find(_error_rates_filename, key)

    
class Parameteriser:
   """A class for determining filenames and making dictionaries for
   template file expansion.  Clients can extend the class to customise
   any of the methods.  The parameteriser generally needs to be kept up
   to date, via the update() method, before methods are called."""

   def __init__(self, simulation_options_filename_extension,
                domain_extents=(1., 1., 1.), finish_time=1.,
                reference_mesh_res=10, reference_time_step_number=10):
       self.__simulation_options_filename_extension = \
           simulation_options_filename_extension
       self.__domain_extents = domain_extents
       self.__finish_time = finish_time
       self.__reference_mesh_res = reference_mesh_res
       self.__reference_time_step_number = reference_time_step_number
       # see update_dict() for an explanation of how the following
       # variables are used
       self.__master_dict = {}
       for d in _built_in_dict_names:
           self.register_dict(d)
       self.__num_expansions = {}

   def update(self, test_command):
       self.__test_command = test_command

   def get(self, what):
        """If what=='level', returns current level name.  If what=='node',
        returns current node name.  Otherwise assumes the argument is a
        level name from which a node name is to be returned.  Return
        value defaults to None."""
        # delegate
        return self.__test_command.get(what)

   def at(self, what):
        """Returns True if where=='root' and we are at the base of the
        handler tree, if where=='leaf' and we are at the end of a
        branch, or if where is the current level_name.
        """
        # delegate
        return self.__test_command.at(where)

   def make_key(self, level_names=None, exclude=[], substitute={}):
        """Forms a string made up of node names corresponding to the
        supplied list of level names.
        """ 
        # delegate
        return self.__test_command.make_key(level_names, exclude, substitute)

   def register_dict(self, dict_name, user_defined_dict={}):
       """Registers a custom named dictionary for use by lookup() and
       expand_template()."""
       try:
           # try updating first so that we do not clobber an existing
           # dictionary
           self.update_dict(dict_name, user_defined_dict)
       except:
           self.__master_dict[dict_name] = user_defined_dict
   
   def update_dict(self, dict_name, user_defined_dict={}):
       """Updates a named dictionary for use by lookup() and
       expand_template().  Note that the passed user-defined dict may
       get merged with a corresponding built-in dict (e.g. 'geometry',
       'simulation').  In fact, the dict argument is optional, in which
       case any built-in dict will get updated according to the current
       position in the handler."""
       # in special cases, start with a built-in dict
       if dict_name=='filenames':
           result_dict = self.filenames_dict()
       elif dict_name=='filename_stems':
           result_dict = self.filename_stems_dict()
       elif dict_name=='filename_extensions':
           result_dict = self.filename_extensions_dict()
       elif dict_name=='geometry':
           result_dict = self.geometry_dict()
       elif dict_name=='simulation_options':
           result_dict = self.simulation_options_dict()
       else:
           result_dict = {}
       # update/merge with the optional dictionary in the argument list
       result_dict.update(user_defined_dict)
       # update/merge the corresponding sub-dictionary in master
       self.__master_dict[dict_name].update(result_dict)
       
   def lookup(self, dict_name, key):
       """Looks up a variable in the dictionary corresponding to
       dict_name."""
       # update the embedded dict, for safety
       self.update_dict(dict_name)
       sub_dict = self.__master_dict[dict_name]
       return sub_dict[key]

   
   def expand_template(self, buffer, dict_name):
       """In a given buffer, replaces all instances of strings prefixed
       '$' with the corresponding values in the dictionary corresponding
       to dict_name."""
       max_num_loops = 5
       # update the embedded dict, for safety
       self.update_dict(dict_name)
       D = self.__master_dict[dict_name]
       if len(D) == 0:
           warn("\nNo entries found in {0} dictionary.".format(dict_name))
       # substitute until there are no more placeholders
       n = 0
       while '$' in buffer and n < max_num_loops:
           buffer = Template(buffer)
           buffer = buffer.safe_substitute(D)
           n += 1
       if n == max_num_loops:
           warn("\nMaximum number of substitutions ({0}) reached".format(n))
       return buffer

   def filename_stems_dict(self):
      """Built-in dictionary for looking up filename stems for options,
      geometries and meshes."""
      return {
          # assume just one case name
          'simulation_options_template':
              make_string(['template', self.get('root')]),
          'simulation_options': self.make_key(),
          'geometry_template':
              make_string(['template', self.make_key(['dim', 'mesh_type'])]),
          'geometry': self.make_key(['dim', 'mesh_type', 'mesh_res']),
          'mesh': self.make_key(['dim', 'mesh_type', 'mesh_res']) }

   def filename_extensions_dict(self):
      """Built-in dictionary for looking up filename extensions for
      options, geometries and meshes."""
      geo_ext = '.geo'
      mesh_ext_list = ['.msh']
      try:
          if self.dim()=="1":
              geo_ext = '.sh'
              mesh_ext_list = ['.bound', '.ele', '.node']
      except KeyError:
          pass
      return {
          'simulation_options_template':
              self.__simulation_options_filename_extension,
          'simulation_options':
              self.__simulation_options_filename_extension,
          'geometry': geo_ext,
          'geometry_template': geo_ext,
          'mesh': mesh_ext_list }
           
   def filenames_dict(self):
      """Built-in dictionary for looking up filenames for options,
      geometries and meshes."""
      S = self.filename_stems_dict()
      E = self.filename_extensions_dict()
      for k in S.keys():
          try:
              # join stem with extension
              S[k] = S[k] + E[k]
          except:
              # could end up with a list
              S[k] = [S[k] + e for e in E[k]]
      return S
   
   def geometry_dict(self):
      """Built-in dictionary for expanding a geometry template.
      Currently minimal, with meshing support for a rectilinear domain
      only.  The client may want to override this method to define new
      entries"""
      # initialise result
      result = {'MESH_NAME': self.lookup('filename_stems', 'mesh')}
      for i, dim_str in enumerate(('X', 'Y', 'Z')):
          try:
              L = self.domain_extents()[i]
          except:
              continue
          n = int(self.element_number(L))
          result.update({
                  'DOMAIN_LENGTH_{0}'.format(dim_str) : str(L),
                  'EL_NUM_{0}'.format(dim_str) : str(n),
                  'EL_SIZE_{0}'.format(dim_str) : str(L/n) })
      return result
   
   def simulation_options_dict(self):
      """Built-in dictionary for expanding a simulation options template.
      Currently minimal.  The client may want to override this method to
      define new entries.
      """
      mesh_name = self.lookup('filename_stems', 'mesh')
      if self.dim()=="1":
         mesh_format = 'triangle'
      else:
         mesh_format = 'gmsh'
      return {
          'MESH_NAME': mesh_name,
          'MESH_FORMAT': mesh_format,
          'DOMAIN_DIM': self.dim(),
          'TIME_STEP': str(self.time_step()),
          'FINISH_TIME': str(self.finish_time()) }

   def element_number(self, domain_length):
      """Helper method.  Calculate element number along domain edge,
      keeping dx, dy and dz approximately the same but scaling
      consistently to higher mesh resolutions."""
      ref_n = numpy.round(self.__reference_mesh_res * domain_length / \
                              self.domain_extents()[0])
      return int(ref_n * int(self.mesh_res()) / self.__reference_mesh_res)

   def time_step(self):
      """Helper method.  Calculate an appropriate time step, maintaining
      a constant Courant number for all mesh resolutions."""
      ref_dt = self.finish_time()/self.__reference_time_step_number
      return ref_dt * self.__reference_mesh_res/int(self.mesh_res())

   def dim(self):
       """Clients may want to override this."""
       return self.get('dim')

   def mesh_res(self):
       """Clients may want to override this."""
       return self.get('mesh_res')

   def domain_extents(self):
       """Clients may want to override this."""
       return self.__domain_extents

   def finish_time(self):
       """Clients may want to override this."""
       return self.__finish_time
       
      
class TestCommand(SelfHandlingCommand):
    """This abstract class, to be extended by concrete Commands, has
    some test-related implementation and also takes on the role of
    storing level details from which to make labels and filenames."""

    def __init__(self, requisite_level_names, handler, message, 
                 parameteriser):
        SelfHandlingCommand.__init__(self, requisite_level_names, 
                                     handler, message)
        self.__tracker = Tracker(_global_key_formatting_dict)
        self.__parameteriser = parameteriser
        # check for dimension setting clash
        if 'dim' in requisite_level_names and \
                _global_test_dimensions() is not None:
            warn("""
global_test_dimensions was set, but here 'dim' is included as a requisite
level, suggesting that several domain dimensions may be encountered.""")

    def inform(self, level_name, node_name, at_leaf, indent, verbose):
        """Records level details, including subtle testing-related details."""

        if level_name == 'mesh_res':
            # mesh_res is a special level; we need to update mesh_res_prev
            try:
                res = self.__tracker.get('mesh_res')
                self.__tracker.update('mesh_res_prev', res)
            except:
                pass

        # in any case record level-node details in both the tracker
        # object and the superclass
        self.__tracker.inform(level_name, node_name, at_leaf, indent)
        SelfHandlingCommand.inform(self, level_name, node_name, at_leaf,
                                   indent, verbose)

        # field is also special; it needs to be reformatted/split into
        # useful bits
        if level_name == 'field':
            pattern = 'Phase(.)::(.*)'
            ph = sub(pattern, '\\1', node_name)
            var = sub(pattern, '\\2', node_name)
            # reformat the field name
            self.__tracker.update('field', str.lower(var) + ph)
            # have a separate dictionary so as not to pollute the main one
            self.field_dict = {'field_long':  node_name,
                               'phase': 'Phase' + ph,
                               'var':   var}
           
    def get(self, what):
        """If what=='level', returns current level name.  If what=='node',
        returns current node name.  Otherwise assumes the argument is a
        level name from which a node name is to be returned.  Return
        value defaults to None."""
        # dimension is a special case
        if what=='dim' and _global_test_dimensions:
            return str(_global_test_dimensions)
        try:
            # try consulting the special field_dict
            return self.field_dict[what]
        except:
            # delegate to the tracker
            return self.__tracker.get(what)

    def at(self, where):
        """Returns True if where=='root' and we are at the base of the
        handler tree, if where=='leaf' and we are at the end of a
        branch, or if where is the current level_name.
        """
        # delegate
        return self.__tracker.at(where)

    def make_key(self, level_names=None, exclude=[], substitute={}):
        """Forms a string made up of node names corresponding to the
        supplied list of level names.
        """ 
        # delegate
        return self.__tracker.make_key(level_names, exclude, substitute)

    def register_dict(self, dict_name, dict_value):
        """Registers a custom named dictionary for use by lookup() and
        expand_template()."""
        # delegate
        self.__parameteriser.update(self)
        self.__parameteriser.register_dict(dict_name, dict_value)

    def update_dict(self, dict_name, dict_value):
        """Updates a named dictionary for use by lookup() and
        expand_template().  Note that the passed user-defined dict may
        get merged with a corresponding built-in dict (e.g. 'geometry',
        'simulation').  In fact, the dict argument is optional, in which
        case any built-in dict will get updated according to the current
        position in the handler."""
        # delegate
        self.__parameteriser.update(self)
        self.__parameteriser.update_dict(dict_name, dict_value)

    def lookup(self, dict_name, key):
        """Looks up a key in the dictionary stored under dict_name."""
        # delegate
        self.__parameteriser.update(self)
        return self.__parameteriser.lookup(dict_name, key)

    def expand_template(self, buffer, dict_name):
        """In a given buffer, replaces all instances of strings prefixed
        '$' with the corresponding values in the dictionary corresponding
        to dict_name."""
        # delegate
        self.__parameteriser.update(self)
        return self.__parameteriser.expand_template(buffer, dict_name)

   
class WriteXMLFile(TestCommand):
    """Writes out an XML file which will be picked up by Fluidity's test
    harness."""
    
    def __init__(self, handler, parameteriser, command_line,
                 assert_threshold_calculator=None):
        TestCommand.__init__(self, ['field', 'norm', 'mesh_res'], handler, 
                             "Generating XML file", parameteriser)
        self.command_line = command_line
        if assert_threshold_calculator is None:
            self.assert_threshold_calculator = AssertThresholdCalculator()
        else:
            self.assert_threshold_calculator = assert_threshold_calculator
            
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
        norm_threshold = self.assert_threshold_calculator.\
                         calculate_norm_threshold(self)
        if norm_threshold is not None:
            key = self.make_key(exclude=['root'])
            self.write_xml_assert('norm', key, 'lt', norm_threshold)
            self.write_message('  norm < {0}'.format(norm_threshold))
                
        # write rate assert
        rate_threshold = self.assert_threshold_calculator.\
                         calculate_rate_threshold(self)
        if rate_threshold is not None:
            mesh_res = self.get('mesh_res')
            if self.mesh_res_prev is None:
                self.mesh_res_prev = mesh_res
                return
            subs_dict = {'mesh_res': '{0}_{1}'.format(
                self.mesh_res_prev, mesh_res)}
            key = self.make_key(exclude=['root'],
                                substitute=subs_dict)
            self.write_xml_assert('rate', key, 'gt', rate_threshold)
            self.write_message('  rate > {0}'.format(rate_threshold))
            self.mesh_res_prev = mesh_res

    def write_xml_begin(self):
        case_name = self.get('root')
        self.xml_file = open(case_name+'.xml', 'w')
        simulation_options_extension = self.lookup('filename_extensions',
                                                   'simulation_options')
        self.xml_file.write(_xml_skeleton_begin.format(
                case_name, getuser(),
                simulation_options_extension,
                self.command_line))
        
    def write_xml_end(self):
        try:
            self.xml_file.write(_xml_skeleton_end)
            self.xml_file.close()
        except AttributeError:
            pass

    def write_xml_assert(self, metric_type, key, rel_op, threshold):
        self.xml_file.write(_xml_skeleton_assert.format(
            metric_type, key, rel_op, threshold))
        
        
class AssertThresholdCalculator:
    """Can be overriden by clients to define thresholds in assert
    statements.
    """
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
            return self.norm_threshold * \
                test_command.lookup('solution', var_name)
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

        src_filename = self.lookup('filenames',
                                   'simulation_options_template')
        tgt_filename = self.lookup('filenames', 'simulation_options')
        with open(src_filename) as src:
            buf = str(src.read())
        buf = self.expand_template(buf, 'simulation_options')
        with open(tgt_filename, 'w') as tgt:
            tgt.write(buf)
        self.write_message('done')
      
            
class ProcessMesh(TestCommand):
    def __init__(self, handler, parameteriser):
        TestCommand.__init__(self, ['mesh_res'], handler,
                             "Meshing", parameteriser)
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
        src_filename = self.lookup('filenames', 'geometry_template')
        tgt_filename = self.lookup('filenames', 'geometry')
        with open(src_filename) as src:
            buf = str(src.read())
        buf = self.expand_template(buf, 'geometry')
        with open(tgt_filename, 'w') as tgt:
            tgt.write(buf)

        # call interval or gmsh
        # TODO: tell user to 'make geo' if needed
        if self.get('dim') == '1':
            lx_str = self.lookup('geometry', 'DOMAIN_LENGTH_X')
            dx_str = self.lookup('geometry', 'EL_SIZE_X')
            call([self.interval_path, '0.0',
                  lx_str, '--dx='+dx_str,
                  self.lookup('filename_stems', 'geometry')])
        else:
            call(['gmsh', '-'+self.get('dim'),
                  self.lookup('filenames', 'geometry')],
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
        options_filename = self.lookup('filenames', 'simulation_options')
        simulation_name = self.lookup('filename_stems', 'simulation_options')
        if self.simulator_verbosity > 0:
            call([self.simulator_path,
                  '-v{0}'.format(simulator_verbosity),
                  '-l {0}.log'.format(simulation_name), options_filename],
                 stdout=open(os.devnull, 'wb'))
        else:
            call([self.simulator_path, options_filename],
                 stdout=open(os.devnull, 'wb'))
        self.write_message('done')

        
class WriteToReport(TestCommand):
    """Writes error norms and/or convergence rates to files."""

    def __init__(self, handler, parameteriser, norm_calculator=None):
        if norm_calculator is None:
            self.norm_calculator = NormCalculator()
        else:
            self.norm_calculator = norm_calculator
        # ordered merge
        rl = list(OrderedDict.fromkeys(
                self.norm_calculator.get_requisite_levels() + ['mesh_res']))
        TestCommand.__init__(self, rl, handler, 
                             "Writing report(s)", parameteriser)

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
            self.norms_file = open(_error_norms_filename, 'w')
            self.rates_file = open(_error_rates_filename, 'w')

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
            if err is None:
                return
            key = self.make_key(exclude=['root'])
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
            key = self.make_key(exclude=['root'],
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
    """Can be overridden by clients so that they can define norms their
    own way."""
    def get_requisite_levels(self):
        return ['field', 'norm', 'mesh_res']
    
    def calculate_norm(self, test_command):
        """if the client intends to examine errors of variable X, he or
        she needs to represent the errors in a diagnostic field called
        XAbsError."""
        filename = test_command.make_key(exclude=['field','norm']) + '.stat'
        if test_command.get('norm')=='1':
            calc_type = "integral"
        elif test_command.get('norm')=='2':
            calc_type = "l2norm"
        phase = test_command.get('phase')
        var = test_command.get('var')+'AbsError'
        try:
            return stat_parser(filename)[ph][var][calc_type][-1]
        except Exception as e:
            e.args += (
                "NormCalculator expected to find {0}::{1} in the stat file;".\
                    format(phase, var) + \
                    " has this been defined in the options file?",)
            raise
        
            
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

        # try removing input files, geometry files, mesh files
        trash = [self.lookup('filenames', 'simulation_options'),
                 self.lookup('filenames', 'geometry'),
                 self.lookup('filenames', 'mesh')]
        self.remove_recursively(trash, True)

        # try removing results files
        results_pattern = \
            self.lookup('filename_stems', 'simulation_options')+'_*.vtu'
        trash = glob(results_pattern)
        if trash:
            self.remove_recursively(trash, False)
            self.write_message('   removed {0}'.format(results_pattern),
                               with_newline=True)

                
    def remove_recursively(self, filename_or_list, verbose):
        for f in filename_or_list:
            try:
                os.remove(f)
                if verbose:
                    self.write_message('   removed {0}'.format(f),
                                       with_newline=True)
            except TypeError:
                self.remove_recursively(f, verbose)
            except OSError:
                pass


class LumpedTestCommand(TestCommand):
    """Chains multiple Commands together as a single TestCommand.  The
    LumpedTestCommand's handler overrides any handlers associated with the
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
        if self.at('leaf'):
            self.write_message('...')
        for cmd in self.commands:
            cmd.execute()
        if self.at('leaf'):
            self.write_message('done')

    def update_dict(self, dict_name, dict_value):
        """Delegates to wrapped commands in addition to calling the
        original superclass method."""
        TestCommand.update_dict(self, dict_name, dict_value)
        for cmd in self.commands:
            cmd.update_dict(dict_name, dict_value)

    
class TestSuite:
    """A front end for running various TestCommands."""
    def __init__(self, mesh_handler, simulation_handler, 
                 results_handler, parameteriser=None, 
                 command_line_in_xml=None,
                 assert_threshold_calculator=None,
                 simulator_name='fluidity',
                 simulation_options_filename_extension='.flml',
                 simulator_verbosity=0,
                 norm_calculator=None,
                 global_test_dimensions=_global_test_dimensions,
                 global_key_formatting_dict=_global_key_formatting_dict,
                 python_verbosity=None):

        # set module variables
        global _global_test_dimensions, _global_key_formatting_dict
        _global_test_dimensions = global_test_dimensions
        _global_key_formatting_dict = global_key_formatting_dict
        if python_verbosity is not None:
            set_global_verbosity(python_verbosity)

        # set handlers 
        self.mesh_handler = mesh_handler
        self.simulation_handler = simulation_handler
        self.results_handler = results_handler

        # set parameteriser
        if parameteriser is None:
            self.parameteriser = Parameteriser(
                simulation_options_filename_extension)
        else:
            self.parameteriser = parameteriser

        # set other stuff
        self.command_line_in_xml = command_line_in_xml
        self.assert_threshold_calculator = assert_threshold_calculator
        self.simulator_name = simulator_name
        self.simulator_verbosity = simulator_verbosity
        self.norm_calculator = norm_calculator
            

    def do(self, what):
        """General front end method"""
        cmd_str = str.lower(what[0:3])
        # decide which commands to run
        if cmd_str == 'all':
            cmd_list = [
                self.create_generate_symbols_command(),
                self.create_write_xml_file_command(),
                self.create_expand_options_command(),
                self.create_mesh_command(),
                self.create_simulate_command(),
                self.create_write_report_command()]
        elif cmd_str == 'xml':
            cmd_list = [self.create_write_xml_file_command()]
        elif cmd_str == 'gen':
            cmd_list = [self.create_generate_symbols_command()]
        elif cmd_str == 'pre':
            cmd_list = [self.create_expand_options_command(),
                        self.create_mesh_command()]
        elif cmd_str == 'pro':
            cmd_list = [self.create_simulate_command()]
        elif cmd_str == 'pos':
            cmd_list = [self.create_write_report_command()]
        elif cmd_str == 'cle':
            cmd_list = [self.create_clean_command()]
        else:
            raise ValueError(
                what + ' is not recognised as an argument to TestSuite.do')
        # execute in turn
        for cmd in cmd_list:
            cmd.handle()
            
    def create_generate_symbols_command(self):
        # blank method.  Expect this to be overridden
        return DoNothing("No symbols to be generated")
        
    def create_write_xml_file_command(self):
        if self.command_line_in_xml is None:
            return DoNothing("""No XML file to be generated
(TestSuite was initialised without command_line_in_xml)""")
        else:
            return WriteXMLFile(self.results_handler,
                                self.parameteriser,
                                self.command_line_in_xml,
                                self.assert_threshold_calculator)

    def create_mesh_command(self):
        return ProcessMesh(self.mesh_handler, self.parameteriser)
        
    def create_expand_options_command(self):
        return ExpandOptionsTemplate(
            self.simulation_handler, self.parameteriser)
        
    def create_simulate_command(self):
        return RunSimulation(self.simulation_handler,
                             self.parameteriser,
                             self.simulator_name,
                             self.simulator_verbosity)
        
    def create_write_report_command(self):
        return WriteToReport(self.results_handler,
                             self.parameteriser,
                             self.norm_calculator)
        
    def create_clean_command(self):
        return CleanUp(self.simulation_handler, self.parameteriser)
