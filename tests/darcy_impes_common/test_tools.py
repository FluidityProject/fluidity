from batch_tools import Command, CommandList, verbose, HandlerLevel, \
    CompositeHandler, LeafHandler, default_fluidity_path
import os
import re
import subprocess
import numpy
import sys
from importlib import import_module
from getpass import getuser

default_fluidity_path()
from fluidity_tools import stat_parser

binary_verbosity = 0

error_rates_filename = "error_rates.txt"
error_norms_filename = "error_norms.txt"
cfl_max_filename = "cfl_max.txt"

## HELPER FUNCTIONS/CLASSES

def join_with_underscores(strings):
    result = None
    # convert string to tuple if it isn't already
    if isinstance(strings, str):
        strings = (strings, )
    for s in strings:
        # for flexibility, any nil values are ignored
        if s is None: continue
        # convert any non-strings
        if s is not str: s = str(s)
        if result is None:
            result = s
        else:
            result = result + '_' + s
    # for safety (making filenames), replace any points with 'p'
    result = re.sub('\.', 'p', result)
    return result


## FOR PUBLIC USE 
   
def find_in_open(report, key):
    """Searches a report and returns the value associated with key
    (ID). """
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

def find_norm(key): return find(error_norms_filename, key)
def find_rate(key): return find(error_rates_filename, key)
def find_cfl(key): return find(cfl_max_filename, key)


class WriteXMLFile(Command):
    def __init__(self, trigger_level_name,
                 norm_threshold=None, rate_threshold=None):
        self.trigger_level_name = trigger_level_name
        self.norm_threshold = norm_threshold
        self.rate_threshold = rate_threshold
        self.stem = None
        self.model = None
        self.saturation2_scale = None
        self.gravity_magnitude = None
        self.dim = None
        self.mesh_type = None
        self.phase_index = None
        self.var_name = None
        self.field_short = None
        self.norm = None
        self.mesh_suffix_0 = None


    def __del__(self):
        # finalise any outstanding open file
            try:
                self.write_xml_end()
            except:
                pass
            

    def write_xml_begin(self):
        case_name = join_with_underscores((
                self.stem, self.saturation2_scale, self.gravity_magnitude))
        self.xml_file = open(case_name+'.xml', 'w')
        self.xml_file.write("""<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE testproblem SYSTEM "regressiontest.dtd">

<testproblem>
  <name>{0}</name>
  <owner userid="{1}"/>
  <tags>diml</tags>
  <problem_definition length="short" nprocs="1">
    <command_line>
python processing.py pre proc post clean
    </command_line>
  </problem_definition>
  <variables/>
  <pass_tests>
    <test name="Solvers converged" language="python">
import os
files = os.listdir("./")
assert(not "matrixdump" in files and not "matrixdump.info" in files)
    </test>""".format(case_name, getuser()))

        
    def write_xml_end(self):
        self.xml_file.write("""
  </pass_tests>
  <warn_tests>
  </warn_tests>
</testproblem>""")
        self.xml_file.close()

            
    def write_xml_snippet(self, metric_type, key, rel_op, threshold):
        self.xml_file.write("""
    <test name="{1}: expect {0} {2} {3:g}" language="python">
from test_tools import find_{0}
assert(find_{0}("{1}") &{2}; {3:g})
    </test>""".format(metric_type, key, rel_op, 
                      threshold))

        
    def execute(self, level_name, value, indent):
        # store upper level details
        # TODO - remove the dependence on named members;
        #   have a client-specified (generic) list of levels
        if level_name == 'stem':
            self.stem = value
        elif level_name == 'saturation2_scale':
            self.saturation2_scale = value
        elif level_name == 'gravity_magnitude':
            self.gravity_magnitude = value

        # if appropriate, read in generated expressions
        # and initialise XML file
        if level_name == self.trigger_level_name:
            # close any outstanding open file
            try:
                self.write_xml_end()
            except:
                pass
            solution_name = join_with_underscores((
                self.stem, self.saturation2_scale, self.gravity_magnitude))
            solution_expressions = import_module(solution_name)
            self.solution_dict = solution_expressions.py_dict
            self.write_xml_begin()

        # continue to the lower levels
        if level_name == 'dim':
            self.dim = value
        elif level_name == 'mesh_type':
            self.mesh_type = value
        elif level_name == 'field':
            # split up field name into useful bits
            pattern = 'Phase(.)::(.*)'
            self.phase_index = re.sub(pattern, '\\1', value)
            self.var_name = re.sub(pattern, '\\2', value)
            self.field_short = str.lower(self.var_name)+self.phase_index
            # read in field magnitudes for computing the norm properly
            if self.norm_threshold is not None:
                self.rescaled_norm_threshold = self.norm_threshold* \
                     self.solution_dict[str.lower(self.var_name)+'_scale']
            
        elif level_name == 'norm':
            self.norm = value
            # important to signal that this is the first mesh before
            # iterating over mesh suffices (resolutions)
            self.mesh_suffix_0 = None
            
        if level_name == 'mesh_suffix':
            mesh_suffix = value
        
            key_stem = join_with_underscores((
                self.model, str(self.dim)+'d', self.mesh_type,
                self.field_short, 'l'+str(self.norm)))

            # write norm check
            if self.norm_threshold is not None:
                key = join_with_underscores((key_stem, mesh_suffix))
                self.write_xml_snippet('norm', key, 'lt',
                                   self.rescaled_norm_threshold)
                
            # write rate check
            if self.rate_threshold is not None:
                if self.mesh_suffix_0 is None:
                    self.mesh_suffix_0 = mesh_suffix
                    return
                key = join_with_underscores((key_stem,
                                             self.mesh_suffix_0+mesh_suffix))
                self.write_xml_snippet('rate', key, 'gt',
                                   self.rate_threshold)
                self.mesh_suffix_0 = mesh_suffix
        

    
class RunSimulation(Command):

    def __init__(self):
        self.darcy_impes_path = os.environ["FLUIDITYPATH"] + "bin/darcy_impes"
        if not os.path.isfile(self.darcy_impes_path): 
            raise IOError("Cannot find the darcy_impes binary.")
        # can add more levels to this list
        self.stem = None
        self.saturation2_scale = None
        self.gravity_magnitude = None
        self.model = None
        self.dim = None
        self.mesh_type = None

    def execute(self, level_name, value, indent):
        # can add more levels to this list
        if level_name == 'stem':
            self.stem = value
        elif level_name == 'saturation2_scale':
            self.saturation2_scale = value
        elif level_name == 'gravity_magnitude':
            self.gravity_magnitude = value
        elif level_name == 'model':
            self.model = value
        elif level_name == 'dim':
            self.dim = value
        elif level_name == 'mesh_type':
            self.mesh_type = value
        elif level_name == 'mesh_suffix':
            mesh_suffix = value

            casename = join_with_underscores((
                self.stem,
                self.saturation2_scale, self.gravity_magnitude, 
                self.model, str(self.dim)+'d',
                self.mesh_type, mesh_suffix))

            # start simulation (TODO: guard against absent mesh)
            if binary_verbosity > 0:
                subprocess.call([self.darcy_impes_path,
                                 '-v{0}'.format(binary_verbosity),
                                 '-l {0}.log'.format(casename),
                                 casename+'.diml'],
                                stdout=open(os.devnull, 'wb'))
            else:
                subprocess.call([self.darcy_impes_path,
                                 casename+'.diml'],
                                stdout=open(os.devnull, 'wb'))


class WriteToReport(Command):
    """Writes error norms and/or convergence rates to file"""

    def __init__(self, norms_file=None, rates_file=None):
        # both args are optional
        self.norms_file = norms_file
        self.rates_file = rates_file
        # can add more levels to this list
        self.stem = None
        self.model = None
        self.saturation2_scale = None
        self.gravity_magnitude = None
        self.dim = None
        self.mesh_type = None
        self.phase_index = None
        self.var_name = None
        self.field_short = None
        self.norm = None
        self.mesh_suffix = None
        self.mesh_suffix_0 = None
        
    def execute(self, level_name, value, indent):
        if level_name!='mesh_suffix':
            level_str = '{0}{1}: {2}\n'.format(indent, level_name, value)
            if self.norms_file is not None:
                self.norms_file.write(level_str)
            if self.rates_file is not None:
                self.rates_file.write(level_str)
            
        # can add more levels to this list
        if level_name == 'stem':
            self.stem = value
        elif level_name == 'saturation2_scale':
            self.saturation2_scale = value
        elif level_name == 'gravity_magnitude':
            self.gravity_magnitude = value
        elif level_name == 'model':
            self.model = value
        elif level_name == 'dim':
            self.dim = value
        elif level_name == 'mesh_type':
            self.mesh_type = value
        elif level_name == 'field':
            # split up field name into useful bits
            pattern = 'Phase(.)::(.*)'
            self.phase_index = re.sub(pattern, '\\1', value)
            self.var_name = re.sub(pattern, '\\2', value)
            self.field_short = str.lower(self.var_name)+self.phase_index
        elif level_name == 'norm':
            self.norm = value
            # important to signal that this is the first mesh before
            # iterating over mesh suffices (resolutions)
            self.mesh_suffix_0 = None

        elif level_name == 'mesh_suffix':
            self.mesh_suffix = value

            err = self.get_norm()
            key_stem = join_with_underscores((
                self.saturation2_scale, self.gravity_magnitude, 
                self.model, str(self.dim)+'d', self.mesh_type,
                self.field_short, 'l'+str(self.norm)))

            if self.norms_file is not None:
                key = join_with_underscores((key_stem, self.mesh_suffix))
                self.norms_file.write('{0}{1}{2:12.3e}\n'.format(
                    indent, key, err))
                if verbose(): sys.stdout.write ('   err: {0:.3e}'.\
                                                format(err))

            # can only start computing rates from the 2nd mesh
            if self.rates_file is not None and self.mesh_suffix_0 is not None:
                # print a new ID and the rate
                key = join_with_underscores((key_stem, 
                                             self.mesh_suffix_0+self.mesh_suffix))
                rate = numpy.log2(self.err0/err)
                self.rates_file.write('{0}{1}{2:12.6f}\n'.format(
                    indent, key, rate))
                if verbose(): sys.stdout.write ('   rate: {0:.6f}'.\
                                                format(rate))

            self.err0 = err
            self.mesh_suffix_0 = self.mesh_suffix


    def get_norm(self):
        # Clients may want to override this method
        filename = join_with_underscores((
            self.stem, self.saturation2_scale, self.gravity_magnitude, 
            self.model, str(self.dim)+'d', self.mesh_type,
            self.mesh_suffix)) + '.stat'
        if self.norm==1:
            self.calc_type = "integral"
        elif self.norm==2:
            self.calc_type = "l2norm"
        return stat_parser(filename)\
            ['Phase'+self.phase_index]\
            [self.var_name+'AbsError']\
            [self.calc_type][-1]



class WriteCFLToReport(Command):
    """Writes max CFL to file"""

    def __init__(self, cfl_max_file):
        self.cfl_max_file = cfl_max_file
        # can add more levels to this list
        self.stem = None
        self.model = None
        self.saturation2_scale = None
        self.gravity_magnitude = None
        self.dim = None
        self.mesh_type = None
        self.mesh_suffix = None
        
    def execute(self, level_name, value, indent):
        if level_name!='mesh_suffix':
            level_str = '{0}{1}: {2}\n'.format(indent, level_name, value)
            self.cfl_max_file.write(level_str)
            
        # can add more levels to this list
        if level_name == 'stem':
            self.stem = value
        elif level_name == 'saturation2_scale':
            self.saturation2_scale = value
        elif level_name == 'gravity_magnitude':
            self.gravity_magnitude = value
        elif level_name == 'model':
            self.model = value
        elif level_name == 'dim':
            self.dim = value
        elif level_name == 'mesh_type':
            self.mesh_type = value
        elif level_name == 'mesh_suffix':
            self.mesh_suffix = value

            cfl_max = max(self.get_cfl_max(1), self.get_cfl_max(2))
            key = join_with_underscores((
                self.saturation2_scale, self.gravity_magnitude, 
                self.model, str(self.dim)+'d', self.mesh_type,
                self.mesh_suffix))
            self.cfl_max_file.write('{0}{1}{2:12.3e}\n'.format(
                    indent, key, cfl_max))
            if verbose(): sys.stdout.write ('   cfl_max: {0:.3e}'.\
                                                format(cfl_max))

    def get_cfl_max(self, phase):
        filename = join_with_underscores((
            self.stem, self.saturation2_scale, self.gravity_magnitude, 
            self.model, str(self.dim)+'d', self.mesh_type,
            self.mesh_suffix)) + '.stat'
        return stat_parser(filename)\
            ['Phase'+str(phase)]\
            ['DarcyVelocityOverPorosityCFL']\
            ['max'][-1]

               
# class Dummy(Command):
#    def __init__(self):
#       pass
#    def execute(self, level_name, value, indent):
#       pass
