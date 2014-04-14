from batch_tools import Command, CommandList, verbose, HandlerList, CompositeHandler
import os
import re
import subprocess
import numpy
import fluidity_tools

def join_with_underscores(strings):
    result = None
    for s in strings:
        # for flexibility, any nil values are ignored
        if s is None: continue
        if result is None:
            result = s
        else:
            result = result + '_' + s
    return result

         
class Dummy(Command):
   def __init__(self):
      pass
   def execute(self, level_name, value, indent):
      pass      
            

class RunSimulation(Command):

    def __init__(self):
        self.binary_path = "../../bin/darcy_impes"
        if not os.path.isfile(self.binary_path): 
            raise IOError("Cannot find the binary.")
        # can add more levels to this list
        self.stem = None
        self.model = None
        self.dim = None
        self.mesh_type = None

    def execute(self, level_name, value, indent):
        # can add more levels to this list
        if level_name == 'stem':
            self.stem = value
        elif level_name == 'model':
            self.model = value
        elif level_name == 'dim':
            self.dim = value
        elif level_name == 'mesh_type':
            self.mesh_type = value
        elif level_name == 'mesh_suffix':
            mesh_suffix = value

            filename = join_with_underscores((
                self.stem, self.model, str(self.dim)+'d',
                self.mesh_type, mesh_suffix)) + '.diml'
            
            # start simulation (TODO: guard against absent mesh)
            subprocess.call([self.binary_path, filename,
                             '-v3'], stdout=open(os.devnull, 'wb'))


class WriteToReport(Command):

    def __init__(self, norms_file=None, rates_file=None):
        # both args are optional
        self.norms_file = norms_file
        self.rates_file = rates_file
        # can add more levels to this list
        self.stem = None
        self.model = None
        self.dim = None
        self.mesh_type = None
        self.phase_index = None
        self.var_name = None
        self.field_short = None
        self.norm = None
        self.mesh_suffix = None
        
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
        elif level_name == 'mesh_suffix':
            self.mesh_suffix = value

            err = self.get_norm()

            if self.norms_file is not None:
                key = join_with_underscores((
                    self.model, str(self.dim)+'d', self.mesh_type,
                    self.field_short, 'l'+str(self.norm),
                    self.mesh_suffix))
                self.norms_file.write('{0}{1}{2:12.3e}\n'.format(
                    indent, key, err))

            # can only start computing rates from the 2nd mesh
            if self.rates_file is not None and self.mesh_suffix!='A':
                # print a new ID and the rate
                key = join_with_underscores((
                    self.model, str(self.dim)+'d', self.mesh_type,
                    self.field_short, 'l'+str(self.norm),
                    self.mesh_suffix0+self.mesh_suffix))
                self.rates_file.write('{0}{1}{2:12.6f}\n'.format(
                    indent, key, numpy.log2(self.err0/err)))
            self.err0 = err
            self.mesh_suffix0 = self.mesh_suffix

    def get_norm(self):
        filename = join_with_underscores((
            self.stem, self.model, str(self.dim)+'d',
            self.mesh_type, self.mesh_suffix)) + '.stat'
        if self.norm==1:
            self.calc_type = "integral"
        elif self.norm==2:
            self.calc_type = "l2norm"
        return fluidity_tools.stat_parser(filename)\
            ['Phase'+self.phase_index]\
            [self.var_name+'AbsError']\
            [self.calc_type][-1]
