#!/usr/bin/env python

"""manufactured_solution_test_tools.py

Documentation to come."""

import sys
sys.path.append('../../python')

import os
import string
import subprocess
import re
import numpy
import glob
from test_tools import Command, CommandList, HandlerLevel, \
   CompositeHandler, LeafHandler, WriteXMLSnippet, RunSimulation, \
   WriteToReport, error_norms_filename, error_rates_filename, verbose
from fluidity_tools import stat_parser
from getpass import getuser

verbose(True)
debug = False

# constants
mesh_template_prefix = 'template'
domain_shape_dict = {1:'line', 2:'rectangle', 3:'cuboid'}
mesh_res_dict = {'A':5, 'B':10, 'C':20, 'D':40,
                            'E':80, 'F':160}


class Parameterisation:
   """Parameters determined from a small, fixed set of
   variables can be determined here.
   """

   def __init__(self, case, dim, mesh_type, mesh_suffix):
      
      # basic parameterisations
      self.dim = dim
      self.domain_shape = domain_shape_dict[dim]
      self.mesh_res = mesh_res_dict[mesh_suffix]

      # filenames
      if dim==1:
         m_t = ''
         g_e = '.sh'
         m_e_l = ['.bound', '.ele', '.node']
      else:
         m_t = '_'+mesh_type
         g_e = '.geo'
         m_e_l = ['.msh']
      self.mesh_name = self.domain_shape+m_t+'_'+mesh_suffix
      self.options_template_filename = "template_"+case+".diml"
      self.options_name = case+'_'+str(dim)+'d'+m_t+'_'+mesh_suffix
      self.options_filename = self.options_name+'.diml'
      self.geo_template_filename = mesh_template_prefix+'_'+\
                                   self.domain_shape+m_t+g_e
      self.geo_filename = self.mesh_name+'.geo'
      self.mesh_filenames = []
      for m_e in m_e_l:
         self.mesh_filenames.append(self.mesh_name+m_e)

      
   def compute_more_options(self):
      """Options template expansion typically requires more stuff.
      It is placed here so as not to burden other clients.
      """
      
      # dictionary entries
      if self.dim==1:
         self.mesh_format = 'triangle'
         # for collapsing inappropriate wall BCs
         self.begin_rm_if_1D = '<!--'
         self.end_rm_if_1D = '-->'
      else:
         self.mesh_format = 'gmsh'
         self.begin_rm_if_1D = ''
         self.end_rm_if_1D = ''
         
      # more dictionary entries
      if self.dim==1:
         self.outlet_ID = '2'
         self.inlet_ID = '1'
         self.wall_IDs = ''
         self.wall_num = '0'
      elif self.dim==2:
         self.outlet_ID = '11'
         self.inlet_ID = '13'
         self.wall_IDs = '12 14'
         self.wall_num = '2'
      elif self.dim==3:
         self.outlet_ID = '28'
         self.inlet_ID = '29'
         self.wall_IDs = '30 31 32 33'
         self.wall_num = '4'


class MSTestToolsCommand(Command):
   """This abstract class, to be extended by concrete Commands, has a
method for doing the boring job of storing level details, i.e. recording
where the Command is in the Handler tree.
   """
   
   def __init__(self):
      self.stem = None
      self.dim = None
      self.mesh_type = None
      self.mesh_suffix = None
      
   def store_level_details(self, level_name, value):
      if level_name == 'stem':
         self.stem = value
      if level_name == 'dim':
         self.dim = value
      if level_name == 'mesh_type':
         self.mesh_type = value
      if level_name == 'mesh_suffix':
         self.mesh_suffix = value
      
   
      
class ExpandOptionsTemplate(MSTestToolsCommand):
   
   def __init__(self, finish_time, solution_dict=None):
      MSTestToolsCommand.__init__(self)
      self.finish_time = finish_time
      # solution_dict may come from an external file generated as part of a
      # manufactured solution.
      self.solution_dict = solution_dict
      
   def execute(self, level_name, value, indent):
      # all levels record level details
      self.store_level_details(level_name, value)
      if level_name == 'mesh_suffix':
         # the last level does the important stuff
         
         # compute useful parameters
         param = Parameterisation(self.stem, self.dim,
                                  self.mesh_type, self.mesh_suffix)
         param.compute_more_options()

         # compose the dictionary
         options_dict = {
            'MESH_FILE': param.mesh_name,
            'MESH_FORMAT': param.mesh_format,
            'DOMAIN_DIM': str(self.dim),

   'TIME_STEP': str(self.compute_time_step()),

   'FINISH_TIME': self.finish_time,
            'INLET_ID': param.inlet_ID,
            'OUTLET_ID': param.outlet_ID,
            'WALL_IDS': param.wall_IDs,
            'WALL_NUM': param.wall_num,
            'BEGIN_RM_IF_1D': param.begin_rm_if_1D,
            'END_RM_IF_1D': param.end_rm_if_1D }
         
         # and use Python's string.Template to do the expansion
         with open(param.options_template_filename) as src:
            buf = string.Template(src.read())
            buf = buf.safe_substitute(options_dict)
         # repeat with the solution dict
         try:
            buf = string.Template(buf)
            buf = buf.safe_substitute(self.solution_dict)
         except Attribute_error:
            pass
         with open(param.options_filename, 'w') as tgt:
            tgt.write(buf)

   def compute_time_step(self):
      """Extracted so as to be overrideable"""
      
      # let the coarsest mesh have one time step and maintain a
      # constant Courant number for all meshes
      scale_factor = float(mesh_res_dict['A'])/ \
                     float(mesh_res_dict[self.mesh_suffix])
      return scale_factor * self.finish_time
      
            
class ProcessMesh(MSTestToolsCommand):

   def __init__(self, domain_extents):
      MSTestToolsCommand.__init__(self)
      self.binary_path = "../../bin/darcy_impes"
      self.domain_extents = domain_extents
      
   def execute(self, level_name, value, indent):
      # all levels record level details
      self.store_level_details(level_name, value)
      if level_name == 'mesh_suffix':
         # the last level does the important stuff

         # compute useful parameters
         param = Parameterisation(self.stem, self.dim, self.mesh_type,
                                  self.mesh_suffix)
         
         # compose the dictionary
         geo_dict = {'MESH_NAME': str(param.mesh_name)}
         for i, dim_str in enumerate(('X', 'Y', 'Z')):
            # calculate element number along domain edge, keeping dx, dy
            # and dz approximately the same but scaling consistently to
            # higher mesh resolutions
            D = self.domain_extents[i]
             # reference number along x-edge
            nx_A = mesh_res_dict['A']
            # reference number along edge in this dimension
            n_A = numpy.round(nx_A * D / self.domain_extents[0])
            # number scaled up for this mesh resolution
            n = int(n_A * param.mesh_res / nx_A)
            geo_dict.update({
               'DOMAIN_LENGTH_'+dim_str : str(D),
               'EL_NUM_'+dim_str : str(n),
               'EL_SIZE_'+dim_str : str(D/n) })
            
         # expand mesh template
         with open(param.geo_template_filename) as src:
            buf = string.Template( src.read() )
         buf = buf.safe_substitute(geo_dict)
         with open(param.geo_filename, 'w') as tgt:
            tgt.write(buf)

         # call interval or gmsh
         if self.dim==1:
            lx=self.domain_extents[0]
            dx=lx/float(param.mesh_res)
            subprocess.call(['../../bin/interval', '0.0',
                             str(lx), '--dx='+str(dx),
                             param.mesh_name])
         else:
            subprocess.call(['gmsh', '-'+str(self.dim),
                             param.geo_filename],
                            stdout=open(os.devnull, 'wb'))

            
class CleanUp(MSTestToolsCommand):
   """Removes generated geometry, mesh, and options files, but leaves
    results intact.  This ensures templates do not get deleted
    (compared with rm -f *.geo, etc)
   """

   def __init__(self):
      MSTestToolsCommand.__init__(self)
      self.stem = None
      self.dim = None
      self.mesh_type = None
      
   def execute(self, level_name, value, indent):
      # all levels record level details
      self.store_level_details(level_name, value)
      if level_name == 'mesh_suffix':
         # the last level does the important stuff

         # compute useful parameters
         param = Parameterisation(self.stem, self.dim, self.mesh_type,
                                  self.mesh_suffix)

         # try removing input files, mesh files, results files
         for f in [param.options_filename, param.geo_filename] + \
             param.mesh_filenames:
            try:
               os.remove(f)
               if verbose(): sys.stdout.write('\n'+indent+'   removed '+f)
            except OSError:
               pass

         
class ManufacturedSolutionTestSuite:
   def __init__(self, case_name, solution_dict, finish_time,
                mesh_type_list, mesh_suffix_list_per_dimension,
                field_name_list, norm_list):

      self.case_name = case_name
      self.solution_dict = solution_dict
      self.finish_time = finish_time

      # needed for e.g. creating 1D meshes, but do not rely on this
      # for true mesh dimensions otherwise
      self.domain_extents = solution_dict['domain_extents']
      
      # build handlers
      mesh_type_handler_level = HandlerLevel('mesh_type', mesh_type_list)
      field_handler_level = HandlerLevel('field', field_name_list)
      norm_handler_level = HandlerLevel('norm', norm_list)
      # maintain one tree for running simulations etc. and one
      # for postprocessing 
      for hdlr_type in ('main', 'post', 'post_onemesh'):
         dim_handlers = []
         
         for dim in (1, 2, 3):
            # start with mesh suffix list.  If 'onemesh', choose only
            # the first mesh (the last mesh was originally chosen, but higher domain dimensions might finish with lower mesh resolutions)
            if hdlr_type=='post_onemesh':
               if len(mesh_suffix_list_per_dimension[dim-1]) > 0:
                  children = HandlerLevel(
                     'mesh_suffix', [mesh_suffix_list_per_dimension[dim-1][0]])
               else:
                  children = []
            else:
               children = HandlerLevel(
                  'mesh_suffix', mesh_suffix_list_per_dimension[dim-1])
            # if postprocessing, first need to treat fields and norms
            if hdlr_type[0:4]=='post':
               children = field_handler_level.add_sub(
                  norm_handler_level.add_sub(children))
            # if dim > 1, also treat irreg and reg meshes 
            if dim > 1:
               children = mesh_type_handler_level.add_sub(children)
            dim_handlers.append(CompositeHandler('dim', dim, children))

         hdlr = CompositeHandler('stem', case_name, dim_handlers)
         if hdlr_type=='post':
            self.post_handler = hdlr
         elif hdlr_type=='post_onemesh':
            self.post_onemesh_handler = hdlr
         else:
            self.main_handler = hdlr

            
   def write_xml_begin(self):
      self.xml_file = open(self.case_name+'.xml', 'w')
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
  <pass_tests>
    <test name="Solvers converged" language="python">
import os
files = os.listdir("./")
assert(not "matrixdump" in files and not "matrixdump.info" in files)
    </test>""".format(self.case_name, getuser()))

      
   def write_xml_end(self):
      self.xml_file.write("""
  </pass_tests>
  <warn_tests>
  </warn_tests>
</testproblem>""".format(self.case_name, getuser()))
      self.xml_file.close()


   def write_xml(self):
      self.write_xml_begin()
      self.post_onemesh_handler.handle(
         # expect results for the first mesh to be good to within 5%
         WriteXMLSnippet(self.xml_file, 'norm', 0.05,
                         self.solution_dict))
      self.post_handler.handle(
         # expect rates to be better than 0.6
         WriteXMLSnippet(self.xml_file, 'rate', 0.6,
                         self.solution_dict))
      self.write_xml_end()


   def preprocess(self):
      cmds = CommandList([
         ExpandOptionsTemplate(self.finish_time, self.solution_dict),
         ProcessMesh(self.domain_extents)])
      self.main_handler.handle( cmds )


   def process(self):
      self.main_handler.handle( RunSimulation() )


   def postprocess(self):
      with open(error_norms_filename, 'w') as f_norms:
         with open(error_rates_filename, 'w') as f_rates:
            self.post_handler.handle(
               WriteToReport(f_norms, f_rates))


   def clean_up(self):
      self.main_handler.handle( CleanUp() )
      
            
   def do(self, what):
      """Front end for various methods"""
      what = str.lower(what[0:3])
      if what=='xml' or what=='all':
         if verbose(): print '\nGenerating XML file'
         self.write_xml()
      if what=='pre' or what=='all':
         if verbose(): print '\nPreprocessing'
         self.preprocess()
      if what=='pro' or what=='all':
         if verbose(): print '\nRunning simulations'
         self.process()
      if what=='pos' or what=='all':
         if verbose(): print '\nAnalysing'
         self.postprocess()
      if what=='cle':
         if verbose(): print '\nCleaning up'
         self.clean_up()
