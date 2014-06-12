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
from test_tools import Command, CommandList, HandlerLevel, \
   CompositeHandler, LeafHandler, WriteXMLFile, RunSimulation, \
   WriteToReport, error_norms_filename, error_rates_filename, verbose, \
   join_with_underscores
from fluidity_tools import stat_parser
from solution_generator import generate
from importlib import import_module
from vtktools import vtu
from copy import deepcopy

verbose(True)
debug = False

# constants
mesh_template_prefix = 'template'
domain_shape_dict = {1:'line', 2:'rectangle', 3:'cuboid'}
mesh_res_dict = {'A':5, 'B':10, 'C':20, 'D':40,
                 'E':80, 'F':160, 'G':320, 'H':640}


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
         m_t = None
         g_e = '.sh'
         m_e_l = ['.bound', '.ele', '.node']
      else:
         m_t = mesh_type
         g_e = '.geo'
         m_e_l = ['.msh']
      self.mesh_name = join_with_underscores((self.domain_shape, 
                                              m_t, mesh_suffix))
      self.options_template_filename = "template_"+case+".diml"
      self.options_name = join_with_underscores((case, 
         str(dim)+'d', m_t, mesh_suffix))
      self.options_filename = self.options_name+'.diml'
      self.geo_template_filename = join_with_underscores((
         mesh_template_prefix, self.domain_shape, m_t))+g_e
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
         self.outlet_ID = '12'
         self.inlet_ID = '14'
         self.wall_IDs = '11 13'
         self.wall_num = '2'
      elif self.dim==3:
         self.outlet_ID = '32'
         self.inlet_ID = '30'
         self.wall_IDs = '28 29 31 33'
         self.wall_num = '4'


class MMSCommand(Command):
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

         
class GenerateSymbols(MMSCommand):
   
   def __init__(self):
      MMSCommand.__init__(self)
      
   def execute(self, level_name, value, indent):
      # upper levels just record level details; the last level does all
      # the important stuff
      self.store_level_details(level_name, value)
      if level_name == 'stem':
      
         # write dictionary
         solution_name = join_with_underscores((self.stem))
         if verbose():
            div_plus_src = generate(solution_name, check_solution=True)
            sys.stdout.write('\n'+indent+'   Reality check')
            for i in range(3):
               dim = i + 1
               sys.stdout.write(
                  '\n'+indent+'   {0}D: divergence - source = {1}'.\
                     format(dim, str(div_plus_src[i])))
         else:
            generate(solution_name)
      
class ExpandOptionsTemplate(MMSCommand):
   
   def __init__(self, mesh_A_number_of_timesteps=1):
      MMSCommand.__init__(self)
      self.mesh_A_number_of_timesteps = mesh_A_number_of_timesteps
      
   def execute(self, level_name, value, indent):
      # all levels record level details
      self.store_level_details(level_name, value)
      
      if level_name == 'stem':
         # at this level it is appropriate to read in the expressions
         # generated from GenerateSymbols (in this instance, once only)
         solution_name = join_with_underscores((self.stem))
         solution_expressions = import_module(solution_name)
         self.finish_time = solution_expressions.py_dict['finish_time']
         self.text_dict = solution_expressions.text_dict
         
      if level_name == 'mesh_suffix':
         # this last level does most of the important stuff
            
         # compute useful stuff
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
            buf = buf.safe_substitute(self.text_dict)
         except AttributeError:
            pass
         with open(param.options_filename, 'w') as tgt:
            tgt.write(buf)

   def compute_time_step(self):
      """Extracted so as to be overrideable"""
      
      # maintain a constant Courant number for all meshes
      scale_factor = float(mesh_res_dict['A'])/ \
                     float(mesh_res_dict[self.mesh_suffix])
      return scale_factor * self.finish_time / self.mesh_A_number_of_timesteps
      
            
class ProcessMesh(MMSCommand):

   def __init__(self):
      MMSCommand.__init__(self)
      self.binary_path = "../../bin/darcy_impes"
      
   def execute(self, level_name, value, indent):
      # all levels record level details
      self.store_level_details(level_name, value)
      
      if level_name == 'stem':
         # at this level it is appropriate to read in the expressions
         # generated from GenerateSymbols (in this instance, once only)
         solution_name = join_with_underscores((self.stem))
         solution_expressions = import_module(solution_name)

         # needed for e.g. creating 1D meshes, but do not rely on this
         # for true mesh dimensions otherwise
         self.domain_extents = solution_expressions.py_dict['domain_extents']

      if level_name == 'mesh_suffix':
         # this last level does most of the important stuff

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
         # TODO tell user to 'make geo' if needed
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

            
class CleanUp(MMSCommand):
   """Removes generated geometry, mesh, and options files, but leaves
    results intact.  This ensures templates do not get deleted
    (compared with rm -f *.geo, etc)
   """

   def __init__(self):
      MMSCommand.__init__(self)
      
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


class Diagnose(MMSCommand):
   """This class was created to help figure out why the pressure field is
   nonconvergent with mesh resolution.  Analysis is limited to 1D."""

   def __init__(self):
      MMSCommand.__init__(self)
      self.use_analytic_vel_and_src = False
      
   def execute(self, level_name, value, indent):
      # all levels record level details
      self.store_level_details(level_name, value)

      if level_name == 'stem':
         # at this level it is appropriate to read in the expressions
         # generated from GenerateSymbols (in this instance, once only)
         solution_name = join_with_underscores((self.stem))
         solution_expressions = import_module(solution_name)
         py_dict = solution_expressions.py_dict
         self.g = py_dict['gravity_magnitude']
         self.mu = (py_dict['viscosity1'],
                    py_dict['viscosity2'])
         self.rho = (py_dict['density1'],
                     py_dict['density2'])
         self.K = (py_dict['permeability1'],
                   py_dict['permeability2'])
         self.u = (py_dict['darcy_velocity1_x_1D'],
                   py_dict['darcy_velocity2_x_1D'])
         self.q = (py_dict['source_saturation1_1D'],
                   py_dict['source_saturation2_1D'])
         

      if level_name == 'mesh_suffix':
         # the last level does the important stuff

         # 2D or 3D => do nothing
         if self.dim > 1: return

         # compute useful parameters
         param = Parameterisation(self.stem, self.dim, self.mesh_type,
                                  self.mesh_suffix)

         data = vtu(param.options_name+'_1.vtu')

         p = data.GetScalarField('Phase1::Pressure')
         s = (data.GetScalarField('Phase1::Saturation'),
              data.GetScalarField('Phase2::Saturation'))
         x = data.GetLocations()[:,0]
         t = data.GetScalarField('Phase1::Time')
         if self.use_analytic_vel_and_src:
            u = (self.u[0](x, t), self.u[1](x, t))
            q = (self.q[0](x, t), self.q[1](x, t))
         else:
            u = (data.GetVectorField('Phase1::DarcyVelocity')[:, 0],
                 data.GetVectorField('Phase2::DarcyVelocity')[:, 0])
            q = (data.GetScalarField('Phase1::SaturationSource'),
                 data.GetScalarField('Phase2::SaturationSource'))
            
         dx = numpy.diff(x[0:2])
         n = len(x)

         # what does the pressure equation evaluate to, with and without
         # source term?  Evaluate using (i) the pressure and saturation
         # fields and (ii) Darcy velocities.  N.B.  central differencing
         # is introduced in (i), so the two methods may not match
         # exactly.
         grad_p = numpy.gradient(p, dx)
         A, B, C, div_p, div_u, div_minus_src_p, div_minus_src_u, \
             numsrc_minus_anasrc = \
             [[deepcopy([]) for i in range(3)] for j in range(8)]
         rms = lambda v: numpy.sqrt(sum(v**2)/float(len(v)))
         for i in [0, 1]:
            A[i] = self.K[i](s)/self.mu[i]
            B[i] = self.rho[i]*self.K[i](s)/self.mu[i]
            div_p[i] = -(numpy.gradient(A[i]*grad_p - B[i]*self.g, dx))
            div_u[i] = numpy.gradient(u[i], dx)
            div_minus_src_p[i] = div_p[i] - q[i]
            div_minus_src_u[i] = div_u[i] - q[i]
            numsrc_minus_anasrc[i] = \
                data.GetScalarField('Phase{0}::SaturationSource'.\
                                       format(i+1)) - self.q[i](x, t)

         # add phases together.  Sum of div should match
         # Phase1::DivergenceTotalDarcyVelocity.  Sum of div_minus_src
         # should be (near) zero.
         div_p[2] = div_p[0] + div_p[1]
         div_u[2] = div_u[0] + div_u[1]
         div_minus_src_p[2] = div_minus_src_p[0] + div_minus_src_p[1] 
         div_minus_src_u[2] = div_minus_src_u[0] + div_minus_src_u[1] 
         numsrc_minus_anasrc[2] = numsrc_minus_anasrc[0] + numsrc_minus_anasrc[1] 

         val_cf = stat_parser(param.options_name+'.stat')\
             ['Phase1']['DivergenceTotalDarcyVelocity']['l2norm'][-1]
         sys.stdout.write('\n'+indent+'  divergence (c.f. stat file, {0:.4e})'.\
                             format(val_cf))
         sys.stdout.write('\n'+indent+' '*6+'p-based calc '+' '*10+'u-based calc')
         for i in [0, 1, 2]:
            sys.stdout.write('\n'+indent+'  ')
            if i == 2:
               i_str = 'both'
            else:
               i_str = 'ph.'+str(i+1)
            sys.stdout.write('       {0}:  {1:.4e}'.\
                                format(i_str, rms(div_p.pop())))
            sys.stdout.write('       {0}:  {1:.4e}'.\
                                format(i_str, rms(div_u.pop())))

         val_cf = stat_parser(param.options_name+'.stat')\
             ['Phase1']['DivergenceTotalDarcyVelocityCorrected']['l2norm'][-1]
         sys.stdout.write(
            '\n'+indent+
            '  divergence minus source term (c.f. stat file, {0:.4e})'.\
               format(val_cf))
         for i in [0, 1, 2]:
            sys.stdout.write('\n'+indent+'  ')
            if i == 2:
               i_str = 'both'
            else:
               i_str = 'ph.'+str(i+1)
            sys.stdout.write('       {0}:  {1:.4e}'.\
                                format(i_str, rms(div_minus_src_p.pop())))
            sys.stdout.write('       {0}:  {1:.4e}'.\
                                format(i_str, rms(div_minus_src_u.pop())))

         # sys.stdout.write('\n'+indent+'  numerical source minus analytic source:')
         # for i in [0, 1, 2]:
         #    sys.stdout.write('\n'+indent+'  ')
         #    if i == 2:
         #       i_str = 'both'
         #    else:
         #       i_str = 'ph.'+str(i+1)
         #    sys.stdout.write('       {0}:  {1:.4e}'.\
         #                        format(i_str, rms(numsrc_minus_anasrc.pop())))

         # sys.stdout.write('\n'+indent+'  source1 minus source2: {1:.4e}'.\
         #                        format(i_str, rms(q[0] - q[1])))

         sys.stdout.write('\n')

               
class ManufacturedSolutionTestSuite:
   def __init__(self, case_name, mesh_type_list,
                mesh_suffix_list_per_dimension, field_name_list,
                norm_list, mesh_A_number_of_timesteps=1):

      self.case_name = case_name
      
      # build handlers
      mesh_type_handler_level = HandlerLevel('mesh_type', mesh_type_list)
      field_handler_level = HandlerLevel('field', field_name_list)
      norm_handler_level = HandlerLevel('norm', norm_list)

      # maintain a (trivial) tree for generating symbols
      self.shallow_handler = LeafHandler(
         'stem', case_name)

      # maintain one tree for running simulations etc. and a deeper one
      # for configuring and computing norms and convergence rates
      for hdlr_type in ('main', 'deep'):
         dim_handlers = []
         
         for dim in (1, 2, 3):
            # start with mesh suffix list
            children = HandlerLevel(
               'mesh_suffix', mesh_suffix_list_per_dimension[dim-1])
            # if postprocessing, first need to treat fields and norms
            if hdlr_type[0:4]=='deep':
               children = field_handler_level.add_sub(
                  norm_handler_level.add_sub(children))
            # if dim > 1, also treat irreg and reg meshes 
            if dim > 1:
               children = mesh_type_handler_level.add_sub(children)
            dim_handlers.append(CompositeHandler('dim', dim, children))

         hdlr = CompositeHandler('stem', case_name, dim_handlers)
         if hdlr_type=='deep':
            self.deep_handler = hdlr
         else:
            self.main_handler = hdlr

      self.mesh_A_number_of_timesteps = mesh_A_number_of_timesteps


   def generate_symbols(self):
      self.shallow_handler.handle( GenerateSymbols() )


   def write_xml(self):
      self.deep_handler.handle( WriteXMLFile('stem', 0.05, 0.6) )


   def preprocess(self):
      cmds = CommandList([
         ExpandOptionsTemplate(self.mesh_A_number_of_timesteps),
         ProcessMesh()])
      self.main_handler.handle( cmds )


   def process(self):
      self.main_handler.handle( RunSimulation() )


   def postprocess(self):
      with open(error_norms_filename, 'w') as f_norms:
         with open(error_rates_filename, 'w') as f_rates:
            self.deep_handler.handle(
               WriteToReport(f_norms, f_rates))


   def clean_up(self):
      self.main_handler.handle( CleanUp() )


   def diagnose(self):
      self.main_handler.handle( Diagnose() )
      
            
   def do(self, what):
      """Front end for various methods"""
      what = str.lower(what[0:3])
      if what=='gen' or what=='all':
         if verbose(): print '\nGenerating symbols'
         self.generate_symbols()
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
      if what=='dia':
         if verbose(): print '\nStats (RMS)'
         self.diagnose()
