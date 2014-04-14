#!/usr/bin/env python

"""manufactured_solution_test_tools.py

"""

import sys
sys.path.append('../../python')

import os
import string
import subprocess
import re
import numpy
import glob
# from batch_tools import Command, CommandList, HandlerList, CompositeHandler, verbose
from test_tools import Command, CommandList, HandlerList, CompositeHandler, RunSimulation, WriteToReport, verbose

error_rates_filename = "error_rates.txt"
error_norms_filename = "error_norms.txt"

debug = True
verbose(True)

class Parameterisation:
   """Parameters determined from a small, fixed set of
   variables can be determined here.
   """

   def __init__(self, case, dim, mesh_type, mesh_suffix):
      # constants
      self.options_template_filename = "template_"+case+".diml"
      self.mesh_template_prefix = 'template'
      self.domain_shape_dict = {1:'line', 2:'rectangle', 3:'cuboid'}
      self.mesh_res_dict = {'A':5, 'B':10, 'C':20, 'D':40,
                            'E':80, 'F':160}
      
      # basic parameterisations
      self.dim = dim
      self.domain_shape = self.domain_shape_dict[dim]
      self.mesh_res = self.mesh_res_dict[mesh_suffix]

      # filenames
      if dim==1:
         mt = ''
         # me = '.ele'
      else:
         mt = '_'+mesh_type
         # me = '.msh'
      self.mesh_name = self.domain_shape+mt+'_'+mesh_suffix
      self.options_name = case+'_'+str(dim)+'d'+mt+'_'+mesh_suffix
      self.geo_template_filename = self.mesh_template_prefix+'_'+\
                                   self.domain_shape+mt+'.geo'
      self.geo_filename = self.mesh_name+'.geo'
      # self.mesh_filename = self.mesh_name+me
      self.options_filename = self.options_name+'.diml'
      
   def compute_more_options(self):
      """Options template expansion typically requires more stuff.
      It is placed here so as not to burden other clients.
      """

      self.timestep_scale_factor = self.mesh_res/self.mesh_res_dict['A']
      
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
         self.gravity_direction = '-1'
         self.outlet_ID = '2'
         self.inlet_ID = '1'
         self.wall_IDs = ''
         self.wall_num = '0'
      elif self.dim==2:
         self.gravity_direction = '-1 0'
         self.outlet_ID = '11'
         self.inlet_ID = '13'
         self.wall_IDs = '12 14'
         self.wall_num = '2'
      elif self.dim==3:
         self.gravity_direction = '-1 0 0'
         self.outlet_ID = '28'
         self.inlet_ID = '29'
         self.wall_IDs = '30 31 32 33'
         self.wall_num = '4'

      
class ExpandOptionsTemplate(Command):
   
   def __init__(self, finish_time, solution_dict=None):
      self.finish_time = finish_time
      # solution_dict may come from an external file generated as part of a
      # manufactured solution.
      self.solution_dict = solution_dict
      self.stem = None
      self.dim = None
      self.mesh_type = None
      
   def execute(self, level_name, value, indent):
      # upper levels just record level names; the last level does all
      # the important stuff
      if level_name == 'stem':
         self.stem = value
      if level_name == 'dim':
         self.dim = value
      if level_name == 'mesh_type':
         self.mesh_type = value
      if level_name == 'mesh_suffix':
         mesh_suffix = value
         # compute useful stuff
         param = Parameterisation(self.stem, self.dim,
                                  self.mesh_type, mesh_suffix)
         param.compute_more_options()
         
         # let the coarsest mesh have one time step and maintain a
         # constant Courant number for all meshes
         time_step = self.finish_time * param.timestep_scale_factor

         # compose the dictionary
         options_dict = {
            'MESH_FILE': param.mesh_name,
            'MESH_FORMAT': param.mesh_format,
            'DOMAIN_DIM': str(self.dim), 
            'TIME_STEP': str(time_step),
            'FINISH_TIME': self.finish_time,
            'GRAVITY_DIRECTION': param.gravity_direction,
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

         

            
class ProcessMesh(Command):

   def __init__(self, domain_length_1D):
      self.binary_path = "../../bin/darcy_impes"
      self.stem = None
      self.dim = None
      self.mesh_type = None
      self.domain_length_1D = domain_length_1D
      
   def execute(self, level_name, value, indent):
      if level_name == 'stem':
         self.stem = value
      if level_name == 'dim':
         self.dim = value
      if level_name == 'mesh_type':
         self.mesh_type = value
      if level_name == 'mesh_suffix':
         mesh_suffix = value

         # compute useful stuff
         param = Parameterisation(self.stem, self.dim,
                                  self.mesh_type, mesh_suffix)

         # expand mesh template
         geo_dict = {
            'EL_NUM': str(param.mesh_res),
            'MESH_NAME': str(param.mesh_name) }
         with open(param.geo_template_filename) as src:
            buf = string.Template( src.read() )
         buf = buf.safe_substitute(geo_dict)
         with open(param.geo_filename, 'w') as tgt:
            tgt.write(buf)

         # call interval or gmsh
         if self.dim==1:
            lx=self.domain_length_1D
            dx=lx/param.mesh_res
            subprocess.call(['../../bin/interval', '0.0',
                             str(lx), '--dx='+str(dx),
                             param.mesh_name])
         else:
            subprocess.call(['gmsh', '-'+str(self.dim),
                             param.geo_filename],
                            stdout=open(os.devnull, 'wb'))

        
# class Analyse(Command):

#    def __init__(self):
#       self.stem = None
#       self.dim = None
#       self.mesh_type = None

#    def execute(self, level_name, value, indent):
#       if level_name == 'stem':
#          self.stem = value
#       if level_name == 'dim':
#          self.dim = value
#       if level_name == 'mesh_type':
#          self.mesh_type = value
#       if level_name == 'mesh_suffix':
#          mesh_suffix = value

      
class ManufacturedSolutionTestSuite:
   def __init__(self, case_name, solution_dict, finish_time,
                domain_length_1D, mesh_type_list,
                mesh_suffix_list_per_dimension,
                field_name_list):

      self.case_name = case_name
      self.solution_dict = solution_dict
      self.finish_time = finish_time

      # domain_length_1D is only needed to create 1D meshes;
      # it can be initialised with a dummy value otherwise
      self.domain_length_1D = domain_length_1D
      
      # build handlers
      mesh_type_handler_list = HandlerList('mesh_type', mesh_type_list)
      field_handler_list = HandlerList('field', field_name_list)
      norm_handler_list = HandlerList('norm', (1, 2))
      # maintain one tree for running simulations etc. and one
      # for postprocessing 
      for hdlr_type in ('main', 'post'):
         dim_handlers = []
         for dim in (1, 2, 3):
            # start with mesh suffix list
            children = HandlerList(
               'mesh_suffix', mesh_suffix_list_per_dimension[dim-1])

            # if postprocessing, first need to treat fields and norms
            if hdlr_type=='post':
               children = field_handler_list.expand(
                  norm_handler_list.expand(children))

            # if dim > 1, also treat irreg and reg meshes 
            if dim > 1:
               children = mesh_type_handler_list.expand(children)
               
            dim_handlers.append(CompositeHandler('dim', dim, children))

         hdlr = CompositeHandler('stem', case_name, dim_handlers)
         if hdlr_type=='post':
            self.post_handler = hdlr
         else:
            self.main_handler = hdlr
            
   def do(self, what):
      what = str.lower(what[0:3])
      
      if what=='pre' or what=='all':
         if verbose(): print '\nPreprocessing'
         cmds = CommandList([
            ExpandOptionsTemplate(self.finish_time, self.solution_dict),
            ProcessMesh(self.domain_length_1D)])
         self.main_handler.handle( cmds )
         
      if what=='pro' or what=='run' or what=='all':
         if verbose(): print '\nRunning simulations'
         self.main_handler.handle( RunSimulation() )
         
      if what=='pos' or what=='ana' or what=='all':
         if verbose(): print '\nAnalysing'
         with open(error_norms_filename, 'w') as f_norms:
            with open(error_rates_filename, 'w') as f_rates:
               self.post_handler.handle(
                  WriteToReport(f_norms, f_rates))


   # def analyse(self):
   #    print '\n{0:s}::{1:s} at ElapsedTime={2:g}'.\
   #       format(self.options.phase_name, self.options.var_name, \
   #              self.options.end_time)
   #    for domain_dim in self.domain_dim_list:
   #       print '   '+str(domain_dim)+'d:'
   #       for domain_size in self.domain_size_list:
   #          print '      scale='+domain_size+':'
   #          for mesh_type in self.mesh_type_list:
   #             print '         mesh='+str(mesh_type)+':'
   #             self.analyse_mesh_series(str(domain_dim)+'D_'+\
   #                                    domain_size+'_'+\
   #                                    str(mesh_type))

   # def clean(self):
   #    print '\nCleaning'
   #    for domain_dim in self.domain_dim_list:
   #       print '   '+str(domain_dim)+'D:'
   #       for domain_size in self.domain_size_list:
   #          print '      scale='+domain_size+':'
   #          for mesh_type in self.mesh_type_list:
   #             print '         mesh='+str(mesh_type)+':'
   #             for mesh_res in self.mesh_res_list:
   #                print '            res='+str(mesh_res)
   #                s = Sim(domain_dim, domain_size, mesh_type,
   #                        mesh_res, self.options)
   #                s.clean()

   def analyse_mesh_series(self, filename_stem):
      sys.path.append(self.options.fluidity_path + "python/")
      from fluidity_tools import stat_parser as stat

      if self.options.read_errors:
         metric_name = 'l2_norm'
      else:
         metric_name = 'max_abs'
     
      n_mesh = len(self.mesh_res_list)
      metric = range(n_mesh)
      delta = range(n_mesh-1)
      order = range(n_mesh-1)

      # for each mesh, take samples within the specified time window.
      # Can use a nearest neighbour lookup.  Note that the time abscissa
      # may vary from mesh to mesh if e.g. adaptive timestepping is
      # used.  
      first = True
      for i_mesh, mesh_res in enumerate(self.mesh_res_list):

         # check stat file exists
         fname = filename_stem+'_'+self.options.mesh_suffix_dict[mesh_res]+\
                 '.stat'
         if not os.path.isfile(fname): continue
         
         # get the time history of the metric of interest
         if self.options.read_errors:
            vals = numpy.abs(stat(fname)[self.options.phase_name]\
                             [self.options.var_name]['l2norm'])
         else:
            min_vals = numpy.abs(stat(fname)[self.options.phase_name]\
                                 [self.options.var_name]['min'])
            max_vals = numpy.abs(stat(fname)[self.options.phase_name]\
                                 [self.options.var_name]['max'])
            vals = numpy.amax(numpy.abs(numpy.vstack((min_vals, max_vals))), axis=0)

         # look up the proper end index by end time
         # (for some reason the last value is sometimes spurious)
         t = stat(fname)['ElapsedTime']['value']
         ii = range(len(t))
         i1 = int(numpy.interp(self.options.end_time, t, ii))
         metric[i_mesh] = vals[i1]

         # print stuff
         if first:
            print 12*' '+'res={0:3g}: {1:s}={2:11.4e}'.\
               format(mesh_res, metric_name, metric[i_mesh])
            first = False
         elif self.options.read_errors:
            arg1 = metric[i_mesh-1]/metric[i_mesh]
            if arg1 == 0.:
               order[i_mesh-1] = numpy.nan
               order_str = 'N/A'
            else:
               arg2 = self.mesh_res_list[i_mesh] \
                      /self.mesh_res_list[i_mesh-1]
               order[i_mesh-1] = numpy.log(arg1)/numpy.log(arg2)
               order_str = '{0:7.4f}'.format(order[i_mesh-1])
               
            print 12*' '+'res={0:3g}: {1:s}={2:11.4e}, order={3:s}'.\
               format(mesh_res, metric_name, metric[i_mesh], \
                      order_str)
         else:
            delta[i_mesh-1] = metric[i_mesh] - metric[i_mesh-1]
            print 12*' '+'res={0:3g}: {1:s}={2:11.4e}, delta={3:11.4e}'.\
               format(mesh_res, metric_name, metric[i_mesh], \
                      delta[i_mesh-1])

         # write more stuff to files if requested
         if self.options.summarise_history:
            dt = stat(fname)['dt']['value']
            fname = filename_stem+'_'+\
                    self.options.mesh_suffix_dict[mesh_res]+'_'+\
                    self.options.var_name+'_'+\
                    self.options.phase_name+'.hist'
            with open(fname, 'w') as f:
               f.write('t' + 10*' ' + 'dt' + 10*' ' + \
                       metric_name + '\n')
               for i in range(len(t)):
                  f.write('{0:8.4e} {1:8.4e} {2:12.4e}\n'.\
                          format(t[i], dt[i], vals[i]))
               
               
