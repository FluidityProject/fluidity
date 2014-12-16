#!/usr/bin/env python

"""Extension of test_tools.py for MMS-based tests associated with the
Darcy IMPES code.

More documentation to come."""

import sys
sys.path.append('../../python')

import os
import string
import subprocess
import re
import numpy
from solution_generator import generate
from importlib import import_module
from copy import deepcopy
from batch_tools import Command, Validate, new_handler as nh
from test_tools import Parameteriser, TestCommand, TestSuite, \
    AssertThresholdCalculator

simulator_name = 'darcy_impes'
simulation_options_filename_extension = '.diml'


class ManufacturedSolutionParameteriser(Parameteriser):
   def __init__(self, simulation_options_filename_extension,
                reference_mesh_res, reference_time_step_number):
      # don't bother initialising the superclass with domain size and
      # finish time; we will override its getters via generated
      # solution expressions.
      Parameteriser.__init__(
         self, simulation_options_filename_extension,
         reference_mesh_res=reference_mesh_res,
         reference_time_step_number=reference_time_step_number)

   def domain_extents(self):
      return self.lookup('solution', 'domain_extents')
   
   def finish_time(self):
      return self.lookup('solution', 'finish_time')

   def geometry_dict(self):
      # initialise result
      result = { 'MESH_NAME': self.lookup('filename_stems', 'mesh') }
      # make entries for each dimension
      for i, dim_str in enumerate(('X', 'Y', 'Z')):
         L = self.domain_extents()[i]
         n = self.element_number(L)
         result.update({
               'DOMAIN_LENGTH_{0}'.format(dim_str) : str(L),
               'EL_NUM_{0}'.format(dim_str) : str(n),
               'EL_SIZE_{0}'.format(dim_str) : str(L/n) })
      return result

   def simulation_options_dict(self):
      dim = self.get('dim')
      # some dictionary entries depend on whether 1D or multi-D
      if dim=="1":
         mesh_format = 'triangle'
         # these markers are for collapsing inappropriate wall BCs
         begin_rm_if_1D = '<!--'
         end_rm_if_1D = '-->'
      else:
         mesh_format = 'gmsh'
         begin_rm_if_1D = ''
         end_rm_if_1D = ''
      # these boundary IDs are dependent on the geometry templates and
      # are different for each dimensional space
      if dim=='1':
         outlet_ID = '2'
         inlet_ID = '1'
         wall_IDs = ''
         wall_num = '0'
      elif dim=='2':
         outlet_ID = '12'
         inlet_ID = '14'
         wall_IDs = '11 13'
         wall_num = '2'
      elif dim=='3':
         outlet_ID = '32'
         inlet_ID = '30'
         wall_IDs = '28 29 31 33'
         wall_num = '4'
      # now compile the dictionary 
      return {
         'MESH_NAME': self.lookup('filename_stems', 'mesh'),
         'MESH_FORMAT': mesh_format,
         'DOMAIN_DIM': dim,
         'TIME_STEP': str(self.time_step()),
         'FINISH_TIME': str(self.finish_time()),
         'INLET_ID': inlet_ID,
         'OUTLET_ID': outlet_ID,
         'WALL_IDS': wall_IDs,
         'WALL_NUM': wall_num,
         'BEGIN_RM_IF_1D': begin_rm_if_1D,
         'END_RM_IF_1D': end_rm_if_1D }
   
      
class GenerateSymbols(TestCommand):
    def __init__(self, handler, parameteriser):
        TestCommand.__init__(self, ['case'], handler,
                            "Generating symbols", parameteriser)
        
    def execute(self):
        # the last level does everything.  Exit if we are not there yet
        if not self.at('leaf'):
            return
            
        # write dictionary
        solution_name = self.make_key()
        if self.is_verbose():
            div_plus_src = generate(solution_name, check_solution=True)
            self.write_message(
               '   Sanity check; expect sum(div(u) - src) = 0...',
               with_newline=True)
            for i in range(3):
                dim = i + 1
                self.write_message('      {0}D: {1}'.\
                                      format(dim, str(div_plus_src[i])),
                                   with_newline=True)
        else:
            generate(solution_name)


class TestCommandDecorator(Command):
    """Wraps a TestCommand, giving it access to the generated solution"""
    def __init__(self, wrapped_test_command, solution_level_names):
       Command.__init__(self, verbosity_threshold=1)
       # make sure the solution level names exist in the wrapped
       # object's handler
       cmd_name = 'TestCommandDecorator\'s wrapped_test_command'
       val = Validate(solution_level_names, cmd_name)
       if val.is_verbose():
            msg = 'Validating '+cmd_name
       else:
            msg = ''
       wrapped_test_command.handle(val, msg)
       self.__wrapped_test_command = wrapped_test_command
       self.__solution_level_names = solution_level_names
       
    def inform(self, level_name, node_name, at_leaf, indent, verbose):
       """Delegation."""
       self.__wrapped_test_command.inform(level_name, node_name,
                                        at_leaf, indent, verbose)
       
    def execute(self):
        """Imports generated expressions, if appropriate, before delegating
        to wrapped TestCommand"""
        # when this is the last of the solution levels, import generated
        # expressions
        if self.__wrapped_test_command.get('level') == \
               self.__solution_level_names[-1]:
           solution_name = self.__wrapped_test_command.make_key(
              self.__solution_level_names)
           self.__solution_expressions = import_module(solution_name)
        
        # when this is the very last level, update the solution- and
        # level_dependent dictionaries
        if self.__wrapped_test_command.at('leaf'):
           self.__wrapped_test_command.register_dict(
              'solution', self.__solution_expressions.py_dict)
           self.__wrapped_test_command.update_dict(
              'simulation_options', self.__solution_expressions.text_dict)
        
        # continue to wrapped execute()
        self.__wrapped_test_command.execute()

        
    def finalise(self):
        """Delegates to wrapped TestCommand"""
        self.__wrapped_test_command.finalise()

        
    def handle(self, other=None, message=None):
       """Exploits the handler of the wrapped class."""
       self.__wrapped_test_command.handle(self)

        
               
class ManufacturedSolutionTestSuite(TestSuite):
    """A front end for running various MMS-tailored
    TestCommands."""
   
    def __init__(self, case_name, mesh_type_list, mesh_res_lists,
                 field_list, norm_list, extra_levels=None, 
                 reference_mesh_res=5, reference_time_step_number=1,
                 command_line_in_xml=None,
                 norm_threshold=0.05, rate_threshold=0.8,
                 simulator_verbosity=0,
                 python_verbosity=1):
        """In contrast to the inherited TestSuite, the present class builds its
        own handlers from lists supplied by the client.  The keyword
        args at the end can be used to make extra handlers at the
        solution level.  For example one might want to iterate over
        different values of pressure and saturation for a given test
        case.  Then one might specify: 'pressure', ['1000', '2000',
        '3000'], 'saturation', ['0.1', '0.2', '0.3'].

        """

        # make a handler for iteration over parameterised solutions
        self.__solution_handler = nh('case', case_name)
        try:
           for level_name, node_names in extra_levels:
              self.__solution_handler = self.__solution_handler * \
                  nh(level_name, node_names)
        except:
           pass

        # extract a corresponding list of level names.  Continuing the
        # example given in the class description, one might end up with
        # ['case', 'pressure', 'saturation'].
        self.__solution_level_names = self.__solution_handler.get_level_names()

        # make another handler for running simulations etc.
        simulation_handler = deepcopy(self.__solution_handler) * \
            nh('dim',
               [nh('1') * \
                    nh('mesh_res', mesh_res_lists[0]),
                nh('2') * \
                    nh('mesh_type', mesh_type_list) * \
                    nh('mesh_res', mesh_res_lists[1]),
                nh('3') * \
                    nh('mesh_type', mesh_type_list) * \
                    nh('mesh_res', mesh_res_lists[2])])

        # ideally the mesh handler wouldn't be dependent on the solution
        # handler, but in this module the solution dictionary contains
        # variables which are linked to the geometry dictionary.  So the
        # mesh handler will have to be as big as the simulation handler.
        mesh_handler = simulation_handler
                
        # and a handler for configuring and computing norms and rates.
        # Note that mesh resolution always comes last.
        results_handler = deepcopy(self.__solution_handler) * \
            nh('dim', [nh('1') * \
                           nh('field', field_list) * \
                           nh('norm', norm_list) * \
                           nh('mesh_res', mesh_res_lists[0]),
                       nh('2') * \
                           nh('mesh_type', mesh_type_list) * \
                           nh('field', field_list) * \
                           nh('norm', norm_list) * \
                           nh('mesh_res', mesh_res_lists[1]),
                       nh('3') * \
                           nh('mesh_type', mesh_type_list) * \
                           nh('field', field_list) * \
                           nh('norm', norm_list) * \
                           nh('mesh_res', mesh_res_lists[2])])

        # create an object for MMS-tailored parameterisations
        self.__parameteriser = ManufacturedSolutionParameteriser(
           simulation_options_filename_extension, reference_mesh_res, 
           reference_time_step_number)
                       
        # using the above objects, initialise the superclass
        TestSuite.__init__(
           self, mesh_handler, simulation_handler, results_handler,
           parameteriser = self.__parameteriser,
           command_line_in_xml = command_line_in_xml,
           assert_threshold_calculator = AssertThresholdCalculator(
              norm_threshold, rate_threshold),
           simulator_name = simulator_name,
           simulation_options_filename_extension = \
              simulation_options_filename_extension,
           simulator_verbosity = simulator_verbosity,
           python_verbosity = python_verbosity)
        
        
    def create_generate_symbols_command(self):
       """New implementation."""
       return GenerateSymbols(self.__solution_handler, self.__parameteriser)

    def create_write_xml_file_command(self):
       """Method override.  Decorates the object created by the original
       TestSuite.
       """
       return TestCommandDecorator(
          TestSuite.create_write_xml_file_command(self),
          self.__solution_level_names)

    def create_mesh_command(self):
       """Method override.  Decorates the object created by the original
       TestSuite.
       """
       return TestCommandDecorator(
          TestSuite.create_mesh_command(self),
          self.__solution_level_names)

    def create_expand_options_command(self):
       """Method override.  Decorates the object created by the original
       TestSuite.
       """
       return TestCommandDecorator(
          TestSuite.create_expand_options_command(self),
          self.__solution_level_names)
