#!/usr/bin/env python

"""buckley_leverett_test_tools.py

For interpolating between simulated fields and quasi-analytic solutions,
reading/writing reports of error norms and convergence rates, and batch
processing the darcy_impes tests.

About the interpolations: it is straightforward to compare the the 1D
simulations to the expected 1D quasi-analytical solution at the latter's
point locations; a piecewise-linear interpolant is constructed between
the simulation's points (see compute_abs_errors).  This is difficult to
extend to 2D and 3D, but not if the interpolation is done the other way
round - the analytical solution is used to interpolate at each
simulation point (interp_using_analytic).  The jump discontinuity of the
analytical solution will be ill-represented, but we can still expect
first order convergence of the l1-errors.  However, l2-errors are
deprecated.  Their reported convergence rates can be seen to be much
less than their l1 counterparts, presumably as a result of the squaring
of errors at the discontinuity.

"""

import sys
sys.path.append('../../python')

import fluidity_tools
import vtktools
import numpy
import re
import subprocess
import os
from batch_tools import new_handler as nh
from test_tools import Parameteriser, TestCommand, TestSuite, \
    AssertThresholdCalculator

error_norms_1D_filename = "error_norms_1d.txt"

# set test_tools verbosity
verbosity(1)
debug = False


## HELPER FUNCTIONS/CLASSES
    
def read_analytic_solution(filename):
    """Converts a file containing a column of x-coordinates and a column of
    values into corresponding lists."""
    x = []; v = [];
    with open(filename, "r") as f:
       f.seek(0)
       for l in f:
          if not l.strip(): continue
          x.append(float(l.split()[0]))
          v.append(float(l.split()[1]))
    if not numpy.all(numpy.diff(x) > 0): 
        print ("Issue with analytic mesh: "
               "x-coordinates not monotonically increasing.")
        sys.exit()
    return x, v

    
def read_numerical_solution(filename, field_name):
    """Converts a vtu file into a list of x-coordinates and field values."""
    f = vtktools.vtu(filename)
    x = f.GetLocations()[:,0]
    v = f.GetScalarField(field_name)
    isort = sorted(range(len(x)), key=lambda i: x[i])
    x = numpy.array([x[i] for i in isort])
    v = numpy.array([v[i] for i in isort])
    if not numpy.all(numpy.diff(x) > 0): 
        print ("Issue with numerical mesh: "
               "x-coordinates not monotonically increasing.")
        sys.exit()
    return x, v


class BLWriteToReport(WriteToReport):

    def __init__(self, handler, parameteriser):
        WriteToReport.__init__(handler, parameteriser,
                               norm_calculator=BLNormCalculator())

    def __del__(self):
        WriteToReport.__del__(self)
        try:
            self.norms_file_1D.close()
        except AttributeError:
            pass

    def execute(self):
        if self.at('root'):
            # if this is the very first level, create the file objects
            self.norms_file = open(_error_norms_filename, 'w')
            self.rates_file = open(_error_rates_filename, 'w')

                
class WriteErrorNorm1D(TestCommand):
    """Compares a 1D numerical solution with the appropriate analytical
    solution, writing the error norm to a report.  The output report
    has two columns: the first has identifiers in the format
    1d_mesh_res_field_ID_normtype; the second has the values
    of the error norms.
    """

    def __init__(self, analytic_stem, norms_file):
        self.analytic_stem = analytic_stem
        self.norms_file = norms_file

    def execute(self):
        if self.at('field'):
            # store the errors and dx for this field.  In the ID, change
            # the dimension to 1 to signify the 1D solution
            num_suf = self.make_key(exclude=['case', 'field'],
                                    substitute={'dim': '1'})
            [self.errs, self.dx] = self.compute_abs_errors(
                num_suf, self.get('field'), self.get['field_long'])
        
        if self.at('leaf'):
            # print the ID and value.  N.B. key is different to num_suf
            # passed to compute_abs_errors.
            key = '{0}_1d_{1}_l{2:g}_{3}'.format(
                self.model, self.field_short, self.norm,
                self.mesh_suffix)
            v = self.compute_norm(self.errs, self.dx)
            self.norms_file.write(
                '{0}{1:24.12e}\n'.format(key,v))
            
    def compute_abs_errors(self, numerical_suffix, analytic_suffix, 
                           field_name):
        """Computes errors between interpolated numerical values and analytic
        point values.  Also returns mean mesh spacing."""
        s_a = self.analytic_stem+'_'+analytic_suffix+'.txt'
        s_n = self.stem+'_'+numerical_suffix+'_1.vtu'
        [x_a, v_a] = read_analytic_solution(s_a)
        [x_n, v_n] = read_numerical_solution(s_n, field_name)
        v_nInterp = numpy.interp(x_a, x_n, v_n)
        return numpy.abs(v_nInterp - v_a), (x_a[-1] - x_a[0])/(len(x_a) - 1)

    def compute_norm(self, errors, mesh_spacing):
        """Computes a norm.  Uniform, linear elements are assumed, so the
        integral over the domain can be approximated by the trapezium rule
        (for the L1 norm at least).
        """
        E = numpy.array(errors)
        if self.norm == 1: 
            I = (E[0]/2 + sum(E[1:-1]) + E[-1]/2)*mesh_spacing
        elif self.norm == 2:
            # does not make much difference, but previously I = (E[0]**2/2 +
            # sum(E[1:-1]**2) + E[-1]**2/2)*mesh_spacing
            I = (E[0]**2 + sum(E[:-1]*E[1:]) + 2*sum(E[1:-1]**2) + \
                     E[-1]**2)*mesh_spacing/3
        else: 
            print("Only l1 or l2 norms may be computed; i.e. self.norm "
                  "must be 1 or 2.")
            sys.exit()
        return I**(1/float(self.norm))
   

# class BLWriteToReport(WriteToReport):
#     """Overrides WriteToReport and its get_norm method
#     in order to allow customised computation of 1D errors.
#     """

#     def __init__(self, handler):
#         WriteToReport.__init__(self, handler)

#     def get_norm(self):
#         if self.tracker.get('dim')==1:
            
#             key = '{0}_1d_{1}_l{2:g}_{3}'.format(
#                 self.model, self.field_short, self.norm,
#                 self.mesh_suffix)
#             v = find(error_norms_1D_filename, key)
#         else:
#             filename = self.stem+'_'+self.model+'_'+\
#                        str(self.dim)+'d_'+self.mesh_suffix+'.stat'
#             if self.norm==1:
#                 self.calc_type = "integral"
#             elif self.norm==2:
#                 self.calc_type = "l2norm"
#             v = fluidity_tools.stat_parser(filename)\
#                 ['Phase'+self.phase_index]\
#                 ['Analytic'+self.var_name+'Error']\
#                 [self.calc_type][-1]
#         return v

    
class BLNormCalculator(NormCalculator):
    def get_requisite_levels(self):
        return ['field', 'mesh_res']
    
    def calculate_norm(self, test_command):
        """There are two ways of computing the norm; see module
        docstring."""
        if test_command.get('dim')==1:
            key = test_command.make_key()
            return find(error_norms_1D_filename, key)
        else:
            return NormCalculator.calculate_norm(self, test_command)
                
            
## FOR PUBLIC USE 

def interp_using_analytic(x, analytic_filename):
    """Interpolates at point x using the profile given by
    analytic_filename."""
    [x_a, v_a] = read_analytic_solution(analytic_filename)
    return numpy.interp(x, x_a, v_a)

       
class BuckleyLeverettTestSuite:
    def __init__(self, case_name, mesh_type_list, mesh_res_lists,
                 field_list, analytic_filename_stem, extra_levels=None, 
                 reference_mesh_res=5, reference_time_step_number=1,
                 command_line_in_xml=None, norm_threshold=0.05,
                 rate_threshold=0.8, simulator_verbosity=0,
                 python_verbosity=1):
        """analytic_filename_stem is a string to be appended with
        shorthand names for the fields in question, which in turn are
        translated from the elements in field_name_list.  If a
        field_name is 'Phase2::Saturation', then the shorthand is
        saturation2."""
        
        self.analytic_filename_stem = analytic_filename_stem

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
        
        # make a handler for meshing, simulations etc
        dim_handlers = []
        for i, dim in ('1', '2', '3'):
            dim_handlers.append(
                nh(dim) * nh('mesh_res', mesh_res_lists[i]))
        sim_handler = deepcopy(shallow_handler) * dim_handlers
        
        # and a handler for configuring and computing norms and rates.
        # Note that mesh resolution comes last.
        norm_list = ['1', '2']
        dim_handlers = []
        for i, dim in ('1', '2', '3'):
            dim_handlers.append(
                nh(dim) * \
                    nh('field', field_list) * \
                    nh('norm', norm_list) * \
                    nh('mesh_res', mesh_res_lists[i]))
        results_handler = deepcopy(shallow_handler) * dim_handlers
                
        # the computation of 1D error norms is a special case where mesh
        # resolution is handled before fields and norms.
        results_1D_handler = \
            deepcopy(shallow_handler) * \
            nh('mesh_res', mesh_res_lists[0]) * \
            nh('field', field_list) * \
            nh('norm', norm_list)

        # create an object for BL-tailored parameterisations
        self.__parameteriser = BuckleyLeverettParameteriser(
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

        self.preproc_cmd = 
        self.sim_handler.handle(CommandList((Preprocess(sim_handler),
                                             RunSimulation(sim_handler),
                                             Preprocess(sim_handler, clean=True))))
 
        
    def generate_error_report_1D(self):
        if verbose(): print '\nComputing 1D error norms'
        # open write file
        with open(error_norms_1D_filename, "w") as f_norms:
            # build handler structure
            handler = self.new_handler(
                self.model_handler_level.add_sub(
                    # 1D, so just add_sub mesh_suffix_handler_levels[0]
                    self.mesh_suffix_handler_levels[0].add_sub(
                        self.field_handler_level.add_sub(
                            self.norm_handler_level))))
            # pass it a command 
            handler.handle(
                WriteErrorNorm1D(self.analytic_filename_stem, f_norms))
    
    
    def generate_convergence_report(self):
        if verbose(): print '\nComputing convergence rates'
        # open write file(s)
        with open(error_rates_filename, "w") as f_rates:
            
            # build handler structure
            dim_handlers = []
            for i, dim in enumerate((1, 2, 3)):
                dim_handlers.append(CompositeHandler('dim', dim,
                        self.field_handler_level.add_sub(
                            self.norm_handler_level.add_sub(
                                self.mesh_suffix_handler_levels[i]))))
            handler = self.new_handler(
                self.model_handler_level.add_sub(dim_handlers))

            # also write ND norms?
            if debug:
                f_norms = open(error_norms_filename, "w")
            else:
                f_norms = None
                
            # pass the handler a command
            try:
                handler.handle(
                    BLWriteToReport(f_norms, f_rates))
            finally:
                if debug: f_norms.close()
    
    def generate_reports(self):
        self.generate_error_report_1D()
        self.generate_convergence_report()

    def create_write_report_command(self):
        return WriteToReport(self.results_handler,
                             self.parameteriser,
                             self.norm_calculator)

                              
