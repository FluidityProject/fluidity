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
from test_tools import Command, HandlerList, CompositeHandler, set_verbose
import pdb

verbose = True
set_verbose(verbose)

error_norms_filename = "error_norms.txt"
error_rates_filename = "error_rates.txt"


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
    if not numpy.all(numpy.diff(x) > 0): 
        print ("Issue with numerical mesh: "
               "x-coordinates not monotonically increasing.")
        sys.exit()
    return x, v
   
   
def find(report, key):
    """Searches a report and returns the value associated with key
    (ID). """
    report.seek(0)
    v = 0.
    for l in report:
        if not l.strip(): continue
        if l.split()[0] == key: 
            v = float(l.split()[1])
            break
    return v

        
class RunSimulation(Command):
    """Preprocesses an input file before running the respective
    simulation."""
        
    def __init__(self):
        self.binary_path = "../../bin/darcy_impes"

    def execute(self, level_name, value):
        # most levels just record level names
        if level_name == 'stem':
            self.stem = value
        if level_name == 'model':
            self.model = value
        if level_name == 'dim':
            self.dim = value
        if level_name == 'grid':
            # last level
            grid = value
            
            # for the most coarse grid, an options file should exist.
            # For the other grids, adapt it
            case = self.stem+'_'+self.model+'_'+self.dim
            filename_new = case+'_'+grid+'.diml'
            if grid!='A':
                filename_A = case+'_A.diml'
                with open(filename_A, 'r') as f_A:
                    f = open(filename_new, 'w')
                    for line in f_A:
                        f.write(line.replace('_A', '_'+grid))
                    f.close()

            # process and clean up 
            if not os.path.isfile(self.binary_path):
                raise IOError("The darcy_impes binary hasn't been built yet.")
            subprocess.call([self.binary_path, filename_new])
            if grid!='A': os.remove(filename_new)
        
                
class WriteErrorNorm1D(Command):
    """Compares a 1D numerical solution with the appropriate analytical
    solution, writing the error norm to a report.  The output report
    has two columns: the first has identifiers in the format
    model_name_1D_grid_name_field_ID_normtype; the second has the values
    of the error norms.
    """

    def __init__(self, analytic_stem, norms_file):
        self.analytic_stem = analytic_stem
        self.norms_file = norms_file

    def execute(self, level_name, value):
        if level_name == 'stem':
            self.stem = value
        if level_name == 'model':
            self.model = value
        if level_name == 'grid':
            self.suffix = self.model+'_1d_'+value
        if level_name == 'field':
            # translate the field_name to shorthand form
            field = value
            self.field_short = str.lower(re.sub('Phase(.)::(.*)', '\\2\\1',
                                                field))
            # store the errors and dx for this field
            [self.errs, self.dx] = self.compute_abs_errors(
                self.suffix, self.field_short, field)
            
        if level_name == 'norm':
            self.norm = value
            
            # print the ID and value
            key = '{0}_{1}_l{2:g}'.format(
                self.suffix, self.field_short, self.norm)
            v = self.compute_norm(self.errs, self.dx)
            self.norms_file.write(
                '{0}{1:24.12e}\n'.format(key,v))
            
    def compute_abs_errors(self, numerical_suffix, analytic_suffix, 
                           field_name):
        """Computes errors between interpolated numerical values and analytic
        point values.  Also returns mean grid spacing."""
        s_a = self.analytic_stem+'_'+analytic_suffix+'.txt'
        s_n = self.stem+'_'+numerical_suffix+'_1.vtu'
        [x_a, v_a] = read_analytic_solution(s_a)
        [x_n, v_n] = read_numerical_solution(s_n, field_name)
        v_nInterp = numpy.interp(x_a, x_n, v_n)
        return numpy.abs(v_nInterp - v_a), (x_a[-1] - x_a[0])/(len(x_a) - 1)

    def compute_norm(self, errors, grid_spacing):
        """Computes a norm.  Uniform, linear elements are assumed, so the
        integral over the domain can be approximated by the trapezium rule
        (for the L1 norm at least).
        """
        E = numpy.array(errors)
        if self.norm == 1: 
            I = (E[0]/2 + sum(E[1:-1]) + E[-1]/2)*grid_spacing
        elif self.norm == 2:
            # does not make much difference, but previously I = (E[0]**2/2 +
            # sum(E[1:-1]**2) + E[-1]**2/2)*grid_spacing
            I = (E[0]**2 + sum(E[:-1]*E[1:]) + 2*sum(E[1:-1]**2) + \
                     E[-1]**2)*grid_spacing/3
        else: 
            print("Only l1 or l2 norms may be computed; i.e. self.norm "
                  "must be 1 or 2.")
            sys.exit()
        return I**(1/float(self.norm))
   

class WriteConvergenceRate(Command):
    """Analyses the rate at which an error norm is converging to zero.  For
    this to work, the previously analysed grid needs to be coarser by a
    factor of two.  The output report has two columns: the first has
    identifiers in the format
    modelname_dimension_fieldID_normtype_gridpair; the second has the
    values of the convergence rates.

    """

    def __init__(self, rates_file, dimension):
        self.rates_file = rates_file
        self.dim = str(dimension)

    def execute(self, level_name, value):
        if level_name == 'stem':
            self.stem = value
        if level_name == 'model':
            self.model = value
        if level_name == 'field':
            # split up field name into useful bits
            self.phase_index = re.sub('Phase(.)::(.*)', '\\1', value)
            self.var_name = re.sub('Phase(.)::(.*)', '\\2', value)
            self.field_short = str.lower(self.var_name)+self.phase_index
        if level_name == 'norm':
            self.norm = value
        if level_name == 'grid':
            grid = value
            v = self.get_norm(grid)
            
            # can only start computing rates from the 2nd grid
            if grid!='A':
                # print a new ID and the rate
                key = self.model+'_'+self.dim+'d_'+self.field_short+\
                      '_l'+str(self.norm)+'_'+self.grid0+grid
                self.rates_file.write('{0}{1:12.6f}\n'.format(
                    key,numpy.log2(self.v0/v)))
            self.v0 = v
            self.grid0 = grid

    def get_norm(self): 
        # defer to subclasses
        pass


class Dim1_WriteConvergenceRate(WriteConvergenceRate):
    """Can be run after WriteErrorNorm1D has been handled.  N.B. Has a
    different ordering and formatting of the IDs on the LHS.
    """
    def __init__(self, rates_file, norms_file):
        WriteConvergenceRate.__init__(self, rates_file, 1)
        self.norms_file = norms_file

    def get_norm(self, grid): 
        suffix = self.model+'_1d_'+grid
        key = '{0}_{1}_l{2:g}'.format(suffix, self.field_short,
                                      self.norm)
        return find(self.norms_file, key)
        

class DimN_WriteConvergenceRate(WriteConvergenceRate):
    """Produces a similar output to sibling Dim1, but uses
    information associated with multidimensional VTU files incorporating
    interp_using_analytic.  This is a more versatile approach as
    explained at the top of the file.
    """
    def __init__(self, rates_file, dimension):
        WriteConvergenceRate.__init__(self, rates_file, dimension)

    def get_norm(self, grid): 
        filename = self.stem+'_'+self.model+'_'+self.dim+'d_'+\
                   grid+'.stat'
        if self.norm==1:
            self.calc_type = "integral"
        elif self.norm==2:
            self.calc_type = "l2norm"
        return fluidity_tools.stat_parser(filename)\
            ['Phase'+self.phase_index]\
            ['Analytic'+self.var_name+'Error']\
            [self.calc_type][-1]
                
            

## FOR PUBLIC USE 

def interp_using_analytic(x, analytic_filename):
    """Interpolates at point x using the profile given by
    analytic_filename."""
    [x_a, v_a] = read_analytic_solution(analytic_filename)
    return numpy.interp(x, x_a, v_a)

       
def find_rate(key):
    """Opens and searches the rates report and returns the value associated
    with key.
    """
    with open(error_rates_filename, "r") as f:
        v = find(f, key)
        return v

       
class BuckleyLeverettTestSuite:
    """Class to help run and postprocess tests.  numerical_filename_stem is
    a string to be appended with the elements in model_name_list and
    grid_name_list, where the latter is one of three elements in the
    list grid_name_list_per_dimension.  analytic_filename_stem is a
    string to be appended with shorthand names for the fields in
    question, which in turn are translated from the elements in
    field_name_list.  If a field_name is 'Phase2::Saturation', then the
    shorthand is saturation2.

    """
    def __init__(self, numerical_filename_stem, model_name_list,
                 grid_name_list_per_dimension, analytic_filename_stem,
                 field_name_list):
        
        self.numerical_filename_stem = numerical_filename_stem
        self.analytic_filename_stem = analytic_filename_stem

        # convert lists
        self.grid_handler_list = []
        for grid_name_list in grid_name_list_per_dimension:
            self.grid_handler_list.append(
                HandlerList('grid', grid_name_list))
        self.dim_handler_list = HandlerList('dim', ('1d', '2d', '3d'))
        self.model_handler_list = HandlerList('model', model_name_list)
        self.field_handler_list = HandlerList('field', field_name_list)
        self.norm_handler_list = HandlerList('norm', (1, 2))
        
            
    def new_handler(self, handlers):
        """Short wrapper for CompositeHandler."""
        return CompositeHandler('stem', self.numerical_filename_stem,
                                handlers )
 
 
    def process_folder(self):
        if verbose: print '\nProcessing folder '
        # build handler structure
        handler = self.new_handler(
            self.model_handler_list.expand(
                self.dim_handler_list.splice(
                    # n.b. splice lists here
                    self.grid_handler_list)))
        # pass it a command 
        handler.handle( RunSimulation() )

        
    def generate_error_report_1D(self):
        if verbose: print '\nComputing 1D error norms'
        # open write file
        with open(error_norms_filename, "w") as f_norms:
            # build handler structure
            handler = self.new_handler(
                self.model_handler_list.expand(
                    # 1D, so just expand grid_handler_list[0]
                    self.grid_handler_list[0].expand(
                        self.field_handler_list.expand(
                            self.norm_handler_list))))
            # pass it a command 
            handler.handle(
                WriteErrorNorm1D(self.analytic_filename_stem, f_norms))
    
    
    def generate_convergence_report(self):
        if verbose: print '\nComputing convergence rates'
        # open write file and start looping over dimensions
        with open(error_rates_filename, "w") as f_rates:
            for i, dim in enumerate((1, 2, 3)):
                # build handler structure
                handler = self.new_handler(
                    self.model_handler_list.expand(
                        self.field_handler_list.expand(
                            self.norm_handler_list.expand(
                                self.grid_handler_list[i]))))
                # pass it a command
                if dim==1:
                    with open(error_norms_filename, "r") as f_norms:
                        handler.handle(
                            Dim1_WriteConvergenceRate(f_rates, f_norms))
                else:
                    handler.handle(
                        DimN_WriteConvergenceRate(f_rates, dim))
        
    
    def generate_reports(self):
        self.generate_error_report_1D()
        self.generate_convergence_report()
