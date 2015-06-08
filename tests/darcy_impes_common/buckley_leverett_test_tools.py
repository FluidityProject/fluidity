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
from test_tools import Command, CommandList, HandlerLevel, CompositeHandler, WriteToReport, RunSimulation, error_norms_filename, error_rates_filename, find, find_rate, verbose

error_norms_1D_filename = "error_norms_1d.txt"

verbose(True)
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

        
class Preprocess(Command):
    def __init__(self, clean=False):
        # alternative role
        self.clean = clean

    def execute(self, level_name, value, indent):
        # most levels just record level names
        if level_name == 'stem':
            self.stem = value
        elif level_name == 'model':
            self.model = value
        elif level_name == 'dim':
            self.dim = value
        elif level_name == 'mesh_suffix':
            # last level
            mesh_suffix = value
            
            case = self.stem+'_'+self.model+'_'+str(self.dim)+'d'
            filename_new = case+'_'+mesh_suffix+'.diml'
            if mesh_suffix!='A':
                if not self.clean:
                    # for the most coarse mesh, an options file should
                    # exist.  For the other meshes, adapt it
                    filename_A = case+'_A.diml'
                    with open(filename_A, 'r') as f_A:
                        f = open(filename_new, 'w')
                        for line in f_A:
                            f.write(line.replace('_A', '_'+mesh_suffix))
                        f.close()
                else:
                    # undo the above
                    os.remove(filename_new)

                
class WriteErrorNorm1D(Command):
    """Compares a 1D numerical solution with the appropriate analytical
    solution, writing the error norm to a report.  The output report
    has two columns: the first has identifiers in the format
    model_name_1d_mesh_suffix_field_ID_normtype; the second has the values
    of the error norms.
    """

    def __init__(self, analytic_stem, norms_file):
        self.analytic_stem = analytic_stem
        self.norms_file = norms_file

    def execute(self, level_name, value, indent):
        if level_name == 'stem':
            self.stem = value
        elif level_name == 'model':
            self.model = value
        elif level_name == 'mesh_suffix':
            self.mesh_suffix = value
        elif level_name == 'field':
            # translate the field_name to shorthand form
            field = value
            self.field_short = str.lower(re.sub('Phase(.)::(.*)', '\\2\\1',
                                                field))
            # store the errors and dx for this field
            num_suf = self.model+'_1d_'+self.mesh_suffix
            [self.errs, self.dx] = self.compute_abs_errors(
                num_suf, self.field_short, field)
            
        elif level_name == 'norm':
            self.norm = value
            
            # print the ID and value.  N.B. key is different to 
            # num_suf passed to compute_abs_errors.
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
   

class BLWriteToReport(WriteToReport):
    """Overrides WriteToReport and its get_norm method
    in order to allow customised computation of 1D errors.
    """

    def __init__(self, norms_file, rates_file):
        WriteToReport.__init__(self, norms_file, rates_file)

    def get_norm(self):
        if self.dim==1:
            key = '{0}_1d_{1}_l{2:g}_{3}'.format(
                self.model, self.field_short, self.norm,
                self.mesh_suffix)
            v = find(error_norms_1D_filename, key)
        else:
            filename = self.stem+'_'+self.model+'_'+\
                       str(self.dim)+'d_'+self.mesh_suffix+'.stat'
            if self.norm==1:
                self.calc_type = "integral"
            elif self.norm==2:
                self.calc_type = "l2norm"
            v = fluidity_tools.stat_parser(filename)\
                ['Phase'+self.phase_index]\
                ['Analytic'+self.var_name+'Error']\
                [self.calc_type][-1]
        return v
                
            

## FOR PUBLIC USE 
    

def interp_using_analytic(x, analytic_filename):
    """Interpolates at point x using the profile given by
    analytic_filename."""
    [x_a, v_a] = read_analytic_solution(analytic_filename)
    return numpy.interp(x, x_a, v_a)

       
class BuckleyLeverettTestSuite:
    """
    Class to help run and postprocess tests.  numerical_filename_stem is
    a string to be appended with the elements in model_name_list and
    mesh_suffix_list, where the latter is one of three elements in the
    list mesh_suffix_list_per_dimension.  analytic_filename_stem is a
    string to be appended with shorthand names for the fields in
    question, which in turn are translated from the elements in
    field_name_list.  If a field_name is 'Phase2::Saturation', then the
    shorthand is saturation2.
    """
    def __init__(self, numerical_filename_stem, model_name_list,
                 mesh_suffix_list_per_dimension, analytic_filename_stem,
                 field_name_list):
        
        self.numerical_filename_stem = numerical_filename_stem
        self.analytic_filename_stem = analytic_filename_stem

        # convert lists
        self.mesh_suffix_handler_levels = []
        for mesh_suffix_list in mesh_suffix_list_per_dimension:
            self.mesh_suffix_handler_levels.append(
                HandlerLevel('mesh_suffix', mesh_suffix_list))
        self.dim_handler_level = HandlerLevel('dim', (1, 2, 3))
        self.model_handler_level = HandlerLevel('model', model_name_list)
        self.field_handler_level = HandlerLevel('field', field_name_list)
        self.norm_handler_level = HandlerLevel('norm', (1, 2))
        
            
    def new_handler(self, handler_level):
        """Short wrapper for CompositeHandler."""
        return CompositeHandler('stem', self.numerical_filename_stem,
                                handler_level )
 
 
    def process_folder(self):
        if verbose(): print '\nProcessing folder '
        # build handler structure
        handler = self.new_handler(
            self.model_handler_level.add_sub(
                self.dim_handler_level.splice(
                    # n.b. splice lists here
                    self.mesh_suffix_handler_levels)))
        # pass it a command 
        handler.handle(CommandList((Preprocess(),
                                    RunSimulation(),
                                    Preprocess(clean=True))))

        
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
