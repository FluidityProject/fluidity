#!/usr/bin/env python

"""Expression generation via sympy for MMS-based tests associated with
the Darcy IMPES code.

More documentation to come."""

from numpy import array, ndarray, linspace, meshgrid, shape, vstack, \
    min, max, zeros
from sympy import Symbol, Function, diff, integrate, sin, cos, pi, exp, sqrt
from sympy.utilities.lambdify import lambdify
import re
import fileinput

## HELPERS

spacetime = (Symbol('x'), Symbol('y'), Symbol('z'), Symbol('t'))
saturations = (Symbol('s1'), Symbol('s2'))

class Vector:
    def __init__(self, components_or_dim):
        if isinstance(components_or_dim, int):
            self.components = array([Symbol('0')]*components_or_dim)
        else:
            self.components = array(components_or_dim)
            # protect against failure of array to wrap symbols properly
            try:
                len(self.components)
            except:
                self.components = array((components_or_dim,))
    def clone(self):
        return Vector(self.components)
    def __setitem__(self, index, value):
        self.components[index] = value
    def __getitem__(self, index):
        return self.components[index]
    def __neg__(self):
        return Vector(-self.components)
    # addition and subtraction can be done elementwise; allow other
    # collection types to be encountered
    def __add__(self, other):
        try:
            return Vector(self.components + other.components)
        except:
            return Vector(self.components + array(other))
    def __radd__(self, other):
        try:
            return Vector(other.components + self.components)
        except:
            return Vector(array(other) + self.components)
    def __sub__(self, other):
        try:
            return Vector(self.components - other.components)
        except:
            return Vector(self.components - array(other))
    def __rsub__(self, other):
        try:
            return Vector(other.components - self.components)
        except:
            return Vector(array(other) - self.components)
    # expect multiplication and division with scalars only.  Note that
    # explicitly looping over vector components seems to be more stable
    # than letting the sub-objects handle it
    def __mul__(self, other):
        result = self.clone()
        for i in range(result.dim()):
            result[i] *= other
        return result
    def __div__(self, other):
        result = self.clone()
        for i in range(result.dim()):
            result[i] /= other
        return result
    def __rmul__(self, other):
        return self * other
    def __str__(self):
        return 'Vector({0})'.format(self.components)
    def __abs__(self):
        result = 0.
        for u in self.components:
            result += u**2
        return sqrt(result)
    def dim(self):
        return len(self.components)
    def divergence(self):
        result = 0.
        for i, u in enumerate(self.components):
            result += diff(u, spacetime[i])
        return result

def grad(scalar, dim):
    result = Vector(dim)
    for i in range(dim):
        result[i] = diff(scalar, spacetime[i])
    return result

def div(vector):
    try:
        return vector.divergence()
    except TypeError:
        return diff(vector, spacetime[0])

def spacetime_lambdify(expression, dim):
    Xt = spacetime[:dim] + spacetime[3:]
    return lambdify(Xt, expression)

def write_expr(genfile, is_text, key, value, args=None, collect_args=False):
    if is_text:
        # capitalise if in text mode
        key = str.upper(key)
    genfile.write('    "'+key+'" : ')

    if is_text:
        # text format
        genfile.write('"')
        try:
            # write a vector 
            for component in value:
                genfile.write(str(component)+' ')
        except:
            # write a scalar 
            genfile.write(str(value))
        genfile.write('"')

    else:
        # python format
        if args is None:
            # write a number (or numbers)
            if len(shape(value)) > 0:
                # protect against bad formatting of numpy.ndarray
                genfile.write('(')
                for v in value:
                    genfile.write('{0}, '.format(v))
                genfile.write(')')
            else:
                genfile.write('{0}'.format(value))
        else:
            # write a lambda
            if not collect_args:
                genfile.write('lambda')
                for i, a in enumerate(args):
                    if i > 0:
                        genfile.write(','.format(a))
                    genfile.write(' {0}'.format(a))
                genfile.write(': ')
            else:
                genfile.write('lambda {0}: '.format(args))
            genfile.write(str(value))
            
    # end
    genfile.write(',\n\n')
    
        
## FOR PUBLIC USE

def generate_coords(extents, point_num):
    """ Generates a Cartesian array of points for checking
    samples of e.g. CFL number"""
    
    # extents can be [[xmin,...,tmin],[xmax,...,tmax]] or just
    # [xmax,...,tmax]
    if len(shape(extents))==1:
        n = len(extents)
        extents = [tuple([0.]*n), extents]
    else:
        n = len(extents[0])
    # point_num can be [nx,...,nt] or just n
    try:
        len(point_num)
    except TypeError:
        point_num = tuple([point_num]*n)
    points = []
    for i in range(n):
        points.append(linspace(extents[0][i], extents[1][i], point_num[i]))
    X = meshgrid(*points)
    for i in range(n):
        x = X[i].flatten()
        if i==0:
            coords = x
        else:
            coords = vstack((coords, x))
    return coords.transpose()


class GenericScalar:
    def __init__(self, symbol, phase_num, expression, diffusivity=Symbol('0'),
                 external_source=Symbol('0'), scale=1.):
        self.symbol = symbol
        self.expression = expression
        self.phase_num = phase_num
        self.diffusivity = diffusivity
        self.external_source = external_source
        self.scale = scale
        
    def compute(self, phase_num, phi, sat, darcy_vel):
        if self.phase_num != phase_num:
            return
        t = Symbol('t')
        dim = darcy_vel.dim()
        self.manuf_source = diff(phi*sat*self.expression, t) \
            + div(darcy_vel*self.expression) \
            - div(phi*sat*self.diffusivity*grad(self.expression, dim)) \
            - self.external_source.subs(self.symbol, self.expression)
        self.external_source_grad = diff(self.external_source, self.symbol)
    
    def write_expressions(self, genfile, is_text, phase_num, dim):
        if self.phase_num != phase_num:
            return
        # append domain dimension to each variable name
        suf = '_' + str(dim) + 'D'
        if dim == 1:
            write_expr(genfile, is_text,
                       str(self.symbol)+str(self.phase_num)+'_scale',
                       self.scale, spacetime[:dim] + spacetime[3:])
        write_expr(genfile, is_text, str(self.symbol)+str(self.phase_num)+suf,
                   self.expression, spacetime[:dim] + spacetime[3:])
        if self.diffusivity != Symbol('0'):
            write_expr(genfile, is_text,
                       'diffusivity_'+str(self.symbol)+str(self.phase_num)+suf,
                       self.diffusivity, spacetime[:dim] + spacetime[3:])
        write_expr(genfile, is_text,
                       'external_source_'+str(self.symbol)+str(self.phase_num)+suf,
                       self.external_source, spacetime[:dim] + spacetime[3:])
        write_expr(genfile, is_text,
                       'external_source_gradient_'+str(self.symbol)+str(self.phase_num)+suf,
                       self.external_source_grad, spacetime[:dim] + spacetime[3:])
        write_expr(genfile, is_text,
                   'manufactured_source_'+str(self.symbol)+str(self.phase_num)+suf,
                   self.manuf_source, spacetime[:dim] + spacetime[3:])

        

class ManufacturedSolution:
    def __init__(self, dim, g_dir, p, s2, generic_scalars=[]):
        self.dim = dim
        self.g_dir = array(g_dir)
        self.p = []
        self.grad_p = []
        for i in range(2):
            try:
                # p is given as multiple fields
                self.p.append(p[i])
            except TypeError:
                # p is given as one common field
                self.p.append(p)
            self.grad_p.append(grad(self.p[i], dim))
        self.s = (1 - s2, s2)
        self.u = [Vector(dim) for i in range(2)]
        self.q = [Symbol('0') for i in range(2)]
        self.generic_scalars = generic_scalars
        
    def compute_phase(self, phase_num, phi, K, mu, rho, g_mag):
        t = Symbol('t')
        g = g_mag*self.g_dir
        i = phase_num - 1
        grad_p_eff = self.grad_p[i] - rho*g
        try:
            grad_p_eff = self.grad_p[i] - rho*g
        except TypeError:
            grad_p_eff = self.grad_p[i] - array(rho*g)
        self.u[i] = -K(self.s)/mu*grad_p_eff
        self.q[i] = diff(phi*self.s[i], t) + div(self.u[i])
        for c in self.generic_scalars:
            c.compute(phase_num, phi, self.s[i], self.u[i])


    def check_max_velocity(self, phase_num, coords_and_time):
        i = phase_num-1
        U = spacetime_lambdify(self.u_mag[i], self.dim)
        u_max = 0.
        if len(shape(coords_and_time))==1:
            coords_and_time = (coords_and_time,)
        for Xt in coords_and_time:
            u_max = max((u_max, U(*Xt)))
        return u_max

    def compute_divergence_minus_source(self):
        result = Symbol('0')
        for i in [0, 1]:
            result += div(self.u[i]) - self.q[i]
        return result

    def write_main_expressions(self, genfile, is_text):
        # append domain dimension to each variable name
        suf = '_' + str(self.dim) + 'D'
        write_expr(genfile, is_text, 'gravity_direction'+suf,
                         self.g_dir)
        for i in range(2):
            write_expr(genfile, is_text, 'pressure'+str(i+1)+suf,
                       self.p[i], spacetime[:self.dim] + spacetime[3:])
            write_expr(genfile, is_text, 'saturation'+str(i+1)+suf, 
                       self.s[i], spacetime[:self.dim] + spacetime[3:])
        
    def write_phase_expressions(self, genfile, is_text, phase_num):
        # append domain dimension to each variable name
        suf = '_' + str(self.dim) + 'D'
        i = phase_num-1
        for j in range(self.dim):
            xi = spacetime[j]
            try:
                uij = self.u[i][j]
            except TypeError:
                uij = self.u[i]
            write_expr(
                genfile, is_text, 'darcy_velocity'+str(i+1)+'_'+str(xi)+suf,
                uij, spacetime[:self.dim] + spacetime[3:])
        write_expr(genfile, is_text,
                   'darcy_velocity'+str(i+1)+'_magnitude'+suf,
                   abs(self.u[i]), spacetime[:self.dim] + spacetime[3:])
        write_expr(genfile, is_text,
                   'source_saturation'+str(i+1)+suf,
                   self.q[i], spacetime[:self.dim] + spacetime[3:])
        for c in self.generic_scalars:
            c.write_expressions(genfile, is_text, phase_num, self.dim)

        

class SolutionHarness:
    def __init__(self, D, T, g_mag, ka, phi, mu, rho, K, p_scale, s_scale,
                 solution_name='solution_expressions'):
        # parameters
        self.genfilename = solution_name+'.py'
        self.D = D                  # domain extents
        self.T = T                  # finish time 
        self.g_mag = g_mag          # gravity magnitude
        self.ka = ka                # absolute permeability
        self.phi = phi              # porosity
        # the following should be two elements long
        self.mu = mu            # viscosity
        self.rho = rho;         # density
        self.K = K              # relperm (functions of s)
        # pressure field magnitudes, for scaling error norms appropriately
        self.p_scale = [Symbol('0') for i in range(2)]
        for i in range(2):
            try:
                # p_scale is given as multiple fields
                self.p_scale[i] = p_scale[i]
            except TypeError:
                # p_scale is given as one common field
                self.p_scale[i] = p_scale
        # repeat for saturation
        self.s_scale = [Symbol('0') for i in range(2)]
        for i in range(2):
            try:
                # s_scale is given as multiple fields
                self.s_scale[i] = s_scale[i]
            except TypeError:
                # s_scale is given as one common field
                self.s_scale[i] = s_scale

    # def compute_reference_p_grad( 
    #         self, phase_num, reference_saturation,
    #         reference_darcy_velocity ):
    #     i = phase_num-1
    #     return self.rho[i]*self.g_mag - \
    #         reference_darcy_velocity*self.mu[i]/ \
    #         self.K(reference_saturation)
        
    def write_dict(self, manufactured_solution_list):
        """Writes two dictionaries - one returning pure Python code as
        numbers and lambdas, the other returning text-formatted Python
        code for substitution into XML options files.
        """
        
        f = open(self.genfilename, 'w')
        f.write('# THIS FILE HAS BEEN AUTOMATICALLY GENERATED.\n\n')

        # although sympy has been used to generate the following
        # dictionary, clients of the dictionary need only depend on
        # numpy.
        f.write('from numpy import sin, cos, pi, exp, sqrt \n\n')

        dict_names = ['py_dict', 'text_dict']
        for dn in dict_names:
            is_text = (dn=='text_dict')
            # write LHS of dictionary assignment 
            f.write(dn + ' = {\n')

            write_expr(f, is_text, 'domain_extents', self.D)
            write_expr(f, is_text, 'finish_time', self.T)
            for i, phase_num in enumerate((1, 2)):
                write_expr(f, is_text, 'pressure{0}_scale'.format(phase_num),
                           self.p_scale[i])
                write_expr(f, is_text, 'saturation{0}_scale'.format(phase_num),
                           self.s_scale[i])
            write_expr(f, is_text, 'gravity_magnitude', self.g_mag)
            write_expr(f, is_text, 'viscosity1', self.mu[0])
            write_expr(f, is_text, 'viscosity2', self.mu[1])
            write_expr(f, is_text, 'density1', self.rho[0])
            write_expr(f, is_text, 'density2', self.rho[1])
            write_expr(f, is_text, 'permeability1',
                       self.K[0](saturations), saturations, collect_args=True)
            write_expr(f, is_text, 'permeability2', 
                       self.K[1](saturations), saturations, collect_args=True)
            write_expr(f, is_text, 'porosity', self.phi)
            write_expr(f, is_text, 'absolute_permeability', self.ka)
            # solutions (one per domain dimension)
            for ms in manufactured_solution_list:
                ms.write_main_expressions(f, is_text)
                for i, phase_num in enumerate((1, 2)):
                    ms.compute_phase(phase_num, self.phi,
                                     self.K[i], self.mu[i],
                                     self.rho[i], self.g_mag)
                    ms.write_phase_expressions(f, is_text, phase_num)

            # finish writing assingnment 
            f.write('    }\n\n')

        # finish writing file
        f.close()

        # now need to reformat expressions with r.e.-based search and
        # replace
        flag_text = False
        for line in fileinput.input(self.genfilename, inplace=1):
            # sympy writes fractions with integers, e.g. 1/2, so need to
            # append decimal points to all integers that aren't acting
            # as exponents, array subscripts or labels.  "1/2." is fine
            # though
            line = re.sub('([^e*][ /*\-+()][0-9]+)([ /*\-+()"])', '\\1.\\2',
                          line.rstrip())
            # a second pass is needed due to that pattern absorbing
            # neighbouring candidates
            line = re.sub('([^e*][ /*\-+()][0-9]+)([ /*\-+()"])', '\\1.\\2',
                          line.rstrip())
            # another important conversion
            line = re.sub('Abs', 'abs', line.rstrip())
            
            # if we are in text_dict, replace x with X[0], etc.
            if flag_text:
                # disregard words
                line = re.sub('([^a-z])x([^a-z])', '\\1X[0]\\2', line.rstrip())
                line = re.sub('([^a-z])y([^a-z])', '\\1X[1]\\2', line.rstrip())
                line = re.sub('([^a-z])z([^a-z])', '\\1X[2]\\2', line.rstrip())
            else:
                if line.find('text_dict') != -1:
                    flag_text = True

            # write the processed line
            print(line)
            

    def reality_check(self, manufactured_solution_list):
        results = []
        for ms in manufactured_solution_list:
            results.append( ms.compute_divergence_minus_source() )
        return results
    
