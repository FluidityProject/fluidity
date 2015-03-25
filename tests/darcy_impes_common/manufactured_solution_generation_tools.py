from numpy import array, ndarray, linspace, meshgrid, shape, vstack, min, max, zeros
from sympy import Symbol, Function, diff, integrate, sin, cos, pi, exp, sqrt
from sympy.utilities.lambdify import lambdify
import re
import fileinput

## HELPERS

spacetime = (Symbol('x'), Symbol('y'), Symbol('z'), Symbol('t'))
saturations = (Symbol('s1'), Symbol('s2'))

def grad(scalar, dim):
    vector = array([Symbol('0')]*dim)
    for i in range(dim):
        vector[i] = diff(scalar, spacetime[i])
    return vector

def div(vector):
    try:
        scalar = 0.
        for i in range(len(vector)):
            scalar = scalar + diff(vector[i], spacetime[i])
        return scalar
    except TypeError:
        return diff(vector, spacetime[0])

def mag(vector):
    try:
        scalar = 0.
        for i in range(len(vector)):
            scalar = scalar + vector[i]**2
        return sqrt(scalar)
    except TypeError:
        return abs(vector)

def spacetime_lambdify(expression, dim):
    Xt = spacetime[:dim] + spacetime[3:]
    return lambdify(Xt, expression)

def write_expr(genfile, is_text, key, value, args=None, separate_args=False):
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
            if separate_args:
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
    
    
class ManufacturedSolution:
    def __init__(self, dim, g_dir, p, s2):
        self.dim = dim
        self.g_dir = array(g_dir)
        self.p = [Symbol('0') for i in range(2)]
        self.grad_p = [[Symbol('0') for j in range(dim)] for i in range(2)]
        for i in range(2):
            try:
                # p is given as multiple fields
                self.p[i] = p[i]
            except TypeError:
                # p is given as one common field
                self.p[i] = p
            self.grad_p[i] = grad(self.p[i], dim)
        self.s = (1 - s2, s2)
        self.u = [[Symbol('0') for j in range(dim)] for i in range(2)]
        self.q = [Symbol('0') for i in range(2)]
        
    def compute_phase(self, phase_num, phi, K, mu, rho, g_mag):
        t = Symbol('t')
        g = g_mag*self.g_dir
        i = phase_num - 1
        try:
            grad_p_eff = self.grad_p[i] - rho*g
        except TypeError:
            grad_p_eff = self.grad_p[i] - array(rho*g)
        # explicitly looping over vector components here seems to be
        # more stable than letting the objects handle it
        for j, gpe in enumerate(grad_p_eff):
            # n.b. for flexibility, K is a function of both sats
            self.u[i][j] = -K(self.s)/mu*gpe
        self.q[i] = diff(phi*self.s[i], t) + div(self.u[i])

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
                       self.p[i], spacetime[:self.dim] + spacetime[3:],
                       separate_args=True)
            write_expr(genfile, is_text, 'saturation'+str(i+1)+suf, 
                       self.s[i], spacetime[:self.dim] + spacetime[3:],
                       separate_args=True)
        
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
            write_expr(genfile, is_text, 'darcy_velocity'+str(i+1)+'_'+str(xi)+suf,
                   uij, spacetime[:self.dim] + spacetime[3:], separate_args=True)
        write_expr(genfile, is_text,
                   'darcy_velocity'+str(i+1)+'_magnitude'+suf,
                   mag(self.u[i]), spacetime[:self.dim] + spacetime[3:],
                   separate_args=True)
        write_expr(genfile, is_text,
                   'source_saturation'+str(i+1)+suf,
                   self.q[i], spacetime[:self.dim] + spacetime[3:],
                   separate_args=True)
        

class SolutionHarness:
    def __init__(self, D, T, g_mag, ka, phi, mu, rho, K, p_scale, s2_scale,
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
        # field magnitudes, for scaling error norms appropriately
        self.p_scale = p_scale 
        self.s2_scale = s2_scale 

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
            write_expr(f, is_text, 'pressure_scale', self.p_scale)
            write_expr(f, is_text, 'saturation_scale', self.s2_scale)
            write_expr(f, is_text, 'gravity_magnitude', self.g_mag)
            write_expr(f, is_text, 'viscosity1', self.mu[0])
            write_expr(f, is_text, 'viscosity2', self.mu[1])
            write_expr(f, is_text, 'density1', self.rho[0])
            write_expr(f, is_text, 'density2', self.rho[1])
            write_expr(f, is_text, 'permeability1',
                       self.K[0](saturations), saturations)
            write_expr(f, is_text, 'permeability2', 
                       self.K[1](saturations), saturations)
            write_expr(f, is_text, 'porosity', self.phi)
            write_expr(f, is_text, 'absolute_permeability', self.ka)
            # solutions (one per domain dimension)
            for ms in manufactured_solution_list:
                ms.write_main_expressions(f, is_text)
                for i, phase_num in enumerate((1, 2)):
                    # TEMP
                    # write_expr(f, is_text, 
                    #               'grad_p_eff'+str(phase_num),
                    #               ms.grad_p[i] -
                    #               array(self.rho[i]*self.g_mag*ms.g_dir),
                    #               spacetime[:ms.dim] + spacetime[3:],
                    #               separate_args=True)
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
    
