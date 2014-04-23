from numpy import array, linspace, meshgrid, shape, vstack, min, max, zeros
from sympy import Symbol, Function, diff, integrate, sin, cos, pi, exp, sqrt
from sympy.utilities.lambdify import lambdify
import re
import fileinput


## HELPERS

spacetime = (Symbol('x'), Symbol('y'), Symbol('z'), Symbol('t'))

def grad(scalar, dim):
    vector = array([Symbol('0')]*dim)
    for i in range(dim):
        vector[i] = diff(scalar, spacetime[i])
    return vector

def div(vector):
    scalar = 0.
    for i in range(len(vector)):
        scalar = scalar + diff(vector[i], spacetime[i])
    return scalar

def mag(vector):
    scalar = 0.
    for i in range(len(vector)):
        scalar = scalar + vector[i]**2
    return sqrt(scalar)

def spacetime_lambdify(expression, dim):
    Xt = spacetime[:dim] + spacetime[3:]
    return lambdify(Xt, expression)

def write_expression_py(genfile, key, value):
    genfile.write('    "'+key+'" : {0}'.format(value))
    genfile.write(',\n\n')

def write_expression_txt(genfile, key, value):
    genfile.write('    "'+key+'" : "')
    try:
        for component in value:
            genfile.write(str(component)+' ')
    except:
        genfile.write(str(value))
    genfile.write('",\n\n')
    
        
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
        self.p = p
        self.s = (1 - s2, s2)
        self.grad_p = grad(p, dim)
        self.u_mag = [Symbol('0')]*2
        self.q = [Symbol('0')]*2
        
    def compute_phase(self, phase_num, phi, K, mu, rho, g_mag):
        t = Symbol('t')
        g = g_mag*self.g_dir
        i = phase_num - 1
        # n.b. K is a function of s[i]
        u = -K(self.s[i])/mu*(self.grad_p - rho*g)
        self.u_mag[i] = mag(u)
        self.q[i] = -diff(phi*self.s[i], t) + div(u)

    def check_max_velocity(self, phase_num, coords_and_time):
        i = phase_num-1
        U = spacetime_lambdify(self.u_mag[i], self.dim)
        u_max = 0.
        if len(shape(coords_and_time))==1:
            coords_and_time = (coords_and_time,)
        for Xt in coords_and_time:
            u_max = max((u_max, U(*Xt)))
        return u_max

    def write_main_expressions(self, genfile):
        # append domain dimension to each variable name
        suf = '_' + str(self.dim) + 'D'
        write_expression_txt(genfile, 'GRAVITY_DIRECTION'+suf, self.g_dir)
        write_expression_txt(genfile, 'PRESSURE'+suf, self.p)
        write_expression_txt(genfile, 'SATURATION2'+suf, self.s[1])
        
    def write_phase_expressions(self, genfile, phase_num):
        # append domain dimension to each variable name
        suf = '_' + str(self.dim) + 'D'
        i = phase_num-1
        write_expression_txt(
            genfile, 'DARCY_VELOCITY'+str(i+1)+'_MAGNITUDE'+suf,
            self.u_mag[i])
        write_expression_txt(
            genfile, 'SOURCE_SATURATION'+str(i+1)+suf, self.q[i])
        

class SolutionHarness:
    def __init__(self, D, g_mag, ka, phi, mu, rho, K, p_scale, s2_scale,
                 solution_name='solution_expressions'):
        # parameters
        self.genfilename = solution_name+'.py'
        self.D = D                  # domain extents
        self.g_mag = g_mag          # gravity magnitude
        self.ka = ka                # absolute permeability
        self.phi = phi              # porosity
        # the following should be two elements long
        self.mu = mu            # viscosity
        self.rho = rho;         # density
        # the following should be a lambda func
        self.K = K              # relperm
        # field magnitudes, for scaling error norms appropriately
        self.p_scale = p_scale 
        self.s2_scale = s2_scale 

    def compute_reference_p_grad( 
            self, phase_num, reference_saturation,
            reference_darcy_velocity ):
        i = phase_num-1
        return self.rho[i]*self.g_mag - \
            reference_darcy_velocity*self.mu[i]/ \
            self.K(reference_saturation)
        
    def write_dict(self, manufactured_solution_list):

        f = open(self.genfilename, 'w')
        f.write('# THIS FILE HAS BEEN AUTOMATICALLY GENERATED.\n\n')
                
        # start writing dictionary
        f.write('solution_dict = {\n\n')
        
        # stuff that will be used again in Python code
        write_expression_py(f, 'domain_extents', self.D)
        
        # constants
        write_expression_py(f, 'pressure1_scale', self.p_scale)
        write_expression_py(f, 'saturation2_scale', self.s2_scale)
        write_expression_txt(f, 'GRAVITY_MAGNITUDE', self.g_mag)
        write_expression_txt(f, 'VISCOSITY1', self.mu[0])
        write_expression_txt(f, 'VISCOSITY2', self.mu[1])
        write_expression_txt(f, 'DENSITY1', self.rho[0])
        write_expression_txt(f, 'DENSITY2', self.rho[1])
        write_expression_txt(f, 'POROSITY', self.phi)
        write_expression_txt(f, 'ABSOLUTE_PERMEABILITY', self.ka)

        # solutions (one per domain dimension)
        for ms in manufactured_solution_list:
            ms.write_main_expressions(f)
            for i, phase_num in enumerate((1, 2)):
                ms.compute_phase(phase_num, self.phi,
                                 self.K, self.mu[i],
                                 self.rho[i], self.g_mag)
                ms.write_phase_expressions(f, phase_num)
            
        # finish writing dictionary
        f.write('    }')
        f.close()

        # Now need to reformat expressions with r.e.-based search and
        # replace
        for line in fileinput.input(self.genfilename, inplace=1):
            # replace x with X[0] (careful with the exp function), etc.
            line = re.sub('([^e])x([^p])', '\\1X[0]\\2', line.rstrip())
            line = re.sub('y', 'X[1]', line.rstrip())
            line = re.sub('z', 'X[2]', line.rstrip())
            # sympy writes fractions with integers, e.g. 1/2, so
            # need to append decimal points
            line = re.sub('([^e*][ /\*-+()][0-9]+)([ /\*-+()])', '\\1.\\2', line.rstrip())
            # a second pass is needed here
            line = re.sub('([^e*][ /\*-+()][0-9]+)([ /\*-+()])', '\\1.\\2', line.rstrip())
            print(line)
