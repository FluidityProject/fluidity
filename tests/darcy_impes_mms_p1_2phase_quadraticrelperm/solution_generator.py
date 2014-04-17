from sympy import Symbol, Function, diff, integrate, sin, cos, pi, exp
from manufactured_solution_generation_tools import \
    ManufacturedSolution, SolutionHarness, generate_coords

def generate():
    print '\nGenerating expressions...'

    # N.B. mesh elements will be sized such that there are mesh_res
    # elements along domain edges in the x-direction, where mesh_res is
    # 5 for mesh 'A', 10 for 'B', etc.  Irregular meshes will try to
    # fill the domain with uniformly sized elements.  This means that
    # the domain should be sized 'nicely' in each dimension if there is
    # to be good convergence on these meshes.
    
    D = (1.0, 1.2, 0.8)         # domain size
    g_mag = 9.8                 # gravity magnitude
    ka = 1.567346939e-9         # absolute permeability
    phi = .4                    # porosity
    
    # regard first phase as air; second as water
    mu = (1.725e-5, 1.e-3)        # viscosity
    rho = (1.284, 1000.);         # density
    K = lambda s: ka * s**2       # relperm

    sh = SolutionHarness(D, g_mag, ka, phi, mu, rho, K)

    # arbitrary p and s2 scales (may affect stability)
    p = 1.
    s2 = 0.01

    # helper 1D functions
    fs = lambda x: 3*(1. - x)*(1.5*x)**2
    fp = lambda x: sin(pi*x)**2

    # 1D
    x = Symbol('x')
    g_dir = (-1)
    s2 = s2*fs(x/D[0])
    p = p*fp(x/D[0])
    ms1d = ManufacturedSolution(1, g_dir, p, s2)
    
    # 2D
    y = Symbol('y')
    g_dir = (-1, 0)
    s2 = s2*fs(y/D[1])
    p = p*fp(y/D[1])
    ms2d = ManufacturedSolution(2, g_dir, p, s2)
    
    # 3D
    z = Symbol('z')
    g_dir = (-1, 0, 0)
    s2 = s2*fs(z/D[2])
    p = p*fp(z/D[2])
    ms3d = ManufacturedSolution(3, g_dir, p, s2)

    # generate expressions for the manufactured solution
    sh.write_dict([ms1d, ms2d, ms3d])

    print 'done.'
