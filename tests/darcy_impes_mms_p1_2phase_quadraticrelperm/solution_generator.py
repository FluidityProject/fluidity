from sympy import Symbol, Function, diff, integrate, sin, cos, pi, exp
import math
from manufactured_solution_generation_tools import \
    ManufacturedSolution, SolutionHarness, generate_coords

def generate(solution_name, check_solution=False):
    
    D = (1.0, 1.2, 0.8)         # domain size [see note 1]
    T = 1.0                     # finish time
    g_mag = 1.                 # gravity magnitude [see note 2]
    ka = 1.567346939e-9         # absolute permeability
    phi = .4                    # porosity

    # arbitrary p and s2 scales [see note 2]
    p1_scale = 1000.
    s2_scale = 0.5

    # regard first phase as air; second as water
    mu = (1.725e-5, 1.e-3)            # viscosity
    rho = (1.284, 1000.)              # density

    # relperms - these do not have placeholders in the options template;
    # any changes made here also need to be made in the template
    sr = (0.2, 0.3)                      # residual sats
    K = (lambda s: ka * (s[0]-sr[0])**2, # quadratic relperm relation
         lambda s: ka * (s[1]-sr[1])**2)
    
    # initialise object containing solution parameters
    sh = SolutionHarness(D, T, g_mag, ka, phi, mu, rho, K,
                         p1_scale, s2_scale, solution_name)

    # initialise our prognostic variables
    p1 = p1_scale
    s2 = s2_scale
    
    # introduce time dependence
    t = Symbol('t')
    p1 = p1 * (3 + cos(pi*t/T))/4
    s2 = s2 * 1/(1 + t/T)

    # helper functions for space dependence
    fs_lon = lambda xi: exp(-xi)
    fs_lat = lambda xi: 3*(1. - xi)*(1.5*xi)**2
    fp_lon = lambda xi: cos(pi*xi)
    fp_lat = lambda xi: sin(pi*xi)
    
    # 1D
    x = Symbol('x')
    g_dir = (1)
    s2 = s2*fs_lon(x/D[0])
    p1 = p1*fp_lon(x/D[0])
    ms1d = ManufacturedSolution(1, g_dir, p1, s2)
    
    # 2D
    y = Symbol('y')
    g_dir = (1, 0)
    s2 = s2*fs_lat(y/D[1])
    p1 = p1*fp_lat(y/D[1])
    ms2d = ManufacturedSolution(2, g_dir, p1, s2)
    
    # 3D
    z = Symbol('z')
    g_dir = (1, 0, 0)
    s2 = s2*fs_lat(z/D[2])
    p1 = p1*fp_lat(z/D[2])
    ms3d = ManufacturedSolution(3, g_dir, p1, s2)

    # generate expressions for the manufactured solution
    sh.write_dict([ms1d, ms2d, ms3d])

    # check that divergence + source term = 0
    if check_solution:
        return sh.reality_check([ms1d, ms2d, ms3d])
    
    # Notes:
    # 
    # [1] mesh elements will be sized such that there are mesh_res
    # elements along domain edges in the x-direction, where mesh_res is
    # 5 for mesh 'A', 10 for 'B', etc.  Irregular meshes will try to
    # fill the domain with uniformly sized elements.  This means that
    # the domain probably needs to be sized 'nicely' in each dimension
    # if there is to be good convergence on these meshes.
    # 
    # [2] make pressure level high enough to have an influence (> 100).
    # High levels of saturation and rho*g_mag/mu have been found to
    # cause numerical instability well below the expected CFL limit.
    # This may be caused by having a highly nonlinear relative
    # permeability term and forcing an unnatural pressure field; it only
    # happens with MMS.

