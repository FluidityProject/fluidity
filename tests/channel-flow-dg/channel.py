import os
from fluidity_tools import stat_parser
from sympy import *

meshtemplate='''
Point(1) = {0, 0, 0, <dx>};
Extrude {0, 1, 0} {
  Point{1};Layers{<layers>};
}
Point(3) = {1, 0, 0, <dx>};
Extrude {0, 1, 0} {
  Point{3};Layers{<layers>};
}
Line(3)={1,3};
Line(4)={2,4};
Line Loop(5) = {4, -2, -3, 1};
Plane Surface(6) = {5};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {4, 3};
Physical Surface(1) = {6};
'''

def generate_meshfile(name,layers):


    file(name,'w').write(
        meshtemplate.replace('<dx>',str(1./layers)
                 ).replace('<layers>',str(layers)))


def run_test(layers, binary):
    '''run_test(layers, binary)

    Run a single test of the channel problem. Layers is the number of mesh
    points in the cross-channel direction. The mesh is unstructured and
    isotropic. binary is a string containing the fluidity command to run.
    The return value is the error in u and p at the end of the simulation.'''

    generate_meshfile("channel.geo",layers)

    os.system("gmsh -2 channel.geo")
    os.system("../../scripts/gmsh2triangle --2d channel.msh")

    os.system(binary+" channel.flml")

    s=stat_parser("channel-flow-dg.stat")
    return (s["Water"]['AnalyticUVelocitySolutionError']['l2norm'][-1],
        s["Water"]['AnalyticPressureSolutionError']['l2norm'][-1])

def forcing(X):
    '''Forcing function. Must be an analytic function of X[1] only'''

    return (X[1]**3,0)

# Numeric verision of the forcing function, for efficiency.
def numeric_forcing(X):
    '''Forcing function. Must be an analytic function of X[1] only'''
    from math import sin, pi

    return (X[1]**3,0)

#Viscosity
mu=1.0

#Note that because Coriolis can't be set from Python, the user has to ensure
#that this matches what it in the flml.
coriolis=1.0
#coriolis=0.0

def analytic_solution(forcing):
    '''Solve the ode d^2u/dx^2 = F/mu subject to u(0)=0, u(1)=0'''

    x=Symbol('x')
    # Constants of integration.
    c1=Symbol('c_1')
    c2=Symbol('c_2')
    
    general=integrate(integrate(-forcing((0,x))[0]/mu,x)+c1,x)+c2

    constants = solve((Eq(general.subs(x,0),0),
                  Eq(general.subs(x,1),0)), c1,c2)

    specific=general.subs(constants)

    return specific

def solution(forcing):
    '''Return a function which is the solution to:
    d^2u/dx^2 = F/mu subject to u(0)=0, u(1)=0'''

    def sol(sx):
        return analytic_solution(forcing).subs(Symbol('x'),sx[1])

    return sol

absorption=0.5

def analytic_solution(forcing):
    '''Return the steady state of the ode du/dt = F - Au'''

    x=Symbol('x')

    u=forcing((0.0,x))[0]/absorption

    return u

def solution(forcing):
    '''Return a function which is the solution to:
    ode du/dt = F - Au'''

    def sol(sx):
        return analytic_solution(forcing).subs(Symbol('x'),sx[1])

    return sol


def analytic_pressure_solution(forcing):
    
    u=analytic_solution(forcing)

    return integrate(-coriolis*u+forcing((0,Symbol('x')))[1], Symbol('x'))

def pressure_solution(forcing):
    '''Return a function which is the solution to:
    dp/dx = f x u The constant of integration is set to 0.'''

    def sol(sx):
        return analytic_pressure_solution(forcing).subs(Symbol('x'),sx[1])

    return sol



