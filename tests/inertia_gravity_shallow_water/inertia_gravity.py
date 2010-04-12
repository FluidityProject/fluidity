import numpy
import os
from fluidity_tools import stat_parser

meshtemplate='''
Point(1) = {0, 0, 0};
Extrude {0, <size>, 0} {
  Point{1}; Layers{<layers>};
}
Extrude {<size>, 0, 0} {
  Line{1}; Layers{<layers>};
}
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Surface(1) = {5};
'''
def generate_meshfile(name,layers):


    file(name+".geo",'w').write(
        meshtemplate.replace('<layers>',str(layers)).replace(
            '<size>',str(2*numpy.pi)))

    os.system("gmsh -2 "+name+".geo")
    os.system("../../scripts/gmsh2triangle --2d "+name+".msh")

def convergence(dx):

    error = numpy.zeros(dx.size)

    binary="../../bin/shallow_water"

    dt=0.0125
    for x in enumerate(dx):
        generate_meshfile("square",x[1])

        last_error=1.e100
        while True:
            print "Simulation with dx="+`x[1]`+" and dt="+`dt`
            
            os.system("spud-set inertia_gravity.swml /timestepping/timestep "+str(dt))
            
            os.system(binary+" inertia_gravity.swml")
            
            s=stat_parser("inertia_gravity.stat")
            this_error=s["Fluid"]['LayerThicknessError']['l2norm'][-1]

            if numpy.abs((this_error-last_error)/this_error) <0.01:
                # Set dt to twice the current level so that the next stage
                # has some chance to converge at the same dt.
                dt=dt*2

                print "converged! dx="+`x[1]`+" and dt="+`dt`
                print "this error="+`this_error`+" last_error="+`last_error`

                error[x[0]]=this_error
                break

            last_error=this_error
            dt=dt*.5
            
    return error



g=1
H=1
c=(g*H)**0.5
f=1

k=numpy.array([1,2])
phi=0.01

def perp(vec):
    '''Return the perp of a 2-vector.'''

    return numpy.array([-vec[1],vec[0]])

def dispersion(k):
    """Return omega as a function of the wave number k for the two
    dimensional shallow water equations"""

    return numpy.sqrt((f**2+c**2*(k[0]**2+k[1]**2)))
    

def velocity_solution(k):
    '''Return inertia gravity wave solution to the shallow water
    equations for a given wave number vector k. The result is a
    function of x and t.'''

    omega=dispersion(k)

    def sol(x,t):
        return numpy.real(
            (1j*k+perp(k)*f/omega)*phi*
            numpy.exp(1j*(numpy.dot(k,x)-omega*t)))
    return sol

def height_solution(k):
    '''Return inertia gravity wave solution to the shallow water
    equations for a given wave number vector k. The result is a
    function of x and t.'''

    omega=dispersion(k)

    def sol(x,t):
        return H*(1.+numpy.real(
            (1j*numpy.dot(k,k)/omega)*phi*
            numpy.exp(1j*(numpy.dot(k,x)-omega*t))))
    return sol
