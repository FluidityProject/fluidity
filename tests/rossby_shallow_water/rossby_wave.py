import numpy
import os
from fluidity_tools import stat_parser

ee = 1.0e-2
beta = 1.0

k=numpy.array([1,2])
h=1.0

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

    os.system("spud-set rossby_wave.swml /physical_parameters/gravity/magnitude "+str(1.0/ee))
    os.system("spud-set rossby_wave.swml /physical_parameters/coriolis/beta_plane/f_0 "+str(1.0/ee))
    os.system("spud-set rossby_wave.swml /material_phase[0]/scalar_field::LayerThickness/prognostic/mean_layer_thickness "+str(1.0/ee))
              
    error = numpy.zeros(dx.size)

    binary="../../bin/shallow_water"

    dt=0.01
    for x in enumerate(dx):
        generate_meshfile("square",x[1])
        
        last_error=1.e100
        while True:
            print "Simulation with dx="+`x[1]`+" and dt="+`dt`
            
            os.system("spud-set rossby_wave.swml /timestepping/timestep "+str(dt))
            
            os.system(binary+" rossby_wave.swml")
            
            s=stat_parser("rossby_wave.stat")
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

def convergence(dx):

    error = numpy.zeros(dx.size)

    binary="../../bin/shallow_water"
    project = "rossby_wave"

    dt=0.0125
    for x in enumerate(dx):
        generate_meshfile("square",x[1])

        last_error=1.e100
        while True:
            print "Simulation with dx="+`x[1]`+" and dt="+`dt`
            
            os.system("spud-set "+project+".swml /timestepping/timestep "+str(dt))
            
            os.system(binary+" "+project+".swml")
            
            s=stat_parser(project+".stat")
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

def convergence_constant_dt(dx):

    error = numpy.zeros(dx.size)

    binary="../../bin/shallow_water"

    dt=1.0e-3
    for x in enumerate(dx):
        generate_meshfile("square",x[1])

        last_error=1.e100
        print "Simulation with dx="+`x[1]`+" and dt="+`dt`
        
        os.system("spud-set rossby_wave.swml /timestepping/timestep "+str(dt))
            
        os.system(binary+" rossby_wave.swml")
            
        s=stat_parser("rossby_wave.stat")
        this_error=s["Fluid"]['LayerThicknessError']['l2norm'][-1]
        
        error[x[0]]=this_error
            
    return error

def perp(vec):
    '''Return the perp of a 2-vector.'''

    return numpy.array([-vec[1],vec[0]])

def dispersion(k):
    """Return omega as a function of the wave number k for the two
    dimensional shallow water equations"""

    return beta*k[0]/(1 + k[0]**2 + k[1]**2)

def velocity_solution(k):
    '''Return inertia gravity wave solution to the shallow water
    equations for a given wave number vector k. The result is a
    function of x and t.'''

    omega=dispersion(k)

    def sol(x,t):
        return numpy.real(h*perp([1j*k[0]*numpy.sin(k[1]*x[1]), \
                                     k[1]*numpy.cos(k[1]*x[1])])* \
                              numpy.exp(1j*k[0]*x[0]+omega*t))
    return sol

def height_solution(k):
    '''Return inertia gravity wave solution to the shallow water
    equations for a given wave number vector k. The result is a
    function of x and t.'''

    omega=dispersion(k)
    
    def sol(x,t):
        return numpy.real(h*numpy.sin(k[1]*x[1])* \
                              numpy.exp(1j*(k[0]*x[0]+omega*t)))
    return sol
