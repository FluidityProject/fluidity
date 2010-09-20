import os
from fluidity_tools import stat_parser
from sympy import *
from numpy import array,max,abs

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

    file(name+".geo",'w').write(
        meshtemplate.replace('<dx>',str(1./layers)
                 ).replace('<layers>',str(layers)))

    os.system("gmsh -2 "+name+".geo")


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


def analytic_pressure_solution(forcing):
    
    u=analytic_solution(forcing)

    return integrate(-coriolis*u+forcing((0,Symbol('x')))[1], Symbol('x'))


def pressure_solution(forcing):
    '''Return a function which is the solution to:
    dp/dx = f x u The constant of integration is set to 0.'''

    def sol(sx):
        return analytic_pressure_solution(forcing).subs(Symbol('x'),sx[1])

    return sol


def plot_theory():
    '''Produce a plot showing the forcing, analytic velocity solution and
    analytic pressure solution'''
    from pylab import \
    plot,figure,quiver,frange,subplot,xticks,yticks,axis,xlabel,ylabel, \
    subplots_adjust 


    figure()

    y=frange(0.0,1,0.05)

    psol=pressure_solution(forcing)

    usol=solution(forcing)

    v=0*y

    x=0*y

    us=array([float(usol(pos)) for pos in zip(x,y)])

    ps=array([float(psol(pos)) for pos in zip(x,y)])

    uf=array([forcing(pos) for pos in zip(x,y)])[:,0]

    subplots_adjust(wspace=0.25)
    subplot(1,3,1)

    quiver(x[1:-1],y[1:-1],uf[1:-1],v[1:-1], scale=1)
    plot(uf,y)
    xticks([0,0.5,1],map(str,[0,0.5,1]))
    yticks([ 0 ,  0.2,  0.4,  0.6,  0.8,  1 ],map(str,[ 0 ,  0.2,  0.4,  0.6,  0.8,  1 ]))
    ylabel("y")
    xlabel("u source")

    subplot(1,3,2)
    plot(us,y)
    quiver(x[1:-1],y[1:-1],us[1:-1],v[1:-1], scale=.03)
    xticks([0,0.01,0.02,0.03],map(str,[0,0.01,0.02,0.03]))
    yticks([])
    xlabel("u solution")

    subplot(1,3,3)
    plot(ps,y)
    xticks([-0.02,-0.01,0],map(str,[-0.02,-0.01,0]))
    yticks([])
    xlabel("p solution")
    
    return uf,us,ps

def plot_results(dx, error):
    '''plot_results(error)

    Produce a plot of the actual errors provided in the argument
    "error". Error should be a two column matrix with the first column being
    the velocity error and the second column the pressure error.
    '''
    from pylab import \
    figure,xticks,yticks,axis,xlabel,ylabel,loglog,legend,title

    figure()

    loglog(dx,error)
    loglog(dx,0.03*dx**2)
    yticks(yticks()[0], map(lambda x: "%3.1e"%x, yticks()[0]))
    xticks(xticks()[0], map(lambda x: "%3.1e"%x, xticks()[0]))

    xlabel("dx")

    title("Convergence of the rotating channel")

    legend(("u error","p error","O(dx^2)"))


def retrieve_results(layers):
    '''retrieve_results(layers)

    For each layer count in the sequence layers, retrieve the velocity and
    pressure error from the simulation results in appropriate channel-n
    directory.

    The first column of the result is the l2 norm of the error in the u
    component of velocity. The second is the l2 norm in the pressure.
    '''
    from numpy import zeros

    error=zeros((len(layers),2))

    for i,layer in enumerate(layers):
        s=stat_parser("channel-%d/rotating_channel.stat"%layer)

        error[i,0]=s["Water"]['AnalyticUVelocitySolutionError']['l2norm'][-1]
        error[i,1]=s["Water"]['AnalyticPressureSolutionError']['l2norm'][-1]
    
    return error
