from math import pi, exp, sin, cos
import pylab
import numpy

A=1.;B=-1.;C=0.
a=1.;b=1.;c=1.

A1 = b/4*(A**2*b*(a**2/(b**2+c**2))+A*B*a)
A2 = c/4*(A**2*c*(b**2/(c**2+a**2))+A*B*b)
A3 = a/4*(A**2*a*(c**2/(a**2+b**2))+A*B*c)

theta=a**2+b**2+c**2
nu=0.01

def initial_cond(XX):
   '''Taylor Green Vortex'''
   x = XX[0]; y = XX[1]; z = XX[2]
   u = A*cos(a*x) * sin(b*y) * sin(c*z)
   v = B*sin(a*x) * cos(b*y) * sin(c*z)
   w = C*sin(a*x) * sin(b*y) * cos(c*z)
   return [u,v,w]

def velocity(XX,t):
   '''Taylor's analytical solution. Double-check my cyclic permutations.'''
   x = XX[0]; y = XX[1]; z = XX[2]
   print 'A1, A2, A3: ', A1, A2, A3
   u = A*(1.-theta*nu*t)*cos(a*x) * sin(b*y) * sin(c*z) + A3/a*t*sin(2*a*x)*cos(2*b*y) - A2/a*t*sin(2*a*x)*cos(2*c*z)
   v = B*(1.-theta*nu*t)*sin(a*x) * cos(b*y) * sin(c*z) + A1/b*t*sin(2*b*y)*cos(2*c*z) - A3/b*t*sin(2*b*y)*cos(2*a*x)
   w = C*(1.-theta*nu*t)*sin(a*x) * sin(b*y) * cos(c*z) + A2/c*t*sin(2*c*z)*cos(2*a*x) - A1/c*t*sin(2*c*z)*cos(2*b*y)
   return [u,v,w]

def dissipation_taylor(t):
   '''Taylor's 1st-order Taylor expansion of mean-square vorticity'''
   omega_t=numpy.zeros(t.size)
   for i in range(len(t)):
      omega_t[i] = nu*0.75*(1.-6*t[i]*nu + (5./48+18.*nu**2)*t[i]**2 - (5./3+36.*nu**2)*nu*t[i]**3)
   return omega_t

def dissipation_goldstein(t):
   '''Goldstein's expansion of mean-square vorticity'''
   omega_g=numpy.zeros(t.size)
   for i in range(len(t)):
      omega_g[i] = nu*0.75*(exp(-6.*t[i]*nu) - 1./384./nu**2*(exp(-6.*t[i]*nu) - 20.*exp(-12.*t[i]*nu) + 35.*exp(-14.*t[i]*nu) - 16.*exp(-16.*t[i]*nu)))
   return omega_g

def plot_dissipation(vort_12,time_12,vort_24,time_24,vort_48,time_48,vort_les12,time_les12):
  ##### Plot time series of turbulent kinetic energy dissipation rate (epsilon)
  # vorticity L2 norm = sqrt(vort_x**2+vort_y**2+vort_z**2)
  # volume-average dissipation rate = nu*(vort_x**2+vort_y**2+vort_z**2)/volume
  vort_12 = nu*vort_12**2/(2*pi)**3
  vort_24 = nu*vort_24**2/(2*pi)**3
  vort_48 = nu*vort_48**2/(2*pi)**3
  vort_les12 = nu*vort_les12**2/(2*pi)**3

  pylab.figure(1)
  pylab.title('Simulated vs. analytical dissipation')
  pylab.xlabel('Time (s)')
  pylab.ylabel('turbulent kinetic energy dissipation rate')
  pylab.plot(time_12, vort_12, linestyle="solid")
  pylab.plot(time_24, vort_24, linestyle="solid")
  pylab.plot(time_48, vort_48, linestyle="solid")
  pylab.plot(time_les12, vort_les12, linestyle="solid")
  pylab.plot(time_12, dissipation_taylor(time_12), linestyle="dashed")
  pylab.plot(time_12, dissipation_goldstein(time_12), linestyle="dashed")
  pylab.legend(("Fluidity-12","Fluidity-24","Fluidity-48","Fluidity-les-12","Taylor approx","Goldstein approx"), loc="upper left")
  pylab.axis([0.0,time_12[-1],0.0,max(dissipation_goldstein(time_12)[-1],vort_12[-1],vort_24[-1],vort_48[-1])])
  pylab.savefig("vortex_dissipation.pdf")
  return


