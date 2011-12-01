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
   # See Orszag et al, 1976
   '''Taylor's 1st-order Taylor expansion of mean-square vorticity'''
   omega_t=numpy.zeros(t.size)
   for i in range(len(t)):
      omega_t[i] = nu*0.75*(1.-6*t[i]*nu + (5./48+18.*nu**2)*t[i]**2 - (5./3+36.*nu**2)*nu*t[i]**3)
      #omega_t[i] = nu*0.75*(1.-6*t[i]*nu + (5./48+18.*nu**2)*t[i]**2 - (5./3+36.*nu**2)*nu*t[i]**3 + (50/99.64+1835/9.16*nu**2+54.*nu**4)*t[i]**4 - (361/44.32+761/12.*nu**2+324/5.*nu**4)*nu*t[i]**5)
   return omega_t

def dissipation_goldstein(t):
   # See Orszag et al, 1976
   '''Goldstein's expansion of mean-square vorticity'''
   omega_g=numpy.zeros(t.size)
   for i in range(len(t)):
      #omega_g[i] = nu*0.75*(exp(-6.*t[i]*nu) - 1./384./nu**2*(exp(-6.*t[i]*nu) - 20.*exp(-12.*t[i]*nu) + 35.*exp(-14.*t[i]*nu) - 16.*exp(-16.*t[i]*nu)))
      omega_g[i] = nu*0.75*(exp(-6.*t[i]*nu) - 1./384./nu**2*(exp(-6.*t[i]*nu) - 20.*exp(-12.*t[i]*nu) + 35.*exp(-14.*t[i]*nu) - 16.*exp(-16.*t[i]*nu)) + 1./nu**4*(4105/635830272*exp(-6.*t[i]*nu)))
   return omega_g

def plot_dissipation(vort_12,time_12,vort_24,time_24,vort_48,time_48,vort_les12,time_les12,vort_les24,time_les24,vort_les48,time_les48,
vort_dynles12,time_dynles12,vort_dynles24,time_dynles24,vort_dynlesa24,time_dynlesa24,vort_dynles48,time_dynles48,
vort_dynlesa48,time_dynlesa48):
  ##### Plot time series of turbulent kinetic energy dissipation rate (epsilon)

  pylab.figure(1)
  pylab.title('Simulated vs. analytical dissipation')
  pylab.xlabel('Time (s)')
  pylab.ylabel('turbulent kinetic energy dissipation rate')
  pylab.plot(time_12, vort_12, linestyle="dashed",color='blue')
  pylab.plot(time_24, vort_24, linestyle="dashed",color='green')
  pylab.plot(time_48, vort_48, linestyle="dashed",color='red')
  pylab.plot(time_les12, vort_les12, linestyle="dashdot",color='blue')
  pylab.plot(time_les24, vort_les24, linestyle="dashdot",color='green')
  pylab.plot(time_les48, vort_les48, linestyle="dashdot",color='red')
  pylab.plot(time_dynles12, vort_dynles12, linestyle="solid",color='blue')
  pylab.plot(time_dynles24, vort_dynles24, linestyle="solid",color='green')
  pylab.plot(time_dynlesa24, vort_dynlesa24, linestyle="solid",color='pink')
  pylab.plot(time_dynles48, vort_dynles48, linestyle="solid",color='red')
  pylab.plot(time_dynlesa48, vort_dynlesa48, linestyle="solid",color='purple')
  #pylab.plot(time_12, dissipation_taylor(time_12), linestyle="dashed")
  pylab.plot(time_12, dissipation_goldstein(time_12), linestyle="solid",color='black')
  pylab.legend(("12","24","48","LES-12","LES-24","LES-48","DYNLES-12","DYNLES-24","DYNLESA-24","DYNLES-48","DYNLESA-48","Goldstein"), loc="best")
  pylab.axis([0.0,time_12[-1],0.0,max(dissipation_goldstein(time_12)[-1],vort_48[-1],vort_les48[-1],vort_dynles48[-1])])
  pylab.savefig("vortex_dissipation.pdf")
  return


