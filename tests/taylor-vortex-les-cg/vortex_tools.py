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

def plot_dissipation(vort_12,time_12,vort_24,time_24,vort_48,time_48,
vort_les_12,time_les_12,vort_les_24,time_les_24,vort_les_48,time_les_48,
vort_dynles_12,time_dynles_12,vort_dynles_24,time_dynles_24,vort_dynles_48,time_dynles_48,
vort_dynlesa_12,time_dynlesa_12,vort_dynlesa_24,time_dynlesa_24,vort_dynlesa_48,time_dynlesa_48,
vort_dynlesa1_12,time_dynlesa1_12,vort_dynlesa1_24,time_dynlesa1_24,vort_dynlesa1_48,time_dynlesa1_48,
vort_dynlesa1a_12,time_dynlesa1a_12,vort_dynlesa1a_24,time_dynlesa1a_24,
vort_dynlesa1b1_12,time_dynlesa1b1_12,vort_dynlesa1b1_24,time_dynlesa1b1_24,vort_dynlesa1b1_48,time_dynlesa1b1_48,
vort_dynlesn_12,time_dynlesn_12):

  ##### Plot time series of turbulent kinetic energy dissipation rate (epsilon)
  pylab.figure(1)
  pylab.title('Simulated vs. analytical dissipation')
  pylab.xlabel('Time (s)')
  pylab.ylabel('turbulent kinetic energy dissipation rate')
  pylab.plot(time_12, vort_12, linestyle="dashdot",color='blue')
  pylab.plot(time_24, vort_24, linestyle="dashed",color='blue')
  pylab.plot(time_48, vort_48, linestyle="solid",color='blue')
  pylab.plot(time_les_12, vort_les_12, linestyle="dashdot",color='green')
  pylab.plot(time_les_24, vort_les_24, linestyle="dashed",color='green')
  pylab.plot(time_les_48, vort_les_48, linestyle="solid",color='green')
  pylab.plot(time_dynles_12, vort_dynles_12, linestyle="dashdot",color='red')
  pylab.plot(time_dynles_24, vort_dynles_24, linestyle="dashed",color='red')
  pylab.plot(time_dynles_48, vort_dynles_48, linestyle="solid",color='red')
  pylab.plot(time_dynlesa_12, vort_dynlesa_12, linestyle="dashdot",color='purple')
  pylab.plot(time_dynlesa_24, vort_dynlesa_24, linestyle="dashed",color='purple')
  pylab.plot(time_dynlesa_48, vort_dynlesa_48, linestyle="solid",color='purple')
  pylab.plot(time_dynlesa1_12, vort_dynlesa1_12, linestyle="dashdot",color='brown')
  pylab.plot(time_dynlesa1_24, vort_dynlesa1_24, linestyle="dashed",color='brown')
  pylab.plot(time_dynlesa1_48, vort_dynlesa1_48, linestyle="solid",color='brown')
  pylab.plot(time_dynlesa1a_12, vort_dynlesa1a_12, linestyle="dashdot",color='magenta')
  pylab.plot(time_dynlesa1a_24, vort_dynlesa1a_24, linestyle="dashed",color='magenta')
  pylab.plot(time_dynlesa1b1_12, vort_dynlesa1b1_12, linestyle="dashdot",color='turquoise')
  pylab.plot(time_dynlesa1b1_24, vort_dynlesa1b1_24, linestyle="dashed",color='turquoise')
  pylab.plot(time_dynlesa1b1_48, vort_dynlesa1b1_48, linestyle="solid",color='turquoise')
  pylab.plot(time_dynlesn_12, vort_dynlesn_12, linestyle="dashdot",color='yellow')

  #pylab.plot(time_12, dissipation_taylor(time_12), linestyle="dashed",color='')
  pylab.plot(time_dynles_12, dissipation_goldstein(time_dynles_12), linestyle="solid",color='black')
  pylab.legend(('dns_12','dns_24','dns_48','vort_les_12','vort_les_24','vort_les_48',
'dynles-12','dynles-24','dynles-48','dynles-aniso-12','dynles-aniso-24','dynles-aniso-48',
'dynles-a1-12','dynles-a1-24','dynles-a1-48','dynles-aniso-a1-12','dynles-aniso-a1-24','dynles-a1b1-12',
'dynles-a1b1-24','dynles-a1b1-48','dynles-nolimc-12','goldstein'),loc='best')

  ltext = pylab.gca().get_legend().get_texts()
  pylab.setp(ltext, fontsize = 8, color = 'b')
  pylab.axis([0,30,0,0.03])
  #pylab.axis([0.0,time_12[-1],0.0,max(dissipation_goldstein(time_12)[-1],vort_48[-1],vort_les48[-1],vort_dynles48[-1])])
  pylab.savefig("vortex_dissipation.pdf")
  return


