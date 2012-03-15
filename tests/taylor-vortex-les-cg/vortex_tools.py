from math import pi, exp, sin, cos
import pylab
import numpy
import sys

A=1.;B=-1.;C=0.
a=1.;b=1.;c=1.

A1 = b/4*(A**2*b*(a**2/(b**2+c**2))+A*B*a)
A2 = c/4*(A**2*c*(b**2/(c**2+a**2))+A*B*b)
A3 = a/4*(A**2*a*(c**2/(a**2+b**2))+A*B*c)

theta=a**2+b**2+c**2

def initial_cond(XX):
   '''Taylor Green Vortex'''
   x = XX[0]; y = XX[1]; z = XX[2]
   u = A*cos(a*x) * sin(b*y) * cos(c*z)
   v = B*sin(a*x) * cos(b*y) * cos(c*z)
   w = C*sin(a*x) * sin(b*y) * cos(c*z)
   return [u,v,w]

def velocity(XX,t):
   '''Taylor's analytical solution. Double-check my cyclic permutations.'''
   x = XX[0]; y = XX[1]; z = XX[2];nu=0.01
   u = A*(1.-theta*nu*t)*cos(a*x) * sin(b*y) * sin(c*z) + A3/a*t*sin(2*a*x)*cos(2*b*y) - A2/a*t*sin(2*a*x)*cos(2*c*z)
   v = B*(1.-theta*nu*t)*sin(a*x) * cos(b*y) * sin(c*z) + A1/b*t*sin(2*b*y)*cos(2*c*z) - A3/b*t*sin(2*b*y)*cos(2*a*x)
   w = C*(1.-theta*nu*t)*sin(a*x) * sin(b*y) * cos(c*z) + A2/c*t*sin(2*c*z)*cos(2*a*x) - A1/c*t*sin(2*c*z)*cos(2*b*y)
   return [u,v,w]


def dissipation_taylor(t,nu):
   # See Orszag et al, 1976
   '''Taylor's 1st-order Taylor expansion of mean-square vorticity'''
   omega_t=numpy.zeros(t.size)
   omega_0=nu*0.75
   for i in range(len(t)):
      omega_t[i] = nu*0.75*(1.-6.*t[i]*nu + (5./48.+18.*nu**2.)*t[i]**2. - (5./3.+36.*nu**2.)*nu*t[i]**3.)/omega_0
      #omega_t[i] = nu*0.75*(1.-6.*t[i]*nu + (5./48.+18.*nu**2.)*t[i]**2. - (5./3.+36.*nu**2.)*nu*t[i]**3. + (50./99.64+1835./9.16*nu**2.+54.*nu**4.)*t[i]**4. - (361./44.32+761./12.*nu**2.+324./5.*nu**4.)*nu*t[i]**5.)/omega_0
   return omega_t

def dissipation_goldstein(t,nu):
   # See Orszag et al, 1976
   '''Goldstein's expansion of mean-square vorticity'''
   omega_g=numpy.zeros(t.size)
   for i in range(len(t)):
      #omega_g[i] = nu*0.75*(exp(-6.*t[i]*nu) - 1./384./nu**2*(exp(-6.*t[i]*nu) - 20.*exp(-12.*t[i]*nu) + 35.*exp(-14.*t[i]*nu) - 16.*exp(-16.*t[i]*nu)) + 1./nu**4*(4105./635830272.*exp(-6.*t[i]*nu)))
      omega_g[i] = nu*0.75*(exp(-6.*t[i]*nu) - 1./384./nu**2.*(exp(-6.*t[i]*nu) - 20.*exp(-12.*t[i]*nu) + 35.*exp(-14.*t[i]*nu) - 16.*exp(-16.*t[i]*nu)))
   return omega_g

def dissipation_goldstein_norm(t,nu):
   # See Orszag et al, 1976
   '''Goldstein's expansion of mean-square vorticity'''
   omega_g_norm=numpy.zeros(t.size)
   omega_0=nu*0.75*(1 - 1./384./nu**2*(1. - 20. + 35. - 16.))
   for i in range(len(t)):
      #omega_g_norm[i] = nu*0.75*(exp(-6.*t[i]*nu) - 1./384./nu**2*(exp(-6.*t[i]*nu) - 20.*exp(-12.*t[i]*nu) + 35.*exp(-14.*t[i]*nu) - 16.*exp(-16.*t[i]*nu)) + 1./nu**4*(4105./635830272.*exp(-6.*t[i]*nu)))/omega_0
      omega_g_norm[i] = nu*0.75*(exp(-6.*t[i]*nu) - 1./384./nu**2.*(exp(-6.*t[i]*nu) - 20.*exp(-12.*t[i]*nu) + 35.*exp(-14.*t[i]*nu) - 16.*exp(-16.*t[i]*nu)))/omega_0
   return omega_g_norm

def plot_velo(XX,t):
  vel=[]
  for i in range(len(t)):
    vel.append(velocity(XX,t[i]))
  vel=numpy.array(vel)

  pylab.figure(2)
  pylab.title('Velocity at point '+str(XX[0])+','+str(XX[1])+','+str(XX[2]))
  pylab.xlabel('Time (s)')
  pylab.ylabel('u')
  pylab.plot(t, vel[:,0])
  pylab.plot(t, vel[:,1])
  pylab.plot(t, vel[:,2])
  pylab.savefig("analytical_velocity.pdf")
  return

def plot_dissipation(arrays):
  ##### Plot time series of turbulent kinetic energy dissipation rate (epsilon)
  pylab.figure(1)
  pylab.title('Simulated vs. analytical dissipation')
  pylab.xlabel('Time (s)')
  pylab.ylabel('turbulent kinetic energy dissipation rate')
  pylab.plot(arrays[0], arrays[1], linestyle="dashdot",color='blue') # dns-12
  pylab.plot(arrays[0], arrays[2], linestyle="dashed",color='blue') # dns-24
  pylab.plot(arrays[0], arrays[3], linestyle="solid",color='blue') # dns-48
  pylab.plot(arrays[0], arrays[4], linestyle="dashdot",color='green') # les-12
  pylab.plot(arrays[0], arrays[5], linestyle="dashed",color='green') # les-24
  pylab.plot(arrays[0], arrays[6], linestyle="solid",color='green') # les-48
  pylab.plot(arrays[0], arrays[7], linestyle="dashdot",color='red') # dynles-12
  pylab.plot(arrays[0], arrays[8], linestyle="dashed",color='red') # dynles-24
  pylab.plot(arrays[0], arrays[9], linestyle="solid",color='red') # dynles-48
  pylab.plot(arrays[0], arrays[10], linestyle="dashdot",color='purple') # dynlesa-12
  pylab.plot(arrays[0], arrays[11], linestyle="dashed",color='purple') # dynlesa-24
  pylab.plot(arrays[0], arrays[12], linestyle="solid",color='purple') # dynlesa-48
  pylab.plot(arrays[0], arrays[13], linestyle="dashdot",color='brown') # dynlesa1-12
  #pylab.plot(arrays[0], arrays[14], linestyle="dashed",color='brown') # dynlesa1-24
  pylab.plot(arrays[0], arrays[15], linestyle="dashdot",color='magenta') # dynlesa1a-12
  pylab.plot(arrays[0], arrays[16], linestyle="dashed",color='magenta') # dynlesa1a-24
  pylab.plot(arrays[0], arrays[17], linestyle="dashdot",color='turquoise') # dynlesa1b1-12
  #pylab.plot(arrays[0], arrays[18], linestyle="dashed",color='turquoise') # dynlesa1b1-24
  pylab.plot(arrays[0], arrays[19], linestyle="dashdot",color='yellow') # dynlesn-12

  pylab.plot(arrays[0], arrays[20], linestyle="solid",color='pink') # analytical-12

  pylab.plot(arrays[0], dissipation_goldstein_norm(arrays[0],1./50.), linestyle="solid",color='black')
  pylab.plot(arrays[0], dissipation_goldstein_norm(arrays[0],1./75.), linestyle="solid",color='black')
  pylab.plot(arrays[0], dissipation_goldstein_norm(arrays[0],1./100.), linestyle="solid",color='black')

  #pylab.plot(arrays[0], dissipation_taylor(arrays[0],1./50.), linestyle="dashdot",color='cyan')
  #pylab.plot(arrays[0], dissipation_taylor(arrays[0],1./75.), linestyle="dashed",color='cyan')
  #pylab.plot(arrays[0], dissipation_taylor(arrays[0],1./100.), linestyle="solid",color='cyan')

  pylab.legend(('dns_12','dns_24','dns_48','vort_les_12','vort_les_24','vort_les_48',
'dynles-12','dynles-24','dynles-48','dynles-aniso-12','dynles-aniso-24','dynles-aniso-48',
'dynles-a1-12','dynles-aniso-a1-12','dynles-aniso-a1-24',
'dynles-a1b1-12','dynles-nolimc-12','analytical-12',
'goldstein-0.02','goldstein-0.0133','goldstein-0.01'),loc='best')
  #,'taylor-0.1','taylor-0.05','taylor-0.01'
  ltext = pylab.gca().get_legend().get_texts()
  pylab.setp(ltext, fontsize = 5, color = 'b')
  pylab.axis([0,30,0,3])
  #pylab.axis([0.0,arrays[0][-1],0.0,max(dissipation_goldstein(arrays[0])[-1],vort_48[-1],vort_les48[-1],vort_dynles48[-1])])
  pylab.savefig("vortex_dissipation.pdf")
  pylab.close()
  return

def print_dims(arrays):
  for i in arrays:
    print len(i)
  return

def normalise(arrays):
  count=0
  for i in arrays:
    if not count==0:
      norm=i[0]
      i[:]=i[:]/norm
    count+=1
  return

