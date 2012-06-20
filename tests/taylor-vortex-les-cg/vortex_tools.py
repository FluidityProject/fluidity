from math import pi, exp, sin, cos
import pylab
import numpy
import sys

A=1.;B=-1.;C=0.
a=1.;b=1.;c=1.

# A1=-1/8, A2=1/8, A3=0
A1 = b/4*(A**2*b*a**2/(b**2+c**2)+A*B*a)
A2 = c/4*(B**2*c*b**2/(c**2+a**2)+B*C*b)
A3 = a/4*(C**2*a*c**2/(a**2+b**2)+C*A*c)

theta=a**2+b**2+c**2

def initial_cond(XX):
   '''Taylor Green Vortex'''
   x = XX[0]; y = XX[1]; z = XX[2]
   u = A*cos(a*x) * sin(b*y) * cos(c*z)
   v = B*sin(a*x) * cos(b*y) * cos(c*z)
   w = 0.
   return [u,v,w]

def velocity(XX,t):
   '''Taylor's analytical solution. Double-check my cyclic permutations.'''
   x = XX[0]; y = XX[1]; z = XX[2];nu=0.01
   u = A*(1.-theta*nu*t)*cos(a*x) * sin(b*y) * sin(c*z) + A3/a*t*sin(2*a*x)*cos(2*b*y) - A2/a*t*sin(2*a*x)*cos(2*c*z)
   v = B*(1.-theta*nu*t)*cos(b*y) * sin(c*z) * sin(a*x) + A1/b*t*sin(2*b*y)*cos(2*c*z) - A3/b*t*sin(2*b*y)*cos(2*a*x)
   w = C*(1.-theta*nu*t)*cos(c*z) * sin(a*x) *sin(b*y) + A2/c*t*sin(2*c*z)*cos(2*a*x) - A1/c*t*sin(2*c*z)*cos(2*b*y)
   #w = C*(1.-theta*nu*t)*cos(c*z) * sin(a*x) *sin(b*y) + A3/c*t*sin(2*c*z)*cos(2*a*x) - A2/c*t*sin(2*c*z)*cos(2*b*y)
   return [u,v,w]

def secondvelocity(XX,t):
   '''Taylor's second approximation of the analytical solution (special case). Double-check my cyclic permutations.'''
   x = XX[0]; y = XX[1]; z = XX[2];nu=0.01
   d1=A(1.-3.*a**2*nu*t + 9./2.*a**4*nu**2*t**2) + 0.5*A*A1*(0.5*t**2-a**2*nu*t**3)
   d2=-d1
   e1=A1*(t - 7.*a**2*nu*t**2 + 3.*a**4*nu**2*t**3)
   e2=-e1
   b3=-0.5*A*A1*(0.5*t**2-a**2*nu*t**3)
   a1=-15./22.*b3
   a2=-a1
   b1=2*a1
   b2=2*a2
   g1=b1
   g2=b2
   g3=-b3
   L3=4./9.*A1**2/a*t**3
   M3=L3
   N3=L3
   y2=g3*a
   y1=g3*a

   u=d1*cos(a*x)*sin(b*y)*sin(c*z) - e2/a*sin(2.*a*x)*sin(2.*c*z) + a1/3.*cos(3.*a*x)*sin(b*y)*sin(c*z) + g2*cos(a*x)*sin(3.*b*y)*sin(c*z) + b3*cos(a*x)*sin(b*y)*sin(3.*c*z) - y2/a*cos(3.*a*x)*sin(b*y)*sin(3.*c*z) + M3*cos(2.*a*x)*sin(2.*b*y)*sin(4.*c*z)
   print u
   v=0
   w=0
   return [u,v,w]

def alternative_initial_cond(XX):
   '''Taylor Green Vortex'''
   x = XX[0]; y = XX[1]; z = XX[2]
   u = A*sin(a*x) * cos(b*y) * cos(c*z)
   v = B*cos(a*x) * sin(b*y) * cos(c*z)
   w = 0.
   return [u,v,w]

def alternative_velocity(XX,t):
   '''Taylor's analytical solution. Double-check my cyclic permutations.'''
   x = XX[0]; y = XX[1]; z = XX[2];nu=0.01
   u = A*(1.-theta*nu*t)*cos(a*x) * sin(b*y) * sin(c*z) + A3/a*t*sin(2*a*x)*cos(2*b*y) - A2/a*t*sin(2*a*x)*cos(2*c*z)
   v = B*(1.-theta*nu*t)*sin(a*x) * cos(b*y) * sin(c*z) + A1/b*t*sin(2*b*y)*cos(2*c*z) - A3/b*t*sin(2*b*y)*cos(2*a*x)
   w = C*(1.-theta*nu*t)*sin(a*x) * sin(b*y) * cos(c*z) + A2/c*t*sin(2*c*z)*cos(2*a*x) - A1/c*t*sin(2*c*z)*cos(2*b*y)
   return [u,v,w]

###############################################

def plot_velo(XX,t):
  u=[]
  for i in range(len(t)):
    u.append(velocity(XX,t[i]))

  pylab.figure(1)
  pylab.title('Velocity at point '+str(XX[0])+','+str(XX[1])+','+str(XX[2]))
  pylab.xlabel('Time (s)')
  pylab.ylabel('u')
  pylab.plot(t, u)
  pylab.savefig('analytical_velocity.pdf')
  return

def dissipation_taylor(t,nu):
   # See Orszag et al, 1976
   '''Taylor's 1st-order Taylor expansion of mean-square vorticity'''
   omega_t=numpy.zeros(t.size)
   omega_0=nu*0.75
   for i in range(len(t)):
      omega_t[i] = nu*0.75*(1.-6.*t[i]*nu + (5./48.+18.*nu**2.)*t[i]**2. - (5./3.+36.*nu**2.)*nu*t[i]**3.)/omega_0
      #omega_t[i] = nu*0.75*(1.-6.*t[i]*nu + (5./48.+18.*nu**2.)*t[i]**2. - (5./3.+36.*nu**2.)*nu*t[i]**3. + (50./99.64+1835./9.16*nu**2.+54.*nu**4.)*t[i]**4. - (361./44.32+761./12.*nu**2.+324./5.*nu**4.)*nu*t[i]**5.)/omega_0
   return omega_t

# beware of difference between goldstein and normalised goldstein:
# these will give different convergence rates.
# normalised is for plotting, absolute is for convergence.

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
  pylab.savefig('analytical_velocity.pdf')
  return

def print_dims(arrays):
  for i in arrays:
    print len(i)
  return

def normalise(arrays):
  for i in arrays:
    norm=i[0]
    i[:]=i[:]/norm
  return

def integrate(tarrays,arrays):
  print 'performing numerical time-integration'
  integrals=[]

  print type(arrays)
  temp=arrays if type(arrays)==tuple else tuple(arrays)

  for i in range(len(arrays)):
    temp=0.0
    result=numpy.array(arrays[i])
    for j in range(len(result)-1):
      dt=tarrays[i][j+1]-tarrays[i][j]
      # trapezium rule of numerical integration
      temp += (result[j+1]+result[j])*dt/2.0
    integrals.append(temp)

  return integrals

def convergence_rate(integrals,benchmark):
  convrates=[]
  for i in range(len(integrals)-1):
    conv = abs((integrals[i]-benchmark)/(integrals[i+1]-benchmark))
    convrates.append(conv)
  return convrates

def plot_dissipation(tarrays,varrays):

  ##### Plot time series of turbulent kinetic energy dissipation rate (epsilon)
  plot1 = pylab.figure(figsize = (12, 10))
  size = 15
  ax = pylab.subplot(111)
  ax.set_title('Simulated vs. analytical dissipation')
  ax.set_xlabel('Time (s)')
  ax.set_ylabel('turbulent kinetic energy dissipation rate')

  ax.plot(tarrays[0], varrays[0], linestyle='dotted',color='blue') # dns-periodic-12
  ax.plot(tarrays[1], varrays[1], linestyle='dashdot',color='blue') # dns-periodic-24
  ax.plot(tarrays[2], varrays[2], linestyle='dashed',color='blue') # dns-periodic-48
  ax.plot(tarrays[3], varrays[3], linestyle='solid',color='blue') # dns-periodic-64
  ax.plot(tarrays[4], varrays[4], linestyle='dotted',color='green') # dynles-periodic-12
  ax.plot(tarrays[5], varrays[5], linestyle='dashdot',color='green') # dynles-periodic-24
  ax.plot(tarrays[6], varrays[6], linestyle='dashed',color='green') # dynles-periodic-48
  ax.plot(tarrays[7], varrays[7], linestyle='solid',color='green') # dynles-periodic-64
  ax.plot(tarrays[8], varrays[8], linestyle='dotted',color='red') # dynles-aniso-periodic-12
  ax.plot(tarrays[9], varrays[9], linestyle='dashdot',color='red') # dynles-aniso-periodic-24
  ax.plot(tarrays[10], varrays[10], linestyle='solid',color='purple') # analytical-12

  ax.plot(tarrays[0], dissipation_goldstein_norm(tarrays[0],0.005), linestyle='none',color='black',marker='s',markerfacecolor='none') # Goldstein-0.005
  ax.plot(tarrays[0], dissipation_goldstein_norm(tarrays[0],0.01), linestyle='none',color='black',marker='x') # Goldstein-0.01
  ax.plot(tarrays[0], dissipation_goldstein_norm(tarrays[0],0.01333), linestyle='none',color='black',marker='o',markerfacecolor='none') # Goldstein-0.01333
  ax.plot(tarrays[0], dissipation_goldstein_norm(tarrays[0],0.015), linestyle='none',color='black',marker='+') # Goldstein-0.015
  #ax.plot(tarrays[], dissipation_taylor(tarrays[],1./50.), linestyle='dashdot',color='cyan')

  pylab.legend(('vort_p_12','vort_p_24','vort_p_48','vort_p_64',
'vort_dynlesp_12','vort_dynlesp_24','vort_dynlesp_48','vort_dynlesp_64',
'vort_dynlesap_12','vort_dynlesap_24',
'vort_analytp_12','vort_gold_0.005','vort_gold_0.01','vort_gold_0.0133','vort_gold_0.015'),loc='best')

  ltext = pylab.gca().get_legend().get_texts()
  pylab.setp(ltext, fontsize = 10, color = 'b')
  pylab.axis([0,20,0,6])

  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  pylab.savefig('vortex_dissipation.pdf')
  return

