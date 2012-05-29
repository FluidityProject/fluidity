from math import pi, exp, sin, cos
import pylab
import numpy
import sys

# molecular viscosity
nu=0.01

def decay(t):
   return exp(-2.*nu*t)

def initial_cond(XX):
   '''2D Taylor Green Vortex'''
   x = XX[0]; y = XX[1]
   u = sin(x) * cos(y)
   v = -cos(x) * sin(y)
   return [u,v]

def velocity(XX,t):
   '''analytical solution'''
   x = XX[0]; y = XX[1]
   u = sin(x) * cos(y) * decay(t)
   v = -cos(x) * sin(y) * decay(t)
   return [u,v]

def pressure(XX,t):
   x = XX[0]; y = XX[1];rho=1.
   return rho/4.*(cos(2.*x)+cos(2.*y))*(decay(t))**2.

def exactfilter(field,width):
   print field, width, field[:]*(1.-width/24.+width/576.)
   # 4th order accurate Taylor expansion of Helmholtz filter
   return field[:]*(1.-w/24.+w/576.)

def exactstrain(XX,t):
   x = XX[0]; y = XX[1]
   sij = cos(x)*cos(y)*decay(t)
   return [[sij, 0.0],[0.0, -sij]]

def strainmagnitude(XX,t):
   x = XX[0]; y = XX[1]
   # strain magnitude is sqrt(2SijSij) and is therefore always positive.
   return abs(2.*cos(x)*cos(y)*(decay(t))**2.)

def dissipation_rate():
  # calculate -d(TKE)/dt
  return

def plot_velo(XX,t):
  u=[]
  for i in range(len(t)):
    u.append(velocity(XX,t[i]))

  pylab.figure(1)
  pylab.title('Velocity at point '+str(XX[0])+','+str(XX[1])+','+str(XX[2]))
  pylab.xlabel('Time (s)')
  pylab.ylabel('u')
  pylab.plot(t, u)
  pylab.savefig("analytical_velocity.pdf")
  return

def normalise(arrays):
  for i in arrays:
    norm=i[0]
    i[:]=i[:]/norm
  return

def differentiate(times,arrays):
  print 'performing numerical time-differentiation'
  derivs=[]
  for i in range(len(arrays)):
    temp=[]
    t=numpy.array(times[i])
    result=numpy.array(arrays[i])
    print len(t)
    print len(result)
    if len(t)==len(result):
      for j in range(len(result)-1):
        dt = t[j+1]-t[j]
        # -d/dt(tke)
        temp.append(-(result[j+1]+result[j])/dt)
      derivs.append(temp)
  return derivs

def integrate(times,arrays):
  print 'performing numerical time-integration'
  integrals=[]
  arrays=arrays if type(arrays)==tuple else tuple(arrays)
  times=times if type(times)==tuple else tuple(times)
  for i in range(len(arrays)):
    temp=0.0
    t=numpy.array(times[i])
    result=numpy.array(arrays[i])
    print len(t);print len(result)
    if len(t)==len(result):
      for j in range(len(result)-1):
        dt = t[j+1]-t[j]
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

def plot_tke(times,arrays):

  ##### Plot time series of turbulent kinetic energy dissipation rate (epsilon)
  plot1 = pylab.figure(figsize = (12, 6))
  size = 15
  ax = pylab.subplot(111)
  ax.set_title('Error in simulated TKE')
  ax.set_xlabel('Time (s)')
  ax.set_ylabel('turbulent kinetic energy error')
  ax.plot(times[0], arrays[0], linestyle="dashdot",color='blue') # dns-12
  ax.plot(times[1], arrays[1], linestyle="dashed",color='blue') # dns-24
  ax.plot(times[2], arrays[2], linestyle="solid",color='blue') # dns-48
  ax.plot(times[3], arrays[3], linestyle="dashdot",color='green') # les-12
  ax.plot(times[4], arrays[4], linestyle="dashed",color='green') # les-24
  ax.plot(times[5], arrays[5], linestyle="solid",color='green') # les-48
  ax.plot(times[6], arrays[6], linestyle="dashdot",color='red') # les-aniso-12
  ax.plot(times[7], arrays[7], linestyle="dashed",color='red') # les-aniso-24
  ax.plot(times[8], arrays[8], linestyle="solid",color='red') # les-aniso-48
  ax.plot(times[9], arrays[9], linestyle="dashdot",color='orange') # dynles-12
  ax.plot(times[10], arrays[10], linestyle="dashed",color='orange') # dynles-24
  ax.plot(times[11], arrays[11], linestyle="solid",color='orange') # dynles-48
  ax.plot(times[12], arrays[12], linestyle="dashdot",color='purple') # dynles-aniso-12
  ax.plot(times[13], arrays[13], linestyle="dashed",color='purple') # dynles-aniso-24
  ax.plot(times[14], arrays[14], linestyle="solid",color='purple') # dynles-aniso-48

  pylab.legend(('dns-12','dns-24','dns-48',
'les-12','les-24','les-48','les-aniso-12','les-aniso-24','les-aniso-48',
'dynles-12','dynles-24','dynles-48','dynles-aniso-12','dynles-aniso-24','dynles-aniso-48'),loc='best')
  ltext = pylab.gca().get_legend().get_texts()
  pylab.setp(ltext, fontsize = 10, color = 'b')
  pylab.axis([0,300,0,0.2])

  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  pylab.savefig("vortex_tke_error.pdf")
  return

def plot_dissipation(times,arrays):

  ##### Plot time series of turbulent kinetic energy dissipation rate (epsilon)
  plot1 = pylab.figure(figsize = (12, 8))
  size = 15
  ax = pylab.subplot(111)
  ax.set_title('Simulated vs. analytical dissipation')
  ax.set_xlabel('Time (s)')
  ax.set_ylabel('turbulent kinetic energy dissipation rate -d/dt(tke)')
  print len(times[0][0:-1]), len(arrays[0])
  ax.plot(times[0][0:-1], arrays[0], linestyle="dashdot",color='blue') # dns-12
  ax.plot(times[1][0:-1], arrays[1], linestyle="dashed",color='blue') # dns-24
  ax.plot(times[2][0:-1], arrays[2], linestyle="solid",color='blue') # dns-48
  ax.plot(times[3][0:-1], arrays[3], linestyle="dashdot",color='green') # les-12
  ax.plot(times[4][0:-1], arrays[4], linestyle="dashed",color='green') # les-24
  ax.plot(times[5][0:-1], arrays[5], linestyle="solid",color='green') # les-48
  ax.plot(times[6][0:-1], arrays[6], linestyle="dashdot",color='red') # les-aniso-12
  ax.plot(times[7][0:-1], arrays[7], linestyle="dashed",color='red') # les-aniso-24
  ax.plot(times[8][0:-1], arrays[8], linestyle="solid",color='red') # les-aniso-48
  ax.plot(times[9][0:-1], arrays[9], linestyle="dashdot",color='orange') # dynles-12
  ax.plot(times[10][0:-1], arrays[10], linestyle="dashed",color='orange') # dynles-24
  ax.plot(times[11][0:-1], arrays[11], linestyle="solid",color='orange') # dynles-48
  ax.plot(times[12][0:-1], arrays[12], linestyle="dashdot",color='purple') # dynles-aniso-12
  ax.plot(times[13][0:-1], arrays[13], linestyle="dashed",color='purple') # dynles-aniso-24
  ax.plot(times[14][0:-1], arrays[14], linestyle="solid",color='purple') # dynles-aniso-48
  ax.plot(times[2][0:-1], arrays[-1], linestyle="solid",color='black') # analytical-12

  pylab.legend(('dns-12','dns-24','dns-48',
'les-12','les-24','les-48','les-aniso-12','les-aniso-24','les-aniso-48',
'dynles-12','dynles-24','dynles-48','dynles-aniso-12','dynles-aniso-24','dynles-aniso-48',
'exact'),loc='best')
  ltext = pylab.gca().get_legend().get_texts()
  pylab.setp(ltext, fontsize = 10, color = 'b')
  pylab.axis([0,150,-35,1])

  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(size)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(size)

  pylab.savefig("vortex_dissipation.pdf")
  return

