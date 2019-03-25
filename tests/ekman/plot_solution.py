#!/usr/bin/env python
from numpy import sqrt,exp, pi, cos, sin
from pylab import find, plot, show, xlabel, ylabel, legend
import vtktools

u=vtktools.vtu('ekman_2.vtu')
xyz=u.GetLocations()

rho_air=1.3
rho_water=1025.34
C_s=1.4e-3
nu_H=100
nu_V=1.4e-2
u10=0
v10=10
f=1.03124e-4
# derived quantities:
BigTau=rho_air/rho_water*C_s*sqrt(u10**2+v10**2)*v10
print('BigTau:', BigTau)
D=pi*sqrt(2*nu_V/f)
print('D:', D)
u0=BigTau*D/sqrt(2.0)/pi/nu_V
print('u0:', u0)

# only look at point (x,y)=(0.0, 0.0)
ind=find( (abs(xyz[:,0]-50)<1.0) & (abs(xyz[:,1]-50)<1.0) )
z=xyz[ind,2]-max(xyz[ind,2])

u_ex=cos(pi/4.0+pi*z/D)
u_ex=u0*exp(pi*z/D)*cos(pi/4.0+pi*z/D)
v_ex=u0*exp(pi*z/D)*sin(pi/4.0+pi*z/D)

uvw=u.GetField('Velocity')

xlabel('U-component')
ylabel('V-component')
plot( u_ex, v_ex , '-o')
plot( uvw[ind,0], uvw[ind,1], '.k')
legend(['analytical solution', 'ICOM result'], loc='upper left')
show()