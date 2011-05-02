from math import exp, pi, log

def velocity(XX):
   '''Velocity'''
   x = XX[0]; y = XX[1]
   sigma=4; eta = sigma*y/x
   # erf(eta)
   u = 2./pi**0.5*(exp(-eta**2)-1.)
   v = 1./sigma/pi**0.5*(1.-exp(-eta**2))
   return [u, v]

def pressure(XX):
   '''Pressure'''
   x = XX[0]; y = XX[1]
   rho = 1.0; uref=1.0
   # Cp = p/rho/uref**2
   p = 0.5*rho*uref**2*log(2.*x-x**2+1./4.)*log(4.*y**3-3.*y**2+5./4.)
   return p

def tke(XX):
   '''Turbulent kinetic energy'''
   x = XX[0]; y = XX[1]
   sigma=4; sigmav = 2.5*sigma; kmax = 0.01; etav = sigmav*y/x
   return kmax*etav**2*exp(1.-etav**2)

def eps(XX):
   '''Turbulent dissipation'''
   x = XX[0]; y = XX[1]
   sigma=4; sigmav = 2.5*sigma; kmax = 0.01; etav = sigmav*y/x
   uref=1.0; L=0.5; nu = uref*L/1.e6; numax = 1000.*nu
   return 0.36*kmax**2/numax*exp(-etav**2)

def grad_u(XX):
   '''grad u'''
   x = XX[0]; y = XX[1]
   sigma=4; eta = sigma*y/x
   u_x = -2./pi**0.5*eta/x*exp(-eta**2)
   u_y = 2./pi**0.5*sigma/x*exp(-eta**2)
   return [u_x, u_y]

def grad_v(XX):
   '''grad v'''
   x = XX[0]; y = XX[1]
   sigma=4; eta = sigma*y/x
   v_x = -2./pi**0.5*eta*y/x**2*exp(-eta**2)
   v_y = 2./pi**0.5*eta/x*exp(-eta**2)
   return [v_x, v_y]

def grad_tke(XX):
   '''grad turbulent kinetic energy'''
   x = XX[0]; y = XX[1]
   sigma=4; sigmav = 2.5*sigma; kmax = 0.01; etav = sigmav*y/x
   k_x = 2.*kmax*sigmav**2*y**2/x**3*(etav**2-1.)*exp(-etav**2)
   k_y = 2.*kmax*sigmav**2*y/x**2*(1.-etav**2)*exp(-etav**2)
   return [k_x, k_y]

def grad_eps(XX):
   '''grad turbulent dissipation'''
   x = XX[0]; y = XX[1]
   sigma=4; sigmav = 2.5*sigma; kmax = 0.01; etav = sigmav*y/x
   uref=1.0; L=0.5; nu = uref*L/1.e6; numax = 1000.*nu
   eps_x = 0.72*kmax**2/numax*sigmav**2*y**2/x**4*exp(-etav**2)
   eps_y = -0.72*kmax**2/numax*sigmav**2/x**2*exp(-etav**2)
   return [eps_x, eps_y]

def forcing_mom(XX):
   '''Forcing function: momentum'''
   x = XX[0]; y = XX[1]
   adv = 1.0; beta = 0.0; mass = 0.0; rho = 1.0
   sigma=4; sigmav = 2.5*sigma; kmax = 0.01; cmu = 0.09
   uref=1.0; L=0.5
   eta = sigma*y/x; etav = sigmav*y/x
   nu = uref*L/1.e6
   numax = 1000.*nu

   u = velocity(XX)[0]
   v = velocity(XX)[1]
   p = pressure(XX)

   p_x = rho*uref**2*(1.-x)/(2.*x-x**2+1./4.)*log(4.*y**3-3.*y**2+5./4.)
   p_y = rho*uref**2*(1.-x)/(4.*y**3-3.*y**2+5./4.)*log(2.*x-x**2+1./4.)
   u_x = grad_u(XX)[0]
   u_y = grad_u(XX)[1]
   v_x = grad_v(XX)[0]
   v_y = grad_v(XX)[1]
   
   u_xx = 4./pi**0.5*eta/x**2*exp(-eta**2)*(1.-eta**2)
   u_yy = -4./pi**0.5*(sigma/x)**2*eta*exp(-eta**2)
   v_xx = 2./pi**0.5*sigma*y**2/x**4*exp(-eta**2)*(3.-2.*eta**2)
   v_yy = 2./pi**0.5*sigma/x**4*exp(-eta**2)*(1.-2.*eta**2)
   
   nut = 1./4.*numax*etav**4*exp(2.-etav**2)
   nut_x = 1./2.*numax*sigmav**4*y**4/x**5*(etav**2-2.)*exp(2.-etav**2)
   nut_y = 1./2.*numax*sigmav**4*y**3/x**4*(2.-etav**2)*exp(2.-etav**2)

   Su = mass*rho*0.0 + adv*rho*(u*u_x+v*u_y) + p_x - (nu+nut)*(u_xx+u_yy) - 2*nut_x*u_x - nut_y*(u_y+v_x)
   Sv = mass*rho*0.0 + adv*rho*(u*v_x+v*v_y) + p_y - (nu+nut)*(v_xx+v_yy) - 2*nut_y*v_y - nut_x*(u_y+v_x)
   return (Su, Sv)

def forcing_k(XX):
   '''Forcing function: k'''
   x = XX[0]; y = XX[1]
   sigma=4; sigmav = 2.5*sigma
   uref=1.0; L=0.5; kmax = 0.01; sigmak = 1.0; cmu = 0.09
   eta = sigma*y/x; etav = sigmav*y/x
   nu = uref*L/1.e6
   numax = 1000.*nu

   u = velocity(XX)[0]
   v = velocity(XX)[1]

   nut = 1./4.*numax*etav**4*exp(2.-etav**2)
   k = tke(XX)
   e = eps(XX)

   k_x = grad_tke(XX)[0]
   k_y = grad_tke(XX)[1]
   k_xx = 2.*kmax*sigmav**2*y**2/x**4*(2.*etav**4-7.*etav**2+3.)*exp(-etav**2)
   k_yy = 2.*kmax*sigmav**2/x**2*(2.*etav**4-5.*etav**2+1.)*exp(-etav**2)
   nut_x = 1./2.*numax*sigmav**4*y**4/x**5*(etav**2-2.)*exp(2.-etav**2)
   nut_y = 1./2.*numax*sigmav**4*y**3/x**4*(2.-etav**2)*exp(2.-etav**2)
   u_x = grad_u(XX)[0]
   u_y = grad_u(XX)[1]
   v_x = grad_v(XX)[0]
   v_y = grad_v(XX)[1]
   S = 2.*(u_x**2+v_y**2)+(u_y+v_x)**2

   Sk = u*k_x + v*k_y - (nu+nut/sigmak)*((k_xx+k_yy) + 1./sigmak*(nut_x*k_x+nut_y*k_y)) - nut*S + e

   return Sk

def forcing_eps(XX):
   '''Forcing function: epsilon'''
   x = XX[0]; y = XX[1]
   sigma=4; sigmav = 2.5*sigma
   uref=1.0; L=0.5; kmax = 0.01
   sigmae = 1.3; ce1 = 1.44; ce2 = 1.92; cmu = 0.09
   eta = sigma*y/x; etav = sigmav*y/x
   nu = uref*L/1.e6
   numax = 1000.*nu

   u = velocity(XX)[0]
   v = velocity(XX)[1]
   nut = 1./4.*numax*etav**4*exp(2.-etav**2)

   # Need to prevent k=0:
   kmin = 1.e-3
   k = max(tke(XX), kmin)
   e = eps(XX)

   eps_x = grad_eps(XX)[0]
   eps_y = grad_eps(XX)[1]
   eps_xx = 0.72*kmax**2/numax*sigmav**2*y**2/x**4*(2.*etav**2-3.)*exp(-etav**2)
   eps_yy = -0.72*kmax**2/numax*sigmav**2/x**2*(2.*etav**2-1.)*exp(-etav**2)

   nut_x = 1./2.*numax*sigmav**4*y**4/x**5*(etav**2-2.)*exp(2.-etav**2)
   nut_y = 1./2.*numax*sigmav**4*y**3/x**4*(2.-etav**2)*exp(2.-etav**2)
   u_x = grad_u(XX)[0]
   u_y = grad_u(XX)[1]
   v_x = grad_v(XX)[0]
   v_y = grad_v(XX)[1]
   S = 2.*(u_x**2+v_y**2)+(u_y+v_x)**2

   Se = u*eps_x + v*eps_y - (nu+nut/sigmae)*((eps_xx+eps_yy) + 1./sigmae*(nut_x*eps_x+nut_y*eps_y)) - ce1*e/k*nut*S + ce2*e**2/k

   return Se

