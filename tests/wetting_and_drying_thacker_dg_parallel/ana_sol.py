# This script is used in diamond

# Set this to the d0 value in diamond!

h0=50.0
eta0=2.0
Rb=430620.0
Rmesh=440000.0
g=9.81

# should be copied from the diamond extrude function. Returns the depth as a positive number.
def bathymetry_function(r):
        from math import cos
        # The middle term is an offset term while the last constant is to ensure positive waterlevel at the boundary domain
        return h0*(Rb**2-r**2)/(Rb**2)-h0*(Rb**2-Rmesh**2)/(Rb**2)+0.5
  
# returns the depth of the analytical fs as a positive number but takes into account the shape of the cup
def analytic_solution(r, t, d0):
        eta=analytic_solution_pure(r, t)
        b=bathymetry_function(r)-d0
        # check if the free surface is lower than the cup and if, then set it to the cup heigth
        return min(eta,b)  

# returns the depth of the analyitical fs as a positive number
def analytic_solution_pure(r, t):
        from math import cos
        A=((h0+eta0)**2-h0**2)/((h0+eta0)**2+h0**2)
        omega=(8.0*g*h0/(Rb**2))**0.5
        eta=(1.0-A**2)**0.5/(1.0-A*cos(omega*t))
        eta=eta-1
        eta=eta-r**2/(Rb**2)*((1.0-A**2)/((1.0-A*cos(omega*t))**2)-1.0)
        eta=h0*eta
        # Add cup offset
        offset=bathymetry_function(Rb)
        return -eta+offset



# returns the velocity components of the analytical fs but takes into account the shape of the cup (that is in dry areas, 0 is returned)
# X is the 3 dimensional position
# t is the time 
# d0 is the wetting and drying parameter
def analytic_solution_vel(X, t, d0):
        r=(X[0]**2+X[1]**2)**0.5
#        if analytic_solution_pure(r, t)>=bathymetry_function(r)-d0:
        if analytic_solution_pure(r, t)>=bathymetry_function(r):  # No d0 correction is necessart here,  since we expect flow up to the point where the analyitical free surface meets the bathymetry
                return [0.0, 0.0, 0.0]
        else:        
                return analytic_solution_vel_pure(X, t)

# returns the velocity components of the analytical fs
def analytic_solution_vel_pure(X, t):
        from math import cos, sin
        A=((h0+eta0)**2-h0**2)/((h0+eta0)**2+h0**2)
        omega=(8.0*g*h0/(Rb**2))**0.5
        u=1.0/(2*(1-A*cos(omega*t)))*omega*X[0]*A*sin(omega*t)
        v=1.0/(2*(1-A*cos(omega*t)))*omega*X[1]*A*sin(omega*t)
        w=omega*A*sin(omega*t)/(1-A*cos(omega*t))*(2*h0*(X[0]**2+X[1]**2)/Rb**2-h0-X[2])

        return [u,v,w]
