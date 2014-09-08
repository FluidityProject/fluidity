# THIS FILE HAS BEEN AUTOMATICALLY GENERATED.

from numpy import sin, cos, pi, exp, sqrt 

py_dict = {
    "domain_extents" : (1.0, 1.2, 0.8, ),

    "finish_time" : 1.0,

    "pressure1_scale" : 1000.0,

    "saturation1_scale" : 0.5,

    "pressure2_scale" : 1000.0,

    "saturation2_scale" : 0.5,

    "gravity_magnitude" : 1.0,

    "viscosity1" : 1.725e-05,

    "viscosity2" : 0.001,

    "density1" : 1.284,

    "density2" : 1000.0,

    "permeability1" : lambda (s1, s2): 1.567346939e-9*(s1 - 0.2)**2,

    "permeability2" : lambda (s1, s2): 1.567346939e-9*(s2 - 0.3)**2,

    "porosity" : 0.4,

    "absolute_permeability" : 1.567346939e-09,

    "gravity_direction_1D" : 1,

    "pressure1_1D" : lambda x, t: (250.0*cos(1.0*pi*t) + 750.0)*cos(1.0*pi*x),

    "saturation1_1D" : lambda x, t: 1 - 0.5*exp(-1.0*x)/(1.0*t + 1),

    "pressure2_1D" : lambda x, t: (250.0*cos(1.0*pi*t) + 750.0)*cos(1.0*pi*x),

    "saturation2_1D" : lambda x, t: 0.5*exp(-1.0*x)/(1.0*t + 1),

    "darcy_velocity1_x_1D" : lambda x, t: -9.0860692115942e-5*(0.8 - 0.5*exp(-1.0*x)/(1.0*t + 1))**2*(-1.0*pi*(250.0*cos(1.0*pi*t) + 750.0)*sin(1.0*pi*x) - 1.284),

    "darcy_velocity1_magnitude_1D" : lambda x, t: 9.0860692115942e-5*((0.8 - 0.5*exp(-1.0*x)/(1.0*t + 1))**4)**(1/2)*Abs(-1.0*pi*(250.0*cos(1.0*pi*t) + 750.0)*sin(1.0*pi*x) - 1.284),

    "source_saturation1_1D" : lambda x, t: 9.0860692115942e-5*pi**2*(0.8 - 0.5*exp(-1.0*x)/(1.0*t + 1))**2*(250.0*cos(1.0*pi*t) + 750.0)*cos(1.0*pi*x) - 9.0860692115942e-5*(0.8 - 0.5*exp(-1.0*x)/(1.0*t + 1))*(-1.0*pi*(250.0*cos(1.0*pi*t) + 750.0)*sin(1.0*pi*x) - 1.284)*exp(-1.0*x)/(1.0*t + 1) + 0.2*exp(-1.0*x)/(1.0*t + 1)**2,

    "darcy_velocity2_x_1D" : lambda x, t: -1.567346939e-6*(-0.3 + 0.5*exp(-1.0*x)/(1.0*t + 1))**2*(-1.0*pi*(250.0*cos(1.0*pi*t) + 750.0)*sin(1.0*pi*x) - 1000.0),

    "darcy_velocity2_magnitude_1D" : lambda x, t: 1.567346939e-6*((-0.3 + 0.5*exp(-1.0*x)/(1.0*t + 1))**4)**(1/2)*Abs(-1.0*pi*(250.0*cos(1.0*pi*t) + 750.0)*sin(1.0*pi*x) - 1000.0),

    "source_saturation2_1D" : lambda x, t: 1.567346939e-6*pi**2*(-0.3 + 0.5*exp(-1.0*x)/(1.0*t + 1))**2*(250.0*cos(1.0*pi*t) + 750.0)*cos(1.0*pi*x) + 1.567346939e-6*(-0.3 + 0.5*exp(-1.0*x)/(1.0*t + 1))*(-1.0*pi*(250.0*cos(1.0*pi*t) + 750.0)*sin(1.0*pi*x) - 1000.0)*exp(-1.0*x)/(1.0*t + 1) - 0.2*exp(-1.0*x)/(1.0*t + 1)**2,

