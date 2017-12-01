#!/usr/bin/env python

# Setup values for the Poiseuille flow simulation

# Please note that the values listed below are only read by Fluidity through
# the .flml file and that the rectangular domain itself is defined in the .geo
# file (loacted in the mesh directory). Please ensure that the domain sizes are
# the same in both files.

# The reference density is set to 1 in the .flml file, but any value would do
# as the velocity field of a Poiseuille flow is independent on the density
# (however density is assumed to be constant).

# Dynamic viscosity of the fluid
visc = 1e-3
# Half-height of the rectangular domain
h = 1.
# Driving pressure, the pressure applied at the inflow boundary
drvP = 1e6
# Pressure at the outflow boundary
outP = 5e5
