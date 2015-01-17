from math import exp

Lx = 7000.0
Ly = 7000.0
R = 200.0 # Vent RADIUS
dx = 50 # Characteristic element length

print "BC_x_w(1) = 0.0"
print "BC_x_e(1) = %f" % R
print "BC_y_s(1) = 0.0"
print "BC_y_n(1) = 0.0"
print "BC_TYPE(1) = 'MASS_INFLOW'"
print "BC_P_g(1) = 101325.0"
print "\n"

print "BC_x_w(2) = 0.0"
print "BC_x_e(2) = %f" % Lx
print "BC_y_s(2) = %f" % Ly
print "BC_y_n(2) = %f" % Ly
print "BC_TYPE(2) = 'P_OUTFLOW'"
print "BC_P_g(2) = %f" % (101325.0*exp(-(9.8*28.42*Ly)/(8314.56*288.15)))
print "\n"

print "BC_x_w(3) = %f" % R
print "BC_x_e(3) = %f" % Lx
print "BC_y_s(3) = 0.0"
print "BC_y_n(3) = 0.0"
print "BC_TYPE(3) = 'FREE_SLIP_WALL'"
print "BC_Tw_g(3) = 0.0"
print "BC_Tw_s(3,1) = 0.0"
print "\n"

print "BC_x_w(4) = 0.0"
print "BC_x_e(4) = 0.0"
print "BC_y_s(4) = 0.0"
print "BC_y_n(4) = %f" % Ly
print "BC_TYPE(4) = 'FREE_SLIP_WALL'"
print "BC_Tw_g(4) = 0.0"
print "BC_Tw_s(4,1) = 0.0"
print "\n"

n = 5
for bc in range(0, int(Ly), dx):
   print "BC_x_e(%d) = %f" % (n, Lx)
   print "BC_x_w(%d) = %f" % (n, Lx)
   print "BC_y_n(%d) = %f" % (n, bc+dx)
   print "BC_y_s(%d) = %f" % (n, bc)
   print "BC_TYPE(%d) = 'P_OUTFLOW'" % (n)
   print "BC_P_g(%d) = %f" % (n, 101325.0*exp(-(9.8*28.42*(bc+(dx/2.0)))/(8314.56*288.15)))
   print "\n"
   n = n+1
