from math import exp

print "BC_x_w(1) = 0.0"
print "BC_x_e(1) = 200.0"
print "BC_y_s(1) = 0.0"
print "BC_y_n(1) = 0.0"
print "BC_TYPE(1) = 'MASS_INFLOW'"
print "BC_P_g(1) = 101325.0"
print "\n"

print "BC_x_w(2) = 0.0"
print "BC_x_e(2) = 30000.0"
print "BC_y_s(2) = 10000.0"
print "BC_y_n(2) = 10000.0"
print "BC_TYPE(2) = 'P_OUTFLOW'"
print "BC_P_g(2) = 31684.736532047"
print "\n"

print "BC_x_w(3) = 200.0"
print "BC_x_e(3) = 30000.0"
print "BC_y_s(3) = 0.0"
print "BC_y_n(3) = 0.0"
print "BC_TYPE(3) = 'FREE_SLIP_WALL'"
print "BC_Tw_g(3) = 0.0"
print "BC_Tw_s(3,1) = 0.0"
print "\n"

print "BC_x_w(4) = 0.0"
print "BC_x_e(4) = 0.0"
print "BC_y_s(4) = 0.0"
print "BC_y_n(4) = 10000.0"
print "BC_TYPE(4) = 'FREE_SLIP_WALL'"
print "BC_Tw_g(4) = 0.0"
print "BC_Tw_s(4,1) = 0.0"
print "\n"

n = 5
for bc in range(0, 10000, 100):
   print "BC_x_e(%d) = %d" % (n, 30000)
   print "BC_x_w(%d) = %d" % (n, 30000)
   print "BC_y_n(%d) = %d" % (n, bc+100)
   print "BC_y_s(%d) = %d" % (n, bc)
   print "BC_TYPE(%d) = 'P_OUTFLOW'" % (n)
   print "BC_P_g(%d) = %f" % (n, 101325.0*exp(-(9.8*28.42*(bc+50.0))/(8314.56*288.15)))
   print "\n"
   n = n+1
