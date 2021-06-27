import time_integrated_errors_mod
from math import log


####################################################################################
error1=time_integrated_errors_mod.get_time_integrated_error('thacker_A.stat')
error2=time_integrated_errors_mod.get_time_integrated_error('thacker_B.stat')
rate=[]
for i in range(0, len(error1)):
  rate.append(log(abs(error1[i]/error2[i]), 2))
 
print 'A_B_rates Free surface rates: ' rate

error1=time_integrated_errors_mod.get_time_integrated_vel_error('thacker_A.stat')
error2=time_integrated_errors_mod.get_time_integrated_vel_error('thacker_B.stat')
rate=[]
for i in range(0, len(error1)):
  rate.append(log(abs(error1[i]/error2[i]), 2))
  
print 'A_B_rates Velocity rates: ' rate
####################################################################################
error1=time_integrated_errors_mod.get_time_integrated_error('thacker_B.stat')
error2=time_integrated_errors_mod.get_time_integrated_error('thacker_C.stat')
rate=[]
for i in range(0, len(error1)):
  rate.append(log(abs(error1[i]/error2[i]), 2))
 
print 'B_C_rates Free surface rates: ' rate

error1=time_integrated_errors_mod.get_time_integrated_vel_error('thacker_B.stat')
error2=time_integrated_errors_mod.get_time_integrated_vel_error('thacker_C.stat')
rate=[]
for i in range(0, len(error1)):
  rate.append(log(abs(error1[i]/error2[i]), 2))
  
print 'B_C_rates Velocity rates: ' rate

####################################################################################
error1=time_integrated_errors_mod.get_time_integrated_error('thacker_C.stat')
error2=time_integrated_errors_mod.get_time_integrated_error('thacker_D.stat')
rate=[]
for i in range(0, len(error1)):
  rate.append(log(abs(error1[i]/error2[i]), 2))
 
print 'C_D_rates Free surface rates: ' rate

error1=time_integrated_errors_mod.get_time_integrated_vel_error('thacker_C.stat')
error2=time_integrated_errors_mod.get_time_integrated_vel_error('thacker_D.stat')
rate=[]
for i in range(0, len(error1)):
  rate.append(log(abs(error1[i]/error2[i]), 2))
  
print 'C_D_rates Velocity rates: ' rate
####################################################################################

error1=time_integrated_errors_mod.get_time_integrated_error('thacker_D.stat')
error2=time_integrated_errors_mod.get_time_integrated_error('thacker_E.stat')
rate=[]
for i in range(0, len(error1)):
  rate.append(log(abs(error1[i]/error2[i]), 2))
 
print 'D_E_rates Free surface rates: ' rate

error1=time_integrated_errors_mod.get_time_integrated_vel_error('thacker_D.stat')
error2=time_integrated_errors_mod.get_time_integrated_vel_error('thacker_E.stat')
rate=[]
for i in range(0, len(error1)):
  rate.append(log(abs(error1[i]/error2[i]), 2))
  
print 'D_E_rates Velocity rates: ' rate
