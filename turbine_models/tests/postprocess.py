import numpy as np
import matplotlib.pyplot as plt
from math import pi
print 'hello from postprocessing\n'

T=2*pi/25.07
StaticData=np.genfromtxt('FoilData/NACA_0012.dat',skip_header=11,skip_footer=51)

Static_AOA=StaticData[:,0]
Static_CL=StaticData[:,1]


Dynamic_AOA=[];
Dynamic_CL=[];
time=[];

for i in range (1,5000):
    DynamicData=np.genfromtxt('Test_'+str(i)+'.dat',skip_header=1,delimiter=', ')
    Dynamic_AOA.append(DynamicData[11,1]*180.0/pi)
    Dynamic_CL.append(DynamicData[11,2])
    time.append(i*0.001/T)

plt.figure(1)
plt.plot(Static_AOA,Static_CL,'-bo',label='Static loading Exp.')
plt.plot(Dynamic_AOA,Dynamic_CL,'-r',label='Dynamic loading ALM-model')
plt.title('Lift response of NACA_0012 under harmonic pitch motion')
#plt.xlim((-10,30))
#plt.ylim((-1,3.5))
plt.xlabel('angle of attack [degrees]')
plt.ylabel('$C_L$')
plt.legend(loc='upper left')

plt.figure(2)
plt.subplot(2,1,1)
plt.plot(time,Dynamic_CL,'-r',label='Dynamic Lift - ALM')
plt.title('Lift timeseries')
#plt.xlim((0,7))
#plt.ylim((0,1.2))
plt.xlabel('Number of Oscillations')
plt.ylabel('$C_L$')
#plt.legend(loc='upper right')
plt.subplot(2,1,2)
plt.plot(time,Dynamic_AOA)
#plt.xlim((0,7))
#plt.ylim((0,10))
plt.xlabel('Number of Oscillations')
plt.ylabel('angle of attack [degrees]')


plt.show()
print 'Goodbye from postprocessing\n'
