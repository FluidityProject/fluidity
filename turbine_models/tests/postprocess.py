import numpy as np
import matplotlib.pyplot as plt
from math import pi
print 'hello from postprocessing\n'

StaticData=np.genfromtxt('FoilData/NACA_0012.dat',skip_header=11,skip_footer=51)

Static_AOA=StaticData[:,0]
Static_CL=StaticData[:,1]

Dynamic_AOA=[];
Dynamic_CL=[];

for i in range (1,100):
    DynamicData=np.genfromtxt('Test_'+str(i)+'.dat',skip_header=1,delimiter=' , ')
    Dynamic_AOA.append(DynamicData[11,1]*180.0/pi)
    Dynamic_CL.append(DynamicData[11,3])

plt.plot(Static_AOA,Static_CL,'-bo',label='Static loading response')
plt.plot(Dynamic_AOA,Dynamic_CL,'-r.',label='Dynamic loading response')
plt.title('Lift response of NACA_0012 under harmonic pitch motion')
plt.xlim((-10,30))
plt.ylim((-1,2.5))
plt.xlabel('angle of attack [degrees]')
plt.ylabel('$C_L$')
plt.legend(loc='upper left')
plt.show()
print 'Goodbye from postprocessing\n'
