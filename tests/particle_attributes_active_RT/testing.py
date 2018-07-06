#def val(t):
from numpy import arange,zeros,reshape,concatenate
x=zeros(250*500,float)
y=zeros(250*500,float)
for i in range(0,500):
  for j in range(0,250):
    x[(i*250)+j]=i/500.
    y[(i*250)+j]=0.5+j/500.
    print x[(i*250)+j]
    print y[(i*250)+j]

  
  #x = 1.0*arange(0,500.)/500.
  #y = 0.5 + arange(0,250.)/250.
#  return reshape(concatenate((x,y)),(2,125000)).T
