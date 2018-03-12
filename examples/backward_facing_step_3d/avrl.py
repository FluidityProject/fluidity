#!/usr/bin/env python
import sys
import numpy

def moving_average(name):
  npy = numpy.load(str(name))
  av = numpy.zeros([len(npy[:,0])])
  av[0]=npy[0,0]
  init=200.
  for i in range(1,len(av)):
    t = npy[i,-1]
    dt = t-npy[i-1,-1]
    # a = (t-dt)/t, b=dt/t
    a = (t-init-dt)/(t-init)
    b = dt/(t-init)
    # av = a*old_av + b*new_val
    # npy[:,1] is the reattachment length of AverageVelocity
    rl = av[i-1]*a + npy[i,1]*b
    av[i] = rl

  return npy,av

#########################################################################

def main():
  print('Only run inside numpy_data directory')
  try:
    name = sys.argv[1]
  except:
    print('argument is name of numpy data file')

  npy, av = moving_average(name)
  print(npy[:,0])
  print(av)
  numpy.save('av_'+str(name), [av,npy[:,-1]])

if __name__ == "__main__":
  sys.exit(main())

