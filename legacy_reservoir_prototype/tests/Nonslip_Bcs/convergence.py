from fluidity_tools import stat_parser as stat
import pylab as p

def stats(fname):
    error=stat(fname+".stat")["Fluid"]["Error%magnitude"]["l2norm"][:]
    t=stat(fname+".stat")["ElapsedTime"]["value"][:]
    return t,error



h=(0.2,0.1,0.05,0.025)


f=('A','B','C','D')

error=[]

for H,F in zip(h,f):
    t,e= stats('newtonian_chris_'+F)
    
    error.append(e[-1])

h0=[error[0]/(2**i) for i in range(len(h))]

p.figure()
p.xscale('log',basex=2)
p.yscale('log',basey=2)
p.scatter(h,error)
p.plot(h,h0,'k--',lw=2)

