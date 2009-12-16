import matplotlib
#matplotlib.interactive(True)
#matplotlib.use('WX')
from pylab import *
import scipy
from scipy import linalg
from scipy.optimize import *
from scipy.integrate import *
from scipy.linalg import *

def get_A_rigid(x,d,rho):
    N=len(rho)
    m=scipy.zeros(N)
    y=scipy.zeros(N)
    A=scipy.zeros([N,N])
    y[1:N]=x[0:N-1]

    m[0]=rho[0]*d[0]*d[0]/(3.0*(d[0]-scipy.sum(y[1:N])))

    for j in range(1,N):
        m[j]=rho[j]*d[j]*d[j]/(3.0*(y[j]+d[j]))

    for j in range(N):
        A[j,j]=m[0]+m[j]+3.0*scipy.sum(m[1:j])
        for k in range(0,j):
            A[j,k]=A[k,j]
        for k in range(j+1,N):
            A[j,k]=m[0]+1.5*m[j]+3.0*scipy.sum(m[1:j])

    return A[1:N,1:N]

def get_A(y,d,rho):
    N=len(rho)
    m=scipy.zeros(N)
    A=scipy.zeros([N,N])


    for j in range(N):
        m[j]=rho[j]*d[j]*d[j]/(3.0*(y[j]+d[j]))

    for j in range(N):
        A[j,j]=m[j]+3.0*scipy.sum(m[0:j])
        for k in range(0,j):
            A[j,k]=A[k,j]
        for k in range(j+1,N):
            A[j,k]=1.5*m[j]+3.0*scipy.sum(m[0:j])

    return A

def hamiltonian(y,c,d,rho):

    if y.ndim==1:
        y=reshape(y,[1,y.shape[0]])
    h=scipy.zeros(y.shape[0])
    N=y.shape[1]/2
    Q=scipy.zeros(N)
    D=scipy.zeros(N)
    m=scipy.zeros(N)
    g=10
    
    for i in range(len(h)):

        AI=scipy.linalg.inv(get_A_rigid(y[i,1:N].reshape(2*N),d,rho))
        AI=0.5*(AI+AI.T)
        h[i]=0.5*scipy.dot(scipy.dot(y[i,N+1:2*N].T,AI),y[i,N+1:2*N])\
          +potential(c,y[i,0:N].reshape(N),d,rho)
         
    return h

def potential(c,x,d,rho):
    N=len(rho)
    D=x+d
    g=10.0

    V=-0.5*scipy.sum(rho*D*(1.0-d/D)*(1.0-d/D))
    for j in range(N): 
            V+=0.5*g/(c*c)*(rho[j]*(scipy.sum(x[j:N])*scipy.sum(x[j:N])\
                   -scipy.sum(x[j+1:N])*scipy.sum(x[j+1:N])))   

    return V

def trav_h(y,t,c,d,rho):

    N=len(rho)
    f=scipy.zeros(2*N-2)
    x=scipy.zeros(2*N)
    A=get_A_rigid(y,d,rho)
    AI=scipy.linalg.inv(A)
    AI=0.5*(AI+AI.T)
    g=10
    x[1:N]=y[0:N-1]
    x[N+1:2*N]=y[N-1:2*N-2]

    Aq=scipy.zeros([N,N,N])
    
    for k in range(N):
        Aq[:,:,k]=rho[0]*d[0]*d[0]\
                       /(3.0*(d[0]-scipy.sum(x[1:N]))*(d[0]-scipy.sum(x[1:N])))
        Aq[k,k,k]+=-rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k,k+1:N,k]+=-1.5*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k+1:N,k,k]+=-1.5*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k+1:N,k+1:N,k]+=-3.0*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))

    f[0:N-1]=scipy.dot(AI,x[N+1:2*N])
        

    for k in range(1,N):
        
        f[N+k-2]=0.5*scipy.dot(scipy.dot(scipy.dot(\
                    scipy.dot(x[N+1:2*N].T,AI),Aq[1:N,1:N,k]),AI),x[N+1:2*N])
        f[N+k-2]+=0.5*rho[k]*(1.0-d[k]*d[k]/((d[k]+x[k])*(d[k]+x[k])))\
                         +rho[0]*scipy.sum(x[1:N])/(d[0]-scipy.sum(x[1:N]))\
                         +0.5*rho[0]*scipy.sum(x[1:N])*scipy.sum(x[1:N])\
                         /((d[0]-scipy.sum(x[1:N]))*(d[0]-scipy.sum(x[1:N])))\
                          -g/(c*c)*((rho[k]-rho[0])*scipy.sum(x[k:N])\
                          +scipy.sum((rho[1:k]-rho[0])*x[1:k]))
    

    return f

def J_trav_h(y,t,c,d,rho):
    N=len(rho)
    J=scipy.zeros([2*N-2,2*N-2])
    x=scipy.zeros(2*N)

    A=get_A_rigid(y,d,rho)
    B=scipy.zeros([N,N])
    C=scipy.zeros([N,N])
    D=scipy.zeros([N,N])
    g=10
    x[1:N]=y[0:N-1]
    x[N+1:2*N]=y[N-1:2*N-2]

    AI=scipy.linalg.inv(A)
    AI=0.5*(AI+AI.T)
    Aq=scipy.zeros([N,N,N])
    Aqq=scipy.zeros([N,N,N])
    
    for k in range(N):

        Aq[:,:,k]=rho[0]*d[0]*d[0]\
                       /(3.0*(d[0]-scipy.sum(x[1:N]))*(d[0]-scipy.sum(x[1:N])))
        Aq[k,k,k]+=-rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k,k+1:N,k]+=-1.5*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k+1:N,k,k]+=-1.5*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k+1:N,k+1:N,k]+=-3.0*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        

        Aqq[:,:,k]=2.0*rho[0]*d[0]*d[0]\
                       /(3.0*(d[0]-scipy.sum(x[1:N]))*(d[0]-scipy.sum(x[1:N]))\
                                                           *(d[0]-scipy.sum(x[1:N])))
        Aqq[k,k,k]+=2.0*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k])*(x[k]+d[k]))
        Aqq[k,k+1:N,k]+=3.0*rho[k]*d[k]*d[k]\
            /(3.0*(x[k]+d[k])*(x[k]+d[k])*(x[k]+d[k]))
        Aqq[k+1:N,k,k]+=3.0*rho[k]*d[k]*d[k]\
            /(3.0*(x[k]+d[k])*(x[k]+d[k])*(x[k]+d[k]))
        Aqq[k+1:N,k+1:N,k]+=6.0*rho[k]*d[k]*d[k]\
            /(3.0*(x[k]+d[k])*(x[k]+d[k])*(x[k]+d[k]))

    for k in range(1,N):
        B[1:N,k]=scipy.dot(scipy.dot(\
                    scipy.dot(AI,Aq[1:N,1:N,k]),AI),x[N+1:2*N])
        for j in range(1,N):
            C[j,k]=-0.5*scipy.dot(x[N+1:2*N].T,scipy.dot(AI,\
                               scipy.dot(Aq[1:N,1:N,j],\
                               scipy.dot(AI,scipy.dot(Aq[1:N,1:N,k],\
                               scipy.dot(AI,x[N+1:2*N]))))
                               +scipy.dot(Aq[1:N,1:N,k],\
                               scipy.dot(AI,scipy.dot(Aq[1:N,1:N,j],\
                               scipy.dot(AI,x[N+1:2*N]))))))

        C[k,k]+=0.5*scipy.dot(x[N+1:2*N].T,scipy.dot(AI,\
                               scipy.dot(Aqq[1:N,1:N,k],\
                                             scipy.dot(AI,x[N+1:2*N]))))
        for j in range(1,N):
            if (j<k):
                D[k,j]=-(rho[j]-rho[0])*g/(c*c)+rho[0]/(d[0]-scipy.sum(x[1:N]))\
                    +rho[0]*scipy.sum(x[1:N])/((d[0]-scipy.sum(x[1:N]))*\
                                                   (d[0]-scipy.sum(x[1:N])))\
                         +rho[0]*scipy.sum(x[1:N])*scipy.sum(x[1:N])\
                         /((d[0]-scipy.sum(x[1:N]))*(d[0]-scipy.sum(x[1:N]))\
                               *(d[0]-scipy.sum(x[1:N])))
            else:
                D[k,j]=-(rho[k]-rho[0])*g/(c*c)+rho[0]/(d[0]-scipy.sum(x[1:N]))\
                    +rho[0]*scipy.sum(x[1:N])/((d[0]-scipy.sum(x[1:N]))*\
                                                   (d[0]-scipy.sum(x[1:N])))\
                         +rho[0]*scipy.sum(x[1:N])*scipy.sum(x[1:N])\
                         /((d[0]-scipy.sum(x[1:N]))*(d[0]-scipy.sum(x[1:N]))\
                               *(d[0]-scipy.sum(x[1:N])))
        D[k,k]+=rho[k]*d[k]*d[k]/((d[k]+x[k])*(d[k]+x[k])*(d[k]+x[k]))


#  J=[ dd H/dqdp   -ddH/dqdq]
#    [ dd H/dpdp   -ddH/dpdq]

   
        J[0:N-1,0:N-1]=-B[1:N,1:N].T
        J[0:N-1,N-1:2*N-2]=C[1:N,1:N]+D[1:N,1:N]
        J[N-1:2*N-2,0:N-1]=AI
        J[N-1:2*N-2,N-1:2*N-2]=B[1:N,1:N]

    return J.T

def eigenvector(c,rho,d):
    N=len(rho)
    M=scipy.zeros([2*N-2,2*N-2])

    M=J_trav_h(scipy.zeros(2*N-2),[0.0],c,d,rho)
    (lamb,v)=linalg.eig(M)
    print lamb
    unstable=[]
    rotate=[]
    for i in range(2*N-2):
        if (scipy.real(lamb[i])>1e-8):
            print 'eigenvalue',i,'=',lamb[i]
            print 'eigen vector=',v[:,i]
            unstable.append(i)
            if (abs(scipy.real(lamb[i]))<1e-8):
                print 'eigenvalue',i,'=',lamb[i]
                print 'eigen vector=',v[:,i]
                rotate.append(i)

    vv=scipy.zeros([2*N-2,2*N-2])
    order=lamb[unstable].argsort()
    forder=scipy.real(lamb).argsort()

    if (len(unstable)==0):
        print 'No unstable manifold'
        return vv, lamb,len(unstable)
    if (len(unstable)==1):
        print 'unstable manifold dimension 1'
        qx=moms_to_change(v[:,unstable].reshape([1,2*N-2]),c,d,rho)
        vv[:,0]=-scipy.real(scipy.sign(qx[0,0])*v[:,unstable[0]])
        if N>1:
            vv[:,1:2*N-2]=scipy.real(v[:,forder[0:2*N-3]])
        print vv
        return vv, scipy.sort(lamb), len(unstable)
    if (len(unstable)==2):
        if (scipy.real(lamb[unstable[0]])==scipy.real(lamb[unstable[1]])):
            print 'unstable manifold dimension 2'
            qx=moms_to_change(v[:,unstable[order[0]]].reshape([1,2*N-2]),c,d,rho)
            vv[:,0]=scipy.sign(qx[0,0])*\
                (v[:,unstable[0]])
            qx=moms_to_change(v[:,unstable[order[1]]].reshape([1,2*N-2]),c,d,rho)
            vv[:,1]=scipy.sign(qx[0,0])*\
                (v[:,unstable[1]])
            if N>2:
                print forder
                vv[:,2:2*N-2]=scipy.real(v[:,forder[0:2*N-4]])
        else:
             print 'unstable manifold dimension 2'
             qx=moms_to_change(v[:,unstable[order[0]]].reshape([1,2*N-2]),c,d,rho)
             vv[:,0]=scipy.sign(qx[0,0])*(v[:,unstable[0]])
             qx=moms_to_change(v[:,unstable[order[1]]].reshape([1,2*N-2]),c,d,rho)
             vv[:,1]=scipy.sign(qx[0,0])*(v[:,unstable[1]])
             if N>2:
                  vv[:,2:2*N-2]=scipy.real(v[:,forder[0:2*N-4]])
        return vv, scipy.sort(lamb),len(unstable)
    else:
        vv=scipy.zeros([2*N-2,len(unstable)])
        for i in range(len(unstable)):
            qx=moms_to_change(v[:,unstable[order[i]]].reshape([1,2*N-2]),c,d,rho)
            vv[:,i]=scipy.sign(qx[0,0])\
                *((v[:,unstable[order[i]]]))
            if N>len(lamb):
                  vv[:,len(lamb):2*N-2]=scipy.real(v[:,forder[0:2*N-len(lamb)]])
        return vv, scipy.sort(lamb),len(unstable)     

def p_half_step(pn,pb,qb,dt,c,d,rho):
    N=len(rho)
    x=scipy.zeros(2*N-2)
    x[0:N-1]=qb
    x[N-1:2*N-2]=pn

    fd=trav_h(x,[0.0],c,d,rho)


    return pn-dt/2.0*fd[N-1:2*N-2]-pb

def J_p_half_step(pn,pb,qb,dt,c,d,rho):
    N=len(rho)
    x=scipy.zeros(2*N-2)
    x[0:N-1]=qb
    x[N-1:2*N-2]=pn

    fd=J_trav_h(x,[0.0],c,d,rho)

    return scipy.eye(N-1,N-1)-dt/2.0*fd[N-1:2*N-2,N-1:2*N-2]

def p_step(pb,qb,dt,c,d,rho):
    N=len(rho)
    x=scipy.zeros(2*N-2)
    x[0:N-1]=qb
    x[N-1:2*N-2]=pb

    fd=trav_h(x,[0.0],c,d,rho)


    return pb+dt/2.0*fd[N-1:2*N-2]

def q_step(qn,pb,qb,dt,c,d,rho):
    N=len(rho)
    x=scipy.zeros(2*N-2)
    x[0:N-1]=qb
    x[N-1:2*N-2]=pb
    fd1=trav_h(x,[0.0],c,d,rho)
    x[0:N-1]=qn
    fd2=trav_h(x,[0.0],c,d,rho)

    return qn-dt/2.0*(fd1[0:N-1]+fd2[0:N-1])-qb

def J_q_step(qn,pb,qb,dt,c,d,rho):
    N=len(rho)
    x=scipy.zeros(2*N-2)
    x[0:N-1]=qn
    x[N-1:2*N-2]=pb

    fd=J_trav_h(x,[0.0],c,d,rho)

    return scipy.eye(N-1,N-1)-dt/2.0*fd[0:N-1,0:N-1]

def gen_leap(y,t,dt,c,d,rho):
    N=len(rho)
    p0=scipy.zeros([N-1])
    q0=scipy.zeros([N-1])
    y1=scipy.zeros([2*N])
    yi=scipy.zeros([2*N-2])
    q0[:]=y[1:N]
    p0[:]=y[N+1:2*N]

    yi[0:N-1]=y[1:N]
    yi[N-1:2*N-2]=y[N+1:2*N]

    yf=yi+dt/2.0*trav_h(yi,0.0,c,d,rho)
#    yf=y

#    print 'ph'
#    ph=fsolve(p_half_step,yf[0:N],args=(p0,q0,dt,c,d,rho)\
    ph=fsolve(p_half_step,yf[N-1:2*N-2],args=(p0,q0,dt,c,d,rho),
              fprime=J_p_half_step,
              xtol=1.0e-10)
    yf[N-1:2*N-2]=ph
    yf=yi+dt*trav_h(yf,0.0,c,d,rho)
    

#    print 'q1'
    q1=scipy.array(fsolve(q_step,yf[0:N-1],args=(ph,q0,dt,c,d,rho),\
                              fprime=J_q_step,xtol=1.0e-10))
#                              xtol=1.0e-8))

    p1=scipy.array(p_step(ph,q1,dt,c,d,rho))

    print 'q1=',q1
    print 'p1=',p1
    
    y1[1:N]=q1
    y1[N+1:2*N]=p1

    y1[0]=-scipy.sum(q1)

#    qx=moms_to_change(y1.reshape([1,2*N-2]),c,d,rho)
#    print 'q_x=', qx

    return y1

def to_thickness(y,d):
    N=y.shape[1]
    nx=y.shape[0]
    D=scipy.zeros([nx,N/2])

    for i in range(nx):
        D[i,0]=d[0]-scipy.sum(y[i,1:N/2])
        for j in range(1,N/2):
            D[i,j]=d[j]+y[i,j]
    

    return D

def moms_to_change(y,c,d,rho):

    nx=y.shape[0]
    N=len(rho)
    py=scipy.zeros([nx,N-1])


    for i in  range(nx):
        A=get_A_rigid(y[i,:].reshape(2*N-2),d,rho)
        AI=scipy.linalg.inv(A)
        AI=0.5*(AI+AI.T)
        py[i,:]=scipy.dot(AI,y[i,N-1:2*N-2])

    return py

def potMap(c,rho,d):
    N=len(d)
    fact=scipy.arange(-0.5,1.0,0.1)
    nx=len(fact)
    V=scipy.zeros([nx,nx])
    X=scipy.zeros([nx,nx])
    Y=scipy.zeros([nx,nx])

    for i in range(N):
        for j in  range(i+1,N):
            for k1 in range(nx):
                for k2 in range(nx):
                    D=scipy.zeros(N)
                    D[i]=fact[k1]*d[i]
                    D[j]=fact[k2]*d[j]
                    X[k1,k2]=D[i]+D[j]+scipy.sum(d[i:N])
                    Y[k1,k2]=D[j]+scipy.sum(d[j:N])
                    V[k1,k2]=potential(c,D,d,rho)

#            print V
            figure()
            pcolor(V)
#            pcolor(V,
#                   shading='flat')

    draw()

    return  


def bifurcationDiagram(c,rho,d):
    N=len(rho)
    c_0=0.1
    lamb=scipy.ndarray([5000,2*N-2])
    cm=scipy.ndarray([5000,2*N-2])
    for i in range(5000):
        cc=c_0+i/250.0
        (v,l,NL)=eigenvector(cc,rho,d)
        cm[i,:]=cc
        lamb[i,:]=scipy.sort(l[:])
#lamb=scipy.reshape(lamb,[len(lamb)/(2*N),2*N])
#cm=scipy.reshape(cm,[len(cm)/(2*N),2*N])

    figure()
    for i in range(lamb.shape[1]):
        plot(cm[:,i],scipy.real(lamb[:,i]))

    (v,l,NL)=eigenvector(c,rho,d)

    plot([c],[l[0]],'ro')
    if N>1:
        plot([c],[l[1]],'bo')

    hold(True)
    title(' Real part of eigenvalues of linearized system')
    xlabel('Wavespeed [m/s]')
    ylabel(r'$Re[\lambda_i]$')
    
    draw()

    return


def trav_maker(c,t_0,t_f,dt,pt,epsilon,rho,d):
 
#c=1.495

#d=scipy.array([50.0,950.0])
#rho=scipy.array([1023,1026])
    N=len(d)
#bifurcationDiagram(c,rho,d)
    (v,l,NL)=eigenvector(c,rho,d)

    print d
    print rho
    theta=2.0*3.141592*(0.0)
    print v.shape
    vv=scipy.zeros(2*N-2)
    if (NL==1):
        vv[:]=1.0*v[:,0]
    if (NL>1):
        for i in range(2*N):
            vv[i]=scipy.array([v[i,0]*scipy.sin(theta)+v[i,0]*scipy.cos(theta)])
    print 'theta=', theta

    DD=scipy.zeros(2*N)
    print vv
    DD[1:N]=epsilon*vv[0:N-1]
    DD[N+1:2*N]=epsilon*vv[N-1:2*N-2]
    print DD
    if True:
        t=arange(0.0,t_f+pt,pt)
        tb=arange(t_0,0.0,pt)
        y=scipy.zeros([len(t)+len(tb),2*N])
        y[len(tb),:]=DD
        for i in range(len(tb)):
            if (any(y[i,0:N]+d<0.0)):
                y[i,:]=0
                break
            yf=y[len(tb)-i,:]
            for k in range(int(pt/dt)):
                yf=gen_leap(yf,tb[len(tb)-i-1],-dt,c,d,rho)
                y[len(tb)-i-1,:]=yf
        for i in range(len(tb),len(tb)+len(t)-1):
            if (any(y[i,0:N]+d<0.0)):
                y[i,:]=0
                break
            yf=y[i,:]
            for k in range(int(pt/dt)):
                yf=gen_leap(yf,t[i-len(tb)],dt,c,d,rho)
            y[i+1,:]=yf
        t=scipy.append(tb,t)    

    yy=to_thickness(y,d)

    if (True):
        figure()
        for i in range (N):
            plot(t,scipy.sum(yy[:,i:N],axis=1),'b')
        hold(True)
        title('Approximate rigid lid travelling wave solution, c='+ "%4.3e "%(c))
        ylabel('Interface height [m]')
        xlabel('X:=x-ct [m]')
        for i in range (N):
            text(pt*len(t)/2,0.8*yy[len(t)/2,i]+scipy.sum(yy[len(t)/2,i+1:N]),
                 '%4.0f'%rho[i])
        figure()
        for i in range (N):
            subplot(N,1,i+1)
            plot(t,scipy.sum(yy[:,i:N],axis=1),'b')
            ylabel('Interface height [m]')
        hold(True)
        title('Approximate rigid lid travelling wave solution, c='+ "%4.3e "%(c))
        xlabel('X:=x-ct [m]')
        figure() 
        pyf=scipy.zeros([len(t),2*N-2])
        pyf[:,0:N-1]=y[:,1:N]
        pyf[:,N-1:2*N-2]=y[:,N+1:2*N]
        py=moms_to_change(pyf,c,d,rho)
        for i in range (N-1):
            plot(t,py[:,i])
        xlabel('X:=x-ct [m]')
        ylabel(r'(D_i)_x')
        title('rates of change, c='+ "%4.3e "%(c))
        figure()
        plot(t,scipy.sum(py*py,axis=1))
        xlabel('X:=x-ct [m]')
        ylabel(r'sum of (D_1)_x^2')
        title('total rate of change, c='+ "%4.3e "%(c))
        figure()    
        plot(t,hamiltonian(y,c,d,rho))
        xlabel('X:=x-ct [m]')
        ylabel(r'H(p,q)')
        title('Hamiltonian, c='+ "%4.3e "%(c))

    file=open('ham_d2.dat','w')
    for j in range(N):
        for i in range(t.size):
            out= "%20.19e "%(yy[i,j])
            file.write(out)
        file.write("\n")
    file.close

    print 'File length =', N*t.size
    draw()
