#import matplotlib
#matplotlib.interactive(True)
#matplotlib.use('WX')
import matplotlib.axes3d as p3
from pylab import *
import scipy
from scipy import linalg
from scipy.optimize import *
from scipy.integrate import *
from scipy.linalg import *
from scipy.fftpack import *
import interrupt
import sys
import curses

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

        AI=scipy.linalg.inv(get_A(y[i,:].reshape(2*N),d,rho))
        AI=0.5*(AI+AI.T)
        h[i]=0.5*scipy.dot(scipy.dot(y[i,N:2*N].T,AI),y[i,N:2*N])\
          +potential(c,y[i,0:N].reshape(N),d,rho)
         
    return h


def get_B(rho):

    N=len(rho)
    B=scipy.zeros([N,N])
    for i in range(N):
        B[i,i:]=rho[i]
        B[i+1:,i]=rho[i]
    return B
        

def potential(c,x,d,rho):
    N=len(rho)
    D=x+d
    g=10.0
    B=get_B(rho)

    V=-0.5*scipy.sum(rho*D*(1.0-d/D)*(1.0-d/D))+0.5*g/(c*c)*scipy.dot(x.T,scipy.dot(B,x))
#    for j in range(N): 
#            V+=0.5*g/(c*c)*(rho[j]*(scipy.sum(x[j:N])*scipy.sum(x[j:N])\
#                   -scipy.sum(x[j+1:N])*scipy.sum(x[j+1:N])))   

    return V

def trav_h(x,t,c,d,rho):

    N=len(x)/2
    f=scipy.zeros(2*N)
    A=get_A(x,d,rho)
    B=get_B(rho)
    AI=scipy.linalg.inv(A)
    AI=0.5*(AI+AI.T)
    g=10

    Aq=scipy.zeros([N,N,N])    
    for k in range(N):
        Aq[k,k,k]=-rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k,k+1:N,k]=-1.5*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k+1:N,k,k]=-1.5*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k+1:N,k+1:N,k]=-3.0*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))

    f[0:N]=scipy.dot(AI,x[N:2*N])
        


    for k in range(N):
        
        f[N+k]=0.5*scipy.dot(scipy.dot(scipy.dot(\
                    scipy.dot(x[N:2*N].T,AI),Aq[:,:,k]),AI),x[N:2*N])
    f[N:]+=0.5*rho*(1.0-d*d/((d+x[0:N])*(d+x[0:N])))\
                          -g/(c*c)*scipy.dot(B,x[0:N])
                        
    

    return f

def J_trav_h(x,t,c,d,rho):
    N=len(x)/2
    J=scipy.zeros([2*N,2*N])
    dH=scipy.zeros([2*N,2*N])

    A=get_A(x,d,rho)
    B=scipy.zeros([N,N])
    C=scipy.zeros([N,N])
    D=scipy.zeros([N,N])
    g=10
    

    AI=scipy.linalg.inv(A)
    AI=0.5*(AI+AI.T)
    Aq=scipy.zeros([N,N,N])
    Aqq=scipy.zeros([N,N,N])
    
    for k in range(N):
        Aq[k,k,k]=-rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k,k+1:N,k]=-1.5*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k+1:N,k,k]=-1.5*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k+1:N,k+1:N,k]=-3.0*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k]))
        
        Aq[k,k,k]=2.0*rho[k]*d[k]*d[k]/(3.0*(x[k]+d[k])*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k,k+1:N,k]=3.0*rho[k]*d[k]*d[k]\
            /(3.0*(x[k]+d[k])*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k+1:N,k,k]=3.0*rho[k]*d[k]*d[k]\
            /(3.0*(x[k]+d[k])*(x[k]+d[k])*(x[k]+d[k]))
        Aq[k+1:N,k+1:N,k]=6.0*rho[k]*d[k]*d[k]\
            /(3.0*(x[k]+d[k])*(x[k]+d[k])*(x[k]+d[k]))

    D=-g/(c*c)*get_B(rho)

    for k in range(N):
        B[:,k]=scipy.dot(scipy.dot(\
                    scipy.dot(AI,Aq[:,:,k]),AI),x[N:2*N])
        for j in range(N):
            C[j,k]=-0.5*scipy.dot(x[N:2*N].T,scipy.dot(AI,\
                               scipy.dot(Aq[:,:,j],\
                               scipy.dot(AI,scipy.dot(Aq[:,:,k],\
                               scipy.dot(AI,x[N:2*N]))))
                               +scipy.dot(Aq[:,:,k],\
                               scipy.dot(AI,scipy.dot(Aq[:,:,j],\
                               scipy.dot(AI,x[N:2*N]))))))

        C[k,k]+=0.5*scipy.dot(x[N:2*N].T,scipy.dot(AI,\
                               scipy.dot(Aqq[:,:,k],\
                                             scipy.dot(AI,x[N:2*N]))))
        D[k,k]+=rho[k]*d[k]*d[k]/((d[k]+x[k])*(d[k]+x[k])*(d[k]+x[k]))


#  J=[ dd H/dqdp   -ddH/dqdq]
#    [ dd H/dpdp   -ddH/dpdq]

   
        J[0:N,0:N]=-B.T
        J[0:N,N:2*N]=C+D
        J[N:2*N,0:N]=AI
        J[N:2*N,N:2*N]=B

        dH[0:N,0:N]=-C-D
        dH[0:N,N:2*N]=-B.T
        dH[N:2*N,0:N]=B
        dH[N:2*N,N:2*N]=AI

    return J

def eigenvector(c,rho,d):
    N=len(rho)
    M=scipy.zeros([2*N,2*N])

    M=J_trav_h(scipy.zeros(2*N),[0.0],c,d,rho)
    (lamb,v)=linalg.eig(M.T)
    for i in range(2*N):
        if (abs(scipy.imag(lamb[i])>1e-8)):
            if (scipy.imag(lamb[i])>0.0):
                v[:,i]=scipy.real(v[:,i])
            else:
                v[:,i]=scipy.imag(v[:,i])
        v[:,i]=v[:,i]/scipy.sqrt(scipy.sum(v[:,i]*v[:,i]))
    print lamb
    unstable=[]
    rotate=[]
    for i in range(2*N):
        if (scipy.real(lamb[i])>1e-8):
            print 'eigenvalue',i,'=',lamb[i]
            print 'eigen vector=',v[:,i]
            unstable.append(i)
            if (abs(scipy.real(lamb[i]))<1e-8):
                print 'eigenvalue',i,'=',lamb[i]
                print 'eigen vector=',v[:,i]
                rotate.append(i)

    vv=scipy.zeros([2*N,2*N])
    order=lamb[unstable].argsort()
    forder=scipy.real(lamb).argsort()
    

    if (len(unstable)==0):
        print 'No unstable manifold'
        return vv, lamb,len(unstable)
    if (len(unstable)==1):
        print 'unstable manifold dimension 1'
        qx=moms_to_change(v[:,unstable].reshape([1,2*N]),c,d,rho)
        vv[:,0]=-scipy.real(scipy.sign(qx[0,0])*v[:,unstable[0]])
        if N>1:
            vv[:,1:2*N]=scipy.real(v[:,forder[0:2*N-1]])
        print vv
        return vv, scipy.sort(lamb), len(unstable)
    if (len(unstable)==2):
        if (scipy.real(lamb[unstable[0]])==scipy.real(lamb[unstable[1]])):
            print 'unstable manifold dimension 2'
            qx=moms_to_change(v[:,unstable[order[0]]].reshape([1,2*N]),c,d,rho)
            vv[:,0]=scipy.sign(qx[0,0])*\
                (v[:,unstable[0]])
            qx=moms_to_change(v[:,unstable[order[1]]].reshape([1,2*N]),c,d,rho)
            vv[:,1]=scipy.sign(qx[0])*\
                (v[:,unstable[1]])
            if N>2:
                print forder
                vv[:,2:2*N]=scipy.real(v[:,forder[0:2*N-2]])
        else:
             print 'unstable manifold dimension 2'
             qx=moms_to_change(v[:,unstable[order[0]]].reshape([1,2*N]),c,d,rho)
             vv[:,0]=scipy.sign(qx[0,0])*(v[:,unstable[0]])
             qx=moms_to_change(v[:,unstable[order[1]]].reshape([1,2*N]),c,d,rho)
             vv[:,1]=scipy.sign(qx[0,0])*(v[:,unstable[1]])
             if N>2:
                  vv[:,2:2*N]=scipy.real(v[:,forder[0:2*N-2]])
        return vv, scipy.sort(lamb),len(unstable)
    else:
        vv=scipy.zeros([2*N,len(unstable)])
        for i in range(len(unstable)):
            qx=moms_to_change(v[:,unstable[order[i]]].reshape([1,2*N]),c,d,rho)
            vv[:,i]=scipy.sign(qx[0,0])\
                *((v[:,unstable[order[i]]]))
            if N>len(lamb):
                  vv[:,len(lamb):2*N]=scipy.real(v[:,forder[0:2*N-len(lamb)]])
        return vv, scipy.sort(lamb),len(unstable)     

def p_half_step(pn,pb,qb,dt,c,d,rho):
    N=len(pb)
    x=scipy.zeros(2*N)
    x[0:N]=qb
    x[N:2*N]=pn

    fd=trav_h(x,[0.0],c,d,rho)


    return pn-dt/2.0*fd[N:2*N]-pb

def J_p_half_step(pn,pb,qb,dt,c,d,rho):
    N=len(rho)
    x=scipy.zeros(2*N)
    x[0:N]=qb
    x[N:2*N]=pn

    fd=J_trav_h(x,[0.0],c,d,rho)

    return scipy.eye(N,N)-dt/2.0*fd[N:2*N,N:2*N]

def p_step(pb,qb,dt,c,d,rho):
    N=len(rho)
    x=scipy.zeros(2*N)
    x[0:N]=qb
    x[N:2*N]=pb

    fd=trav_h(x,[0.0],c,d,rho)


    return pb+dt/2.0*fd[N:2*N]

def q_step(qn,pb,qb,dt,c,d,rho):
    N=len(rho)
    x=scipy.zeros(2*N)
    x[0:N]=qb
    x[N:2*N]=pb
    fd1=trav_h(x,[0.0],c,d,rho)
    x[0:N]=qn
    fd2=trav_h(x,[0.0],c,d,rho)

    return qn-dt/2.0*(fd1[0:N]+fd2[0:N])-qb

def J_q_step(qn,pb,qb,dt,c,d,rho):
    N=len(rho)
    x=scipy.zeros(2*N)
    x[0:N]=qn
    x[N:2*N]=pb

    fd=J_trav_h(x,[0.0],c,d,rho)

    return scipy.eye(N,N)-dt/2.0*fd[0:N,0:N]

def gen_leap(y,t,dt,c,d,rho):
    N=len(rho)
    p0=scipy.zeros([N])
    q0=scipy.zeros([N])
    y1=scipy.zeros([2*N])
    q0[:]=y[0:N]
    p0[:]=y[N:2*N]
    

    if (t<0.0) :
#        tau=1000.0/((1000.0+t)*(1000.0+t))
        tau=(100.0/(t+100.0))/10.0
    else:
        tau=1.0

    yf=y+dt*tau/2.0*trav_h(y,0.0,c,d,rho)
#    yf=y

#    print 'ph'
#    ph=fsolve(p_half_step,yf[0:N],args=(p0,q0,dt,c,d,rho)\
    ph=fsolve(p_half_step,scipy.array(yf[N:2*N]),args=(p0,q0,dt*tau,c,d,rho),
              fprime=J_p_half_step,
              xtol=1.0e-8)
    yf[N:2*N]=ph
    yf=y+dt*tau*trav_h(yf,0.0,c,d,rho)
    

#    print 'q1'
    q1=scipy.array(fsolve(q_step,scipy.array(yf[0:N]),
                          args=(ph,q0,dt*tau,c,d,rho),
                          fprime=J_q_step,
                          xtol=1.0e-8))
#                              xtol=1.0e-8))

    p1=scipy.array(p_step(ph,scipy.array(q1),dt*tau,c,d,rho))

#    print 'q1=',q1
#    print 'p1=',p1
    
    if (N>1):
        y1[0:N]=q1[0:N]
        y1[N:2*N]=p1[0:N]
    else:
        y1[0]=q1
        y1[1]=p1

#    qx=moms_to_change(y.reshape([1,2*N]),c,d,rho)
#    print 'q_x=', qx

    return y1

def to_thickness(y,d):
    N=y.shape[1]
    nx=y.shape[0]
    D=scipy.zeros([nx,N/2])

    for i in range(nx):
        for j in range(N/2):
            D[i,j]=d[j]+y[i,j]
    

    return D

def moms_to_change(y,c,d,rho):

    nx=y.shape[0]
    N=y.shape[1]/2
    h=scipy.zeros(N)
    m=scipy.zeros(N)
    mq=scipy.zeros([N,N])
    f=scipy.zeros(2*N)
    py=scipy.zeros([nx,N])


    for i in  range(nx):
        A=get_A(y[i,:].reshape(2*N),d,rho)
        AI=scipy.linalg.inv(A)
        AI=0.5*(AI+AI.T)
   

        py[i,:]=scipy.dot(AI,y[i,N:2*N])

    return py

def potMap(c,rho,d,ev=True,db=[]):
    N=len(d)
    fact=scipy.arange(0.0,1.01,0.01)
    nx=len(fact)
    if len(db)==0:
        print 'db=[]'
        db=scipy.zeros([2,N])
        db[:,0]=[-0.005*scipy.sum(d),0.005*scipy.sum(d)]
        for i in range(1,N):
            db[:,i]=[-0.4*scipy.sum(d[i:N]),0.3*scipy.sum(d[i:N])]
    V=scipy.zeros([nx,nx])
    X=scipy.zeros([nx,nx])
    Y=scipy.zeros([nx,nx])
    df=scipy.zeros(N)
#    df[0]=0.005
#    df[1:N]=0.1
    (v,l,NL)=eigenvector(c,rho,d)
    v=v/scipy.sqrt(scipy.sum(v[0:N,0]*v[0:N,0]))
    fig=[]

    if (N==3):
        V3=scipy.zeros([nx,nx])
        X3=scipy.zeros([nx,nx])
        Y3=scipy.zeros([nx,nx])
        fig3d=figure()
        ax=p3.Axes3D(fig3d)
        for k3 in range(nx):
            for k1 in range(nx):
                for k2 in range(nx):
                    D=scipy.zeros(N)
                    D[0]=0.3*(db[0,0]+(db[1,0]-db[0,0])*fact[k1])
                    D[1]=0.3*(db[0,1]+(db[1,1]-db[0,1])*fact[k2])
                    D[2]=1.0*(db[0,2]+(db[1,2]-db[0,2])*fact[k3])
                    Y3[k1,k2]=D[0]+scipy.sum(d[0:N])
                    X3[k1,k2]=D[1]+scipy.sum(d[1:N])
                    D[0]=D[0]-D[1]
                    D[1]=D[1]-D[2]
                    V3[k1,k2]=potential(c,D,d,rho)


            ax.contour3D(X3,Y3,V3+D[2]+d[2],[D[2]+d[2],D[2]+d[2]])
        ax.set_zlim(0.5*db[0,2]+d[2],0.5*db[1,2]+d[2])

               

    for i in range(N):
        for j in  range(i+1,N):
            for k1 in range(nx):
                for k2 in range(nx):
                    D=scipy.zeros(N)
                    D[i]=db[0,i]+(db[1,i]-db[0,i])*fact[k1]
                    D[j]=db[0,j]+(db[1,j]-db[0,j])*fact[k2]
                    Y[k1,k2]=D[i]+scipy.sum(d[i:N])
                    X[k1,k2]=D[j]+scipy.sum(d[j:N])
                    D[i]=D[i]-D[j]
                    V[k1,k2]=potential(c,D,d,rho)

#            print V
            figure()
            contour(X,Y,V,[0,0],colors='k',linewidth=1.0) 
            CP=contour(X,Y,V,colors='k') 
            clabel(CP, colors = 'k', fmt = '%2.1f', fontsize=14)

            xlabel('h_%i'%(j+1))
            ylabel('h_%i'%(i+1))
            fig.append(gca())
            hold(True)
            if (ev):
                if (NL>0):
                    plot([scipy.sum(d[j:N]-1.e0*v[j:N,0])],
                         [ scipy.sum(d[i:N]-1.e0*v[i:N,0])],'ro')
#            pcolor(V,
#                   shading='flat')

    if (N==3):
        fig.append(ax)
    draw()

    return fig


def bifurcationDiagram(c,rho,d):
    N=len(rho)
    c_0=0.1
    lamb=scipy.ndarray([5000,2*N])
    cm=scipy.ndarray([5000,2*N])
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
    ylabel(r'Re[\lambda_i]')
    
    draw()

    return

def get_uwp(y,c,d,x,dx):
    N=shape(y)[1]
    nx=shape(y)[0]
    nz=10
    z=scipy.zeros(nz+1)
    z[0:nz]=scipy.arange(0.0,1.0,1.0/nz)
    z[nz]=1.0

    u=scipy.zeros([nx,N])
    w=scipy.zeros([nx,N,nz+1])
    p=scipy.zeros([nx,N,nz+1])
    X=scipy.zeros([nx,N,nz+1])
    Z=scipy.zeros([nx,N,nz+1])

    for i in range(nx):
        X[i,:,:]=x[i]
        for j in range(N):
            u[i,j]=c*(1.0-d[j]/y[i,j])
            for k in range(nz+1):
                p[i,j,k]=-(y[i,j]*z[k]+scipy.sum(y[i,j+1:N]))*u[i,j]
                Z[i,j,k]=y[i,j]*z[k]+scipy.sum(y[i,j+1:N])
            for l in range(j,N):
                p[i,j,:]+=y[i,l]*(u[i,l]-u[i,j])
    for i in range(1,nx-1):
        for j in range(N):
            for k in range(nz+1):
                w[i,j,k]=-(y[i,j]*z[k]\
                               +scipy.sum(y[i,j+1:N]))*(u[i+1,j]-u[i-1,j])\
                               /(2.0*dx)
            for l in range(j,N):
                w[i,j,:]+=(y[i+1,l]*(u[i+1,l]-u[i+1,j])\
                               -y[i+1,l]*(u[i+1,l]-u[i+1,j]))/(2.0*dx)

    return u, w, p, X, Z

def shooting(c,t_0,t_f,dt,pt,epsilon,rho,d,th=0.0):
    if (t_0<0.0):
        th_out=fmin(shooting_wrapper,th,args=(c,t_0,t_f,dt,pt,epsilon,rho,d))
    print 'fmin=',th_out
    return th_out

def shooting_wrapper(theta,c,t_0,t_f,dt,pt,epsilon,rho,d):
#     yy,y,t=trav_maker(c,t_f,0.0,dt,pt,epsilon,rho,d,th=theta,pic=False)
     yy,y,t=trav_maker(c,0.0,t_f,dt,pt,epsilon,rho,d,th=theta,pic=False)
     nx=yy.shape[0]
     N=len(rho)
#     return scipy.sum(scipy.sum(y[0:nx/10,:]*y[0:nx/10,:],axis=1))
#    AI=scipy.linalg.inv(get_A(y[0,:].reshape(2*N),d,rho))
     py=moms_to_change(y,c,d,rho)
#     return -min((y[:,0:N]*py[:,:]).reshape(N*nx))
#     return (scipy.sum(y[0,:]*y[0,:]))
     ndx=(scipy.sign(py[:,0])!=scipy.sign(py[0,0])).argmax()
     print ndx
     rf=min(sum(y[ndx:,0:N]*y[ndx:,0:N]+py[ndx:,:]*py[ndx:,:],axis=1))
     print rf
     return scipy.sqrt(rf)

def progBarOut(screen,progstr,theta,t,rf):
    screen.box()
    screen.addstr(1,30,"trying theta=%f"%theta[0])
    screen.addstr(2,30,"integration progress bar")
    st='%s'%progstr
    screen.addstr(3,1,st)
    screen.addstr(4,1,'t=%f, rf=%f '%(t,rf))
    screen.refresh()

def trav_maker(c,t_0,t_f,dt,pt,epsilon,rho,d,th=scipy.zeros(19),pic=True):
    stdscr = curses.initscr()
    curses.noecho()
    curses.cbreak()
    s=curses.newwin(8,77)
    while True:
        try:
            N=len(d)
            (v,l,NL)=eigenvector(c,rho,d)
            theta=2.0*3.141592*th
            vv=scipy.zeros(2*N)

            if (NL>0):
                vv[:]=cos(theta[0])*v[:,0]
                for i in range(2*N-1):
                    vv[:]+=sin(theta[i])*v[:,i+1]
                vv=vv/scipy.sqrt(scipy.sum(vv*vv))

            DD=scipy.zeros(2*N)
            for i in range(2*N):
                DD[i]=-epsilon*vv[i]
            if True:
                s.clear()
                t=arange(0.0,t_f+pt,pt)
                prog=interrupt.progressBar(0,100,75)
                tb=arange(t_0,0.0,pt)
                y=scipy.zeros([len(t)+len(tb),2*N])
                y[len(tb),:]=DD
                for i in range(len(tb)):
#                    print 't=%f'%tb[len(tb)-i-1]
                    prog.updateAmount(100*i/len(t))
                    if prog.new:  progBarOut(s,prog,th,0.0,0.0)
                    if (any(y[i,0:N]+d<0.0)):
                        y[i,:]=0.0
                        break
                    yf=y[len(tb)-i,:]
                    for k in range(int(pt/dt)):
                        yf=gen_leap(yf,tb[len(tb)-i-1],-dt,c,d,rho)
                    y[len(tb)-i-1,:]=yf
                for i in range(len(tb),len(tb)+len(t)-1):
#                    print 't=%f'%t[i-len(tb)]
#                    print 'M=%f'%scipy.sum(rho*(y[i,0:N]+d))
                    prog.updateAmount(100*(i+len(tb)-1)/len(t))
                    rf=(y[i,N-1]*y[i,N-1]+y[i,N]*y[i,N])
                    if prog.new: progBarOut(s,prog,th,
                                            t[i-len(tb)],scipy.sqrt(rf))
                    if (any(y[i,0:N]+d<0.0)):
                        y[i,:]=0
                        break
                    yf=y[i,:]
                    for k in range(int(pt/dt)):
                        yf=gen_leap(yf,t[i-len(tb)],dt,c,d,rho)
                    y[i+1,:]=yf
            t=scipy.append(tb,t)    

            yy=to_thickness(y,d)
            if pic:
                figure()
                print scipy.sum(yy[:,0:N],axis=1)-scipy.sum(d[0:N])
                for i in range (N):
                    subplot(N,1,i+1)
                    plot(t,scipy.sum(yy[:,i:N],axis=1)-scipy.sum(d[0:N]),'b')
    #            plot(t,[-scipy.sum(d[0:i])]*scipy.ones(len(t)),'k:')
                    ylabel('Interface height [m]')


                file=open('ham_d2.dat','w')
                for j in range(N):
                    for i in range(t.size):
                        out= "%20.19e "%(yy[i,j])
                        file.write(out)
                    file.write("\n")
                file.close

                print 'File length =', N*t.size
                draw()
            
        except interrupt.MyError:
            print "Hello everybody"
            return scipy.zeros(0),scipy.zeros(0),scipy.zeros(0)
        finally:
            curses.echo()
            curses.nocbreak()
            curses.endwin()
            print d
            print rho
            print th
        return yy,y,t

def pic_maker(yy,y,t,c,d,rho):
    N=len(d)
    if (N==0):
        return
    if len(t)<2:
        return
    pt=t[1]-t[0]
    if (True):
        if (len(t)>1001):
            figure()
            fft_out=scipy.fftpack.fft(yy[1000:len(t),N-1])
            plot(fft_out[1:])
        figure()
        for i in range (N):
            plot(t,scipy.sum(yy[:,i:N],axis=1)-scipy.sum(d[0:N]),'k')
        hold(True)
#        title('Approximate travelling wave solution, c='+ "%4.3e "%(c))
        ylabel('Interface height [m]')
        xlabel('X:=x-ct [m]')
        ylim(ymin=-scipy.sum(d))
        for i in range (N):
            text(pt*len(t)/2,0.8*yy[len(t)/2,i]+scipy.sum(yy[len(t)/2,i+1:N]),
                 '%4.0f'%rho[i])
        figure()
        for i in range (N):
            subplot(N,1,i+1)
            plot(t,scipy.sum(yy[:,i:N],axis=1)-scipy.sum(d[0:N]),'k')
            hold(True)
            plot(t,[-scipy.sum(d[0:i])]*scipy.ones(len(t)),'k:')
            ylabel('Interface height [m]')
        hold(True)
#        title('Approximate travelling wave solution, c='+ "%4.3e "%(c))
        xlabel('X:=x-ct[m]')
        figure()
        for i in range (N):
            subplot(N,1,i+1)
            plot(t,yy[:,i],'k')
            hold(True)
            plot(t,d[i]*scipy.ones(len(t)),'k:')
            ylabel('Interface height [m]')
        hold(True)
#        title('Approximate travelling wave solution, c='+ "%4.3e "%(c))
        xlabel('X:=x-ct[m]')
        figure() 
        py=moms_to_change(y,c,d,rho)
        for i in range (N):
            plot(t,py[:,i])
        xlabel('X:=x-ct [m]')
        ylabel(r'(D_i)_x')
        title('rates of change, c='+ "%4.3e "%(c))
        figure()
        plot(t,scipy.sum(py*py,axis=1))
        xlabel('X:=x-ct [m]')
        ylabel(r'sum of (D_1)_x^2')
        title('total rate of change, c='+ "%4.3e "%(c))
        for i in range(N):
            for j in range(N):
                figure()
                plot(yy[:,i],py[:,j])
                plot([yy[0,i]],[py[0,j]],'ko',
                     [yy[len(t)-1,i]],[py[len(t)-1,j]],'kx')
                xlabel('D_%i'%(i+1))
                ylabel(r'dD_%i'%(j+1))
                title('Projection plot, c=%4.3e'%(c))
                if (i!=j):
                    figure()
                    plot(yy[:,i],yy[:,j])
                    xlabel('D_%i'%(i+1))
                    ylabel(r'D_%i'%(j+1))
                    title('Projection plot, c=%4.3e'%(c))
                    figure()
                    plot(py[:,i],py[:,j])
                    xlabel('dD_%i'%(i+1))
                    ylabel(r'dD_%i'%(j+1))
                    title('Projection plot, c=%4.3e'%(c))
        figure()
        plot(t,scipy.sum(y*y,axis=1))
        xlabel('X:=x-ct [m]')
        ylabel(r'sum of (D_1)_x^2')
        title('total solution norm, c='+ "%4.3e "%(c))
        figure()    
        plot(t,hamiltonian(y,c,d,rho))
        xlabel('X:=x-ct [m]')
        ylabel(r'H(p,q)')
        title('Hamiltonian, c='+ "%4.3e "%(c))

        (u,w,p,X,Z)=get_uwp(yy,c,d,t,pt)
        figure()
        hold(True)
        for i in range (N):
            contour(X[:,i,:],Z[:,i,:],p[:,i,:])
        for i in range (N):
            plot(t,scipy.sum(yy[:,i:N],axis=1),'k')
        db=scipy.zeros([2,N])
        db[0,0]=-35.0*abs(scipy.sum(y[:,0:N],axis=1)).max()
        db[1,0]=35.0*abs(scipy.sum(y[:,0:N],axis=1)).max()
        for i in range(1,N-1):
            db[0,i]=2.0*scipy.sum(y[:,i:N],axis=1).min()\
                -1.0*scipy.sum(y[:,i:N],axis=1).max()
            db[1,i]=2.0*scipy.sum(y[:,i:N],axis=1).max()\
                -1.0*scipy.sum(y[:,i:N],axis=1).min()
        db[0,N-1]=2.0*y[:,N-1].min()\
            -1.0*y[:,N-1].max()
        db[1,N-1]=2.0*y[:,N-1].max()\
            -1.0*y[:,N-1].min()

        print db
#        fig=potMap(c,rho,d,ev=False,db=db)
        fig=potMap(c,rho,d,ev=False)
        print fig
        for i in range (N):
            for j in range(i+1,N):
                fig[(i)*(i+1)/2-1+j].plot(scipy.sum(yy[:,j:N],axis=1),
                                          scipy.sum(yy[:,i:N],axis=1),'k:')
        if (N==3):
            fig[len(fig)-1].plot3D(scipy.sum(yy[:,1:N],axis=1),
                                   scipy.sum(yy[:,0:N],axis=1),
                                   scipy.sum(yy[:,2:N],axis=1))
            fig[len(fig)-1].set_zlim(d[2]-abs(y[:,2]).max(),
                                     d[2]+abs(y[:,2]).max())
    draw()

    return
