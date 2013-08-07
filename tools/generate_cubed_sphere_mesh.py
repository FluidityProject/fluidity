#! /usr/bin/env python
import argparse
import sys
import numpy as np
from math import sqrt, pi, tan
from mpi4py import MPI

haloheader='''<?xml version="1.0" encoding="utf-8" ?>
<halos process="<proc>" nprocs="<nprocs>">
'''

halostart='''    <halo level="<level>" n_private_nodes="<n>">
'''

halotemplate='''        <halo_data process="<hp>">
            <send><send_list></send>
            <receive><receive_list></receive>
        </halo_data>
'''

haloend='''    </halo>
'''

halofooter='''</halos>
'''

def get_sizes(panel,xn):

    nxp=(nx/n) # base value
    extra=nx%n # extras to add on later
    nyp=nxp # for now

    # this deals with the case where nx is not exactly divisible by n
    # by putting an extra node in the first 'extra' partitions in the
    # x and y dirns (on panel 0 and corresponding partitions on other panels):
    
    if panel<2:
        if xn[0] < extra:
            nxp+=1
    else:
        if n-xn[0]-1 < extra:
            nxp+=1

    if panel in [0,1,2,5]:
        if xn[1] < extra:
            nyp+=1
    else:
        if n-xn[1]-1 < extra:
            nyp+=1


    #---------------------------------------------------------------------
    # make sure front and back panels own their edges

    # If the partition is on the right or top edge of panel 0 we need to
    # include the right or top edge in this partition...
    if panel==0:
        if xn[0]==n-1:
            nxp+=1
        if xn[1]==n-1:
            nyp+=1

    #...taking it from the bottom of panel 2 or 1:
#    if panel==2 and xn[1]==0:
#        nyp-=1

    # If the partition is at the top or right of panel 3 we need to
    # include the top or right edge in this partition...
    if panel==3:
        if xn[0]==0:
            nxp+=1
        if xn[1]==0:
            nyp+=1

     #...taking it from the left of panel 5 or 4:
    if (panel==5 or panel==4) and xn[0]==0:
        nxp-=1

    #---------------------------------------------------------------------
    # make panels line up 

    # partitions at the top of panel 1 own neither their top nor their
    # bottom edges:
    if panel==1 and xn[1]==n-1:
        nyp-=1

    # partitions on the left of panel 2 own both their left and right
    # edges while those on the right own neither and those on the top
    # do not own their top or bottom edges:
    if panel==2:
        if xn[0]==0:
            nxp+=1
        if xn[0]==n-1:
            nxp-=1
        if xn[1]==n-1:
            nyp-=1

    # partitions at the top of panel 4 own neither of their top or
    # bottom edges, while partitions at the bottom own both:
    if panel==4:
        if xn[1]==0:
            nyp+=1
        if xn[1]==n-1:
            nyp-=1

    return nxp, nyp

def get_start_coords(panel,xn):

    xstart=np.zeros((2))
    # base value of nxp:
    nxp_base=nx/n
    # this deals with the case where nx is not exactly divisible by n
    # by putting an extra node in the first 'extra' partitions in the
    # x and y dirns (on panel 0 and corresponding partitions on other panels):
    extra=nx%n # extras to add on
    if panel<2:
        if xn[0] < extra:
            xstart[0]=(xn[0]*(nxp_base+1))*dx
        else:
            xstart[0]=(extra*(nxp_base+1)+(xn[0]-extra)*nxp_base)*dx
    else:
        if n-xn[0]-1 < extra:
            xstart[0]=(((n-extra)*nxp_base)+(extra-n+xn[0])*(nxp_base+1))*dx
        else:
            xstart[0]=(xn[0]*nxp_base)*dx
        
    if panel in [0,1,2,5]:
        if xn[1] < extra:
            xstart[1]=(xn[1]*(nxp_base+1))*dx
        else:
            xstart[1]=(extra*(nxp_base+1)+(xn[1]-extra)*nxp_base)*dx
    else:
        if n-xn[1]-1 < extra:
            xstart[1]=(((n-extra)*nxp_base)+(extra-n+xn[1])*(nxp_base+1))*dx
        else:
            xstart[1]=(xn[1]*nxp_base)*dx

    #---------------------------------------------------------------------
    # make sure front (0) and back (3) panels own their edges - set an offset
    # for neighbouring panels:

     # offset left of panel 4 and 5:
    if (panel==4 or panel==5) and xn[0]==0:
        xstart[0]+=dx

    #---------------------------------------------------------------------
    # make panels line up 

    # partitions on panel 1 own their left and top edges, so need
    # an offset in the y-direction:
    if panel==1:
        xstart[1]+=dx

    # partitions on panel 2 own their top and right edges, apart
    # from those on the left which also own their left edge (so no
    # offset for these) and those on the right which lose their right
    # edge to panel 4 (this is done by setting nxp accordingly in
    # get_sizes) 

    if panel==2:
        xstart[1]+=dx
        if xn[0]!=0:
            xstart[0]+=dx

    # partitions on panel 3 own their right and top edges, apart from
    # those on the left which also own their left edges (so no offset
    # for these) and those on the bottom which also own their bottom
    # edges (no offset for these either)
    if panel==3:
        if xn[0]!=0:
            xstart[0]+=dx
        if xn[1]!=0:
            xstart[1]+=dx

    # partitions on panel 4 own their left and top edges so need an
    # offset in the y-direction (apart from partitions on the bottom
    # which retain their bottom edge):
    if panel==4 and xn[1]!=0:
        xstart[1]+=dx

    return xstart

def find_neighbour_info(panel,coord,d):

    # number of adjacent panel in 'd' direction
    # d=0 is left, d=1 is under, d=2 is right, d=3 is above
    p=padj[d][panel]
    # local 'coord' coordinate in current panel - this is the
    # coordinate that changes along this edge of the panel
    lc=lx[panel][coord]
    # index of this coordinate in adjacent panel:
    plc=lx[p].index(lc)
    # is this coordinate the same in the adjacent panel?
    same_coord=plc==coord
    # note, if the coordinate is the same, by construction it runs in the
    # same direction, if is not the same, it runs in the other direction:
    adj=p*nsq
    diff=step[coord]
    loc=[i for i in range(2)]
    if same_coord:
        cadj=coord
        loc[coord]=xn[coord]
        if d==0:
            loc[coord-1]=n-1
            adj+=n-1
        elif d==1:
            loc[coord-1]=n-1
            adj+=n*(n-1)
        else:
            loc[coord-1]=0
    else:
        cadj=coord-1
        diff=-step[coord-1]
        loc[coord-1]=n-1-xn[coord]
        if d<=1:
            loc[coord]=n-1
            adj+=nsq-1
        elif d==2:
            loc[coord]=0
            adj+=n-1
        elif d==3:
            loc[coord]=0
            adj+=n*(n-1)

    adj+=diff*xn[coord]

    nxpadj=get_sizes(p,loc)

    return adj, diff, nxpadj[cadj]

comm=MPI.COMM_WORLD
r=comm.Get_rank()
size=comm.Get_size()

parser=argparse.ArgumentParser(description='Creates node, ele and halo files for a cubed sphere mesh.')
parser.add_argument('meshfile', help='Mesh filename, without extension.')
parser.add_argument('n', metavar='n',type=int,help='number of partitions along each side of the cube face')
parser.add_argument('nx', metavar='nx',type=int,help='number of nodes along each side of the cube face')
parser.add_argument('--debug',help='test mesh and write out debugging files',action='store_true')
parser.add_argument('--project',choices=['equidistant','equiangular'],default='None',help='specify projection from cube mesh to sphere')

args=parser.parse_args()
n=args.n
nx=args.nx

# decomp cube into 6*n^2 partitions
nsq=n**2
# each cube face (panel) has nx^2 nodes
nxsq=nx**2
dx=2./float(nx)
if r==0:
    import pylab
    import mpl_toolkits.mplot3d.axes3d as p3d
    fig=pylab.figure()
    ax=p3d.Axes3D(fig)
    allx=np.empty((size*(nx/n+5)**2,3))
else:
    allx=None

# panel numberings:
# 0 is the front, 1 the top.
#    -----
#    | 1 |
#-------------
#| 5 | 0 | 2 |
#-------------
#    | 4 |
#    -----
#    | 3 |
#    -----

# setup mappings to neighbouring panels.
# left, under, right, above
padj=[[5,5,1,1,3,3],[4,0,0,2,2,4],[2,2,4,4,0,0],[1,3,3,5,5,1]]

# starting coordinates for each panel (bottom left corner):
startx=[[-1.,-1.,-1.],[-1.,-1.,1.],[1.,-1.,1.],[1.,1.,1.],[1.,1.,-1.],[-1.,1.,-1.]]
# direction of coordinate change:
lx=[[0,2,1],[0,1,2],[2,1,0],[2,0,1],[1,0,2],[1,2,0]]
dirn=[[1,1],[1,1],[-1,1],[-1,-1],[-1,-1],[-1,1]]

# who owns shared side?
edge_owned=[[-666,1,1,-666,1,1],[0,-666,0,0,-666,1],[0,1,-666,0,0,-666],[-666,1,1,-666,1,1],[0,-666,1,0,-666,0],[0,0,-666,0,1,-666]]

# within each panel, the partition numbering is:
# -----------------------------
# |  n(n-1)  | ... |   n^2-1  |
# |          |     |          |
# -----------------------------
# |     .    |     |    .     |
# |     .    |     |    .     |
# -----------------------------
# |    0     | ... |  (n-1)   |
# |          |     |          |
# -----------------------------

# find which panel we're on:
panel=r/nsq

# find where we are on this panel:
xn=[0,0]
xn[0]=r%n
if panel==0:
    xn[1]=r/n
else:
    xn[1]=(r%(panel*nsq))/n

nxp, nyp=get_sizes(panel,xn)
nxny=nxp*nyp
# number of owned eles #plus number of side halo eles
neles=(nxp-1)*(nyp-1)#+4*(nxp-1)+4*(nyp-1)

xstart=get_start_coords(panel,xn)

# find neighbouring partition numbers and halo node numbers:
# neigh contains the numbers of the neighbouring partitions, numbered
# anti-clockwise starting from the bottom left corner neighbour
# initialise neigh with values that work away from panel edges:
neigh=np.ma.array([r-n-1,r-n,r-n+1,r+1,r+n+1,r+n,r+n-1,r-1])
neigh.mask=np.zeros((8))

# edge_count counts how many panel edges this partition is on. When
# count==2, cn will have the value of the index of neigh that
# corresponds to a non-existent partition:
cn=0
edge_count=0

# First need to find sizes of neighbouring partitions:
step=[1,n]
nxpadj=get_sizes(panel, [xn[0]-1,xn[1]])
nxpl=nxpadj[1]
nxpadj=get_sizes(panel, [xn[0],xn[1]-1])
nxpu=nxpadj[0]
nxpadj=get_sizes(panel, [xn[0]+1,xn[1]])
nxpr=nxpadj[1]
nxpadj=get_sizes(panel, [xn[0],xn[1]+1])
nxpa=nxpadj[0]

# number of bottom left node in neighbour:
bln=np.ma.array([1 for i in range(8)])

# if we're on the left of a panel, adjust neighbours to the left:
if xn[0]==0:
    edge_count+=1
    cn-=1
    adj, diff, nxpl =find_neighbour_info(panel,1,0)
    neigh[0]=adj-diff
    neigh[6:8]=[adj+diff,adj]
    if panel%2==1:
        bln[0]=3
        bln[6:8]=3

# if we're on the right of a panel, adjust neighbours to the right:
if xn[0]==n-1:
    edge_count+=1
    cn+=1
    adj, diff, nxpr=find_neighbour_info(panel,1,2)
    neigh[2:5]=[adj-diff,adj,adj+diff]
    if panel%2==0: bln[2:5]=3

# if we're on the bottom of a panel, adjust neighbours underneath:
if xn[1]==0:
    edge_count+=1
    cn+=1
    adj, diff, nxpu=find_neighbour_info(panel,0,1)
    neigh[0:3]=[adj-diff,adj,adj+diff]
    if panel%2==0: bln[0:3]=2

# if we're on the top of a panel, adjust neighbours above:
if xn[1]==n-1:
    edge_count+=1
    cn+=7
    adj, diff, nxpa=find_neighbour_info(panel,0,3)
    neigh[4:7]=[adj+diff,adj,adj-diff]
    if panel%2==1: bln[4:7]=2

cnmask={}
cnmask[0]=0
cnmask[2]=2
cnmask[6]=6
cnmask[8]=4

# if we're on 2 edges of the panel one of the corner neighbours
# doesn't exist. cn contains the index of that neighbour:
if edge_count==2:
    neigh.mask[cnmask[cn]]=1

# if there is just one partition on this panel, all 4 corner
# neighbours don't exist:
if edge_count==4:
    neigh.mask[0]=1
    neigh.mask[2]=1
    neigh.mask[4]=1
    neigh.mask[6]=1

bln.mask=neigh.mask

mask2=np.row_stack([neigh.mask,neigh.mask])

nrecv=np.ma.zeros((2,8),dtype=np.int)
nrecv[0]=[1,nxpu,1,nxpr,1,nxpa,1,nxpl]
nrecv[1]=[3,nxpu,3,nxpr,3,nxpa,3,nxpl]
nrecv.mask=mask2

# add onto expected number of nodes to receive:
if (panel==1 or panel==2) and xn[0]==0:
    if xn[1]==0:
        nrecv[:,1]+=1
    if xn[1]==n-1:
        nrecv[:,5]+=1
    
if (panel==4 or panel==5) and xn[1]==0:
    if xn[0]==0:
        nrecv[:,7]+=1
    if xn[0]==n-1:
        nrecv[:,3]+=1
    
# work out which nodes we need to send:
# number of halo nodes to send:
n_halo_nodes=np.ma.zeros((2,8),dtype=np.int)
n_halo_nodes[0]=[1,nxp,1,nyp,1,nxp,1,nyp]
n_halo_nodes[1]=[3,nxp,3,nyp,3,nxp,3,nyp]
n_halo_nodes.mask=mask2

# diff[0] walks along the halo (for edge halos)
# diff[1] steps from the first level 1 node to the first level 2 node (for edge halos)
# for corner halos adding both diff[0] and diff[1] gives the corner of the level 2 halo
diff=np.ma.zeros((2,8),dtype=np.int)
diff[0]=[nxp,  1,  -1, nxp,-nxp,  -1,   1,-nxp]
diff[1]=[1,  nxp, nxp,  -1,  -1,-nxp,-nxp,  1]
diff.mask=mask2

first_halo_nodes=np.ma.zeros((2,8),dtype=np.int)
first_halo_nodes[0]=[0,0,nxp-1,nxp-1,nyp*nxp-1,nyp*nxp-1,(nyp-1)*nxp,(nyp-1)*nxp]
first_halo_nodes[1,:]=first_halo_nodes[0,:]+diff[1,:]
first_halo_nodes.mask=mask2

send_nodes=[[],[]]
send_nodes[0]=[[] for i in range(8)]
send_nodes[1]=[[] for i in range(8)]

for i in range(8):
    if not neigh.mask[i]:
        send_nodes[0][i]=range(first_halo_nodes[0,i],first_halo_nodes[0,i]+n_halo_nodes[0,i]*diff[0,i],diff[0,i])

for i in range(8):
    if i%2==1:
        send_nodes[1][i]=range(first_halo_nodes[1,i],first_halo_nodes[1,i]+n_halo_nodes[1,i]*diff[0,i],diff[0,i])
    elif not neigh.mask[i]:
        send_nodes[1][i]=list(np.array([first_halo_nodes[0,i]]*3)+np.array((diff[0,i],diff[0,i]+diff[1,i],diff[1,i])))

loc={}
loc[0]=[7,5,1,3]
loc[3]=[1,7,3,5]
end={}
end[0]=1
end[3]=0
# need this for the case where there is only one partition per panel
indlist=[]

if panel==0 or panel==3:
    if xn[0]==0:
        if xn[1]==0:
            indlist.append(loc[panel][0])
        if xn[1]==n-1:
            indlist.append(loc[panel][1])
    if xn[0]==n-1:
        if xn[1]==0:
            indlist.append(loc[panel][2])
        if xn[1]==n-1:
            indlist.append(loc[panel][3])
    for l in indlist:
        n_halo_nodes[:,l]+=1
        ind=end[panel]*len(send_nodes[0][l])
        send_nodes[0][l].insert(ind,send_nodes[1][l].pop(ind-end[panel]))
        ind=end[panel]*len(send_nodes[1][l])
        send_nodes[1][l].insert(ind,send_nodes[1][l][ind-end[panel]]+diff[1,l])
        ind=end[panel]*len(send_nodes[1][l])
        send_nodes[1][l].insert(ind,send_nodes[0][l][end[panel]*(len(send_nodes[0][l])-1)]+diff[1,l])

# Generate cube coordinates for this partition:
if args.project==None:
    a=1.0
else:
    a=sqrt(3)/3.
startxp=startx[panel]
lx0=startxp[lx[panel][0]]+dirn[panel][0]*xstart[0]
ly0=startxp[lx[panel][1]]+dirn[panel][1]*xstart[1]

x=np.linspace(lx0,lx0+dirn[panel][0]*(nxp-1)*dx,nxp)
y=np.linspace(ly0,ly0+dirn[panel][1]*(nyp-1)*dx,nyp)
xx,yy=np.meshgrid(x,y)
zz=np.zeros((nxny))
zz[:]=startxp[lx[panel][2]]
pos=np.zeros((nxny,3))
pos[:,lx[panel][0]]=xx.reshape(nxny)
pos[:,lx[panel][1]]=yy.reshape(nxny)
pos[:,lx[panel][2]]=zz

# project cube to sphere:
if args.project=='equidistant':
    for i in range(nxny):
        rr=sqrt(pos[i,0]**2+pos[i,1]**2+pos[i,2]**2)
        pos[i,:]=a**2*pos[i,:]/rr
elif args.project=='equiangular':
    pos[:,lx[panel][0]]=np.tan((pi/4.)*pos[:,lx[panel][0]])
    pos[:,lx[panel][1]]=np.tan((pi/4.)*pos[:,lx[panel][1]])
    for i in range(nxny):
        rr=sqrt(pos[i,0]**2+pos[i,1]**2+pos[i,2]**2)
        pos[i,:]=a**2*pos[i,:]/rr

#import pylab
#import mpl_toolkits.mplot3d.axes3d as p3d
#fig=pylab.figure()
#ax=p3d.Axes3D(fig)
#ax.scatter(pos[:,0],pos[:,1],pos[:,2],s=20,marker='x',color='orange')
#pylab.show()

halo_start=np.ma.zeros((2,8),dtype=np.int)
halo_start.mask=np.row_stack([neigh.mask,neigh.mask])

for i in range(8):
    if not neigh.mask[i]:
        print 'process ', r,'is sending ', len(send_nodes[0][i]), 'nodes to process ', neigh[i]
        comm.Send(pos[send_nodes[0][i],:],dest=neigh[i])
d1=[np.empty((nrecv[0][i],3)) for i in range(8) if not neigh.mask[i]]
j=0
for i in range(8):
    if not neigh.mask[i]:
        print 'process ', r, 'is ready to receive ', d1[j].shape[0], 'nodes from process ', neigh[i]
        comm.Recv(d1[j],source=neigh[i])
        j+=1
for i in range(8):
    if not neigh.mask[i]:
#        print 'process ', r,'is sending ', len(send_nodes[1][i]), 'nodes to process ', neigh[i]
        comm.Send(pos[send_nodes[1][i],:],dest=neigh[i])
d2=[np.empty((nrecv[1][i],3)) for i in range(8) if not neigh.mask[i]]
j=0
for i in range(8):
    if not neigh.mask[i]:
#        print 'process ', r, 'is ready to receive ', d2[j].shape[0], 'nodes from process ', neigh[i]
        comm.Recv(d2[j],source=neigh[i])
        j+=1

j=0
for i in range(8):
    if not neigh.mask[i]:
        halo_start[0,i]=pos.shape[0]+1
        pos=np.append(pos,d1[j][::-1],axis=0)
        j+=1
j=0
for i in range(8):
    if not neigh.mask[i]:
        halo_start[1,i]=pos.shape[0]+1
#        if panel==neigh[i]/nsq and i%2==0:
#            pos=np.append(pos,d2[j],axis=0)
#        else:
        pos=np.append(pos,d2[j][::-1],axis=0)
        j+=1

send_list=[[str(),str()] for i in range(size)]
recv_list=[[str(),str()] for i in range(size)]
hp=[i for i in range(size)]
for i in range(8):
    if not neigh.mask[i]:
        send_list[neigh[i]][0]=' '.join(str(n+1) for n in send_nodes[0][i])
        send_list[neigh[i]][1]=' '.join(str(n+1) for n in send_nodes[1][i])
        recv_list[neigh[i]][0]=' '.join(str(n) for n in range(halo_start[0,i]+nrecv[0][i]-1,halo_start[0,i]-1,-1))
#        if panel==neigh[i]/nsq and i%2==0:
#            recv_list[neigh[i]][1]=' '.join(str(n) for n in range(halo_start[1,i],halo_start[1,i]+nrecv[1][i]))
#        else:
        recv_list[neigh[i]][1]=' '.join(str(n) for n in range(halo_start[1,i]+nrecv[1][i]-1,halo_start[1,i]-1,-1))

# write out nodes and coordinates belonging to this partition:
filename=args.meshfile+'_'+str(r)+'.node'
# number of nodes in this partition:
nnods=pos.shape[0]
# create node array with space for node number
nodes=np.zeros((nnods,4))
# insert node numbers into first column:
nodes[:,0]=range(1,nnods+1)
# put coordinates into nodes array:
nodes[:,1:4]=pos
# write array to file:
#np.savetxt(filename,nodes,fmt="%d %.10f %.10f %.10f",header=str(nnods)+" 4 1\n",comments="")
f=open(filename,'w')
f.write(str(nnods)+" 3 0 0\n")
np.savetxt(f,nodes,fmt="%d %.10f %.10f %.10f")
f.close()

# write out ele file for this partition:
filename=args.meshfile+'_'+str(r)+'.ele'
# create ele array:
eles=np.zeros((neles,6))
# the first elements are those that contain only nodes owned by this partition
for i in range(nyp-1):
    i1=i*nxp-i
    i2=(i+1)*nxp-1-i
    eles[i1:i2,1]=range(i*nxp+1,(i+1)*nxp)
    eles[i1:i2,2]=eles[i1:i2,1]+1
    eles[i1:i2,3]=eles[i1:i2,1]+nxp
    eles[i1:i2,4]=eles[i1:i2,1]+nxp+1

# now the elements containing owned nodes and level 1 halo nodes:
inner_corner_nodes=[1,nxp,(nyp-1)*nxp+1,nxp*nyp]
ends=[0,1,3,2]

node_order={}
node_order[1]=[1,2,3,4]
node_order[3]=[2,4,1,3]
node_order[5]=[4,3,2,1]
node_order[7]=[3,1,4,2]
node_mapping=[[3,1,4,2],[2,4,1,3]]
i1=(nxp-1)*(nyp-1)
i2=i1
inner=[[],[]]
outer=[[],[]]
for i in range(1,8,2):
    j=i/2
    pneigh=neigh[i]/nsq
    inner[0]=range(inner_corner_nodes[ends[j]],inner_corner_nodes[ends[(j+1)%4]]+diff[0,i],diff[0,i])
    outer[0]=range(halo_start[0,i],halo_start[0,i]+nrecv[0,i])
    # need a separate copy as we don't always want to change both lists
    inner[1]=list(outer[0])
    outer[1]=range(halo_start[1,i],halo_start[1,i]+nrecv[1,i])  
    # deal with halos to the left and right
    for m in [-1,1]:
        # ind1 is 0 if m is -1 or -1 if m is 1:
        ind1=-(m+1)/2
        # ind2 gives the number of the edge halo to the left/right of i:
        ind2=(i+2*m)%8
        # ind3 gives the number of the corner halo to the left/right of i:
        ind3=(i+m)%8
        n_owned_edges=edge_owned[panel][pneigh]+edge_owned[panel][neigh[ind2]/nsq]
        
        # if corner halo exists:
        if not neigh.mask[ind3]:
            if pneigh==neigh[ind3]/nsq and (i<ind2 or pneigh!=neigh[ind2]/nsq):
                inner[0].insert(-ind1*len(inner[0]),halo_start[0,ind2]+(ind1+1)*(nrecv[0,ind2]-1))
                inner[0].insert(-ind1*len(inner[0]),halo_start[1,ind2]+(ind1+1)*(nrecv[1,ind2]-1))
                outer[0].insert(-ind1*len(outer[0]),halo_start[0,ind3])
                outer[0].insert(-ind1*len(outer[0]),halo_start[1,ind3]-ind1*2)
                inner[1].insert(-ind1*len(inner[1]),halo_start[0,ind3])
                inner[1].insert(-ind1*len(inner[1]),halo_start[1,ind3]-ind1*2)
                outer[1].insert(-ind1*len(outer[1]),halo_start[1,ind3]+2*(ind1+1))
                outer[1].insert(-ind1*len(outer[1]),halo_start[1,ind3]+1)
        # if corner halo doesn't exist and the cube edge between the 2
        # adjacent panels is not owned by the panel 'underneath':
        elif not edge_owned[pneigh][neigh[ind2]/nsq]:
            # if this partition owns both or neither of the 2 cube
            # edges it shares with partitions to the 'left' and
            # 'below':
            if n_owned_edges==2:
                outer[0].insert(-ind1*len(outer[0]),halo_start[0,ind2]+(ind1+1)*(nrecv[0,ind2]-1))
                inner[1].insert(-ind1*len(outer[1]),halo_start[0,ind2]+(ind1+1)*(nrecv[0,ind2]-1))
                outer[1].insert(-ind1*len(outer[1]),halo_start[1,ind2]+(ind1+1)*(nrecv[1,ind2]-1))
            # if this partition only owns the cube edge between it and the
            # panel 'below':
            elif edge_owned[panel][pneigh]:
                inner[0].insert(-ind1*len(inner[0]),halo_start[0,ind2]+1+(ind1+1)*(nrecv[0,ind2]-3))
                outer[0].insert(-ind1*len(outer[0]),halo_start[0,ind2]+(ind1+1)*(nrecv[0,ind2]-1))
                inner[1].insert(-ind1*len(outer[1]),halo_start[0,ind2]+(ind1+1)*(nrecv[0,ind2]-1))
                outer[1].insert(-ind1*len(outer[1]),halo_start[1,ind2]+(ind1+1)*(nrecv[1,ind2]-1))
            elif n_owned_edges==0:
                outer[1].insert(-ind1*len(outer[1]),halo_start[1,ind2]+(ind1+1)*(nrecv[1,ind2]-1))
                inner[1].insert(-ind1*len(inner[1]),halo_start[0,ind2]+(ind1+1)*(nrecv[0,ind2]-1))
        elif n_owned_edges==0:
            inner[0].insert(-ind1*len(inner[0]),halo_start[0,ind2]+(ind1+1)*(nrecv[0,ind2]-1))
        # if this partition only owns one cube edge (and the edge between
        # the 2 adjacent panels is owned by the panel 'underneath'):
        elif n_owned_edges==1:
            inner.append([])
            outer.append([])
            outer[0].pop(ind1)
            outer[-1].insert(-ind1*len(outer[-1]),outer[1].pop(ind1))
            outer[-1].insert(-ind1*len(outer[-1]),outer[1].pop(ind1))
            outer[-1]=outer[-1][::-1]
            inner[-1].insert(-ind1*len(inner[-1]),outer[1][-ind1*(len(outer[1])-1)])
            outer[1].insert(-ind1*len(outer[1]),inner[1][ind1])
            inner[-1].insert(-ind1*len(inner[-1]),inner[1].pop(ind1))
    ele_nodes=[node_order[i]]*3
    if bln[i]>1:
        new_order=map(lambda x:node_mapping[bln[i]-2][x-1], node_order[i])
        ele_nodes[edge_owned[pneigh][panel]:3]=[new_order]*(3-edge_owned[pneigh][panel])
    for l in range(len(inner)):
        n_neweles=len(inner[l])-1
        neles+=n_neweles
        neweles=np.zeros((n_neweles,6))
        neweles[:,ele_nodes[l][0]]=outer[l][0:n_neweles]
        neweles[:,ele_nodes[l][1]]=outer[l][1:n_neweles+1]
        neweles[:,ele_nodes[l][2]]=inner[l][0:n_neweles]
        neweles[:,ele_nodes[l][3]]=inner[l][1:n_neweles+1]
        eles=np.append(eles,neweles,axis=0)

    inner=inner[0:2]
    outer=outer[0:2]

# insert ele numbers into first column:
eles[:,0]=range(1,neles+1)
# insert boundary id of 0 into last column:
eles[:,5]=0
# write eles to file:
f=open(filename,'w')
f.write(str(neles)+" 4 1\n")
np.savetxt(f,eles,fmt="%d")
f.close()

# debugging tests
if args.debug:
    filename='cube_'+str(r)+'.debug'
    f=open(filename,'w')
    for e in range(0,neles):
        # check that element is not twisted
        # get 2 vectors from local nodes (1,2), (1,3) and (1,4)
        v1=pos[eles[e,2]-1]-pos[eles[e,1]-1]
        v2=pos[eles[e,3]-1]-pos[eles[e,1]-1]
        v3=pos[eles[e,4]-1]-pos[eles[e,1]-1]
        # vec1.vec3 and vec2.vec3 should be > 0
        try:
            assert((v1*v3).sum()>0)
            assert((v2*v3).sum()>0)
        except:
            f.write('ele: '+str(e+1)+' is twisted\n')
            f.write(str(eles[e,1])+' '+str(eles[e,2])+' '+str(eles[e,3])+' '+str(eles[e,4])+'\n')
            f.write(str(v1)+' '+str(v2)+' '+str(v3)+'\n')
            f.write(str((v1*v3).sum())+'\n')
            f.write(str((v2*v3).sum())+'\n')
        # check that distances between nodes in element are <=dx
        v3=pos[eles[e,4]-1]-pos[eles[e,2]-1]
        v4=pos[eles[e,4]-1]-pos[eles[e,3]-1]
        try:
            assert(sqrt((v1*v1).sum())<=2./(nx-1.))
            assert(sqrt((v2*v2).sum())<=2./(nx-1.))
            assert(sqrt((v3*v3).sum())<=2./(nx-1.))
            assert(sqrt((v4*v4).sum())<=2./(nx-1.))
        except:
            f.write('ele: '+str(e+1)+' has too great a distance between its nodes\n')
            f.write(str(eles[e,1])+' '+str(eles[e,2])+' '+str(eles[e,3])+' '+str(eles[e,4])+'\n')
            f.write(str(v1)+' '+str(v2)+' '+str(v3)+' '+str(v4)+'\n')
            f.write(str(sqrt((v1*v1).sum()))+'\n')
            f.write(str(sqrt((v2*v2).sum()))+'\n')
            f.write(str(sqrt((v3*v3).sum()))+'\n')
            f.write(str(sqrt((v4*v4).sum()))+'\n')
    f.close()

halofile=args.meshfile+'_'+str(r)+'.halo'
h=open(halofile,'w')
h.write(haloheader.replace('<proc>',str(r)).replace('<nprocs>',str(size)))
h.write(halostart.replace('<level>',str(1)).replace('<n>',str(nxp*nyp)))
for i in range(size):
    h.write(halotemplate.replace('<hp>',str(hp[i])).replace('<send_list>',send_list[i][0]).replace('<receive_list>',recv_list[i][0]))
h.write(haloend)
h.write(halostart.replace('<level>',str(2)).replace('<n>',str(nxp*nyp)))
for i in range(size):
    h.write(halotemplate.replace('<hp>',str(hp[i])).replace('<send_list>',send_list[i][0]+' '+send_list[i][1]).replace('<receive_list>',recv_list[i][0]+' '+recv_list[i][1]))
h.write(haloend)
h.write(halofooter)
h.close()
