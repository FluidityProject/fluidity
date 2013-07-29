#! /usr/bin/env python
import argparse
import sys
import numpy as np
from math import sqrt

def find(x):
    for k, nods in enumerate(halonodes):
        try:
            l = nods.index(x)
            return k, l
        except:
            continue

def node_mapping(xarr):
    ans=[]
    for x in xarr:
        if x<=nnods[i]:
            ans.append(x+nnods[0:i].sum())
        else:
            ii,jj=find(x)
            ans.append(sent_nodes[ii][i][jj]+nnods[0:ii].sum())
    return ans

# get base filename and number of files
parser=argparse.ArgumentParser(description='Creates single .node and .ele files for a cubed sphere mesh generated in parallel.')
parser.add_argument('basename', help='base filename without number or extension')
parser.add_argument('n', metavar='n', type=int, help='number of partitions origianlly generated')
args=parser.parse_args()
n=args.n

# open new node and ele files
nodefile=open(args.basename+'_serial.node','w')
elefile=open(args.basename+'_serial.ele','w')

# array for storing number of owned nodes for each process
nnods=np.zeros((n,),dtype=int)
# list of lists of sent nodes
sent_nodes=[[[] for i in range(n)] for i in range(n)]
# initial nodes and elements from first process
allnodes=np.loadtxt(args.basename+'_0.node',skiprows=1)
alleles=np.loadtxt(args.basename+'_0.ele',skiprows=1)
# get number of owned nodes and sent nodes from first process
h=open(args.basename+'_'+str(0)+'.halo','r')
h.readline()
h.readline()
nnods[0]=int(h.readline().split('\"')[3])
for j in range(n):
    h.readline()
    sent_nodes[0][j]=[int(x) for x in h.readline().split('>')[1].split('<')[0].split()]
    h.readline()
    h.readline()
h.close()
# set of local node numbers owned by this process
nodes=set(range(1,nnods[0]+1))
# find any eles containing nodes from processors of higher rank and remove them
e=[]
for k in range(alleles.shape[0]):
    if not set(alleles[k,1:5]).issubset(nodes):
        e.append(k)
alleles=alleles[np.setdiff1d(np.arange(alleles.shape[0]),e),:]

# loop over remaining files
for i in range(1,n):
    # read in nodes and adjust node number
    total_nnods=nnods.sum()
    allnodes=np.append(allnodes[0:total_nnods,:],np.loadtxt(args.basename+'_'+str(i)+'.node',skiprows=1),axis=0)
    allnodes[total_nnods:,0]+=total_nnods

    # read in next file of elements
    neweles=np.loadtxt(args.basename+'_'+str(i)+'.ele',skiprows=1)

    # read in level 1 halo nodes and number of owned nodes for this process
    halonodes=[[] for j in range(n)]
    h=open(args.basename+'_'+str(i)+'.halo','r')
    h.readline()
    h.readline()
    nnods[i]=int(h.readline().split('\"')[3])
    # set of local node numbers owned by this process
    nodes=set(range(1,nnods[i]+1))
    for j in range(n):
        h.readline()
        sent_nodes[i][j]=[int(x) for x in h.readline().split('>')[1].split('<')[0].split()]
        halonodes[j]=[int(x) for x in h.readline().split('>')[1].split('<')[0].split()]
        h.readline()
        if j<i:
            nodes=nodes.union(set(halonodes[j]))
    # find any eles containing nodes from processors of higher rank and remove them
    e=[]
    for k in range(neweles.shape[0]):
        if not set(neweles[k,1:5]).issubset(nodes):
            e.append(k)
    neweles=neweles[np.setdiff1d(np.arange(neweles.shape[0]),e),:]
    for e in neweles[:,1:5]:
        alleles=np.append(alleles,np.array([0]+node_mapping(e)+[0]).reshape(1,6), axis=0)

# add on nodes from final process
h=open(args.basename+'_'+str(n-1)+'.halo','r')
h.readline()
h.readline()
total_nnods+=int(h.readline().split('\"')[3])
h.close()

neles=alleles.shape[0]

nx=sqrt((total_nnods-2)/6.)
f=open('serial.debug','w')
for e in range(0,neles):
    # check that element is not twisted
    # get 2 vectors from local nodes (1,2), (1,3) and (1,4)
    v1=allnodes[alleles[e,2]-1,1:4]-allnodes[alleles[e,1]-1,1:4]
    v2=allnodes[alleles[e,3]-1,1:4]-allnodes[alleles[e,1]-1,1:4]
    v3=allnodes[alleles[e,4]-1,1:4]-allnodes[alleles[e,1]-1,1:4]
    # vec1.vec3 and vec2.vec3 should be > 0
    try:
        assert((v1*v3).sum()>0)
        assert((v2*v3).sum()>0)
    except:
        f.write('ele: '+str(e+1)+' is twisted\n')
        f.write(str(alleles[e,1])+' '+str(alleles[e,2])+' '+str(alleles[e,3])+' '+str(alleles[e,4])+'\n')
        f.write(str(v1)+' '+str(v2)+' '+str(v3)+'\n')
        f.write(str((v1*v3).sum())+'\n')
        f.write(str((v2*v3).sum())+'\n')
    # check that distances between nodes in element are <=dx
    v3=allnodes[alleles[e,4]-1,1:4]-allnodes[alleles[e,2]-1,1:4]
    v4=allnodes[alleles[e,4]-1,1:4]-allnodes[alleles[e,3]-1,1:4]
    try:
        assert(sqrt((v1*v1).sum())<=2./(nx-1.))
        assert(sqrt((v2*v2).sum())<=2./(nx-1.))
        assert(sqrt((v3*v3).sum())<=2./(nx-1.))
        assert(sqrt((v4*v4).sum())<=2./(nx-1.))
    except:
        f.write('ele: '+str(e+1)+' has too great a distance between its nodes\n')
        f.write(str(alleles[e,1])+' '+str(alleles[e,2])+' '+str(alleles[e,3])+' '+str(alleles[e,4])+'\n')
        f.write(str(v1)+' '+str(v2)+' '+str(v3)+' '+str(v4)+'\n')
        f.write(str(sqrt((v1*v1).sum()))+'\n')
        f.write(str(sqrt((v2*v2).sum()))+'\n')
        f.write(str(sqrt((v3*v3).sum()))+'\n')
        f.write(str(sqrt((v4*v4).sum()))+'\n')
f.close()

# write out nodes
nodefile.write(str(total_nnods)+" 3 0 0\n")
np.savetxt(nodefile,allnodes[0:total_nnods,:],fmt="%d %.10f %.10f %.10f")
nodefile.close()

# write out elements
alleles[:,0]=range(1,alleles.shape[0]+1)
alleles[:,5]=0
elefile.write(str(neles)+" 4 1\n")
np.savetxt(elefile,alleles[0:neles,:],fmt="%d")
elefile.close()

