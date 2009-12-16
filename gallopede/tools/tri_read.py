import sys
import string
import array
import gallopede_functions_trav
import gallopede_functions
import os

os.system('cp '+sys.argv[1]+".node Positions.node")
os.system('cp '+sys.argv[1]+".ele Positions.ele")

nodeFile=open(sys.argv[1] + ".node")


data = nodeFile.readline().strip()
data = data.split()
nnodes=int(data[0])
nattribs=int(data[2])
nbnd=int(data[3])

nodeXDictionary={}
nodeYDictionary={}
nodeBndDictionary={}
nodePer=[]
for line in nodeFile.readlines():
       if not line:
              continue
       line = line.strip()
       if len(line) == 0:
              continue
       data = string.splitfields(line)
       if (data[0][0] in string.digits):
              nodeXDictionary[int(data[0])] = string.atof(data[1]) 
              nodeYDictionary[int(data[0])] = string.atof(data[2]) 
              if (nbnd==1):
                     nodeBndDictionary[int(data[0])] = int(data[3+nattribs])
                     if(not(int(data[3+nattribs])==0) and not (int(data[3+nattribs])==1)):
                             nodePer.append(int(data[0]))
              else:
                     nodeBndDictionary[int(data[0])] = 0

nodeFile.close()

polyN1Dictionary={}
polyN2Dictionary={}
polyBndDictionary={}
polyPerDictionary={}
polyFile=open(sys.argv[1] + ".poly")
data = polyFile.readline().strip()
npnodes=int(data[0])
for n in range(npnodes):       
       data = string.splitfields(polyFile.readline())
       if (int(data[0])== npnodes) :
              continue
data = string.splitfields(polyFile.readline().strip())
npoly=int(data[0])
nbnd=int(data[1])
for n in range(npoly):
       line=polyFile.readline()
       if not line:
              continue
       line = line.strip()
       if len(line) == 0:
              continue
       data = string.splitfields(line)
       if (data[0][0] in string.digits):
              polyN1Dictionary[int(data[0])] = int(data[1])
              polyN2Dictionary[int(data[0])] = int(data[2]) 
              if (nbnd==1):
                     polyBndDictionary[int(data[0])] = int(data[3])
              else:
                     polyBndDictionary[int(data[0])] = 1
              polyper=array.array('i')
              polyper.fromlist([int(data[1]),int(data[2])])
              if(int(data[3]) in polyPerDictionary):                     
                     polyper.fromlist(polyPerDictionary[int(data[3])].tolist())       
              polyPerDictionary[int(data[3])]= polyper
data = string.splitfields(polyFile.readline().strip())

holeXDictionary={}
holeYDictionary={}
maxbnd=polyPerDictionary.keys()
for x in maxbnd:
       if x<2:
              maxbnd.remove(x)
nholes=int(data[0])

for n in range(nholes):
       line=polyFile.readline()
       if not line:
              continue
       line = line.strip()
       if len(line) == 0:
              continue       
       data = string.splitfields(line)
       holeXDictionary[int(data[0])]=string.atof(data[1])
       holeYDictionary[int(data[0])]=string.atof(data[2])       
polyFile.close()
              
elelookup={}
for n in range(nnodes):
       elelookup[n+1]=n+1
for n in maxbnd:
       for m in polyPerDictionary[n]:
              c=len(polyPerDictionary[n])
              p_i=polyPerDictionary[n].index(m)
              mm=polyPerDictionary[-n][c-p_i-1]
              elelookup[m]=min(elelookup[m],elelookup[mm])
              elelookup[mm]=elelookup[m]


eleFile=open(sys.argv[1] + ".ele")
data=string.splitfields(eleFile.readline())
nele=int(data[0])
nloc=int(data[1])
neattribs=data[2]
eleDictionary={}

for n in range(nele):
       data=string.splitfields(eleFile.readline())
       for nn in (range(nloc)):
              data[nn+1]=elelookup[int(data[nn+1])] 
       eleDictionary[n+1]=data

ghostDictionary={}
for n in range(nele):
       for nn in eleDictionary[n+1][1:nloc+1]:
              if nn not in ghostDictionary:
                        ghostDictionary[nn]=[]
              ghostDictionary[nn].append(n+1)

print ghostDictionary

if (nloc==6):
       for n in maxbnd:
              m=polyPerDictionary[n]
              for nn in range(len(m)/2):
                     mm=m[2*nn:2*nn+2]
                     pair=[]
                     for x in ghostDictionary[elelookup[mm[0]]]:
                            if elelookup[mm[1]] in \
                                  eleDictionary[x][1:nloc+1]:
                                   p1=1+eleDictionary[x][1:nloc+1].index(elelookup[mm[0]])
                                   p2=1+eleDictionary[x][1:nloc+1].index(elelookup[mm[1]])
                                   print x,p1, p2
                                   print gallopede_functions.middle(p1,p2)
                                   pair.append(eleDictionary[x][gallopede_functions.middle(p1,p2)])
                     elelookup[pair[0]]=min(pair)
                     elelookup[pair[1]]=min(pair)
       for n in range(nele):
              data=eleDictionary[n+1]
              for nn in (range(nloc)):
                     data[nn+1]=elelookup[int(data[nn+1])] 
                     eleDictionary[n+1]=data
       ghostDictionary={}
       for n in range(nele):
              for nn in eleDictionary[n+1][1:nloc+1]:
                     if nn not in ghostDictionary:
                            ghostDictionary[nn]=[]
                     ghostDictionary[nn].append(n)


ghostlist=[]
for n in range(nnodes):
       if( not (n+1) in ghostDictionary):
              ghostlist.append(n+1)
nghosts=len(ghostlist)
real_num=range(1,nnodes-nghosts+1)
old_num=ghostDictionary.keys()

X=nodeXDictionary.values()
Y=nodeYDictionary.values()
b=gallopede_functions_trav.bottom(X,Y)
u=gallopede_functions_trav.U(nodeXDictionary.values(),nodeYDictionary.values())
v=gallopede_functions_trav.V(nodeXDictionary.values(),nodeYDictionary.values())
D=gallopede_functions_trav.height(nodeXDictionary.values(),nodeYDictionary.values())


if (sys.argv[2]=="variables"):
       nodeFile=open("Variables.node",'w')
       outline= "%d %d %d %d\n" % (nnodes-nghosts,2,3*len(D[0])+1,1)
       nodeFile.write(outline)
       for n in real_num:
              outline= "%d %f %f" %(n,nodeXDictionary[old_num[n-1]], \
                                         nodeYDictionary[old_num[n-1]])
              for l in range(len(D[0])):
                     outline+= " %8.7e %8.7e %8.7e" % (u[old_num[n-1]-1][l],\
                             v[old_num[n-1]-1][l],D[old_num[n-1]-1][l])
              outline += " %8.7e"%(b[old_num[n-1]-1])
              if (abs(nodeBndDictionary[old_num[n-1]])==1):
                     outline += " %d\n" %(nodeBndDictionary[old_num[n-1]])       
              else:
                     outline += " %d\n" %(0)
              nodeFile.write(outline)       
              
       eleFile=open("Variables.ele",'w')
       outline= "%d %d %d \n" % (nele,nloc,0)
       eleFile.write(outline)

       for n in range(nele):
              outline=eleDictionary[n+1][0]
              for nn in range(nloc):
                     outline +=" %d"%(old_num.\
                                    index(eleDictionary[n+1][nn+1])+1)
              outline += "\n"
              eleFile.write(outline)
       eleFile.close()

if (sys.argv[2]=="clean"):

       nodeFile=open("Connectivity.node",'w')
       outline= "%d %d %d %d\n" % (nnodes-nghosts,2,0,1)
       nodeFile.write(outline)
       for n in real_num:
              outline= "%d %f %f" %(n,nodeXDictionary[old_num[n-1]], \
                                         nodeYDictionary[old_num[n-1]])
              if (abs(nodeBndDictionary[old_num[n-1]])==1):
                     outline += " %d\n" %(nodeBndDictionary[old_num[n-1]])       
              else:
                     outline += " %d\n" %(0)
              nodeFile.write(outline)       

       eleFile=open("Connectivity.ele",'w')
       outline= "%d %d %d \n" % (nele,nloc,0)
       eleFile.write(outline)

       for n in range(nele):
              outline=eleDictionary[n+1][0]
              for nn in range(nloc):
                     outline +=\
                         " %d"%(old_num.index(eleDictionary[n+1][nn+1])+1)
              outline += "\n"
              eleFile.write(outline)
       eleFile.close()
