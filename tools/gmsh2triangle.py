#!/usr/bin/env python3

from optparse import OptionParser
import re
import sys
import os.path

#####################################################################
# Script starts here.
optparser=OptionParser(usage='usage: %prog [options] <filename>',
                       add_help_option=True,
                       description="""This takes a Gmsh 2.0 .msh ascii file """ + 
                       """and produces .node, .ele and .edge or .face files.""")

optparser.add_option("--2D", "--2d", "-2",
                  help="discard 3rd coordinate of node positions",
                  action="store_const", const=2, dest="dim", default=3)

optparser.add_option("--shell", "-s",
                  help="interpret a 3d input mesh as a 2d mesh on a spherical shell",
                  action="store_const", const=True, dest="spherical_shell", default=False)

optparser.add_option("--internal-boundary", "-i", 
                  help="mesh contains internal faces - this option is required if you have assigned " +
                       "a physical boundary id to lines (2D) or surfaces (3D) that are not on the domain boundary",
                  action="store_const", const=True, dest="internal_faces", default=False)

(options, argv) = optparser.parse_args()

if len(argv)<1:
    optparser.print_help()
    sys.exit(1)

if argv[0][-4:]!=".msh":
    sys.stderr.write("Mesh filename must end in .msh\n")
    optparser.print_help()
    sys.exit(1)
    

basename=os.path.basename(argv[0][:-4])

mshfile=open(argv[0], 'r')

# Header section
assert(mshfile.readline().strip()=="$MeshFormat")
assert(mshfile.readline().strip()in["2 0 8", "2.1 0 8", "2.2 0 8"])
assert(mshfile.readline().strip()=="$EndMeshFormat")

# Nodes section
while mshfile.readline().strip() !="$Nodes":
    pass
nodecount=int(mshfile.readline())

if nodecount==0:
  sys.stderr.write("ERROR: No nodes found in mesh.\n")
  sys.exit(1)

if nodecount<0:
  sys.stderr.write("ERROR: Negative number of nodes found in mesh.\n")
  sys.exit(1)

dim=options.dim
nodefile_dim = dim
if options.spherical_shell:
  nodefile_dim -= 1

nodefile_linelist = []
for i in range(nodecount):
    # Node syntax
    line = mshfile.readline().split()
    # compare node id assigned by gmsh to consecutive node id (assumed by fluidity)
    if eval(line[0])!=i+1:
      print(line[0], i+1)
      sys.stderr.write("ERROR: Nodes in gmsh .msh file must be numbered consecutively.")
    nodefile_linelist.append( line[1:dim+1] )

assert(mshfile.readline().strip()=="$EndNodes")

# Elements section
assert(mshfile.readline().strip()=="$Elements")
elementcount=int(mshfile.readline())

# Now loop over the elements placing them in the appropriate buckets.
edges=[]
triangles=[]
tets=[]
quads=[]
hexes=[]

for i in range(elementcount):
    
    element=mshfile.readline().split()

    if (element[1]=="1"):
        edges.append(element[-2:]+[element[3]])
    elif (element[1]=="2"):
        triangles.append(element[-3:]+[element[3]])
    elif (element[1]=="3"):
        quads.append(element[-4:]+[element[3]])
    elif (element[1]=="4"):
        tets.append(element[-4:]+[element[3]])
    elif (element[1]=="5"):
        hexes.append(element[-8:]+[element[3]])
    elif(element[1]=="15"):
        # Ignore point elements
        pass
    else:
        sys.stderr.write("Unknown element type "+str(element[1])+'\n')
        sys.exit(1)

if len(tets) > 0:
  if len(hexes) > 0:
    sys.stderr.write("Warning: Mixed tet/hex mesh encountered - discarding hexes")
  if len(quads) > 0:
    sys.stderr.write("Warning: Mixed tet/quad mesh encountered - discarding quads")
elif len(triangles) > 0:
  if len(hexes) > 0:
    sys.stderr.write("Warning: Mixed triangle/hex mesh encountered - discarding hexes")
  if len(quads) > 0:
    sys.stderr.write("Warning: Mixed triangle/quad mesh encountered - discarding quads")

if len(tets)>0:
    dim=3
    loc=4
    node_order=[1, 2, 3, 4]
    elements=tets
    faces=triangles
    elefile=open(basename+".ele", "w")
    facefile=open(basename+".face", "w")

elif len(triangles)>0:
    dim=2
    loc=3
    node_order=[1, 2, 3]
    elements=triangles
    faces=edges
    elefile=open(basename+".ele", "w")
    facefile=open(basename+".edge", "w")

elif len(hexes)>0:
    dim=3
    loc=8
    node_order=[1, 2, 4, 3, 5, 6, 8, 7]
    elements=hexes
    faces=quads
    elefile=open(basename+".ele", "w")
    facefile=open(basename+".face", "w")

elif len(quads)>0:
    dim=2
    loc=4
    node_order=[1, 2, 4, 3]   # don't really know if this is right
    elements=quads
    faces=edges
    elefile=open(basename+".ele", "w")
    facefile=open(basename+".edge", "w")

else:
    sys.stderr.write("Unable to determine dimension of problem\n")
    sys.exit(1)

nodefile=open(basename+".node", 'w')
nodefile.write(str(nodecount)+" "+str(nodefile_dim)+" 0 0\n")
j=0
for i in range(nodecount):    
  j=j+1
  nodefile.write(" ".join( [str(j)] + nodefile_linelist[i] )+"\n")

nodefile.write("# Produced by: "+" ".join(argv)+"\n")
nodefile.close()

# Output ele file
elefile.write(str(len(elements))+" "+str(loc)+" 1\n")

for i, element in enumerate(elements):
    elefile.write(str(i+1)+" ")
    for j in node_order:
      elefile.write(" ".join(element[j-1:j])+" ")
    elefile.write(" ".join(element[-1:]))
    elefile.write(" "+"\n")

elefile.write("# Produced by: "+" ".join(sys.argv)+"\n")
elefile.close()

# Output ele or face file
if options.internal_faces:
  # make node element list
  ne_list = [set() for i in range(nodecount)]
  for i, element in enumerate(elements):
      element=[eval(element[j-1]) for j in node_order]
      for node in element:
        ne_list[node-1].add(i)

  # make face list, containing: face_nodes, surface_id, element_owner
  facelist=[]
  for face in faces:
    # last entry of face is surface-id
    face_nodes=[eval(node) for node in face[:-1]]
    # loop through elements around node face_nodes[0]
    for ele in ne_list[face_nodes[0]-1]:
      element=[eval(elements[ele][j-1]) for j in node_order]
      if set(face_nodes) < set(element):
        facelist.append(face+[str(ele+1)])
  
  facefile.write(str(len(facelist))+" 2\n")
  faces=facelist

else:
  facefile.write(str(len(faces))+" 1\n")

for i,face in enumerate(faces):
    facefile.write(str(i+1)+" "+" ".join(face)+"\n")

facefile.write("# Produced by: "+" ".join(sys.argv)+"\n")
facefile.close()

