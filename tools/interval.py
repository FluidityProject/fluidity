#!/usr/bin/env python3
#    Copyright (C) 2007 Imperial College London and others.
#
#    Please see the AUTHORS file in the main source directory for a full list
#    of copyright holders.
#
#    Prof. C Pain
#    Applied Modelling and Computation Group
#    Department of Earth Science and Engineering
#    Imperial College London
#
#    amcgsoftware@imperial.ac.uk
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#    USA

import sys
import os
from numpy import *
from optparse import OptionParser


parser=OptionParser(usage='usage: %prog [options] left right name',
                    add_help_option=True)

parser.add_option("--dx", dest="dx",type="float",
                  help="constant interval spacing.")

parser.add_option("--variable_dx", dest="variable_dx",type="string",
                  help="""Name of a file containing a python function val(X)
which defines the spacing (dx) at each point in the domain.""")

parser.add_option("--region_ids", dest="region_ids",type="string",
                  help="""Name of a file containing a python function val(X)
which defines the region_id at each point in the domain.""")

parser.add_option("--reverse", action="store_true", dest="reverse",
                  help="Reverse order of mesh.")

parser.add_option("--3d", action="store_const", const=3, dest="dim", default=1,
                  help="output 3d coordinates.")

(options, argv) = parser.parse_args()

try:
    left=float(argv[0])
    right=float(argv[1])
    dx=options.dx
    name=argv[2]
except:
    parser.print_help()
    sys.exit(1)

if(right<left):
    print("Error: right should be greater than left")
    parser.print_help()
    sys.exist(1)

if options.dx:
    dx=options.dx
    # This ensures the rightmost point is actually present.
    right=right+0.01*dx
    nodes=arange(left, right, dx)
elif options.variable_dx:
    exec(open(options.variable_dx).read())
    nodes=[left]
    while nodes[-1]<(right-float(str(val(right)))):
        nodes.append(nodes[-1]+ float(str(val(nodes[-1]))))
    #force last node to be equal to right as specified by user
    nodes.append(right)
else:
    dx=right-left
    # This ensures the rightmost point is actually present.
    right=right+0.01*dx
    nodes=arange(left, right, dx)

if(options.reverse):
    nodes = nodes[::-1]

dim=options.dim

eles=[(i+1,i+2) for i in range(len(nodes)-1)]

# Write the mesh file in Gmsh format
meshfile=open(name+".msh", "w")

meshfile.write("""$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
""")
meshfile.write("%d\n" % (len(nodes)))

for n in range(len(nodes)):
    meshfile.write("%d %.15f 0 0 \n" % (n + 1, nodes[n]))
meshfile.write("""$EndNodes
$Elements
""")
meshfile.write("%d\n" % (len(eles)+2))
# Mark first and last node
meshfile.write("1 15 2 1 1 1\n")
meshfile.write("2 15 2 2 %d %d\n" % (len(eles)+1, len(eles)+1))

if options.region_ids:
    exec(open(options.region_ids).read())
else:
    def val(X):
        return 1

for e in range(len(eles)):
    ele = eles[e]
    bid = val(0.5*(nodes[ele[0]-1]+nodes[ele[1]-1]))
    meshfile.write("%d 1 2 %d 1 %d %d\n" % (e + 3, bid, ele[0], ele[1]))

meshfile.write("$EndElements\n")
meshfile.close()
