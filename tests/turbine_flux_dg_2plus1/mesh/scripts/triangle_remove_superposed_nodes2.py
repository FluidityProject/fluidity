#!/usr/bin/env python3
import sys
import triangle
import copy
import numpy
from sets import Set
#input surface_id, filename



# This scripts removes superposed nodes, including elements and edges, in a triangle file set. 
# It should be called before seperate_internal_boundary2.py to avoid empty elements
# 
# Changelog: 
# 4.4.10 Initial version
# 5.5.10 Some more detailed help
# 5.5.10 Version 2 of this script can handle .edge with (more than just the boundary id) attributes.
# 


# Here is an examle mesh4.geo file for this script:

# Point(1) = {0, 0, 0, 2};
# Point(2) = {1, 0, 0, 2};
# Point(3) = {1, 1, 0, 2};
# Point(4) = {0, 1, 0, 2};
# Point(5) = {0.5, 0, 0, 2};
# Point(6) = {0.5, 1, 0, 2};
# Point(7) = {0.500001, 0, 0, 2};
# Point(8) = {0.500001, 1, 0, 2};
# Point(9) = {0.4, -0.1, 0, 2};
# Point(10) = {0.4, 1.1, 0, 2};
# 
# 
# Line(1) = {4, 1};
# Line(2) = {1, 9};
# Line(3) = {9, 5};
# Line(4) = {5, 6};
# Line(9) = {6, 10};
# Line(10) = {10, 4};
# 
# Line(5) = {8, 7};
# Line(6) = {7, 2};
# Line(7) = {2, 3};
# Line(8) = {3, 8};
# 
# Physical Line(20) = {1};
# Physical Line(21) = {2};
# Physical Line(22) = {3};
# Physical Line(23) = {4};
# Physical Line(28) = {9};
# Physical Line(29) = {10};
# 
# Physical Line(24) = {5};
# Physical Line(25) = {6};
# Physical Line(26) = {7};
# Physical Line(27) = {8};
# 
# Line Loop(10) = {4, 9, 10, 1, 2, 3};
# Line Loop(11) = {8, 5, 6, 7};
# 
# Plane Surface(11) = {10};
# Plane Surface(12) = {11};
# Physical Surface(12) = {11, 12};

# Example command line:

#  triangle_add_edgeowner.py mesh4
#  triangle_remove_superposed_nodes2.py mesh4_edgow 23

########################################################################################################
if not len(sys.argv)==3:
        print "Usage: seperate_internal_boundary.py file boundary_id"        
        print ""
        print "output fixed .ele and .node file with removed superposed nodes on the given boundary"
        print ""
        print "The nodes have to have a distance less than 1e-8 to be considered as superposed"
        print "The outout files will be have the suffix nosup"

        exit()

filename=sys.argv[1]
boundary_id=int(sys.argv[2])
debug=2

triangle.read_nodefile(filename+'.node')
if triangle.dim!=2:
        print "Only 2 dim meshes supported so far"
triangle.read_edgefile(filename+'.edge')
triangle.read_elefile(filename+'.ele')

edgein_to_test=triangle.edges_with_bid(boundary_id)
nodeid_to_test=Set()
for edge in edgein_to_test:
        nodeid_to_test=nodeid_to_test.union(Set(triangle.edges[edge]))
print "Checking nodes ", nodeid_to_test

# Look for superposed nodes at the position if nodeid
eleid_del=Set()
edgeid_del=Set()
nodeid_del=Set()
for nodeid in nodeid_to_test:
        for k in range(0,len(triangle.nodes)):
                if k+1==nodeid:
                        continue
                if ((triangle.nodes[k][0]-triangle.nodes[nodeid-1][0])**2+(triangle.nodes[k][1]-triangle.nodes[nodeid-1][1])**2)<1e-8:
                        print "\tFound two nodes at the same position with ids: ", nodeid, " ", triangle.nodes[nodeid-1] , " and ", k+1, triangle.nodes[k], " ", ". The latter node will be deleted."
                        nodeid_del.add(k+1) # remove node with id k
                        print "Remove unnecessary node ", k+1
                        # replace nodes k by nodeid in the element list
                        for index, ele in enumerate(triangle.eles):
                                 triangle.eles[index]=[nodeid if x==k+1 else x for x in triangle.eles[index]]
                                 if (len(Set(triangle.eles[index]))!=len(triangle.eles[index])): # We dont need elements with nodes k+1 and nodeid
                                        eleid_del.add(index+1)
                                        print "Remove unnecessary element ", index+1
                        # replace nodes k by nodeid in the edge list
                        for index, edge in enumerate(triangle.edges):
                                 triangle.edges[index]=[nodeid if x==k+1 else x for x in triangle.edges[index]]
                                 if (len(Set(triangle.edges[index]))!=len(triangle.edges[index])): # We dont need edges with nodes k+1 and nodeid
                                        edgesid_del.add(index+1)
                                        print "Remove unnecessary edge ", index+1
                        
nodeid_delar=list(nodeid_del)
edgeid_delar=list(edgeid_del)
eleid_delar=list(eleid_del)
nodeid_delar.sort()
edgeid_delar.sort()
eleid_delar.sort()
nodeid_delar.reverse()
edgeid_delar.reverse()
eleid_delar.reverse()
print "Nodes to delete: ", nodeid_delar
print "Edges to delete: ", edgeid_delar
print "Elements to delete: ", eleid_delar

for e in edgeid_delar:
        triangle.delete_edgeid(e)
for e in eleid_delar:
        triangle.delete_eleid(e)
# Note: Delete nodes at last to make sure that they are not used by any elements or edges
for n in nodeid_delar:
        triangle.delete_nodeid(n)

if debug>0: 
        print "save node file as ", filename, "_nosup2.node"
triangle.save_nodefile(triangle.nodes, 2, filename+"_nosup2.node")
if debug>0: 
        print "save ele file as ", filename, "_nosup2.ele"
triangle.save_elefile(triangle.eles, triangle.region_ids, filename+"_nosup2.ele")
if debug>0: 
        print "save edge file as ", filename, "_nosup2.edge"
triangle.save_edgefile2(triangle.edges, triangle.surface_ids, triangle.edge_attribute1, filename+"_nosup2.edge")






