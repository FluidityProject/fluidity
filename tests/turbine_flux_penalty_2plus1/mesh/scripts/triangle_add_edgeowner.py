#!/usr/bin/env python3
import sys
import triangle
import copy
import numpy
from sets import Set
#input surface_id, filename

# 5.5.2010: this script adds a new attribute to the .edge file which holds the "owner" element number of this edge


# Here is an examle geo file for this script:

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




########################################################################################################
def nodedupl_recursion(elein, edgein, nodeid, boundary_id):
     global copy_eles, copy_edges, copy_nodes, debug, copy_surface_ids, copy_surface_id, copy_surfaceowner_ids, copy_region_ids

     next_edgein=triangle.get_partner_edge(edgein, nodeid, boundary_id)

     if next_edgein==None:
        print "Reached one end of the surface boundary."
        return
     if debug>1: 
        print "Lets loop around nodeid", nodeid, " starting with ele", elein+1, " with boundary edge ", edgein+1, " until we reach the next surface edge with id ", next_edgein+1
     next_elein_list=triangle.get_eles_on_ele_side(elein, nodeid, edgein, boundary_id)

     if debug>1: 
             print "Duplicate edge ", next_edgein +1
     copy_edges.append(triangle.edges[next_edgein])
     copy_surface_ids.append(new_surface_id)
     copy_surfaceowner_ids.append(next_elein_list[len(next_elein_list)-1]+1) # update copy_surfaceowner_ids for the new edge
     # update copy_surfaceowner_ids for the old edge
     if triangle.ele_with_edgeids(next_edgein)[0]==next_elein_list[len(next_elein_list)-1]:
        copy_surfaceowner_ids[next_edgein]=triangle.ele_with_edgeids(next_edgein)[1]+1
     else:
        copy_surfaceowner_ids[next_edgein]=triangle.ele_with_edgeids(next_edgein)[0]+1

     if (triangle.edges[next_edgein][0]==nodeid):
        next_nodeid=triangle.edges[next_edgein][1]
     else:
        next_nodeid=triangle.edges[next_edgein][0]
     nodedupl_recursion(next_elein_list[len(next_elein_list)-1], next_edgein, next_nodeid, boundary_id)
     


########################################################################################################
if not len(sys.argv)==2:
        print "Usage: seperate_internal_boundary.py file"        
        print ""
        print "output fixed .edge, .ele and .node file with new edge attribute holding the element owner of the edge. "
        print ""
        print "The outout files will be have the suffix edgow"
        exit()

filename=sys.argv[1]
debug=2

triangle.read_nodefile(filename+'.node')
if triangle.dim!=2:
        print "Only 2 dim meshes supported so far"
triangle.read_edgefile(filename+'.edge')
triangle.read_elefile(filename+'.ele')

copy_eles=copy.deepcopy(triangle.eles)
copy_region_ids=copy.deepcopy(triangle.region_ids)
copy_edges=copy.deepcopy(triangle.edges)
copy_surface_ids=copy.deepcopy(triangle.surface_ids)
copy_surfaceowner_ids=[-1 for i in range(0,len(triangle.surface_ids))] # Will store the elemed id for each surface edge
copy_nodes=copy.deepcopy(triangle.nodes)

# Now assign the surfaceowner_id to the external boundaries
for e in range(0,len(copy_surfaceowner_ids)):
        if copy_surfaceowner_ids[e]>=0:
                print "Internal Error. Ask simon.funke@gmail.com!"
                exit()                 
        if len(triangle.ele_with_edgeids(e))!=1:
                print "Error Found internal boundary!"
                exit()
        copy_surfaceowner_ids[e]=triangle.ele_with_edgeids(e)[0]+1

if debug>0: 
        print "save node file as ", filename, "_edgow.node"
triangle.save_nodefile(copy_nodes, 2, filename+"_edgow.node")
if debug>0: 
        print "save ele file as ", filename, "_edgow.ele"
triangle.save_elefile(copy_eles, copy_region_ids, filename+"_edgow.ele")
if debug>0: 
        print "save edge file as ", filename, "_edgow.edge"
triangle.save_edgefile2(copy_edges, copy_surface_ids, copy_surfaceowner_ids, filename+"_edgow.edge")

