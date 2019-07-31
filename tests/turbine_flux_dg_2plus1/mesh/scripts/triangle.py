#!/usr/bin/env python3
import sys
from sets import Set

edges=[]
surface_ids=[] # surface ids of the edges
nodes=[]
eles=[]
region_ids=[] # region ids of the elements
edge_attribute1=[] # first edge attribute 
dim=-1
debug=0


def read_nodefile(filename):
   f = open(filename, 'r')
   global nodes, dim
   [nodes, dim]=pnode(f)
   f.close()

def read_elefile(filename):
   f = open(filename, 'r')
   global eles
   eles=peles(f)
   f.close()

def read_edgefile(filename):
   f = open(filename, 'r')
   global edges, surface_ids, edge_attribute1
   [edges, surface_ids, edge_attribute1]=pedges(f)
   f.close()

def save_nodefile(nodes, dim, filename):
   f = open(filename, 'w')

   f.write(`len(nodes)`+' '+`dim`+' 0 0\n')
   k=1
   for node in nodes:
        f.write(`k`+' ' + ''.join([`n`+' ' for n in node])+'\n')
        k=k+1
   f.write('# Created with triangle.py')
   f.close()

def save_elefile(eles, par1, filename):
   f = open(filename, 'w')
   f.write(`len(eles)`+' '+`len(eles[0])`+' 1\n')
   if len(par1)!=len(eles):
        print "Error in save_elefile(): Length of element list is ", len(eles), "  but length of region_id list is: ", len(par1)
        exit()
   k=1
   for ele in eles:
        f.write(`k`+' ' + ''.join([`e`+' ' for e in ele])+`par1[k-1]`+'\n')
        k=k+1
   f.write('# Created with triangle.py')
   f.close()


def save_edgefile(edges, atr1, filename):
   f = open(filename, 'w')
   f.write(`len(edges)`+' 1\n')
   if len(atr1)!=len(edges):
        print "Error in save_edgefile(): Length of edge list is ", len(edges), "  but length of atr1 list is: ", len(par1)
        exit()
   k=1
   for edge in edges:
        f.write(`k`+' ' + ''.join([`e`+' ' for e in edge])+`atr1[k-1]`+'\n')
        k=k+1
   f.write('# Create with triangly.py')
   f.close()

def save_edgefile2(edges, atr1, atr2, filename):
   f = open(filename, 'w')
   f.write(`len(edges)`+' 2\n')
   if len(atr1)!=len(edges) or len(atr2)!=len(edges):
        print "Error in save_edgefile(): Length of edge list is ", len(edges), "  but length of atr1, atr2 list is: ", len(atr1), ', ', len(atr2)
        exit()
   k=1
   for edge in edges:
        f.write(`k`+' ' + ''.join([`e`+' ' for e in edge])+`atr1[k-1]`+' '+ `atr2[k-1]`+'\n')
        k=k+1
   f.write('# Create with triangle.py')
   f.close()

def pnode(fnode):
    """ Parse the .node file
    and returns the correspondind
    list of the nodes
    """

    nodes = []
    count = 0
    dim = 0
    nbnodes = 0

    for line in fnode:

        if line[0]=='#':
            continue
        else:
            parts = line.split()
            if len(parts)==0:
                continue
            if count == 0:
                nbnodes = int(parts[0])
                dim = int(parts[1])
                count += 1
                if dim not in [2,3]:
                    print "Error dim in .node file!"

            else:
                if dim == 3:
                    nodes.append(([float(x) for x in parts[1:]]))
                else:
                    nodes.append(([float(x) for x in parts[1:3]]))
    if (len(nodes) != nbnodes):
        print "Error: actual number of nodes does not corresponds " + str(len(nodes)), nbnodes
        sys.exit(2)
                                 
    return nodes, dim

def pedges(fedge):
    """ Parse the .edge file and returns the corresponding
    list of edges
    """

    edges = []
    surface_ids=[]
    attribute1=[]
    count = 0
    nbedges = 0
    boundary = 0
    for line in fedge:
        if line[0]=='#':
            continue
        else:
            parts = line.split()
            if len(parts)==0:
                continue
        
            if len(parts)==0:
                continue
            if count == 0:
                nbedges = int(parts[0])
                boundary = int(parts[1])
                count += 1
                if boundary>2:
                    print "Error: too many attributes in .edge file!"

            else:
                 edges.append(([int(x) for x in parts[1:3]]))
                 if boundary >= 1:
                   surface_ids.append(int(parts[3]))
                 if boundary >= 2:
                   attribute1.append(int(parts[4]))

    if (len(edges) != nbedges):
        print "Error: actual number of edges does not corresponds " + str(len(edges)), str(nbedges)
        sys.exit(2)

    return [edges, surface_ids, attribute1]

def peles(fele):
    """ Parse the .ele file and returns the corresponding
    list of elements
    """

    eles = []
    count = 0
    nbeles = 0
    for line in fele:

        if line[0]=='#':
            continue
        else:
            parts = line.split()

            if count == 0:
                nbeles = int(parts[0])
                nodesperele = int(parts[1])
                region = int(parts[2])
                count += 1

            else:
                eles.append(([int(x) for x in parts[1:nodesperele+1]]))
                if region==1:
                        region_ids.append(int(parts[len(parts)-1]))

    if (len(eles) != nbeles):
        print "Error: actual number of eles does not corresponds " + str(len(eles)), str(nbeles)
        sys.exit(2)

    return eles

# Returns array of edge indices whose boundary marker ids equal bid
def edges_with_bid(bid):
   global edges
   ret=[]
   
   for i in range(0,len(surface_ids)):
        if surface_ids[i] == bid:
                ret.append(i)
   return ret

# Returns array of ele indices which have the given edge id
def ele_with_edgeids(edgeids):
   global edges, eles, nodes
   ret=[]
   nodeids=edges[edgeids][0:len(edges[edgeids])] # get nodeids of edgeid
   # loop through all elements and check if they include all these nodes
   i=0
   for ele in eles:
        if Set(nodeids).issubset(Set(ele[0:len(ele)])):
                ret.append(i)
        i=i+1
   if len(ret)>2 :
        print "There are ", len(ret), " elements adjacent to face with id ", edgeids, ". Only 2 were expected."
        exit()
   return ret

# Returns array of ele indices which have the given node ids
def ele_with_nodeids(nodeids):
   global eles, nodes
   ret=[]
   # loop through all elements and check if they include all these nodes
   i=0
   for ele in eles:
        if Set(nodeids).issubset(Set(ele[0:len(ele)])):
                ret.append(i)
        i=i+1
   if len(ret)!=2 :
        print "There are ", len(ret), " elements adjacent to nodes with ids ", nodeids, ". Only 2 were expected."
        exit()
   return ret

# returns the edges with nodeids
def edgeid_from_nodeids(nodeids):
        global edges
        ret=[]
        i=0
        for edge in edges:
                if Set(nodeids).issubset(Set(edge[0:len(edge)])):
                        ret.append(i)
                i=i+1
        return ret

# returns the edge with boundary_id and wit node nodeid but not with id edgeid
def get_partner_edge(edgein, nodeid, boundary_id):
        global edges
        cands=edges_with_bid(boundary_id)
        for edge in cands:
                if edge!=edgein and (nodeid in edges[edge]):
                        return edge

def get_boundaryid_from_edgeid(edgeid):
        global surface_ids
        return surface_ids[edgeid]

def has_ele_edge_on_boundaryid(elein, nodeid, boundary_id):
        global eles, debug
        for n in eles[elein]:
                if n==nodeid:
                        continue
                edge=edgeid_from_nodeids([n,nodeid])
                if edge!=[] and get_boundaryid_from_edgeid(edge[0])==boundary_id:
                        if debug>3:
                                print "Found boundary edge with id", edgeid_from_nodeids([n,nodeid])[0]+1
                        return True
        return False
                
def get_eles_on_ele_side(elein, nodeid, edgein, boundary_id):
        global eles, edges
        if not has_ele_edge_on_boundaryid(elein, nodeid, boundary_id):
                print "Error in get_eles_on_ele_side: Given Element has no edge with secified boundary id"
                sys.exit()
        ret=[elein]
        thirdnode=Set(eles[elein]).difference(Set(edges[edgein])).pop()
        while True: 
                if debug>3:
                        print "Found third node of current element wit index: ", thirdnode
                        print "Looking for neighbour with nodes ", thirdnode, nodeid
                partner_eleins=ele_with_nodeids([thirdnode, nodeid])
                if partner_eleins[0]==elein:
                        partner_elein=partner_eleins[1]
                else:
                        partner_elein=partner_eleins[0]
                if debug>3:
                        print "Found neighbour with eleid", partner_elein+1
                ret.append(partner_elein)
                if has_ele_edge_on_boundaryid(partner_elein, nodeid, boundary_id):
                        if debug>3:
                                print "Reached the end of the element loop"
                        return ret
                if debug>3:
                        print "New elemnt has no desired boundary edge, loop again"
                elein=partner_elein
                thirdnode=Set(eles[elein]).difference(Set([nodeid, thirdnode])).pop()


# deletes the node with id nodeid and updates the elelist and edgelist 
def delete_nodeid(nodeid):
        global eles, edges, nodes
        for index, ele in enumerate(eles):
                if nodeid in ele:
                        print "Error in delete_nodeid(", nodeid, "). This node is still in use in element ", index+1, ": ", ele, "."
                        exit()
                eles[index]=[x if x<nodeid else x-1 for x in eles[index]]
        for index, edge in enumerate(edges):
                if nodeid in edge:
                        print "Error in delete_nodeid(", nodeid, "). This node is still in use in edge ", index+1, "."
                        exit()
                edges[index]=[x if x<nodeid else x-1 for x in edges[index]]
        nodes.pop(nodeid-1)

# deletes the ele with id eleid 
def delete_eleid(eleid):
        global eles, region_ids
        eles.pop(eleid-1)
        region_ids.pop(eleid-1)

# deletes the edge with id edgeid
def delete_edgeid(edgeid):
        global edges, surface_ids
        edges.pop(edgeid-1)
        surface_ids.pop(nodeid-1)
        edge_attribute1.pop(nodeid-1)
