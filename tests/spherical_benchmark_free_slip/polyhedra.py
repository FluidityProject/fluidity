from numpy import array
from numpy.linalg import norm
from math import sqrt

# 8-3-2017 s.kramer@imperial.ac.uk:
# the icosahedron() code is partly based on http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
# the subdivide() routine is based on http://prideout.net/blog/?p=44 , but completely rewritten to actually produce a valid mesh
# with unique vertices and not depend on euclid

def icosahedron():
    """Construct a 20-sided polyhedron"""
    t = (1.0 + sqrt(5.0)) / 2.0;
    verts = [(-1,  t,  0),
             ( 1,  t,  0),
             (-1, -t,  0),
             ( 1, -t,  0),

             ( 0, -1,  t),
             ( 0,  1,  t),
             ( 0, -1, -t),
             ( 0,  1, -t),

             ( t,  0, -1),
             ( t,  0,  1),
             (-t,  0, -1),
             (-t,  0,  1)]
    # now normalize them:
    s = sqrt(1 + t**2)
    verts = [tuple(x/s for x in vertex) for vertex in verts]

    faces = [(0, 11, 5),
             (0, 5, 1),
             (0, 1, 7),
             (0, 7, 10),
             (0, 10, 11),
             (1, 5, 9),
             (5, 11, 4),
             (11, 10, 2),
             (10, 7, 6),
             (7, 1, 8),
             (3, 9, 4),
             (3, 4, 2),
             (3, 2, 6),
             (3, 6, 8),
             (3, 8, 9),
             (4, 9, 5),
             (2, 4, 11),
             (6, 2, 10),
             (8, 6, 7),
             (9, 8, 1)]
    return verts, faces


def create_new_node(face, i, j, verts, new_nodes):
    a, b = face[i], face[j]
    if (a,b) in new_nodes:
        c = new_nodes[(a,b)]
    elif (b,a) in new_nodes:
        c = new_nodes[(b,a)]
    else:
        s = array(verts[a]) + array(verts[b])
        verts.append(s/norm(s))
        c = len(verts)-1
        new_nodes[(a,b)] = c
    return c

def subdivide(verts, faces):

    """Subdivide each triangle into four triangles, pushing verts to the unit sphere"""


    # map between pairs of existing verts (indices) and the new node (number) created inbetween
    new_nodes = {}

    # we instantiate the enumerated list first, as we only want to loop over already existing faces
    for fi, face in list(enumerate(faces)):

        # Create three new verts at the midpoints of each edge:
        i = create_new_node(face, 0, 1, verts, new_nodes)
        j = create_new_node(face, 1, 2, verts, new_nodes)
        k = create_new_node(face, 0, 2, verts, new_nodes)

        # Split the current triangle into four smaller triangles:
        faces.append((i, j, k))
        faces.append((face[0], i, k))
        faces.append((i, face[1], j))
        faces[fi] = (k, j, face[2])
