from numpy import array
from numpy.linalg import norm

# 8-3-2017 s.kramer@imperial.ac.uk:
# the icosahedron() code is copied from http://prideout.net/blog/?p=44 - "I consider this code to be on the public domain" (sic.)
# the subdivide() routine is based on it, but completely rewritten to actually produce a valid mesh
# with unique vertices and not depend on euclid

def icosahedron():

    """Construct a 20-sided polyhedron"""
    faces = [ \
        (0,1,2),
        (0,2,3),
        (0,3,4),
        (0,4,5),
        (0,5,1),
        (11,6,7),
        (11,7,8),
        (11,8,9),
        (11,9,10),
        (11,10,6),
        (1,2,6),
        (2,3,7),
        (3,4,8),
        (4,5,9),
        (5,1,10),
        (6,7,2),
        (7,8,3),
        (8,9,4),
        (9,10,5),
        (10,6,1) ]
    verts = [ \
        ( 0.000,  0.000,  1.000 ),
        ( 0.894,  0.000,  0.447 ),
        ( 0.276,  0.851,  0.447 ),
        (-0.724,  0.526,  0.447 ),
        (-0.724, -0.526,  0.447 ),
        ( 0.276, -0.851,  0.447 ),
        ( 0.724,  0.526, -0.447 ),
        (-0.276,  0.851, -0.447 ),
        (-0.894,  0.000, -0.447 ),
        (-0.276, -0.851, -0.447 ),
        ( 0.724, -0.526, -0.447 ),
        ( 0.000,  0.000, -1.000 ) ]

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
