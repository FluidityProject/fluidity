#!/usr/bin/python

import sys
import re

pts = []
tet = []
tets = []
tet_no = 0
fh = open ( sys.argv[1] )

def append_tet ( tet ):
    if tet == []:
        return
    if len(tet) != 4:
        print tet
        print ( "Wrong number of points in tet %d!" % len ( tets ) )
        sys.exit ( 1 )
    tets.append ( tet )

def write_vtk_header ():
    print "# vtk DataFile Version 2.0"
    print "VTK LEGACY FILE FORMAT ASCII - tetrahedral sample file for H5FED test"
    print "ASCII"
    print "DATASET UNSTRUCTURED_GRID"
    print

def write_vtk_points_double ( pts ):
    print "POINTS %d DOUBLE" % ( len(pts) )
    for pt in pts:
        print "%s %s %s" % ( pt[0], pt[1], pt[2] )
    print

def write_vtk_cells ( cells ):
    num_cells = len ( cells )
    size = num_cells
    for cell in cells:
        size += len ( cell )
    print "CELLS %d %d" % ( num_cells, size )
    for cell in cells:
        print "%d\t" % (len(cell)),
        for el in cell:
            print "%d " % (el),
        print
    print

def write_vtk_cell_types ( cells ):
    print "CELL_TYPES %d" % ( len(cells) )
    for cell in cells:
        print "10"
    print

def write_vtk_fooder ( cells ):
    print "CELL_DATA %d" % ( len(cells) )
    print "SCALARS cell_attribute_data float 1"
    print "LOOKUP_TABLE default"
    for i in range(len(cells)):
        print "%d" % (i)

for line in fh:
    line = line.rstrip('\n')
    if line[0:3] == "TET":
        append_tet ( tet )
        tet = []
    else:
        pt = re.sub ( ' ', '', line).split(',')
        if not pt in pts:
            pts.append ( pt )
        tet.append ( pts.index ( pt ) )

append_tet ( tet )
write_vtk_header()
write_vtk_points_double ( pts )
write_vtk_cells ( tets )
write_vtk_cell_types ( tets )
write_vtk_fooder ( tets )
