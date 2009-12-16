# Author: Daryl Harrison

# Enthought library imports.
from enthought.traits.api import HasTraits, Instance, String
from enthought.tvtk.api import tvtk

# Local imports.
from os import path

######################################################################
# `TriangleWriter` class
######################################################################
class TriangleWriter(HasTraits):
    """
    
    < Description here >
    """

    # The version of this class.  Used for persistence.
    __version__ = 0

    # Grid to write.
    grid = Instance(tvtk.UnstructuredGrid, allow_none=False)
    basename = String

    def __init__(self, input_grid, output_basename):
        # Expects basename without file extension. Appropriate .node, .face, .ele files will be created.

        self.grid = input_grid
        self.basename = output_basename

    def write(self):
        self.write_node_file()

        cell_type = self.grid.get_cell(0).cell_type

        if (cell_type == 5):
            # triangles - write .face file
            self.write_face_file()
        else:
            # tetrahedra - write .ele file
            self.write_ele_file()

    ########################################
    # File writing methods.

    def write_node_file(self):
        # Data arrays: [attribute 0] ... [attribute n] [boundary marker]
        f = open(self.basename+'.node', 'w')

        points = self.grid.number_of_points
        dimensions = len(self.grid.points[0])

        number_of_arrays = self.grid.point_data.number_of_arrays
        if (number_of_arrays > 0):
            # True if last of arrays in point_data is an IntArray
            boundary_marker = (self.grid.point_data.get_array(number_of_arrays-1).data_type == 6)
            # Number of arrays in point_data (subtract 1 if boundary_marker is True)
            attributes = (number_of_arrays-1 if boundary_marker else number_of_arrays)
        else:
            boundary_marker = False
            attributes = 0

        f.write(`points`+" "+`dimensions`+" "+`attributes`+ " "+`int(boundary_marker)`)

        arrays = []
        for i in range(number_of_arrays):
            arrays.append(self.grid.point_data.get_array(i))

        for i in range(points):
            f.write("\n"+`i+1`+" ") # write point number (starting from 1)

            point = self.grid.get_point(i)
            for j in range(dimensions):
                f.write(`point[j]`+ " ")

            for j in range(number_of_arrays):
                f.write(`arrays[j].get_value(i)`+" ")

        f.write("\n# Produced by: TriangleWriter in MayaVi2")


    def write_face_file(self):
        # .face file contains triangles
        # Data arrays: [boundary marker]
        f = open(self.basename+'.face', 'w')

        faces = self.grid.number_of_cells
        # Assuming that the first of any arrays is the boundary marker, and it must be an IntArray
        boundary_marker = (self.grid.cell_data.number_of_arrays > 0 and self.grid.cell_data.get_array(0).data_type == 6)

        f.write(`faces`+" "+`int(boundary_marker)`)

        if (boundary_marker):
            boundary_marker_array = self.grid.cell_data.get_array(0)

        for i in range(faces):
            f.write("\n"+`i+1`+" ") # write face number (starting from 1)

            cell = self.grid.get_cell(i)
            for j in range(3): # nodes per triangle
                point = cell.get_point_id(j)
                f.write(`point+1`+ " ")

            if (boundary_marker):
                f.write(`boundary_marker_array.get_value(i)`)

        f.write("\n# Produced by: TriangleWriter in MayaVi2")


    def write_ele_file(self):
        # .ele file contains tetrahedrons
        # Data arrays: [attribute 0] ... [attribute n]
        f = open(self.basename+'.ele', 'w')

        tetrahedra = self.grid.number_of_cells
        nodes_per_tetrahedron = len(self.grid.get_cell(0).point_ids)
        attributes = self.grid.cell_data.number_of_arrays

        f.write(`tetrahedra`+" "+`nodes_per_tetrahedron`+" "+`attributes`)

        attribute_arrays = []
        for i in range(attributes):
            attribute_arrays.append(self.grid.cell_data.get_array(i))

        for i in range(tetrahedra):
            f.write("\n"+`i+1`+" ") # write tetrahedron number (starting from 1)

            cell = self.grid.get_cell(i)
            for j in range(nodes_per_tetrahedron):
                point = cell.get_point_id(j)
                f.write(`point+1`+ " ")

            for j in range(attributes):
                f.write(`attribute_arrays[j].get_value(i)`+" ")

        f.write("\n# Produced by: TriangleWriter in MayaVi2")