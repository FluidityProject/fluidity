import numpy
from ufl import *
from pyop2 import ffc_interface, op2
from pyop2.utils import uniquify

import state_types as fluidity_state
from state_types import *
from ufl_expr import *
from bcs import *

valuetype = numpy.float64
op2.init()

class FieldDict(dict):

    def __getitem__(self, key):
        return super(FieldDict, self).__getitem__(key)

# from http://www.toofishes.net/blog/python-cached-property-decorator/
class cached_property(object):
    '''A read-only @property that is only evaluated once. The value is cached
    on the object itself rather than the function or class; this should prevent
    memory leakage.'''
    def __init__(self, fget, doc=None):
        self.fget = fget
        self.__doc__ = doc or fget.__doc__
        self.__name__ = fget.__name__
        self.__module__ = fget.__module__

    def __get__(self, obj, cls):
        if obj is None:
            return self
        obj.__dict__[self.__name__] = result = self.fget(obj)
        return result

field2rank = {'ScalarField': 0,
              'VectorField': 1,
              'TensorField': 2}

field2element = {'ScalarField': FiniteElement,
                 'VectorField': VectorElement,
                 'TensorField': TensorElement}

dimloc2cell = {1: {1: 'vertex', 2: 'interval'},
               2: {3: 'triangle', 4: 'quadrilateral'},
               3: {4: 'tetrahedron', 6: 'hexahedron'}}

type2family = {'lagrangian': 'Lagrange'}

def ufl_cell(element):
    try:
        return dimloc2cell[element.dimension][element.quadrature.loc]
    except KeyError:
        raise RuntimeError("Elements of dimension %d and loc %d are not suppported" \
                % (element.dimension, element.loc))

def ufl_family(field):
    try:
        return "%s%s" % ("Discontinuous " if field.mesh.continuity < 0 else "", \
                         type2family[field.fluidity_element.type])
    except KeyError:
        raise RuntimeError("Elements of type %s are not supported" % field.fluidity_element.type)

def ufl_element(field):
    e = field.fluidity_element
    return field2element[field.description](ufl_family(field), ufl_cell(e), e.degree)

class Mesh(fluidity_state.Mesh):

    def __init__(self, *args, **kwargs):
        fluidity_state.Mesh.__init__(self, *args, **kwargs)
        self._boundaries_computed = False

    @cached_property
    def element_set(self):
        return op2.Set(self.element_count, "%s_elements" % self.name)

    @cached_property
    def node_set(self):
        return op2.Set(self.node_count, "%s_nodes" % self.name)

    @cached_property
    def element_node_map(self):
        return op2.Map(self.element_set, self.node_set, self.shape.loc, \
                self.ndglno - 1, "%s_elem_node" % self.name)

    def compute_boundaries(self):
        # FIXME: This check should use a @once decorator. However, I was confused about how to
        # implement it as there doesn't seem to be a canonical, safe example of once's implementation.
        if self._boundaries_computed:
            return
        self.faces.compute_boundaries(self)
        self._boundaries_computed = True

class Faces(fluidity_state.Faces):

    def __init__(self, surface_node_list, face_element_list, boundary_ids):
        fluidity_state.Faces.__init__(self, surface_node_list, face_element_list, boundary_ids)
        # The sets of elements on each boundary
        self.boundary_elem_sets = {}
        # The facet id of the facet on the boundary for a given element
        self.boundary_facets = {}
        # The mapping from elements on each boundary to the global element number
        self.boundary_maps = {}
        # The sets of nodes on each boundary
        self.surface_elem_sets = {}
        # A mapping from the elements on a boundary to the global node numbers
        self.surface_elements_to_nodes_maps = {}

    @cached_property 
    def boundary_list(self):
        return uniquify(self.boundary_ids)
    
    def compute_boundaries(self, mesh):
        # Get face_list CSR data structures. It's not really a CSR matrix.
        row_ptr = self.face_list.indptr
        col_idx = self.face_list.indices
        val     = self.face_list.data
        
        # Compute the list of (element, surface_element, local_facet) 
        # tuples on the boundary
        boundary_elements = []
        for element in xrange(len(row_ptr)-1):
            begin, end = row_ptr[element], row_ptr[element+1]
            neighbour_elements = col_idx[begin:end]
            facet_ids = val[begin:end]
            for local_facet, neighbour in enumerate(neighbour_elements):
                if neighbour < 0:
                    # Facet is on the mesh boundary
                    surface_element = facet_ids[local_facet]-1 # Facet numbers offset by 1
                    boundary_elements.append((element, surface_element, local_facet))
        
        # Initialise structures to hold data underlying the OP2 data structures
        element_count = {}
        boundary_faces = {}
        elem_node_mapping = {}
        for boundary in self.boundary_list:
            element_count[boundary] = 0
            boundary_faces[boundary] = []
            elem_node_mapping[boundary] = []

        # Count the number of elements in each set, create the element->node mapping, 
        # and store the local facet id for each element on the boundary
        for element, surface_element, local_facet in boundary_elements:
            boundary = self.boundary_ids[surface_element]
            boundary_faces[boundary].append(local_facet) # May need transforming for UFC
            elem_node_mapping[boundary] += mesh.ele_nodes(element)
            element_count[boundary] += 1

        # Find the surface elements in each boundary and their global nodes
        surface_element_count = {}
        surface_elements = {}
        for boundary in self.boundary_list:
            surface_element_count[boundary] = 0
            surface_elements[boundary] = []
        for surface_element, boundary in enumerate(mesh.faces.boundary_ids):
            surface_element_count[boundary] += 1
            surface_elements[boundary].append(surface_element)
        
        surface_elements_to_nodes = {}
        for boundary, surface_elements in surface_elements.iteritems():
            surface_elements_to_nodes[boundary] = []
            for surface_element in surface_elements:
                local_nodes = mesh.faces.surface_mesh.ele_nodes(surface_element)
                # Global node number is in Fortran numbering
                global_nodes = [ mesh.faces.surface_node_list[n]-1 for n in local_nodes ]
                surface_elements_to_nodes[boundary] += global_nodes

        # Construct OP2 data structures on top of this information
        for boundary in self.boundary_list:
            self.boundary_elem_sets[boundary] = \
                op2.Set(element_count[boundary])
            self.surface_elem_sets[boundary] = \
                op2.Set(surface_element_count[boundary])
            self.boundary_maps[boundary] = \
                op2.Map(self.boundary_elem_sets[boundary], 
                        mesh.node_set, mesh.shape.loc, 
                        elem_node_mapping[boundary])
            self.surface_elements_to_nodes_maps[boundary] = \
                op2.Map(self.surface_elem_sets[boundary], 
                        mesh.node_set, mesh.faces.surface_mesh.shape.loc, 
                        surface_elements_to_nodes[boundary])
            self.boundary_facets[boundary] = \
                op2.Dat(self.boundary_elem_sets[boundary], 1, 
                        boundary_faces[boundary], numpy.int32)

    def check_boundary(self, boundary):
        if boundary not in self.boundary_list:
            raise ValueError("Boundary %d not in boundary list." % boundary)

class FieldCoefficient(Coefficient):
    """Coefficient derived from a Fluidity field."""

    @property
    def value_shape(self):
        return (self.mesh.shape.dimension,)*self.rank() or 1

    @property
    def element_set(self):
        return self.mesh.element_set

    @cached_property
    def node_set(self):
        return op2.Set(self.node_count, "%s_nodes" % self.name)

    @cached_property
    def dat(self):
        return op2.Dat(self.node_set, self.value_shape, \
                self.val, valuetype, self.name)

    @cached_property
    def element_node_map(self):
        return op2.Map(self.mesh.element_set, self.node_set, self.mesh.shape.loc, \
                self.mesh.ndglno - 1, "%s_elem_node" % self.name)

    def temporary_dat(self, name=None):
        return op2.Dat(self.node_set, self.value_shape, \
                numpy.zeros(self.node_count), valuetype, name)

class ScalarField(FieldCoefficient, fluidity_state.ScalarField):

    def __init__(self,n,v,ft,op,uid,mesh=None,element=None,count=None):
        fluidity_state.ScalarField.__init__(self, n, v, ft, op, uid, mesh)
        FieldCoefficient.__init__(self, element or ufl_element(self), count)

    @property
    def fluidity_element(self):
        return fluidity_state.ScalarField.shape(self)

    def _reconstruct(self, element, count):
        # This code is class specific
        return ScalarField(self.name, self.val, self.field_type, self.option_path,
                self.uid, self.mesh, element, count)

class VectorField(FieldCoefficient, fluidity_state.VectorField):

    def __init__(self,n,v,ft,op,dim,uid,mesh=None,element=None,count=None):
        fluidity_state.VectorField.__init__(self, n, v, ft, op, dim, uid, mesh)
        FieldCoefficient.__init__(self, element or ufl_element(self), count)

    @property
    def fluidity_element(self):
        return fluidity_state.VectorField.shape(self)

    def _reconstruct(self, element, count):
        # This code is class specific
        return VectorField(self.name, self.val, self.field_type, self.option_path,
                self.dimension, self.uid, self.mesh, element, count)

class TensorField(FieldCoefficient, fluidity_state.TensorField):

    def __init__(self,n,v,ft,op,dim0,dim1,uid,mesh=None,element=None,count=None):
        fluidity_state.TensorField.__init__(self, n, v, ft, op, dim0, dim1, uid, mesh)
        FieldCoefficient.__init__(self, element or ufl_element(self), count)

    @property
    def fluidity_element(self):
        return fluidity_state.TensorField.shape(self)

    def _reconstruct(self, element, count):
        # This code is class specific
        return TensorField(self.name, self.val, self.field_type, self.option_path,
                self.dimension[0], self.dimension[1], self.uid,
                self.mesh, element, count)

from solving import solve
