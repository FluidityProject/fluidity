import numpy
from ufl import *
from pyop2 import ffc_interface, op2

import state_types as fluidity_state
from state_types import *
from ufl_expr import *

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

def ufl_cell(shape):
    try:
        return dimloc2cell[shape.dimension][shape.quadrature.loc]
    except KeyError:
        raise RuntimeError("Elements of dimension %d and loc %d are not suppported" \
                % (shape.dimension, shape.loc))

def ufl_family(field):
    try:
        return "%s%s" % ("Discontinuous " if field.mesh.continuity < 0 else "", \
                         type2family[field.shape().type])
    except KeyError:
        raise RuntimeError("Elements of type %s are not supported" % field.shape().type)

def ufl_element(field):
    s = field.shape()
    return field2element[field.description](ufl_family(field), ufl_cell(s), s.degree)

class Mesh(fluidity_state.Mesh):

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

class FieldCoefficient(Coefficient):
    """Coefficient derived from a Fluidity field."""

    def __init__(self, field, element=None, count=None):
        field.rank = field2rank[field.description]
        self._field = field
        super(FieldCoefficient, self).__init__(element or ufl_element(field), count)

    @property
    def field(self):
        return self._field

    @property
    def mesh(self):
        return self.field.mesh

    @property
    def name(self):
        return self.field.name

    @property
    def value_shape(self):
        return (self.field.mesh.shape.dimension,)*self.field.rank or 1

    @property
    def element_set(self):
        return self.field.mesh.element_set

    @cached_property
    def node_set(self):
        return op2.Set(self.field.node_count, "%s_nodes" % self.name)

    @cached_property
    def dat(self):
        return op2.Dat(self.node_set, self.value_shape, \
                self.field.val, valuetype, self.name)

    @cached_property
    def element_node_map(self):
        return op2.Map(self.field.mesh.element_set, self.node_set, self.field.mesh.shape.loc, \
                self.field.mesh.ndglno - 1, "%s_elem_node" % self.name)

    def temporary_dat(self, name):
        return op2.Dat(self.node_set, self.value_shape, \
                numpy.zeros(self.field.node_count), valuetype, name)

    def _reconstruct(self, element, count):
        # This code is class specific
        return FieldCoefficient(self._field, element, count)

class ScalarField(FieldCoefficient):

    def __init__(self,n,v,ft,op,uid,mesh=None):
        field = fluidity_state.ScalarField(n, v, ft, op, uid, mesh)
        super(ScalarField, self).__init__(field)

class VectorField(FieldCoefficient):

    def __init__(self,n,v,ft,op,dim,uid,mesh=None):
        field = fluidity_state.VectorField(n, v, ft, op, dim, uid, mesh)
        super(VectorField, self).__init__(field)

class TensorField(FieldCoefficient):

    def __init__(self,n,v,ft,op,dim0,dim1,uid,mesh=None):
        field = fluidity_state.TensorField(n, v, ft, op, dim0, dim1, uid, mesh)
        super(TensorField, self).__init__(field)

from solving import solve
