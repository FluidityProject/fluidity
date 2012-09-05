import numpy
from ufl import *
from pyop2 import ffc_interface, op2

import state_types as fluidity_state
from state_types import *

valuetype = numpy.float64

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
        return dimloc2cell[shape.dimension][shape.loc]
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

class FieldCoefficient(Coefficient):
    """Coefficient derived from a Fluidity field."""

    def __init__(self, field):
        field.rank = field2rank[field.description]
        self._field = field
        super(FieldCoefficient, self).__init__(ufl_element(field))

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

def solve(a, L):
    """Solve LSE a*x = L for x."""

    # FIXME do this only once
    op2.init(backend='sequential')

    # Generate code for mass and rhs assembly.

    a_code = ffc_interface.compile_form(a, "a")
    L_code = ffc_interface.compile_form(L, "L")

    a_kernel = op2.Kernel(a_code, "a_cell_integral_0_0" )
    L_kernel = op2.Kernel(L_code, "l_cell_integral_0_0")

    f = L.compute_form_data().coefficients[0]

    # Set up simulation data structures

    b = f.temporary_dat("b")
    x = f.temporary_dat("x")

    # Create sparsity and matrix

    sparsity = op2.Sparsity((f.element_node_map, f.element_node_map), \
            f.value_shape, "%s_sparsity" % f.name)
    mat = op2.Mat(sparsity, valuetype, "%s_mat" % f.name)

    # Assemble and solve
    
    elem_node = f.element_node_map
    op2.par_loop(a_kernel, f.element_set(3,3),
                 mat((elem_node[op2.i[0]], elem_node[op2.i[1]]), op2.INC),
                 a.measures()[0].domain_data().dat(elem_node, op2.READ))

    op2.par_loop(L_kernel, f.element_set,
                 b(elem_node, op2.INC),
                 L.measures()[0].domain_data().dat(elem_node, op2.READ),
                 f.dat(elem_node, op2.READ))

    op2.solve(mat, b, x)

    return x
