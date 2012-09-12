import ufl.argument
from ufl.assertions import ufl_assert
from ufl.finiteelement import FiniteElementBase

class Argument(ufl.argument.Argument):
    
    def __init__(self, element, mesh, count=None):
        super(Argument, self).__init__(element, count)
        self._mesh = mesh

    @property
    def mesh(self):
        return self._mesh

    def reconstruct(self, element=None, mesh=None, count=None):
        if mesh is None or mesh == self._mesh:
            mesh = self._mesh
        if element is None or element == self._element:
            element = self._element
        if count is None or count == self._count:
            count = self._count
        if count is self._count and element is self._element:
            return self
        ufl_assert(isinstance(element, FiniteElementBase),
                   "Expecting an element, not %s" % element)
        ufl_assert(isinstance(count, int),
                   "Expecting an int, not %s" % count)
        ufl_assert(element.value_shape() == self._element.value_shape(),
                   "Cannot reconstruct an Argument with a different value shape.")
        return Argument(element, mesh, count)

def TestFunction(field):
    return Argument(field.element(), field.mesh, -2)

def TrialFunction(field):
    return Argument(field.element(), field.mesh, -1)
