# Copyright (C) 2009-2011 Anders Logg
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Garth N. Wells, 2012
#
# First added:  2008-10-22
# Last changed: 2012-08-18

__all__ = ["SubDomain", "DirichletBC", "PeriodicBC", "homogenize", "near"]

import types
import ufl
import sys
from pyop2.utils import as_tuple

# From Dolfin
EPS = 3.0e-16

def near(x, x0):
    return x0-EPS < x < x0+EPS

class SubDomain(object):
    
    def mark(self, mesh, id):
        coordinates = getattr(sys.modules['__main__'], 'coordinates')
        mesh.faces.boundaries[id] = []
        for node in mesh.faces.surface_node_list:
            # If node inside boundary then mark with that boundary id.
            # Need the coordinate field to do this.
            if self.inside(coordinates.val[node], 0):
                mesh.faces.boundaries[id].append(node)

    def inside(self, x, on_boundary):
        raise NotImplementedError("inside must be implemented.")

class BoundaryCondition(object):
    pass

class DirichletBC(BoundaryCondition):

    def __init__(self, *args, **kwargs):
        "Create Dirichlet boundary condition"

        # Copy constructor
        if len(args) == 1:
            if not isinstance(args[0], DirichletBC):
                log.error("bcs.py, create DirichletBC: " \
                          "Expecting a DirichletBC as only argument"\
                          " for copy constructor")

            other = args[0]
            self.mesh_arg = other.mesh_arg
            self.function_arg = other.function_arg
            self.domain_args = other.domain_args
            return

        # Only support setting BC value as float or tuple of floats initially.
        self.function_arg = as_tuple(args[1], float)
        self.domain_args = args[2:]

        # Add method argument if it's given
        if "method" in kwargs:
            args = tuple(list(args) + [kwargs["method"]])

# Creattion of Python class to avoid issue of SWIG directors going out
# of scope
class PeriodicBC(BoundaryCondition):

    def __init__(self, *args, **kwargs):
        "Create Periodic boundary condition"

        # Copy constructor
        if len(args) == 1:
            if not isinstance(args[0], cpp.PeriodicBC):
                cpp.dolfin_error("bcs.py",
                                 "create PeriodicBC",
                                 "Expecting a DirichleBC as only argument"\
                                 " for copy constructor")

            # Initialize base class
            cpp.PeriodicBC.__init__(self, args[0])

        elif len(args) == 2:
            if not isinstance(args[0], cpp.FunctionSpace):
                cpp.dolfin_error("bcs.py",
                                 "create PeriodicBC",
                                 "Expecting a FunctionSpace as first"\
                                 " constructor argument")

            if not isinstance(args[1], cpp.SubDomain):
                cpp.dolfin_error("bcs.py",
                                 "create PeriodicBC",
                                 "Expecting a SubDomain as second"\
                                 " constructor argument")

            # Store SubDomain to avoid scoping issue with SWIG directors
            self.domain_args = args[1:]

            # Initialize base class
            cpp.PeriodicBC.__init__(self, *args)

        else:
            cpp.dolfin_error("bcs.py",
                             "create PeriodicBC",
                             "Too many arguments passed to constructor")


def homogenize(bc):
    """
    Return a homogeneous version of the given boundary condition.

    *Arguments*
        bc
            a :py:class:`DirichletBC <dolfin.fem.bcs.DirichletBC>` instance,
            or a list/tuple of
            :py:class:`DirichletBC <dolfin.fem.bcs.DirichletBC>` instances.

    Other types of boundary conditions, like periodic, are ignored.

    If the given boundary condition is a list of boundary conditions,
    then a list of homogeneous boundary conditions is returned.

    """

    # Handle case when boundary condition is a list
    if isinstance(bc, (list, tuple)):
        bcs = bc
        return [homogenize(bc) for bc in bcs]

    # Only consider Dirichlet boundary conditions
    if not isinstance(bc, cpp.DirichletBC):
        return bc

    # Create zero function
    V = bc.function_space()
    if V.element().value_rank() == 0:
        zero = Constant(0)
    elif V.element().value_rank() == 1:
        zero = Constant([0]*V.element().value_dimension(0))
    else:
        cpp.dolfin_error("bcs.py",
                         "homogenize boundary condition",
                         "Unhandled value rank %d for homogenization of boundary conditions" % \
                             V.element().value_rank())

    # Create homogeneous boundary condition
    if len(bc.domain_args) == 1:
        new_bc = cpp.DirichletBC(V, zero, bc.domain_args[0])
    elif len(bc.domain_args) == 2:
        new_bc = cpp.DirichletBC(V, zero, bc.domain_args[0], bc.domain_args[1])
    else:
        cpp.dolfin_error("bcs.py",
                         "homogenize boundary condition",
                         "Unknown type of boundary specification")
    new_bc.domain_args = bc.domain_args

    return new_bc
