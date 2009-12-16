#!/usr/bin/env python

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

"""
Tools for dealing with gmsh mesh files
"""

import unittest

import fluidity.diagnostics.elements as elements
import fluidity.diagnostics.utils as utils

GMSH_UNKNOWN = None
GMSH_LINE = 1
GMSH_TRIANGLE = 2
GMSH_TETRAHEDRON = 4
GMSH_QUAD = 3
GMSH_HEXAHEDRON = 5
 
gmshElementTypeIds = ( \
    GMSH_UNKNOWN, \
    GMSH_LINE, GMSH_TRIANGLE, GMSH_QUAD, \
    GMSH_TETRAHEDRON, GMSH_HEXAHEDRON \
  )
 
class GmshElementType(elements.ElementType):
  """
  Class defining a Gmsh element type
  """

  _gmshElementTypeIdToElementTypeId = { \
      GMSH_UNKNOWN:elements.ELEMENT_UNKNOWN, \
      GMSH_LINE:elements.ELEMENT_LINE, \
      GMSH_TRIANGLE:elements.ELEMENT_TRIANGLE, GMSH_QUAD:elements.ELEMENT_QUAD, \
      GMSH_TETRAHEDRON:elements.ELEMENT_TETRAHEDRON, GMSH_HEXAHEDRON:elements.ELEMENT_HEXAHEDRON \
    }
  _elementTypeIdToGmshElementTypeId = utils.DictInverse(_gmshElementTypeIdToElementTypeId)
  
  def __init__(self, dim = None, nodeCount = None, gmshElementTypeId = None):
    if gmshElementTypeId is None:
      elements.ElementType.__init__(self, dim = dim, nodeCount = nodeCount)
    else:
      elements.ElementType.__init__(self, elementTypeId = self._gmshElementTypeIdToElementTypeId[gmshElementTypeId])
      
    self._UpdateGmshElementTypeId()
    self.RegisterEventHandler("elementTypeIdChange", self._UpdateGmshElementTypeId)
        
    return
    
  def _UpdateGmshElementTypeId(self):
    """
    Update the Gmsh type ID to reflect the element type ID
    """
    
    self._gmshElementTypeId = self._elementTypeIdToGmshElementTypeId[self._elementTypeId]
    
    return
    
  def GetGmshElementTypeId(self):
    return self._gmshElementTypeId
    
  def SetGmshElementTypeId(self, gmshElementTypeId):
    self.SetElementTypeId(self._gmshElementTypeIdToElementTypeId[gmshElementTypeId])
  
    return
    
class gmshtoolsUnittests(unittest.TestCase):    
  def testGmshElementType(self):
    type = GmshElementType(dim = 2, nodeCount = 4)
    self.assertEquals(type.GetGmshElementTypeId(), GMSH_QUAD)
    type.SetDim(3)
    self.assertEquals(type.GetGmshElementTypeId(), GMSH_TETRAHEDRON)
    type.SetNodeCount(8)
    self.assertEquals(type.GetGmshElementTypeId(), GMSH_HEXAHEDRON)
    type.SetGmshElementTypeId(GMSH_LINE)
    self.assertEquals(type.GetDim(), 1)
    self.assertEquals(type.GetNodeCount(), 2)
    self.assertRaises(KeyError, type.SetGmshElementTypeId, GMSH_UNKNOWN)
    self.assertRaises(AssertionError, type.SetDim, -1)
    self.assertRaises(AssertionError, type.SetNodeCount, -1)
    
    return
