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
Field classes
"""

import math
import unittest

import fluidity.diagnostics.debug as debug
import fluidity.diagnostics.elements as elements
import fluidity.diagnostics.meshes as meshes

try:
  import numpy
except:
  debug.deprint("Warning: Failed to import numpy module")

import fluidity.diagnostics.optimise as optimise
import fluidity.diagnostics.utils as utils

class Field:
  """
  A generic field
  """

  def __init__(self, mesh, type = None, shape = None, data = None, name = None):  
    if type is None:
      assert(shape is None)
  
    self.SetName(name)
    
    self._mesh = mesh
    
    self._type = type
    self._shape = shape
      
    self._NewData()
    if not data is None:
      self.SetData(data)
      
    mesh.RegisterEventHandler("nodesAdded", self._OnNodesAdded)
  
    return
    
  def __add__(self, data):
    assert(len(data) == self.NodeCoordsCount())
  
    for i, datum in enumerate(data):
      self._data[i] += datum
    
    return self
    
  def __div__(self, factor):
    for i, datum in enumerate(self._data):
      self._data[i] /= factor
      
    return self
    
  def _DataLen(self):
    assert(not self._shape is None)
    if len(self._shape) == 0:
      return 0
      
    if not hasattr(self, "_dataLen"):
      self._dataLen = self._shape[0]
      for length in self._shape[1:]:
        self._dataLen *= length
    
    return self._dataLen
    
  def _NewData(self):
    if self._shape is None:
      self._data = [None for i in range(self.NodeCoordsCount())]
    else:
      self._data = []
      for i in range(self.NodeCoordsCount()):
        self._data.append(numpy.array([self._type() for i in range(self._DataLen())]))
        self._data[-1].shape = self._shape
    
    return
    
  def _OnNodesAdded(self):
    while len(self._data) < self.NodeCoordsCount():
      self._data.append(numpy.array([self._type() for i in range(self._DataLen())]))
      self._data[-1].shape = self._shape
      
    return
    
  def GetName(self):
    return self._name
    
  def SetName(self, name):
    self._name = name
    
    return
    
  def GetMesh(self):
    return self._mesh
    
  def NodeCount(self):
    return self._mesh.NodeCount()
    
  def NodeCoordsCount(self):
    return self.NodeCount()
    
  def GetType(self):
    return self._type
    
  def GetShape(self):
    return self._shape
    
  def GetVal(self, index):
    return self._data[index]
    
  def SetVal(self, index, val):
    if self._shape is None:
      self._data[index] = self._type(val)
    else:
      self._data[index] = numpy.array(val)
      self._data[index].shape = self._shape
    
    return
    
  def GetData(self):
    return self._data
    
  def SetData(self, data):
    assert(len(data) <= self.NodeCoordsCount())
    
    self._NewData()
    
    for i, datum in enumerate(data):
      self.SetVal(i, datum)
      
    return
    
  def NearbyVal(self, coord):
    """
    Return data value for a point in the mesh near to the supplied coordinate
    (not necessarily the nearest neighbour).
    """
  
    assert(self.NodeCoordsCount() > 0)
    assert(len(coord) == self._mesh.GetDim())
    
    nodeCoords, data = utils.KeyedSort(self._mesh.GetNodeCoords(), self._data, returnSortedKeys = True)
    
    if len(nodeCoords) == 1:
      return data[0]
    elif coord <= nodeCoords[0]:
      return data[0]
    elif coord >= nodeCoords[-1]:
      return data[-1]
      
    lBound = 0
    uBound = len(nodeCoords) - 1
    while True:
      index = (uBound + lBound) / 2
      
      if coord == nodeCoords[index]:
        break
      elif coord < nodeCoords[index]:
        uBound = index - 1
      else:
        lBound = index + 1
        
      if lBound >= uBound:
        index = min(lBound, uBound)
        break
        
    diff1 = 0.0
    diff2 = 0.0
    for i, val in enumerate(coord):
      diff1 += val - nodeCoords[index][i] ** 2.0
      diff2 += val - nodeCoords[index + 1][i] ** 2.0
    diff1 = math.sqrt(diff1)
    diff2 = math.sqrt(diff2)
    
    if diff2 >= diff1:
      return data[index]
    else:
      return data[index + 1]
      
  def IsTensorOfRank(self, rank):
    if not self._type is float:
      return False
    elif len(self._shape) < rank:
      return False
    
    for length in self._shape[rank:]:
      if not length == 1:
        return False
    
    return True
   
  def IsScalarField(self):    
    return self.IsTensorOfRank(0)
    
  def IsVectorField(self):
    return self.IsTensorOfRank(1)
    
  def IsTensorField(self):
    return self.IsTensorOfRank(2)
    
  def InsertInVtu(self, vtu):
    if self.IsScalarField():
      data = numpy.array(self._data)
      data.shape = (self.NodeCoordsCount(), 1)
      vtu.AddScalarField(self._name, data)
    elif self.IsVectorField():
      data = numpy.array(self._data)
      data.shape = tuple([self.NodeCoordsCount()] + list(self._shape))
      vtu.AddVectorField(self._name, data)
    elif self.IsTensorField():
      data = numpy.array(self._data)
      data.shape = tuple([self.NodeCoordsCount()] + list(self._shape))
      vtu.AddField(self._name, data)
    else:
      raise Exception("Invalid field type")
    
    return
    
  def ToVtu(self):
    vtu = self._mesh.ToVtu()
    self.InsertInVtu(vtu)
    
    return vtu
    
class fieldsUnittests(unittest.TestCase):
  def testField(self):
    # 2D triangle mesh
    mesh = meshes.Mesh(2)
    mesh.AddNodeCoord([0.0, 0.0])
    mesh.AddNodeCoord([1.0, 0.0])
    mesh.AddNodeCoord([0.0, 1.0])
    mesh.AddVolumeElement(elements.Element(nodes = [0, 1, 2]))
    mesh.AddSurfaceElement(elements.Element(nodes = [0, 1]))
    mesh.AddSurfaceElement(elements.Element(nodes = [0, 2]))
    mesh.AddSurfaceElement(elements.Element(nodes = [1, 2]))
    
    field = Field(mesh, name = "TestField", type = float, shape = [2])
    
    self.assertFalse(field.IsScalarField())
    self.assertTrue(field.IsVectorField())
    self.assertFalse(field.IsTensorField())
    
    self.assertEquals(field.NodeCount(), 3)
    for vec in field.GetData():
      self.assertEquals(len(vec), 2)
      for comp in vec:
        self.assertAlmostEquals(comp, 0.0)
        
    field.SetVal(2, [1.0, 2.0])
    vec = field.GetVal(2)
    self.assertEquals(len(vec), 2)
    self.assertAlmostEquals(vec[0], 1.0)
    self.assertAlmostEquals(vec[1], 2.0)
    
    if optimise.DebuggingEnabled():
      self.assertRaises(ValueError, field.SetVal, 2, [1.0, 2.0, 3.0])
    
    return
