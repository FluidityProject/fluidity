#!/usr/bin/env python3

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

import array
import copy
import ctypes
import os
import tempfile
import unittest

import fluidity.diagnostics.bounds as bounds
import fluidity.diagnostics.calc as calc
import fluidity.diagnostics.debug as debug
import fluidity.diagnostics.elements as elements
import fluidity.diagnostics.filehandling as filehandling
import fluidity.diagnostics.meshes as meshes
import fluidity.diagnostics.mesh_halos as mesh_halos
import fluidity.diagnostics.utils as utils

GMSH_UNKNOWN = None
GMSH_LINE = 1
GMSH_TRIANGLE = 2
GMSH_TETRAHEDRON = 4
GMSH_QUAD = 3
GMSH_HEXAHEDRON = 5
GMSH_POINT = 15

gmshElementTypeIds = ( \
    GMSH_UNKNOWN, \
    GMSH_LINE, GMSH_TRIANGLE, GMSH_QUAD, \
    GMSH_TETRAHEDRON, GMSH_HEXAHEDRON, \
    GMSH_POINT \
  )

def FromGmshNodeOrder(nodes, type):
  """
  Permute Gmsh node ordering into default node ordering
  """

  newNodes = nodes

  if type.GetElementTypeId() == elements.ELEMENT_QUAD:
    newNodes = copy.deepcopy(nodes)
    newNodes[3] = nodes[2]
    newNodes[2] = nodes[3]

  return newNodes

def ToGmshNodeOrder(nodes, type):
  """
  Permute Gmsh node ordering into default node ordering
  """

  newNodes = nodes

  if type.GetElementTypeId() == elements.ELEMENT_QUAD:
    newNodes = copy.deepcopy(nodes)
    newNodes[2] = nodes[3]
    newNodes[3] = nodes[2]

  return newNodes

class GmshElementType(elements.ElementType):
  """
  Class defining a Gmsh element type
  """

  _gmshElementTypeIdToElementTypeId = { \
      GMSH_UNKNOWN:elements.ELEMENT_UNKNOWN, \
      GMSH_LINE:elements.ELEMENT_LINE, \
      GMSH_TRIANGLE:elements.ELEMENT_TRIANGLE, GMSH_QUAD:elements.ELEMENT_QUAD, \
      GMSH_TETRAHEDRON:elements.ELEMENT_TETRAHEDRON, GMSH_HEXAHEDRON:elements.ELEMENT_HEXAHEDRON, \
      GMSH_POINT:elements.ELEMENT_VERTEX
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

def ByteSwap(fileHandle):
  """
  Given a file handle, read a single 32-bit int, which should have the value 1.
  Determine whether this matches the endianness of the current architecture, or a swap is required.
  """

  iArr = array.array("i")
  iArr.fromfile(fileHandle, 1)
  if iArr[0] == 1:
    swap = False
  else:
    iArr.byteswap()
    if iArr[0] == 1:
      swap = True
    else:
      raise Exception("Invalid one byte")

  return swap


def ReadNonCommentLine(fileHandle):
  line = fileHandle.readline().decode("utf8")
  while len(line) > 0:
    line = line.strip()
    if len(line) > 0:
      return line
    line = fileHandle.readline().decode("utf8")

  return line


def ReadBinaryMshV2(fileHandle, dataSize):
  if dataSize == 4:
    realFormat = "f"
  elif dataSize == 8:
    realFormat = "d"
  else:
    raise Exception("Unrecognised real size " + str(dataSize))

  swap = ByteSwap(fileHandle)

  line = ReadNonCommentLine(fileHandle)
  assert(line == "$EndMeshFormat")

  # Read the Nodes section (no PhysicalNames section in binary)
  line = ReadNonCommentLine(fileHandle)
  assert(line == "$Nodes")

  # number-of-nodes
  # node-number x-coord y-coord z-coord
  # ...

  line = ReadNonCommentLine(fileHandle)
  nNodes = int(line)
  # Assume dense node IDs, but not necessarily ordered
  seenNode = [False for i in range(nNodes)]
  nodeIds = []
  nodes = []
  lbound = [calc.Inf() for i in range(3)]
  ubound = [-calc.Inf() for i in range(3)]
  for i in range(nNodes):
    iArr = array.array("i")
    rArr = array.array(realFormat)
    iArr.fromfile(fileHandle, 1)
    rArr.fromfile(fileHandle, 3)
    if swap:
      iArr.byteswap()
      rArr.byteswap()
    nodeId = iArr[0]
    coord = rArr
    assert(nodeId > 0)
    assert(not seenNode[nodeId - 1])
    seenNode[nodeId - 1] = True
    nodeIds.append(nodeId)
    nodes.append(coord)
    for j in range(3):
      lbound[j] = min(lbound[j], coord[j])
      ubound[j] = max(ubound[j], coord[j])

  line = ReadNonCommentLine(fileHandle)
  assert(line == "$EndNodes")

  nodes = utils.KeyedSort(nodeIds, nodes)
  bound = bounds.BoundingBox(lbound, ubound)
  indices = bound.UsedDimIndices()
  dim = len(indices)
  if dim < 3:
    nodes = [[coord[index] for index in indices] for coord in nodes]

  mesh = meshes.Mesh(dim)
  mesh.AddNodeCoords(nodes)

  # Read the Elements section
  line = ReadNonCommentLine(fileHandle)
  assert(line == "$Elements")

  # number-of-elements
  # element-header-binary
  # element-binary
  # ...

  # where element-header-binary is: elm-type num-elm num-tags
  # where element-binary is:
  #   num-elm * (4 + num-tags*4 + node-num*4)
  # node-num physical-tag elementary-tag node-nums ...

  line = ReadNonCommentLine(fileHandle)
  nEles = int(line)
  i = 0
  while i < nEles:
    iArr = array.array("i")
    iArr.fromfile(fileHandle, 3)
    if swap:
      iArr.byteswap()
    typeId = iArr[0]
    nSubEles = iArr[1]
    nIds = iArr[2]

    type = GmshElementType(gmshElementTypeId = typeId)

    for j in range(nSubEles):
      iArr = array.array("i")
      iArr.fromfile(fileHandle, 1 + nIds + type.GetNodeCount())
      if swap:
        iArr.byteswap()
      eleId = iArr[0]
      assert(eleId > 0)
      ids = iArr[1:1 + nIds]
      nodes = FromGmshNodeOrder(utils.OffsetList(iArr[-type.GetNodeCount():], -1), type)

      element = elements.Element(nodes, ids)

      if type.GetDim() == dim - 1:
        mesh.AddSurfaceElement(element)
      elif type.GetDim() == dim:
        mesh.AddVolumeElement(element)
      else:
        debug.deprint("Warning: Element of type " + str(type) + " encountered in " + str(dim) + " dimensions")

    i += nSubEles
  assert(i == nEles)

  line = ReadNonCommentLine(fileHandle)
  assert(line == "$EndElements")

  return mesh


def ReadBinaryMshV4(fileHandle, dataSize):
  if dataSize == 4:
    sizeFormat = "i"
  elif dataSize == 8:
    sizeFormat = "l"
  else:
    raise Exception("Unrecognised size_t size " + str(dataSize))

  swap = ByteSwap(fileHandle)

  line = ReadNonCommentLine(fileHandle)
  assert(line == "$EndMeshFormat")

  # skip ahead to Nodes section (possibly bypassing Entities)
  while line != "$Nodes":
    line = ReadNonCommentLine(fileHandle)
  assert(line == "$Nodes")

  # numEntityBlock(size_t) numNodes(size_t)
  #   minNodeTag(size_t) maxNodeTag(size_t)
  #
  # entityDim(int) entityTag(int) parametric(int) numNodes(size_t)
  #
  #  nodeTag(size_t) ...
  #  x(double) y(double) z(double) ...
  # ...
  sArr = array.array(sizeFormat)
  sArr.fromfile(fileHandle, 4)
  if swap:
    sArr.byteswap()
  numBlocks = sArr[0]
  numNodes = sArr[1]
  # assume dense nodes (we can check using the min/max id fields)
  seenNode = [False] * numNodes
  nodeIds = []
  nodes = []
  lbound = [calc.Inf() for i in range(3)]
  ubound = [-calc.Inf() for i in range(3)]
  for b in range(numBlocks):
    iArr = array.array("i")
    sArr = array.array(sizeFormat)
    iArr.fromfile(fileHandle, 3)
    sArr.fromfile(fileHandle, 1)
    if swap:
      iArr.byteswap()
      sArr.byteswap()

    subNodes = sArr[0]
    tagArr = array.array(sizeFormat)
    tagArr.fromfile(fileHandle, subNodes)
    if swap:
      tagArr.byteswap()

    for i in range(subNodes):
      rArr = array.array("d")
      rArr.fromfile(fileHandle, 3)
      if swap:
        rArr.byteswap()

      nodeId = tagArr[i]
      coord = rArr
      assert(nodeId > 0)
      assert(not seenNode[nodeId - 1])
      seenNode[nodeId - 1] = True
      nodeIds.append(nodeId)
      nodes.append(coord)
      for j in range(3):
        lbound[j] = min(lbound[j], coord[j])
        ubound[j] = max(ubound[j], coord[j])

  line = ReadNonCommentLine(fileHandle)
  assert(line == "$EndNodes")

  nodes = utils.KeyedSort(nodeIds, nodes)
  bound = bounds.BoundingBox(lbound, ubound)
  indices = bound.UsedDimIndices()
  dim = len(indices)
  if dim < 3:
    nodes = [[coord[index] for index in indices] for coord in nodes]

  mesh = meshes.Mesh(dim)
  mesh.AddNodeCoords(nodes)

  # read the Elements section
  line = ReadNonCommentLine(fileHandle)
  assert(line == "$Elements")

  # numEntityBlocks(size_t) numElements(size_t)
  #   minElementTag(size_t) maxElementTag(size_t)
  #
  # entityDim(int) entityTag(int) elementType(int) numElems(size_t)
  #   elementTag(size_t) nodeTag(size_t) ...
  #   ...
  # ...
  sArr = array.array(sizeFormat)
  sArr.fromfile(fileHandle, 4)
  if swap:
    sArr.byteswap()
  numEntities = sArr[0]
  numElems = sArr[1]
  for i in range(numEntities):
    iArr = array.array("i")
    sArr = array.array(sizeFormat)
    iArr.fromfile(fileHandle, 3)
    sArr.fromfile(fileHandle, 1)

    if swap:
      iArr.byteswap()
      sArr.byteswap()

    entityDim = iArr[0]
    entityTag = iArr[1]
    elementType = iArr[2]
    elemsInBlock = sArr[0]

    type = GmshElementType(gmshElementTypeId=elementType)

    for j in range(elemsInBlock):
      sArr = array.array(sizeFormat)
      sArr.fromfile(1 + type.GetNodeCount())
      if swap:
        sArr.byteswap()

      ids = [entityTag, sArr[0]]
      nodes = FromGmshNodeOrder(utils.OffsetList(sArr[1:], -1), type)
      element = elements.Element(nodes, ids)

      if type.GetDim() == dim - 1:
        mesh.AddSurfaceElement(element)
      elif type.GetDim() == dim:
        mesh.AddVolumeElement(element)
      else:
        debug.deprint("Warning: Element of type " + str(type) + " encountered in " + str(dim) + " dimensions")

    line = ReadNonCommentLine(fileHandle)
    assert(line == "$EndElements")

    return mesh

def ReadAsciiMshV2(fileHandle):
  line = ReadNonCommentLine(fileHandle)
  assert(line == "$EndMeshFormat")

  # Read the Nodes section
  line = ReadNonCommentLine(fileHandle)
  assert(line == "$Nodes")

  line = ReadNonCommentLine(fileHandle)
  nNodes = int(line)
  # Assume dense node IDs, but not necessarily ordered
  seenNode = [False for i in range(nNodes)]
  nodeIds = []
  nodes = []
  lbound = [calc.Inf() for i in range(3)]
  ubound = [-calc.Inf() for i in range(3)]
  for i in range(nNodes):
    line = ReadNonCommentLine(fileHandle)
    lineSplit = line.split()
    assert(len(lineSplit) == 4)
    nodeId = int(lineSplit[0])
    coord = [float(comp) for comp in lineSplit[1:]]
    assert(nodeId > 0)
    assert(not seenNode[nodeId - 1])
    seenNode[nodeId - 1] = True
    nodeIds.append(nodeId)
    nodes.append(coord)
    for j in range(3):
      lbound[j] = min(lbound[j], coord[j])
      ubound[j] = max(ubound[j], coord[j])

  line = ReadNonCommentLine(fileHandle)
  assert(line == "$EndNodes")

  nodes = utils.KeyedSort(nodeIds, nodes)
  bound = bounds.BoundingBox(lbound, ubound)
  indices = bound.UsedDimIndices()
  dim = len(indices)
  if dim < 3:
    nodes = [[coord[index] for index in indices] for coord in nodes]

  mesh = meshes.Mesh(dim)
  mesh.AddNodeCoords(nodes)

  # Read the Elements section

  line = ReadNonCommentLine(fileHandle)
  assert(line == "$Elements")

  line = ReadNonCommentLine(fileHandle)
  nEles = int(line)
  for i in range(nEles):
    line = ReadNonCommentLine(fileHandle)
    lineSplit = line.split()
    assert(len(lineSplit) > 3)
    eleId = int(lineSplit[0])
    assert(eleId > 0)
    typeId = int(lineSplit[1])
    nIds = int(lineSplit[2])

    type = GmshElementType(gmshElementTypeId = typeId)
    ids = [int(id) for id in lineSplit[3:3 + nIds]]
    nodes = FromGmshNodeOrder([int(node) - 1 for node in lineSplit[-type.GetNodeCount():]], type)
    element = elements.Element(nodes, ids)
    if type.GetDim() == dim - 1:
      mesh.AddSurfaceElement(element)
    elif type.GetDim() == dim:
      mesh.AddVolumeElement(element)
    else:
      debug.deprint("Warning: Element of type " + str(type) + " encountered in " + str(dim) + " dimensions")

  line = ReadNonCommentLine(fileHandle)
  assert(line == "$EndElements")

  return mesh


def ReadAsciiMshV4(fileHandle):
  line = ReadNonCommentLine(fileHandle)
  assert(line == "$EndMeshFormat")

  # skip to the Nodes section
  while line != "$Nodes":
    line = ReadNonCommentLine(fileHandle)
  assert(line == "$Nodes")

  line = ReadNonCommentLine(fileHandle)
  lineSplit = line.split()
  numBlocks = int(lineSplit[0])
  numNodes = int(lineSplit[1])
  seenNode = [False] * numNodes
  nodeIds = []
  nodes = []
  lbound = [calc.Inf() for i in range(3)]
  ubound = [-calc.Inf() for i in range(3)]
  for b in range(numBlocks):
    line = ReadNonCommentLine(fileHandle)
    lineSplit = line.split()
    subNodes = int(lineSplit[3])
    tagArr = [int(ReadNonCommentLine(fileHandle)) for i in range(subNodes)]
    for i in range(subNodes):
      line = ReadNonCommentLine(fileHandle)
      lineSplit = line.split()
      assert(len(lineSplit) == 3)

      nodeId = tagArr[i]
      coord = [float(comp) for comp in lineSplit]
      assert(nodeId > 0)
      assert(not seenNode[nodeId - 1])

      seenNode[nodeId - 1] = True
      nodeIds.append(nodeId)
      nodes.append(coord)

      for j in range(3):
        lbound[j] = min(lbound[j], coord[j])
        ubound[j] = max(ubound[j], coord[j])

  line = ReadNonCommentLine(fileHandle)
  assert(line == "$EndNodes")

  nodes = utils.KeyedSort(nodeIds, nodes)
  bound = bounds.BoundingBox(lbound, ubound)
  indices = bound.UsedDimIndices()
  dim = len(indices)
  if dim < 3:
    nodes = [[coord[index] for index in indices] for coord in nodes]

  mesh = meshes.Mesh(dim)
  mesh.AddNodeCoords(nodes)

  # read the Elements section
  line = ReadNonCommentLine(fileHandle)
  assert(line == "$Elements")

  line = ReadNonCommentLine(fileHandle)
  lineSplit = line.split()
  numEntities = int(lineSplit[0])
  numElems = int(lineSplit[1])
  for i in range(numEntities):
    line = ReadNonCommentLine(fileHandle)
    lineSplit = line.split()

    entityDim = int(lineSplit[0])
    entityTag = int(lineSplit[1])
    elementType = int(lineSplit[2])
    elemsInBlock = int(lineSplit[3])

    type = GmshElementType(gmshElementTypeId=elementType)

    for j in range(elemsInBlock):
      line = ReadNonCommentLine(fileHandle)
      lineSplit = line.split()
      assert(len(lineSplit) == 1 + type.GetNodeCount())

      ids = [entityTag, int(lineSplit[0])]
      nodes = FromGmshNodeOrder([int(node) - 1 for node in lineSplit[1:]], type)
      element = elements.Element(nodes, ids)

      if type.GetDim() == dim - 1:
        mesh.AddSurfaceElement(element)
      elif type.GetDim() == dim:
        mesh.AddVolumeElement(element)
      else:
        debug.deprint("Warning: Element of type " + str(type) + " encountered in " + str(dim) + " dimensions")

  line = ReadNonCommentLine(fileHandle)
  assert(line == "$EndElements")

  return mesh


def ReadMsh(filename):
  """
  Read a Gmsh msh file
  """

  fileHandle = open(filename, "rb")

  basename = filename.split(".")[0]
  hasHalo = filehandling.FileExists(basename + ".halo")

  # Read the MeshFormat section

  line = ReadNonCommentLine(fileHandle)
  assert(line == "$MeshFormat")

  line = ReadNonCommentLine(fileHandle)
  lineSplit = line.split()
  assert(len(lineSplit) == 3)
  version = lineSplit[0]
  fileType = int(lineSplit[1])
  dataSize = int(lineSplit[2])

  if version[0] == "4" and version < "4.1":
    raise Exception("gmshtools doesn't handle msh4 with minor version 0")

  if fileType == 1:
    # Binary format
    if version[0] == "4":
      mesh = ReadBinaryMshV4(fileHandle, dataSize)
    elif version[0] == "2":
      mesh = ReadBinaryMshV2(fileHandle, dataSize)
    else:
      raise Exception("Unknown gmsh major version")
  elif fileType == 0:
    # ASCII format
    if version[0] == "4":
      mesh = ReadAsciiMshV4(fileHandle)
    elif version[0] == "2":
      mesh = ReadAsciiMshV2(fileHandle)
    else:
      raise Exception("Unknown gmsh major version")
  else:
    raise Exception("File type " + str(fileType) + " not recognised")

  fileHandle.close()

  if hasHalo:
    # Read the .halo file
    debug.dprint("Reading .halo file")

    if mesh_halos.HaloIOSupport():
      halos = mesh_halos.ReadHalos(basename + ".halo")
      mesh.SetHalos(halos)
    else:
      debug.deprint("Warning: No .halo I/O support")

  return mesh

def WriteMsh(mesh, filename, binary = True):
  """
  Write a Gmsh msh file
  """

  if binary:
    # Binary format

    fileHandle = open(filename, "wb")

    # Write the MeshFormat section
    fileHandle.write("$MeshFormat\n".encode("utf8"))
    version = 2.1
    fileType = 1
    dataSize = ctypes.sizeof(ctypes.c_double)
    fileHandle.write(utils.FormLine([version, fileType, dataSize]).encode("utf8"))

    iArr = array.array("i", [1])
    iArr.tofile(fileHandle)
    fileHandle.write("\n".encode("utf8"))

    fileHandle.write("$EndMeshFormat\n".encode("utf8"))

    # Write the Nodes section

    fileHandle.write("$Nodes\n".encode("utf8"))
    fileHandle.write(utils.FormLine([mesh.NodeCoordsCount()]).encode("utf8"))

    for i, nodeCoord in enumerate(mesh.GetNodeCoords()):
      nodeCoord = list(nodeCoord)
      while len(nodeCoord) < 3:
        nodeCoord.append(0.0)
      assert(len(nodeCoord) == 3)

      iArr = array.array("i", [i + 1])
      rArr = array.array("d", nodeCoord)
      iArr.tofile(fileHandle)
      rArr.tofile(fileHandle)
    fileHandle.write("\n".encode("utf8"))

    fileHandle.write("$EndNodes\n".encode("utf8"))

    # Write the Elements section

    fileHandle.write("$Elements\n".encode("utf8"))
    fileHandle.write(utils.FormLine([mesh.SurfaceElementCount() + mesh.VolumeElementCount()]).encode("utf8"))

    eleSort = {}
    for ele in mesh.GetSurfaceElements() + mesh.GetVolumeElements():
      eleType = ele.GetType()
      gmshType = GmshElementType(dim = eleType.GetDim(), nodeCount = eleType.GetNodeCount())

      key = (gmshType.GetGmshElementTypeId(), len(ele.GetIds()))
      if key in eleSort:
        eleSort[key].append(ele)
      else:
        eleSort[key] = [ele]

    index = 1
    for gmshEleId, nIds in eleSort:
      eles = eleSort[(gmshEleId, nIds)]
      iArr = array.array("i", [gmshEleId, len(eles), nIds])
      iArr.tofile(fileHandle)
      for ele in eles:
        iArr = array.array("i", [index] + list(ele.GetIds()) + utils.OffsetList(ToGmshNodeOrder(ele.GetNodes(), ele.GetType()), 1))
        iArr.tofile(fileHandle)
        index += 1
    assert(index == mesh.SurfaceElementCount() + mesh.VolumeElementCount() + 1)
    fileHandle.write("\n".encode("utf8"))

    fileHandle.write("$EndElements\n".encode("utf8"))
  else:
    # ASCII format

    fileHandle = open(filename, "w")

    # Write the MeshFormat section
    fileHandle.write("$MeshFormat\n")
    version = 2.1
    fileType = 0
    dataSize = ctypes.sizeof(ctypes.c_double)
    fileHandle.write(utils.FormLine([version, fileType, dataSize]))
    fileHandle.write("$EndMeshFormat\n")

    # Write the Nodes section

    fileHandle.write("$Nodes\n")
    fileHandle.write(utils.FormLine([mesh.NodeCoordsCount()]))
    for i, nodeCoord in enumerate(mesh.GetNodeCoords()):
      nodeCoord = list(nodeCoord)
      while len(nodeCoord) < 3:
        nodeCoord.append(0.0)
      assert(len(nodeCoord) == 3)
      fileHandle.write(utils.FormLine([i + 1, nodeCoord]))
    fileHandle.write("$EndNodes\n")

    # Write the Elements section

    fileHandle.write("$Elements\n")
    fileHandle.write(utils.FormLine([mesh.SurfaceElementCount() + mesh.VolumeElementCount()]))
    for i, ele in enumerate(mesh.GetSurfaceElements() + mesh.GetVolumeElements()):
      eleType = ele.GetType()
      gmshType = GmshElementType(dim = eleType.GetDim(), nodeCount = eleType.GetNodeCount())
      ids = ele.GetIds()
      fileHandle.write(utils.FormLine([i + 1, gmshType.GetGmshElementTypeId(), len(ids), ids, utils.OffsetList(ToGmshNodeOrder(ele.GetNodes(), eleType), 1)]))
    fileHandle.write("$EndElements\n")

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

  def testMshIo(self):
    tempDir = tempfile.mkdtemp()
    oldMesh = meshes.Mesh(2)
    oldMesh.AddNodeCoord([0.0, 0.0])
    oldMesh.AddNodeCoord([1.0, 0.0])
    oldMesh.AddNodeCoord([0.0, 1.0])
    oldMesh.AddVolumeElement(elements.Element(nodes = [0, 1, 2]))
    oldMesh.AddSurfaceElement(elements.Element(nodes = [0, 1]))
    oldMesh.AddSurfaceElement(elements.Element(nodes = [0, 2]))
    oldMesh.AddSurfaceElement(elements.Element(nodes = [1, 2]))
    filename = os.path.join(tempDir, "temp")
    WriteMsh(oldMesh, filename, binary = False)
    newMesh = ReadMsh(filename)
    filehandling.Rmdir(tempDir, force = True)
    self.assertEquals(oldMesh.GetDim(), newMesh.GetDim())
    self.assertEquals(oldMesh.NodeCount(), newMesh.NodeCount())
    self.assertEquals(oldMesh.SurfaceElementCount(), newMesh.SurfaceElementCount())
    self.assertEquals(oldMesh.VolumeElementCount(), newMesh.VolumeElementCount())

    tempDir = tempfile.mkdtemp()
    oldMesh = meshes.Mesh(2)
    oldMesh.AddNodeCoord([0.0, 0.0])
    oldMesh.AddNodeCoord([1.0, 0.0])
    oldMesh.AddNodeCoord([0.0, 1.0])
    oldMesh.AddVolumeElement(elements.Element(nodes = [0, 1, 2]))
    oldMesh.AddSurfaceElement(elements.Element(nodes = [0, 1]))
    oldMesh.AddSurfaceElement(elements.Element(nodes = [0, 2]))
    oldMesh.AddSurfaceElement(elements.Element(nodes = [1, 2]))
    filename = os.path.join(tempDir, "temp")
    WriteMsh(oldMesh, filename, binary = True)
    newMesh = ReadMsh(filename)
    filehandling.Rmdir(tempDir, force = True)
    self.assertEquals(oldMesh.GetDim(), newMesh.GetDim())
    self.assertEquals(oldMesh.NodeCount(), newMesh.NodeCount())
    self.assertEquals(oldMesh.SurfaceElementCount(), newMesh.SurfaceElementCount())
    self.assertEquals(oldMesh.VolumeElementCount(), newMesh.VolumeElementCount())

    return
