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
Fluidity related tools
"""

import glob
import os
import subprocess
import tempfile
import unittest

import fluidity.diagnostics.debug as debug

try:
  from fluidity_tools import *
except:
  debug.deprint("Warning: Failed to import fluidity_tools module")
  
import fluidity.diagnostics.filehandling as filehandling
import fluidity.diagnostics.utils as utils
  
def FluidityBinary(binaries = ["dfluidity-debug", "dfluidity"]):
  """
  Return the command used to call Fluidity
  """

  binary = None
  
  for possibleBinary in binaries:
    process = subprocess.Popen(["which", possibleBinary], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    process.wait()
    if process.returncode == 0:
      binary = possibleBinary
      break
  
  if binary is None:
    raise Exception("Failed to find Fluidity binary")
  
  return binary
  
class Stat:
  """
  Class for handling .stat files. Similiar to the dictionary returned by
  stat_parser, but with some annoying features fixed.
  """

  def __init__(self, filename = None, delimiter = "%"):
    self.SetDelimiter(delimiter)
    self._s = {}
    if not filename is None:
      self.Read(filename)
      
    return
    
  def __getitem__(self, key):
    """
    Index into the .stat with the given key (or path)
    """
  
    def SItem(s, key, delimiter):
      if key in s:
        return s[key]
      else:
        keySplit = key.split(delimiter)
        for i in range(len(keySplit)):
          key = utils.FormLine(keySplit[:i], delimiter = delimiter, newline = False)
          if key in s:
            if isinstance(s[key], dict):
              try:
                # Tolerate a failure when recursing, as the key may have been
                # eroneously split
                return SItem(s[key], utils.FormLine(keySplit[i:], delimiter = delimiter, newline = False), delimiter)
              except Exception:
                pass
            else:
              return s[key]

        raise Exception("Key not found")
    
    item = SItem(self._s, key, self._delimiter)
    if isinstance(item, dict):
      subS = Stat(delimiter = self._delimiter)
      subS._s = item
      return subS
    else:
      return item
    
  def __str__(self):
    paths = self.Paths()
    string = "Stat file:\n"
    for i, path in enumerate(paths):
      string += path
      if i < len(paths):
        string += "\n"
        
    return string
    
  def _PathSplit(self, path):
    """
    Return a list of keys into the stat dictionary for the supplied path
    """
    
    def SPathSplit(s, delimiter, path):
      pathSplit = path.split(delimiter)
      index = 0
      newPath = pathSplit[index]
      while not newPath in s.keys():
        index += 1
        newPath += delimiter + pathSplit[index]
        
      paths = []
      paths.append(newPath)
      if isinstance(s[newPath], dict):
        paths += SPathSplit(s[newPath], delimiter, path[len(newPath) + len(delimiter):])
        
      return paths
  
    return SPathSplit(self._s, self._delimiter, path)
    
  def keys(self):
    return self.Paths()
    
  def GetDelimiter(self):
    return self._delimiter 
   
  def SetDelimiter(self, delimiter):
    self._delimiter = delimiter
    
    return
    
  def Paths(self):
    """
    Return all valid paths
    """
  
    def SPaths(s, delimiter, base = ""):
      if len(base) > 0:
        base += delimiter
    
      paths = []
      for key in s.keys():
        if isinstance(s[key], dict):
          paths += SPaths(s[key], delimiter, base = base + key)
        else:
          paths.append(base + key)
          
      return paths
    
    return SPaths(self._s, self._delimiter)
    
  def PathLists(self):
    """
    Return all valid paths as a series of key lists
    """
    
    def SPathLists(s, delimiter, base = []):    
      paths = []
      for key in s.keys():
        if isinstance(s[key], dict):
          paths += SPathLists(s[key], delimiter, base = base + [key])
        else:
          paths.append(base + [key])
          
      return paths
    
    return SPathLists(self._s, self._delimiter)
    
  def HasPath(self, path):
    """
    Return whether the supplied path is valid for this Stat
    """
    
    try:
      self[path]
      return True
    except Exception:
      return False
    
  def FormPath(self, *args):
    path = ""
    for i, arg in enumerate(args):
      path += arg
      if i < len(args) - 1:
        path += self._delimiter
    
    return path
    
  def FormPathFromList(self, pathList):
    path = ""
    for i, entry in enumerate(pathList):
      path += entry
      if i < len(pathList) - 1:
        path += self._delimiter
        
    return path
    
  def SplitPath(self, path):
    return path.split(self._delimiter)
    
  def Read(self, filename):
    """
    Read a .stat file
    """
    
    def ParseRawS(s, delimiter):    
      newS = {}
      for key1 in s.keys():
        assert(not key1 in ["val", "value"])
        if isinstance(s[key1], dict):
          if len(s[key1].keys()) == 1 and s[key1].keys()[0] in ["val", "value"]:
            newS[str(key1)] = [float(val) for val in s[key1][s[key1].keys()[0]]]
          else:
            subS = ParseRawS(s[key1], delimiter)
            newS[str(key1)] = {}
            for key2 in subS.keys():
              newS[str(key1)][str(key2)] = subS[key2]
        else:        
          rank = len(s[key1].shape)
          if rank > 1:
            assert(rank == 2)
            for i in range(len(s[key1])):
              # I'm not very happy with this hack - this should be fixed
              # properly at some point
              newS[str(key1) + delimiter + str(i + 1)] = [val for val in s[key1][i]]
          else:
            try:
              newS[str(key1)] = [float(val) for val in s[key1]]
            except TypeError:
              debug.deprint("Type error for data " + str(s[key1]), 0)
              raise Exception("ParseRawS failure")
            except ValueError:
              debug.deprint("Value error for data " + str(s[key1]), 0)
              raise Exception("ParseRawS failure")
          
      return newS
      
    self._s = ParseRawS(stat_parser(filename), self._delimiter)
    
    return
    
  def Write(self, filename):
    """
    Write to a .stat file
    """
    
    def FieldTag(column, name, statistic, materialPhase = None, components = None):
      tagList = []
      tagList.append("<field")
      tagList.append("column=\"" + str(column) + "\"")
      tagList.append("name=\"" + name + "\"")
      tagList.append("statistic=\"" + statistic + "\"")
      if not materialPhase is None:
        tagList.append("material_phase=\"" + materialPhase + "\"")
      if not components is None or components == 1:
        tagList.append("components=" + str(components) + "\"")
      tagList.append("/>")
      
      return utils.FormLine(tagList)
      
    def EntryComponents(path):
      try:
        components = len(self[path])
      except TypeError:
        components = 1
      
      return
    
    materialPhaseTags = []
    materialPhaseData = []
    otherTags = []
    otherData = []
    for path in self.Paths():
      pathSplit = self._PathSplit(path)
      
      if len(pathSplit) == 1:
        otherTags.append(pathSplit + [EntryComponents(path)])
        otherData.append(self[path])
      elif len(pathSplit) == 2:
        otherTags.append(pathSplit + [EntryComponents(path)])
        otherData.append(self[path])
      elif len(pathSplit) == 3:
        materialPhaseTags.append(pathSplit + [EntryComponents(path)])
        materialPhaseData.append(self[path])
      else:
        assert(len(pathSplit) > 0)
        otherTags.append(pathSplit[:2] + [utils.FormLine(pathSplit[2:], delimiter = self._delimiter, newline = False)] + [EntryComponents(path)])
        otherData.append(self[path])
        
    statHandle = open(filename, "w")
    statHandle.write("<header>\n")
     
    column = 1
    for tag in otherTags:
      assert(len(tag) == 2)
      statHandle.write(FieldTag(column, tag[0], tag[1], components = tag[2]))
      column += 1
    for tag in materialPhaseTags:
      assert(len(tag) == 3)
      statHandle.write(FieldTag(column, tag[1], tag[2], materialPhase = tag[0], components = tag[3]))
      column += 1
     
    statHandle.write("</header>\n")
    
    entries = 0
    if len(materialPhaseData) > 0:
      entries = len(materialPhaseData[0])
    elif len(otherData) > 0:
      entries = len(otherData[0])
      
    for i in range(entries):
      for data in otherData:
        statHandle.write(str(data[i]) + " ")
      for data in materialPhaseData:
        statHandle.write(str(data[i]) + " ")
      statHandle.write("\n")
     
    statHandle.flush()
    statHandle.close()
      
    return
    
def DetectorArrays(stat):
  """
  Return a dictionary of detector array lists contained in the supplied stat
  """
  
  arrayNames = []
  arrayIndices = []
  notArrayNames = []
  
  # First, find the names of the arrays
  for path in stat.PathLists():  
    if isinstance(stat[stat.FormPathFromList(path)], Stat):
      # This isn't a leaf node
      continue
      
    finalKey = path[-1]
    keySplit = finalKey.split("_")
    if len(keySplit) == 0 or not len(stat.SplitPath(keySplit[-1])) in [1, 2] or not utils.IsIntString(stat.SplitPath(keySplit[-1])[0]):
      # This definitely can't be an entry in a detector array
      continue
    
    arrayName = utils.FormLine(path[:-1] + [utils.FormLine(keySplit[:-1], delimiter = "_", newline = False)], delimiter = "%", newline = False)
    if len(stat.SplitPath(keySplit[-1])) > 1:
      arrayName = stat.FormPath(arrayName, stat.SplitPath(keySplit[-1])[1])
      index = stat.SplitPath(keySplit[-1])[0]
    else:
      index = keySplit[-1]
      
    if arrayName in notArrayNames:
      # We've already discovered that this isn't an entry in a detector array
      continue
          
    if arrayName in arrayNames:
      arrayIndex = arrayNames.index(arrayName)
    
    if index <= 0:
      # This isn't a valid index
      notArrayNames.append(arrayName)
      if arrayName in arrayNames:
        arrayIndices.remove(arrayIndices[arrayIndex])
        arrayNames.remove(arrayName)
      continue
    if arrayName in arrayNames and index in arrayIndices[arrayIndex]:
      # We've already seen this index
      notArrayNames.append(arrayName)
      arrayIndices.remove(arrayIndices[arrayIndex])
      arrayNames.remove(arrayName)
      continue
        
    if arrayName in arrayNames:
      arrayIndices[arrayIndex].append(index)
    else: 
      arrayNames.append(arrayName)
      arrayIndices.append([index])

  assert(len(arrayNames) == len(arrayIndices))
  for i, name in enumerate(arrayNames):
    intIndices = [int(index) for index in arrayIndices[i]]
    intIndices.sort()
    for i, index in enumerate(intIndices):
      if not i + 1 == index:
        # The indices are not consecutive from one
        arrayNames.remove(name)
        arrayIndices.remove(arrayIndices[i])
        break
  
  # OK, we have a list of valid detector array names. Let's form the dictionary
  # of data.
  
  data = {}
  for i, name in enumerate(arrayNames):
    indices = utils.KeyedSort([int(index) for index in arrayIndices[i]], arrayIndices[i])
    
    data[name] = []
    for index in indices:
      if stat.HasPath(name + "_" + index):
        data[name].append(stat[name + "_" + index])
      else:
        nameSplit = stat.SplitPath(name)
        data[name].append(stat[stat.FormPath(stat.FormPathFromList(nameSplit[:-1]) + "_" + index, nameSplit[-1])])
    
  return data
  
def VtuFilenames(project, firstId, lastId = None, extension = ".vtu"):
  """
  Return vtu filenames for a Fluidity simulation, in the supplied range of IDs
  """
  
  if lastId is None:
    lastId = firstId
  assert(lastId >= firstId)
  
  filenames = []
  for id in range(firstId, lastId + 1):
    filenames.append(project + "_" + str(id) + extension)
  
  return filenames

def PVtuFilenames(project, firstId, lastId = None, extension = ".pvtu"):
  """
  Return pvtu filenames for a Fluidity simulation, in the supplied range of IDs
  """
  
  return VtuFilenames(project, firstId, lastId = lastId, extension = extension)
  
def FindVtuFilenames(project, firstId, lastId = None, extensions = [".vtu", ".pvtu"]):
  """
  Find vtu filenames for a Fluidity simulation, in the supplied range of IDs
  """
  
  if lastId is None:
    lastId = firstId
  assert(lastId >= firstId)
        
  filenames = []
  for id in range(firstId, lastId + 1):
    filename = project + "_" + str(id)
    for i, ext in enumerate(extensions):
      try:
        os.stat(filename + ext)
        filename += ext
        break
      except OSError:
        pass
      if i == len(extensions) - 1:
        raise Exception("Failed to find input file with ID " + str(id))
    filenames.append(filename)  
  
  return filenames

def FindPvtuVtuFilenames(project, firstId, lastId = None, pvtuExtension = ".pvtu", vtuExtension = ".vtu"):
  """
  Find vtu filenames for a Fluidity simulation, in the supplied range of IDs,
  retuning the vtus for each given pvtu filename
  """
  
  pvtuFilenames = FindVtuFilenames(project, firstId, lastId = lastId, extensions = [pvtuExtension])
  
  filenames = []
  for pvtuFilename in pvtuFilenames:
    i = 0
    while True:
      filename = pvtuFilename[:-len(pvtuExtension)] + "_" + str(i) + vtuExtension
      try:
        os.stat(filename)
        filenames.append(filename)
        i += 1
      except OSError:
        break
    if i == 0:
      raise Exception("Failed to find vtus for pvtu " + pvtuFilename)
    
  return filenames
  
def FindMaxVtuId(project, extensions = [".vtu", ".pvtu"]):
  """
  Find the maximum Fluidity vtu ID for the supplied project name
  """
  
  filenames = []
  for ext in extensions:
    filenames += glob.glob(project + "_?*" + ext)
  
  maxId = None
  for filename in filenames:
    id = filename[len(project) + 1:-len(filename.split(".")[-1]) - 1]
    try:
      id = int(id)
    except ValueError:
      continue
      
    if maxId is None:
      maxId = id
    else:
      maxId = max(maxId, id)
  
  return maxId
  
def FindFinalVtu(project, extensions = [".vtu", ".pvtu"]):
  """
  Final the final Fluidity vtu for the supplied project name
  """
  
  return FindVtuFilenames(project, FindMaxVtuId(project, extensions = extensions), extensions = extensions)[0]
  
def FindPFilenames(basename, extension):
  filenames = []
  i = 0
  while True:
    filename = basename + "_" + str(i)
    if filehandling.FileExists(filename + extension):
      filenames.append(filename)
    else:
      break
    i += 1
  
  return filenames

class fluiditytoolsUnittests(unittest.TestCase):
  def testFluidityToolsSupport(self):
    import fluidity_tools
    
    return
    
  def testStatIO(self):
    tempDir = tempfile.mkdtemp()
    filename = os.path.join(tempDir, "temp.stat")
    stat = Stat()
    stat.Write(filename)
    stat = Stat(filename)
    filehandling.Rmdir(tempDir, force = True)
    
    return
  
  def testFindVtuFilenames(self):
    tempDir = tempfile.mkdtemp()
    project = os.path.join(tempDir, "project")
    for i in range(2, 4):
      filehandling.Touch(project + "_" + str(i) + ".vtu")
    filenames = FindVtuFilenames(project, firstId = 2, lastId = 3)
    self.assertEquals(len(filenames), 2)
    filehandling.Touch(project + "_4.pvtu")
    filenames = FindVtuFilenames(project, firstId = 2, lastId = 4)
    self.assertEquals(len(filenames), 3)
    filehandling.Rmdir(tempDir, force = True)
    
    return
  
  def testFindPvtuVtuFilenames(self):
    tempDir = tempfile.mkdtemp()
    project = os.path.join(tempDir, "project")
    for i in range(2, 4):
      for j in range(3):
        filehandling.Touch(project + "_" + str(i) + "_" + str(j) + ".vtu")
      filehandling.Touch(project + "_" + str(i) + ".pvtu")
    filenames = FindPvtuVtuFilenames(project, firstId = 2, lastId = 3)
    self.assertEquals(len(filenames), 6)
    filehandling.Rmdir(tempDir, force = True)
    
    return
