#!/usr/bin/env python

import array
import exceptions
import os
import math
import fileinput
import re
import numpy
from xml.dom.minidom import parseString

try:
  import vtktools
except:
  pass

def parse_s(str):
    """Parse a .s file. Makes a dict vals that you use like:
    vals['maxp'][0] etc."""

    vars = []
    vals = {}
    fieldsize = 13 # how big is each number in the list?
    fields_per_line = 12 # how many variables printed out per line
    no_lines = 1 # how many lines correspond to each timestep
    if os.stat(str)[6] == 0:
      raise Exception, "Error: %s must not be empty!" % str

    f = open(str, "r")
    for line in f:
        if line.startswith("@@"):
            var_count = int(line[2:])
    f.close()

    # OK. Fix for bizarre behaviour when compiled with Sun compiler.
    # I'm not going near study.F to try to find what causes this:
    # ugly hack time!
    expected_line_pattern = []
    if var_count <= fields_per_line: expected_line_pattern = [var_count]
    else:
      tmpvarcount = var_count
      while tmpvarcount > fields_per_line:
        expected_line_pattern.append(fields_per_line)
        tmpvarcount = tmpvarcount - fields_per_line
      if tmpvarcount > 0: expected_line_pattern.append(tmpvarcount)

    # now let's see if it matches the expected pattern
    linecount = 0
    f = open(str, "r")
    actual_line_pattern = []
    for line in f:
      if not line.startswith("@"):
        linecount = linecount + 1
        actual_line_pattern.append(len(line[:-1]) / fieldsize)

    f.close()

    expected_line_pattern = expected_line_pattern * (linecount / len(expected_line_pattern))

    no_lines = int(math.ceil(var_count / float(fields_per_line)))
    lines = []

    f = open(str, "r")

    for line in f:
        if line.startswith("@("):
            numnames = line.split() # (NUM):NAME
            for numname in numnames:
                name = numname.strip().split(':')[-1].lower() # I love python
                vars.append(name)
                vals[name] = []
        elif not line.startswith("@"):
            lines.append(line)

    i = 0
    while i < len(lines):
        line = ""
        if actual_line_pattern[i] < expected_line_pattern[i]:
          print "Warning: .s file is not formatted as advertised. Skipping line %s." % i
          i = i + 1
          continue

        for j in range(i, i+no_lines):
            line += lines[j][:-1]

        line += '\n'
        for var in vars:
            k = vars.index(var)
            try:
              vals[var].append(float(line[fieldsize*k:fieldsize*(k+1)]))
            except ValueError:
              num = line[fieldsize*k:fieldsize*(k+1)]
              reverse_num = num[::-1]
              reverse_num = reverse_num.replace('-', '-e', 1).replace('+', '+e', 1)
              vals[var].append(float(reverse_num[::-1]))

        i += no_lines

    return vals

if __name__ == "__main__":
    import sys
    var = parse_s(sys.argv[1])
    print `var`

def compare_variables(reference, current, error, zerotol=1.0e-14):
    """This takes in an array for a particular variable
    (e.g. kinetic energy) containing the values of that
    variable at the timesteps. It compares the output
    of the current run against a reference run,
    checking that the relative error is within the bound
    specified."""

    assert len(reference) == len(current)
    relerrs = []

    for i in range(len(current)):
       if abs(reference[i]) > zerotol: # decide if reference[i] is "0.0" or not
           diff = abs(reference[i] - current[i])
           relerr = (diff / abs(reference[i]))
           relerrs.append(relerr)
       else:
           relerrs.append(abs(current[i])) # not really a relative error but however

    maxerr = max(relerrs)
    print "Asserting max relative error is smaller than", error
    print "max relative error: %s; index: %s" % (maxerr, relerrs.index(maxerr))
    assert maxerr < error

def compare_variable(reference, current, error, zerotol=1.0e-14):
    """Compares current value with reference value. Relative error
    should be smaller than 'error'. If the reference is within
    zerotol however, an absolute error will be used."""

    diff = abs(reference - current)
    if abs(reference) > zerotol: # decide if reference is "0.0" or not
       relerr = (diff / abs(reference))
       print "Asserting relative error is smaller than", error
       print "relative error: " + `relerr`
       assert relerr < error
    else:
       print "Asserting absolute error is smaller than", error
       print "absolute error: "+ `diff`
       assert diff < error


def tsunami_hit(fieldname, sfile, tol=0.00025):
    """This is for automation of tsunami modelling validation.
    Given a fieldname corresponding to the free surface height
    at a detector, it finds the first point of inflection
    in the values of that field (in time) then returns the acctim
    of that timestep."""

    s = parse_s(sfile)
    d = s[fieldname]

    for i in range(len(d)):
      if abs(d[i]) < tol: d[i] = 0.0

    oldval = d[0]

    # the parser works in 2 phases.
    # phase 1, mode == 0: the value is constant.
    # phase 2, mode == 1: the value has decreased (going_down = True)
    #                     or increased (going_down = False).
    #                     we stay in this mode until
    #                     the value increases (decreases);
    #                     this is the timestep we're looking for.

    mode = 0

    timestep = 0

    for i in range(len(d)):
        if mode == 0:
            if d[i] > oldval:
                mode = 1
                going_down = False
            if d[i] < oldval:
                mode = 1
                going_down = True
        if mode == 1:
            if going_down and d[i] > oldval: # we were going down, but are now gone up
                timestep = i
                break
            if not going_down and d[i] < oldval: # we were going up, but are now gone down
                timestep = i
                break

        oldval = d[i]

    return s['acctim'][timestep]


def getDistanceMeshDensity(file):
  v = vtktools.vtu(file)
  l = [0.0] * v.ugrid.GetNumberOfPoints()
  a = vtktools.arr(l)
  for i in range(v.ugrid.GetNumberOfPoints()):
    neighbours = v.GetPointPoints(i)
    sum = 0.0
    for neighbour in neighbours:
      sum = sum + v.GetDistance(i, neighbour)
    a[i] = sum / len(neighbours)
  return a

def getElementMeshDensity(file):
  v = vtktools.vtu(file)
  l = [0.0] * v.ugrid.GetNumberOfPoints()
  a = vtktools.arr(l)
  c = v.ugrid.GetCell(1)

  for i in range(v.ugrid.GetNumberOfPoints()):
    eles = v.GetPointCells(i)
    sum = 0.0
    for ele in eles:
      points = v.ugrid.GetCell(ele).GetPoints().GetData()
      sum = sum + c.ComputeVolume(points.GetTuple3(1), points.GetTuple3(2),
                                  points.GetTuple3(3), points.GetTuple3(4))
    a[i] = sum / len(eles)
  return a

class stat_parser(dict):
    """Parse a .stat file. The resulting mapping object is a heirarchy
of dictionaries. Most entries are of the form:

   parser[material_phase][field][statistic].

for example:

   p=stat_parser(filename)
   p['Material1']['Speed']['max']
"""

    def __init__(self, filename, subsample = 1):
    
      assert(subsample > 0)

      statfile=file(filename, "r")
      header_re=re.compile(r"</header>")
      xml="" # xml header.

      # extract the xml header stopping when </header> is reached.
      while 1:
        line=statfile.readline()
        if line=="":
          raise Exception("Unable to read .stat file header")
        xml=xml+line
        if re.search(header_re, line):
          break

      # now parse the xml.
      parsed=parseString(xml)
      
      binaryFormat = False
      constantEles = parsed.getElementsByTagName("constant")
      for ele in constantEles:
        name = ele.getAttribute("name")
        type = ele.getAttribute("type")
        value = ele.getAttribute("value")
        if name == "format":
          assert(type == "string")
          if value == "binary":
            binaryFormat = True
       
      nColumns = 0  
      for field in parsed.getElementsByTagName("field"):
        components = field.getAttribute("components")
        if components:
          nColumns += int(components)
        else:
          nColumns += 1

      if binaryFormat:           
        for ele in constantEles:
          name = ele.getAttribute("name")
          type = ele.getAttribute("type")
          value = ele.getAttribute("value")
          if name == "real_size":
            assert(type == "integer")            
            real_size = int(value)            
            if real_size == 4:
              realFormat = 'f'
            elif real_size == 8:
              realFormat = 'd'
            else:
              raise Exception("Unexpected real size: " + str(real_size))
          elif name == "integer_size":
            assert(type == "integer")            
            integer_size = int(value)            
            if not integer_size == 4:
              raise Exception("Unexpected integer size: " + str(real_size))
       
        nOutput = (os.path.getsize(filename + ".dat") / (nColumns * real_size)) / subsample
        
        columns = numpy.empty((nColumns, nOutput))
        statDatFile = file(filename + ".dat", "rb")   
        index = 0     
        while True:
          values = array.array(realFormat)
          try:
            values.fromfile(statDatFile, nColumns)
          except EOFError:
            break
            
          for i, value in enumerate(values):
            columns[i][index] = value
            
          index += 1
          if index > nOutput:
            # Ignore incomplete lines
            break
          if subsample > 1:
            # Ignore non-sampled lines
            statDatFile.seek(real_size * (subsample - 1) * nColumns, 1)        
        statDatFile.close()
        assert(index == nOutput)
      else:
        columns = [[] for i in range(nColumns)]
        lineNo = 0
        for line in statfile:
          entries = map(float, line.split())
             # Ignore incomplete lines        # Ignore non-sampled lines
          if len(entries) == len(columns) and lineNo == 0:
            map(list.append, columns, entries)
          lineNo = (lineNo + 1) % subsample
        columns = numpy.array(columns)
              
      for field in parsed.getElementsByTagName("field"):
        material_phase=field.getAttribute("material_phase")
        name=field.getAttribute("name")
        column=field.getAttribute("column")
        statistic=field.getAttribute("statistic")
        components=field.getAttribute("components")

        if material_phase:
          if not self.has_key(material_phase):
            self[material_phase]={}
          current_dict=self[material_phase]
        else:
          current_dict=self

        if not current_dict.has_key(name):
          current_dict[name]={}

        if components:
            column=int(column)
            components=int(components)
            current_dict[name][statistic]=columns[column-1:column-1+components]
        else:
            current_dict[name][statistic]=columns[int(column)-1]

def test_steady(vals, error, test_count = 1):
  """
  Test the test_count elements before the last element of vals against the
  last element of vals. If they are not within error of the last element of
  vals, raise an exception. Otherwise return. Can be used to test that a
  simulation has reached a steady state.
  """

  last_val = vals[len(vals) - 1]
  max_difference = 0.0
  index = None
  for i in range(len(vals) - test_count - 1, len(vals) - 1):
    difference = abs(vals[i] - last_val)
    if difference > max_difference:
      max_difference = difference
      index = i
  print "max difference: %s; index %s" % (max_difference, index)
  assert max_difference < error

  return

if __name__ == "__main__":
    import sys
    var = parse_s(sys.argv[1])
    print `var`

def shell():
  '''
  shell()

  Return ipython shell. To actually start the shell, invoke the function
  returned by this function.

  This is particularly useful for debugging embedded
  python or for crawling over the data when something has gone wrong.
  '''
  import sys
  
  if not hasattr(sys,"argv"):
    sys.argv=[]

  try:
    from IPython.Shell import IPShellEmbed
  except ImportError:
    sys.stderr.write(
      """
      *****************************************************
      *** Failed to import IPython. This probably means ***
      *** you don't have it installed. Please install   ***
      *** IPython and try again.                        ***
      *****************************************************
      """)
    raise

  banner = """
  This is an IPython shell embedded in Fluidity. You can use it to examine
  or even set variables. Press CTRL+d to exit and return to Fluidity.
  """

  ipshell = IPShellEmbed(banner=banner)

  return ipshell
  
