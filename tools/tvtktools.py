from enthought.tvtk.api import tvtk
import math
import numpy
import scipy.linalg

# All returned arrays are cast into either numpy or numarray arrays
arr=numpy.array

class vtu:
  """Unstructured grid object to deal with VTK unstructured grids."""
  def __init__(self, filename):
    """Creates a vtu object by reading the specified file."""
    gridreader=tvtk.XMLUnstructuredGridReader(file_name=filename)
    gridreader.update()
    self.ugrid=gridreader.output
    self.filename=filename
    
  def GetScalarField(self, name):
    """Returns an array with the values of the specified scalar field."""
    return arr(self.ugrid.point_data.get_array(name)) 
 
  def GetVectorField(self, name):
    """Returns an array with the values of the specified vector field."""
    return arr(self.ugrid.point_data.get_array(name)) 

  def GetVectorNorm(self, name):
    """Return the field with the norm of the specified vector field."""
    v = self.GetVectorField(name)
    n = []
    norm = scipy.linalg.norm
    for node in range(self.ugrid.number_of_points):
      n.append(norm(v[node]))

    return arr(n)
  
  def GetField(self,name):
    """Returns an array with the values of the specified field."""
    pointdata=self.ugrid.point_data
    vtkdata=pointdata.get_array(name)
    nc=vtkdata.number_of_components
    nt=vtkdata.number_of_tuples
    array=arr(vtkdata)
    if nc==9:
      return array.reshape(nt,3,3)
    elif nc==4:
      return array.reshape(nt,2,2)
    else:
      return array.reshape(nt,nc)
    
  def Write(self, filename=[]):
    """Writes the grid to a vtu file.
    
    If no filename is specified it will use the name of the file originally 
    read in, thus overwriting it!
    """
    if filename==[]: 
      filename=self.filename
    gridwriter=tvtk.XMLUnstructuredGridWriter(file_name=filename, input=self.ugrid)
    gridwriter.write()
    
  def AddScalarField(self, name, array):
    """Adds a scalar field with the specified name using the values from the array."""
    # In vtktools.py the following used SetNumberOfValues=len(array)
    data = tvtk.FloatArray(number_of_tuples=len(array), name=name)
    for i in range(len(array)):
      data.set_value(i, array[i])

    pointdata=self.ugrid.point_data
    pointdata.add_array(data)
    pointdata.set_active_scalars(name)

  def AddVectorField(self, name, array):
    """Adds a vector field with the specified name using the values from the array."""
    n=array.size
    # In vtktools.py the following used SetNumberOfValues=n
    data = tvtk.FloatArray(number_of_components=array.shape[1], number_of_tuples=n, name=name)
    for i in range(n):
      data.set_value(i, array.reshape(n)[i])

    pointdata=self.ugrid.point_data
    pointdata.add_array(data)
    pointdata.set_active_vectors(name)

  def AddField(self, name, array):
    """Adds a field with arbitrary number of components under the specified name using."""
    n=array.size
    sh=arr(array.shape)
    # number of tuples is sh[0]
    # number of components is the product of the rest of sh
    data = vtk.vtkFloatArray(number_of_components=sh[1:].prod(), number_of_tuples=n, name=name)
    flatarray=array.reshape(n)
    for i in range(n):
      data.set_value(i, flatarray[i])

    pointdata=self.ugrid.point_data
    pointdata.add_array(data)

  def ApplyProjection(self, projection_x, projection_y, projection_z):
    """Applys a projection to the grid coordinates. This overwrites the existing values."""
    npoints = self.ugrid.number_of_points
    for i in range (npoints):
      (x,y,z) = self.ugrid.get_point(i)
      new_x = eval (projection_x)
      new_y = eval (projection_y)
      new_z = eval (projection_z)
      self.ugrid.points.set_point(i, new_x, new_y, new_z)
                
  def ProbeData(self, coordinates, name):
    """Interpolate field values at these coordinates."""

    # Initialise locator
    bbox = self.ugrid.bounds
    locator = tvtk.PointLocator(data_set=self.ugrid, tolerance=10.0)
    locator.update()

    # Initialise probe
    points = tvtk.Points()
    ilen, jlen = coordinates.shape
    for i in range(ilen):
      points.insert_next_point(coordinates[i][0], coordinates[i][1], coordinates[i][2])
    polydata = tvtk.PolyData(points=points)
    probe = tvtk.ProbeFilter(input=polydata, source=self.ugrid)
    probe.update()

    # Reposition invalid nodes at nearest mesh vertices
    alid_ids = probe.valid_points
    
    valid_points = tvtk.Points()
    valid_loc = 0
    for i in range(ilen):
      if valid_ids.get_tuple1(valid_loc) == i:
        valid_points.insert_next_point(coordinates[i][0], coordinates[i][1], coordinates[i][2])
        valid_loc = valid_loc + 1
      else:
        nearest = locator.find_closest_point([coordinates[i][0], coordinates[i][1], coordinates[i][2]])
        point = self.ugrid.points.get_point(nearest)
        valid_points.insert_next_point(point[0], point[1], point[2])
    polydata.points=valid_points
    probe.input=polydata
    probe.update()
        
    # Get final updated values
    pointdata=probe.output.point_data
    vtkdata=pointdata.get_array(name)
    nc=vtkdata.number_of_components()
    nt=vtkdata.number_of_tuples()
    array = arr(vtkdata)
    array.shape = (nt,nc)
    return array
    
  def RemoveField(self, name):
    """Removes said field from the unstructured grid."""
    self.ugrid.point_data.remove_array(name)

  def GetLocations(self):
    """Returns an array with the locations of the nodes."""
    return arr(self.ugrid.points.data)
    
  def GetCellPoints(self, id):
    """Returns an array with the node numbers of each cell (ndglno)."""
    idlist=tvtk.IdList()
    self.ugrid.get_cell_points(id, idlist)
    return arr([idlist.get_id(i) for i in range(idlist.number_of_ids)])
    
  def GetFieldNames(self):
    """Returns the names of the available fields."""
    pointdata=self.ugrid.point_data
    return [pointdata.get_array_name(i) for i in range(pointdata.number_of_arrays)]

  def GetPointCells(self, id):
    """Return an array with the elements which contain a node."""
    idlist=tvtk.IdList()
    self.ugrid.get_point_cells(id, idlist)
    return arr([idlist.get_id(i) for i in range(idlist.number_of_ids())])

  def GetPointPoints(self, id):
    """Return the nodes connecting to a given node."""
    cells = self.GetPointCells(id)
    lst = []
    for cell in cells:
      lst = lst + list(self.GetCellPoints(cell))

    s = set(lst) # remove duplicates
    return arr(list(s)) # make into a list again

  def GetDistance(self, x, y):
    """Return the distance in physical space between x and y."""
    posx = self.ugrid.get_point(x)
    posy = self.ugrid.get_point(y)
    return math.sqrt(sum([(posx[i] - posy[i])**2 for i in range(len(posx))]))

  def Crop(self, min_x, max_x, min_y, max_y, min_z, max_z):
    """Trim off the edges defined by a bounding box."""
    trimmer = tvtk.ExtractUnstructuredGrid(input=self.ugrid, extent=(min_x, max_x, min_y, max_y, min_z, max_z))
    trimmer.update()
    trimmed_ug = trimmer.output

    self.ugrid = trimmed_ug
    
  def StructuredPointProbe(self, nx, ny, nz, bounding_box=None):
    """ Probe the unstructured grid dataset using a structured points dataset. """
       
    bbox = [0.0,0.0, 0.0,0.0, 0.0,0.0]
    if bounding_box==None:
      bbox = self.ugrid.bounds
    else:
      bbox = bounding_box

    spacing = [0.0, 0.0, 0.0]

    if nx>1: spacing[0] = (bbox[1]-bbox[0])/(nx-1.0)
    if ny>1: spacing[1] = (bbox[3]-bbox[2])/(ny-1.0)
    if nz>1: spacing[2] = (bbox[5]-bbox[4])/(nz-1.0)

    sgrid = tvtk.StructuredPoints(dimensions=(nx, ny, nz), origin=[bbox[0],bbox[2],bbox[4]], spacing=spacing)

    probe = tvtk.ProbeFilter (source=self.ugrid, input=sgrid)
    
    probe.update()
    
    return probe.output
