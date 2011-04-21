import numpy,sys,copy,operator

class State:
  def __init__(self,n=""):
    self.scalar_fields = {}
    self.vector_fields = {}
    self.tensor_fields = {}
    self.meshes = {}
    self.name = n
  def __repr__(self):
    return '(State) %s' % self.name
  def print_fields(self):
    print "scalar: ",self.scalar_fields,"\nvector:",self.vector_fields,"\ntensor:",self.tensor_fields,"\n"

class Field:
  def __init__(self,n,ft,op,description):
    self.name = n
    self.field_type = ft
    self.option_path = op
    self.description = description

  def __repr__(self):
    return '(%s) %s' % (self.description,self.name)

  def set_mesh(self,mesh):
    self.mesh = mesh
    self.element_count = self.mesh.element_count
    self.ele_count = self.mesh.element_count

  def shape(self):
    return self.mesh.shape

  def ele_loc(self,ele_number):
    # Return the number of nodes of element ele_number
    return self.mesh.ele_loc(ele_number)

  def ele_nodes(self,ele_number):
    # Return a pointer to a vector containing the global node numbers of
    # element ele_number, i.e. all node numbers associated with the element ele_number
    return self.mesh.ele_nodes(ele_number)

  def ele_shape(self,ele_number):
    return self.mesh.shape

  def addto(self, node, val):
      '''Add val to node of this field. If node is a scalar then val must
      have the shape of one data item in this field. If node is a sequence then
      the leading dimension of val must match the length of node and the
      remaining dimensions must match the shape of a data item in this field.'''
      try:
          for ii,vv in zip(node,val):
              self.val[ii]=self.val[ii]+vv
      except TypeError:
          # In this case it's presumably scalars.
          self.val[node]=self.val[node]+val

  def set(self, node, val):
      '''Set node of this field to val.  If node is a scalar then val must
      have the shape of one data item in this field. If node is a sequence then
      the leading dimension of val must match the length of node and the
      remaining dimensions must match the shape of a data item in this
      field.'''
      
      try:
          for ii,vv in zip(node,val):
              self.val[ii]=vv
      except TypeError:
          self.val[node]=val

  def __setitem__(self, node, val):
    self.set(node, val)
     
  def node_val(self,node):
    return self.val[node]

  def __getitem__(self, node):
    if node != int(node):
      raise TypeError
    if node > self.node_count:
      raise IndexError
    return self.node_val(node)

  def ele_val(self,ele_number):
    # Return the values of field at the nodes of ele_number
    return numpy.array(map(self.node_val,self.ele_nodes(ele_number)))

  def ele_val_at_quad(self,ele_number):
    # Return the values of field at the quadrature points of ele_number
    shape_n = numpy.matrix(self.ele_shape(ele_number).n)
    ele_val = numpy.matrix(self.ele_val(ele_number))
    #print "ele_val:",ele_val.shape
    #print "shape_n:",shape_n.shape
    return numpy.array(ele_val*shape_n)

  def ele_region_id(self,ele_number):
    return self.mesh.ele_region_id(ele_number)

  def remap_ele(self, ele_number, mesh):
    assert self.mesh.continuity >= mesh.continuity
    if mesh.continuity >= 0: # if we are CG
      assert self.mesh.shape.degree <= mesh.shape.degree
    assert not (self.mesh.shape.type == "bubble" and mesh.shape.type == "lagrangian")

    # we should check for periodic/nonperiodic here, but
    # our mesh type doesn't know whether it's periodic or not ...

    # we really ought to cache locweight, as it's constant
    # for each element, and only depends on the target mesh --
    # but we currently don't, sorry
    locweight = numpy.zeros((mesh.shape.loc, self.mesh.shape.loc))
    for i in range(mesh.shape.loc):
      for j in range(self.mesh.shape.loc):
        locweight[i,j] = self.mesh.shape.eval_shape(j, mesh.shape.local_coords(i))

    return numpy.dot(locweight, self.ele_val(ele_number))

class ScalarField(Field):
  "A scalar field"  
  description = "ScalarField"
  def __init__(self,n,v,ft,op):
    Field.__init__(self,n,ft,op,self.description)
    self.val = v
    self.node_count = self.val.shape[0]


class VectorField(Field):
  "A vector field"
  description = "VectorField"
  def __init__(self,n,v,ft,op,dim):
    Field.__init__(self,n,ft,op,self.description)
    self.val = v
    self.dimension = dim
    self.node_count=self.val.shape[0]
    
class TensorField(Field):
  "A tensor field"
  description = "VectorField"
  def __init__(self,n,v,ft,op,dim0,dim1):
    Field.__init__(self,n,ft,op,self.description)
    self.val = v
    self.dimension = numpy.array([dim0,dim1])
    self.node_count=self.val.shape[0]


class Mesh:
  "A mesh"
  def __init__(self,ndglno,elements,nodes,continuity,name,option_path,region_ids):
    self.ndglno = ndglno
    self.element_count = elements
    self.node_count = nodes
    self.continuity = continuity
    self.name = name
    self.option_path = option_path
    self.region_ids = region_ids
    self.shape = Element(0,0,0,0,[],[], [], 0,0,0,0, "unknown", "unknown")

  def __repr__(self):
    return '(Mesh) %s' % self.name

  def ele_shape(self,ele_number):
    # Returns the shape of this mesh
    return self.shape
  
  def ele_loc(self,ele_number):
    # Returns the loc of the shape of this mesh
    return self.shape.loc

  def ele_nodes(self,ele_number):
    # Return all nodes associated with the element ele_number
    base = self.shape.loc*ele_number
    nodes = []
    for i in range(self.shape.loc):
      # Subtract 1, since the nodes are numbered from 1 in ndglno
      nodes.append(self.ndglno[base+i]-1)
    return nodes

  def ele_region_id(ele_number):
    # Return the region_id of element ele_number
    return self.mesh.region_ids[ele_number]

class Element:
  "An element"
  def __init__(self,dim,loc,ngi,degree,n,dn, coords, size_spoly_x,size_spoly_y,size_dspoly_x,size_dspoly_y, family, type):
    self.dimension = dim  # 2d or 3d?
    self.loc = loc  # Number of nodes
    self.ngi = ngi  # Number of gauss points
    self.degree = degree  # Polynomial degree of element
    # Shape functions: n is for the primitive function, dn is for partial derivatives, dn_s is for partial derivatives on surfaces
    # n is loc x ngi, dn is loc x ngi x dim
    self.n = n
    self.dn = dn
    self.coords = coords
    self.family = family
    self.type = type

    # Initialize spoly and dspoly, make sure to transpose due to Fortran silliness
    self.spoly = [ [ Polynomial([],0) for inner in range(size_spoly_y) ] for outer in range(size_spoly_x)]
    self.dspoly = [ [ Polynomial([],0) for inner in range(size_dspoly_y) ] for outer in range(size_dspoly_x)]

  # Methods for setting up attributes
  def set_quadrature(self,quadrature):
    self.quadrature = quadrature
  def set_surface_quadrature(self,quadrature):
    self.surface_quadrature = quadrature
  def set_polynomial_s(self,poly,x,y):
      self.spoly[x-1][y-1] = poly
  def set_polynomial_ds(self,poly,x,y):
    self.dspoly[x-1][y-1] = poly
  def eval_shape(self, node, coords):
    result = 1.0
    for i in range(len(self.spoly)):
      result = result * self.spoly[i][node].eval(coords[i])
    return result
  def local_coords(self, node):
    return self.coords[node,:]

class Quadrature:
  "Quadrature"
  def __init__(self,w,locations,dim,degree,loc,ngi):
    self.dimension = dim  # Dimension of the elements for which quadrature is required.
    self.degree = degree  # Degree of accuracy of quadrature
    self.loc = loc        # Number of vertices of the element.
    self.ngi = ngi        # Number of quadrature points
    self.weights = w       # Quadrature weights
    self.locations = locations  # Locations of quadrature points


class Polynomial:
  "Polynomial"
  def __init__(self,coefs,degree):
    self.coefficients = coefs
    self.degree = degree
  def __repr__(self):
    return '(Polynomial)'+str(self.coefficients)
  def eval(self, x):
    val = 0.0
    for i in range(self.degree+1):
      val = val + self.coefficients[i] * x**i
    return val

class Transform:
  "Transform with information about the detwei and Jacobian"
  # Note that so far only the dim == ldim == (2||3) have been tested
  
  def __init__(self,ele_num,field):
    self.ele_num = ele_num
    self.element = field.mesh.shape
    self.field = field
    # Jacobian matrix and its inverse at each quadrature point (dim x dim x field.mesh.shape.ngi)
    # Facilitates access to this information externally
    self.J = [ [] for i in range(field.mesh.shape.ngi) ]
    self.invJ = [ [] for i in range(field.mesh.shape.ngi) ]
    
    # Calculate detwei, i.e. the gauss weights transformed by the coordinate transform
    # from real to computational space
    self.transform_to_physical_detwei(field)

  def set_J(self,J,gi):
    # Set J for the specified quadrature point and calculate its inverse
    self.J[gi] = numpy.matrix(J)
    if max(J.shape) != min(J.shape):
      # we don't have a square Jacobian, i.e. mesh_dim < coord_dim
      if hasattr(self, "invJ"):
        del self.invJ
    else:
      self.invJ[gi] = numpy.matrix(numpy.linalg.inv(self.J[gi]))

  def transform_to_physical_detwei(self,field):
     # field is the Coordinate Field

    # Column n of X is the position of the nth node. (dim x n%loc)
    # only need position of n nodes since Jacobian is only calculated once

    # element is the referenced velocity Element
    ele_num = self.ele_num
    X = numpy.transpose(numpy.matrix(field.ele_val(ele_num)))
    element = field.mesh.shape 

    # Quadrature weights for physical coordinates.
    self.detwei = numpy.zeros(element.ngi)
    dim = field.dimension     # Dimension of space
    ldim = element.dimension  # Dimension of element
    if(dim == ldim):
      if(dim==1):
        for gi in range(element.ngi):
          J = numpy.zeros([dim, ldim])
          J[0,0] = numpy.dot(X[0,:],element.dn[:,gi,0])  # Still wrong probably!
          self.detwei[gi] = abs(J[0,0])*element.quadrature.weights[gi]
          # The Jacobian is the transpose of the J that was calculated
          self.set_J(numpy.transpose(J),gi)
      elif(dim==2 or dim == 3):
        for gi in range(element.ngi):
          J = numpy.dot(X, element.dn[:,gi,:])
          self.detwei[gi] = abs( numpy.linalg.det(J)) * element.quadrature.weights[gi]
          # The Jacobian is the transpose of the J that was calculated
          self.set_J(numpy.transpose(J),gi)
      else:
        sys.exit("More than 3 dimensions!")

    # Lower dimensional element (ldim) embedded in higher dimensional space (dim)
    elif(ldim<dim):
      # 1-dim element embedded in 'dim'-dimensional space:
      if(ldim==1):
        for gi in range(element.ngi):
          J = numpy.zeros([dim, ldim])
          J[:,0] = numpy.dot(X, element.dn[:,gi,0])
          det = numpy.linalg.norm(J[:,0])
          self.detwei[gi] = det * element.quadrature.weights[gi]
          self.set_J(J,gi)

      # 1-dim element embedded in 'dim'-dimensional space:
      elif(ldim==2):
        # J is 2 columns of 2 'dim'-dimensional vectors:
        for gi in range(element.ngi):
          J = numpy.dot(X, element.dn[:,gi,:])
          # Outer product times quad. weight
          self.detwei[gi] = abs( J[1,0]*J[2,1]-J[2,0]*J[1,1]
                          - J[2,0]*J[0,1]+J[0,0]*J[2,1]
                          + J[0,0]*J[1,1]-J[1,0]*J[0,1]) * element.quadrature.weights[gi]
          self.set_J(numpy.transpose(J),gi)

    else:
      sys.exit("Dimension of shape exceeds dimension of coordinate field.")
    numpy.matrix(self.detwei)


  def grad(self,shape):
    # Evaluate derivatives in physical space

    # dn is loc x ngi x dim

    # Derivatives (dm_t and dn_t in Fortran)

    # Create a new object with the same values of 
    newshape = copy.copy(shape)
    newshape.dn = numpy.zeros((shape.loc,shape.ngi,self.field.dimension))

    for gi in range(self.field.mesh.shape.ngi):
      if (shape.dimension == self.field.dimension):
        for i in range(shape.loc):
          a = numpy.array(numpy.dot(shape.dn[i,gi], self.invJ[gi]))
          newshape.dn[i,gi,:] = a
      else:
        if shape.dimension == 1:
          invdetJ = 1.0/numpy.linalg.norm(self.J[gi][0,:])
          for i in range(shape.loc):
            newshape.dn[i,gi,:] = numpy.dot(self.J[gi], shape.dn[i, gi, :] * invdetJ**2)
        else:
          raise Exception, "Sorry, haven't worked out how to do this yet"
    return newshape


  def shape_shape(self,shape1,shape2):
    # For each node in each element shape1, shape2 calculate the
    # coefficient of the integral int(shape1shape2)dV.
    #
    # In effect, this calculates a mass matrix.
    m = numpy.zeros((shape1.loc,shape2.loc))

    for i in range(shape1.loc):
      for j in range(shape2.loc):
        m[i,j] = numpy.dot(shape1.n[i]*shape2.n[j], self.detwei)
    return m


  def shape_dshape(self,shape,dshape):
    # For each node in element shape and transformed gradient dshape
    # calculate the coefficient of the integral int(shape dshape)dV

    # The dimensions of dshape are: (nodes, gauss points, dimensions)

    # dshape is usually the element returned by calling element.grad()

    # dshape_loc = size(self.gradient)
    dshape_loc = len(dshape.dn)
    dshape_dim = dshape.dn.shape[2]
    field = self.field

    shape_dshape = numpy.zeros((shape.loc, dshape_loc, dshape_dim))
    for i in range(shape.loc):
      for j in range(dshape_loc):
        shape_dshape[i,j,:] = ( numpy.matrix(self.detwei) * numpy.matrix(self.spread(shape.n[i], dshape_dim)*dshape.dn[j])).reshape((dshape_dim,))
    return shape_dshape

  def spread(self,arr,dim):
    # Simple version of the Fortran spread function which returns 
    # an array of a higher dimension with the same values as the original one
    # in the new dimension
    if dim==1:
      return arr
    elif dim==2:
      return numpy.column_stack((arr,arr))
    elif dim==3:
      return numpy.column_stack((arr,arr,arr))



def test_shape_dshape(state):
  # Tests shape_dshape (surprised?) - pass in the state with the coordiate field
  cf = state.vector_fields['Coordinate']
  lumpmass = ScalarField("LumpMass", numpy.zeros(cf.node_count),"","")
  lumpmass.set_mesh(cf.mesh)
  
  # Set our linear function f to the x values of the Coordinate Mesh
  f = ScalarField("Linear",cf.val1,"","")
  f.set_mesh(cf.mesh)

  psi = VectorField("psi",numpy.zeros(cf.node_count),numpy.zeros(cf.node_count),numpy.zeros(cf.node_count),"","",2)
  psi.set_mesh(cf.mesh)

  for el in range(lumpmass.element_count()):
    el_f = lumpmass.ele_nodes(el)
    t = Transform(el,cf)
    masslump_local = t.shape_shape(cf.mesh.shape,cf.mesh.shape).sum(1)
    lumpmass.addto(el_f,masslump_local)
    rhs_local = t.shape_dshape(f.mesh.shape,t.grad(f.mesh.shape))
    rhs_f = 0
    for i in range(rhs_local.shape[1]):
      rhs_f = rhs_f + numpy.matrix(rhs_local[:,i,:]) * f.node_val(el_f[i])
    psi.addto(el_f,rhs_f)

  print psi.val[0]/lumpmass.val

  return 0
