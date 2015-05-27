subroutine test_jacobian
  !!< test computation of jacobian for a 2D triangle embedded in 3D space

  use transform_elements
  use elements
  use quadrature
  use shape_functions
  use fields
  use unittest_tools

  implicit none

  integer, parameter :: dim=2, loc=3, quad_degree=4
  integer, parameter :: xdim=3

  type(mesh_type):: mesh
  type(vector_field):: X
  type(element_type), pointer :: shape
  type(quadrature_type), pointer :: quad
  real, allocatable, dimension(:,:,:) :: J
  real, allocatable, dimension(:) :: detwei
  logical :: fail

  allocate(quad)
  allocate(shape)

  quad=make_quadrature(loc, dim, quad_degree)
  shape=make_element_shape(loc, dim, 1, quad)
  allocate(J(dim,xdim,shape%ngi))
  allocate(detwei(shape%ngi))

  ! create single triangle mesh
  call allocate(mesh, loc, 1, shape, "OneElementMesh")
  call set_ele_nodes(mesh, 1, (/1,2,3/))

  ! and 3D positions field on it
  call allocate(X, 3, mesh, "Coordinate")
  call set(X, 1, (/0.0, 0.0, 0.0/))
  call set(X, 2, (/2.0, 0.0, 0.0/))
  call set(X, 3, (/0.0, 1.0, 1.0/))

  ! compute jacobian
  call compute_jacobian(X, 1, J, detwei=detwei)

  fail = abs(sum(detwei)-sqrt(2.0))>1e-10
  call report_test("[compute_jacobian]", fail, .false., "Incorrect Jacobian")

end subroutine test_jacobian
