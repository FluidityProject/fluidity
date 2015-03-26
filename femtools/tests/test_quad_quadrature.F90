subroutine test_quad_quadrature
  
  use unittest_tools
  use elements
  use fetools
  use shape_functions
  use fields
  implicit none

  logical :: fail
  type(mesh_type) :: tri_mesh, quad_mesh
  type(vector_field) :: tri_X, quad_X
  type(element_type) :: quad_shape, tri_shape
  type(quadrature_type) :: quad_quadrature, tri_quadrature
  real, dimension(:,:,:), allocatable :: J_quad, J_tri
  real, dimension(:), allocatable :: detwei_quad, detwei_tri
  real, dimension(4,4) :: quad_mass
  real, dimension(6,6) :: l_tri_mass
  real, dimension(4,4) :: global_tri_mass
  real, dimension(4,6) :: local2global
  integer :: i

  quad_quadrature = make_quadrature(vertices=4,dim=2,degree=2)
  tri_quadrature = make_quadrature(vertices=3,dim=2,degree=4)

  allocate(J_quad(2,2,quad_quadrature%ngi))
  allocate(J_tri(2,2,tri_quadrature%ngi))
  allocate(detwei_quad(quad_quadrature%ngi))
  allocate(detwei_tri(tri_quadrature%ngi))
 
  quad_shape=make_element_shape(vertices=4, dim=2, degree=1, &
       &quad= quad_quadrature)
  tri_shape=make_element_shape(vertices=3, dim=2, degree=2, &
       &quad= tri_quadrature)

  !This unit test is based on the fact that the Q1 space on a single
  !quadrilateral can be represented exactly by the P2 space on the same
  !quadrilateral subdivided into two quadratically-mapped triangles.  The
  !dividing line between the triangles is quadratic and passes through the
  !mean of the four vertices.

  ! create single quad mesh
  ! numbering: 3 4
  !            1 2
  call allocate(quad_mesh, 4, 1, quad_shape, "OneElementMesh")
  call set_ele_nodes(quad_mesh, 1, (/1,2,3,4/))

  ! and 2D positions field on it
  call allocate(quad_X, 2, quad_mesh, "Coordinate")
  call set(quad_X, 1, (/0.0, 0.0/))
  call set(quad_X, 2, (/1.0, 0.0/))
  call set(quad_X, 3, (/0.2, 1.2/))
  call set(quad_X, 4, (/1.3, 1.5/))

  !compute the mass matrix using quadrilateral quadrature
  call compute_jacobian(quad_X, 1, J=J_quad,detwei=detwei_quad)
  quad_mass = shape_shape(quad_shape,quad_shape,detwei_quad)

  ! create 2 P2 triangles mesh
  !numbering: 6 5 3         3  4
  !           4 2 8
  !           1 7 9         1  2
  call allocate(tri_mesh, 9, 2, tri_shape, "TwoElementMesh")
  call set_ele_nodes(tri_mesh, 1, (/ 1,2,3,4,5,6 /))
  call set_ele_nodes(tri_mesh, 2, (/ 3,2,1,8,7,9 /))

  ! and 2D positions field on it
  call allocate(tri_X, 2, tri_mesh, "Coordinate")
  call set(tri_X, 1, node_val(quad_X, 1))
  call set(tri_X, 2, sum(quad_X%val, 2)/4.0)
  call set(tri_X, 3, node_val(quad_X, 4))
  call set(tri_X, 4, (node_val(quad_X, 1)+node_val(quad_X,3))/2.0)
  call set(tri_X, 5, (node_val(quad_X, 3)+node_val(quad_X,4))/2.0)
  call set(tri_X, 6, node_val(quad_X, 3))
  call set(tri_X, 7, (node_val(quad_X, 1)+node_val(quad_X,2))/2.0)
  call set(tri_X, 8, (node_val(quad_X, 2)+node_val(quad_X,4))/2.0)
  call set(tri_X, 9, node_val(quad_X, 2))

  !compute the mass matrix using triangular quadrature
  call compute_jacobian(tri_X, 1, J=J_tri,detwei=detwei_tri)
  l_tri_mass = shape_shape(tri_shape,tri_shape,detwei_tri)
  !local2global(i,:) gives coefficients of expansion of Q1 basis function i
  !into P2 basis functions in this triangle
  local2global(1,:) = (/1.,0.25,0.,0.5,0.,0./)
  local2global(2,:) = (/0.,0.25,0.,0.,0.,0./)
  local2global(3,:) = (/0.,0.25,0.,0.5,0.5,1./)
  local2global(4,:) = (/0.,0.25,1.,0.,0.5,0./)
  global_tri_mass = matmul(local2global,matmul(l_tri_mass,transpose(local2global)))

  call compute_jacobian(tri_X, 2, J=J_tri,detwei=detwei_tri)
  l_tri_mass = shape_shape(tri_shape,tri_shape,detwei_tri)
  !local2global(i,:) gives coefficients of expansion of Q1 basis function i
  !into P2 basis functions in this triangle
  local2global(1,:) = (/0.,0.25,1.,0.,0.5,0./)
  local2global(2,:) = (/0.,0.25,0.,0.5,0.5,1./)
  local2global(3,:) = (/0.,0.25,0.,0.,0.,0./)
  local2global(4,:) = (/1.,0.25,0.,0.5,0.,0./)
  global_tri_mass = global_tri_mass + &
       matmul(local2global,matmul(l_tri_mass,transpose(local2global)))

  fail = any(abs(quad_mass - global_tri_mass) > 1e-12)
  call report_test("[quad_quadrature]", fail, .false., "matrices not the same")

end subroutine test_quad_quadrature
