
subroutine test_simplex_tensor
  use fields
  use mesh_files
  use unittest_tools
  use conformity_measurement
  use metric_tools
  implicit none

  type(vector_field) :: X_1, X_2
  real, dimension(1, 1) :: m_1
  real, dimension(2, 2) :: m_2, c_2
  logical :: fail
  X_1=read_mesh_files("data/interval-0.0-1.0-1.0", quad_degree=4, format="gmsh")

  m_1 = simplex_tensor(X_1, 1)

  fail = (m_1(1, 1) /= 1.0)
  call report_test("[test_simplex_tensor]", fail, .false., "")

  call deallocate(X_1)
  X_2=read_mesh_files("data/triangle.1", quad_degree=4, format="gmsh")

  m_2 = simplex_tensor(X_2, 1)
  c_2(1,:) = (/0.75, 0.0/)
  c_2(2,:) = (/0.0, 1.0/)

  fail = (m_2 .fne. c_2)
  call report_test("[test_simplex_tensor]", fail, .false., "")
end subroutine test_simplex_tensor
