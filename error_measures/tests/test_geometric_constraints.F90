
subroutine test_geometric_constraints
  use fields
  use mesh_files
  use unittest_tools
  use geometric_constraints_metric
  use vtk_interfaces
  implicit none

  type(vector_field), target :: X
  type(mesh_type), pointer :: mesh
  type(tensor_field) :: metric
  real, dimension(3, 3) :: init
  logical :: fail
  real :: a, b
  integer :: stat

  X=read_mesh_files("data/tet", quad_degree=4, format="gmsh")
  mesh => X%mesh

  call allocate(metric, mesh, "ErrorMetric")
  init = get_matrix_identity(3) / 100**2
  call set(metric, init)

  use_geometric_constraints_metric = .true.

  call vtk_write_fields("data/geometric_constraints", 0, X, mesh, tfields=(/metric/))

  call form_geometric_constraints_metric(X, metric)

  call vtk_write_fields("data/geometric_constraints", 1, X, mesh, tfields=(/metric/))

end subroutine test_geometric_constraints
