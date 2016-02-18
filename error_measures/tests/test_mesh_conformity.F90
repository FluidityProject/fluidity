
subroutine test_mesh_conformity
  use fields
  use mesh_files
  use unittest_tools
  use conformity_measurement
  implicit none

  type(vector_field), target :: X
  type(mesh_type), pointer :: mesh
  type(tensor_field) :: metric
  type(scalar_field) :: conformity
  logical :: fail
  real :: a, b

  X=read_mesh_files("data/eqtriangle.1", format="gmsh", quad_degree=4)
  mesh => X%mesh

  call allocate(metric, mesh, "ErrorMetric", field_type=FIELD_TYPE_CONSTANT)
  call set(metric, get_matrix_identity(X%dim))

  conformity = piecewise_constant_field(mesh, "MeshConformity")

  call compute_mesh_conformity(metric, X, conformity)
  
  fail = (node_val(conformity, 1) > 1e-7)
  call report_test("[test_mesh_conformity]", fail, .false., "Perfect meshes give zero.")

  call set(metric, 2 * get_matrix_identity(X%dim))
  call compute_mesh_conformity(metric, X, conformity)
  fail = (node_val(conformity, 1) .fne. (1 / sqrt(2.0)))
  call report_test("[test_mesh_conformity]", fail, .false., "Double the length scale: 1/sqrt(2.0)")

  call set(metric, 10 * get_matrix_identity(X%dim))
  call compute_mesh_conformity(metric, X, conformity)
  a = node_val(conformity, 1)
  call set(metric, 0.1 * get_matrix_identity(X%dim))
  call compute_mesh_conformity(metric, X, conformity)
  b = node_val(conformity, 1)

  fail = (a .fne. b)
  call report_test("[test_mesh_conformity]", fail, .false., "Too big and too short should be the same")
end subroutine test_mesh_conformity
