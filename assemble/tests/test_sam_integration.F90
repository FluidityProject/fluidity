subroutine test_sam_integration

  use fields
  use mesh_files
  use conservative_interpolation_module
  use unittest_tools
  use state_module
  use vtk_interfaces
  use sam_integration
  use spud
  use reference_counting
  use populate_state_module
  implicit none

  type(vector_field) :: test_field_v
  type(scalar_field), pointer :: ptr_field_s
  type(vector_field), pointer :: ptr_field_v
  type(mesh_type), pointer :: mesh_ptr
  type(state_type), dimension(:), pointer :: states
  logical :: fail
  integer, dimension(10) :: options
  integer :: stat

  call load_options("data/sam_integration.flml")
  call populate_state(states)

  mesh_ptr => extract_mesh(states(1), "VelocityMesh")
  call allocate(test_field_v, 3, mesh_ptr, "TestVector")
  call set(test_field_v, (/2.0, 3.0, 4.0/))
  test_field_v%option_path = "/material_phase[0]/vector_field::TestVector"
  call set_option("/geometry/dimension", 3, stat=stat)
  call set_option("/material_phase[0]/vector_field::TestVector/name", "TestVector", stat=stat)
  call set_option("/material_phase[0]/vector_field::TestVector/prognostic", "1.0", stat=stat)
  call set_option("/material_phase[0]/vector_field::TestVector/prognostic/mesh[0]/name", "CoordinateMesh", stat=stat)
  call insert(states(1), test_field_v, "TestVector")
  call deallocate(test_field_v)

  call vtk_write_state("data/sam_integration", 0, state=states)

  options = (/0, 4, 1, 1, 1, 1, 0, 0, 0, 0/)
  call sam_drive(states, options)

  call vtk_write_state("data/sam_integration", 1, state=states)

  ptr_field_s => extract_scalar_field(states(1), "TestScalar")
  fail = (node_val(ptr_field_s, 1) /= 1.0)
  call report_test("[sam_integration scalar]", fail, .false., "Keep a constant scalar field.")

  ptr_field_v => extract_vector_field(states(1), "TestVector")
  fail = any(node_val(ptr_field_v, 1) /= (/2.0, 3.0, 4.0/))
  call report_test("[sam_integration vector]", fail, .false., "Keep a constant vector field.")
end subroutine test_sam_integration
