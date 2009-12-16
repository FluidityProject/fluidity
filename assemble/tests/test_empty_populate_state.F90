subroutine test_empty_populate_state
  use spud
  use populate_state_module
  use global_parameters
  use state_module
  use unittest_tools
  use vtk_interfaces
  implicit none

  type(state_type), dimension(:), pointer :: states => null()
  logical :: fail

  is_active_process = .false.
  call load_options("data/empty-mesh.flml")
  call populate_state(states)

  call report_test("[empty populate state]", .false., .false., "Crashing is failing")

  fail = (scalar_field_count(states(1)) == 0)
  call report_test("[empty populate state]", fail, .false., "Should have scalar fields")

  fail = (vector_field_count(states(1)) == 0)
  call report_test("[empty populate state]", fail, .false., "Should have vector fields")

  call vtk_write_state("data/empty_populate_state", state=states)
end subroutine test_empty_populate_state
