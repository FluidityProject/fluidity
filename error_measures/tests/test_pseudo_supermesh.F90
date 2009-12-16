
subroutine test_pseudo_supermesh
  use fields
  use unittest_tools
  use pseudo_supermesh
  use vtk_interfaces
  use reference_counting
  implicit none

  type(vector_field), target :: X_init, X_supermesh
  type(state_type) :: vtk_state
  logical :: fail
  character(len=255), dimension(2) :: strings

  strings(1) = "data/pseudo_supermesh_0.vtu"
  strings(2) = "data/pseudo_supermesh_1.vtu"

  call vtk_read_state(strings(1), vtk_state)
  X_init = extract_vector_field(vtk_state, "Coordinate")

  call add_faces(X_init%mesh)

  call compute_pseudo_supermesh(strings, &
                              & X_init, X_supermesh)

  call vtk_write_fields("data/pseudo_supermesh", 2, X_supermesh, X_supermesh%mesh)

  fail = (X_supermesh%refcount%count /= 1)
  call report_test("[supermesh refcount]", fail, .false., "")

  fail = (X_init%refcount%count /= 1)
  call report_test("[initial mesh refcount]", fail, .false., "")

  call deallocate(vtk_state)
  call deallocate(X_supermesh)
  call print_references(-1)

end subroutine test_pseudo_supermesh
