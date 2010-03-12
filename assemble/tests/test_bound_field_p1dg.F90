subroutine test_bound_field_p1dg

  use populate_state_module
  use fields
  use state_module
  use spud
  use state_fields_module
  use sparse_tools
  use bound_field_module
  use vtk_interfaces
  use unittest_tools

  implicit none

  type(state_type), dimension(:), pointer :: states
  type(scalar_field), pointer :: u, max_bound, min_bound, lumped_mass
  type(scalar_field) :: inverse_lumped_mass
  type(csr_matrix), pointer :: mass
  type(vector_field), pointer :: x
  integer :: ele
  integer, dimension(:), pointer :: nodes
  integer :: node
  logical :: fail

  call load_options("data/bound_field.flml")
  call populate_state(states)

  u => extract_scalar_field(states(1), "Unbounded")
  max_bound => extract_scalar_field(states(1), "MaxBound")
  min_bound => extract_scalar_field(states(1), "MinBound")

  do ele=1,ele_count(u)
    nodes => ele_nodes(u, ele)
    call set(u, nodes(1), 1.2)
  end do

  mass => get_mass_matrix(states, u%mesh)
  lumped_mass => get_lumped_mass(states, u%mesh)

  call allocate(inverse_lumped_mass, lumped_mass%mesh, "InverseLumpedMass")
  inverse_lumped_mass%val = 1.0/lumped_mass%val

  x => extract_vector_field(states(1), "Coordinate")

  call vtk_write_fields("data/bounding", 0, x, u%mesh, sfields=(/u/))
  call bound_field_diffuse(u, max_bound, min_bound, mass, lumped_mass, inverse_lumped_mass)
  call vtk_write_fields("data/bounding", 1, x, u%mesh, sfields=(/u/))

  fail = .false.
  do node=1,node_count(u)
    if (node_val(u, node) > 1.0) then
      fail = .true.
    end if
  end do

  call report_test("[bound_field_p1dg]", fail, .false., "")

end subroutine test_bound_field_p1dg
