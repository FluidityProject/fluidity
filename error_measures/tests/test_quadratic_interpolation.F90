subroutine test_quadratic_interpolation

  use fields
  use state_module
  use interpolation_module
  use vtk_interfaces
  use unittest_tools
  implicit none

  type(state_type) :: state_3, state_5
  type(mesh_type), pointer :: old_mesh, new_mesh
  type(vector_field), pointer :: old_position, new_position
  type(scalar_field), dimension(3) :: old_fields, new_fields
  integer :: node, i
  real :: x, y, z
  real, dimension(3) :: pos
  logical :: fail
  character(len=20) :: buf

  call vtk_read_state("data/cube-itv3.vtu", state_3)
  call vtk_read_state("data/cube-itv5.vtu", state_5)

  old_mesh => extract_mesh(state_3, "Mesh")
  new_mesh => extract_mesh(state_5, "Mesh")
  old_position => extract_vector_field(state_3, "Coordinate")
  new_position => extract_vector_field(state_5, "Coordinate")

  do i=1,3
    write(buf, '(i0)') i
    call allocate(old_fields(i), old_mesh, "Temperature" // trim(buf))
    call allocate(new_fields(i), new_mesh, "Temperature" // trim(buf))
  end do

  do node=1,node_count(old_mesh)
    pos = node_val(old_position, node)
    x = pos(1) ; y = pos(2) ; z = pos(3)
    old_fields(1)%val(node) = x**2
    old_fields(2)%val(node) = y**2
    old_fields(3)%val(node) = x**2 + y**2
  end do

  call quadratic_interpolation(old_fields, old_position, new_fields, new_position)

  call vtk_write_fields("data/quadratic_interpolation", 0, old_position, old_mesh, sfields=old_fields, vfields=(/old_position/))
  call vtk_write_fields("data/quadratic_interpolation", 1, new_position, new_mesh, sfields=new_fields, vfields=(/new_position/))

  fail = .false.
  do node=1,node_count(new_mesh)
    pos = node_val(new_position, node)
    x = pos(1) ; y = pos(2) ; z = pos(3)
    if (.not. fequals(node_val(new_fields(1), node), x**2, 0.01)) then
      write(0,*) "node == ", node
      write(0,*) "node_val(new_fields(1), node) == ", node_val(new_fields(1), node)
      write(0,*) "x**2 == ", x**2
      fail = .true.
    end if
    if (.not. fequals(node_val(new_fields(2), node), y**2, 0.01)) then
      write(0,*) "node == ", node
      write(0,*) "node_val(new_fields(2), node) == ", node_val(new_fields(2), node)
      write(0,*) "y**2 == ", y**2
      fail = .true.
    end if
    if (.not. fequals(node_val(new_fields(3), node), x**2 + y**2, 0.01)) then
      write(0,*) "node == ", node
      write(0,*) "node_val(new_fields(3), node) == ", node_val(new_fields(3), node)
      write(0,*) "x**2 + y**2 == ", x**2 + y**2
      fail = .true.
    end if
  end do
  call report_test("[quadratic interpolation]", fail, .false., "All nodal values should be exact.")

end subroutine test_quadratic_interpolation
