subroutine test_linear_interpolation

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
  type(tensor_field) :: old_tensor, new_tensor
  type(vector_field) :: old_vector, new_vector
  integer :: node, i
  real :: x, y, z
  real, dimension(3) :: pos
  real, dimension(3, 3) :: matrix
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
    call insert(state_3, old_fields(i), old_fields(i)%name)
    call deallocate(old_fields(i))

    call allocate(new_fields(i), new_mesh, "Temperature" // trim(buf))
    call insert(state_5, new_fields(i), new_fields(i)%name)
    call deallocate(new_fields(i))
  end do

  call allocate(old_tensor, old_mesh, "Viscosity")
  call insert(state_3, old_tensor, "Viscosity")
  call deallocate(old_tensor)

  call allocate(new_tensor, new_mesh, "Viscosity")
  call insert(state_5, new_tensor, "Viscosity")
  call deallocate(new_tensor)

  call allocate(old_vector, 3, old_mesh, "Velocity")
  call insert(state_3, old_vector, "Velocity")
  call deallocate(old_vector)

  call allocate(new_vector, 3, new_mesh, "Velocity")
  call insert(state_5, new_vector, "Velocity")
  call deallocate(new_vector)

  do node=1,node_count(old_mesh)
    pos = node_val(old_position, node)
    x = pos(1) ; y = pos(2) ; z = pos(3)
    old_fields(1)%val(node) = x
    old_fields(2)%val(node) = y
    old_fields(3)%val(node) = x + y
    old_tensor%val(:, :, node) = 0.0
    old_tensor%val(1, 1, node) = x
    old_tensor%val(2, 2, node) = y
    old_tensor%val(3, 3, node) = z
    old_vector%val(1,node) = x
    old_vector%val(2,node) = y
    old_vector%val(3,node) = z
  end do

  call linear_interpolation(state_3, state_5)

  call vtk_write_fields("data/linear_interpolation", 0, old_position, old_mesh, sfields=old_fields, vfields=(/old_position/))
  call vtk_write_fields("data/linear_interpolation", 1, new_position, new_mesh, sfields=new_fields, vfields=(/new_position/))

  fail = .false.
  do node=1,node_count(new_mesh)
    pos = node_val(new_position, node)
    x = pos(1) ; y = pos(2) ; z = pos(3)
    if (.not. fequals(node_val(new_fields(1), node), x, 0.01)) then
      write(0,*) "node == ", node
      write(0,*) "node_val(new_fields(1), node) == ", node_val(new_fields(1), node)
      write(0,*) "x**2 == ", x**2
      fail = .true.
    end if
    if (.not. fequals(node_val(new_fields(2), node), y, 0.01)) then
      write(0,*) "node == ", node
      write(0,*) "node_val(new_fields(2), node) == ", node_val(new_fields(2), node)
      write(0,*) "y**2 == ", y**2
      fail = .true.
    end if
    if (.not. fequals(node_val(new_fields(3), node), x + y, 0.01)) then
      write(0,*) "node == ", node
      write(0,*) "node_val(new_fields(3), node) == ", node_val(new_fields(3), node)
      write(0,*) "x**2 + y**2 == ", x**2 + y**2
      fail = .true.
    end if
    if (.not. fequals(node_val(new_vector, 1, node), x, 0.01)) then
      write(0,*) "node == ", node
      fail = .true.
    end if
    if (.not. fequals(node_val(new_vector, 2, node), y, 0.01)) then
      write(0,*) "node == ", node
      fail = .true.
    end if
    if (.not. fequals(node_val(new_vector, 3, node), z, 0.01)) then
      write(0,*) "node == ", node
      fail = .true.
    end if
    matrix = node_val(new_tensor, node)
    if (.not. fequals(matrix(1, 1), x, 0.01)) then
      write(0,*) "node == ", node
      write(0,*) "node_val(new_tensor, node) == ", node_val(new_tensor, node)
      write(0,*) "x == ", x
      fail = .true.
    end if
    if (.not. fequals(matrix(2, 2), y, 0.01)) then
      write(0,*) "node == ", node
      write(0,*) "node_val(new_tensor, node) == ", node_val(new_tensor, node)
      write(0,*) "y == ", y
      fail = .true.
    end if
    if (.not. fequals(matrix(3, 3), z, 0.01)) then
      write(0,*) "node == ", node
      write(0,*) "node_val(new_tensor, node) == ", node_val(new_tensor, node)
      write(0,*) "z == ", z
      fail = .true.
    end if
  end do
  call report_test("[linear interpolation]", fail, .false., "All nodal values should be exact.")

end subroutine test_linear_interpolation
