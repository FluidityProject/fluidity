subroutine test_mba_supermesh

  use fields
  use populate_state_module
  use spud
  use state_module
  use form_metric_field
  use metric_assemble
  use adapt_state_module
  use field_derivatives
  use vtk_interfaces
  use solenoidal_interpolation_module
  use global_parameters
  use supermesh_construction
  implicit none

  interface
    function id(X) result(m)
      real, dimension(:), intent(in) :: X
      real, dimension(size(X), size(X)) :: m
    end function id
  end interface

  interface
    function nid(X) result(m)
      real, dimension(:), intent(in) :: X
      real, dimension(size(X), size(X)) :: m
    end function nid
  end interface

  type(state_type), dimension(:), pointer :: states_A => null()
  type(state_type), dimension(:), pointer :: states_B => null()
  type(tensor_field) :: metric
  type(mesh_type), pointer :: mesh_A, mesh_B
  type(vector_field), pointer :: x_A, x_B
  type(supermesh) :: sup

  call load_options("data/solenoidal_interpolation_A.flml")
  call populate_state(states_A)
  call populate_state(states_B)

  mesh_A => extract_mesh(states_A(1), "VelocityMesh")
  x_A => extract_vector_field(states_A(1), "Coordinate")

  call allocate(metric, mesh_A, "Metric")

  call set_from_function(metric, id, x_A)
  call adapt_state(states_A, metric)

  mesh_B => extract_mesh(states_B(1), "VelocityMesh")
  x_B => extract_vector_field(states_B(1), "Coordinate")

  call allocate(metric, mesh_B, "Metric")

  call set_from_function(metric, nid, x_B)
  call adapt_state(states_B, metric)

  call vtk_write_state("data/mba_supermesh", 0, state=states_A)
  call vtk_write_state("data/mba_supermesh", 1, state=states_B)

  x_A => extract_vector_field(states_A(1), "Coordinate")
  x_B => extract_vector_field(states_B(1), "Coordinate")
  sup = construct_supermesh(x_A, x_B)
  call vtk_write_fields("data/mba_supermesh", 2, sup%positions, sup%positions%mesh)

  call deallocate(states_A(1))
  call deallocate(states_B(1))
  call deallocate(sup)

end subroutine test_mba_supermesh

function id(X) result(m)
  real, dimension(:), intent(in) :: X
  real, dimension(size(X), size(X)) :: m
  integer :: i

  m = 0.0
  do i=1,size(X)
    m(i,i) = 1.0
  end do
end function id

function nid(X) result(m)
  real, dimension(:), intent(in) :: X
  real, dimension(size(X), size(X)) :: m
  integer :: i

  m = 0.0
  do i=1,size(X)
    m(i,i) = i
  end do
end function nid
