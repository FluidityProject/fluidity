subroutine test_trace
  
  use fields
  use field_derivatives
  use vtk_interfaces
  use state_module
  use unittest_tools
  implicit none

  type(scalar_field) :: tensor_trace
  type(tensor_field) :: tensor
  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  logical :: fail

  interface
    function solution(pos)
      real, dimension(:) :: pos
      real, dimension(size(pos), size(pos)) :: solution
    end function
  end interface

  call vtk_read_state("data/pseudo2d.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  call allocate(tensor, mesh, "Tensor")
  call allocate(tensor_trace, mesh, "Trace")

  call set_from_function(tensor, solution, positions)
  call trace(tensor, tensor_trace)

  fail = any(tensor_trace%val /= 3.0)

  call report_test("[trace]", fail, .false., "trace(matrix of 1s) == dim")

end subroutine test_trace

function solution(pos)
  real, dimension(:) :: pos
  real, dimension(size(pos), size(pos)) :: solution
  real :: x,y,z
  x = pos(1); y = pos(2); z = pos(3)

  solution = 1.0
end function solution
