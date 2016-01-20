subroutine test_compute_g_hat

  use fields
  use sparse_tools
  use adjacency_lists
  use mesh_files
  use anisotropic_zz_module
  use unittest_tools
  implicit none

  type(scalar_field) :: field
  type(vector_field) :: positions
  type(mesh_type) :: pwc_mesh
  integer :: ele
  real, dimension(2,2) :: g_hat
  logical :: fail

  interface
    function linear_func(pos) result(val)
      real, dimension(:), intent(in) :: pos
      real :: val
    end function linear_func
  end interface

  positions = read_mesh_files("data/pslgA", quad_degree=4, format="gmsh")
  call add_nelist(positions%mesh)
  pwc_mesh = piecewise_constant_mesh(positions%mesh, "PiecewiseConstantPlease")
  call allocate(field, positions%mesh, "Field")
  call set_from_function(field, linear_func, positions)

  fail = .false.
  do ele=1,ele_count(positions)
    g_hat = compute_g_hat(field, positions, ele_shape(pwc_mesh, ele), ele)
    fail = fail .or. (any(abs(g_hat) > epsilon(0.0)))
  end do

  call report_test("[compute g_hat]", fail, .false., "For a linear function, recovered gradient should be the same")


end subroutine test_compute_g_hat

function linear_func(pos) result(val)
  real, dimension(:), intent(in) :: pos
  real :: val

  real :: x, y
  x = pos(1); y = pos(2)
  val = x + y
end function linear_func
