subroutine compute_interpolation_error_noadapt

  use metric_assemble
  use adapt_state_module
  use form_metric_field
  use state_module
  use vtk_interfaces
  use vector_tools
  use unittest_tools
  use edge_length_module
  use gradation_metric
  use mpi
  use interpolation_error
  implicit none
  
  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field) :: field
  real :: inf, l2, h1

  interface
    function solution(pos)
      real, dimension(:) :: pos
      real :: solution
    end function
  end interface
  interface
    function gradsoln(pos)
      real, dimension(:) :: pos
      real, dimension(size(pos)) :: gradsoln
    end function
  end interface

  call vtk_read_state("data/1x1square.vtu", state)

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  call allocate(field, mesh, "Field")

  call set_from_function(field, solution, positions)

  inf = compute_interpolation_error_inf(solution, field, positions)
  l2 =  compute_interpolation_error_l2(solution, field, positions)
  h1  = compute_interpolation_error_h1(gradsoln, field, positions)
  write(0,*) "inf: (", mesh%nodes, ", ", inf, ")"
  write(0,*) "l2: (", mesh%nodes, ", ", l2, ")"
  write(0,*) "h1: (", mesh%nodes, ", ", h1, ")"
end subroutine compute_interpolation_error_noadapt

function solution(pos)
  real :: solution
  real, dimension(:) :: pos
  x = pos(1); y = pos(2); z = pos(3)

  solution = y * x**2  + y**3 + tanh(10.0 * (sin(5.0*y) - 2.0*x)) + 4.0
end function solution

function gradsoln(pos)
  real, dimension(:), intent(in) :: pos
  real, dimension(size(pos)) :: gradsoln
  real :: x, y, z, dx, dy, dz

  x = pos(1); y = pos(2); z = pos(3)
  dx = 2*x*y - 20.0*sech(10.0*(sin(5.0*y) - 2.0*x))**2
  dy = 50.0*cos(5.0*y)*sech(10.0*(sin(5.0*y) - 2.0*x))**2 + 3*y**2 +  x**2
  dz = 0.0
  gradsoln(1) = dx; gradsoln(2) = dy; gradsoln(3) = dz
end function gradsoln

function sech(x)
  real, intent(in) :: x
  real :: sech

  sech = 1.0 / cosh(x)
end function sech
