#include "fdebug.h"

subroutine test_anisotropic_zz_3d
  use fields
  use anisotropic_zz_module
  use mesh_files
  use vtk_interfaces
  use solvers
  use edge_length_module
  use bounding_box_metric
  use field_options
  use state_module
  use form_metric_field
  use adapt_state_unittest_module, only : adapt_state => adapt_state_unittest
  use interpolation_error
  implicit none

  type(vector_field) :: positions
  type(scalar_field) :: u
  type(state_type) :: state
  type(tensor_field) :: metric
  integer :: loop, stat
  character(len=255) :: path
  logical :: fail
  real :: tau, eta, h1
  interface
    function solution(pos) result(val)
      real, dimension(:), intent(in) :: pos
      real :: val
    end function solution
  end interface
  interface
    function gradsoln(pos)
      real, dimension(:) :: pos
      real, dimension(size(pos)) :: gradsoln
    end function gradsoln
  end interface

  call set_global_debug_level(2)

  positions = read_mesh_files("data/cube_anisotropic", quad_degree=4, format="gmsh")
  call initialise_bounding_box_metric(positions)
  call insert(state, positions, "Coordinate")
  call insert(state, positions%mesh, "Mesh")
  call deallocate(positions)

  tau = 1.0

  do loop=1,10

    positions = extract_vector_field(state, "Coordinate")
    call insert(state, positions%mesh, "Mesh")
    call allocate(u, positions%mesh, "U")
    call set_from_function(u, solution, positions)
    u%option_path = "/fields/u"
    call set_option("/fields/u/prognostic/adaptivity_options/anisotropic_zienkiewicz_zhu/tau", tau, stat=stat)
    path = "/fields/u/prognostic/adaptivity_options/anisotropic_zz"
    call set_solver_options(path, ksptype='cg', pctype='sor', rtol=1.0e-8, max_its=10000)

    call allocate(metric, positions%mesh, "Metric")
    call zero(metric)

    call compute_anisotropic_zz_metric(u, positions, metric, eta_estimate=eta)
    h1 = compute_interpolation_error_h1(gradsoln, u, positions)
    ewrite(2,*) "h1: ", h1
    if (eta > 0.75*tau .and. eta < 1.25*tau .and. loop /= 1) then
      exit
    end if
    call adaptivity_bounds(state, 0.000001, 1.0)
    call bound_metric(metric, state)

    call insert(state, u, "U")
    call deallocate(u)
    call set_global_debug_level(0)
    call adapt_state(state, metric)
    call set_global_debug_level(2)
    call deallocate(metric)
    u = extract_scalar_field(state, "U")
    positions = extract_vector_field(state, "Coordinate")
    call set_from_function(u, solution, positions)
    call vtk_write_state("data/anisotropic_zz_3d", loop, state=(/state/))
  end do

  fail = (ele_count(positions) > 300)
  call report_test("[anisotropic_zz]", fail, .false., "")

end subroutine test_anisotropic_zz_3d

function solution(pos) result(val)
  real, dimension(:), intent(in) :: pos
  real :: val

  real :: x, y, z
  real, parameter :: eps=0.01

  x = pos(1); y = pos(2); z = pos(3)
  !val = (x-0.5)**2 + (y-0.5)**2 + (z-0.5)**2

  val = exp(-x/eps) + exp(-y/eps) + exp(-z/eps)
end function solution

function gradsoln(pos) result(val)
  real, dimension(:), intent(in) :: pos
  real, dimension(size(pos)) :: val
  real :: x, y, z
  real, parameter :: eps=0.01
  x = pos(1); y = pos(2); z = pos(3)

  !val(1) = 2*(x-0.5)
  !val(2) = 2*(y-0.5)
  !val(3) = 2*(z-0.5)
  val(1) = -exp(-(x/eps))/eps
  val(2) = -exp(-(y/eps))/eps
  val(3) = -exp(-(z/eps))/eps
end function gradsoln
