#include "fdebug.h"

subroutine test_anisotropic_zz
  use fields
  use anisotropic_zz_module
  use mba_adapt_module
  use mesh_files
  use vtk_interfaces
  use edge_length_module
  use bounding_box_metric
  use field_options
  use state_module
  use form_metric_field
  use interpolation_error
  use huang_metric_module
  use populate_state_module
  use adapt_state_module
  use global_parameters
  implicit none

  type(vector_field) :: positions
  type(scalar_field) :: u
  type(state_type) :: state
  type(tensor_field) :: metric
  integer :: loop, stat
  logical :: fail
  real :: h1
  real :: eta, tau
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

  positions = read_mesh_files("data/unit_metric", quad_degree=4, format="gmsh")
  call initialise_bounding_box_metric(positions)
  call insert(state, positions, "Coordinate")
  call deallocate(positions)

  adaptivity_mesh_name = "Mesh"

  call compute_domain_statistics((/state/))

  tau = 0.8

  do loop=1,10

    state%name = "StateNameHere"
    positions = extract_vector_field(state, "Coordinate")
    if (.not. has_faces(positions%mesh)) then
      call add_faces(positions%mesh)
    end if
    call insert(state, positions%mesh, "Mesh")
    call insert(state, positions%mesh, "CoordinateMesh")
    call allocate(u, positions%mesh, "U")
    call set_from_function(u, solution, positions)
    u%option_path = "/fields/u"
    call set_option("/fields/u/prognostic/adaptivity_options/anisotropic_zienkiewicz_zhu/tau", tau, stat=stat)
    call set_option("/fields/u/prognostic/adaptivity_options/huang_metric/seminorm", (/1, 2/), stat=stat)
    call set_option("/mesh_adaptivity/hr_adaptivity/maximum_number_of_nodes", 100000, stat=stat)

    call allocate(metric, positions%mesh, "Metric")
    call zero(metric)

    call compute_anisotropic_zz_metric(u, positions, metric, eta_estimate=eta)
!    call compute_hessian(u, positions, metric)
!    call form_huang_metric(metric, u, positions, tau)
    h1 = compute_interpolation_error_h1(gradsoln, u, positions)
    ewrite(2,*) "h1: ", h1
    call adaptivity_bounds(state, 0.000001, 10.0)
    call bound_metric(metric, state)

    call deallocate(u)

    call adapt_state(state, metric)

    positions = extract_vector_field(state, "Coordinate")
    call allocate(u, positions%mesh, "U")
    call set_from_function(u, solution, positions)
    call vtk_write_fields("anisotropic_zz", loop, positions, positions%mesh, sfields=(/u/))
    call deallocate(u)

  end do

  fail = ((eta - tau)/tau > 0.01)
  call report_test("[anisotropic_zz]", fail, .false., "")

end subroutine test_anisotropic_zz

function solution(pos) result(val)
  real, dimension(:), intent(in) :: pos
  real :: val

  real :: x, y

  x = pos(1); y = pos(2)
  !val = (4-4*exp(-100.*x)-4*(1-exp(-100.))*x)*y*(1-y)
  val = (x-0.5)**2 + (y-0.5)**2
end function solution

function gradsoln(pos) result(val)
  real, dimension(:), intent(in) :: pos
  real, dimension(size(pos)) :: val

  real :: x, y
  x = pos(1); y = pos(2)
  val(1) = 2*(x - 0.5)
  val(2) = 2*(y - 0.5)
  !val(1) = (400.0*exp(-(100.0*x)) - 4.0)*(1.0 - y)*y
  !val(2) = (-4*exp(-(100.0*x)) - 4.0*x + 4)*(1 - y) - (-4*exp(-(100.0*x)) - 4.0*x + 4)*y
end function gradsoln
