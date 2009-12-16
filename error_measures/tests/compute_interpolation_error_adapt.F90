subroutine compute_interpolation_error_adapt
#define ERROR 0.025

  use global_parameters, only: current_debug_level, pseudo2d_coord
  use metric_assemble
  use adapt_state_unittest_module, only : adapt_state => adapt_state_unittest
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
  type(scalar_field), pointer :: field_ptr
  type(tensor_field) :: metric, hessian
  type(metric_options) :: opts
  type(anisotropic_edge_options), target :: edge_opts
  real, dimension(4), target :: xhsamp, yhsamp, zhsamp
  real, dimension(4), target :: hminxx, hminxy, hminxz, hminyy, hminyz, hminzz
  real, dimension(4), target :: hmaxxx, hmaxxy, hmaxxz, hmaxyy, hmaxyz, hmaxzz
  integer :: nhsamp
  real :: inf, l2, h1
  integer :: i

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

  pseudo2d_coord = 3
  call vtk_read_state("data/mymesh.vtu", state)

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  call allocate(field, mesh, "Field")
  call zero(field)
  call set_from_function(field, solution, positions)
  call insert(state, field, "Field")
  field_ptr => extract_scalar_field(state, "Field")
  

  field_ptr%options%relative = .false.
  field_ptr%options%error = ERROR
  field_ptr%options%min_psi = 1e-10
  field_ptr%options%square_eigs = .false.

  call allocate(metric, mesh, "Metric")
  call allocate(hessian, mesh, "Hessian")

  opts%min_edge_length = 0.0005
  opts%max_edge_length = 2.0
  opts%use_anisotropic_edge_length = .true.
  xhsamp(1) = 0.0; yhsamp(1) = 0.0; zhsamp(1) = 0.0
  xhsamp(2) = 1.0; yhsamp(2) = 0.0; zhsamp(2) = 1.0
  xhsamp(3) = 1.0; yhsamp(3) = 1.0; zhsamp(3) = 0.0
  xhsamp(4) = 0.0; yhsamp(4) = 1.0; zhsamp(4) = 1.0
  hminxx = 0.0005; hminxy = 0.0; hminxz = 0.0; hminyy = 0.0005; hminyz = 0.0; hminzz = 0.00025
  hmaxxx = 1.0; hmaxxy = 0.0; hmaxxz = 0.0; hmaxyy = 1.0; hmaxyz = 0.0; hmaxzz = 0.002
  nhsamp = 4
  edge_opts%no_samp = NHSAMP
  edge_opts%x => XHSAMP(1:NHSAMP)
  edge_opts%y => YHSAMP(1:NHSAMP)
  edge_opts%z => ZHSAMP(1:NHSAMP)
  edge_opts%hminxx => HMINXX(1:NHSAMP); edge_opts%hmaxxx => HMAXXX(1:NHSAMP)
  edge_opts%hminxy => HMINXY(1:NHSAMP); edge_opts%hmaxxy => HMAXXY(1:NHSAMP)
  edge_opts%hminxz => HMINXZ(1:NHSAMP); edge_opts%hmaxxz => HMAXXZ(1:NHSAMP)
  edge_opts%hminyy => HMINYY(1:NHSAMP); edge_opts%hmaxyy => HMAXYY(1:NHSAMP)
  edge_opts%hminyz => HMINYZ(1:NHSAMP); edge_opts%hmaxyz => HMAXYZ(1:NHSAMP)
  edge_opts%hminzz => HMINZZ(1:NHSAMP); edge_opts%hmaxzz => HMAXZZ(1:NHSAMP)
  opts%anisotropic_edge_opts => edge_opts

  gamma0 = 2.5
  use_gradation_metric = .true.
  gradation_initialised = .true.
  call assemble_metric((/state/), metric, opts)
  call vtk_write_fields("data/interpolation_error_adapted", 0, positions, mesh, &
                        sfields=(/field/), tfields=(/metric/))
  call adapt_state(state, metric)

  do i=1,3
    mesh => extract_mesh(state, "Mesh")
    positions => extract_vector_field(state, "Coordinate")
    field_ptr => extract_scalar_field(state, "Field")
    call set_from_function(field_ptr, solution, positions)
    field_ptr%options%relative = .false.
    field_ptr%options%error = ERROR
    field_ptr%options%min_psi = 1e-10
    field_ptr%options%square_eigs = .false.
    call deallocate(metric); call allocate(metric, mesh, "Metric")
    call assemble_metric((/state/), metric, opts)
    call vtk_write_fields("data/interpolation_error_adapted", i, positions, mesh, &
                        sfields=(/field_ptr/))
    call adapt_state(state, metric) 
  end do

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  field_ptr => extract_scalar_field(state, "Field")

  write(0,*) ele_ngi(field_ptr, 1)
  call set_from_function(field_ptr, solution, positions)
  inf = compute_interpolation_error_inf(solution, field_ptr, positions)
  l2  = compute_interpolation_error_l2(solution, field_ptr, positions)
  h1  = compute_interpolation_error_h1(gradsoln, field_ptr, positions)
  call vtk_write_fields("data/interpolation_error_adapted", 4, positions,  mesh, sfields=(/field_ptr/)) 
  write(0,*) "field%options%error == ", field%options%error
  write(0,*) "inf: (", mesh%nodes, ", ", inf, ")"
  write(0,*) "l2: (", mesh%nodes, ", ", l2, ")"
  write(0,*) "h1: (", mesh%nodes, ", ", h1, ")"
end subroutine compute_interpolation_error_adapt

function solution(pos)
  real :: solution
  real, dimension(:) :: pos
  real :: x,y,z
  real, parameter :: PI=4.0*atan(1.0)
  x = pos(1); y = pos(2); z = pos(3)

  solution = cos(11*PI*x) + sin(PI*(2*y + 1)/2)/PI
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
