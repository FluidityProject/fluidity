subroutine test_gradation

  use vtk_interfaces
  use metric_assemble
  use unittest_tools
  use adapt_state_unittest_module, only : adapt_state => adapt_state_unittest
  use state_module
  use mpi
  use edge_length_module
  use gradation_metric
  use field_derivatives
  use form_metric_field
  use fields
  use global_parameters, only: pseudo2d_coord
  use surfacelabels
  implicit none

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field) :: field
  type(scalar_field), pointer :: field_ptr
  type(tensor_field) :: metric
  type(metric_options) :: opts
  type(anisotropic_edge_options), target :: edge_opts
  real, dimension(4), target :: xhsamp, yhsamp, zhsamp
  real, dimension(4), target :: hminxx, hminxy, hminxz, hminyy, hminyz, hminzz
  real, dimension(4), target :: hmaxxx, hmaxxy, hmaxxz, hmaxyy, hmaxyz, hmaxzz
  integer :: nhsamp
  logical :: fail
  integer :: i

  interface
    function func(pos) result(val)
      real, dimension(:), intent(in) :: pos
      real :: val
    end function func
  end interface

  pseudo2d_coord = 3

  fail = .false.

  call vtk_read_state("data/gradation_input.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  
  call allocate(field, mesh, "Field")
  call allocate(metric, mesh, "Metric")
  call set_from_function(field, func, positions)

  field%options%relative = .false.
  field%options%error = 0.01
  field%options%min_psi = 0.1

  opts%min_edge_length = 0.0001
  opts%max_edge_length = 2.0
  opts%use_anisotropic_edge_length = .true.
  xhsamp(1) = 0.0; yhsamp(1) = 0.0; zhsamp(1) = 0.0
  xhsamp(2) = 1.0; yhsamp(2) = 0.0; zhsamp(2) = 1.0
  xhsamp(3) = 1.0; yhsamp(3) = 1.0; zhsamp(3) = 0.0
  xhsamp(4) = 0.0; yhsamp(4) = 1.0; zhsamp(4) = 1.0
  hminxx = 0.0001; hminxy = 0.0; hminxz = 0.0; hminyy = 0.0001; hminyz = 0.0; hminzz = 0.01
  hmaxxx = 2.0; hmaxxy = 0.0; hmaxxz = 0.0; hmaxyy = 2.0; hmaxyz = 0.0; hmaxzz = 0.02
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

  call insert(state, field, "Field")

  use_gradation_metric = .true. ; gamma0 = 1.1; gradation_initialised = .true.
  call assemble_metric((/state/), metric, opts)
  call vtk_write_fields("data/gradation_adapt", 0, positions, mesh, &
                        sfields=(/field/))
  call adapt_state(state, metric)

  do i=1,6
    mesh => extract_mesh(state, "Mesh")
    positions => extract_vector_field(state, "Coordinate")
    field_ptr => extract_scalar_field(state, "Field")
    field_ptr%options%relative = .false.
    field_ptr%options%error = 0.01
    field_ptr%options%min_psi = 0.1

    call set_from_function(field_ptr, func, positions)
    call deallocate(metric); call allocate(metric, mesh, "Metric")
    call assemble_metric((/state/), metric, opts)
    call adapt_state(state, metric)
  end do

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  field_ptr => extract_scalar_field(state, "Field")
  call vtk_write_fields("data/gradation_adapt", 1, positions,  mesh, sfields=(/field_ptr/))

  call report_test("[gradation]", .false., .false., "Gradation runs. That's enough for now.")
end subroutine test_gradation

function func(pos) result(val)
  real, dimension(:), intent(in) :: pos
  real :: val

  val = tanh(25.0 * (pos(1) + 1))
end function func
