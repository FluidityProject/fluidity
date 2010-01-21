subroutine test_anisotropically_bounded_metric

  use metric_assemble
  use form_metric_field
  use state_module
  use vtk_interfaces
  use vector_tools
  use unittest_tools
  use edge_length_module
  use mpi
  implicit none

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: position_field
  type(scalar_field), pointer :: pressure_field
  type(tensor_field) :: metric
  type(metric_options) :: opts
  type(anisotropic_edge_options), target :: edge_opts
  real, dimension(4), target :: xhsamp, yhsamp, zhsamp
  real, dimension(4), target :: hminxx, hminxy, hminxz, hminyy, hminyz, hminzz
  real, dimension(4), target :: hmaxxx, hmaxxy, hmaxxz, hmaxyy, hmaxyz, hmaxzz
  logical :: fail = .false., warn = .false.
  integer :: i, nhsamp
  real :: x, y, z 
  real :: max_eigenbound, min_eigenbound
  integer :: ierr

  call MPI_init(ierr)

  call vtk_read_state("data/test_spr.vtu", state)

  mesh => extract_mesh(state, "Mesh")
  position_field => extract_vector_field(state, "Coordinate")
  pressure_field => extract_scalar_field(state, "Pressure")

  do i=1,mesh%nodes
    x = position_field%val(1)%ptr(i)
    y = position_field%val(2)%ptr(i)
    z = position_field%val(3)%ptr(i)
    pressure_field%val(i) = 10.0 * x * x + 0.5 * y * y
  end do

  call insert(state, pressure_field, "Pressure")

  call allocate(metric, mesh, "Metric")

  opts%min_edge_length = 0.1
  opts%max_edge_length = 1.0
  opts%use_anisotropic_edge_length = .true.
  xhsamp(1) = 0.0; yhsamp(1) = 0.0; zhsamp(1) = 0.0
  xhsamp(2) = 1.0; yhsamp(2) = 0.0; zhsamp(2) = 1.0
  xhsamp(3) = 1.0; yhsamp(3) = 1.0; zhsamp(3) = 0.0
  xhsamp(4) = 0.0; yhsamp(4) = 1.0; zhsamp(4) = 1.0
  hminxx = 20.0; hminxy = 0.0; hminxz = 0.0; hminyy = 50.0; hminyz = 0.0; hminzz = 10.0
  hmaxxx = 5000.0; hmaxxy = 0.0; hmaxxz = 0.0; hmaxyy = 150.0; hmaxyz = 0.0; hmaxzz = 500.0
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

  ! the corresponding eigenbounds are:
  max_eigenbound = 1.0
  min_eigenbound = 0.25

  metric%val = 0.0
  call bound_metric(metric, pressure_field, position_field, opts)

  fail = .false.
  do i=1,mesh%nodes
    if (metric%val(:, :, i) .fne. get_mat_diag(eigenvalue_from_edge_length((/hmaxxx(1), hmaxyy(1), hmaxzz(1)/)))) fail = .true.
  end do


  call report_test("[anisotropic bounds deal with zero]", fail, warn, "Zero hessian should give max edge length bounds.")

  call MPI_finalize(ierr)
end subroutine test_anisotropically_bounded_metric
