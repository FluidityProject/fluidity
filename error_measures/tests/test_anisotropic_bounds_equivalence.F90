subroutine test_anisotropic_bounds_equivalence

  use global_parameters, only: current_debug_level, pseudo2d_coord
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
  implicit none

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field) :: field
  type(tensor_field) :: metric_iso, metric_aniso
  type(metric_options) :: opts
  type(anisotropic_edge_options), target :: edge_opts
  real, dimension(4), target :: xhsamp, yhsamp, zhsamp
  real, dimension(4), target :: hminxx, hminxy, hminxz, hminyy, hminyz, hminzz
  real, dimension(4), target :: hmaxxx, hmaxxy, hmaxxz, hmaxyy, hmaxyz, hmaxzz
  logical :: fail
  integer :: i, nhsamp
  real :: x, y, z 

  pseudo2d_coord = 3
  call vtk_read_state("data/1x1square-delaunay.vtu", state)

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  call allocate(field, mesh, "Field") 

  do i=1,mesh%nodes
    x = positions%val(1)%ptr(i)
    y = positions%val(2)%ptr(i)
    z = positions%val(3)%ptr(i)
    field%val(i) = (y-0.5)**3
  end do
  call insert(state, field, "Field")
  call allocate(metric_iso, mesh, "Metric")
  call allocate(metric_aniso, mesh, "Metric")

  opts%min_edge_length = 0.001
  opts%max_edge_length = 1.0
  opts%use_anisotropic_edge_length = .false.
  xhsamp(1) = 0.0; yhsamp(1) = 0.0; zhsamp(1) = 0.0
  xhsamp(2) = 1.0; yhsamp(2) = 0.0; zhsamp(2) = 1.0
  xhsamp(3) = 1.0; yhsamp(3) = 1.0; zhsamp(3) = 0.0
  xhsamp(4) = 0.0; yhsamp(4) = 1.0; zhsamp(4) = 1.0
  hminxx = 0.001; hminxy = 0.0; hminxz = 0.0; hminyy = 0.001; hminyz = 0.0; hminzz = 0.001
  hmaxxx = 1.0; hmaxxy = 0.0; hmaxxz = 0.0; hmaxyy = 1.0; hmaxyz = 0.0; hmaxzz = 1.0
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

  gamma0 = 2.0
  use_gradation_metric = .false.
  gradation_initialised = .true.
  call assemble_metric((/state/), metric_iso, opts)
  opts%use_anisotropic_edge_length = .true.
  call assemble_metric((/state/), metric_aniso, opts)

  fail = .false.
  do i=1,mesh%nodes
    if (metric_iso%val(:, :, i) .fne. metric_aniso%val(:, :, i)) then
      fail = .true.
    end if
  end do

  call report_test("[anisotropic bounds equivalence]", fail, .false., &
  & "Forming the metric with isotropic and anisotropic bounds &
  & should give the same results in this case.")
end subroutine test_anisotropic_bounds_equivalence
