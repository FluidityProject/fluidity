subroutine test_boundary_layer_adaptivity_two

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
  type(vector_field) :: gradient
  type(scalar_field), pointer :: ptr_field
  type(tensor_field) :: metric, hessian
  type(metric_options) :: opts
  type(anisotropic_edge_options), target :: edge_opts
  real, dimension(4), target :: xhsamp, yhsamp, zhsamp
  real, dimension(4), target :: hminxx, hminxy, hminxz, hminyy, hminyz, hminzz
  real, dimension(4), target :: hmaxxx, hmaxxy, hmaxxz, hmaxyy, hmaxyz, hmaxzz
  logical :: fail = .false., warn = .false.
  integer :: i, nhsamp
  real :: x, y, z 
  integer :: ierr

  pseudo2d_coord = 3
  call vtk_read_state("data/1x1square-delaunay.vtu", state)

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  call allocate(field, mesh, "Field") 

  do i=1,mesh%nodes
    x = positions%val(1)%ptr(i)
    y = positions%val(2)%ptr(i)
    z = positions%val(3)%ptr(i)
    field%val(i) = tanh(50 * (x-0.5))
  end do
  field%options%relative = .false.
  field%options%error = 0.01
  field%options%min_psi = 1e-10
  field%options%square_eigs = .false.

  call insert(state, field, "Field")
  call allocate(metric, mesh, "Metric")

  opts%min_edge_length = 0.001
  opts%max_edge_length = 1.0
  opts%use_anisotropic_edge_length = .false.
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

  gamma0 = 2.0
  use_gradation_metric = .false.
  gradation_initialised = .true.
  call assemble_metric((/state/), metric, opts)
  call vtk_write_fields("data/boundary_layer_adapt_two", 0, positions, mesh, &
                        sfields=(/field/), tfields=(/metric/))
  call adapt_state(state, metric)

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  ptr_field => extract_scalar_field(state, "Field")

  call vtk_write_fields("data/boundary_layer_adapt_two", 1, positions,  mesh, sfields=(/ptr_field/)) 
end subroutine test_boundary_layer_adaptivity_two
