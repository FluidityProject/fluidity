#include "confdefs.h"

subroutine compute_driven_cavity_adapt

  use unittest_tools
  use metric_assemble
  use fields
  use state_module
  use vtk_interfaces
  use adapt_state_unittest_module, only : adapt_state => adapt_state_unittest
  use form_metric_field
  use field_options
  implicit none

  type(state_type) :: state(1)
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(vector_field), pointer :: velocity
  type(tensor_field) :: metric
  type(metric_options) :: opts
  type(anisotropic_edge_options), target :: edge_opts
  real, dimension(4), target :: xhsamp, yhsamp, zhsamp
  real, dimension(4), target :: hminxx, hminxy, hminxz, hminyy, hminyz, hminzz
  real, dimension(4), target :: hmaxxx, hmaxxy, hmaxxz, hmaxyy, hmaxyz, hmaxzz
  integer :: nhsamp

  call vtk_read_state("data/driven_cavity_58.vtu", state(1))
  mesh => extract_mesh(state(1), "Mesh")
  positions => extract_vector_field(state(1), "Coordinate")
  velocity => extract_vector_field(state(1), "Velocity")

  call adaptivity_options(state(1), velocity, (/0.00125, 0.00125, 0.0/), .false.)
  call allocate(metric, mesh, "ErrorMetric")

  opts%min_edge_length = 0.01
  opts%max_edge_length = 1.0
  opts%use_anisotropic_edge_length = .true.
  xhsamp(1) = 0.0; yhsamp(1) = 0.0; zhsamp(1) = 0.0
  xhsamp(2) = 1.0; yhsamp(2) = 0.0; zhsamp(2) = 1.0
  xhsamp(3) = 1.0; yhsamp(3) = 1.0; zhsamp(3) = 0.0
  xhsamp(4) = 0.0; yhsamp(4) = 1.0; zhsamp(4) = 1.0
  hminxx = 0.002; hminxy = 0.0; hminxz = 0.0; hminyy = 0.002; hminyz = 0.0; hminzz = 0.005
  hmaxxx = 0.1; hmaxxy = 0.0; hmaxxz = 0.0; hmaxyy = 0.1; hmaxyz = 0.0; hmaxzz = 0.02
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

  call assemble_metric(state, metric, opts)
  call vtk_write_state("data/driven_cavity_adapt", 0, state=state)
  call adapt_state(state(1), metric)
  call vtk_write_state("data/driven_cavity_adapt", 1, state=state)

  call deallocate(metric)
  call deallocate(state(1))

end subroutine compute_driven_cavity_adapt
