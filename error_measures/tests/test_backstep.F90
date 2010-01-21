#include "confdefs.h"

subroutine test_backstep

  use global_parameters, only: current_debug_level
  use unittest_tools
  use metric_assemble
  use edge_length_module
  use fields
  use state_module
  use vtk_interfaces
  use adapt_state_unittest_module, only : adapt_state => adapt_state_unittest
  use form_metric_field
#ifdef HAVE_MPI
  use mpi
#endif
  implicit none

  type(state_type) :: state
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field), pointer :: field
  type(tensor_field) :: metric
  type(metric_options) :: opts
  type(anisotropic_edge_options), target :: edge_opts
  real, dimension(4), target :: xhsamp, yhsamp, zhsamp
  real, dimension(4), target :: hminxx, hminxy, hminxz, hminyy, hminyz, hminzz
  real, dimension(4), target :: hmaxxx, hmaxxy, hmaxxz, hmaxyy, hmaxyz, hmaxzz
  integer :: i, nhsamp
  integer :: ierr

  call vtk_read_state("data/backstep.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  field => extract_scalar_field(state, "Temperature")
  call allocate(metric, mesh, "Metric")

  opts%min_edge_length = 0.01
  opts%max_edge_length = 1.0
  opts%use_anisotropic_edge_length = .true.
  xhsamp(1) = 0.0; yhsamp(1) = 0.0; zhsamp(1) = 0.0
  xhsamp(2) = 1.0; yhsamp(2) = 0.0; zhsamp(2) = 1.0
  xhsamp(3) = 1.0; yhsamp(3) = 1.0; zhsamp(3) = 0.0
  xhsamp(4) = 0.0; yhsamp(4) = 1.0; zhsamp(4) = 1.0
  hminxx = 0.01; hminxy = 0.0; hminxz = 0.0; hminyy = 0.2; hminyz = 0.0; hminzz = 0.01
  hmaxxx = 7.5; hmaxxy = 0.0; hmaxxz = 0.0; hmaxyy = 2.5; hmaxyz = 0.0; hmaxzz = 2.0
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

  call dprintf(-1, "%p\n$", state%scalar_fields(2))
  call dprintf(-1, "%p\n$", state%scalar_fields(2)%val)
  call assemble_metric((/state/), metric, opts)
  call vtk_write_fields("data/backstep", 0, positions, mesh, sfields=(/field/), tfields=(/metric/))
  call adapt_state(state, metric)

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  field => extract_scalar_field(state, "Temperature")
  call vtk_write_fields("data/backstep", 1, positions, mesh, sfields=(/field/)) 

  call deallocate(metric)
  call deallocate(state)

end subroutine test_backstep
