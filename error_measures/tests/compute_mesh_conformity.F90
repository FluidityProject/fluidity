#include "confdefs.h"

subroutine compute_mesh_conformity

  use global_parameters, only: current_debug_level, pseudo2d_coord
  use unittest_tools
  use metric_assemble
  use edge_length_module
  use fields
  use state_module
  use vtk_interfaces
  use adapt_state_module
  use form_metric_field
  use conformity_measurement
  use field_options
#ifdef HAVE_MPI
  use mpi
#endif
  implicit none

  type(state_type) :: state, state_array(1)
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions
  type(scalar_field) :: pressure, edgelen
  type(tensor_field) :: metric
  type(metric_options) :: opts

  integer :: i
  real :: x, y, z

  current_debug_level = 0
  pseudo2d_coord = 3

  call vtk_read_state("data/pseudo2d.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  call allocate(pressure, mesh, "Pressure")
  call allocate(edgelen, mesh, "Edge lengths")
  call allocate(metric, mesh, "Metric")

  do i=1,mesh%nodes
    x = positions%val(1)%ptr(i)
    y = positions%val(2)%ptr(i)
    z = positions%val(3)%ptr(i)
!    if (x == 0.0) pressure%val(i) = 0.0
!    if (x <= 3.0) pressure%val(i) = x   * 1.0e6  ! ramp up .. rapidly
!    if (x >  3.0) pressure%val(i) = 3.0 * 1.0e6
    pressure%val(i) = x * x
  end do
  
  call adaptivity_options(state, pressure, 0.05, .false.)

  call insert(state, pressure, "Pressure")
  positions => extract_vector_field(state, "Coordinate")

  opts%min_edge_length = 0.01
  opts%max_edge_length = 1.0
  opts%use_anisotropic_edge_length = .false.

  state_array(1) = state
  call assemble_metric(state_array, metric, opts)
  call insert_mesh_conformity(state_array, metric)
  call remove_scalar_field(state_array(1), "Pressure")
  call remove_scalar_field(state_array(1), "PressureInterpolationErrorBound")
  call vtk_write_state("data/mesh_conformity", 0, state=state_array)


  call deallocate(metric)
  call deallocate(state)

end subroutine compute_mesh_conformity
