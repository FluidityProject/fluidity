#include "confdefs.h"

subroutine test_mba_adapt

  use unittest_tools
  use metric_assemble
  use edge_length_module
  use fields
  use state_module
  use vtk_interfaces
  use mba_adapt_module
  use form_metric_field
  use field_options
  use spud
  implicit none
#ifdef HAVE_MPI
  include "mpif.h"
#endif

  type(state_type) :: state, state_array(1)
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions, velocity_pointer
  type(vector_field) :: velocity
  type(scalar_field) :: pressure, edgelen
  type(tensor_field) :: metric

  integer :: i, stat
  real :: x, y

#ifdef HAVE_MBA

  call vtk_read_state("data/mms_a.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  call add_faces(mesh)
  positions => extract_vector_field(state, "Coordinate")
  call allocate(pressure, mesh, "Pressure")
  call allocate(velocity, 2, mesh, "Velocity")
  call allocate(edgelen, mesh, "Edge lengths")
  call allocate(metric, mesh, "Metric")

  do i=1,mesh%nodes
    x = positions%val(1)%ptr(i)
    y = positions%val(2)%ptr(i)
    pressure%val(i) = x * x
    velocity%val(1)%ptr(i) = x
    velocity%val(2)%ptr(i) = y
  end do
  
  call adaptivity_options(state, pressure, 0.001, .false.)

  call insert(state, pressure, "Pressure")
  call insert(state, velocity, "Velocity")
  positions => extract_vector_field(state, "Coordinate")

  call set_option("/mesh_adaptivity/hr_adaptivity/constant_size_constraint/minimum_edge_length", 0.01, stat=stat)
  call set_option("/mesh_adaptivity/hr_adaptivity/constant_size_constraint/maximum_edge_length", 1.0, stat=stat)

  state_array(1) = state
  !call assemble_metric(state_array, metric)
  metric%val = 0.0
  metric%val(1, 1, :) = 1.0; metric%val(2, 2, :) = 1.0
  call get_edge_lengths(metric, edgelen)
  call vtk_write_state("data/2d_adapt", 0, state=(/state/))
  call mba_adapt(state, metric)

  call report_test("[adaptivity runs]", .false., .false., "Congratulations! &
                   & You didn't crash.")

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  velocity_pointer => extract_vector_field(state, "Velocity")
  call vtk_write_state("data/2d_adapt", 1, state=(/state/))
#endif

  call report_test("[adaptivity output]", .false., .false., "Congratulations! &
                   & The output from adaptivity might even be OK if you get this far.")

#ifdef HAVE_MBA
  call deallocate(metric)
  call deallocate(state)
#endif

end subroutine test_mba_adapt
