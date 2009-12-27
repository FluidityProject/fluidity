#include "confdefs.h"

subroutine test_adaptivity

  use global_parameters, only: current_debug_level
  use node_boundary, only: pseudo2d_coord
  use unittest_tools
  use metric_assemble
  use edge_length_module
  use fields
  use state_module
  use vtk_interfaces
  use adapt_state_unittest_module, only : adapt_state => adapt_state_unittest
  use form_metric_field
  use field_options
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

  integer :: i
  real :: x, y, z

  current_debug_level = 0
  pseudo2d_coord = 3

  call vtk_read_state("data/pseudo2d.vtu", state)
  mesh => extract_mesh(state, "Mesh")
  call add_faces(mesh)
  positions => extract_vector_field(state, "Coordinate")
  ! Update mesh descriptor on positons
  positions%mesh=mesh

  call allocate(pressure, mesh, "Pressure")
  call allocate(velocity, 3, mesh, "Velocity")
  
  call allocate(metric, mesh, "Metric")

  do i=1,mesh%nodes
    x = positions%val(1)%ptr(i)
    y = positions%val(2)%ptr(i)
    z = positions%val(3)%ptr(i)
!    if (x == 0.0) pressure%val(i) = 0.0
!    if (x <= 3.0) pressure%val(i) = x   * 1.0e6  ! ramp up .. rapidly
!    if (x >  3.0) pressure%val(i) = 3.0 * 1.0e6
    pressure%val(i) = x * x
    velocity%val(1)%ptr(i) = x
    velocity%val(2)%ptr(i) = y
    velocity%val(3)%ptr(i) = z
  end do
  
  call adaptivity_options(state, pressure, 1.0, .false.)

  call insert(state, pressure, "Pressure")
  call insert(state, velocity, "Velocity")
  call deallocate(pressure)
  call deallocate(velocity)
  
  positions => extract_vector_field(state, "Coordinate")

  call adaptivity_bounds(state, 0.01, 1.0)

  state_array(1) = state
  call assemble_metric(state_array, metric)
  
  call allocate(edgelen, mesh, "Edge lengths")
  call get_edge_lengths(metric, edgelen)
  call vtk_write_fields("data/adapt", 0, positions, mesh, sfields=(/edgelen/), vfields=(/velocity/), tfields=(/metric/))
  call deallocate(edgelen)
  
  call adapt_state(state, metric)

  call report_test("[adaptivity runs]", .false., .false., "Congratulations! &
                   & You didn't crash.")

  mesh => extract_mesh(state, "Mesh")
  positions => extract_vector_field(state, "Coordinate")
  velocity_pointer => extract_vector_field(state, "Velocity")
  call vtk_write_fields("data/adapt", 1, positions,  mesh, vfields=(/velocity_pointer/)) 

  call report_test("[adaptivity output]", .false., .false., "Congratulations! &
                   & The output from adaptivity might even be OK if you get this far.")

  call deallocate(metric)
  call deallocate(state)
  
  call report_test_no_references()

end subroutine test_adaptivity
