#include "confdefs.h"

subroutine test_vector_field_adaptivity

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
  use spud
  implicit none
#ifdef HAVE_MPI
  include "mpif.h"
#endif

  type(state_type) :: state(1)
  type(mesh_type), pointer :: mesh
  type(vector_field), pointer :: positions, vfield_pointer
  type(vector_field) :: vfield, vfield_weight
  type(scalar_field) :: edgelen
  type(tensor_field) :: metric

  integer :: i, stat
  real :: x, y, z

  current_debug_level = 0
  pseudo2d_coord = 3

  call vtk_read_state("data/pseudo2d.vtu", state(1))
  mesh => extract_mesh(state(1), "Mesh")
  call add_faces(mesh)
  positions => extract_vector_field(state(1), "Coordinate")
  ! Update mesh descriptor on positons
  positions%mesh=mesh

  call allocate(vfield, 3, mesh, "VField")
  call allocate(vfield_weight, 3, mesh, "VFieldInterpolationErrorBound", field_type=FIELD_TYPE_CONSTANT)
  call allocate(edgelen, mesh, "Edge lengths")
  call allocate(metric, mesh, "Metric")

  do i=1,mesh%nodes
    x = positions%val(1,i)
    y = positions%val(2,i)
    z = positions%val(3,i)

    vfield%val(1,i) = x
    vfield%val(2,i) = y
    vfield%val(3,i) = z
  end do

  vfield%option_path = "/whatever"
  call add_option("/whatever/virtual/adaptivity_options/absolute_measure", stat=stat)

  call set(vfield_weight, (/1.0, 1.0, 1.0/))
  
  call insert(state(1), vfield, "VField")
  call insert(state(1), vfield_weight, "VFieldInterpolationErrorBound")
  positions => extract_vector_field(state(1), "Coordinate")

  call adaptivity_bounds(state(1), 0.01, 1.0)

  call assemble_metric(state, metric)
  call get_edge_lengths(metric, edgelen)
  call vtk_write_fields("data/vfield_adapt", 0, positions, mesh, sfields=(/edgelen/), vfields=(/vfield, vfield_weight/), tfields=(/metric/))
  call adapt_state(state(1), metric)

  call report_test("[adaptivity runs]", .false., .false., "Congratulations! &
                   & You didn't crash.")

  mesh => extract_mesh(state(1), "Mesh")
  positions => extract_vector_field(state(1), "Coordinate")
  vfield_pointer => extract_vector_field(state(1), "VField")
  call vtk_write_fields("data/vfield_adapt", 1, positions,  mesh, vfields=(/vfield_pointer/)) 

  call report_test("[adaptivity output]", .false., .false., "Congratulations! &
                   & The output from adaptivity might even be OK if you get this far.")
end subroutine test_vector_field_adaptivity
