#include "fdebug.h"

module pseudo_supermesh
!!< Compute a pseudo-supermesh of a given set of meshes.
!!< A supermesh is a mesh such that if a node
!!< exists at x in any of the input meshes,
!!< a node exists at x in the supermesh.
!!< In a pseudo-supermesh, this requirement
!!< is relaxed; instead we demand that 
!!< local mesh density in a ball around 
!!< a given point x in the supermesh is
!!< the densest of all the densities 
!!< around x in the input meshes -- i.e.,
!!< we substitute guarantees about exact
!!< nodal placement for heuristic statements
!!< about local node density.

  use adapt_state_unittest_module, only : adapt_state => adapt_state_unittest
  use fields
  use interpolation_module
  use vtk_interfaces
  use conformity_measurement
  use merge_tensors
  use edge_length_module
  use limit_metric_module
  implicit none

  contains

  subroutine compute_pseudo_supermesh(snapshots, starting_positions, super_positions, no_its, mxnods)
    !!< snapshots is a list of VTUs containing the meshes we want
    !!< to merge.
    !!< starting_positions is the initial mesh + positions to interpolate the
    !!< metric tensor fields describing the snapshot meshes onto.
    !!< super_positions is the output -- a positions field on a mesh.
    character(len=255), dimension(:), intent(in) :: snapshots
    type(vector_field), intent(in) :: starting_positions
    type(vector_field), intent(out) :: super_positions
    integer, intent(in), optional :: no_its, mxnods

    integer :: lno_its
    integer :: it, i

    type(mesh_type) :: current_mesh, vtk_mesh
    type(vector_field) :: current_pos, vtk_pos
    type(state_type) :: vtk_state, temp_state
    type(state_type) :: interpolation_input, interpolation_output

    type(tensor_field) :: merged_metric, interpolated_metric
    type(tensor_field) :: vtk_metric

    if (present(no_its)) then
      lno_its = no_its
    else
      lno_its = 3
    end if

    current_pos = starting_positions
    call incref(starting_positions)
    current_mesh = starting_positions%mesh
    call incref(current_mesh)

    do it=1,lno_its
      call allocate(merged_metric, current_mesh, "MergedMetric")
      call zero(merged_metric)
      call allocate(interpolated_metric, current_mesh, "InterpolatedMetric")
      call zero(interpolated_metric)
      call insert(interpolation_output, interpolated_metric, "InterpolatedMetric")
      call insert(interpolation_output, current_mesh, "Mesh")
      call insert(interpolation_output, current_pos, "Coordinate")

!      call allocate(edgelen, current_mesh, "EdgeLengths")

      do i=1,size(snapshots)
        call zero(interpolated_metric)
        call vtk_read_state(trim(snapshots(i)), vtk_state)
        vtk_mesh = extract_mesh(vtk_state, "Mesh")
        vtk_pos  = extract_vector_field(vtk_state, "Coordinate")
        call compute_mesh_metric(vtk_pos, vtk_metric)
        call insert(interpolation_input, vtk_metric, "InterpolatedMetric")
        call insert(interpolation_input, vtk_mesh, "Mesh")
        call insert(interpolation_input, vtk_pos, "Coordinate")
        call vtk_write_state("interpolation_input", i, state=(/interpolation_input/))
        call linear_interpolation(interpolation_input, interpolation_output)
        call vtk_write_state("interpolation_output", i, state=(/interpolation_output/))
        call merge_tensor_fields(merged_metric, interpolated_metric)
        call deallocate(vtk_metric)
        call deallocate(interpolation_input)
        call deallocate(vtk_state)
      end do

      call deallocate(interpolated_metric)
      call deallocate(interpolation_output)

      call insert(temp_state, current_mesh, "Mesh")
      call insert(temp_state, current_pos, "Coordinate")
      ! Assuming current_mesh had a refcount of one,
      ! it now has a refcount of two.
!      call get_edge_lengths(merged_metric, edgelen)
!      call vtk_write_fields("supermesh_before_adapt", it, current_pos, current_mesh, sfields=(/edgelen/), tfields=(/merged_metric/))
      if (present(mxnods)) then
        call limit_metric(current_pos, merged_metric, min_nodes=1, max_nodes=mxnods)
      end if
      call adapt_state(temp_state, merged_metric)
!      call vtk_write_state("supermesh_after_adapt", it, state=(/temp_state/))
      call deallocate(merged_metric)
!      call deallocate(edgelen)
      ! Now it has a refcount of one, as adapt_state
      ! has destroyed the old one and created a new mesh
      ! with refcount one.

      ! We're finished with the current_mesh, so let it be
      ! deallocated if no one else is using it.
      call deallocate(current_mesh)
      call deallocate(current_pos)

      current_mesh = extract_mesh(temp_state, "Mesh")
      current_pos  = extract_vector_field(temp_state, "Coordinate")
      call incref(current_mesh)
      call incref(current_pos)
      call deallocate(temp_state)
    end do

    call deallocate(current_mesh)
    super_positions = current_pos
  end subroutine compute_pseudo_supermesh

end module pseudo_supermesh
