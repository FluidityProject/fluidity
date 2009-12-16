#include "fdebug.h"

module geometric_constraints_metric
!!< This module wraps Gerard's geometric constraints
!!< code and applied it during metric formation.

  use vtk_interfaces
  use metric_tools
  use merge_tensors
  use edge_length_module
  use fields
  use node_boundary
  use form_metric_field
  use gradation_metric
  use parallel_tools
  use halos
  use state_module
  implicit none

  logical :: use_geometric_constraints_metric = .false.
  logical :: geometric_constraints_initialised = .false.

  contains

  subroutine initialise_geometric_constraints_metric
    use_geometric_constraints_metric = .false.
    if (have_option("/mesh_adaptivity/hr_adaptivity/geometric_constraints")) then
      use_geometric_constraints_metric = .true.
    end if
    geometric_constraints_initialised = .true.
  end subroutine initialise_geometric_constraints_metric

  subroutine form_geometric_constraints_metric(positions, error_metric, state)
    type(tensor_field), intent(inout) :: error_metric !!< The metric formed so far
    type(vector_field), intent(in) :: positions
    type(state_type), intent(in) :: state

    integer :: dim
    integer :: stat, stat2
    real :: gradation_parameter

    real, dimension(error_metric%dim * error_metric%dim * node_count(error_metric)) :: geometric_edge_lengths_raw
    type(tensor_field) :: geometric_edge_lengths
    integer :: snloc, nselements
    integer :: noits, grad_count

    type(scalar_field) :: edgelen, geometric_edgelen
    integer, save :: adaptcnt = 0
#ifdef HAVE_MPI
    include 'mpif.h'
    integer::noits_max, ierr
#endif

    integer, dimension(:), pointer :: lsenlist
    logical :: debug_metric
    
    if(.not.use_geometric_constraints_metric) then
        return
    end if

    debug_metric = have_option("/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages")

    snloc = face_loc(error_metric%mesh, 1)
    nselements = surface_element_count(error_metric%mesh)
    allocate(lsenlist(1:nselements*snloc))
    call getsndgln(error_metric%mesh, lsenlist)
    
    dim = error_metric%dim
    ewrite(2,*) "++: Applying geometric constraints"

    call FindGeometryConstraints(positions, geometric_edge_lengths_raw)
    
    deallocate(lsenlist)

    geometric_edge_lengths = wrap_tensor_field(error_metric%mesh, geometric_edge_lengths_raw, "GeometricEdgeLengths")
    call get_option("/mesh_adaptivity/hr_adaptivity/enable_gradation/gradation_parameter", gradation_parameter, stat=stat)
    call set_option("/mesh_adaptivity/hr_adaptivity/enable_gradation/gradation_parameter", 1.1, stat=stat2)

    call bound_metric(geometric_edge_lengths, state)
    if (.not. isparallel()) then
      call form_gradation_metric(positions, geometric_edge_lengths)
    else
      noits = 2
      grad_count = 0
      do while ((noits.gt.1).and.(grad_count.le.3))
         call form_gradation_metric(positions, geometric_edge_lengths, noits)
         call halo_update(error_metric)
         grad_count = grad_count + 1
#ifdef HAVE_MPI
         CALL MPI_Allreduce (noits, noits_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
         noits = noits_max
#endif
      end do
    end if

    if (stat == 0) then
      call set_option("/mesh_adaptivity/hr_adaptivity/enable_gradation/gradation_parameter", gradation_parameter)
    else
      call delete_option("/mesh_adaptivity/hr_adaptivity/enable_gradation/gradation_parameter")
    end if

    if (debug_metric) then
      call allocate(edgelen, error_metric%mesh, "Desired edge lengths")
      call allocate(geometric_edgelen, error_metric%mesh, "Geometric edge lengths")
      call get_edge_lengths(geometric_edge_lengths, geometric_edgelen)
      call get_edge_lengths(error_metric, edgelen)
      call vtk_write_fields(trim("geometric_constraints_metric"), adaptcnt, positions, positions%mesh, &
                             sfields=(/edgelen, geometric_edgelen/), tfields=(/error_metric, geometric_edge_lengths/))
      call deallocate(edgelen)
      call deallocate(geometric_edgelen)
      adaptcnt = adaptcnt + 1
    end if
    call merge_tensor_fields(error_metric, geometric_edge_lengths)

    call deallocate(geometric_edge_lengths)
  end subroutine form_geometric_constraints_metric

end module geometric_constraints_metric
