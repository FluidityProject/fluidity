#include "fdebug.h"

module metric_assemble

  use VTK_interfaces
  use surfacelabels
!  use fields
  use field_derivatives
  use form_metric_field
  use merge_tensors
  use edge_length_module
  use interpolation_metric
  use goals
  use goal_metric
  use convective_metric
  use gradation_metric
  use bounding_box_metric
  use boundary_metric
  use geometric_constraints_metric
  use limit_metric_module
  use metric_tools
  use state_module
  use halos
  use spud
  use parallel_tools
  use metric_advection
  use anisotropic_gradation
  use richardson_metric_module
  use anisotropic_zz_module
  use reference_meshes
  
  implicit none

  private
  public assemble_metric
  
  contains

  subroutine assemble_metric(state, error_metric)
    !!< This routine drives the metric assembly logic, which is contained in
    !!< other modules.

    type(state_type), dimension(:), intent(inout) :: state
    type(tensor_field), intent(inout) :: error_metric
    integer :: noits, grad_count

    integer :: i, stat

    type(vector_field), pointer :: positions
    integer, save :: adaptcnt = 0
    character(len=20) :: buf
    logical :: debug_metric, vertically_structured_adaptivity
    type(tensor_field), pointer :: max_tensor

#ifdef HAVE_MPI
    include 'mpif.h'
    integer::noits_max, ierr
#endif

    ewrite(2,*) "+: Assembling metric"
    debug_metric = have_option("/mesh_adaptivity/hr_adaptivity/debug/write_metric_stages")
    ! is this metric going to be collapsed in the vertical to do horizontal adaptivity with it?
    vertically_structured_adaptivity = have_option("/mesh_adaptivity/hr_adaptivity/vertically_structured_adaptivity")

    do i=1,size(state)
      positions => extract_vector_field(state(i), "Coordinate", stat=stat)
      if (stat == 0) exit
    end do

    max_tensor => extract_tensor_field(state(1), "MinMetricEigenbound")

    if (debug_metric) then
      do i=1,size(state)
        write(buf, '(i0)') i 

        call vtk_write_state(trim("metric_input_") // trim(buf),&
             & adaptcnt, state=(/state(i)/)) 
      end do
      adaptcnt = adaptcnt + 1
    end if

    call zero(error_metric)

    call initialise_boundcount(error_metric%mesh, positions)

    call initialise_interpolation_metric
    call initialise_goal_metric
    call initialise_convective_metric
    call initialise_gradation_metric
    call initialise_anisotropic_gradation
    call initialise_bounding_box_metric(positions)
    call initialise_boundary_metric
    call initialise_geometric_constraints_metric
    call initialise_metric_advection
    call initialise_richardson_number_metric

    if (use_goal_metric) then
      call form_goal_metric(state, error_metric)
      call halo_update(error_metric)
    else if (use_interpolation_metric) then
      call form_interpolation_metric(state, error_metric)
      call halo_update(error_metric)
      call form_anisotropic_zz_metric(state, error_metric)
    end if
    
    if (use_richardson_number_metric) then
      call form_richardson_number_metric(state, error_metric)
      call halo_update(error_metric)
    end if 
    
    if (use_boundary_metric) then
      call form_boundary_metric(error_metric, positions)
      call halo_update(error_metric)
    end if

    if (use_convective_metric) then
       call form_convective_metric(state(1), error_metric)
       call halo_update(error_metric)
    end if
    
    call enforce_reference_meshes(state, positions, error_metric)

    if (use_geometric_constraints_metric) then
      call form_geometric_constraints_metric(positions, error_metric, state(1))
      call halo_update(error_metric)
    end if

    if (use_anisotropic_gradation) then
      call form_anisotropic_gradation_metric(error_metric, positions, state(1))
    else if (use_gradation_metric) then
      if (.not. isparallel()) then
        call form_gradation_metric(positions, error_metric)
      else
        noits = 2
        grad_count = 0
        do while ((noits.gt.1).and.(grad_count.le.3))
           call form_gradation_metric(positions, error_metric, noits)
           call halo_update(error_metric)
           grad_count = grad_count + 1
#ifdef HAVE_MPI
           CALL MPI_Allreduce (noits, noits_max, 1, getpinteger(), MPI_MAX, MPI_COMM_WORLD, ierr)
           assert(ierr == MPI_SUCCESS)
           noits = noits_max
#endif
        end do
      end if
    end if

    if (use_metric_advection) then
      call form_advection_metric(error_metric, state(1))
      call halo_update(error_metric)
    end if

    if (use_bounding_box_metric) then
      call form_bounding_box_metric(positions, error_metric, max_tensor)
      call halo_update(error_metric)
    end if

    ! for vertically structured, the limiting should happen after the vertical collapsing
    if (.not. vertically_structured_adaptivity) then
      call limit_metric(positions, error_metric)
    end if
    call halo_update(error_metric)
    
  end subroutine assemble_metric

end module metric_assemble
