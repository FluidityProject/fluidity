#include "confdefs.h"
#include "fdebug.h"

module zoltan_detectors

#ifdef HAVE_ZOLTAN

  use zoltan, only: zoltan_int
  use zoltan_global_variables, only: zoltan_global_uen_to_new_local_numbering, &
                                     zoltan_global_old_local_numbering_to_uen, &
                                     zoltan_global_new_positions
  use data_structures, only: has_key, fetch
  use detector_data_types
  use detector_tools
  use detector_parallel
  use fields

  implicit none

  public :: zoltan_detectors_migrate

  private

contains

  subroutine zoltan_detectors_migrate()
    type(detector_list_ptr), dimension(:), pointer :: detector_lists
    type(detector_linked_list), dimension(:), allocatable :: send_lists, recv_lists
    type(vector_field), pointer :: xfield
    integer :: l, nprocs, dim

    if (get_num_detector_lists() <= 0) return

    nprocs = getnprocs()
    allocate( send_lists(nprocs) )
    allocate( recv_lists(nprocs) )
    xfield => zoltan_global_new_positions
    dim = xfield%dim

    ! Migrate detectors from all lists
    call get_registered_detector_lists(detector_lists)
    do l=1, get_num_detector_lists()
       
       call zoltan_detectors_prepare_migration( detector_lists(l)%ptr, send_lists )

       call detectors_send_all_to_all(detector_lists(l)%ptr, send_lists, recv_lists, dim, tracking=.false.)

       call zoltan_detectors_finalise_migration( detector_lists(l)%ptr, recv_lists )

       call distribute_detectors( detector_lists(l)%ptr, xfield )
    end do
    
  end subroutine zoltan_detectors_migrate

  subroutine zoltan_detectors_prepare_migration(detector_list, send_lists)
    type(detector_linked_list), pointer, intent(inout) :: detector_list
    type(detector_linked_list), dimension(:), intent(inout) :: send_lists

    type(detector_type), pointer :: detector, detector_to_move
    integer :: j, nprocs, local_ele, new_owner, detector_uen

    ewrite(2,*) "Preparing detector migration for list: "//trim(detector_list%name)
    ewrite(2,*) "Found "//int2str( detector_list%length )//" local and "//int2str( detector_list%total_num_det )//" global detectors"
    nprocs = getnprocs()

    detector => detector_list%first
    do while(associated(detector))

       ! Translate detector element to UEN
       if (.not. has_key(zoltan_global_old_local_numbering_to_uen, detector%element)) then
          ewrite(-1,*) "No uen found in Zoltan for detector ", detector%id_number, " with local element number ", detector%element
          FLAbort("No universal element number for detector found before Zoltan migration")
       end if
       detector%element = fetch(zoltan_global_old_local_numbering_to_uen, detector%element)

       if (has_key(zoltan_global_uen_to_new_local_numbering, detector%element)) then
          local_ele = fetch(zoltan_global_uen_to_new_local_numbering, detector%element)
          new_owner = element_owner(zoltan_global_new_positions%mesh, local_ele)

          if (new_owner == getprocno()) then
             ! Keep this detector
             detector => detector%next
          else
             ! Move detector to send list
             detector_to_move => detector
             detector => detector%next
             call move(detector_to_move, detector_list, send_lists(new_owner))
             detector_to_move => null()
          end if

       else
          ! Can't determine destination, so flag detector for later broadcast
          detector%element = -1
          detector => detector%next
       end if
    end do

  end subroutine zoltan_detectors_prepare_migration

  subroutine zoltan_detectors_finalise_migration(detector_list, recv_lists)
    type(detector_linked_list), pointer, intent(inout) :: detector_list
    type(detector_linked_list), dimension(:), intent(inout) :: recv_lists

    type(detector_type), pointer :: detector, detector_to_move
    integer :: j, p, dim

    dim = zoltan_global_new_positions%dim

    ! First update the detectors we didn't send
    detector => detector_list%first
    do while(associated(detector))

       if (detector%element > 0) then

          ! Update the element number for the detector
          if(has_key(zoltan_global_uen_to_new_local_numbering, detector%element)) then
             detector%element = fetch(zoltan_global_uen_to_new_local_numbering, detector%element)
          else
             ewrite(-1,*) "No new local element found by Zoltan for retained detector ", detector%id_number, " with uen ", detector%element
             FLAbort("No local element number for detector found after Zoltan migration")
          end if
       end if

       detector => detector%next       
    end do

    do p=1, getnprocs()
       detector => recv_lists(p)%first
       do while(associated(detector))

          ! Update the element number for the detector
          if(has_key(zoltan_global_uen_to_new_local_numbering, detector%element)) then
             detector%element = fetch(zoltan_global_uen_to_new_local_numbering, detector%element)
          else
             ewrite(-1,*) "No local element number found in Zoltan for detector ", detector%id_number, " with uen ", detector%element
             FLAbort("No local element number for detector found after Zoltan migration")
          end if

          ! Establish local coordinates
          if (.not. allocated(detector%local_coords)) allocate( detector%local_coords(dim+1) )
          detector%local_coords = local_coords(zoltan_global_new_positions, detector%element, detector%position)

          ! If there is a list of detector names, use it, otherwise set ID as name
          if (allocated(detector_list%detector_names)) then
             detector%name=detector_list%detector_names(detector%id_number)
          else
             detector%name=int2str(detector%id_number)
          end if

          ! Move received detector back into parent list
          detector_to_move => detector
          detector => detector%next       
          call move(detector_to_move, recv_lists(p), detector_list)
       end do
    end do

  end subroutine zoltan_detectors_finalise_migration

#endif
  
end module zoltan_detectors
