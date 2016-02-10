#include "confdefs.h"
#include "fdebug.h"

module zoltan_detectors

#ifdef HAVE_ZOLTAN

  use zoltan, only: zoltan_int
  use data_structures, only: has_key, fetch
  use detector_data_types
  use fields
  use zoltan_global_variables, only: zoltan_global_uen_to_new_local_numbering, zoltan_global_old_local_numbering_to_uen, zoltan_global_new_positions
  use detector_tools
  use detector_parallel


  implicit none

  public :: prepare_detectors_for_packing
  private

  contains

  subroutine prepare_detectors_for_packing(ndets_in_ele, to_pack_detector_lists, num_ids, global_ids)
    ! Goes through all local detectors and moves the ones we want to send to the according
    ! to_pack_detector_list for the element we found the detector in
    integer, intent(out), dimension(:) :: ndets_in_ele
    type(detector_linked_list), dimension(:), intent(inout), target :: to_pack_detector_lists
    integer(zoltan_int), intent(in) :: num_ids 
    integer(zoltan_int), intent(in), dimension(*) :: global_ids
    
    integer :: i, j, det_list, new_ele_owner, total_det_to_pack, detector_uen
    integer :: new_local_element_number
    
    type(detector_list_ptr), dimension(:), pointer :: detector_list_array => null()
    type(detector_type), pointer :: detector => null(), detector_to_move => null()
    logical :: found_det_element
    
    ewrite(1,*) "In prepare_detectors_for_packing"    
    
    assert(num_ids == size(ndets_in_ele))
    assert(num_ids == size(to_pack_detector_lists))

    ! loop through all registered detector lists
    call get_registered_detector_lists(detector_list_array)
    do det_list = 1, size(detector_list_array)

       ! search through all the local detectors in this list
       detector => detector_list_array(det_list)%ptr%first
       detector_loop: do while (associated(detector))
          ! store the list ID with the detector, so we can map the detector back when receiving it
          detector%list_id=det_list

          ! translate detector element to uen
          if (.not. has_key(zoltan_global_old_local_numbering_to_uen, detector%element)) then
             ewrite(-1,*) "No uen found in Zoltan for detector ", detector%id_number, " with local element number ", detector%element
             FLAbort("No universal element number for detector found in Zoltan")
          end if
          detector_uen = fetch(zoltan_global_old_local_numbering_to_uen, detector%element)

          ! loop over all the elements we're interested in
          found_det_element=.false.
          element_loop: do j=1, num_ids

             ! check whether detector is in this element
             if (detector_uen == global_ids(j)) then
                found_det_element=.true.

                ! work out new owner
                if (has_key(zoltan_global_uen_to_new_local_numbering, detector_uen)) then
                   new_local_element_number = fetch(zoltan_global_uen_to_new_local_numbering, detector_uen)
                   new_ele_owner = element_owner(zoltan_global_new_positions%mesh, new_local_element_number)
                else
                   new_ele_owner = -1
                end if

                ! check whether old owner is new owner
                if (new_ele_owner == getprocno()) then
                   ndets_in_ele(j) = 0
                   detector => detector%next
                else
                   ! If not, move detector to the pack_list for this element and
                   ! increment the number of detectors in that element
                   ndets_in_ele(j) = ndets_in_ele(j) + 1
                
                   detector_to_move => detector
                   detector => detector%next

                   ! Update detector%element to be universal element number
                   ! so we can unpack to new element number
                   detector_to_move%element = detector_uen
               
                   ! Move detector to list of detectors we need to pack
                   call move(detector_to_move, detector_list_array(det_list)%ptr, to_pack_detector_lists(j))
                   detector_to_move => null()
                end if

                ! We found the right element, so we can skip the others
                exit element_loop
             end if          
          end do element_loop

          ! If we didn't find an element for the detector, we have to advance it
          if (.not.found_det_element) detector => detector%next
       end do detector_loop
    end do

    ! Sanity checks and logging
    total_det_to_pack=0
    do i=1, num_ids
       assert(ndets_in_ele(i) == to_pack_detector_lists(i)%length)
       total_det_to_pack=total_det_to_pack+to_pack_detector_lists(i)%length
    end do
    ewrite(2,*) "Moved", total_det_to_pack, "detectors to to_pack_detector_lists"
    ewrite(1,*) "Exiting prepare_detectors_for_packing"
    
  end subroutine prepare_detectors_for_packing

#endif
  
end module zoltan_detectors
