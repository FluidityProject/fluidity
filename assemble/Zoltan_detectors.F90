#include "confdefs.h"
#include "fdebug.h"

module zoltan_detectors

#ifdef HAVE_ZOLTAN

  use zoltan, only: zoltan_int
  use zoltan_global_variables, only: zoltan_global_uen_to_new_local_numbering, zoltan_global_uen_to_old_local_numbering, zoltan_global_old_local_numbering_to_uen, zoltan_global_new_positions, zoltan_global_new_positions_mesh_nhalos
  use data_structures, only: has_key, fetch
  use halos_derivation, only: ele_owner

  use detector_data_types
  use detector_tools
  use diagnostic_variables, only: default_stat, remove_det_from_current_det_list

  use fields

  implicit none

  public :: pack_detector, unpack_detector, update_detector, prepare_detectors_for_packing
  private

  contains

  subroutine pack_detector(detector, buff, ndims, ndata_per_det)
    ! Packs the element, position, id_number and type of the
    ! detector into buff
    type(detector_type), pointer, intent(in) :: detector
    real, intent(out) :: buff(ndata_per_det)
    integer, intent(in) :: ndims, ndata_per_det

    buff(1) = detector%element
    buff(2:ndims+1) = detector%position
    buff(ndims+2) = detector%id_number
    buff(ndims+3) = detector%type
    
  end subroutine pack_detector

  subroutine unpack_detector(detector, buff, ndims, ndata_per_det)
    ! Unpacks the element, position, id_number and type of the
    ! detector from buff
    type(detector_type), pointer, intent(inout) :: detector
    real, intent(in) :: buff(ndata_per_det)
    integer, intent(in) :: ndims, ndata_per_det

    detector%element = buff(1)
    detector%position = reshape(buff(2:ndims+1),(/ndims/))
    detector%id_number = buff(ndims+2)
    detector%type = buff(ndims+3)
    
  end subroutine unpack_detector

  subroutine update_detector(detector, positions)
    ! Updates the local, inital_owner and local_coords for a
    ! detector based off other information from the detector
    ! and the vector field of positions
    type(detector_type), pointer, intent(inout) :: detector
    type(vector_field), intent(in) :: positions

    detector%local = .true.
    if (detector%type == STATIC_DETECTOR) then
       detector%initial_owner=-1
    else
       detector%initial_owner=getprocno()
    end if
    detector%local_coords=local_coords(positions,detector%element,detector%position)
    detector%name=default_stat%name_of_detector_in_read_order(detector%id_number)    
    
  end subroutine update_detector

  subroutine prepare_detectors_for_packing(ndets_in_ele, detector_list, to_pack_detectors_list, num_ids, global_ids)
    ! Plan is to first count how many detectors are in each element we're sending
    ! While doing this if we find a detector which will need to be sent, we remove
    ! it from detectors_list and add it to the end of the to_pack_detectors_list
    ! It's important it goes on the end of this list as we're relying on the 
    ! counting being done in the same order as we'll do the packing

    integer, intent(out), dimension(:) :: ndets_in_ele
    type(detector_linked_list), intent(inout), target :: detector_list
    type(detector_linked_list), intent(out), target :: to_pack_detectors_list
    integer(zoltan_int), intent(in) :: num_ids 
    integer(zoltan_int), intent(in), dimension(*) :: global_ids
    
    integer :: i, j
    integer :: old_universal_element_number, old_local_element_number, new_local_element_number
    integer :: new_ele_owner
    
    type(detector_type), pointer :: detector => null(), detector_to_move => null()
    
    integer :: original_detector_list_length
    
    ewrite(1,*) "In prepare_detectors_for_packing"
    
    ewrite(3,*) "Length of detector list BEFORE prepare_detectors_for_packing: ", detector_list%length
    
    assert(num_ids == size(ndets_in_ele))
    
    ! loop over all the elements we're interested in
    do i=1, num_ids
       
       ! work out some details about the current element
       old_universal_element_number = global_ids(i)
       old_local_element_number = fetch(zoltan_global_uen_to_old_local_numbering, old_universal_element_number)
       
       if (has_key(zoltan_global_uen_to_new_local_numbering, old_universal_element_number)) then
          new_local_element_number = fetch(zoltan_global_uen_to_new_local_numbering, old_universal_element_number)
          new_ele_owner = ele_owner(new_local_element_number, zoltan_global_new_positions%mesh, zoltan_global_new_positions%mesh%halos(zoltan_global_new_positions_mesh_nhalos))
       else
          new_ele_owner = -1
       end if
       
       if (new_ele_owner == getprocno()) then
          ndets_in_ele(i) = 0
       else
          ! start at the beginning of the detector list
          detector => detector_list%firstnode
          
          original_detector_list_length = detector_list%length
          
          ! search through all the detectors on this processor
          do j=1, original_detector_list_length
             ! check whether detector is in this element
             if (detector%element == old_local_element_number) then
                ! increment the number of detectors in that element
                ndets_in_ele(i) = ndets_in_ele(i) + 1
                
                detector_to_move => detector
                detector => detector%next

                ! Update detector%element to be universal element number
                ! so we can unpack to new element number
                detector_to_move%element = old_universal_element_number
               
                ! remove the detector from detector_list
                call remove_det_from_current_det_list(detector_list, detector_to_move)

                ! add detector to list of detectors we need to pack
                call insert(to_pack_detectors_list, detector_to_move)

                ewrite(3,*) "Detector ", detector_to_move%id_number, "(ID) removed from detector_list and added to to_pack_detectors_list."

                detector_to_move => null()

             end if
          end do
          assert(detector_list%length == (original_detector_list_length - ndets_in_ele(i)))
       end if
       
    end do

    assert(sum(ndets_in_ele(:)) == to_pack_detectors_list%length)
    
    ewrite(3,*) "Length of detector list AFTER prepare_detectors_for_packing: ", detector_list%length
    
    ewrite(1,*) "Exiting prepare_detectors_for_packing"
    
  end subroutine prepare_detectors_for_packing

#endif
  
end module zoltan_detectors
