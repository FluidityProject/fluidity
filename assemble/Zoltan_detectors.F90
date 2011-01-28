#include "confdefs.h"

module zoltan_detectors

#ifdef HAVE_ZOLTAN

  use detector_data_types
  use detector_tools
  ! needed for getting detector%name
  use diagnostic_variables

  use fields

  implicit none

  public :: pack_detector, unpack_detector, update_detector
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
    detector%name=name_of_detector_in_read_order(detector%id_number)    
    
  end subroutine update_detector

#endif
  
end module zoltan_detectors
