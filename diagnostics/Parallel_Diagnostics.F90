#include "fdebug.h"

module parallel_diagnostics

  use fields
  use fldebug
  use halos

  implicit none
  
  private
  
  public :: calculate_node_halo, calculate_element_halo

contains

  subroutine calculate_node_halo(s_field)
    type(scalar_field), intent(inout) :: s_field
    
    integer :: i, j
    
    call zero(s_field)
    do i = halo_count(s_field), 1, -1
      do j = 1, halo_proc_count(s_field%mesh%halos(i))
        call set(s_field, halo_sends(s_field%mesh%halos(i), j), spread(float(i), 1, halo_send_count(s_field%mesh%halos(i), j)))
        call set(s_field, halo_receives(s_field%mesh%halos(i), j), spread(-float(i), 1, halo_receive_count(s_field%mesh%halos(i), j)))
      end do
    end do
    
  end subroutine calculate_node_halo

  subroutine calculate_element_halo(s_field)
    type(scalar_field), intent(inout) :: s_field
    
    integer :: i, j
    type(element_type), pointer :: shape

    assert(ele_count(s_field) > 0)
    shape => ele_shape(s_field, 1)
    if(shape%degree /= 0) then
      FLAbort("element_halo diagnostic requires a degree 0 mesh")
    end if
    
    call zero(s_field)
    do i = element_halo_count(s_field), 1, -1
      do j = 1, halo_proc_count(s_field%mesh%element_halos(i))
        call set(s_field, halo_sends(s_field%mesh%element_halos(i), j), spread(float(i), 1, halo_send_count(s_field%mesh%element_halos(i), j)))
        call set(s_field, halo_receives(s_field%mesh%element_halos(i), j), spread(-float(i), 1, halo_receive_count(s_field%mesh%element_halos(i), j)))
      end do
    end do
    
  end subroutine calculate_element_halo

end module parallel_diagnostics
