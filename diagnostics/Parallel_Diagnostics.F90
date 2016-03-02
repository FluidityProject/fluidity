#include "fdebug.h"

module parallel_diagnostics

  use fldebug
  use fields
  use halos

  implicit none
  
  private
  
  public :: calculate_node_halo, calculate_universal_numbering, &
    & calculate_element_halo, calculate_element_ownership, &
    & calculate_element_universal_numbering

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
  
  subroutine calculate_universal_numbering(s_field)
    type(scalar_field), intent(inout) :: s_field
    
    integer :: i, nhalos
    type(halo_type), pointer :: halo
    
    nhalos = halo_count(s_field)
    if(nhalos > 0) then
      halo => s_field%mesh%halos(nhalos)
      do i = 1, node_count(s_field)
        call set(s_field, i, float(halo_universal_number(halo, i)))
      end do
    else
      do i = 1, node_count(s_field)
        call set(s_field, i, float(i))
      end do
    end if
  
  end subroutine calculate_universal_numbering

  subroutine calculate_element_halo(s_field)
    type(scalar_field), intent(inout) :: s_field
    
    integer :: i, j
    type(element_type), pointer :: shape

    assert(ele_count(s_field) > 0)
    shape => ele_shape(s_field, 1)
    if(shape%degree /= 0) then
      FLExit("element_halo diagnostic requires a degree 0 mesh")
    end if
    
    call zero(s_field)
    do i = element_halo_count(s_field), 1, -1
      do j = 1, halo_proc_count(s_field%mesh%element_halos(i))
        call set(s_field, halo_sends(s_field%mesh%element_halos(i), j), spread(float(i), 1, halo_send_count(s_field%mesh%element_halos(i), j)))
        call set(s_field, halo_receives(s_field%mesh%element_halos(i), j), spread(-float(i), 1, halo_receive_count(s_field%mesh%element_halos(i), j)))
      end do
    end do
    
  end subroutine calculate_element_halo

  subroutine calculate_element_ownership(s_field)
    type(scalar_field), intent(inout) :: s_field

    integer :: i, nhalos
    type(element_type), pointer :: shape
    type(halo_type), pointer :: ele_halo

    assert(ele_count(s_field) > 0)
    shape => ele_shape(s_field, 1)
    if(shape%degree /= 0) then
      FLExit("element_halo_ownership diagnostic requires a degree 0 mesh")
    end if
    assert(node_count(s_field) == ele_count(s_field))

    nhalos = element_halo_count(s_field)
    if(nhalos > 0) then
      ele_halo => s_field%mesh%element_halos(nhalos)
      do i = 1, node_count(s_field)
        call set(s_field, i, float(halo_node_owner(ele_halo, i)))
      end do
    else
      call set(s_field, 1.0)
    end if
    
  end subroutine calculate_element_ownership

  subroutine calculate_element_universal_numbering(s_field)
    type(scalar_field), intent(inout) :: s_field
    
    integer :: i, nhalos
    type(element_type), pointer :: shape
    type(halo_type), pointer :: ele_halo

    assert(ele_count(s_field) > 0)
    shape => ele_shape(s_field, 1)
    if(shape%degree /= 0) then
      FLExit("element_universal_numbering diagnostic requires a degree 0 mesh")
    end if
    assert(node_count(s_field) == ele_count(s_field))

    nhalos = element_halo_count(s_field)
    if(nhalos > 0) then
      ele_halo => s_field%mesh%element_halos(nhalos)
      do i = 1, node_count(s_field)
        call set(s_field, i, float(halo_universal_number(ele_halo, i)))
      end do
    else
      do i = 1, node_count(s_field)
        call set(s_field, i, float(i))
      end do
    end if
    
  end subroutine calculate_element_universal_numbering

end module parallel_diagnostics
