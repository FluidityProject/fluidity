module element_set

  implicit none

  interface

     subroutine ele_get_ele(i,ele)
       integer, intent(in)  :: i
       integer, intent(out) :: ele
     end subroutine ele_get_ele

     subroutine ele_fetch_list(arr)
       use iso_c_binding
       integer(c_int), dimension(*), intent(out) :: arr
      end subroutine ele_fetch_list

      subroutine ele_get_size(size)
        integer, intent(out) :: size
      end subroutine ele_get_size

      subroutine ele_add_to_set(i)
        integer, intent(in) :: i
      end subroutine ele_add_to_set

   end interface

  contains

  subroutine eleset_add(i)
    !!< Add i to the set of elements to be considered.
    integer, intent(in) :: i

    call ele_add_to_set(i)
  end subroutine

  subroutine eleset_get_size(size)
    integer, intent(out) :: size

    call ele_get_size(size)
  end subroutine 

  subroutine eleset_fetch_list(arr)
    !!< Fetch the list and clear it.
    integer, dimension(:), intent(out) :: arr
    call ele_fetch_list(arr)
  end subroutine

  subroutine eleset_get_ele(i, ele)
    !!< Get element i in the set.
    integer, intent(in) :: i
    integer, intent(out) :: ele
    call ele_get_ele(i, ele)
  end subroutine eleset_get_ele


end module element_set
