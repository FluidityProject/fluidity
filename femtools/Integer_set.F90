module integer_set_module
  ! Don't use this directly, use data_structures
  use iso_c_binding, only: c_ptr
  type integer_set
    type(c_ptr) :: address
  end type integer_set

  interface
    subroutine integer_set_create_c(i)
      use iso_c_binding, only: c_ptr
      type(c_ptr), intent(out) :: i
    end subroutine integer_set_create_c

    subroutine integer_set_delete_c(i)
      use iso_c_binding, only: c_ptr
      type(c_ptr), intent(inout) :: i
    end subroutine integer_set_delete_c

    subroutine integer_set_insert_c(i, v, c)
      use iso_c_binding, only: c_ptr
      type(c_ptr), intent(inout) :: i
      integer, intent(in) :: v
      integer, intent(out) :: c
    end subroutine integer_set_insert_c

    pure subroutine integer_set_length_c(i, l)
      use iso_c_binding, only: c_ptr
      type(c_ptr), intent(in) :: i
      integer, intent(out) :: l
    end subroutine integer_set_length_c

    subroutine integer_set_fetch_c(i, idx, val)
      use iso_c_binding, only: c_ptr
      type(c_ptr), intent(in) :: i
      integer, intent(in) :: idx
      integer, intent(out) :: val
    end subroutine integer_set_fetch_c

    subroutine integer_set_has_value_c(i, val, bool)
      use iso_c_binding, only: c_ptr
      type(c_ptr), intent(in) :: i
      integer, intent(in) :: val
      integer, intent(out) :: bool
    end subroutine integer_set_has_value_c
  end interface

  interface allocate
    module procedure integer_set_allocate_single, integer_set_allocate_vector
  end interface

  interface insert
    module procedure integer_set_insert
  end interface

  interface deallocate
    module procedure integer_set_delete_single, integer_set_delete_vector
  end interface

  interface has_value
    module procedure integer_set_has_value
  end interface

  interface key_count
    module procedure integer_set_length_single, integer_set_length_vector
  end interface

  interface fetch
    module procedure integer_set_fetch
  end interface

  private
  public :: integer_set, allocate, deallocate, has_value, key_count, fetch, insert, &
          & set_complement, set2vector, set_intersection

  contains 

  subroutine integer_set_allocate_single(iset)
    type(integer_set), intent(out) :: iset
    iset = integer_set_create()
  end subroutine integer_set_allocate_single
  
  subroutine integer_set_allocate_vector(iset)
    type(integer_set), dimension(:), intent(out) :: iset
    
    integer :: i
    
    do i = 1, size(iset)
      call allocate(iset(i))
    end do
  
  end subroutine integer_set_allocate_vector

  function integer_set_create() result(iset)
    type(integer_set) :: iset
    call integer_set_create_c(iset%address)
  end function integer_set_create

  subroutine integer_set_delete_single(iset)
    type(integer_set), intent(inout) :: iset
    call integer_set_delete_c(iset%address)
  end subroutine integer_set_delete_single
  
  subroutine integer_set_delete_vector(iset)
    type(integer_set), dimension(:), intent(inout) :: iset
    
    integer :: i
    
    do i = 1, size(iset)
      call deallocate(iset(i))
    end do
    
  end subroutine integer_set_delete_vector

  subroutine integer_set_insert(iset, val, changed)
    type(integer_set), intent(inout) :: iset
    integer, intent(in) :: val
    logical, intent(out), optional :: changed
    integer :: lchanged

    call integer_set_insert_c(iset%address, val, lchanged)

    if (present(changed)) then
      changed = (lchanged == 1)
    end if
  end subroutine integer_set_insert

  pure function integer_set_length_single(iset) result(len)
    type(integer_set), intent(in) :: iset
    integer :: len

    call integer_set_length_c(iset%address, len)
  end function integer_set_length_single
  
  pure function integer_set_length_vector(iset) result(len)
    type(integer_set), dimension(:), intent(in) :: iset

    integer, dimension(size(iset)) :: len
    
    integer :: i
    
    do i = 1, size(iset)
      len(i) = key_count(iset(i))
    end do
  
  end function integer_set_length_vector

  function integer_set_fetch(iset, idx) result(val)
    type(integer_set), intent(in) :: iset
    integer, intent(in) :: idx
    integer :: val

    call integer_set_fetch_c(iset%address, idx, val)
  end function integer_set_fetch

  function integer_set_has_value(iset, val) result(bool)
    type(integer_set), intent(in) :: iset
    integer, intent(in) :: val
    logical :: bool

    integer :: lbool
    call integer_set_has_value_c(iset%address, val, lbool)
    bool = (lbool == 1)
  end function integer_set_has_value

  subroutine set_complement(complement, universe, current)
    ! complement = universe \ current
    type(integer_set), intent(out) :: complement
    type(integer_set), intent(in) :: universe, current
    integer :: i, val

    call allocate(complement)
    do i=1,key_count(universe)
      val = fetch(universe, i)
      if (.not. has_value(current, val)) then
        call insert(complement, val)
      end if
    end do
  end subroutine set_complement

  subroutine set_intersection(intersection, A, B)
    ! intersection = A n B
    type(integer_set), intent(out) :: intersection
    type(integer_set), intent(in) :: A, B
    integer :: i, val

    call allocate(intersection)
    do i=1,key_count(A)
      val = fetch(A, i)
      if (has_value(B, val)) then
        call insert(intersection, val)
      end if
    end do
  end subroutine set_intersection

  function set2vector(iset) result(vec)
    type(integer_set), intent(in) :: iset
    integer, dimension(key_count(iset)) :: vec
    integer :: i

    do i=1,key_count(iset)
      vec(i) = fetch(iset, i)
    end do
  end function set2vector

end module integer_set_module
