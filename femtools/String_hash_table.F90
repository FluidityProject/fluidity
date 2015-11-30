#include "fdebug.h"

module string_hash_table_module
  ! Don't use this directly, use data_structures
  use iso_c_binding, only: c_ptr
  use fldebug

  implicit none

  type string_hash_table
    type(c_ptr) :: address
  end type string_hash_table

  integer :: MAX_LEN=254

  interface
    subroutine string_hash_table_create_c(i) bind(c)
      use iso_c_binding, only: c_ptr
      type(c_ptr), intent(out) :: i
    end subroutine string_hash_table_create_c

    subroutine string_hash_table_delete_c(i) bind(c)
      use iso_c_binding, only: c_ptr
      type(c_ptr), intent(in) :: i
    end subroutine string_hash_table_delete_c

    subroutine string_hash_table_insert_c(i, k, v) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int
      type(c_ptr), intent(in) :: i
      character(kind=c_char), intent(in) :: k(*)
      integer(c_int), intent(in) :: v
    end subroutine string_hash_table_insert_c

    subroutine string_hash_table_insert_pointer_c(i, k, p) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int
      type(c_ptr), intent(in) :: i
      character(kind=c_char), intent(in) :: k(*)
      type(c_ptr), intent(in) :: p
    end subroutine string_hash_table_insert_pointer_c

    subroutine string_hash_table_fetch_c(i, key, val) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int
      type(c_ptr), intent(in) :: i
      character(kind=c_char), intent(in) :: key(*)
      integer(c_int), intent(out) :: val
    end subroutine string_hash_table_fetch_c

    subroutine string_hash_table_fetch_pointer_c(i, key, ptr) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int
      type(c_ptr), intent(in) :: i
      character(kind=c_char), intent(in) :: key(*)
      type(c_ptr), intent(out) :: ptr
    end subroutine string_hash_table_fetch_pointer_c

    subroutine string_hash_table_remove_c(i, key, stat) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int
      type(c_ptr), intent(in) :: i
      character(kind=c_char), intent(in) :: key(*)
      integer(c_int), intent(out) :: stat
    end subroutine string_hash_table_remove_c

    subroutine string_hash_table_has_key_c(i, key, bool) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int
      type(c_ptr), intent(in) :: i
      character(kind=c_char), intent(in) :: key(*)
      integer(c_int), intent(out) :: bool
    end subroutine string_hash_table_has_key_c

    subroutine string_hash_table_get_first_c(i, key, key_len, val,status) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int
      type(c_ptr), intent(in) :: i
      character(kind=c_char), intent(in) :: key(*)
      integer(c_int), intent(out) :: key_len, status, val
    end subroutine string_hash_table_get_first_c

    subroutine string_hash_table_get_next_c(i, key, key_len, val,status) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int
      type(c_ptr), intent(in) :: i
      character(kind=c_char), intent(in) :: key(*)
      integer(c_int), intent(out) :: key_len, status, val
    end subroutine string_hash_table_get_next_c

    subroutine string_hash_table_get_first_pointer_c(i, key, key_len, ptr, status) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int
      type(c_ptr), intent(in) :: i
      character(kind=c_char), intent(in) :: key(*)
      integer(c_int), intent(out) :: key_len, status
      type(c_ptr), intent(out) :: ptr
    end subroutine string_hash_table_get_first_pointer_c

    subroutine string_hash_table_get_next_pointer_c(i, key, key_len, ptr, status) bind(c)
      use iso_c_binding, only: c_ptr, c_char, c_int
      type(c_ptr), intent(in) :: i
      character(kind=c_char), intent(in) :: key(*)
      integer(c_int), intent(out) :: key_len, status
      type(c_ptr), intent(out) :: ptr
    end subroutine string_hash_table_get_next_pointer_c

  end interface

  interface allocate
    module procedure string_hash_table_allocate
  end interface

  interface insert
    module procedure string_hash_table_insert, string_hash_table_insert_pointer
  end interface

  interface remove
    module procedure string_hash_table_remove
  end interface

  interface deallocate
    module procedure string_hash_table_delete
  end interface

  interface has_key
    module procedure string_hash_table_has_key
  end interface

  interface fetch
    module procedure string_hash_table_fetch, string_hash_table_fetch_v, string_hash_table_fetch_integer_pointer
  end interface

  interface print
    module procedure print_hash_table
  end interface

  interface copy
    module procedure string_hash_table_copy
  end interface

  interface get_first
    module procedure string_hash_table_get_first, string_hash_table_get_first_pointer
  end interface

  interface get_next
    module procedure string_hash_table_get_next, string_hash_table_get_next_pointer
  end interface

  private
  public :: string_hash_table, allocate, deallocate, has_key, fetch, insert, &
       print, remove, copy, get_first, get_next

  contains 

    function c_wrap(str)
      use iso_c_binding, only: C_CHAR, C_NULL_CHAR
      character(len=*) :: str
      character(kind=c_char, len=MAX_LEN+1) :: c_wrap

      assert(len(str)<max_len)      
      c_wrap=str//C_NULL_CHAR
    end function c_wrap

  subroutine string_hash_table_copy(shash_copy, shash)
    type(string_hash_table), intent(out) :: shash_copy
    type(string_hash_table), intent(inout) :: shash
    
    integer :: key_val,stat
    character(len=MAX_LEN+1) :: key

    !!! note key is full length here to allow for null termination

    call allocate(shash_copy)

    stat=1
    call get_first(shash,key,key_val, stat)
    do while (stat/=0)
      call insert(shash_copy, key, key_val)
      call get_next(shash,key,key_val,stat)
    end do

  end subroutine string_hash_table_copy

  subroutine string_hash_table_allocate(shash)
    type(string_hash_table), intent(out) :: shash
    shash = string_hash_table_create()
  end subroutine string_hash_table_allocate

  function string_hash_table_create() result(shash)
    type(string_hash_table) :: shash
    call string_hash_table_create_c(shash%address)
  end function string_hash_table_create

  subroutine string_hash_table_delete(shash)
    type(string_hash_table), intent(in) :: shash
    call string_hash_table_delete_c(shash%address)
  end subroutine string_hash_table_delete

  subroutine string_hash_table_insert(shash, key, val)
    use iso_c_binding
    type(string_hash_table), intent(in) :: shash
    character(len=*) :: key
    integer, intent(in) :: val

    call string_hash_table_insert_c(shash%address, c_wrap(key), val)
  end subroutine string_hash_table_insert

  subroutine string_hash_table_insert_pointer(shash, key, ptr)
    use iso_c_binding
    type(string_hash_table), intent(in) :: shash
    character(len=*) :: key
    type(c_ptr), intent(in) :: ptr

    call string_hash_table_insert_pointer_c(shash%address, c_wrap(key), ptr)
  end subroutine string_hash_table_insert_pointer

  function string_hash_table_fetch(shash, key) result(val)
    type(string_hash_table), intent(in) :: shash
    character(len=*), intent(in) :: key
    integer :: val

    call string_hash_table_fetch_c(shash%address, c_wrap(key), val)
  end function string_hash_table_fetch

  subroutine string_hash_table_remove(shash, key)
    type(string_hash_table), intent(in) :: shash
    character(len=*), intent(in) :: key
    integer :: stat

    call string_hash_table_remove_c(shash%address, c_wrap(key), stat)
    assert(stat == 1)
  end subroutine string_hash_table_remove

  function string_hash_table_fetch_v(shash, keys) result(vals)
    type(string_hash_table), intent(inout) :: shash
    character(len=*), intent(in), dimension(:) :: keys
    integer, dimension(size(keys)) :: vals
    integer :: i

    do i=1,size(keys)
      call string_hash_table_fetch_c(shash%address, c_wrap(keys(i)), vals(i))
    end do
  end function string_hash_table_fetch_v

  function string_hash_table_fetch_integer_pointer(shash, key, shape) result(ptr)
    use iso_c_binding
    type(string_hash_table), intent(in) :: shash
    character(len=*), intent(in) :: key
    integer, dimension(:), intent(in) :: shape
    integer, dimension(:), pointer :: ptr

    type(c_ptr) :: lptr

    call string_hash_table_fetch_pointer_c(shash%address, c_wrap(key), lptr)
    call c_f_pointer(lptr,ptr,shape)

  end function string_hash_table_fetch_integer_pointer

  function string_hash_table_has_key(shash, key) result(bool)
    type(string_hash_table), intent(in) :: shash
    character(len=*), intent(in) :: key
    logical :: bool

    integer :: lbool

    call string_hash_table_has_key_c(shash%address, c_wrap(key), lbool)
    bool = (lbool == 1)
  end function string_hash_table_has_key

  subroutine string_hash_table_get_first(shash, key, val, stat)
    type(string_hash_table), intent(inout) :: shash
    character(len=MAX_LEN), intent(out) :: key
    integer, intent(out) :: val
    integer :: stat

    character(len=MAX_LEN+1) :: lkey
    integer :: key_len

    key=''
    call string_hash_table_get_first_c(shash%address, lkey, key_len, val, stat)
    key=lkey(:key_len)

  end subroutine string_hash_table_get_first

  subroutine string_hash_table_get_first_pointer(shash, key, ptr, stat)
    type(string_hash_table), intent(inout) :: shash
    character(len=MAX_LEN), intent(out) :: key
    type(c_ptr), intent(out) :: ptr
    integer :: stat

    character(len=MAX_LEN+1) :: lkey
    integer :: key_len

    key=''
    call string_hash_table_get_first_pointer_c(shash%address, lkey, key_len, ptr, stat)
    key=lkey(:key_len)

  end subroutine string_hash_table_get_first_pointer

  subroutine string_hash_table_get_next(shash, key, val, stat)
    type(string_hash_table), intent(in) :: shash
    character(len=MAX_LEN), intent(out) :: key
    integer, intent(out) :: val
    integer :: stat

    character(len=MAX_LEN+1) :: lkey
    integer :: key_len

    key=''
    call string_hash_table_get_next_c(shash%address, lkey, key_len, val, stat)
    key=lkey(:key_len)

  end subroutine string_hash_table_get_next

  subroutine string_hash_table_get_next_pointer(shash, key, ptr, stat)
    type(string_hash_table), intent(in) :: shash
    character(len=MAX_LEN), intent(out) :: key
    type(c_ptr), intent(out) :: ptr
    integer :: stat

    character(len=MAX_LEN+1) :: lkey
    integer :: key_len

    key=''
    call string_hash_table_get_next_pointer_c(shash%address, lkey, key_len,&
         ptr, stat)
    key=lkey(:key_len)

  end subroutine string_hash_table_get_next_pointer

  subroutine print_hash_table(shash, priority)
    type(string_hash_table), intent(inout) :: shash
    integer, intent(in) :: priority

    integer :: val, stat
    character(len=MAX_LEN) :: key

    ewrite(priority,*) "Writing hash table: "
    call get_first(shash,key,val,stat)
    do while (stat/=0)
      ewrite(priority,*) trim(key), " --> ", val
      call get_next(shash,key,val,stat)
    end do
  end subroutine print_hash_table

end module string_hash_table_module
