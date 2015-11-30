subroutine test_string_hash_table

  use data_structures
  use unittest_tools
  use iso_c_binding, only : c_loc
  implicit none

  type(string_hash_table) :: shash
  integer :: len, i
  logical :: fail

  character(len=3) :: test_keys(3)=["One","Two","Six"]
  integer :: test_vals(3)=[10,20,60]

  integer, dimension(:), pointer :: ptr1, ptr2

  call allocate(shash)
  call insert(shash, "One", 10)
  call insert(shash, "Two", 20)
  call insert(shash, "Six", 60)

!  len = key_count(ihash)
!  fail = (len /= 3)
!  call report_test("[key_count]", fail, .false., "Should be 3")

  do i=1,3
    fail = (fetch(shash, test_keys(i)) /= test_vals(i))
    call report_test("[fetch]", fail, .false., "Should give i*10")
  end do 

  allocate(ptr1(4))
  ptr1=[(i,i=1,4)]

  call insert(shash, "Pointer1", c_loc(ptr1))
  ptr2=>fetch(shash,"Pointer1",[4])

  do i=1,4
     fail = (ptr2(i) /= i)
     call report_test("[fetch_pointer]", fail, .false., "Should give i")
  end do

  fail = has_key(shash, "Nine")
  call report_test("[string_hash_table_has_value]", fail, .false., "Should be .false.!")

  fail = .not. has_key(shash, "Two")
  call report_test("[string_hash_table_has_value]", fail, .false., "Should be .true.!")

  call remove(shash, "Six")
  fail = has_key(shash, "Six")
  call report_test("[string_hash_table_has_value]", fail, .false., "Should be .false.!")

  deallocate(ptr2)
  call deallocate(shash)

end subroutine test_string_hash_table
