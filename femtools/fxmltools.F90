module fxmltools

  use iso_c_binding

  implicit none

  interface c_wrap
     module procedure c_wrap_char_array, c_wrap_string
  end interface c_wrap

  interface asbytes
     module procedure ints_asbytes, longs_asbytes, &
          floats_asbytes, vecfloats_asbytes, &
          doubles_asbytes, vecdoubles_asbytes, &
          c_ptr_asbytes
  end interface asbytes

  contains

    function c_wrap_char_array(string) result(c_wrap)
      character, intent(in) :: string(:)
      character :: c_wrap(size(string)+1) 

      c_wrap(1:size(string)) = string
      c_wrap(size(string)+1:size(string)+1) = C_NULL_CHAR

    end function c_wrap_char_array

    function c_wrap_string(string) result(c_wrap)
      character(len=*), intent(in) :: string
      character :: c_wrap(len(string)+1)
      
      integer :: i

      do i = 1, len(string)
         c_wrap(i) = string(i:i)
      end do
      c_wrap(len(string)+1) = C_NULL_CHAR

    end function c_wrap_string

    function ints_asbytes(i) result(bytes)
      integer(c_int), intent(in), dimension(:) :: i
      character(kind=c_char, len=4*size(i)) :: bytes
      bytes = transfer(i, bytes)
    end function ints_asbytes

    function longs_asbytes(i) result(bytes)
      integer(c_long), intent(in), dimension(:) :: i
      character(kind=c_char, len=8*size(i)) :: bytes
      bytes = transfer(i, bytes)
    end function longs_asbytes

    function floats_asbytes(a) result(bytes)
      real(c_float), intent(in), dimension(:) :: a
      character(kind=c_char, len=4*size(a)) :: bytes
      bytes = transfer(a, bytes)
    end function floats_asbytes

    function doubles_asbytes(a) result(bytes)
      real(c_double), intent(in), dimension(:) :: a
      character(kind=c_char, len=8*size(a)) :: bytes
      bytes = transfer(a, bytes)
    end function doubles_asbytes

    function vecfloats_asbytes(a) result(bytes)
      real(c_float), intent(in), dimension(:,:) :: a
      character(kind=c_char, len=4*size(a)) :: bytes
      bytes = asbytes(reshape(a,[size(a)]))
    end function vecfloats_asbytes

    function vecdoubles_asbytes(a) result(bytes)
      real(c_double), intent(in), dimension(:,:) :: a
      character(kind=c_char, len=8*size(a)) :: bytes
      bytes = asbytes(reshape(a,[size(a)]))
    end function vecdoubles_asbytes

    function c_ptr_asbytes(a, bytesize) result(bytes)
      type(c_ptr), intent(in):: a
      integer(c_int), intent(in) :: bytesize
      integer(c_int) :: tmp
      character(kind=c_char, len=bytesize), target :: bytes

      interface
         subroutine memcpy(dest, src, n) bind(c)
           use iso_c_binding
           type(c_ptr), value :: dest, src
           integer(c_int), value :: n
         end subroutine  memcpy
      end interface
      call memcpy(c_loc(bytes), a, bytesize)
    end function c_ptr_asbytes
    
    function reverse(a)
      character(len=*) :: a
      character(len=len(a)) :: reverse

      integer :: i, j, n
      n = len(a)
      do i=1, n
         j=n+1-i
         reverse(j:j) = a(i:i)
      end do
    end function reverse

    function str(a)
      character, dimension(:) :: a
      character(len=size(a)) :: str

      integer :: i
      do i=1, size(a)
         str(i:i) = a(i)
      end do
    end function str

    subroutine copy2chars(s, a)
      character (len=*), intent(in) :: s
      character, allocatable, intent(out) :: a(:)

      integer :: i

      allocate(a(len(s)))

      do i=1, len(s)
         a(i) = s(i:i)
      end do
    end subroutine copy2chars

    pure function int2str_len(i)

      !!< Count number of digits in i.

      integer, intent(in) :: i
      integer :: int2str_len 
      
      int2str_len = 1
      if (i/=0) int2str_len = int2str_len + floor(log10(abs(real(i))))
      if (i<0) int2str_len = int2str_len + 1

    end function int2str_len
    
    function int2str (i)

      !!< Convert integer i into a fortran string.

      integer, intent(in) :: i
      character(len=int2str_len(i)) :: int2str

      write(int2str,"(i0)") i

    end function int2str

    pure function intarray2str_len(i)

      integer, intent(in), dimension(:) :: i
      integer :: intarray2str_len, j 

      intarray2str_len = 0

      do j=1,size(i)
         intarray2str_len = intarray2str_len + int2str_len(i(j))
      end do
      if (size(i)>0) intarray2str_len = intarray2str_len+size(i)-1

    end function intarray2str_len

    function intarray2str(i)

      !!< Convert integer array i into a c string string.

      integer, intent(in), dimension(:) :: i
      character(len=intarray2str_len(i)) :: intarray2str
      integer :: j, k, n, lens(size(i))

      n = size(i)

      do j=1, n
         lens(j) = int2str_len(i(j))
      end do

      k=0
      do j=1, n-1
         write(intarray2str(k+1:k+lens(j)), "(i0)") i(j)
         k=k+lens(j)+1
         intarray2str(k:k) = ' '
      end do
      if (n>0) write(intarray2str(k+1:k+lens(n)), "(i0)") i(n)

    end function intarray2str

end module fxmltools
