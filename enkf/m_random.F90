#include "fdebug.h"

module m_random

implicit none

contains

  subroutine random(work1,n)
    !  Returns a vector of random values N(variance=1,mean=0)
    implicit none
    integer, intent(in) :: n
    real,   intent(out) :: work1(n)
    real,   allocatable :: work2(:)
    real, parameter   ::  pi=3.141592653589

    allocate (work2(n))

    call random_number(work1)
    call random_number(work2)
    work1= sqrt(-2.0*log(work1))*cos(2.0*pi*work2)

    deallocate(work2)
  end subroutine random


  subroutine randn(n, vect)
    implicit none
    integer, intent(in) :: n
    real, intent(out) :: vect(n)

    integer :: i
    real :: a(2), r

    i = 0
    do while (i < n)
       call random_number(a)
       a = 2.0 * a - 1.0
       r = a(1) * a(1) + a(2) * a(2)
       if (r > 1.0) then
          cycle
       end if
       i = i + 1
       ! assume that r is never equal to 0 - PS
       r = sqrt(-2.0 * log(r) / r);
       vect(i) = r * a(1);
       if (i == n) then
          exit
       end if
       i = i + 1
       vect(i) = r * a(2);
    end do
  end subroutine randn

end module m_random
