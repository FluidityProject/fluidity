!! test matrix-free PETSc solve 
#include "fdebug.h"
subroutine test_limit_hull

  use unittest_tools
  use slope_limiters_dg
  use vector_tools
  use fldebug
  implicit none
  
  logical :: fail=.false., warn=.false.
  real, dimension(:,:), allocatable :: hull
  real, dimension(2) :: Ubar, dU
  real :: alpha

  !Set up a square hull
  allocate(hull(2,4))
  hull(1,:) = (/0.0,0.0,1.0,1.0/)
  hull(2,:) = (/0.0,1.0,1.0,0.0/)
  Ubar = (/0.0,0.0/)
  dU = (/1.0,1.0/)
  alpha = limit_hull(Ubar,dU,hull)
  fail = abs(alpha-1.0)>1.0e-8
  call report_test("[HullCorner]", fail, .false., &
    "Point is inside hull. No limiting needed.")

  Ubar = (/0.0,0.0/)
  dU = (/1.0,1.0/)
  dU = dU/0.9
  alpha = limit_hull(Ubar,dU,hull)
  print *, alpha
  print *, 'THIS IS BROKEN!'
  fail = abs(alpha-0.9)>1.0e-8  
  call report_test("[HullOutsideCorner]", fail, .false., &
    "Limiter value should be 0.9")

  Ubar = (/0.0,0.0/)
  dU = (/1.0,1.0/)
  dU = dU/1.1
  assert(.not.outside_hull(Ubar+dU,hull))
  alpha = limit_hull(Ubar,dU,hull)
  fail = abs(alpha-1.0)>1.0e-8  
  call report_test("[HullCompletelyCorner]", fail, .false., &
    "Point is well inside hull.")

  deallocate(hull)
  
end subroutine test_limit_hull
