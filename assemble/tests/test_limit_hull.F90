!! test matrix-free PETSc solve 
#include "fdebug.h"
subroutine test_limit_hull

  use unittest_tools
  use slope_limiters_dg
  use vector_tools
  use fldebug
  implicit none
  
  logical :: fail=.false., warn=.false.
  real, dimension(:,:), allocatable :: hull, vec
  real, dimension(:,:), pointer :: hull_p
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
  
  allocate(vec(2,30))
  vec(1,:) = (/0.55999695,  0.19157609,  0.67084366,  0.53064358,  0.42953486,&
       0.17604345,  0.08961044,  0.54511186,  0.19793959,  0.73359114,&
       0.86645983,  0.61999456,  0.94238344,  0.06399593,  0.71469049,&
       0.63856371,  0.09311535,  0.30548422,  0.01811958,  0.87038398,&
       0.64970014,  0.24037835,  0.20051147,  0.49689274,  0.7896298,&
       0.71427819,  0.09383192,  0.50716164,  0.50429313,  0.40448067/)
  vec(2,:) = (/0.10077083,  0.12588343,  0.24690842,  0.71961343,  0.81142182,&
       0.97042683,  0.75487892,  0.48052802,  0.7666026 ,  0.00432256,&
       0.10182016,  0.70007533,  0.64223023,  0.17429445,  0.85179439,&
       0.45244942,  0.76442205,  0.22065827,  0.92986088,  0.5579516,&
       0.45782496,  0.90416389,  0.45292503,  0.57262274,  0.28157226,&
       0.87835304,  0.04570724,  0.71090654,  0.9608029 ,  0.70510874/)
  call convex_hull(Vec,Hull_p)
  Ubar = Hull_p(:,1)
  dU = (/2.0,2.0/)
  alpha = limit_hull(Ubar,dU,hull_p)
  fail = alpha<0.0.or.alpha>=1.0
  call report_test('[RandomHullOutside]', fail, .false., &
       & "Point is outside hull.")
  dU = (/-2.0,-2.0/)
  alpha = limit_hull(Ubar,dU,hull_p)
  fail = alpha<0.0.or.alpha>=1.0
  call report_test('[RandomHullOutside2]', fail, .false., &
       & "Point is outside hull.")
  dU = (/-2.0,2.0/)
  alpha = limit_hull(Ubar,dU,hull_p)
  fail = alpha<0.0.or.alpha>=1.0
  call report_test('[RandomHullOutside3]', fail, .false., &
       & "Point is outside hull.")
  dU = (/0.55999695-Hull_p(1,1),0.10077083-Hull_p(2,1)/)
  alpha = limit_hull(Ubar,dU,hull_p)
  fail = abs(alpha-1.0)>1.0e-8
  call report_test('[RandomHullInside]', fail, .false., &
       & "Point is inside hull.")
  dU = (/0.5*(Hull_p(1,5)+Hull_p(1,6))-Hull_p(1,1),&
       0.5*(Hull_p(2,5)+Hull_p(2,6))-Hull_p(2,1)/)/0.9
  alpha = limit_hull(Ubar,dU,hull_p)
  fail = abs(alpha-0.9)>1.0e-8
  call report_test('[RandomHull_0.9]', fail, .false., &
       & "Point is outside hull, alpha = 0.9.")

  Ubar = sum(Hull_p,2)/size(Hull_p,2)
  dU = (0.5*(Hull_p(:,5)+Hull_p(:,6))-Ubar)/0.9
  alpha = limit_hull(Ubar,dU,hull_p)
  fail = abs(alpha-0.9)>1.0e-8
  call report_test('[RandomMiddleHull_0.9]', fail, .false., &
       & "Point is outside hull, alpha = 0.9.")
  deallocate(hull_p,vec)
end subroutine test_limit_hull
