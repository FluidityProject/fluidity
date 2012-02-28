subroutine test_convex_hull2d
  use vector_tools
  use unittest_tools
  implicit none  
  real, dimension(:,:), allocatable :: Vec
  real, dimension(:,:), pointer :: Hull
  logical :: fail
  real, dimension(:), allocatable :: s

  allocate(Vec(2,5))
  Vec(1,:) = (/0.,1.,0.,1.,0.5/)
  Vec(2,:) = (/0.,0.,1.,1.,0.0/)
  call convex_hull(Vec,Hull)
  fail = size(Hull,2).ne.4
  call report_test("[SquareHull]", fail, .false., "Convex hull should have 4&
       & points")  
  deallocate(Hull)

  Vec(1,:) = (/0.,1.,0.,1.,0.5/)
  Vec(2,:) = (/0.,0.,1.,1.,-1.0e-7/)
  call convex_hull(Vec,Hull)
  fail = size(Hull,2).ne.5
  call report_test("[HouseHull]", fail, .false., "Convex hull should have 5&
       & points")  
  deallocate(Hull)

  allocate(s(5))
  Vec(1,:) = 1.0-2.0*s
  Vec(2,:) = 3.232*s+3.2
  call convex_hull(Vec,Hull)
  fail = size(Hull,2).ne.2
  call report_test("[LineHull]", fail, .false., "Convex hull should have 2&
       & points")  
  deallocate(s)
  deallocate(Hull)

  Vec(1,:) = 1.0
  Vec(2,:) = 1.0
  call convex_hull(Vec,Hull)
  fail = size(Hull,2).ne.1
  call report_test("[ODHull]", fail, .false., "Convex hull should have 1&
       & point")  
  deallocate(Vec)
  deallocate(Hull)

  allocate(Vec(2,8))
  Vec(1,:) = (/0.,1.,0.,1.,0.,1.,0.9,1./)
  Vec(2,:) = (/0.,0.,1.,1.,1.,1.,0.8,0./)
  call convex_hull(Vec,Hull)
  fail = size(Hull,2).ne.4
  call report_test("[SquareDuplicatesHull]", fail, .false., "Convex hull sho&
       &uld have 4 points")  

  fail = .not.outside_hull((/1.5,1.0/),hull)
  call report_test("[OutsideHull]", fail, .false., "Point is outsid&
       &e hull.")  

  fail = outside_hull((/0.5,0.99/),hull)
  call report_test("[InsideHull]", fail, .false., "Point is inside hull.")  

  deallocate(Hull)

  allocate(Hull(2,1))
  Hull = 0.
  Hull(2,1) = 2.

  fail = .not.outside_hull((/0.5,0.99/),hull)
  call report_test("[Outside0DHull]", fail, .false., "Point is outside hull.")  

  fail = outside_hull((/0.0,2.0/),hull)
  call report_test("[Inside0DHull]", fail, .false., "Point is inside hull.")  

  deallocate(Hull)

  allocate(Hull(2,2))
  !y = 1 - 0.4x
  Hull(:,1) = (/0.,1./)
  Hull(:,2) = (/0.5,0.8/)

  fail = .not.outside_hull((/0.3,0.99/),hull)
  call report_test("[NotColinear1DHull]", fail, .false., "Point is outside h&
       &ull.")

  ! 1 - 0.4*0.3 = 0.88
  fail = outside_hull((/0.3,0.88/),hull)
  call report_test("[InsideColinear1DHull]", fail, .false., "Point is inside hull.")

  ! 1 - 0.4*-0.1 = 1.04
  fail = .not.outside_hull((/-0.1,1.04/),hull)
  call report_test("[OutsideBelowColinear1DHull]", fail, .false., "Point is outside hull.")

  ! 1 - 0.4*0.55 = 0.78
  fail = .not.outside_hull((/0.55,0.78/),hull)
  call report_test("[OutsideAboveColinear1DHull]", fail, .false., "Point is outside hull.")

  deallocate(Hull)

  allocate(Hull(2,4))
  hull(1,:) = (/0.0,1.0,1.0,0.0/)
  hull(2,:) = (/0.0,0.0,1.0,1.0/)
  fail = outside_hull((/0.8,0.8/),hull)
  call report_test("[InsideSquareHull]", fail, .false., "Point is inside hull.")
  deallocate(hull)
end subroutine test_convex_hull2d
