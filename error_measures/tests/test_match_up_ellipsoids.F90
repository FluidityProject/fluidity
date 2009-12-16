subroutine test_match_up_ellipsoids
  
  use gradation_metric
  use unittest_tools
  implicit none

  real, dimension(3, 3) :: vec_P, vec_Q
  real, dimension(3) :: val_P, val_Q
  integer, dimension(3) :: perm_P, perm_Q
  logical :: fail

  vec_P(1, :) = (/0.00000000000000000000E+00, 0.36972895530693085375E-01, -0.99931626875382983943E+00/)
  vec_P(2, :) = (/0.00000000000000000000E+00, -0.99931626875382983943E+00, -0.36972895530693078436E-01/)
  vec_P(3, :) = (/0.10000000000000000000E+01, 0.00000000000000000000E+00, 0.00000000000000000000E+00/)

  vec_Q(1, :) = (/0.36972895530693432320E-01, 0.00000000000000000000E+00, -0.99931626875382983943E+00/)
  vec_Q(2, :) = (/-0.99931626875382972841E+00, 0.00000000000000000000E+00, -0.36972895530693432320E-01/)
  vec_Q(3, :) = (/0.00000000000000000000E+00, 0.10000000000000000000E+01, 0.00000000000000000000E+00/)

  val_P = (/0.16000000000000000000E+02, 0.16000000000000003553E+02, 0.10000000000000000000E+03/)
  val_Q = (/0.15999999999999998224E+02, 0.16000000000000000000E+02, 0.10000000000000000000E+03/)

  call match_up_ellipsoids(vec_P, val_P, perm_P, vec_Q, val_Q, perm_Q)

  fail = .false.
  if (perm_Q(2) == perm_Q(3)) fail = .true.
  call report_test("[match up ellipsoids]", fail, .false., &
  "Previous bugs should not happen again.")

end subroutine test_match_up_ellipsoids
