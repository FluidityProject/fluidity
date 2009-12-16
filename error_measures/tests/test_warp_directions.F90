subroutine test_warp_directions

  use gradation_metric
  use metric_tools
  use unittest_tools
  use vector_tools
  implicit none

  real, dimension(3, 3) :: P, Q, vec_P, vec_Q
  real, dimension(3) :: val_P, val_Q
  logical :: changed_P, changed_Q
  real :: dist

  logical :: fail

  vec_P(1, :) = (/0.10000000000000000000E+01, 0.00000000000000000000E+00, 0.00000000000000000000E+00/)
  vec_P(2, :) = (/0.00000000000000000000E+00, 0.10000000000000000000E+01, 0.00000000000000000000E+00/)
  vec_P(3, :) = (/0.00000000000000000000E+00, 0.00000000000000000000E+00, 0.10000000000000000000E+01/)
  val_P = (/0.16000000000000000000E+02, 0.16000000000000000000E+02, 0.10000000000000000000E+05/)

  vec_Q(1, :) = (/0.10000000000000000000E+01, 0.00000000000000000000E+00, 0.00000000000000000000E+00/)
  vec_Q(2, :) = (/0.00000000000000000000E+00, 0.00000000000000000000E+00, 0.10000000000000000000E+01/)
  vec_Q(3, :) = (/0.00000000000000000000E+00, 0.10000000000000000000E+01, 0.00000000000000000000E+00/)
  val_Q = (/0.16000000000000000000E+02, 0.16000000000000000000E+02, 0.10000000000000000000E+05/)

  dist = 0.320156213733379

  P = vec_P ; Q = vec_Q
  changed_P = .false. ; changed_Q = .false.
  call warp_directions(vec_P, val_P, changed_P, vec_Q, val_Q, changed_Q, dist)

  fail = .false.
  if (changed_P) then
    write(0,*) "P has changed!"
    fail = .true.
  end if

  if (changed_Q) then
    write(0,*) "Q has changed!"
    fail = .true.
  end if

  call report_test("[warp_directions spheroids]", fail, .false., &
                   "Warping these matrices should not change Q.")

end subroutine test_warp_directions
