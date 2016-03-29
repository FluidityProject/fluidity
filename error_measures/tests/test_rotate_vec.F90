subroutine test_rotate_vec

  use gradation_metric
  use metric_tools, only: get_angle
  use unittest_tools
  implicit none

  real, dimension(3, 3) :: vec_P, vec_Q
  integer, dimension(3) :: perm_P, perm_Q
  integer :: idx
  real :: angle
  real :: in_angle, out_angle
  logical :: fail

  vec_P(1, :) = (/0.10000000000000000000E+01, 0.00000000000000000000E+00, 0.00000000000000000000E+00/)
  vec_P(2, :) = (/0.00000000000000000000E+00, 0.00000000000000000000E+00, 0.10000000000000000000E+01/)
  vec_P(3, :) = (/0.00000000000000000000E+00, 0.10000000000000000000E+01, 0.00000000000000000000E+00/)

  vec_Q(1, :) = (/0.00000000000000000000E+00, 0.55115565036679803335E+00, -0.83440245030126314330E+00/)
  vec_Q(2, :) = (/0.00000000000000000000E+00, -0.83440245030126314330E+00, -0.55115565036679803335E+00/)
  vec_Q(3, :) = (/0.10000000000000000000E+01, 0.00000000000000000000E+00, 0.00000000000000000000E+00/)

  perm_P = (/3, 2, 1/)
  perm_Q = (/3, 1, 2/)

  idx = 1
  angle = 0.351727550700944

  in_angle = get_angle(vec_P(:, 3), vec_Q(:, 3))

  call rotate_vec(vec_Q, perm_Q, vec_P, perm_P, idx, -1 * angle)

  out_angle = get_angle(vec_P(:, 3), vec_Q(:, 3))

  fail = .false.
  if (out_angle .fgt. in_angle) then
    fail = .true.
    write(0,*) "in_angle == ", in_angle
    write(0,*) "out_angle == ", out_angle
  end if
  call report_test("[rotate vec]", fail, .false., "Rotating vectors is supposed to bring them closer together!")

end subroutine test_rotate_vec
