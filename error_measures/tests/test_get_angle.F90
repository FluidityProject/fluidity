subroutine test_get_angle

  use metric_tools
  use unittest_tools
  implicit none

  real, dimension(3) :: vecA, vecB
  real :: angle, angle2, pi
  logical :: fail
  integer :: i
  character(len=20) :: buf

  pi = 4.0*atan(1.0)

  vecA = (/1.0, 0.0, 0.0/)
  vecB = (/0.0, 1.0, 0.0/)
  angle = get_angle(vecA, vecB)

  fail = .false.
  if (.not. fequals(angle, pi / 2.0)) fail = .true.
  call report_test("[get angle orthogonal]", fail, .false., "These vectors are orthogonal.")

  vecA = (/1.0, 0.0, 0.0/)
  vecB = (/1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0/)
  angle = get_angle(vecA, vecB)

  fail = .false.
  if (.not. fequals(angle, pi / 4.0)) fail = .true.
  call report_test("[get angle 45 degrees]", fail, .false., "These vectors are at 45 degrees.")

  ! We want get_angle to ignore the sign of the vector, you see.

  vecA = (/1.0, 0.0, 0.0/)
  vecB = (/-1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0/)
  angle = get_angle(vecA, vecB)

  fail = .false.
  if (.not. fequals(angle, pi / 4.0)) fail = .true.
  call report_test("[get angle 45 degrees]", fail, .false., "These vectors are at 45 degrees.")

  do i=1,5
    write(buf, '(i0)') i

    vecA = random_vector(3)
    vecB = random_vector(3)
    angle = get_angle(vecA, vecB)
    angle2 = get_angle(vecA, -1 * vecB)

    fail = .false.
    if (.not. fequals(angle, angle2)) fail = .true.
    call report_test("[get angle ignores sign " // trim(buf) // "]", &
      fail, .false., "The angle computation should lie within [0, Pi/2].")

    angle  = get_angle(vecA, vecB)
    angle2 = get_angle(vecB, vecA)
    fail = .false.
    if (.not. fequals(angle, angle2)) fail = .true.
    call report_test("[get angle is commutative " // trim(buf) // "]", &
      fail, .false., "The angle computation should be commutative.")
  end do


  vecA = (/0.0, 1.0, 0.0/)
  vecB = (/-0.10398997059599775217E+00, -0.99457834584081084017E+00, 0.0/)
  angle = get_angle(vecA, vecB)

  fail = .false.
  if (is_nan(angle)) fail = .true.
  call report_test("[get angle nan]", fail, .false., "Angle computation should never return NaN.")

  vecA = (/-0.24842986328315525002E+00, 0.96862714938199356851E+00, -0.66369050939521093482E-02/)
  vecB = (/-0.24842986326122473706E+00, 0.96862714929648652262E+00, -0.66369183942164558174E-02/)
  angle = get_angle(vecA, vecB)

  fail = .false.
  if (is_nan(angle)) fail = .true.
  call report_test("[get angle nan]", fail, .false., "Angle computation should never return NaN.")

end subroutine test_get_angle
