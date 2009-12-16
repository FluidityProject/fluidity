subroutine test_get_rotation_matrix

  use quadrature
  use metric_tools
  use unittest_tools
  use vector_tools
  implicit none

  real, dimension(3) :: a3, b3, out3
  real, dimension(3, 3) :: mat3

  real, dimension(2) :: a2, b2, out2
  real, dimension(2, 2) :: mat2

  logical :: fail
  integer :: i, j
  character(len=20) :: buf

  do i=1,5
    write(buf,'(i0)') i

    a3 = random_vector(3) ; a3 = a3 / norm(a3)
    b3 = random_vector(3) ; b3 = b3 / norm(b3)
    mat3 = get_rotation_matrix(a3, b3)

    out3 = matmul(mat3, a3)

    fail = .false.
    do j=1,3
      if (.not. fequals(out3(j), b3(j))) fail = .true.
    end do
    call report_test("[rotation matrix 3d]", fail, .false., "If the rotation matrix &
                    & is supposed to map a -> b, it had better map a -> b.")

    fail = .false.
    if (det(mat3) /= 1.0) fail = .false.
    call report_test("[rotation matrix 3d det]", fail, .false., "Rotation matrices &
                    & have unitary determinant.")
  end do

  do i=1,5
    write(buf,'(i0)') i

    a2 = random_vector(2) ; a2 = a2 / norm(a2)
    b2 = random_vector(2) ; b2 = b2 / norm(b2)
    mat2 = get_rotation_matrix(a2, b2)

    out2 = matmul(mat2, a2)

    fail = .false.
    do j=1,2
      if (.not. fequals(out2(j), b2(j))) fail = .true.
    end do
    call report_test("[rotation matrix 2d]", fail, .false., "If the rotation matrix &
                    & is supposed to map a -> b, it had better map a -> b.")

    fail = .false.
    if (det(mat2) /= 1.0) fail = .false.
    call report_test("[rotation matrix 2d det]", fail, .false., "Rotation matrices &
                    & have unitary determinant.")
  end do

  a3 = (/1.0, 0.0, 0.0/)
  b3 = a3
  mat3 = get_rotation_matrix(a3, b3)

  fail = .false.
  if (is_nan(mat3(1,1))) fail = .true.
  call report_test("[rotation matrix NaN]", fail, .false., "Rotation matrices &
                  & should not return NaN.")

  a3 = (/0.647603843384180, 0.605624220172685, -0.462416009643120/)
  b3 = -1 * a3
  mat3 = get_rotation_matrix(a3, b3)

  fail = .false.
  if (is_nan(mat3(1,1))) fail = .true.
  call report_test("[rotation matrix NaN]", fail, .false., "Rotation matrices &
                  & should not return NaN.")

  a3 = (/0.98686556700264593811E+00, -0.16149555124340336798E+00, 0.39420290632705429906E-02/)
  b3 = (/0.98686556608677455937E+00, -0.16149555683264460448E+00, 0.39420293687526817075E-02/)
  mat3 = get_rotation_matrix(a3, b3)
  
  fail = .false.
  if (is_nan(mat3(1,1))) fail = .true.
  call report_test("[rotation matrix NaN]", fail, .false., "Rotation matrices &
                  & should not return NaN.")

  a3 = (/0.0, 1.0, 0.0/)
  b3 = (/0.0, 1.0/sqrt(2.0), 1.0/sqrt(2.0)/)
  mat3 = get_rotation_matrix(a3, b3)
  
  fail = .false.
  if (det(mat3) .fne. 1.0) fail = .true.
  call report_test("[rotation matrix determinant]", fail, .false., "Rotation &
                  & matrices have determinant one.")

end subroutine test_get_rotation_matrix
