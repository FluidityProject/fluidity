subroutine test_rotation_to_vertical
  !!< This subroutine ensures that the rotate_to_vertical function behaves
  !!< itself. 
  use coordinates
  use global_parameters, only: isphere
  use unittest_tools
  implicit none

  !! The rotation matrix
  real, dimension(3,3) :: R
  !! The identity matrix
  real, dimension(3,3) :: I
  !! The local up vector from which to calculate.
  real, dimension(3) :: up
  !! Flags for failure and warning.
  logical :: fail, warn

  fail=.false.
  warn=.false.

  I=reshape((/ 1.0, 0.0, 0.0, &
       &       0.0, 1.0, 0.0, &
       &       0.0, 0.0, 1.0/), (/3,3/))

  up=(/1.0, 1.0, 1.0/)
  
  ! Set isphere=1 so we are on a sphere.
  isphere=1

  R=rotation_to_vertical(up)

  fail=maxval(abs(matmul(transpose(R),R)-I))>100*epsilon(0.0)

  call report_test("[rotation matrix is unitary]", fail, warn, "Rotation mat&
       &rix is not unitary.")
  
  fail=maxval(abs(matmul(R, up)-(/0.0, 0.0, sqrt(3.0)/))) >100*epsilon(0.0)
  
  call report_test("[rotate up vector to vertical]", fail, warn, "Rotation matrix does&
       & not rotate up vector to vertical.")

end subroutine test_rotation_to_vertical
