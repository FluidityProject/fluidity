subroutine test_get_deformation_matrix

  use merge_tensors
  use vector_tools
  use unittest_tools
  implicit none

  logical :: fail, warn
  real, dimension(3, 3) :: F, Finv, M, I, I2, tmp
  integer :: j
  character(len=20) :: buf

  ! I2 is the identity matrix.
  I2 = get_matrix_identity(3)
  
  fail = .false.; warn = .false.
  F = get_deformation_matrix(I2)
  tmp = F - I2
  if (.not. mat_zero(tmp)) fail = .true.
  call report_test("[get_deformation_matrix identity]", fail, warn, "F(I) == I")
  
  do j=1,5
    fail = .false.; warn = .false.
    write(buf,'(i0)') j
    
    M = random_posdef_matrix(3)
    
    F = get_deformation_matrix(M)
    
    Finv = F
    call invert(Finv)

    ! M should == F^T * F
    tmp = matmul(transpose(F), F)
    tmp = tmp - M ! should be 0
    if (.not. mat_zero(tmp)) fail = .true.
    call report_test("[get_deformation_matrix " // trim(buf) // " equality]", fail, warn, "M should equal F^T * F")
    
    ! F^-T * M * F^-1 should be the identity.
    I = matmul(matmul(transpose(Finv), M), Finv)
    tmp = I - I2
    if (.not. mat_zero(tmp)) fail = .true.
    call report_test("[get_deformation_matrix " // trim(buf) // " inverse]", fail, warn, "F^-T * M * F^-1 should be the identity")
  end do

end subroutine test_get_deformation_matrix
