subroutine test_pseudoinverse

  use vector_tools
  use unittest_tools

  implicit none

  real, dimension(4,5) :: A, B 
  real, dimension(5,4) :: A_pinv

  integer :: i,j
  logical :: fail

  A=reshape((/1,0,0,0,0,0,0,4,0,3,0,0,0,0,0,0,2,0,0,0/),(/4,5/))

  ! Calculate A*
  A_pinv=pseudoinverse(A)

  ! check that (A*)*=A
  B=pseudoinverse(A_pinv)
  fail=any(abs(B-A)>1e-10)
  call report_test("[pseudoinverse]", fail, .false., "(A*)*=A")

  ! check that (AA*)A=A
  fail=any(abs(matmul(A,matmul(A_pinv,A))-A)>1e-10)
  call report_test("[pseudoinverse]", fail, .false., "(AA*)A=A")
  
  ! check that (A*A)A*=A*
  fail=any(abs(matmul(A_pinv,matmul(A,A_pinv))-A_pinv)>1e-10)
  call report_test("[pseudoinverse]", fail, .false., "(A*A)A*=A*")

  ! check that (AA*)^T=(AA*)
  fail=any(abs(transpose(matmul(A,A_pinv))-matmul(A,A_pinv))>1e-10)
  call report_test("[pseudoinverse]", fail, .false., "(AA*)^T=(AA*)")

  ! check that (A*A)^T=(A*A)
  fail=any(abs(transpose(matmul(A_pinv,A))-matmul(A_pinv,A))>1e-10)
  call report_test("[pseudoinverse]", fail, .false., "(A*A)^T=(A*A)")

end subroutine test_pseudoinverse
