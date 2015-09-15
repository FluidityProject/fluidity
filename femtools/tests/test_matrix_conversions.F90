#include "fdebug.h"
subroutine test_matrix_conversions

  use sparse_tools
  use petsc_tools
  use unittest_tools
#ifdef HAVE_PETSC_MODULES
  use petsc 
#endif
  implicit none
  
#include "petsc_legacy.h"

  Mat:: M
  type(dynamic_csr_matrix):: R, S
  type(csr_matrix):: A, B
  logical fail

  R=random_sparse_matrix(99, 100, 1001)
  
  A=dcsr2csr(R)
  
  S=csr2dcsr(A)
  
  fail= .not. fequals(R, S, 1e-8)

  call report_test("[dcsr2csr2dcsr]", fail, .false., &
    "Converting from dcsr_matrix to csr_matrix and back failed.")
    
  M=csr2petsc(A)
  
  B=petsc2csr(M)
  
  S=csr2dcsr(B)
  
  fail= .not. fequals(R, S, 1e-8)
  
  call report_test("[csr2petsc2csr]", fail, .false., &
    "Converting from csr_matrix to PETSc Mat and back failed.")

end subroutine test_matrix_conversions

