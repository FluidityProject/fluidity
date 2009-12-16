#include "fdebug.h"
module gallopede_block_solvers

  use sparse_tools
  use fldebug
  use global_parameters_gallopede

contains

  subroutine gallopede_block_solve(X_arr, A_csr, &
       b_arr, ksp_type, pc_type, errtol, noi, zero)
    real, dimension(:), intent(inout) :: X_arr
    type(block_csr_matrix), intent(in) :: A_csr
    type(csr_matrix) :: A11, A12, A21, A22
    real, dimension(:), intent(in) :: b_arr
    real, intent(in), optional :: errtol
    integer, intent(in), optional :: noi !number of iterations
    logical, intent(in), optional :: zero

#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"

    KSPType, intent(in) :: ksp_type
    PCType, intent(in) :: pc_type

    !locals 
    integer :: idx, locrow, loccol, glorow, glocol
    integer, dimension(:), allocatable :: loc, nnv ! nnv == # of nonzeroes per row
    PetscErrorCode :: ierr
    Vec :: x, b ! Approximate solution vector, RHS
    Mat :: A
    PetscInt :: n, i, iterations
    KSP :: ksp ! Krylov subspace context
    PC  :: pc  ! Preconditioner context
    KSPConvergedReason :: reason
    PetscReal :: norm

    real :: lerrtol
    integer :: lnoi

    if(present(errtol)) then
       lerrtol = errtol
    else
       lerrtol = PETSC_DEFAULT_DOUBLE_PRECISION
    end if

    if(present(noi)) then
       lnoi = noi
    else
       lnoi = PETSC_DEFAULT_INTEGER
    end if

    !get blocks
    A11 = block(A_csr,1,1)
    A21 = block(A_csr,2,1)
    A12 = block(A_csr,1,2)      
    A22 = block(A_csr,2,2)

    ierr = 0

    n = size(b_arr)
    allocate(loc(n), nnv(0:n-1))
    CHECK(n)

    ! Create the vectors.

    MSG("Creating vector b")
    call VecCreateSeq(MPI_COMM_SELF, n, b, ierr)
    MSG("Creating vector x")
    call VecCreateSeq(MPI_COMM_SELF, n, x, ierr)
    MSG("setting sizes")
    call VecSetSizes(b, n, n, ierr)
    call VecSetSizes(x, n, n, ierr)

    ! Convert from normal fortran vector to PETSc vector.
    ! PETSc vectors are indexed from zero, which is why
    ! you have to create this index array to tell PETSc
    ! where to put the data.

    MSG("Making list")
    do i=1,n
       loc(i) = i - 1
    end do

    ! VecSetValues is conceptually the same as
    ! do i=1,n
    !    b(loc(i)) = dif_vec(i)
    ! end do

    MSG("Setting values in vecs")
    call VecSetValues(b, n, loc, b_arr, INSERT_VALUES, ierr)
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)

    call VecSetValues(X, n, loc, X_arr, INSERT_VALUES, ierr)
    call VecAssemblyBegin(X, ierr)
    call VecAssemblyEnd(X, ierr)

    ! Now set up the matrix. First get the size.

    !matrix is a block 
    locrow =  2*size(A11,1)
    loccol = locrow

    ! compute the sparsity pattern & construct the matrix

    !n is length of b
    do i=0,n/2-1
       nnv(i) = 2*row_length(A11, i+1)
    end do
    nnv(n/2:n-1) = nnv(0:n/2-1)

    call MatCreateSeqAIJ(MPI_COMM_SELF, locrow, loccol, 0, nnv, A, ierr)
    call MatZeroEntries(A, ierr)
    do i=0,n/2-1
       !A11
       call MatSetValues(A, 1, i, row_length(A11, i+1), &
            row_m_ptr(A11, i+1) - 1, &
            & row_val_ptr(A11, i+1), INSERT_VALUES, ierr)
       !A12
       call MatSetValues(A, 1, i, row_length(A12, i+1), &
            row_m_ptr(A12, i+1) - 1 + n/2, &
            & row_val_ptr(A12, i+1), INSERT_VALUES, ierr)
       !A21
       call MatSetValues(A, 1, i + n/2, row_length(A21, i+1), &
            row_m_ptr(A21, i+1) - 1, &
            & row_val_ptr(A21, i+1), INSERT_VALUES, ierr)
       !A22
       call MatSetValues(A, 1, i + n/2, row_length(A22, i+1), &
            row_m_ptr(A22, i+1) - 1 + n/2, &
            & row_val_ptr(A22, i+1), INSERT_VALUES, ierr)
    end do

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

    call MatNorm(A,NORM_1,norm);
    CHECK(norm)

    ! Set up the solver context
    CHECK(ksp_type)
    CHECK(pc_type)

    call KSPCreate(MPI_COMM_SELF, ksp, ierr); CHKERRQ(ierr)
    call KSPSetType(ksp, ksp_type, ierr); CHKERRQ(ierr)
    call KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN, ierr); CHKERRQ(ierr)
    call KSPGetPC(ksp, pc, ierr); CHKERRQ(ierr)
    call PCSetType(pc, pc_type, ierr); CHKERRQ(ierr)
    call KSPSetTolerances(ksp, lerrtol, lerrtol,     &
         PETSC_DEFAULT_DOUBLE_PRECISION, lnoi, ierr); CHKERRQ(ierr)

    CHECK(lerrtol)

    ! Solve the damn thing already

    call KSPSolve(ksp, b, x, ierr)
    call KSPGetConvergedReason(ksp, reason, ierr)
    call KSPGetIterationNumber(ksp, iterations, ierr)

    CHECK(reason)
    CHECK(iterations)

    call VecGetValues(x, n, loc, X_arr, ierr)
    !write (0,*) "X_arr == ", X_arr

    ! Clean up

    call VecDestroy(x, ierr)
    call VecDestroy(b, ierr)
    call MatDestroy(A, ierr)
    call KSPDestroy(ksp, ierr)

    deallocate(loc, nnv)

  end subroutine gallopede_block_solve

end module gallopede_block_solvers
