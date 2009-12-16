#include "fdebug.h"
module gallopede_solvers

  use sparse_tools
  use global_parameters_gallopede
  use fldebug
  use data_structures
  use multilayer_tools

  implicit none

  public :: gallopede_solve, gallopede_block_solve, gallopede_solve_lift, gallopede_block_solve_lift, Initialize_Petsc

  private  
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscviewer.h"

  integer :: solve_counter = 0
  
  contains

 subroutine Initialize_Petsc()
    PetscErrorCode :: ierr
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    CHKERRQ(ierr)
  end subroutine Initialize_Petsc

    subroutine gallopede_solve(X_arr, A_csr, b_arr, ksp_type, &
         pc_type, errtol, noi)
      real, dimension(:), intent(out) :: X_arr
      type(csr_matrix), intent(in) :: A_csr
      real, dimension(:), intent(in) :: b_arr
      double precision, intent(in), optional :: errtol
      integer, intent(in), optional :: noi !number of iterations

      KSPType, intent(in) :: ksp_type
      PCType, intent(in) :: pc_type

      !locals 
      integer :: idx, locrow, loccol, glorow, glocol
      integer, dimension(:), allocatable :: loc, nnv ! nnv == # of nonzeroes per row
      PetscTruth :: pstrue
      PetscErrorCode :: ierr
      Vec :: x, b ! Approximate solution vector, RHS
      Mat :: A
      PetscInt :: n, i, iterations
      KSP :: ksp ! Krylov subspace context
      PC  :: pc  ! Preconditioner context
      KSPConvergedReason :: reason

      real :: lerrtol
      integer :: lnoi

      ewrite(1,*) 'subroutine gallopede_solve'

      !call set_global_debug_level(3)

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

      ierr = 0

      n = size(b_arr)
      allocate(loc(n), nnv(0:n-1))
      ewrite(2,*)(n)

      ! Create the vectors.

      ewrite(2,*)("Creating vector x")
      call VecCreateSeq(PETSC_COMM_WORLD, n, x,ierr)
      CHKERRQ(ierr)
      ewrite(2,*)("Creating vector b")
      call VecCreateSeq(PETSC_COMM_WORLD, n, b,ierr)

      ewrite(2,*)("setting sizes")
      call VecSetSizes(b, n, n, ierr)
      call VecSetSizes(x, n, n, ierr)

      ! Zero the answer.
      ewrite(2,*)("Zeroing answer")
      call VecZeroEntries(x, ierr)

      ! Convert from normal fortran vector to PETSc vector.
      ! PETSc vectors are indexed from zero, which is why
      ! you have to create this index array to tell PETSc
      ! where to put the data.

      ewrite(2,*)("Making list")
      do i=1,n
         loc(i) = i - 1
      end do

      ! VecSetValues is conceptually the same as
      ! do i=1,n
      !    b(loc(i)) = dif_vec(i)
      ! end do

      ewrite(2,*)("Setting values in vecs")
      call VecSetValues(b, n, loc, b_arr, INSERT_VALUES, ierr)
      call VecAssemblyBegin(b, ierr)
      call VecAssemblyEnd(b, ierr)

      ! Now set up the matrix. First get the size.

      ewrite(2,*)("setting up matrix")

      !matrix is a block 
      locrow =  size(A_csr,1)
      loccol = locrow

      ewrite(2,*) "compute the sparsity pattern"

      !n is length of b
      do i=0,n-1
         nnv(i) = row_length(A_csr, i+1)
      end do

      ewrite(2,*) "construct the matrix"

      call MatCreateSeqAIJ(MPI_COMM_SELF, locrow, loccol, 0, nnv, A, ierr)
      call MatZeroEntries(A, ierr)
      do i=0,n-1
         call MatSetValues(A, 1, i, row_length(A_csr, i+1), &
              row_m_ptr(A_csr, i+1) - 1, &
              & row_val_ptr(A_csr, i+1), INSERT_VALUES, ierr)
      end do

      ewrite(2,*) "assemble the matrix"

      call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
      call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

      ewrite(2,*) "Set up the solver context"

      call KSPCreate(MPI_COMM_SELF, ksp, ierr); CHKERRQ(ierr)
      call KSPSetType(ksp, ksp_type, ierr); CHKERRQ(ierr)
      call KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN, ierr); 
      CHKERRQ(ierr)
      call KSPGetPC(ksp, pc, ierr); CHKERRQ(ierr)
      call PCSetType(pc, pc_type, ierr); CHKERRQ(ierr)
      call KSPSetTolerances(ksp, lerrtol, 0.d-20*lerrtol,     &
           PETSC_DEFAULT_DOUBLE_PRECISION, lnoi, ierr); CHKERRQ(ierr)

      ewrite(2,*) "Solve the system"

      call KSPSolve(ksp, b, x, ierr)
      call KSPGetConvergedReason(ksp, reason, ierr)
      call KSPGetIterationNumber(ksp, iterations, ierr)
      ewrite(2,*)(reason)
      ewrite(2,*)(iterations)

      if(reason<0) then
         if (reason==-3) then
            FLAbort("failed to solve matrix in given iterations")
         else
            FLAbort("failed to solve matrix : diverged")
         end if
      end if

      call VecGetValues(x, n, loc, X_arr, ierr)

      ewrite(2,*) "clean up"

      call VecDestroy(x, ierr)
      call VecDestroy(b, ierr)
      call MatDestroy(A, ierr)
      call KSPDestroy(ksp, ierr)

      deallocate(loc, nnv)

!      call set_global_debug_level(1)

      ewrite(1,*) 'END subroutine gallopede_solve'

    end subroutine gallopede_solve

    subroutine gallopede_block_solve(X_arr, A_csr, &
       b_arr, ksp_type, pc_type, errtol, noi, zero)
#include "finclude/petscksp.h"
    real, dimension(:), intent(inout) :: X_arr
    type(block_csr_matrix), intent(in) :: A_csr
    real, dimension(:), intent(in) :: b_arr
    real, intent(in), optional :: errtol
    integer, intent(in), optional :: noi !number of iterations
    logical, intent(in), optional :: zero

    KSPType, intent(in) :: ksp_type
    PCType, intent(in) :: pc_type

    !locals 
    integer :: idx, locrow, loccol, glorow, glocol
    integer, dimension(:), allocatable :: loc, nnv ! nnv == # of nonzeroes per row
    PetscErrorCode :: ierr
    Vec :: x, b ! Approximate solution vector, RHS
    Mat :: A
    PetscInt :: n, i, iterations, nblock, iB, jB
    PetscViewer :: petview
    KSP :: ksp ! Krylov subspace context
    PC  :: pc  ! Preconditioner context
    KSPConvergedReason :: reason
    PetscReal :: norm, emax, emin
    PetscTruth :: flag

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

    ierr = 0

    n = size(b_arr)
    allocate(loc(n), nnv(0:n-1))
    ewrite(2,*)(n)

    ewrite(2,*)(sum(b_arr))

    ! Create the vectors.

    ewrite(2,*)("Creating vector b")
    call VecCreateSeq(MPI_COMM_SELF, n, b, ierr)
    ewrite(2,*)("Creating vector x")
    call VecCreateSeq(MPI_COMM_SELF, n, x, ierr)
    ewrite(2,*)("setting sizes")
    call VecSetSizes(b, n, n, ierr)
    call VecSetSizes(x, n, n, ierr)

    ! Convert from normal fortran vector to PETSc vector.
    ! PETSc vectors are indexed from zero, which is why
    ! you have to create this index array to tell PETSc
    ! where to put the data.

    ewrite(2,*)("Making list")
    do i=1,n
       loc(i) = i - 1
    end do

    ! VecSetValues is conceptually the same as
    ! do i=1,n
    !    b(loc(i)) = dif_vec(i)
    ! end do

    ewrite(2,*)("Setting values in vecs")
    call VecSetValues(b, n, loc, b_arr, INSERT_VALUES, ierr)
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)

    call VecSetValues(X, n, loc, X_arr, INSERT_VALUES, ierr)
    call VecAssemblyBegin(X, ierr)
    call VecAssemblyEnd(X, ierr)

    ! Now set up the matrix. First get the size.

    !matrix is a block 
    locrow = size(A_csr,1)
    loccol = locrow

    ! compute the sparsity pattern & construct the matrix

    nblock = block_size(A_csr,1)
    ewrite(2,*)(nblock)
    ewrite(2,*)(N_vels)

    ewrite(2,*)("Making list")
    !n is length of b
    nnv = 0
    ewrite(2,*)(size(nnv))
    ewrite(2,*)(A_csr%blocks(1)*nblock)

    do iB = 1, A_csr%blocks(1)
       do i=0,nblock-1
          nnv(i+(iB-1)*nblock) = &
               A_csr%blocks(1)*row_length(A_csr,i+1)
       end do
    end do

    if(.false.) then
       do i = 0, n-1
          ewrite(2,*)(nnv(i))
          if(nnv(i)==0) then
             ewrite(2,*)("?")
             ewrite(2,*)(i)
          end if
       end do
       stop
    end if

    ewrite(2,*)("Creating matrix object")
    call MatCreateSeqAIJ(MPI_COMM_SELF, locrow, loccol, 0, nnv, A, ierr)
    ewrite(2,*)("setting matrix entries")
    call MatZeroEntries(A, ierr)
    if(.false.) then
       do i = 0, n-1
          call MatSetValues(A, 1, i,1,i,1.0, INSERT_VALUES, ierr)
       end do
    else
       do iB = 1, A_csr%blocks(1)
          do jB = 1, A_csr%blocks(1)
             do i=0,nblock-1
                call MatSetValues(A, 1, i + (iB-1)*nblock, &
                     row_length(A_csr, i+1), &
                     (jB-1)*nblock + row_m_ptr(A_csr, i+1) - 1, &
                     & row_val_ptr(A_csr, iB, jB, i+1), INSERT_VALUES, ierr)
             end do
          end do
       end do
    end if
    ewrite(2,*)("assembling matrix")
    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
    
    ewrite(2,*)("checking norm")
    call MatNorm(A,NORM_1,norm, ierr);
    ewrite(2,*)(norm)

    ewrite(1,*)("Is matrix symmetric")
    call MatIsSymmetric(A,1d-8,flag)
    ewrite(1,*)(flag==PETSC_TRUE)

!    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"solmat.m",petview,ierr)
!    call PetscViewerSetFormat(petview,PETSC_VIEWER_ASCII_MATLAB,ierr)
!    call MatView(A,petview,ierr)
!    call PetscViewerDestroy(petview,ierr)
    
!    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"rhs.m",petview,ierr)
!    call PetscViewerSetFormat(petview,PETSC_VIEWER_ASCII_MATLAB,ierr)
!    call VecView(b,petview,ierr)
!    call PetscViewerDestroy(petview,ierr)


    ! Set up the solver context
    ewrite(2,*)(ksp_type)
    ewrite(2,*)(pc_type)

    ewrite(2,*)("setting up solver context")

    call KSPCreate(MPI_COMM_SELF, ksp, ierr); CHKERRQ(ierr)
    call KSPSetType(ksp, ksp_type, ierr); CHKERRQ(ierr)
    call KSPSetOperators(ksp, A, A,&
         DIFFERENT_NONZERO_PATTERN, ierr); CHKERRQ(ierr) 
!    call KSPSetComputeSingularValues(ksp,PETSC_TRUE,ierr)
    
!     ewrite(1,*) "extreme Matrix SVs", emax , emin

    call KSPGetPC(ksp, pc, ierr); CHKERRQ(ierr)
    call PCSetType(pc, pc_type, ierr); CHKERRQ(ierr)
    call KSPSetTolerances(ksp, lerrtol,lerrtol,     &
         PETSC_DEFAULT_DOUBLE_PRECISION, lnoi, ierr); CHKERRQ(ierr)

    ewrite(2,*)(lerrtol)

    ewrite(1,*)("solving equation")
!    call KSPMonitorSet(ksp,KSPMonitorSingularValue,PETSC_NULL_OBJECT,&
!         PETSC_NULL_FUNCTION,ierr)
    call KSPSolve(ksp, b, x, ierr)
    call KSPGetConvergedReason(ksp, reason, ierr)
    call KSPGetIterationNumber(ksp, iterations, ierr)
!    call KSPComputeExtremeSingularValues(ksp,emax,emin,ierr)
    
     ewrite(1,*) "extreme Matrix SVs", emax , emin

    ewrite(1,*)(reason)
    ewrite(1,*)(iterations)

    if(reason<0) then
       if (reason==-3) then
          FLAbort("failed to solve matrix in given iterations")
       else
          FLAbort("failed to solve matrix : diverged")
       end if
    end if

    call VecGetValues(x, n, loc, X_arr, ierr)

    ! Clean up

    call VecDestroy(x, ierr)
    call VecDestroy(b, ierr)
    call MatDestroy(A, ierr)
    call KSPDestroy(ksp, ierr)

    deallocate(loc, nnv)

    ewrite(1,*) '  end subroutine gallopede_block_solve'
 
  end subroutine gallopede_block_solve

  subroutine gallopede_solve_lift(x1,x2, &
       A11, A12, A21, A22, &
       b1, b2, &
       bcs, &
       ksp_type, &
       pc_type, errtol, noi)
    real, dimension(:), intent(inout) :: x1, x2
    type(csr_matrix), intent(in) :: A11,A12,A21,A22
    real, dimension(:), intent(in) :: b1, b2
    real, intent(in), optional :: errtol
    integer, intent(in), optional :: noi !number of iterations
    type(bc_info), intent(in) :: bcs

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

    real :: lerrtol
    integer :: lnoi
    integer, dimension(:), pointer :: rpointer => null()
    real, dimension(:), pointer :: rvalpointer => null()
    integer :: rlength
    logical :: flag

    real :: info(MAT_INFO_SIZE), norm
    real, dimension(:), allocatable :: vec

    ewrite(1,*) 'subroutine gallopede_solve_lift'

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

    ierr = 0

    n = 2*bcs%N_interior + bcs%N_tangents
    
    allocate(loc(n), nnv(0:n-1))
    ewrite(2,*)(n)

     ! Create the vectors.

    ewrite(2,*)("Creating vector b")
    call VecCreateSeq(MPI_COMM_SELF, n, b, ierr)
    ewrite(2,*)("Creating vector x")
    call VecCreateSeq(MPI_COMM_SELF, n, x, ierr)
    ewrite(2,*)("setting sizes")
    call VecSetSizes(b, n, n, ierr)
    call VecSetSizes(x, n, n, ierr)

    ! Zero the answer.
    ewrite(2,*)("Zeroing answer")
    call VecZeroEntries(x, ierr)

    ! Convert from normal fortran vector to PETSc vector.
    ! PETSc vectors are indexed from zero, which is why
    ! you have to create this index array to tell PETSc
    ! where to put the data.

    ewrite(2,*)("Making list")
    do i=1,n
       loc(i) = i - 1
    end do

    ! VecSetValues is conceptually the same as
    ! do i=1,n
    !    b(loc(i)) = dif_vec(i)
    ! end do

    ewrite(2,*)("Setting values in vecs")
    ewrite(2,*) 'interior x-cpts'

    call VecSetValues(b, size(bcs%interior_list), &
         loc(1:size(bcs%interior_list)), &
         b1(bcs%interior_list), INSERT_VALUES, ierr)
    call VecSetValues(b, size(bcs%interior_list), &
         size(bcs%interior_list) + loc(1:size(bcs%interior_list)), &
         b2(bcs%interior_list), INSERT_VALUES, ierr)
    if (bcs%N_tangents>0) then
       call VecSetValues(b, size(bcs%tangent_list), &
            2*size(bcs%interior_list) + loc(1:size(bcs%tangent_list)), &
            b1(bcs%tangent_list)*bcs%tangents(1,:) + &
            b2(bcs%tangent_list)*bcs%tangents(2,:) &
            , INSERT_VALUES, ierr)
    end if

    ewrite(2,*) 'assembling vecs'
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)
    CHKERRQ(ierr)

    ! Now set up the matrix. First get the size.

    ewrite(2,*)("setting up matrix")

    ! compute the sparsity pattern & construct the matrix

    ewrite(2,*) 'computing nnv'
    !n is length of b
    do i=0,bcs%N_interior-1
       nnv(i) = get_lifted_nnv(row_m_ptr(A11,bcs%interior_list(i+1)), bcs)
       nnv(bcs%N_interior+i) = nnv(i)
    end do
    ewrite(2,*) size(bcs%tangent_list)
    ewrite(2,*) bcs%N_tangents
    do i=0,bcs%N_tangents-1
       nnv(2*bcs%N_interior+i) = get_lifted_nnv( &
            row_m_ptr(A11,bcs%tangent_list(i+1)), bcs)
    end do

    ewrite(2,*) 'creating matrix'

    call MatCreateSeqAIJ(MPI_COMM_SELF, n, n, 0, nnv, A, ierr)
    CHKERRQ(ierr)
    call MatZeroEntries(A, ierr);CHKERRQ(ierr)

    ewrite(2,*) 'internal nodes'
    do i=0,bcs%N_interior-1
       call internal_row(i+1,bcs,A11,A12,rpointer,rvalpointer)
       call MatSetValues(A,1,i,size(rpointer),rpointer-1, &
            rvalpointer,INSERT_VALUES, &
               ierr)
       CHKERRQ(ierr)
       call internal_row(i+1,bcs,A21,A22,rpointer,rvalpointer)
       call MatSetValues(A,1,i+bcs%N_interior, &
            size(rpointer),rpointer-1,rvalpointer,INSERT_VALUES,ierr)
       CHKERRQ(ierr)
    end do
    ewrite(2,*) 'boundary nodes'
    if (bcs%N_tangents>0) then
       do i=0,bcs%N_tangents-1
          call tangent_row(i+1,bcs,A11,A12,A21,A22,rpointer, &
               rvalpointer)
          call MatSetValues(A,1,i+2*bcs%N_interior, &
               size(rpointer),rpointer-1,rvalpointer,INSERT_VALUES, &
               ierr)
          CHKERRQ(ierr)
       end do
    end if


    ewrite(2,*) 'assembling matrix'

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

    ewrite(2,*) 'Set up the solver context'
    ewrite(2,*) ksp_type
    ewrite(2,*) pc_type

    ewrite(2,*)("checking norm")
    call MatNorm(A,NORM_1,norm, ierr);
    ewrite(2,*)(norm)

    call KSPCreate(MPI_COMM_SELF, ksp, ierr); CHKERRQ(ierr)
    call KSPSetType(ksp, ksp_type, ierr); CHKERRQ(ierr)
    call KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN, ierr); 
    CHKERRQ(ierr)
    call KSPGetPC(ksp, pc, ierr); CHKERRQ(ierr)
    call PCSetType(pc, pc_type, ierr); CHKERRQ(ierr)
    call KSPSetTolerances(ksp, lerrtol, lerrtol,     &
         PETSC_DEFAULT_DOUBLE_PRECISION, lnoi, ierr); CHKERRQ(ierr)

    ewrite(2,*) 'solve the system'

    call KSPSolve(ksp, b, x, ierr)
    call KSPGetConvergedReason(ksp, reason, ierr)
    call KSPGetIterationNumber(ksp, iterations, ierr)
    ewrite(2,*) reason
    ewrite(2,*) iterations

     if(reason<0) then
         if (reason==-3) then
            FLAbort("failed to solve matrix in given iterations")
         else
            FLAbort("failed to solve matrix : diverged")
         end if
      end if

    ewrite(2,*) 'reassemble the solution vectors'

    ewrite(2,*) '1 interior cpts'
    allocate( vec(n) )
    vec = 0
    x1 = 0.0
    x2 = 0.0
    call VecGetValues(x,bcs%n_interior, loc(1:bcs%N_interior), &
         vec(1:bcs%N_interior), ierr)
    x1(bcs%interior_list) = vec(1:bcs%N_interior)
    call VecGetValues(x,bcs%N_interior, &
         loc(bcs%N_interior+1:2*bcs%N_interior), &
         vec(1:bcs%N_interior), ierr)
    x2(bcs%interior_list) = vec(1:bcs%N_interior)
    if (bcs%N_tangents >0) then
       call VecGetValues(x, bcs%N_tangents, &
            loc(2*bcs%N_interior+1:2*bcs%N_interior+bcs%N_tangents), &
            vec(1:bcs%N_tangents), ierr)
       x1(bcs%tangent_list) = vec(1:bcs%N_tangents)
    end if
    deallocate( vec )

    if (bcs%N_tangents >0) then
       x2(bcs%tangent_list) = x1(bcs%tangent_list)
       x1(bcs%tangent_list) = x1(bcs%tangent_list)*bcs%tangents(1,:)
       x2(bcs%tangent_list) = x2(bcs%tangent_list)*bcs%tangents(2,:)
    end if
       

    ewrite(2,*) 'Clean up'

    call VecDestroy(x, ierr)
    call VecDestroy(b, ierr)
    call MatDestroy(A, ierr)
    call KSPDestroy(ksp, ierr)

    deallocate(loc, nnv)

    if(associated(rpointer)) then
       deallocate(rpointer)
       rpointer=>null()
    end if
    if(associated(rvalpointer)) then
       deallocate(rvalpointer)
       rvalpointer=>null()
    end if

    ewrite(1,*) 'END subroutine gallopede_solve_lift'

  end subroutine gallopede_solve_lift

  subroutine gallopede_block_solve_lift(xvec, &
       b1, b2, M,&
       bcs, &
       ksp_type, &
       pc_type, rerrtol, abserrtol, noi, zero)
!    type(layer), dimension(N_Layers), intent(inout) :: x1, x2
    real, dimension(:), intent(inout):: xvec
!    type(layer), dimension(N_Layers), intent(in) :: b1, b2
    real, dimension(:), intent(in) :: b1,b2
    type(block_csr_matrix), intent(in) :: M
    double precision, intent(in), optional :: rerrtol, abserrtol
    integer, intent(in), optional :: noi !number of iterations
    type(bc_info), intent(in) :: bcs
    logical, intent(in), optional :: zero

    KSPType, intent(in) :: ksp_type
    PCType, intent(in) :: pc_type

    !locals 
    logical :: lzero
    integer :: idx, locrow, loccol, glorow, glocol
    integer, dimension(:), allocatable :: loc, nnv ! nnv == # of nonzeroes per row
    PetscErrorCode :: ierr
    Vec :: x, b ! Approximate solution vector, RHS
    Mat :: A
    PetscInt :: n, i, iterations, layer_i,layer_j
    KSP :: ksp ! Krylov subspace context
    PC  :: pc  ! Preconditioner context
    KSPConvergedReason :: reason

    real :: lrerrtol, labserrtol
    integer :: lnoi
    integer, dimension(:), pointer :: rpointer => null()
    real, dimension(:), pointer :: rvalpointer => null()
    integer :: rlength, layer_slot
    logical :: flag

    real :: info(MAT_INFO_SIZE), norm, sum0
    real, dimension(:), allocatable :: vec

    ewrite(1,*) 'subroutine gallopede_block_solve_lift'

    ewrite(1,*) bcs%N_tangents, 'N_tangents = '

    solve_counter = solve_counter + 1

    assert(size(xvec)==2*size(b1))
     assert(size(xvec)==2*size(b2))

!     call set_global_debug_level(1)

    !if(solve_counter==2) then
    !   ewrite(1,*) maxval(M%val)
    !   FLAbort('asdfasdfasdf')
    !end if

    if(present(zero)) then
       lzero = zero
    else
       lzero = .false.
    end if

    if(present(rerrtol)) then 
       lrerrtol = rerrtol
    else
       lrerrtol = PETSC_DEFAULT_DOUBLE_PRECISION
    end if

    if(present(abserrtol)) then
       labserrtol = abserrtol
    else
       labserrtol = PETSC_DEFAULT_DOUBLE_PRECISION
    end if

    if(present(noi)) then
       lnoi = noi
    else
       lnoi = PETSC_DEFAULT_INTEGER
    end if

    ierr = 0

    n = (2*bcs%N_interior + bcs%N_tangents)*N_Layers

    allocate(loc(n), nnv(0:n-1))
    ewrite(2,*)(n)

    ! Create the vectors.

    ewrite(2,*)("Creating vector b")
    call VecCreateSeq(MPI_COMM_SELF, n, b, ierr)
    ewrite(2,*)("Creating vector x")
    call VecCreateSeq(MPI_COMM_SELF, n, x, ierr)
    ewrite(2,*)("setting sizes")
    call VecSetSizes(b, n, n, ierr)
    call VecSetSizes(x, n, n, ierr)

    ! Zero the answer.
    ewrite(2,*)("Zeroing answer")
    call VecZeroEntries(x, ierr)

    ! Convert from normal fortran vector to PETSc vector.
    ! PETSc vectors are indexed from zero, which is why
    ! you have to create this index array to tell PETSc
    ! where to put the data.

    ewrite(2,*)("Making list")
    do i=1,n
       loc(i) = i - 1
    end do

    ! VecSetValues is conceptually the same as
    ! do i=1,n
    !    b(loc(i)) = dif_vec(i)
    ! end do

    ewrite(2,*)("Setting values in vecs")

    do layer_i = 1, N_Layers
       layer_slot = (layer_i-1)*(2*bcs%N_interior+bcs%N_tangents)

       call VecSetValues(b, size(bcs%interior_list), &
            layer_slot+loc(1:size(bcs%interior_list)), &
            b1((layer_i-1)*n_vels+bcs%interior_list), INSERT_VALUES, ierr)
       call VecSetValues(b, size(bcs%interior_list), &
            layer_slot+size(bcs%interior_list)+loc(1:size(bcs%interior_list)),&
            b2((layer_i-1)*n_vels+bcs%interior_list), INSERT_VALUES, ierr)
       if (bcs%N_tangents >0) then
          call VecSetValues(b, size(bcs%tangent_list), &
               layer_slot+2*size(bcs%interior_list)+ &
               loc(1:size(bcs%tangent_list)), &
               b1((layer_i-1)*n_vels+bcs%tangent_list)*bcs%tangents(1,:) + &
               b2((layer_i-1)*n_vels+bcs%tangent_list)*bcs%tangents(2,:), &
               INSERT_VALUES, ierr)
       end if
    end do
    ewrite(2,*) 'assembling vecs'
    call VecAssemblyBegin(b, ierr)
    call VecAssemblyEnd(b, ierr)
    CHKERRQ(ierr)

    if(.not.lzero) then
       do layer_i = 1, N_Layers
          layer_slot = (layer_i-1)*(2*bcs%N_interior+bcs%N_tangents)
          call VecSetValues(x, size(bcs%interior_list), &
               layer_slot+loc(1:size(bcs%interior_list)), &
               xvec((layer_i-1)*2*n_vels+bcs%interior_list),&
               INSERT_VALUES, ierr)
          call VecSetValues(x, size(bcs%interior_list), &
               layer_slot+size(bcs%interior_list)+&
               loc(1:size(bcs%interior_list)),&
               xvec((layer_i-1)*2*n_vels+n_vels+bcs%interior_list),&
               INSERT_VALUES, ierr)
          if (bcs%N_tangents >0) then
             call VecSetValues(x, size(bcs%tangent_list), &
                  layer_slot+2*size(bcs%interior_list)+ &
                  loc(1:size(bcs%tangent_list)), &
                  xvec((layer_i-1)*2*n_vels+bcs%tangent_list)&
                  *bcs%tangents(1,:) + &
                  xvec((layer_i-1)*2*n_vels+n_vels+bcs%tangent_list)&
                  *bcs%tangents(2,:), &
                  INSERT_VALUES, ierr)
          end if
       end do
       ewrite(2,*) 'assembling vecs'
       call VecAssemblyBegin(x, ierr)
       call VecAssemblyEnd(x, ierr)
       CHKERRQ(ierr)
    end if

    ! Now set up the matrix. First get the size.

    ewrite(2,*)("setting up matrix")

    ! compute the sparsity pattern & construct the matrix

    ewrite(2,*) 'computing nnv'
    !n is length of b
    do i=0,bcs%N_interior-1
       nnv(i) = N_Layers* &
            get_lifted_nnv(row_m_ptr(M,bcs%interior_list(i+1)), bcs)
       nnv(bcs%N_interior+i) = nnv(i)
    end do
    ewrite(2,*) size(bcs%tangent_list)
    ewrite(2,*) bcs%N_tangents
    if (bcs%N_tangents >0) then
       do i=0,bcs%N_tangents-1
          nnv(2*bcs%N_interior+i) = N_Layers* &
               get_lifted_nnv(row_m_ptr(M,bcs%tangent_list(i+1)), bcs)
       end do
    end if

    ewrite(2,*) 'n_layers', n_layers
    ewrite(2,*) 'n nodes', 2*bcs%N_interior+bcs%N_tangents
    ewrite(2,*) 'size(nnv)', size(nnv)
    do i=1,N_layers-1
       nnv(i*(2*bcs%N_interior+bcs%N_tangents): &
            (i+1)*(2*bcs%N_interior+bcs%N_tangents)-1) = &
            nnv(0:(2*bcs%N_interior+bcs%N_tangents)-1)
    end do

    ewrite(2,*) 'creating matrix'

    call MatCreateSeqAIJ(MPI_COMM_SELF, n, n, 0, nnv, A, ierr)
    CHKERRQ(ierr)
    call MatZeroEntries(A, ierr);CHKERRQ(ierr)

    norm = 0.

    assert(bcs%N_interior+bcs%N_tangents==N_vels)
!    assert(bcs%N_tangents>0)

    ewrite(2,*) 'internal nodes'
    do i=0,bcs%N_interior-1
       do layer_i = 1, N_Layers
          do layer_j = 1, N_Layers
             call internal_row_block(i+1,1,layer_i,layer_j, &
                  bcs,M,rpointer,rvalpointer)
             sum0 = sum(rvalpointer)
             norm = max(norm,maxval(abs(rvalpointer)))
             assert(maxval(rpointer)<N_Layers*2*N_vels+1)
             assert(minval(rpointer)>0)
             assert(size(rpointer)>0)
             assert(sum(rvalpointer)<huge(0.0))
             call MatSetValues(A,1, &
                  (layer_i-1)*(2*bcs%N_interior+bcs%N_tangents) + i, &
                  size(rpointer),rpointer-1, &
                  rvalpointer,INSERT_VALUES, &
                  ierr)
             CHKERRQ(ierr)
             call internal_row_block(i+1,2,layer_i,layer_j, &
                  bcs,M,rpointer,rvalpointer)
             norm = max(norm,maxval(abs(rvalpointer)))
             sum0 = sum0 + sum(rvalpointer)
!             assert(maxval(rpointer)<N_Layers*2*N_vels)
             assert(minval(rpointer)>0)
             assert(size(rpointer)>0)
             assert(sum(rvalpointer)<huge(0.0))
             call MatSetValues(A,1, &
                  (layer_i-1)*(2*bcs%N_interior+bcs%N_tangents) + &
                  i+bcs%N_interior, &
                  size(rpointer),rpointer-1,rvalpointer,INSERT_VALUES,ierr)
             CHKERRQ(ierr)
          end do
       end do
    end do
    ewrite(2,*) 'boundary nodes'
    if (bcs%N_tangents >0) then
       do i=0,bcs%N_tangents-1
          do layer_i = 1, N_Layers
             do layer_j = 1, N_Layers
                call tangent_row_block(i+1,layer_i,layer_j, &
                     bcs,M,rpointer,rvalpointer)
                norm = max(norm,maxval(abs(rvalpointer)))
                sum0 = sum0 + sum(rvalpointer)
                assert(maxval(rpointer)<N_Layers*2*N_vels)
                assert(minval(rpointer)>0)
                assert(size(rpointer)>0)
                assert(sum(rvalpointer)<huge(0.0))
                call MatSetValues(A,1, &
                  (layer_i-1)*(2*bcs%N_interior+bcs%N_tangents) + &
                  i+2*bcs%N_interior, &
                  size(rpointer),rpointer-1,rvalpointer,INSERT_VALUES, &
                  ierr)
                CHKERRQ(ierr)
             end do
          end do
       end do
    end if
!    assert(norm>0.0)
!    assert(sum0>0.0)

    ewrite(2,*) 'max vals', norm

    ewrite(2,*) 'assembling matrix'

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr)

    ewrite(2,*)("checking 1 norm")
    call MatNorm(A,NORM_1,norm, ierr);
    ewrite(2,*)(norm)
    ewrite(2,*)("checking infinity norm")
    call MatNorm(A,NORM_INFINITY,norm, ierr);
    ewrite(2,*)(norm)

    !   call MatView(A, PETSC_VIEWER_STDERR_WORLD, ierr)

    ewrite(2,*) 'Set up the solver context'
    ewrite(2,*) ksp_type
    ewrite(2,*) pc_type

    call KSPCreate(MPI_COMM_SELF, ksp, ierr); CHKERRQ(ierr)
    call KSPSetType(ksp, ksp_type, ierr); CHKERRQ(ierr)
    call KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN, ierr); 
    CHKERRQ(ierr)
    if(.not.lzero) then
       call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
    end if
    call KSPGetPC(ksp, pc, ierr); CHKERRQ(ierr)
    call PCSetType(pc, pc_type, ierr); CHKERRQ(ierr)
    call KSPSetTolerances(ksp, lrerrtol, labserrtol, &
         PETSC_DEFAULT_DOUBLE_PRECISION, lnoi, ierr); CHKERRQ(ierr)

    ewrite(2,*) 'solve the system'

    call KSPSolve(ksp, b, x, ierr)
    call KSPGetConvergedReason(ksp, reason, ierr)
    call KSPGetIterationNumber(ksp, iterations, ierr)
    ewrite(2,*) reason
    ewrite(2,*) iterations
          
       if(reason<0) then
          if (reason==-3) then
             FLAbort("failed to solve matrix in given iterations")
          else
             FLAbort("failed to solve matrix : diverged")
          end if
       end if
       


  ewrite(2,*) 'reassemble the solution vectors'

    ewrite(2,*) '1 interior cpts'
    allocate( vec(n) )
    vec = 0
!    call clear(xvec)
    xvec=0.
 !   call clear(x2)

    do layer_i = 1,N_Layers

       call VecGetValues(x,bcs%n_interior, &
            (layer_i-1)*(2*bcs%N_interior+bcs%N_tangents) + &
            loc(1:bcs%N_interior), &
            vec(1:bcs%N_interior), ierr)
       xvec((layer_i-1)*2*n_vels+bcs%interior_list) = vec(1:bcs%N_interior)
       call VecGetValues(x,bcs%N_interior, &
            (layer_i-1)*(2*bcs%N_interior+bcs%N_tangents) + &
            loc(bcs%N_interior+1:2*bcs%N_interior), &
            vec(1:bcs%N_interior), ierr)
       xvec((layer_i-1)*2*n_vels+n_vels+bcs%interior_list) = &
            vec(1:bcs%N_interior)
       if ( bcs%N_tangents >0) then   
          call VecGetValues(x, bcs%N_tangents, &
               (layer_i-1)*(2*bcs%N_interior+bcs%N_tangents) + &
               loc(2*bcs%N_interior+1:2*bcs%N_interior+bcs%N_tangents), &
               vec(1:bcs%N_tangents), ierr)
          xvec((layer_i-1)*2*n_vels+bcs%tangent_list) = vec(1:bcs%N_tangents)
          
          xvec((layer_i-1)*2*n_vels+n_vels+bcs%tangent_list) = &
               vec(1:bcs%N_tangents)
          xvec((layer_i-1)*2*n_vels+bcs%tangent_list) = &
               xvec((layer_i-1)*2*n_vels+bcs%tangent_list)*bcs%tangents(1,:)
          xvec((layer_i-1)*2*n_vels+n_vels+bcs%tangent_list) = &
               xvec((layer_i-1)*2*n_vels+n_vels+bcs%tangent_list)&
               *bcs%tangents(2,:)
       end if

    end do
    deallocate( vec )

    ewrite(2,*) 'Clean up'

    call VecDestroy(x, ierr)
    call VecDestroy(b, ierr)
    call MatDestroy(A, ierr)
    call KSPDestroy(ksp, ierr)

    deallocate(loc, nnv)

    if(associated(rpointer)) then
       deallocate(rpointer)
       rpointer=>null()
    end if
    if(associated(rvalpointer)) then
       deallocate(rvalpointer)
       rvalpointer=>null()
    end if

!    call set_global_debug_level(1)

    ewrite(1,*) 'END subroutine gallopede_block_solve_lift'

  end subroutine gallopede_block_solve_lift

  function get_lifted_nnv(row,bcs) result(nnv)
    integer, dimension(:), target, intent(in) :: row
    type(bc_info) :: bcs
    integer :: nnv

    !locals
    integer :: i

    nnv = 0
    do i = 1, size(row)
       if(bcs%bc_marker(row(i)) == 1) nnv = nnv + 2
       if(bcs%bc_marker(row(i)) == 2) nnv = nnv + 1
    end do

  end function get_lifted_nnv

  subroutine internal_row(i,bcs,A1,A2,rpointer,rvalpointer)
    integer, intent(in) :: i
    type(bc_info), intent(in) :: bcs
    type(csr_matrix), intent(in) :: A1,A2
    integer, dimension(:), pointer :: rpointer
    real, dimension(:), pointer :: rvalpointer

    !locals
    integer :: j,ninternals,ntangents,k,rlength
    integer, dimension(:), pointer :: row

    !ewrite(1,*) 'subroutine internal_row'

    if(associated(rpointer)) then
       deallocate(rpointer)
       rpointer=>null()
    end if
    if(associated(rvalpointer)) then
       deallocate(rvalpointer)
       rvalpointer=>null()
    end if

    ninternals=0
    ntangents=0

    !row = row_m_ptr(A1, bcs%interior_list(i))
    row => A1%colm(A1%findrm(bcs%interior_list(i)): &
         A1%findrm(bcs%interior_list(i)+1)-1)
    do k = 1, size(row)
       if(bcs%bc_marker(row(k))==1) ninternals = ninternals + 1
       if(bcs%bc_marker(row(k))==2) ntangents = ntangents + 1
    end do
    rlength = 2*ninternals + ntangents
    allocate(rpointer(rlength),rvalpointer(rlength))    

    j = 1
    do k = 1,size(row)
       select case(bcs%bc_marker(row(k)))
       case (1)
          rpointer(j) = bcs%lifted_ordering(row(k))
          rvalpointer(j) = val(A1,bcs%interior_list(i),row(k))
          rpointer(j+1) = bcs%N_interior + bcs%lifted_ordering(row(k))
          rvalpointer(j+1) = val(A2,bcs%interior_list(i),row(k))
          j = j+2
       case (2)
          rpointer(j) = 2*bcs%N_interior + bcs%lifted_ordering(row(k))
          rvalpointer(j) = &
               val(A1,bcs%interior_list(i),row(k)) * &
               bcs%tangents(1,bcs%lifted_ordering(row(k))) &
               + val(A2,bcs%interior_list(i),row(k)) * &
               bcs%tangents(2,bcs%lifted_ordering(row(k)))
               j = j+1
       end select
          
    end do

    !ewrite(1,*) 'END subroutine internal_row'

  end subroutine internal_row

  subroutine internal_row_block(i,cartesian_cpt,layer_i,layer_j, &
       bcs,M,rpointer,rvalpointer)
    integer, intent(in) :: i, cartesian_cpt,layer_i,layer_j
    type(bc_info), intent(in) :: bcs
    type(block_csr_matrix), intent(in) :: M
    integer, dimension(:), pointer :: rpointer
    real, dimension(:), pointer :: rvalpointer

    !locals
    integer :: j,ninternals,ntangents,k,rlength
    integer, dimension(:), pointer :: row
    real, dimension(:), pointer :: M1val,M2val

    !ewrite(1,*) 'subroutine internal_row'

    if(associated(rpointer)) then
       deallocate(rpointer)
       rpointer=>null()
    end if
    if(associated(rvalpointer)) then
       deallocate(rvalpointer)
       rvalpointer=>null()
    end if

    ninternals=0
    ntangents=0

    row => M%colm(M%findrm(bcs%interior_list(i)): &
         M%findrm(bcs%interior_list(i)+1)-1)
    M1val => row_val_ptr(M,(layer_i-1)*2 + cartesian_cpt, &
         (layer_j-1)*2 + 1, bcs%interior_list(i))
    M2val => row_val_ptr(M,(layer_i-1)*2 + cartesian_cpt, &
         (layer_j-1)*2 + 2, bcs%interior_list(i))

    do k = 1, size(row)
       if(bcs%bc_marker(row(k))==1) ninternals = ninternals + 1
       if(bcs%bc_marker(row(k))==2) ntangents = ntangents + 1
    end do
    rlength = 2*ninternals + ntangents
    allocate(rpointer(rlength),rvalpointer(rlength))
    
    j = 1
    do k = 1,size(row)
       select case(bcs%bc_marker(row(k)))
       case (1)
          rpointer(j) = (layer_j-1)*(2*bcs%N_interior+bcs%N_tangents) + &
               bcs%lifted_ordering(row(k))
          rvalpointer(j) = m1val(k)
          rpointer(j+1) = (layer_j-1)*(2*bcs%N_interior+bcs%N_tangents) + & 
               bcs%N_interior + bcs%lifted_ordering(row(k))
          rvalpointer(j+1) = m2val(k)
          j = j+2
       case (2)
          rpointer(j) = (layer_j-1)*(2*bcs%N_interior+bcs%N_tangents)& 
               + 2*bcs%N_interior + bcs%lifted_ordering(row(k))
          rvalpointer(j) = m1val(k) * &
               bcs%tangents(1,bcs%lifted_ordering(row(k))) &
               + m2val(k) * &
               bcs%tangents(2,bcs%lifted_ordering(row(k)))
          j = j+1
       end select
    end do

  end subroutine internal_row_block

  subroutine tangent_row(i,bcs,A11,A12,A21,A22,rpointer, &
       rvalpointer)
    integer, intent(in) :: i
    type(bc_info), intent(in) :: bcs
    type(csr_matrix), intent(in) :: A11,A12,A21,A22
    integer, dimension(:), pointer :: rpointer
    real, dimension(:), pointer :: rvalpointer

    !locals
    integer :: j,ninternals,ntangents,k,rlength
    integer, dimension(:), pointer :: row
    real, dimension(:), pointer :: rowval

    !ewrite(1,*) '  subroutine tangent_row'

    if(associated(rpointer)) then
       deallocate(rpointer)
       rpointer=>null()
    end if
    if(associated(rvalpointer)) then
       deallocate(rvalpointer)
       rvalpointer=>null()
    end if
        
    ninternals=0
    ntangents=0

    row => A11%colm(A11%findrm(bcs%tangent_list(i)): &
         A11%findrm(bcs%tangent_list(i)+1)-1)
    do k = 1, size(row)
       if(bcs%bc_marker(row(k))==1) ninternals = ninternals + 1
       if(bcs%bc_marker(row(k))==2) ntangents = ntangents + 1
    end do
    rlength = 2*ninternals + ntangents
    allocate(rpointer(rlength),rvalpointer(rlength))

    j = 1
    do k = 1,row_length(A11,bcs%tangent_list(i))
       select case(bcs%bc_marker(row(k)))
       case (1)
          rpointer(j) = bcs%lifted_ordering(row(k))
          rvalpointer(j) = val(A11,bcs%tangent_list(i),row(k)) * &
               bcs%tangents(1,i) + &
               val(A21,bcs%tangent_list(i),row(k)) * &
               bcs%tangents(2,i)
          rpointer(j+1) = bcs%N_interior + bcs%lifted_ordering(row(k))
          rvalpointer(j+1) = val(A12,bcs%tangent_list(i),row(k)) * &
               bcs%tangents(1,i) + &
               val(A22,bcs%tangent_list(i),row(k)) * &
               bcs%tangents(2,i)
          j = j+2
       case (2)
          rpointer(j) = 2*bcs%N_interior + bcs%lifted_ordering(row(k))
          rvalpointer(j) = &
               val(A11,bcs%tangent_list(i),row(k)) * &
               bcs%tangents(1,i) * &
               bcs%tangents(1,bcs%lifted_ordering(row(k))) &
               + val(A12,bcs%tangent_list(i),row(k)) * &
               bcs%tangents(1,i) * &
               bcs%tangents(2,bcs%lifted_ordering(row(k))) &
               + val(A21,bcs%tangent_list(i),row(k)) * &
               bcs%tangents(2,i) * &
               bcs%tangents(1,bcs%lifted_ordering(row(k))) &
               + val(A22,bcs%tangent_list(i),row(k)) * &
               bcs%tangents(2,i) * &
               bcs%tangents(2,bcs%lifted_ordering(row(k))) 
               j = j+1
       end select
          
    end do

  end subroutine tangent_row

  subroutine tangent_row_block(i,layer_i,layer_j, &
       bcs,M,rpointer,rvalpointer)
    integer, intent(in) :: i,layer_i,layer_j
    type(bc_info), intent(in) :: bcs
    type(block_csr_matrix), intent(in) :: M
    integer, dimension(:), pointer :: rpointer
    real, dimension(:), pointer :: rvalpointer

    !locals
    integer :: j,ninternals,ntangents,k,rlength
    integer, dimension(:), pointer :: row
    real, dimension(:), pointer :: m11val,m12val,m21val,m22val

    !ewrite(1,*) '  subroutine tangent_row'

    if(associated(rpointer)) then
       deallocate(rpointer)
       rpointer=>null()
    end if
    if(associated(rvalpointer)) then
       deallocate(rvalpointer)
       rvalpointer=>null()
    end if
        
    ninternals=0
    ntangents=0

    row => M%colm(M%findrm(bcs%tangent_list(i)): &
         M%findrm(bcs%tangent_list(i)+1)-1)
    M11val => row_val_ptr(M,(layer_i-1)*2 + 1, &
         (layer_j-1)*2 + 1, bcs%tangent_list(i))
    M12val => row_val_ptr(M,(layer_i-1)*2 + 1, &
         (layer_j-1)*2 + 2, bcs%tangent_list(i))
    M21val => row_val_ptr(M,(layer_i-1)*2 + 2, &
         (layer_j-1)*2 + 1, bcs%tangent_list(i))
    M22val => row_val_ptr(M,(layer_i-1)*2 + 2, &
         (layer_j-1)*2 + 2, bcs%tangent_list(i))
    do k = 1, size(row)
       if(bcs%bc_marker(row(k))==1) ninternals = ninternals + 1
       if(bcs%bc_marker(row(k))==2) ntangents = ntangents + 1
    end do
    rlength = 2*ninternals + ntangents
    allocate(rpointer(rlength),rvalpointer(rlength))

    j = 1
    do k = 1,row_length(M,bcs%tangent_list(i))
       select case(bcs%bc_marker(row(k)))
       case (1)
          rpointer(j) = (layer_j-1)*(2*bcs%N_interior+bcs%N_tangents) + &
               bcs%lifted_ordering(row(k))
          rvalpointer(j) = m11val(k)*bcs%tangents(1,i) + &
               m21val(k)*bcs%tangents(2,i)
          rpointer(j+1) = (layer_j-1)*(2*bcs%N_interior+bcs%N_tangents) & 
               + bcs%N_interior + bcs%lifted_ordering(row(k))
          rvalpointer(j+1) = M12val(k)*bcs%tangents(1,i) + &
               m22val(k)*bcs%tangents(2,i)
          j = j+2
       case (2)
          rpointer(j) = (layer_j-1)*(2*bcs%N_interior+bcs%N_tangents)  & 
               + 2*bcs%N_interior + bcs%lifted_ordering(row(k))
          rvalpointer(j) = &
               m11val(k)*bcs%tangents(1,i) * &
               bcs%tangents(1,bcs%lifted_ordering(row(k))) &
               + M12val(k)*bcs%tangents(1,i) * &
               bcs%tangents(2,bcs%lifted_ordering(row(k))) &
               + M21val(k) * &
               bcs%tangents(2,i) * &
               bcs%tangents(1,bcs%lifted_ordering(row(k))) &
               + m22val(k)*bcs%tangents(2,i) * &
               bcs%tangents(2,bcs%lifted_ordering(row(k))) 
          j = j+1
       end select
          
    end do

  end subroutine tangent_row_block

end module gallopede_solvers
