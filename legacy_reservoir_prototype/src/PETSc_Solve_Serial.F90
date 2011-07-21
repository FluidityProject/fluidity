#include "fdebug.h"
module PETScSolve_Serial
  use fldebug


  ! PETSc Include files
#ifdef HAVE_PETSC

  use petsc
#if PETSC_VERSION_MINOR==0
  use petscvec 
  use petscmat 
  use petscksp 
  use petscpc 
  use petscis 
  use petscmg  
#endif
  implicit none
  
#include "petscversion.h"
#if PETSC_VERSION_MINOR==0
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscisdef.h"
#else
#include "finclude/petscdef.h"
#endif



#else
  implicit none


#endif
  ! PetSC
  logical,parameter :: usepetsc = .true.

  public :: PETScMatrixSolve


    ! THERE IS A SERIOUS BUG WITH CALLING THIS PETSC SOLVER WRAPPER
    ! IN THAT THE MEMORY USAGE CONTINUOUSLY GROWS. DO NOT USE IT.


contains

  subroutine PETScMatrixSolve(A,x,b,nonods, &
       & findfe,colfe,number_of_unknowns, max_resid)

    implicit none

    integer :: nonods,number_of_unknowns
    integer ::  findfe(nonods+1), colfe(number_of_unknowns)
    real    :: a(number_of_unknowns),x(nonods),b(nonods),max_resid
#ifdef HAVE_PETSC
    !Locals

  ! Petsc Pointers that need to be global
    Mat             :: petA
    PetscInt        :: peti,petn,petits,PETJ
    PetscScalar     :: petreal
    PetscErrorCode  :: ierr
    PetscMPIInt     :: petsize,rank


    Vec              :: petx,petb
    KSP              :: ksp
    PC               :: pc, subprec, subsubprec
    double precision :: petnorm,pettol, resid_norm
    PetscInt, allocatable :: petindex(:)
    PetscTruth       :: flg
    PetscReal :: rtol, abstol, dtol, residual_norm,rnorm_it
    PetscInt  :: maxits, mgcycles,mg_levels, itnum,max_nzeros,restart_level
    PetscInt  :: krylov_its, numrows, numcols, nzerosa, bs,one,zerop
    KSPConvergedReason                :: reason
    PetscLogDouble :: flops, flops1, flops2
    PetscInt, allocatable :: idxm(:), idxn(:)
    integer :: i, ele,first1,iloc,jloc,inod,jnod, base1, nloc,j, ifree
    integer :: nodi,nodj,element,eloc,acloc,jfree,matpos,ncolfe
    integer, allocatable :: mat_row_nzero(:)
    real :: time1, time2


    ! local variables
    character*80 :: boomeramg,pc_hypre_boomeramg_max_iter,pc_hypre_boomeramg_rtol, &
         vcycle,rel_tol_char,pc_prometheus_rowmax,relax_weight_all,outer_relax_weight_all, &
         pc_hypre_boomeramg_relax_weight_all,pc_hypre_boomeramg_outer_relax_weight_all

    integer :: iter_inner_max = 100,preconditioner = 5,vcycle_number = 1
    real    :: iter_inner_tolerance = 1.0e-6,gs_relax = 1.0
    real    :: vcycle_rel_tol = 1.0e-4,omega = 0.8

    allocate(mat_row_nzero(nonods),idxn(nonods))

    ! Find number of nonzeros for each row of the CSR matrix.

    max_nzeros = 0
    do j = 1, nonods
      ifree = findfe(j+1) - findfe(j)
       mat_row_nzero(j) = ifree
       if(max_nzeros < ifree) max_nzeros = ifree
    end do
    max_nzeros = nonods
    ewrite(1,*) '  Solve Using Petsc Routines ...'
    petsize  = 1; rank = 0 ; zerop = 0
    zerop = 0
    petn = nonods
    ewrite(1,*) '  On Rank:', petsize
    call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-petn', petn, flg, ierr)
    ewrite(1,*) '   constructing matrix:'

    call MatCreateSeqAIJ(PETSC_COMM_SELF,petn,petn,max_nzeros,PETSC_NULL_INTEGER,peta,ierr) 
    call MatCreateSeqAIJ(PETSC_COMM_SELF,petn,petn,zerop,mat_row_nzero,peta,ierr) 
    call MatSetFromOptions(peta,ierr) ! SET OPTIONS FROM COMMAND LINE

    ncolfe = number_of_unknowns ; one = 1 

    do nodi = 1, nonods
      do matpos = findfe(nodi), findfe(nodi+1) - 1
          nodj = colfe(matpos)
          peti = nodi - 1
          petj = nodj - 1
          petreal = a(matpos)
          !ewrite(1,*) petreal,peti,petj
          !call MatSetValues(PetA,one,peti,one,petj,petreal, INSERT_VALUES,ierr)
          call MatSetValue(peta,peti,petj,petreal,INSERT_VALUES,ierr)
       end do ! end of nodj count loop
    end do ! end of nodi loop
    
    
    call MatAssemblyBegin(petA,MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd(petA,MAT_FINAL_ASSEMBLY,ierr)
 !   call MatView(peta,PETSC_VIEWER_DRAW_WORLD,ierr);pause
    call VecCreate(PETSC_COMM_WORLD,petx,ierr)
    call VecSetSizes(petx,PETSC_DECIDE,petn,ierr)
    call VecSetFromOptions(petx,ierr)
    call VecDuplicate(petx,petb,ierr)
    ewrite(1,*) '  Set RHS Vector'
    do peti = 0, petn-1
       i = peti + 1
       petreal = b(i)
       call VecSetValue(petb,peti,petreal,INSERT_VALUES,ierr)
    end do

    boomeramg = 'boomeramg'    
    pc_hypre_boomeramg_max_iter = '-pc_hypre_boomeramg_max_iter'  
    pc_hypre_boomeramg_rtol = '-pc_hypre_boomeramg_rtol'  
    pc_hypre_boomeramg_relax_weight_all = '-pc_hypre_boomeramg_relax_weight_all'
    pc_hypre_boomeramg_outer_relax_weight_all = '-pc_hypre_boomeramg_outer_relax_weight_all'

    pc_prometheus_rowmax = '-pc_prometheus_rowmax'  


    ewrite(1,*)'  Set up Solver:'
    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
    call KSPSetOperators(ksp,petA,petA,DIFFERENT_NONZERO_PATTERN,ierr)
    ! Set Tolerances
    pettol = 1.d-10
    rtol = 1.0e-10    ! Relative Convergence Tolerance - relative decrease in the residual norm
    abstol = 1.0e-10 ! Absolute Convergence Tolerance - absolute size of the residual norm
    dtol = 1.0e2 ! divergence tolerance - amount residual can increase before KSPDefaultConverged() concludes that the method is diverging 
    maxits = 1000 ! Maximum Number of Iterations
    mg_levels = 4!Set multigrid levels
    restart_level = 10 !Number of Krylov Vectors + 1
    call KSPSetTolerances(ksp,PETtol,PETSC_DEFAULT_DOUBLE_PRECISION,  &
         &     PETSC_DEFAULT_DOUBLE_PRECISION,maxits,ierr)

    ! call KSPGMRESSetRestart(ksp,restart_level ,ierr)
    krylov_its = restart_level - 1
    ! call KSPGMRESModifiedGramSchmidtOrthogonalization(ksp,krylov_its,ierr)

    ! Set up Preconditioner
    call KSPGetPC(ksp, pc, ierr); !CHKERRQ(ierr)


    ! Setup the preconditioner
    ! set the preconditioner type 
    if (preconditioner .eq. 1) then
       ! jacobi
       call PCSetType(pc,PCJACOBI,ierr) 
    else if (preconditioner .eq. 2) then
       ! SOR   
       ! set SOR
       call PCSetType(pc,PCSOR,ierr) 
       ! set the relax parameter from options
       omega = gs_relax
       call PCSORSetOmega(pc,omega,ierr)   
    else if (preconditioner .eq. 3) then
       ! SSOR   
       ! set SOR
       call PCSetType(pc,PCSOR,ierr) 
       ! make it SSOR pc
       call PCSORSetSymmetric(pc,SOR_SYMMETRIC_SWEEP,ierr)
       ! set the relax parameter from options
       omega = gs_relax
       call PCSORSetOmega(pc,omega,ierr)
    else if (preconditioner .eq. 4) then
       ! Eigenstat (SSOR improved)   
       call PCSetType(pc,PCEISENSTAT,ierr) 
       ! set the relax parameter from options
       omega = gs_relax
       call PCEisenstatSetOmega(pc,omega,ierr)
    else if (preconditioner .eq. 5) then
       ! Incomplete LU   
       call PCSetType(pc,PCILU,ierr) 
    else if (preconditioner .eq. 6) then
       ! Incomplete Cholesky   
       call PCSetType(pc,PCICC,ierr) 
    else if (preconditioner .eq. 7) then
       ! HYPRE boomeramg 
       call PCSetType(pc,PCHYPRE,ierr)
       ! set to boomeramg
       call PCHYPRESetType(pc,boomeramg,ierr)
       write(unit=vcycle,fmt=*) vcycle_number
       ! set the multi grid cycles     
       !        call PetscOptionsSetValue(pc_hypre_boomeramg_max_iter,trim(vcycle),ierr)  
       write(unit=rel_tol_char,fmt=*) vcycle_rel_tol
       ! set the multigrid relative tolerance (should striclty be 0.0)
       call PetscOptionsSetValue(pc_hypre_boomeramg_rtol,trim(rel_tol_char),ierr)
       write(unit=relax_weight_all,fmt=*) gs_relax
       ! set the relax parameter  
       call PetscOptionsSetValue(pc_hypre_boomeramg_relax_weight_all,trim(relax_weight_all),ierr)         
       write(unit=outer_relax_weight_all,fmt=*) gs_relax
       ! set the outer relax parameter
       call PetscOptionsSetValue(pc_hypre_boomeramg_outer_relax_weight_all,trim(outer_relax_weight_all),ierr)
       !else if (preconditioner .eq. 8) then
       ! Prometheus (i.e. diagonal scaling preconditioning) 
       !call PCSetType(pc,PCPROMETHEUS,ierr) 
       ! use max of row rather than diagonal  
       !call PetscOptionsSetValue(pc_prometheus_rowmax,ierr)  
    else
       ! default use jacobi
       call PCSetType(pc,PCJACOBI,ierr) 
    end if ! if (preconditioner .eq. 1) then     


    ewrite(1,*) '  Finished building Vector'
    call KSPSetType(ksp, KSPFGMRES, ierr); !CHKERRQ(ierr)

    ! Set up KSP Options for PETSc
    call KSPSetFromOptions(ksp,ierr)
    call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)

    ! call KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, PETSC_NULL, &
    !           PETSC_NULL, ierr)

    ewrite(1,*) '  begin solver'
    call VecSet(petx,0.0, ierr)
    call cpu_time(time1)
    call PetscGetFlops(flops1, ierr)
    call KSPSolve(ksp,petb,PETx,ierr)
    call cpu_time(time2) 
    call PetscGetFlops(flops2, ierr)

    ewrite(1,*) '  End Solver' 
    ! call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
    allocate(petindex(nonods))
    do i =1, nonods
       petindex(i) = i-1
    end do
    ! Copy PETSc solution to fortran vector
    call VecGetValues(petx,petn, petindex ,x,ierr)
    ! Vector Norm of Soltution
    call VecNorm(petx,NORM_2,petnorm,ierr)
    ! Get number of iterations and reason for convergence
    call KSPGetIterationNumber(ksp,petits,ierr)
    call KSPGetConvergedReason(ksp, reason, ierr)
    ! Output solver stats if above a certain output.
    if(current_debug_level > 3) call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
    ewrite(3, "( '> PETSc diagnostics reason[', i5, '] iterations[', i5, ']' )") reason, petits
    call get_PETSC_Converged_Reason(reason)
    call KSPGetResidualNorm(ksp, residual_norm,ierr)
!    ewrite(3, "( '> PETSc Residual[', e10.6, ']')") residual_norm
    call PetscGetFlops(flops,ierr)
    ! write(*,*) "Number of MG LEVELS:", mg_levels
    ewrite(1,*) 'CPU time spent in solver:', time2-time1
    ewrite(1,*) 'MFlops counted by Petsc:', (flops2-flops1)/1e6
    ewrite(1,*) 'MFlops/sec:', (flops2-flops1)/((time2-time1)*1e6)
    ewrite(1,*) 'Residual Norm:', residual_norm


    resid_norm = residual_norm
    max_resid = resid_norm
    ! write(3,*) "RESIDUAL NORM:",resid_norm 
    ewrite(3,100) petnorm,petits
    call VecDestroy(petb,ierr)
    call VecDestroy(petx,ierr)
    call MatDestroy(petA,ierr)
    call KSPDestroy(ksp,ierr)

    !Write to file
    !open(unit = 3, file = "resid", status = 'UNKNOWN',FORM = 'formatted')
    !write(3,*) residual_norm, petits
    !write(3,100) petnorm,petits
    !write(3,*) 'CPU time spent in solver:', time2-time1
    !write(3,*) 'MFlops counted by Petsc:', (flops2-flops1)/1e6
    !write(3,*) 'MFlops/sec:', (flops2-flops1)/((time2-time1)*1e6)
    ! close(3) 
    DEALLOCATE(MAT_ROW_NZERO)

100 format('  Norm of solution = ',e10.4,',  Iterations = ',i5)
200 format('Norm of error < 1.e-12,Iterations = ',i5)
300 format('Norm of residual = ',e10.6)

    ! call PetscFinalize(ierr)

#endif
  end subroutine PETScMatrixSolve


  subroutine get_PETSC_Converged_Reason(reason)
    implicit none
    integer :: reason
    select case(reason)
       !! CONVERGED CASES
    case(2)
       ewrite(1,*) "  PETSc Converged Relative Tolerance"
    case(3)
       ewrite(1,*) "  PETSc Converged Absolute Tolerance"
    case(4)
       ewrite(1,*) "  PETSc Converged Iterations"
    case(5)
       ewrite(1,*) "  PETSc Converged CG Neg Curve"
    case(6)
       ewrite(1,*) "  PETSc Converged CG Constrained"
    case(7)
       ewrite(1,*) "  PETSc Converged Step Length"
    case(8)
       ewrite(1,*) "  PETSc Converged Happy Breakdown"
       !! DIVERGED CASES
    case(-2)
       ewrite(-1,*) "  PETSc Diverged Null"
    case(-3)
       ewrite(-1,*) "  PETSc Diverged Iterations"
    case(-4)
       ewrite(-1,*) "  PETSc Diverged Diverged Tolerance Reached"
    case(-5)
       ewrite(-1,*) "  PETSc Diverged Breakdown"
    case(-6)
       ewrite(-1,*) "  PETSc Diverged Breakdown BiCG"
    case(-7)
       ewrite(-1,*) "  PETSc Diverged Non Symmetric"
    case(-8)
       ewrite(-1,*) "  PETSc Diverged Indefinite Preconditioner"
    case(-9)
       ewrite(-1,*) "  PETSc Diverged NaN"
    case(-10)
       ewrite(-1,*) "  PETSc Diverged Indefinite Matrix"


    end select

  end subroutine get_PETSC_Converged_Reason

  subroutine petsc_finish()
    implicit none
#ifdef HAVE_PETSC
    PetscErrorCode   :: ierr
    call PetscFinalize(ierr)
#endif 
  end subroutine petsc_finish

  subroutine petsc_start()
    implicit none
#ifdef HAVE_PETSC
    PetscErrorCode   :: ierr
    call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
#endif 
  end subroutine petsc_start
  

end Module PETScSolve_Serial
