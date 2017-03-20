#include "fdebug.h"
!! This module contains multigrid related subroutines, such as the smoothed
!! aggregation preconditioner.
module multigrid
use FLDebug
use spud
use futils
use parallel_tools
#ifdef HAVE_PETSC_MODULES
  use petsc
#endif
use Sparse_tools
use Petsc_Tools
use sparse_tools_petsc
implicit none
#include "petsc_legacy.h"

!! Some parameters that change the behaviour of 
!! the smoothed aggregation method. All of
!! of these can also be set as PETSC_options:
!! 
!! --mymg_maxlevels, --mymg_coarsesize, --mymg_epsilon and --mymg_omega
!!
!! maximum number of multigrid levels:
integer, public, parameter:: MULTIGRID_MAXLEVELS_DEFAULT=25
!! the maximum number of nodes at the coarsest level 
!! (that is solved by a direct solver):
!! in serial we use a direct solver:
integer, public, parameter:: MULTIGRID_COARSESIZE_DEFAULT_SERIAL=5000
!! in parallel we coarsen a bit further:
integer, public, parameter:: MULTIGRID_COARSESIZE_DEFAULT_PARALLEL=100
!! epsilon determines the relatively strong connections:
PetscReal, public, parameter:: MULTIGRID_EPSILON_DEFAULT=0.01
!! epsilon is divided by epsilon_decay after each coarsening,
!! as the coarser levels are generally less anisotropic and therefore
!! we need a less strong criterium for strong connections
PetscReal, public, parameter:: MULTIGRID_EPSILON_DECAY_DEFAULT=1.0
!! omega in the prolongation smoother:
PetscReal, public, parameter:: MULTIGRID_OMEGA_DEFAULT=2.0/3.0
!! number of smoother iterations going down
integer, public, parameter:: MULTIGRID_NOSMD_DEFAULT=1
!! number of smoother iterations going up
integer, public, parameter:: MULTIGRID_NOSMU_DEFAULT=1
!! max size of clusters (in first round), 0 means follow Vanek'96
integer, public, parameter:: MULTIGRID_CLUSTERSIZE_DEFAULT=0

integer, private, parameter:: ISOLATED=0, COUPLED=-1

integer, public, parameter :: &
     !No internal smoothing
     INTERNAL_SMOOTHING_NONE=0, &
     !No SOR on outer level of the multigrids, whole PC wrapped in 
     !SOR sweeps
     INTERNAL_SMOOTHING_WRAP_SOR=1, &
     !SOR on outer level of the multigrids, no wrapping SOR
     INTERNAL_SMOOTHING_SEPARATE_SOR=2

!=====================================
!!Stuff for internal smoother -- cjc
!=====================================
!PC for internal smoother, it will be MG
PC :: internal_smoother_pc
!Matrix for internal smoother
Mat :: internal_smoother_mat
!list of surface nodes
integer, dimension(:), pointer, save :: surface_node_list => null()
!list of surface values, for copying
PetscReal, dimension(:), pointer, save :: surface_values => null()
!=====================================

private
public SetupSmoothedAggregation, SetupMultigrid, DestroyMultigrid
  
  ! see at the bottom of Petsc_Tools.F90
  ! For some reason "use"ing the interface from the petsc_tools module
  ! - if made public - doesn't work.
  interface
    subroutine myMatGetInfo(A, flag, info, ierr)
       Mat, intent(in):: A
       MatInfoType, intent(in):: flag
       double precision, dimension(:), intent(out):: info
       PetscErrorCode, intent(out):: ierr
    end subroutine myMatGetInfo
  end interface

contains

subroutine SetUpInternalSmoother(surface_node_list_in,matrix,pc, &
     no_top_smoothing)
!!< This subroutine sets up the internal additive smoother as described in
!!< Kramer and Cotter (2009) in preparation
  integer, intent(in), dimension(:) :: surface_node_list_in
  type(csr_matrix), intent(in) :: matrix
  PC, intent(inout) :: pc
  logical, intent(in), optional :: no_top_smoothing
  !
  PetscObject:: myPETSC_NULL_OBJECT
  type(csr_matrix) :: matrix_internal
  integer :: row,i, ierr, nsurface
  integer, dimension(:), pointer :: r_ptr
  
  integer, dimension(:), allocatable :: surface_flag
  logical :: lno_top_smoothing
  ! these are used to point at things in csr_matrices so should be using
  ! real instead of PetscReal
  real, dimension(:), pointer :: r_val_ptr

  lno_top_smoothing = .false.
  if(present(no_top_smoothing)) then
     lno_top_smoothing = no_top_smoothing
  end if
  nsurface = size(surface_node_list_in)
  allocate(surface_node_list(nsurface))
  surface_node_list = surface_node_list_in
  allocate(surface_values(nsurface))
  surface_values = 0.

  allocate(surface_flag(size(matrix,1)))
  surface_flag = 1
  surface_flag(surface_node_list) = 0

  call allocate(matrix_internal,matrix%sparsity)
  call set(matrix_internal,matrix)
  
  !zero all surface columns in matrix_internal
  do row = 1, size(matrix,1)
     r_ptr => row_m_ptr(matrix_internal,row)
     r_val_ptr => row_val_ptr(matrix_internal,row)
     
     if(any(surface_flag(r_ptr)==0)) then
        call set(matrix_internal,row,r_ptr,r_val_ptr*surface_flag(r_ptr))
     end if
  end do
  
  !zero all surface row in matrix_internal
  !and put a 1 on the diagonal
  do i = 1, size(surface_node_list)
     row = surface_node_list(i)
     r_ptr => row_m_ptr(matrix_internal,row)
     r_val_ptr => row_val_ptr(matrix_internal,row)
     
     call set(matrix_internal,row,r_ptr,0.0*r_val_ptr)
     call set(matrix_internal,row,row,1.0)
  end do

  !create matrix
  internal_smoother_mat = csr2petsc(matrix_internal)
  !set up multigrid
  call PCCreate(MPI_COMM_SELF,internal_smoother_pc,ierr)
  call PCSetType(internal_smoother_pc,PCMG,ierr)
  call SetupSmoothedAggregation(internal_smoother_pc, &
       Internal_Smoother_Mat, ierr, no_top_smoothing=lno_top_smoothing)

  ! PCSetOperators needs to be in small caps due to macro hack in include/petsc_legacy.h
  call pcsetoperators(internal_smoother_pc,Internal_Smoother_Mat, Internal_Smoother_Mat, ierr)

  !set up pc to output
  myPETSC_NULL_OBJECT=PETSC_NULL_OBJECT
  call PCShellSetApply(pc,ApplySmoother,PETSC_NULL_OBJECT,ierr)
  if (myPETSC_NULL_OBJECT/=PETSC_NULL_OBJECT) then
    FLAbort("PETSC_NULL_OBJECT has changed please report to skramer")
  end if  

  surface_node_list = surface_node_list - 1

end subroutine SetUpInternalSmoother

subroutine ApplySmoother(dummy,vec_in,vec_out,ierr_out)
  Vec, intent(in) :: vec_in
  Vec, intent(inout) :: vec_out
  PetscErrorCode, intent(inout) :: ierr_out
  integer, intent(in) :: dummy
  
  integer :: ierr
  
  ierr_out=0
  
  call PCApply(internal_smoother_pc,vec_in,vec_out, ierr)

  call VecSetValues(vec_out, size(surface_node_list), &
       surface_node_list, real(spread(0.0, 1, size(surface_values)), kind = PetscScalar_kind), &
       INSERT_VALUES, ierr)

end subroutine ApplySmoother

subroutine DestroyInternalSmoother()
!!< Remove all trace of this preconditioner
  implicit none
  !
  integer :: ierr
  !
  deallocate( surface_node_list )
  surface_node_list => null()
  deallocate( surface_values )
  surface_values => null()
  call PCDestroy(internal_smoother_pc,ierr)
  call MatDestroy(internal_smoother_mat,ierr)

end subroutine DestroyInternalSmoother

subroutine SetupMultigrid(prec, matrix, ierror, &
  external_prolongators, surface_node_list, matrix_csr, &
  internal_smoothing_option)
!!< This subroutine sets up the multigrid preconditioner including
!!< all options (vertical_lumping, internal_smoother)
PC, intent(inout):: prec
Mat, intent(in):: matrix
!! ierror=0 upon succesful return, otherwise ierror=1 and everything 
!! will be deallocated
integer, intent(out):: ierror
!! use external prolongator at the finest level
type(petsc_csr_matrix), dimension(:), optional, intent(in):: external_prolongators
!! if present, use additve smoother that solves the eliptic problem with 
!! the solution of the last multigrid iteration at the top surface as
!! dirichlet boundary condition
integer, optional, dimension(:):: surface_node_list
type(csr_matrix), intent(in), optional :: matrix_csr
integer, optional, intent(in) :: internal_smoothing_option

integer :: linternal_smoothing_option

  PetscErrorCode ierr
  PC subprec, subsubprec
  
  !! Get internal smoothing options
  linternal_smoothing_option = INTERNAL_SMOOTHING_NONE
  if(present(internal_smoothing_option)) then
     linternal_smoothing_option = internal_smoothing_option
  end if

  select case (linternal_smoothing_option)
  case (INTERNAL_SMOOTHING_NONE)
     !Don't apply internal smoothing, just regular mg
     call SetupSmoothedAggregation(prec, matrix, ierror, &
          external_prolongators=external_prolongators)
  case (INTERNAL_SMOOTHING_WRAP_SOR)
     !Apply the internal smoothing with wrapped SOR
     if(.not.present(surface_node_list)) then
        FLAbort('surface_node_list is needed for chosen internal smoothing option')
     end if
     if(.not.present(matrix_csr)) then
        FLAbort('matrix_csr must also be provided for internal smoother')
     end if

     ! The overall preconditioner is compositive multiplicative
     call PCSetType(prec, PCCOMPOSITE, ierr)
     call PCCompositeSetType(prec, PC_COMPOSITE_MULTIPLICATIVE, ierr)
     ! consisting of outer SOR iterations and the composite PC
     call PCCompositeAddPC(prec, PCSOR, ierr)     
     call PCCompositeAddPC(prec, PCCOMPOSITE, ierr)
     call PCCompositeAddPC(prec, PCSOR, ierr)
     
     !set up the forward SOR
     call PCCompositeGetPC(prec, 0, subprec, ierr)
     call PCSORSetSymmetric(subprec,SOR_FORWARD_SWEEP,ierr)
     call PCSORSetIterations(subprec,1,1,ierr)
     call PCSORSetOmega(subprec,real(1.0, kind = PetscReal_kind),ierr)

     ! set up the backward SOR
     call PCCompositeGetPC(prec, 2, subprec, ierr)
     call PCSORSetSymmetric(subprec,SOR_BACKWARD_SWEEP,ierr)
     call PCSORSetIterations(subprec,1,1,ierr)
     call PCSORSetOmega(subprec,real(1.0, kind = PetscReal_kind),ierr)
     
     !set up the middle PC 
     call PCCompositeGetPC(prec, 1, subprec, ierr)
     !It is compositive additive
     call PCSetType(subprec, PCCOMPOSITE, ierr)
     call PCCompositeSetType(subprec, PC_COMPOSITE_ADDITIVE, ierr)
     !consisting of the vertical lumped mg, and the internal smoother 
     !which is a shell
     call PCCompositeAddPC(subprec, PCMG, ierr)
     call PCCompositeAddPC(subprec, PCSHELL, ierr)
     ! set up the vertical_lumped mg
     call PCCompositeGetPC(subprec, 0, subsubprec, ierr)
     call SetupSmoothedAggregation(subsubprec, matrix, ierror, &
          external_prolongators, no_top_smoothing=.true.)
     !set up the "internal" mg shell
     call PCCompositeGetPC(subprec, 1, subsubprec, ierr)
     call SetupInternalSmoother(surface_node_list,matrix_csr,subsubprec, &
          no_top_smoothing=.true.)

  case (INTERNAL_SMOOTHING_SEPARATE_SOR)
     !Apply the internal smoothing with separate SOR for each mg
     if(.not.present(surface_node_list)) then
        FLAbort('surface_node_list is needed for chosen internal smoothing option')
     end if
     if(.not.present(matrix_csr)) then
        FLAbort('matrix_csr must also be provided for internal smoother')
     end if

     !It is compositive additive
     call PCSetType(prec, PCCOMPOSITE, ierr)
     call PCCompositeSetType(prec, PC_COMPOSITE_ADDITIVE, ierr)
     !consisting of the vertical lumped mg, and the internal smoother 
     !which is a shell
     call PCCompositeAddPC(prec, PCMG, ierr)
     call PCCompositeAddPC(prec, PCSHELL, ierr)
     ! set up the vertical_lumped mg
     call PCCompositeGetPC(prec, 0, subprec, ierr)
     call SetupSmoothedAggregation(subprec, matrix, ierror, &
          external_prolongators=external_prolongators)
     ! set up the "internal" mg shell
     call PCCompositeGetPC(prec, 1, subprec, ierr)
     call SetupInternalSmoother(surface_node_list,matrix_csr,subprec)

  case default
     FLAbort('bad internal smoothing option')
  end select

end subroutine SetupMultigrid
  
subroutine DestroyMultigrid(prec)
PC, intent(inout):: prec

  PetscErrorCode ierr
  PCType pctype
  PC subprec
  
  call PCGetType(prec, pctype, ierr)
  if (pctype==PCCOMPOSITE) then
    ! destroy the real mg
    call PCCompositeGetPC(prec, 0, subprec, ierr)
    ! destroy the internal mg
    call PCCompositeGetPC(prec, 1, subprec, ierr)
    call DestroyInternalSmoother()
  end if
  
end subroutine DestroyMultigrid

subroutine SetupSmoothedAggregation(prec, matrix, ierror, &
  external_prolongators,no_top_smoothing)
!!< This subroutine sets up the preconditioner for using the smoothed
!!< aggregation method (as described in Vanek et al. 
!!< Computing 56, 179-196 (1996).
PC, intent(inout):: prec
Mat, intent(in):: matrix
!! ierror=0 upon succesful return, otherwise ierror=1 and everything 
!! will be deallocated
integer, intent(out):: ierror
!! use external prolongator at the finest level
type(petsc_csr_matrix), dimension(:), optional, intent(in):: external_prolongators
!! Don't do smoothing on the top level
logical, intent(in), optional :: no_top_smoothing

  Mat, allocatable, dimension(:):: matrices, prolongators
  KSP ksp_smoother
  PC  prec_smoother
  Vec lvec, rvec
  MatNullSpace nullsp
  PetscErrorCode ierr
  PetscReal epsilon, epsilon_decay, omega
  PetscInt maxlevels, coarsesize
  PetscReal:: eigval
  PetscScalar :: Px2
  Vec:: eigvec, Px
  PetscReal, allocatable, dimension(:):: emin, emax
  PetscObject:: myPETSC_NULL_OBJECT
  integer, allocatable, dimension(:):: contexts
  integer i, j, ri, nolevels, m, n, top_level
  integer nosmd, nosmu, clustersize, no_external_prolongators
  logical forgetlastone
  logical lno_top_smoothing

    ! this might be already done, but it doesn't hurt:
    call PCSetType(prec, PCMG, ierr)

    call PCMGGetLevels(prec, nolevels, ierr)
    if (ierr==0 .and. nolevels>0) then
      ewrite(2,*) "Assuming mg preconditioner is used with the same matrix and same options as before"
      ierror = 0
      return
    end if

    lno_top_smoothing = .false.
    if(present(no_top_smoothing)) then
       lno_top_smoothing = no_top_smoothing
    end if

    call SetSmoothedAggregationOptions(epsilon, epsilon_decay, omega, maxlevels, coarsesize, &
      nosmd, nosmu, clustersize)
      
    ! In the following level i=1 is the original, fine, problem
    ! i=nolevels corresponds to the coarsest problem
    !
    ! Unfortunately PETSc uses a reverse numbering in which ri=0
    ! is the coarsest problem and ri=nolevels-1 the finest problem
    !
    ! maxlevels is the, user specified (default=25) maximum for nolevels
    !   matrices(i)     is PETSc matrix at level i
    !   prolongators(i) is PETSc prolongator between level i+1 and level i
    !       its transpose is the restriction between level i and level i+1
    allocate(matrices(1:maxlevels), prolongators(1:maxlevels-1), &
      contexts(1:maxlevels-1))
      
    if (present(external_prolongators)) then
      no_external_prolongators=size(external_prolongators)
    else
      no_external_prolongators=0
    end if

    forgetlastone=.false.
    matrices(1)=matrix
    do i=1, maxlevels-1
      ewrite(3,*) '---------------------'
      ewrite(3,*) 'coarsening from level',i,' to ',i+1
      if (i<=no_external_prolongators) then
         prolongators(i)=external_prolongators(i)%M
         ewrite(2,*) "Using provided external prolongator"
         ewrite(2,*) "Coarsening from", size(external_prolongators(i),1), &
           "to", size(external_prolongators(i),2), "nodes"
      else
         prolongators(i)=Prolongator(matrices(i), epsilon, omega, clustersize)
         epsilon=epsilon/epsilon_decay
      end if

      if (prolongators(i)==PETSC_NULL_OBJECT) then
        if (IsParallel()) then
          ! in parallel we give up
          ewrite(-1,*) "ERROR: mg preconditioner setup failed"
          ewrite(-1,*) "This may be caused by local partitions being too small"
        else
          ! in serial you may want to try something else automatically
          ewrite(0,*) 'WARNING: mg preconditioner setup failed'          
          ewrite(0,*) 'This probably means the matrix is not suitable for it.'
        end if
        do j=1+no_external_prolongators, i-1
          call MatDestroy(prolongators(j), ierr)
          call MatDestroy(matrices(j), ierr)
        end do
        deallocate(matrices, prolongators, contexts)
        ! Need to set n/o levels (to 1) otherwise PCDestroy will fail:
        myPETSC_NULL_OBJECT=PETSC_NULL_OBJECT
        call PCMGSetLevels(prec, 1, PETSC_NULL_OBJECT, ierr)
        if (myPETSC_NULL_OBJECT/=PETSC_NULL_OBJECT) then
           FLAbort("PETSC_NULL_OBJECT has changed please report to skramer")
        end if
        ierror=1
        return
      end if
      
      ! prolongator between i+1 and i, is restriction between i and i+1
      call MatGetLocalSize(prolongators(i), m, n, ierr)
      call MatPtAP(matrices(i), prolongators(i), MAT_INITIAL_MATRIX, real(1.0, kind = PetscReal_kind), &
        matrices(i+1), ierr)
      
      call allmin(n)
      if (n<coarsesize) exit
    end do
    
    if (forgetlastone) i=i-1
    
    if (i<maxlevels) then
      nolevels=i+1
    else
      nolevels=i
    end if
    
    myPETSC_NULL_OBJECT=PETSC_NULL_OBJECT
    call PCMGSetLevels(prec, nolevels, PETSC_NULL_OBJECT, ierr)
    if (myPETSC_NULL_OBJECT/=PETSC_NULL_OBJECT) then
       FLAbort("PETSC_NULL_OBJECT has changed please report to skramer")
    end if
    
    if (lno_top_smoothing) then
      top_level=nolevels-2
      ! set smoother at finest level to "none"
      call PCMGGetSmootherUp(prec, nolevels-1, ksp_smoother, ierr)
      call SetupNoneSmoother(ksp_smoother, matrices(1))
    else
      top_level=nolevels-1
    end if
    
    if (nosmu<0 .and. nosmd<0) then
      
      allocate(emin(1:nolevels-1), emax(1:nolevels))
      
      call PowerMethod(matrices(1), eigval, eigvec)
      emax(1)=eigval
      call VecDestroy(eigvec, ierr)
      
      ! loop over reverse index where 1 is fine and nolevels coarse:
      do ri=2, nolevels
        call PowerMethod(matrices(ri), eigval, eigvec)
        emax(ri)=eigval
        
        myPETSC_NULL_OBJECT=PETSC_NULL_OBJECT
        call MatCreateVecs(prolongators(ri-1), PETSC_NULL_OBJECT, Px, ierr)
        if (myPETSC_NULL_OBJECT/=PETSC_NULL_OBJECT) then
           FLAbort("PETSC_NULL_OBJECT has changed please report to skramer")
        end if
        call MatMult(prolongators(ri-1), eigvec, Px, ierr)
        call VecNorm(Px, NORM_2, Px2, ierr)
        emin(ri-1)=eigval/Px2**2.
        call VecDestroy(Px, ierr)
        
        call VecDestroy(eigvec, ierr)
      end do
      
          
      ! loop over the 'PETSc' index where 0 is coarse and nolevels-1 is fine:
      do i=1, top_level
        ! reverse index where 1 is fine and nolevels coarse:
        ri=nolevels-i
        ewrite(2,*) "level, emin, emax", ri, emin(ri), emax(ri)
        call PCMGGetSmootherUp(prec, i, ksp_smoother, ierr)
        call SetupChebychevSmoother(ksp_smoother, matrices(ri), &
          emin(ri), emax(ri), -nosmu)
        call PCMGGetSmootherDown(prec, i, ksp_smoother, ierr)
        call SetupChebychevSmoother(ksp_smoother, matrices(ri), &
          emin(ri), emax(ri), -nosmd)
      end do
        
      deallocate(emin, emax)
      
    else if (nosmu>0 .and. nosmd>0) then
    
      ! loop over the 'PETSc' index where 0 is coarse and nolevels-1 is fine:
      do i=1, top_level
        ! reverse index where 1 is fine and nolevels coarse:
        ri=nolevels-i
        call PCMGGetSmootherUp(prec, i, ksp_smoother, ierr)
        if (IsParallel()) then
          call SetupSORSmoother(ksp_smoother, matrices(ri), SOR_LOCAL_SYMMETRIC_SWEEP, nosmu)
        else
          call SetupSORSmoother(ksp_smoother, matrices(ri), SOR_FORWARD_SWEEP, nosmu)
        end if

        call PCMGGetSmootherDown(prec, i, ksp_smoother, ierr)
        if (lno_top_smoothing .and. i==nolevels-1) then
          call SetupNoneSmoother(ksp_smoother, matrices(ri))
        else if (IsParallel()) then
          call SetupSORSmoother(ksp_smoother, matrices(ri), SOR_LOCAL_SYMMETRIC_SWEEP, nosmu)
        else
          call SetupSORSmoother(ksp_smoother, matrices(ri), SOR_BACKWARD_SWEEP, nosmd)
        end if
        
      end do
      
    else
    
      FLAbort("Can't combine chebychev and sor smoothing")
      
    end if
      
    ! loop over the 'PETSc' index where 0 is coarse and nolevels-1 is fine:
    do i=1, nolevels-1
      ! reverse index where 1 is fine and nolevels coarse:
      ri=nolevels-i
      
      call PCMGSetInterpolation(prec, i, prolongators(ri), ierr)
      call PCMGSetRestriction(prec, i, prolongators(ri), ierr)
      ! This won't actually destroy them, as they are still refered to
      ! as interpolators and restrictors. PETSc will destroy them once
      ! these references are destroyed in KSP/PCDestroy:
      if (ri>no_external_prolongators) then
        call MatDestroy(prolongators(ri), ierr)
      end if
      
    end do
    
    ! Create rhs's for coarsest to one but finest level:
    ! This shouldn't be necessary, but PETSc messes up leaving
    !   the vector when destroying the preconditioner. 
    ! (believed fixed in PETSc 3.0.0)
    do i=0, nolevels-2
      ri=nolevels-i
      ! using PETSC_NULL_OBJECT for rvec leaks a reference
      call MatCreateVecs(matrices(ri), lvec, rvec, ierr)
      call PCMGSetRHS(prec, i, lvec, ierr)
      ! Again, this does not yet destroy rhs immediately:
      call VecDestroy(lvec, ierr)
      call VecDestroy(rvec, ierr)
    end do
      
    ! residual needs to be set if PCMG is used with KSPRICHARDSON
    call MatCreateVecs(matrices(1), lvec, rvec, ierr)
    call PCMGSetR(prec, nolevels-1, lvec, ierr)
    call VecDestroy(lvec, ierr)
    call VecDestroy(rvec, ierr)

    ! solver options coarsest level:
    call PCMGGetCoarseSolve(prec, ksp_smoother, ierr)
    call KSPGetPC(ksp_smoother, prec_smoother, ierr)
    call MatGetNullSpace(matrix, nullsp, ierr)
    if (IsParallel() .or. (nullsp/=PETSC_NULL_OBJECT .and. ierr==0)) then
      ! if parallel or if we have a null space: use smoothing instead of direct solve
      call SetupSORSmoother(ksp_smoother, matrices(nolevels), &
        SOR_LOCAL_SYMMETRIC_SWEEP, 20)
    else
      call KSPSetOperators(ksp_smoother, matrices(nolevels), matrices(nolevels), ierr)
      call KSPSetType(ksp_smoother, KSPPREONLY, ierr)
      call PCSetType(prec_smoother, PCLU, ierr)
      call KSPSetTolerances(ksp_smoother, 1.0e-100_PetscReal_kind, 1e-8_PetscReal_kind, 1e10_PetscReal_kind, 300, ierr)
    end if
    
    ! destroy our references of the operators at levels 2 (==one but finest) up to coarsest
    do i=2, nolevels
      call MatDestroy(matrices(i), ierr)
    end do
    
    deallocate(matrices, prolongators, contexts)

    ! succesful return
    ierror=0
    
end subroutine SetupSmoothedAggregation
  
subroutine SetupSORSmoother(ksp, matrix, sortype, iterations)
KSP, intent(in):: ksp
Mat, intent(in):: matrix
MatSORType, intent(in):: sortype
integer, intent(in):: iterations
  
  PC:: pc
  PetscErrorCode:: ierr
  
  call KSPSetType(ksp, KSPRICHARDSON, ierr)
  call KSPSetOperators(ksp, matrix, matrix, ierr)
  ! set 1 richardson iteration, as global iteration inside pcsor might be more efficient
  call KSPSetTolerances(ksp, PETSC_DEFAULT_REAL, &
    PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
    1, ierr)
  call KSPSetNormType(ksp, KSP_NORM_NONE, ierr)
  
  call KSPGetPC(ksp, pc, ierr)
  call PCSetType(pc, PCSOR, ierr)
  call PCSORSetSymmetric(pc, sortype, ierr)
  call PCSORSetOmega(pc, real(1.0, kind = PetscReal_kind), ierr)
  call PCSORSetIterations(pc, iterations, 1, ierr)

end subroutine SetupSORSmoother

subroutine SetupNoneSmoother(ksp, matrix)
KSP, intent(in):: ksp
Mat, intent(in):: matrix
  
  PC:: pc
  PetscErrorCode:: ierr
  
  call KSPSetType(ksp, KSPRICHARDSON, ierr)
  call KSPSetOperators(ksp, matrix, matrix, ierr)
  call KSPSetTolerances(ksp, PETSC_DEFAULT_REAL, &
    PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
    0, ierr)
  call KSPRichardsonSetScale(ksp,real(0.0, kind = PetscReal_kind),ierr)
  call KSPSetNormType(ksp, KSP_NORM_NONE, ierr)
  
  call KSPGetPC(ksp, pc, ierr)
  call PCSetType(pc,PCNONE,ierr)
  
end subroutine SetupNoneSmoother
    
subroutine SetupChebychevSmoother(ksp, matrix, emin, emax, iterations)
KSP, intent(in):: ksp
Mat, intent(in):: matrix
PetscReal, intent(in):: emin, emax
integer, intent(in):: iterations
  
  PC:: pc
  PetscErrorCode:: ierr
  
  call KSPSetType(ksp, KSPCHEBYSHEV, ierr)
  call KSPSetOperators(ksp, matrix, matrix, ierr)
  call KSPSetTolerances(ksp, PETSC_DEFAULT_REAL, &
    PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
    iterations, ierr)
  call KSPChebyshevSetEigenvalues(ksp, emax, emin, ierr)
  call KSPSetNormType(ksp, KSP_NORM_NONE, ierr)

  call KSPGetPC(ksp, pc, ierr)
  call PCSetType(pc, PCNONE, ierr)

end subroutine SetupChebychevSmoother
  
subroutine SetSmoothedAggregationOptions(epsilon, epsilon_decay, omega, maxlevels, &
  coarsesize, nosmd, nosmu, clustersize)
PetscReal, intent(out):: epsilon, epsilon_decay, omega
integer, intent(out):: maxlevels, coarsesize
integer, intent(out):: nosmd, nosmu, clustersize

  PetscBool flag
  PetscErrorCode ierr

    call PetscOptionsGetReal(PETSC_NULL_OBJECT, '', '-mymg_epsilon', epsilon, flag, ierr)
    if (.not. flag) then
      epsilon=MULTIGRID_EPSILON_DEFAULT
    end if
    call PetscOptionsGetReal(PETSC_NULL_OBJECT, '', '-mymg_epsilon_decay', epsilon_decay, flag, ierr)
    if (.not. flag) then
      epsilon_decay=MULTIGRID_EPSILON_DECAY_DEFAULT
    end if
    call PetscOptionsGetReal(PETSC_NULL_OBJECT, '', '-mymg_omega', omega, flag, ierr)
    if (.not. flag) then
      omega=MULTIGRID_OMEGA_DEFAULT
    end if
    call PetscOptionsGetInt(PETSC_NULL_OBJECT, '', '-mymg_maxlevels', maxlevels, flag, ierr)
    if (.not. flag) then
      maxlevels=MULTIGRID_MAXLEVELS_DEFAULT
    end if
    call PetscOptionsGetInt(PETSC_NULL_OBJECT, '', '-mymg_coarsesize', coarsesize, flag, ierr)
    if (.not. flag) then
      if (IsParallel()) then
        coarsesize=MULTIGRID_COARSESIZE_DEFAULT_PARALLEL
      else
        coarsesize=MULTIGRID_COARSESIZE_DEFAULT_SERIAL
      end if
    end if
    call PetscOptionsGetInt(PETSC_NULL_OBJECT, '', '-mymg_nosmd', nosmd, flag, ierr)
    if (.not. flag) then
      nosmd=MULTIGRID_NOSMD_DEFAULT
    end if
    call PetscOptionsGetInt(PETSC_NULL_OBJECT, '', '-mymg_nosmu', nosmu, flag, ierr)
    if (.not. flag) then
      nosmu=MULTIGRID_NOSMU_DEFAULT
    end if
    call PetscOptionsGetInt(PETSC_NULL_OBJECT, '', '-mymg_clustersize', clustersize, flag, ierr)
    if (.not. flag) then
      clustersize=MULTIGRID_CLUSTERSIZE_DEFAULT
    end if

    ewrite(2,*) 'multgrid options- epsilon:', epsilon
    ewrite(2,*) 'multgrid options- epsilon_decay:', epsilon_decay
    ewrite(2,*) 'multgrid options- omega: ', omega
    ewrite(2,*) 'multgrid options- maxlevels:', maxlevels
    ewrite(2,*) 'multgrid options- coarsesize: ', coarsesize
    ewrite(2,*) 'multgrid options- n/o smoother its down (nosmd): ', nosmd
    ewrite(2,*) 'multgrid options- n/o smoother its up (nosmu): ', nosmu
    ewrite(2,*) 'multgrid options- maximum clustersize: ', clustersize
    
end subroutine SetSmoothedAggregationOptions

function Prolongator(A, epsilon, omega, maxclustersize, cluster) result (P)
!!< Constructs coarse grid and prolongator operator between coarse and fine
!!< grid based on the matrix A.
Mat:: P
Mat, intent(in):: A
!! strong connectivity criterion as in Vanek '96
PetscReal, intent(in):: epsilon 
!! overrelaxtion in jacobi-smoothed aggregation
PetscReal, intent(in):: omega
!! maximum size of clusters
integer, intent(in):: maxclustersize 
!! 0 means only clustering as in Vanek '96
!! if supplied returns cluster number each node is assigned to:
integer, optional, dimension(:), intent(out):: cluster

  PetscErrorCode:: ierr
  PetscInt:: diagminloc
  PetscReal:: diagmin
  Vec:: sqrt_diag, inv_sqrt_diag, diag, one
  double precision, dimension(MAT_INFO_SIZE):: matrixinfo
  integer, dimension(:), allocatable:: findN, N, R
  integer:: nrows, nentries, ncols
  integer:: jc, ccnt, base
    
  ! find out basic dimensions of A
  call MatGetLocalSize(A, nrows, ncols, ierr)
  ! use Petsc_Tools's MatGetInfo because of bug in earlier patch levels of petsc 3.0
  call myMatGetInfo(A, MAT_LOCAL, matrixinfo, ierr)
  nentries=matrixinfo(MAT_INFO_NZ_USED)
  call MatGetOwnerShipRange(A, base, PETSC_NULL_INTEGER, ierr)
  ! we decrease by 1, so base+i gives 0-based petsc index if i is the local fortran index:
  base=base-1
  
  allocate(findN(1:nrows+1), N(1:nentries), R(1:nrows))
     
  ! rescale the matrix: a_ij -> a_ij/sqrt(aii*ajj)
  call MatCreateVecs(A, diag, sqrt_diag, ierr)
  call MatGetDiagonal(A, diag, ierr)
  call VecMin(diag, diagminloc, diagmin, ierr)
  if (diagmin<=0.0) then
    ewrite(0,*) 'Multigrid preconditioner "mg" requires strictly positive diagonal'
    FLExit("Zero or negative value on the diagonal")
  end if
  
  !
  call VecCopy(diag, sqrt_diag, ierr)
  call VecSqrtAbs(sqrt_diag, ierr)
  !
  call VecDuplicate(sqrt_diag, inv_sqrt_diag, ierr)
  call VecCopy(sqrt_diag, inv_sqrt_diag, ierr)
  call VecReciprocal(inv_sqrt_diag, ierr)
  call MatDiagonalScale(A, inv_sqrt_diag, inv_sqrt_diag, ierr)
  
  ! construct the strongly coupled neighbourhoods N_i around each node
  ! and use R to register isolated nodes i with N_i={i}
  call Prolongator_init(R, ccnt, findN, N, A, base, epsilon)
  
  if (ccnt==0 .and. nrows>0) then
    ! all nodes are isolated, strongly diagonal dominant matrix
    ! we should solve by other means
    if (present(cluster)) cluster=ISOLATED
    deallocate(findN, N, R)
    ! we return PETSC_NULL; callers of this function should check for this
    P=PETSC_NULL_OBJECT
    return
  else if (100*ccnt<99*nrows .and. .not. IsParallel()) then
    ! more than 1% isolated nodes, give a warning
    ewrite(2,*) "Percentage of isolated nodes: ", (100.0*(nrows-ccnt))/nrows
    ewrite(2,*) "Warning: more than 1 perc. isolated nodes - this may mean mg is not the most suitable preconditioner"
    ewrite(2,*) "On small meshes with a lot of boundary nodes, this is typically fine though."
  end if

  ! Step 1 - Startup aggregation
  ! select some of the coupled neighbourhoods as an initial (incomplete) 
  !   covering
  call Prolongator_step1(R, jc, findN, N, maxclustersize)
  
  ! Step 2 - Enlarging the decomposition sets (aggregates)
  ! add remaining COUPLED but yet uncovered nodes to one of the aggregates
  call Prolongator_step2(R, findN, N, A, base)
  
  ! Step 3 - Handling the remnants
  ! the remaining nodes, that are COUPLED but neither in the original covering
  ! or assigned in step 2, are assigned to new aggregates
  call Prolongator_step3(R, jc, findN, N)
  
  ! jc is now the n/o aggregates, i.e. the n/o coarse nodes
  ! R(i) is now either the coarse node, fine node i is assigned to
  !             or ==ISOLATED
  
  ewrite(3,*) 'Fine nodes: ', nrows
  ewrite(3,*) 'Isolated fine nodes: ', nrows-ccnt
  ewrite(3,*) 'Aggregates: ', jc
  
  ! now scale a_ij -> a_ii^-1/2 * a_ij * ajj^1/2, i.e. starting from the
  ! original matrix: a_ij -> a_ii^-1 * a_ij
  call MatDiagonalScale(A, inv_sqrt_diag, sqrt_diag, ierr)
  
  ! we now have all the stuff to create the prolongator
  call create_prolongator(P, nrows, jc, findN, N, R, A, base, omega)
  
  ! now restore the original matrix
  ! unfortunately MatDiagonalScale is broken for one-sided scaling, i.e.
  ! supplying PETSC_NULL(_OBJECT) for one the vectors
  call VecDuplicate(diag, one, ierr)
  call VecSet(one, 1.0_PetscReal_kind, ierr)
  call MatDiagonalScale(A, diag, one, ierr)
   
  if (present(cluster)) cluster=R
  deallocate(R, N, findN)
  
  call VecDestroy(diag, ierr)
  call VecDestroy(sqrt_diag, ierr)
  call VecDestroy(inv_sqrt_diag, ierr)
  call VecDestroy(one, ierr)
    
end function Prolongator

subroutine create_prolongator(P, nrows, ncols, findN, N, R, A, base, omega)
  
  Mat, intent(out):: P
  integer, intent(in):: nrows ! number of fine nodes
  integer, intent(in):: ncols ! number of clusters
  integer, dimension(:), intent(in):: findN, N, R
  ! A needs to be left rescaled with the inverse diagonal: D^-1 A
  Mat, intent(in):: A
  integer, intent(in):: base
  PetscReal, intent(in):: omega
  
  PetscObject:: myPETSC_NULL_OBJECT
  PetscErrorCode:: ierr
  Vec:: rowsum_vec
  PetscReal, dimension(:), allocatable:: Arowsum
  PetscReal:: aij(1), rowsum
  integer, dimension(:), allocatable:: dnnz, onnz
  integer:: i, j, k, coarse_base
  
  allocate(dnnz(1:nrows), Arowsum(1:nrows))
  
  ! work out nnz in each row of the new prolongator
  dnnz=0
  do i=1, nrows
    ! this is an overestimate it should count the number
    ! of different R(j) values in each rows
    do k=findN(i), findN(i+1)-1
      j=N(k)
      if (R(j)>0) then
        dnnz(i)=dnnz(i)+1
      end if
    end do
    ! since we overestimate, we don't want to get > ncols
    dnnz(i)=min(dnnz(i), ncols)
  end do      

  if (IsParallel()) then
    ! for the moment the prolongator is completely local:
    allocate(onnz(1:nrows))
    onnz=0
    
    call MatCreateAIJ(MPI_COMM_FEMTOOLS, nrows, ncols, PETSC_DECIDE, PETSC_DECIDE, &
      PETSC_NULL_INTEGER, dnnz, PETSC_NULL_INTEGER, onnz, P, ierr)
    call MatSetOption(P, MAT_USE_INODES, PETSC_FALSE, ierr)
      
    ! get base for coarse node/cluster numbering
    call MatGetOwnerShipRangeColumn(P, coarse_base, PETSC_NULL_INTEGER, ierr)
    ! subtract 1 to convert from 1-based fortran to 0 based petsc
    coarse_base=coarse_base-1
  else
    call MatCreateAIJ(MPI_COMM_SELF, nrows, ncols, nrows, ncols, &
      PETSC_NULL_INTEGER, dnnz, 0, PETSC_NULL_INTEGER, P, ierr)
    call MatSetOption(P, MAT_USE_INODES, PETSC_FALSE, ierr)
    ! subtract 1 from each cluster no to get petsc 0-based numbering
    coarse_base=-1
  end if
  call MatSetup(P, ierr)
  
  myPETSC_NULL_OBJECT=PETSC_NULL_OBJECT
  call MatCreateVecs(A, rowsum_vec, PETSC_NULL_OBJECT, ierr)
  if (myPETSC_NULL_OBJECT/=PETSC_NULL_OBJECT) then
    FLAbort("PETSC_NULL_OBJECT has changed please report to skramer")
  end if
  call VecPlaceArray(rowsum_vec, Arowsum, ierr)
  call MatGetRowSum(A, rowsum_vec, ierr)
    
  do i=1, nrows
    rowsum=0.0
    ! the filtered matrix only contains the entries in N_i:
    do k=findN(i), findN(i+1)-1
      j=N(k)
      call MatGetValues(A, 1, (/ base+i /), 1, (/ base+j /),  aij, ierr)
      rowsum=rowsum+aij(1)
      if (R(j)>0) then
        call MatSetValue(P, base+i, coarse_base+R(j), -omega*aij(1), &
          ADD_VALUES, ierr)
      end if
    end do
    if (R(i)>0) then
      call MatSetValue(P, base+i, coarse_base+R(i), &
        1+omega*( rowsum-Arowsum(i) ), ADD_VALUES, ierr)
    end if
  end do
    
  call MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY, ierr)
  
  call VecDestroy(rowsum_vec, ierr)
    
end subroutine create_prolongator

subroutine Prolongator_init(R, ccnt, findN, N, A, base, epsilon)
! construct the strongly coupled neighbourhoods N_i around each node
! and use R to register isolated nodes i with N_i={i}
integer, dimension(:), intent(out):: R, findN, N
integer, intent(out):: ccnt
Mat, intent(in):: A
integer, intent(in):: base
PetscReal, intent(in):: epsilon

  PetscErrorCode:: ierr
  PetscReal, dimension(:), allocatable:: vals(:)
  integer, dimension(:), allocatable:: cols(:)
  PetscReal aij, eps_sqrt
  integer i, j, k, p, ncols
  
  ! workspace for MatGetRow
  allocate( vals(1:size(N)), cols(1:size(n)) )
  
  eps_sqrt=sqrt(epsilon)
  
  ccnt=0 ! counts the coupled nodes, i.e. nodes that are not isolated
  p=1
  do i=1, size(R)
    findN(i)=p
    call MatGetRow(A, base+i, ncols, cols, vals, ierr)
    do k=1, ncols
      j=cols(k)-base
      ! ignore non-local columns
      if (j<1 .or. j>size(R)) cycle
      aij=vals(k)
      if (abs(aij)>eps_sqrt) then
        N(p)=j
        p=p+1
      end if
    end do
    ! N_i always contains i itself, so p has to be increased with at least 2
    if (p>findN(i)+1) then
      R(i)=COUPLED
      ccnt=ccnt+1
    else
      R(i)=ISOLATED
    end if
    call MatRestoreRow(A, base+i, ncols, cols, vals, ierr)
  end do
  findN(i)=p
  
end subroutine Prolongator_init

subroutine Prolongator_step1(R, jc, findN, N, maxclustersize)
integer, dimension(:), intent(inout):: R
integer, intent(out):: jc
integer, dimension(:), intent(in):: findN, N
integer, intent(in):: maxclustersize
  ! Step 1 - Startup aggregation
  ! select some of the coupled neighbourhoods as an initial (incomplete) 
  !   covering
  
  integer, dimension(:), allocatable:: clustersize
  integer i, j
  
  allocate(clustersize(1:size(R)))
  clustersize=0
  jc=0 ! count the covering sets, aka aggregates
  do i=1, size(R)
    ! if the entire N_i is still a subset of R:
    if (R(i)==COUPLED) then
      if (all(R(N(findN(i):findN(i+1)-1))==COUPLED)) then
        jc=jc+1 ! new aggregate
        ! assign all of N_i to aggr. jc and remove from R:
        R(N(findN(i):findN(i+1)-1))=jc
        clustersize(jc)=findN(i+1)-findN(i)
      else
        ! here we deviate from Vanek'96 by allowing this node to be added
        ! to existing clusters, as long as maxclustersize isn't exceeded
        do j=findN(i), findN(i+1)-1
          if (N(j)==i) cycle
          if (R(N(j))>0) then
            if (clustersize(R(N(j)))<maxclustersize) exit
          end if
        end do
        if (j<findN(i+1)) then
          R(i)=R(N(j))
          clustersize(R(N(j)))=clustersize(R(N(j)))+1
        end if
      end if
    end if
  end do
  deallocate(clustersize)
  
end subroutine Prolongator_step1

subroutine Prolongator_step2(R, findN, N, A, base)
! Step 2 - Enlarging he decomposition sets
! add remaining COUPLED but yet uncovered nodes to one of the aggregates
integer, dimension(:), intent(inout):: R
integer, dimension(:), intent(in):: findN, N
Mat, intent(in):: A
integer, intent(in):: base

  PetscErrorCode:: ierr
  PetscReal:: maxc, aij(1)
  integer:: i, j, k, p
  
  do i=1, size(R)
    if (R(i)==COUPLED) then
      ! find the strongest coupling in N_i that is assigned to one of the
      ! aggregates in step 1
      maxc=0
      k=0
      do p=findN(i), findN(i+1)-1
        j=N(p)
        if (R(j)>0) then
          j=N(p)
          call MatGetValues(A, 1, (/ base+i /), 1, (/ base+j /),  aij, ierr)
          if (abs(aij(1))>maxc) then
            maxc=abs(aij(1))
            k=R(j)
          end if
        end if
      end do
      ! remove i from R, register its assignment to aggregate k
      if (k/=0) then
        R(i)=-k+COUPLED ! the assignment is temp. registered with a negative 
          ! number, so as not to actually change the aggregates of step 1
          ! before completion of step 2, the negative numbers should be below
          ! the value of COUPLED
      end if
    end if
  end do

end subroutine Prolongator_step2

subroutine Prolongator_step3(R, jc, findN, N)
! Step 3 - Handling the remnants
! the remaining nodes, that are COUPLED but neither in the original covering
! or assigned in step 2, are assigned to new aggregates
integer, dimension(:), intent(inout):: R
integer, intent(inout):: jc
integer, dimension(:), intent(in):: findN, N

  integer i, j, p
  
  do i=1, size(R)
    if (R(i)==COUPLED) then
      ewrite(3,*) 'step3 action! ',i, ':', N(findN(i):findN(i+1)-1)
      jc=jc+1 ! add new aggregate
      do p=findN(i), findN(i+1)-1
        j=N(p)
        if (R(j)==COUPLED) then
          R(j)=jc
        end if
      end do
    else if (R(i)<COUPLED) then
      ! definitively assign nodes from step 2 to their aggregates:
      R(i)=-(R(i)-COUPLED)
    end if
  end do

end subroutine Prolongator_step3
  
subroutine PowerMethod(matrix, eigval, eigvec)
Mat, intent(in):: matrix
PetscReal, intent(out):: eigval
Vec, intent(out):: eigvec

  integer, parameter:: MAX_ITERATIONS=100
  PetscReal, parameter:: TOLERANCE=0.001
  Vec:: x_k, x_kp1
  PetscErrorCode:: ierr
  PetscReal, dimension(2):: dot_prods
  PetscReal:: rho_k, rho_kp1, norm2
  integer:: i
  PetscRandom:: pr
  PetscObject:: myPETSC_NULL_OBJECT  
  
  call MatCreateVecs(matrix, x_kp1, x_k, ierr)
  
  ! initial guess
  call PetscRandomCreate(PETSC_COMM_WORLD, pr, ierr)
  myPETSC_NULL_OBJECT=PETSC_NULL_OBJECT
  call VecSetRandom(x_k, PETSC_NULL_OBJECT, ierr)
  if (myPETSC_NULL_OBJECT/=PETSC_NULL_OBJECT) then
    FLAbort("PETSC_NULL_OBJECT has changed please report to skramer")
  end if
  call PetscRandomDestroy(pr, ierr)

  rho_k=0.0
  
  do i=1, MAX_ITERATIONS
    call  MatMult(matrix, x_k, x_kp1, ierr)
    ! compute < x_k+1, x_k+1 > and < x_k+1, x_k >
    call VecMDot(x_kp1, 2, (/ x_kp1, x_k /), dot_prods, ierr)
    norm2=sqrt(dot_prods(1))
    rho_kp1=dot_prods(2)
    
    call VecScale(x_kp1, 1.0_PetscReal_kind/norm2, ierr)
    
    ! convergence criterium
    if (abs(rho_kp1-rho_k)<TOLERANCE*rho_kp1) exit    
    
    ! copy for next iteration
    call VecCopy(x_kp1, x_k, ierr)
    rho_k=rho_kp1
  end do
  
  if (i>MAX_ITERATIONS) then
    FLAbort("PowerMethod failed to converge")
  end if
  
  eigval=rho_kp1
  eigvec=x_kp1
  
  call VecDestroy(x_k, ierr)  
  
end subroutine PowerMethod

end module multigrid
