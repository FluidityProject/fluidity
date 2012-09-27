#include "fdebug.h"
!! This module contains multigrid related subroutines, such as the smoothed
!! aggregation preconditioner.
module multigrid_interface
use Petsc_Tools
use Sparse_tools
use sparse_tools_petsc
use FLDebug
use spud
use futils
use parallel_tools
use multigrid
#include "petscversion.h"
#ifdef HAVE_PETSC_MODULES
  use petsc
#if PETSC_VERSION_MINOR==0
  use petscvec
  use petscmat
  use petscksp
  use petscpc
  use petscis
  use petscmg
#endif
#endif
implicit none
#ifdef HAVE_PETSC_MODULES
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscsysdef.h"
#else
#include "finclude/petsc.h"
#if PETSC_VERSION_MINOR==0
#include "finclude/petscmat.h"
#include "finclude/petscvec.h"
#include "finclude/petscviewer.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscsys.h"
#endif
#endif
#if PETSC_VERSION_MINOR==2
#define KSP_NORM_NO KSP_NORM_NONE
#endif

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
  
contains

subroutine SetupMultigrid(prec, matrix, ierror, &
  external_prolongators, surface_node_list, matrix_csr, &
  internal_smoothing_option, has_null_space)
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
!! option to prevent a direct solve at the coarsest level
logical, optional, intent(in) :: has_null_space

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
          external_prolongators=external_prolongators, &
          has_null_space=has_null_space)
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
          external_prolongators, no_top_smoothing=.true., &
          has_null_space=has_null_space)
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
          external_prolongators=external_prolongators, &
          has_null_space=has_null_space)
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

  call PCSetOperators(internal_smoother_pc,Internal_Smoother_Mat, &
       Internal_Smoother_Mat,DIFFERENT_NONZERO_PATTERN,ierr)

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

subroutine SetupSmoothedAggregation(prec, matrix, ierror, &
  external_prolongators,no_top_smoothing, has_null_space)
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
!! option to prevent a direct solve at the coarsest level
logical, optional, intent(in) :: has_null_space


  if (present(external_prolongators)) then
    call setup_mg(prec, matrix, ierror, &
       external_prolongators=external_prolongators%M, &
       no_top_smoothing=no_top_smoothing, &
       has_null_space=has_null_space)
  else
    call setup_mg(prec, matrix, ierror, &
       no_top_smoothing=no_top_smoothing, &
       has_null_space=has_null_space)
  end if

end subroutine SetupSmoothedAggregation

end module multigrid_interface
