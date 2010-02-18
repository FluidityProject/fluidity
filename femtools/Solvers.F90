!  Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    C.Pain@Imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module solvers
  use FLDebug
  use elements
  use Petsc_tools
  use Signal_Vars
  use Multigrid
  use Sparse_Tools
  use sparse_tools_petsc
  use sparse_matrices_fields
  use Fields
  use Global_Parameters
  use spud
  use halos
#ifdef HAVE_PETSC_MODULES
  use petsc 
  use petscvec 
  use petscmat 
  use petscksp 
  use petscpc
#endif
  implicit none
  ! Module to provide explicit interfaces to matrix solvers.

#ifdef HAVE_PETSC
#ifdef HAVE_PETSC_MODULES
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
#else
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#endif
#else
#define PetscReal real
#define PetscInt integer
#define KSPGMRES "gmres"
#define PCSOR "sor"
#endif

#ifdef HAVE_PETSC
  ! stuff used in the PETSc monitor (see petsc_solve_callback_setup() below)
  Vec :: petsc_monitor_exact
  Vec :: petsc_monitor_x
  real, dimension(:), pointer :: petsc_monitor_error => null()
  PetscLogDouble, dimension(:), pointer :: petsc_monitor_flops => null()
  integer :: petsc_monitor_iteration = 0
#endif
    
private

public solcg, gmres, petsc_solve, set_solver_options, &
   complete_solver_option_path

! meant for unit-testing solver code only:
public petsc_solve_setup, petsc_solve_core, petsc_solve_destroy, &
  petsc_solve_copy_vectors_from_scalar_fields, &
  setup_ksp_from_options, SetupKSP

interface petsc_solve
   module procedure petsc_solve_scalar, petsc_solve_vector, &
     petsc_solve_scalar_multiple, petsc_solve_scalar_multiple_rhs_array, &
     petsc_solve_vector_components, &
     petsc_solve_tensor_components, &
     petsc_solve_csr, petsc_solve_block_csr, petsc_solve_block_csr_multiple, &
     petsc_solve_csr_multiple, &
     petsc_solve_scalar_petsc_csr, petsc_solve_vector_petsc_csr
end interface
  
interface set_solver_options
    module procedure set_solver_options_with_path, &
      set_solver_options_scalar, set_solver_options_vector, set_solver_options_tensor
end interface set_solver_options

contains

subroutine solcg(xk,f,nonods,fredop,nnodp,zero,  &
     a,fina,cola,ncola,nbigm,        &
     halo_tag,k, &
     option_path, checkconvergence)
  use parallel_tools
  integer, intent(in):: nonods, fredop, nnodp,  nbigm, ncola
  real, intent(inout):: xk(fredop) ! IN:  initial guess (if ZERO==.false.)
                                   ! OUT: solution vector
  real, intent(in):: f(fredop)  ! rhs of equation
  logical zero ! if zero==.true. use zero initial guess
  real, intent(in):: A(nbigm) ! values of the matrix
  integer, intent(in):: fina(nonods+1) ! row structure of matrix
  integer, intent(in):: cola(ncola)    ! column indices of matrix
  integer, intent(in):: halo_tag ! for parallel computations
  integer, intent(out):: k ! returns number of performed iterations
  character(len=*), optional, intent(in):: option_path
  logical, optional, intent(in):: checkconvergence ! if .true. (default)
    ! fail if not converged
  type(csr_sparsity) :: sparsity
    
#ifdef HAVE_PETSC
  type(block_csr_matrix) matrix
  integer :: idimm, lhalo_tag, stat
  type(halo_type) :: halo
#endif

#ifdef HAVE_PETSC
  idimm = fredop / nonods

  if (IsParallel()) then
    ! Halos are stored on meshes / sparsities. We don't have sufficient
    ! information for a parallel solve.
    FLAbort("Cannot solve using solcg in parallel")
  else
    lhalo_tag=0
  end if
  !call import_halo(lhalo_tag, halo,stat)
  stat = -1

  if (stat==0) then
     sparsity=wrap(fina, colm=cola, name='TemporarySparsity_solcg', &
          row_halo=halo, column_halo=halo)
     matrix=wrap(sparsity, blocks=(/ idimm, idimm /), val=a, &
          name='TemporaryMatrix_solcg')
     call deallocate(halo)
  else
     sparsity=wrap(fina, colm=cola, name='TemporarySparsity_solcg')
     matrix=wrap(sparsity, blocks=(/ idimm, idimm /), val=a, &
          name='TemporaryMatrix_solcg')
  end if
  call deallocate(sparsity)
  call petsc_solve_block_csr(xk, matrix, f, &
    startfromzero=ZERO, checkconvergence=checkconvergence, &
    iterations=k, option_path=option_path)
  call deallocate(matrix)
  ewrite(2,*) 'Out of petsc_(block_)solve_csr.'
#endif
  
end subroutine solcg

subroutine gmres(x,b,nonods,fredop,nnodp,zero, &
     a,fina,cola,ncola,nbigm, &
     halo_tag,k, &
     option_path, checkconvergence)
  use parallel_tools
  integer nnodp
  integer fredop,ncola,nonods,nbigm
  integer fina(nonods+1),cola(ncola)
  real a(nbigm)
  integer halo_tag,k
  real x(fredop),b(fredop)
  logical zero
  character(len=*), optional, intent(in) :: option_path
  logical, optional, intent(in) :: checkconvergence
  type(csr_sparsity) :: sparsity
  
#ifdef HAVE_PETSC
  type(block_csr_matrix) matrix
  integer idimm, lhalo_tag, stat
  type(halo_type) halo

  idimm = fredop / nonods

  if (IsParallel()) then
    ! Halos are stored on meshes / sparsities. We don't have sufficient
    ! information for a parallel solve.
    FLAbort("Cannot solve using gmres in parallel")
  else
    lhalo_tag=0
  end if
  !call import_halo(lhalo_tag, halo, stat)
  stat = -1
  
  if (stat==0) then
     sparsity=wrap(fina, colm=cola, &
          name='TemporarySparsity_gmres', &
          row_halo=halo, column_halo=halo)
     matrix=wrap(sparsity, blocks=(/ idimm, idimm /), val=a, &
          name='TemporaryMatrix_gmres')
     call deallocate(halo)
  else
     sparsity=wrap(fina, colm=cola, &
          name='TemporarySparsity_gmres')
     matrix=wrap(sparsity, blocks=(/ idimm, idimm /), val=a, &
          name='TemporaryMatrix_gmres')
     call deallocate(halo)
  end if

  call deallocate(sparsity)
  call petsc_solve_block_csr(x, matrix, b, &
    startfromzero=zero, checkconvergence=checkconvergence, &
    iterations=k, &
    option_path=option_path)
  call deallocate(matrix)
  ewrite(2,*) 'Out of petsc_(block_)solve_csr.'
#endif

end subroutine gmres

subroutine petsc_solve_scalar(x, matrix, rhs, option_path, &
  preconditioner_matrix, prolongator, surface_node_list,exact, &
  error_filename,internal_smoothing_option)
  !!< Solve a linear system the nice way.
  type(scalar_field), intent(inout) :: x
  type(scalar_field), intent(in) :: rhs
  type(csr_matrix), intent(in) :: matrix
  character(len=*), optional, intent(in) :: option_path
  !! 2 experimental arguments to improve preconditioning with extra outside information
  !! provide approximation the matrix (only to be used in combination with pctype='KSP')
  type(csr_matrix), optional, intent(in) :: preconditioner_matrix
  !! prolongator to be used at the first level of 'mg'
  type(csr_matrix), optional, intent(in) :: prolongator
  !! exact solution, for testing
  type(scalar_field), optional, intent(in) :: exact
  !! file to dump errors and flops
  character(len=*), optional, intent(in) :: error_filename
  !! surface_node_list for internal smoothing
  integer, dimension(:), optional, intent(in) :: surface_node_list
  !! internal smoothing option
  integer, intent(in), optional :: internal_smoothing_option

#ifdef HAVE_PETSC
  KSP ksp
  Mat A
  Vec y, b

  type(petsc_numbering_type) petsc_numbering
  character(len=OPTION_PATH_LEN) solver_option_path, name, option_path_in
  integer literations
  logical lstartfromzero
  
  assert(size(x%val)==size(rhs%val))
  assert(size(x%val)==size(matrix,2))
#ifdef DDEBUG
  if (.not.associated(matrix%sparsity%column_halo)) then
    assert(size(rhs%val)==size(matrix,1))
  else
    ! in parallel we allow for the matrix to not have the halo nodes
    ! whereas the rhs will always contain them (even if not used)
    if (size(matrix,1)/=size(rhs%val)) then
      ! get_nowned_nodes() seems completely buggered for T10 halo
      ! this should be fixed by dham+jrmaddison's mesh+halo integration
      !assert(size(matrix,1)==get_nowned_nodes(matrix%sparsity%halo_tag))
      assert( size(matrix,1)==halo_nowned_nodes(matrix%sparsity%column_halo) )
    end if
  end if
#endif
  
  ! option_path_in may still point to field 
  ! (so we have to add "/prognostic/solver" below)
  if (present(option_path)) then
    option_path_in=option_path
  else
    option_path_in=x%option_path
  end if
  
  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
        name, solver_option_path, lstartfromzero, &
        matrix=matrix, &
        option_path=option_path_in, &
        preconditioner_matrix=preconditioner_matrix, &
        prolongator=prolongator,surface_node_list=surface_node_list, &
        exact=exact,internal_smoothing_option=internal_smoothing_option)
 
  ! copy array into PETSc vecs
  call petsc_solve_copy_vectors_from_scalar_fields(y, b, x, &
       & matrix, rhs, petsc_numbering, lstartfromzero)
     
  ! the solve and convergence check
  call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
        name, solver_option_path, lstartfromzero, &
        literations, x0=x%val)
        
  ! Copy back the result using the petsc numbering:
  call petsc2field(y, petsc_numbering, x, rhs)
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy(y, A, b, ksp, petsc_numbering, exact, &
       & error_filename)
  
#else
  FLAbort("Petsc_solve called while not configured with PETSc")
#endif
  
end subroutine petsc_solve_scalar

subroutine petsc_solve_scalar_multiple(x, matrix, rhs, option_path)
  !!< Solves multiple scalar fields with the same matrix.
  !!< Need to specify an option_path as there's no default
  type(scalar_field), dimension(:), intent(inout) :: x
  type(scalar_field), dimension(:), intent(in) :: rhs
  type(csr_matrix), intent(in) :: matrix
  character(len=*), optional, intent(in) :: option_path

#ifdef HAVE_PETSC
  KSP ksp
  Mat A
  Vec y, b

  type(petsc_numbering_type) petsc_numbering
  character(len=OPTION_PATH_LEN) solver_option_path, name
  integer literations
  logical lstartfromzero
  integer i
  
  assert(size(x)==size(rhs))
  do i=1, size(x)
    assert(size(x(i)%val)==size(rhs(i)%val))
    assert(size(x(i)%val)==size(matrix,2))
    assert(size(rhs(i)%val)==size(matrix,1))
  end do
  
  ewrite(1,*) 'Solving for multiple scalar fields at once'

  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
        name, solver_option_path, lstartfromzero, &
        matrix=matrix, &
        option_path=option_path)
  
  do i=1, size(x)
 
    ! copy array into PETSc vecs
    call petsc_solve_copy_vectors_from_scalar_fields(y, b, &
         x(i), matrix, rhs(i), &
         petsc_numbering, lstartfromzero)
     
    ! the solve and convergence check
    call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
        name, solver_option_path, lstartfromzero, &
        literations, x0=x(i)%val)
        
    ! Copy back the result using the petsc numbering:
    call petsc2field(y, petsc_numbering, x(i), rhs(i))
    
  end do
    
  ewrite(1,*) 'Finished solving all scalar fields'
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy(y, A, b, ksp, petsc_numbering)
  
#else  
  FLAbort("Petsc_solve called while not configured with PETSc")
#endif
  
end subroutine petsc_solve_scalar_multiple

subroutine petsc_solve_scalar_multiple_rhs_array(x, matrix, rhs, option_path)
  !!< Solve a linear system the nice way.
  !!< Solves multiple scalar fields with the same matrix.
  !!< This version allows to specify the rhs as an array
  type(scalar_field), dimension(:), intent(inout) :: x
  real, dimension(:,:), intent(in) :: rhs
  type(csr_matrix), intent(in) :: matrix
  character(len=*), optional, intent(in) :: option_path

#ifdef HAVE_PETSC

  type(scalar_field), dimension(1:size(rhs,2)):: rhs_fields
  integer i
  
  do i=1, size(rhs,2)
    rhs_fields(i)=wrap_scalar_field(x(i)%mesh, rhs(:,i), &
      name='petsc_solve_scalar_multiple_rhs_array'//int2str(i))
  end do
  
  call petsc_solve_scalar_multiple(x, matrix, rhs_fields, &
     option_path=option_path)
  
  do i=1, size(rhs,2)
    call deallocate(rhs_fields(i))
  end do
  
#else
  FLAbort("Petsc_solve called while not configured with PETSc")
#endif
  
end subroutine petsc_solve_scalar_multiple_rhs_array

subroutine petsc_solve_vector(x, matrix, rhs, option_path, deallocate_matrix)
  !!< Solve a linear system the nice way. Options for this
  !!< come via the options mechanism. 
  type(vector_field), intent(inout) :: x
  type(vector_field), intent(in) :: rhs
  type(block_csr_matrix), intent(inout) :: matrix
  character(len=*), optional, intent(in) :: option_path
  !! deallocate the matrix after it's been copied
  logical, intent(in), optional :: deallocate_matrix

#ifdef HAVE_PETSC
  KSP ksp
  Mat A
  Vec y, b

  type(petsc_numbering_type) petsc_numbering
  character(len=OPTION_PATH_LEN) solver_option_path, name, option_path_in
  integer literations
  logical lstartfromzero
  
  type(csr_matrix) :: matrixblock
  type(scalar_field) :: rhsblock, xblock
  integer :: i
  
  assert(x%dim==rhs%dim)
  assert(size(x%val(1)%ptr)==size(rhs%val(1)%ptr))
  assert(size(x%val(1)%ptr)==block_size(matrix,2))
  assert(size(rhs%val(1)%ptr)==block_size(matrix,1))
  assert(x%dim==blocks(matrix,2))
  assert(rhs%dim==blocks(matrix,1))
  
  ! option_path_in may still point to field 
  ! (so we have to add "/prognostic/solver" below)
  if (present(option_path)) then
    option_path_in=option_path
  else
    option_path_in=x%option_path
  end if
  
  if(matrix%diagonal) then
    assert(blocks(matrix,1)==blocks(matrix,2))
    ! only want to solve using the diagonal blocks
    do i = 1, blocks(matrix,1)
      matrixblock=block(matrix,i,i)
      rhsblock = extract_scalar_field(rhs, i)
      xblock = extract_scalar_field(x, i)
  
      ! setup PETSc object and petsc_numbering from options and 
      call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
            name, solver_option_path, lstartfromzero, &
            matrix=matrixblock, &
            option_path=option_path_in)
    
      ! copy array into PETSc vecs
      call petsc_solve_copy_vectors_from_scalar_fields(y, b, xblock, &
          & matrixblock, rhsblock, petsc_numbering, lstartfromzero)
      
      if(present_and_true(deallocate_matrix).and.(i==blocks(matrix,1))) then
        call deallocate(matrix)
      end if
        
      ! the solve and convergence check
      call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
            name, solver_option_path, lstartfromzero, &
            literations, x0=xblock%val)
            
      ! Copy back the result using the petsc numbering:
      call petsc2field(y, petsc_numbering, xblock, rhsblock)
      
      ! destroy all PETSc objects and the petsc_numbering
      call petsc_solve_destroy(y, A, b, ksp, petsc_numbering)
          
    end do
    
  else
  
    ! setup PETSc object and petsc_numbering from options and 
    call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
          name, solver_option_path, lstartfromzero, &
          block_matrix=matrix, &
          option_path=option_path_in)
          
    if(present_and_true(deallocate_matrix)) then
      call deallocate(matrix)
    end if
  
    ! copy array into PETSc vecs
    call petsc_solve_copy_vectors_from_vector_fields(y, b, x, rhs, petsc_numbering, lstartfromzero)
      
    ! the solve and convergence check
    call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
            name, solver_option_path, lstartfromzero, &
            literations, vector_x0=x)
          
    ! Copy back the result using the petsc numbering:
    call petsc2field(y, petsc_numbering, x)
    
    ! destroy all PETSc objects and the petsc_numbering
    call petsc_solve_destroy(y, A, b, ksp, petsc_numbering)
  end if
  
#else
  FLAbort("Petsc_solve called while not configured with PETSc")
#endif
  
end subroutine petsc_solve_vector
  
subroutine petsc_solve_vector_components(x, matrix, rhs, option_path)
  !!< Solve a linear system the nice way. Options for this
  !!< come via the options mechanism. This version solves a linear system
  !!< for each of the components of rhs each time with the same matrix.
  type(vector_field), intent(inout) :: x
  type(vector_field), intent(in) :: rhs
  type(csr_matrix), intent(in) :: matrix
  character(len=*), optional, intent(in) :: option_path

#ifdef HAVE_PETSC
  KSP ksp
  Mat A
  Vec y, b

  type(scalar_field) x_component, rhs_component
  type(petsc_numbering_type) petsc_numbering
  character(len=OPTION_PATH_LEN) solver_option_path, name, option_path_in
  integer literations, i
  logical lstartfromzero
  
  assert(x%dim==rhs%dim)
  assert(size(x%val(1)%ptr)==size(rhs%val(1)%ptr))
  assert(size(x%val(1)%ptr)==size(matrix,2))
  assert(size(rhs%val(1)%ptr)==size(matrix,1))
  
  ! option_path_in may still point to field 
  ! (so we have to add "/prognostic/solver" below)
  if (present(option_path)) then
    option_path_in=option_path
  else
    option_path_in=x%option_path
  end if
  
  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
        name, solver_option_path, lstartfromzero, &
        matrix=matrix, &
        option_path=option_path_in)
 
  ewrite(1,*) 'Solving for multiple components of a vector field'
  
  do i=1, x%dim
    
    ewrite(1, *) 'Now solving for component: ', i
     x_component=extract_scalar_field(x, i)
     rhs_component=extract_scalar_field(rhs, i)
     ! copy array into PETSc vecs
     call petsc_solve_copy_vectors_from_scalar_fields(y, b, x_component, matrix, rhs_component, petsc_numbering, lstartfromzero)
     
     ! the solve and convergence check
     call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
          name, solver_option_path, lstartfromzero, &
          literations, x0=x_component%val)
        
     ! Copy back the result using the petsc numbering:
     call petsc2field(y, petsc_numbering, x_component, rhs_component)
     
  end do
  
  ewrite(1,*) 'Finished solving all components.'
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy(y, A, b, ksp, petsc_numbering)
  
#else
  FLAbort("Petsc_solve called while not configured with PETSc")
#endif
  
end subroutine petsc_solve_vector_components

subroutine petsc_solve_scalar_petsc_csr(x, matrix, rhs, option_path, &
  prolongator, exact, surface_node_list)
  !!< Solve a linear system the nice way. Options for this
  !!< come via the options mechanism. 
  type(scalar_field), intent(inout) :: x
  type(scalar_field), intent(in) :: rhs
  type(petsc_csr_matrix), intent(inout) :: matrix
  character(len=*), optional, intent(in) :: option_path
  !! prolongator to be used at the first level of 'mg'
  type(csr_matrix), optional, intent(in) :: prolongator
  !! exact solution, for testing
  type(scalar_field), optional, intent(in) :: exact
  !! surface_node_list for internal smoothing
  integer, dimension(:), optional, intent(in) :: surface_node_list

#ifdef HAVE_PETSC
  KSP ksp
  Vec y, b

  character(len=OPTION_PATH_LEN) solver_option_path, name, option_path_in
  integer literations
  logical lstartfromzero
  
  assert(size(x%val)==size(rhs%val))
  assert(size(x%val)==size(matrix,2))
  assert(size(rhs%val)==size(matrix,1))
  
  ! option_path_in may still point to field 
  ! (so we have to add "/prognostic/solver" below)
  if (present(option_path)) then
    option_path_in=option_path
  else
    option_path_in=x%option_path
  end if
    
  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup_petsc_csr(y, b, ksp, &
        name, solver_option_path, lstartfromzero, &
        matrix, option_path=option_path_in, &
        prolongator=prolongator,surface_node_list=surface_node_list, &
        exact=exact)
        
  ! copy array into PETSc vecs
  call petsc_solve_copy_vectors_from_scalar_fields(y, b, x, rhs=rhs, &
     petsc_numbering=matrix%row_numbering, startfromzero=lstartfromzero)
    
  ! the solve and convergence check
  call petsc_solve_core(y, matrix%M, b, ksp, matrix%row_numbering, &
          name, solver_option_path, lstartfromzero, &
          literations, x0=x%val)
        
  ! Copy back the result using the petsc numbering:
  call petsc2field(y, matrix%column_numbering, x)
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy_petsc_csr(y, b, ksp)
    
#else
  FLAbort("Petsc_solve called while not configured with PETSc")
#endif
  
end subroutine petsc_solve_scalar_petsc_csr

subroutine petsc_solve_vector_petsc_csr(x, matrix, rhs, option_path)
  !!< Solve a linear system the nice way. Options for this
  !!< come via the options mechanism. 
  type(vector_field), intent(inout) :: x
  type(vector_field), intent(in) :: rhs
  type(petsc_csr_matrix), intent(inout) :: matrix
  character(len=*), optional, intent(in) :: option_path

#ifdef HAVE_PETSC
  KSP ksp
  Vec y, b

  character(len=OPTION_PATH_LEN) solver_option_path, name, option_path_in
  integer literations
  logical lstartfromzero
  
  assert(x%dim==rhs%dim)
  assert(size(x%val(1)%ptr)==size(rhs%val(1)%ptr))
  assert(size(x%val(1)%ptr)==block_size(matrix,2))
  assert(size(rhs%val(1)%ptr)==block_size(matrix,1))
  assert(x%dim==blocks(matrix,2))
  assert(rhs%dim==blocks(matrix,1))
  
  ! option_path_in may still point to field 
  ! (so we have to add "/prognostic/solver" below)
  if (present(option_path)) then
    option_path_in=option_path
  else
    option_path_in=x%option_path
  end if
    
  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup_petsc_csr(y, b, ksp, &
        name, solver_option_path, lstartfromzero, &
        matrix, option_path=option_path_in)
        
  ! copy array into PETSc vecs
  call petsc_solve_copy_vectors_from_vector_fields(y, b, x, rhs, &
     matrix%row_numbering, lstartfromzero)
    
  ! the solve and convergence check
  call petsc_solve_core(y, matrix%M, b, ksp, matrix%row_numbering, &
          name, solver_option_path, lstartfromzero, &
          literations, vector_x0=x)
        
  ! Copy back the result using the petsc numbering:
  call petsc2field(y, matrix%column_numbering, x)
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy_petsc_csr(y, b, ksp)
    
#else
  FLAbort("Petsc_solve called while not configured with PETSc")
#endif
  
end subroutine petsc_solve_vector_petsc_csr

subroutine petsc_solve_tensor_components(x, matrix, rhs, &
  symmetric, option_path)
  !!< Solve a linear system the nice way. Options for this
  !!< come via the options mechanism. This version solves a linear system
  !!< for each of the components of rhs each time with the same matrix.
  type(tensor_field), intent(inout) :: x
  type(tensor_field), intent(in) :: rhs
  type(csr_matrix), intent(in) :: matrix
  ! if .true. assume rhs is symmetric (so we need to solve for fewer components)
  logical, optional, intent(in):: symmetric
  character(len=*), optional, intent(in) :: option_path

#ifdef HAVE_PETSC
  KSP ksp
  Mat A
  Vec y, b

  type(scalar_field) x_component, rhs_component
  type(petsc_numbering_type) petsc_numbering
  character(len=OPTION_PATH_LEN) solver_option_path, name, option_path_in
  integer literations, i, j, startj
  logical lstartfromzero
  
  assert(x%dim==rhs%dim)
  assert(size(x%val,3)==size(rhs%val,3))
  assert(size(x%val,3)==size(matrix,2))
  assert(size(rhs%val,3)==size(matrix,1))
  
  ! option_path_in may still point to field 
  ! (so we have to add "/prognostic/solver" below)
  if (present(option_path)) then
    option_path_in=option_path
  else
    option_path_in=x%option_path
  end if
  
  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
        name, solver_option_path, lstartfromzero, &
        matrix=matrix, &
        option_path=option_path_in)
 
  ewrite(1,*) 'Solving for multiple components of a tensor field'
  
  startj=1
  do i=1, x%dim
     
     if (present(symmetric)) then
       if (symmetric) then
         ! only computes with rhs(i,j) where j>=i
         startj=i
       end if
     end if
     
     do j=startj, x%dim
       
        ewrite(1, *) 'Now solving for component: ', i, j
       
        x_component=extract_scalar_field(x, i, j)
        rhs_component=extract_scalar_field(rhs, i, j)
        ! copy array into PETSc vecs
        call petsc_solve_copy_vectors_from_scalar_fields(y, b, x_component, matrix, rhs_component, petsc_numbering, lstartfromzero)
         
        ! the solve and convergence check
        call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
              name, solver_option_path, lstartfromzero, &
              literations, x0=x_component%val)
            
        ! Copy back the result using the petsc numbering:
        call petsc2field(y, petsc_numbering, x_component, rhs_component)
     
     end do
  end do
  
  ewrite(1,*) 'Finished solving all components.'
  
  if (present(symmetric)) then
     if (symmetric) then
       
       ewrite(2,*) 'This is a symmetric matrix'
       ewrite(2,*) 'Only components (i,j) with j>=i have been solved for.'
       ewrite(2,*) 'Now copying these to components (j,i).'

       ! copy results x(i,j) of equations with rhs(i,j) where j>=i to x(j,i)
       do i=1, x%dim
          do j=i, x%dim
             x%val(j,i,:)=x%val(i,j,:)
          end do
       end do
     end if
  end if
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy(y, A, b, ksp, petsc_numbering)
  
#else
  FLAbort("Petsc_solve called while not configured with PETSc")
#endif
  
end subroutine petsc_solve_tensor_components

subroutine petsc_solve_csr(x, matrix, rhs, option_path, &
  startfromzero, checkconvergence, iterations, prolongator)
!!< Solve a linear system. This is the lower level version of 
!!< petsc_solve_scalar.
!! Equation to be solved is: matrix * x= rhs
real, dimension(:), intent(inout):: x
type(csr_matrix), intent(in):: matrix
real, dimension(:), intent(in):: rhs
!! path in options tree to solver options:
character(len=*), optional, intent(in):: option_path
!! startfromzero=.true. is default. startfromzero==.false. can be 
logical, optional, intent(in):: startfromzero
!! checkconvergence determines whether to check that the solve finished
!! succesfully (if this call is part of a iterative process this might not
!! be what you want)
logical, optional, intent(in):: checkconvergence
!! returns the number of performed iterations:
integer, optional, intent(out):: iterations
!! prolongator to be used at the first level of 'mg'
type(csr_matrix), optional, intent(in) :: prolongator

  type(block_csr_matrix) block_matrix
  
  ! redefine matrix as a 1 by 1 block matrix !
  block_matrix=wrap(matrix%sparsity, (/ 1, 1 /), matrix%val, &
    name='TemporaryMatrix_petsc_solve_csr')
  
  ! now we simply use the block_csr solver as a special case
  call petsc_solve_block_csr(x, block_matrix, rhs, &
    option_path=option_path, &
    startfromzero=startfromzero, checkconvergence=checkconvergence, &
    iterations=iterations, prolongator=prolongator)
    
  call deallocate(block_matrix)

end subroutine petsc_solve_csr

subroutine petsc_solve_block_csr(x, matrix, rhs, option_path, &
  startfromzero, checkconvergence, iterations, prolongator)
!!< Solves a linear system based on a block matrix.
!! The equation to be solved is: matrix * x= rhs
real, dimension(:), intent(inout):: x
type(block_csr_matrix), intent(in):: matrix
real, dimension(:), intent(in):: rhs
!! path in options tree to solver options:
character(len=*), optional, intent(in):: option_path
!! startfromzero=.true. is default. startfromzero==.false. can be 
logical, optional, intent(in):: startfromzero
!! checkconvergence determines whether to check that the solve finished
!! succesfully (if this call is part of a iterative process this might not
!! be what you want)
logical, optional, intent(in):: checkconvergence
!! returns the number of performed iterations per equation:
integer, optional, intent(out):: iterations
!! prolongator to be used at the first level of 'mg'
type(csr_matrix), optional, intent(in) :: prolongator

#ifdef HAVE_PETSC
  KSP ksp
  Mat A
  Vec y, b

  type(petsc_numbering_type) petsc_numbering
  character(len=OPTION_PATH_LEN) solver_option_path, name
  integer literations
  logical lstartfromzero
  
  assert(.not.matrix%diagonal)
  
  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
        name, solver_option_path, lstartfromzero, &
        block_matrix=matrix, &
        option_path=option_path, &
        startfromzero_in=startfromzero, &
        prolongator=prolongator)
  if (name=='return') return
 
  ! copy array into PETSc vecs
  call petsc_solve_copy_vectors_from_array(y, b, x, rhs, petsc_numbering, lstartfromzero)
     
  ! the solve and convergence check
  call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
        name, solver_option_path, lstartfromzero, &
        literations, x0=x, checkconvergence=checkconvergence)
        
  if (present(iterations)) then
     iterations=literations
  end if  
     
  ! Copy back the result using the petsc numbering:
  call petsc2array(y, petsc_numbering, x)
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy(y, A, b, ksp, petsc_numbering)
  
#else
  FLAbort("Petsc_solve called while not configured with PETSc")
#endif

end subroutine petsc_solve_block_csr

subroutine petsc_solve_csr_multiple(x, matrix, rhs, option_path, &
  startfromzero, checkconvergence, iterations)
!!< Solve a linear system. This is the lower level version of 
!!< petsc_solve_scalar. 
!! Equation to be solved is: matrix * x(:,i)= rhs(:,i) 
real, dimension(:,:), intent(inout):: x
type(csr_matrix), intent(in):: matrix
real, dimension(:,:), intent(in):: rhs
!! path in options tree to solver options:
character(len=*), optional, intent(in):: option_path
!! startfromzero=.true. is default.
logical, optional, intent(in):: startfromzero
!! checkconvergence determines whether to check that the solve finished
!! succesfully (if this call is part of a iterative process this might not
!! be what you want)
logical, optional, intent(in):: checkconvergence
!! returns the number of performed iterations:
integer, optional, dimension(:), intent(out):: iterations

  type(block_csr_matrix) block_matrix
  
  ! redefine matrix as a 1 by 1 block matrix !
  block_matrix=wrap(matrix%sparsity, (/ 1, 1 /), matrix%val, &
    name='TemporaryMatrix_petsc_solve_csr')
  
  ! now we simply use the block_csr solver as a special case
  call petsc_solve_block_csr_multiple(x, block_matrix, rhs, &
    option_path=option_path, &
    startfromzero=startfromzero, checkconvergence=checkconvergence, &
    iterations=iterations)
    
  call deallocate(block_matrix)

end subroutine petsc_solve_csr_multiple
  
subroutine petsc_solve_block_csr_multiple(x, matrix, rhs, option_path, &
  startfromzero, checkconvergence, iterations)
!!< Solves the same linear system based on a block matrix a number of times with different righthand sides.
!!< This is the lower level version
!!< of petsc_solve_scalar/petsc_solve_vector. 
!! i-th equation to be solved is: matrix * x(:,i)= rhs(:,i)
real, dimension(:,:), intent(inout):: x
type(block_csr_matrix), intent(in):: matrix
real, dimension(:,:), intent(in):: rhs
!! path in options tree to solver options:
character(len=*), optional, intent(in):: option_path
!! startfromzero=.true. is default.
logical, optional, intent(in):: startfromzero
!! checkconvergence determines whether to check that the solve finished
!! succesfully (if this call is part of a iterative process this might not
!! be what you want)
logical, optional, intent(in):: checkconvergence
!! returns the number of performed iterations per equation:
integer, optional, intent(out):: iterations(:)

#ifdef HAVE_PETSC
  KSP ksp
  Mat A
  Vec y, b

  type(petsc_numbering_type) petsc_numbering
  character(len=OPTION_PATH_LEN) solver_option_path, name
  integer literations, i
  logical lstartfromzero
  
  assert(.not.matrix%diagonal)
  
  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
        name, solver_option_path, lstartfromzero, &
        block_matrix=matrix, &
        option_path=option_path, &
        startfromzero_in=startfromzero)

  ewrite(1, *) 'Solving multiple equations with the same matrix'
  ewrite(1, *) 'Number of equations: ', size(rhs,2)
  
  do i=1, size(rhs,2)
    
     ! copy array into PETSc vecs
     call petsc_solve_copy_vectors_from_array(y, b, x(:,i), rhs(:,i), petsc_numbering, lstartfromzero)
     
     ewrite(1, *) 'Solving equation no.:', i   

     ! the solve and convergence check
     call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
        name, solver_option_path, lstartfromzero, &
        literations, x0=x(:,i), checkconvergence=checkconvergence)
     if (name=='return') return
        
     if (present(iterations)) then
        iterations(i)=literations
     end if  
     
     ! Copy back the result using the petsc numbering:
     call petsc2array(y, petsc_numbering, x(:,i))
     
  end do

  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy(y, A, b, ksp, petsc_numbering)
  
#else
  FLAbort("Petsc_solve called while not configured with PETSc")
#endif

end subroutine petsc_solve_block_csr_multiple
  
function complete_solver_option_path(option_path)
character(len=*), intent(in):: option_path
character(len=OPTION_PATH_LEN):: complete_solver_option_path

  ! at the moment only prognostic fields have a solver options block
  ! under [field_path]/prognostic/solver/. Other cases should be
  ! implemented here:
  if (have_option(trim(option_path)//'/prognostic/solver')) then
    complete_solver_option_path=trim(option_path)//'/prognostic/solver'
  ! some diagnostic cases have solver blocks now
  else if (have_option(trim(option_path)//'/diagnostic/solver')) then
    complete_solver_option_path=trim(option_path)//'/diagnostic/solver'
  else if (have_option(trim(option_path)//'/solver')) then
    complete_solver_option_path=trim(option_path)//'/solver'
  else
    ewrite(-1,*) 'option_path: ', trim(option_path)
    FLAbort("Missing solver element in provided option_path.")
  end if
  
end function complete_solver_option_path

#ifdef HAVE_PETSC
subroutine petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
  name, solver_option_path, startfromzero, &
  matrix, block_matrix, option_path, startfromzero_in, &
  preconditioner_matrix, prolongator, surface_node_list, &
  exact,internal_smoothing_option)
!!< sets up things needed to call petsc_solve_core
!! Stuff that comes out:
!!
!! PETSc solution vector
Vec, intent(out):: y
!! PETSc matrix
Mat, intent(out):: A
!! PETSc rhs vector
Vec, intent(out):: b
!! Solver object
Mat, intent(out):: ksp
!! numbering from local (i.e. fluidity speak: global) to PETSc (fluidity: universal) numbering
type(petsc_numbering_type), intent(out):: petsc_numbering
!! name of solve, to be printed on log output
character(len=*), intent(out):: name
!! returns the option path to solver/ block for new options, otherwise ""
character(len=*), intent(out):: solver_option_path
!! whether to start with zero initial guess
logical, intent(out):: startfromzero

!! Stuff that goes in:
!! 
!! provide either a matrix or block_matrix to be solved
type(csr_matrix), target, optional, intent(in):: matrix
type(block_csr_matrix), target, optional, intent(in):: block_matrix
!! for new options 
character(len=*), optional, intent(in):: option_path
!! whether to start with zero initial guess (as passed in)
logical, optional, intent(in):: startfromzero_in
!! provide approximation the matrix (only to be used in combination with pctype='KSP')
type(csr_matrix), optional, intent(in) :: preconditioner_matrix
!! prolongator to be used at the first level of 'mg'
type(csr_matrix), optional, intent(in) :: prolongator
!! Stuff needed for internal smoother
integer, dimension(:), optional, intent(in) :: surface_node_list
integer, optional, intent(in) :: internal_smoothing_option
!! exact solution for computing errors
type(scalar_field), optional, intent(in) :: exact
  
  logical, dimension(:), pointer:: inactive_mask
  integer, dimension(:), allocatable:: ghost_nodes
  Mat:: pmat
  MatStructure:: matrix_structure
  ! one of the PETSc supplied orderings see
  ! http://www-unix.mcs.anl.gov/petsc/petsc-as/snapshots/petsc-current/docs/manualpages/MatOrderings/MatGetOrdering.html
  MatOrderingType:: ordering_type
  logical:: use_reordering
  real time1, time2
  integer ierr
  logical:: parallel, timing, have_cache
  type(halo_type), pointer ::  halo
  integer i, j
  KSP, pointer:: ksp_pointer

  if (.not. present(option_path)) then
    FLAbort("No option_path provided to solver call.")
  end if
  solver_option_path=complete_solver_option_path(option_path)

  if (have_option(trim(option_path)//'/name')) then
    call get_option(trim(option_path)//'/name', name)
    ewrite(1,*) 'Inside petsc_solve_(block_)csr, solving for: ', trim(name)
  else
    ewrite(1,*) 'Inside petsc_solve_(block_)csr, solving using option_path: ', trim(option_path)
    name=option_path
  end if

  timing=(debug_level()>=2)
  if (timing) then
    call cpu_time(time1)
  end if
  
  startfromzero=have_option(trim(solver_option_path)//'/start_from_zero')
  if (present_and_true(startfromzero_in) .and. .not. startfromzero) then
    ewrite(2,*) 'Note: startfromzero hard-coded to .true.'
    ewrite(2,*) 'Ignoring setting from solver option.'
    startfromzero=.true.
  end if
  
  ksp=PETSC_NULL_OBJECT
  if (present(matrix)) then
    if (associated(matrix%ksp)) then
      ksp=matrix%ksp
    end if
  else if (present(block_matrix)) then
    if (associated(matrix%ksp)) then
      ksp=block_matrix%ksp
    end if
  end if
  
  if (ksp/=PETSC_NULL_OBJECT) then
    ! oh goody, we've been left something useful!
    call KSPGetOperators(ksp, A, Pmat, matrix_structure, ierr)
    have_cache=.true.
    
    if (have_option(trim(solver_option_path)// &
      '/preconditioner::mg/vertical_lumping/internal_smoother')) then
       ! this option is unsafe with caching, as it needs
       ! DestroyMultigrid to be called on top of PCDestroy to destroy
       ! all its associated objects
       FLExit("Sorry, can't combine internal_smoother with cache_solver_context")
    end if
  else
    ! no cache - we just have to do it all over again
    have_cache=.false.
  end if
  
  ewrite(1, *) 'Assembling matrix.'
  
  call PetscOptionsGetString("", "-ordering_type", ordering_type, &
        use_reordering, ierr)
  if (.not. use_reordering) then
    call get_option(trim(solver_option_path)//'/reordering[0]/name', &
      ordering_type, stat=ierr)
    use_reordering= (ierr==0)
  end if

  if (present(matrix)) then
     ewrite(2, *) 'Number of rows == ', size(matrix, 1)

     ! Create the matrix & vectors.
     
     inactive_mask => get_inactive_mask(matrix)
     ! create list of inactive, ghost_nodes
     if (associated(inactive_mask)) then
        allocate( ghost_nodes(1:count(inactive_mask)) )
        j=0
        do i=1, size(matrix,1)
           if (inactive_mask(i)) then
             j=j+1
             ghost_nodes(j)=i
           end if
        end do
     else
        allocate( ghost_nodes(1:0) )
     end if

     ! set up numbering used in PETSc objects:
     ! NOTE: we use size(matrix,2) here as halo rows may be absent
     call allocate(petsc_numbering, &
        nnodes=size(matrix,2), nfields=1, &
        halo=matrix%sparsity%column_halo, &
        ghost_nodes=ghost_nodes)

     if (use_reordering) then
        call reorder(petsc_numbering, matrix%sparsity, ordering_type)
     end if
     
     if (.not. have_cache) then
       ! create PETSc Mat using this numbering:
       A=csr2petsc(matrix, petsc_numbering, petsc_numbering)
     end if
      
     halo=>matrix%sparsity%column_halo
            
  elseif (present(block_matrix)) then
    
     ewrite(2, *) 'Number of rows == ', size(block_matrix, 1)
     ewrite(2, *) 'Number of blocks == ', blocks(block_matrix,1)
     assert(.not.block_matrix%diagonal)

     ! Create the matrix & vectors.

     ! set up numbering used in PETSc objects:
     call allocate(petsc_numbering, &
        nnodes=block_size(block_matrix,2), nfields=blocks(block_matrix,1), &
        halo=block_matrix%sparsity%column_halo)

      if (use_reordering) then
        call reorder(petsc_numbering, block_matrix%sparsity, ordering_type)
      end if
     
      if (.not. have_cache) then
        ! create PETSc Mat using this numbering:
        A=block_csr2petsc(block_matrix, petsc_numbering, petsc_numbering)
      end if
      
      halo=>block_matrix%sparsity%column_halo
      
  else
  
     ewrite(-1,*) "So what am I going to solve???"
     FLAbort("Wake up!")
      
  end if

  ewrite(1, *) 'Matrix assembly completed.'
  

  if (IsParallel()) then
    parallel= (associated(halo))
  else
    parallel=.false.
  end if
  
  if (have_cache) then
    ! write the cached solver options to log:
    call ewrite_ksp_options(ksp)
  else if (present(preconditioner_matrix)) then
    ewrite(2,*)  'Using provided preconditioner matrix'
    pmat=csr2petsc(preconditioner_matrix, petsc_numbering)
    
    ewrite(2, *) 'Using solver options defined at: ', trim(solver_option_path)
    call SetupKSP(ksp, A, pmat, solver_option_path, parallel, &
      startfromzero_in=startfromzero_in, &
      prolongator=prolongator,surface_node_list=surface_node_list, &
      matrix_csr=matrix,exact=exact, &
      internal_smoothing_option=internal_smoothing_option)
  else
    ewrite(2, *) 'Using solver options defined at: ', trim(solver_option_path)
    call SetupKSP(ksp, A, A, solver_option_path, parallel, &
      startfromzero_in=startfromzero_in, &
      prolongator=prolongator,surface_node_list=surface_node_list, &
      matrix_csr=matrix,exact=exact, &
      internal_smoothing_option=internal_smoothing_option)
  end if
  
  if (.not. have_cache .and. have_option(trim(solver_option_path)// &
    &'/cache_solver_context')) then
    
    ! save the ksp solver context for future generations
    ! (hack with pointer to convince intel compiler that it's
    !  really just the pointed-to value I'm changing)
    if (present(matrix)) then
      ksp_pointer => matrix%ksp
    else if (present(block_matrix)) then
      ksp_pointer => block_matrix%ksp
    end if
    if (associated(ksp_pointer)) then
      ksp_pointer = ksp
      
      ! make sure we don't destroy it, the %ksp becomes a separate reference
      call PetscObjectReference(ksp, ierr)
    else
      ! matrices coming from block() can't cache
      FLAbort("User wants to cache solver context, but no proper matrix is provided.")
    end if
    
  else if (have_cache) then
  
    ! ksp is a copy of matrix%ksp, make it a separate reference, 
    ! so we can KSPDestroy it without destroying matrix%ksp
    call PetscObjectReference(ksp, ierr)
    
    ! same for the matrix, kspgetoperators returns the matrix reference
    ! owned by the ksp - make it a separate reference
    call PetscObjectReference(A, ierr)
    
  end if
  
  b=PetscNumberingCreateVec(petsc_numbering)
  call VecDuplicate(b, y, ierr)
  
  if (timing) then
    call cpu_time(time2)
    ewrite(2,*) "Time spent in Petsc setup: ", time2-time1
  end if
      
end subroutine petsc_solve_setup
  
subroutine petsc_solve_setup_petsc_csr(y, b, ksp, &
  name, solver_option_path, startfromzero, &
  matrix, option_path, startfromzero_in, &
  prolongator,surface_node_list, exact)
!!< sets up things needed to call petsc_solve_core
!! Stuff that comes out:
!!
!! PETSc solution vector
Vec, intent(out):: y
!! PETSc rhs vector
Vec, intent(out):: b
!! Solver object
Mat, intent(out):: ksp
!! name of solve, to be printed on log output
character(len=*), intent(out):: name
!! returns the option path to solver/ block for new options, otherwise ""
character(len=*), intent(out):: solver_option_path
!! whether to start with zero initial guess
logical, intent(out):: startfromzero

!! Stuff that goes in:
!! 
!! provide either a matrix or block_matrix to be solved
type(petsc_csr_matrix), intent(inout):: matrix
!! for new options 
character(len=*), intent(in):: option_path
!! whether to start with zero initial guess (as passed in)
logical, optional, intent(in):: startfromzero_in

!! additional info for "mg" preconditioner:
!! prolongator to be used at the first level of 'mg'
type(csr_matrix), optional, intent(in) :: prolongator
!! Stuff needed for internal smoother
integer, dimension(:), optional, intent(in) :: surface_node_list
!! exact solution for computing errors
type(scalar_field), optional, intent(in) :: exact


  real time1, time2
  integer ierr
  logical parallel, timing

  solver_option_path=complete_solver_option_path(option_path)

  if (have_option(trim(option_path)//'/name')) then
    call get_option(trim(option_path)//'/name', name)
    ewrite(1,*) 'Inside petsc_solve_(block_)csr, solving for: ', trim(name)
  else
    ewrite(1,*) 'Inside petsc_solve_(block_)csr, solving using option_path: ', trim(option_path)
    name=option_path
  end if

  timing=(debug_level()>=2)
  if (timing) then
    call cpu_time(time1)
  end if
  
  call assemble(matrix)
  
  startfromzero=have_option(trim(solver_option_path)//'/start_from_zero')
  if (present_and_true(startfromzero_in) .and. .not. startfromzero) then
    ewrite(2,*) 'Note: startfromzero hard-coded to .true.'
    ewrite(2,*) 'Ignoring setting from solver option.'
    startfromzero=.true.
  end if
  
  if (IsParallel()) then
    parallel= associated(matrix%row_halo)
  else
    parallel= .false.
  end if
  
  ewrite(2, *) 'Using solver options defined at: ', trim(solver_option_path)
  call SetupKSP(ksp, matrix%M, matrix%M, solver_option_path, parallel, &
      startfromzero_in=startfromzero_in, &
      prolongator=prolongator,surface_node_list=surface_node_list, &
      exact=exact)
  
  b=PetscNumberingCreateVec(matrix%column_numbering)
  call VecDuplicate(b, y, ierr)
  
  if (timing) then
    call cpu_time(time2)
    ewrite(2,*) "Time spent in Petsc setup: ", time2-time1
  end if
      
end subroutine petsc_solve_setup_petsc_csr

subroutine petsc_solve_copy_vectors_from_array(y, b,  x, rhs,  petsc_numbering, startfromzero)
Vec, intent(inout):: y, b
real, dimension(:), intent(in):: x, rhs
type(petsc_numbering_type), intent(in):: petsc_numbering
logical, intent(in):: startfromzero

  ewrite(1, *) 'Assembling RHS.'
  
  ! create PETSc vec for rhs using above numbering:
  call array2petsc(rhs, petsc_numbering, b)
  
  ewrite(1, *) 'RHS assembly completed.'

  if (.not. startfromzero) then
    
    ewrite(1, *) 'Assembling initial guess.'

    ! create PETSc vec for initial guess and result using above numbering:
    call array2petsc(x, petsc_numbering, y)

    ewrite(1, *) 'Initial guess assembly completed.'
    
  end if
  
end subroutine petsc_solve_copy_vectors_from_array

subroutine petsc_solve_copy_vectors_from_scalar_fields(y, b,  x, matrix, rhs,  petsc_numbering, startfromzero)
Vec, intent(inout):: y, b
type(scalar_field), target, intent(in):: x, rhs
type(csr_matrix), optional, intent(in):: matrix
type(petsc_numbering_type), intent(in):: petsc_numbering
logical, intent(in):: startfromzero

  type(scalar_field):: ghost_rhs, petsc_solve_rhs, tmp_rhs
  type(mesh_type), pointer:: mesh
  logical, dimension(:), pointer:: inactive_mask
  
  ewrite(1, *) 'Assembling RHS.'  
  
  if (present(matrix)) then
    inactive_mask => get_inactive_mask(matrix)
  else
    nullify(inactive_mask)
  end if
  
  if (associated(inactive_mask)) then
    ! this takes care of the actual lifting of ghost columns, i.e.
    ! row that have column indices referring to ghost nodes, need to 
    ! move their coefficient multiplied with the bc value moved to the rhs
    mesh => rhs%mesh
    call allocate(ghost_rhs, mesh, name="GhostRHS")
    call allocate(petsc_solve_rhs, mesh, name="PetscSolveRHS")
    call allocate(tmp_rhs, mesh, name="TempRHS")
    
    where (inactive_mask)
      ghost_rhs%val=rhs%val
    elsewhere
      ghost_rhs%val=0.
    end where
    
    ! not all processes that see a ghost node may have the right
    ! value set for it. This ensures the owner gets to set the value:
    if(associated(matrix%sparsity%column_halo)) call halo_update(matrix%sparsity%column_halo, ghost_rhs)
    
    ! tmp_rhs is the rhs contribution of lifting the ghost columns
    call mult(tmp_rhs, matrix, ghost_rhs)
    
    call set(petsc_solve_rhs, rhs)
    call addto(petsc_solve_rhs, tmp_rhs, scale=-1.0)
    
    ! note that we don't set the rhs value for the ghost rows
    ! the right value will be substituted after we return from the solve
    call field2petsc(petsc_solve_rhs, petsc_numbering, b)
    
    call deallocate(ghost_rhs)
    call deallocate(petsc_solve_rhs)
    call deallocate(tmp_rhs)
    
  else
  
    ! create PETSc vec for rhs using above numbering:
    call field2petsc(rhs, petsc_numbering, b)
    
  end if
  
  ewrite(1, *) 'RHS assembly completed.'

  if (.not. startfromzero) then
    
    ewrite(1, *) 'Assembling initial guess.'

    ! create PETSc vec for initial guess and result using above numbering:
    call field2petsc(x, petsc_numbering, y)

    ewrite(1, *) 'Initial guess assembly completed.'
    
  end if
  
end subroutine petsc_solve_copy_vectors_from_scalar_fields

subroutine petsc_solve_copy_vectors_from_vector_fields(y, b,  x, rhs,  petsc_numbering, startfromzero)
Vec, intent(inout):: y, b
type(vector_field), intent(in):: x, rhs
type(petsc_numbering_type), intent(in):: petsc_numbering
logical, intent(in):: startfromzero

  ewrite(1, *) 'Assembling RHS.'
  
  ! create PETSc vec for rhs using above numbering:
  call field2petsc(rhs, petsc_numbering, b)
  
  ewrite(1, *) 'RHS assembly completed.'

  if (.not. startfromzero) then
    
    ewrite(1, *) 'Assembling initial guess.'

    ! create PETSc vec for initial guess and result using above numbering:
    call field2petsc(x, petsc_numbering, y)

    ewrite(1, *) 'Initial guess assembly completed.'
    
  end if
  
end subroutine petsc_solve_copy_vectors_from_vector_fields

subroutine petsc_solve_core(y, A, b, ksp, petsc_numbering, &
  name, solver_option_path, startfromzero, &
  iterations, x0, vector_x0, checkconvergence)
!!< inner core of matrix solve, called by all versions of petsc_solve
!! IN: inital guess, OUT: solution
Vec, intent(inout):: y
!! PETSc matrix
Mat, intent(in):: A
!! PETSc vector with right hand side of the equation
Vec, intent(in):: b
!! solver object
KSP, intent(inout):: ksp
!! numbering from local (i.e. fluidity speak: global) to PETSc (fluidity: universal) numbering
type(petsc_numbering_type), intent(in):: petsc_numbering
!! name of solve, to be printed on log output
character(len=*), intent(in):: name
!! for new options option path to solver/ block
character(len=*), intent(in):: solver_option_path
!! whether to start with zero initial guess
logical, intent(in):: startfromzero
!! returns number of performed iterations
integer, intent(out):: iterations
!! initial guess (written out in matrixdump after failed solve)
real, dimension(:), optional, intent(in):: x0
!! initial guess (written out in matrixdump after failed solve) in vector_field form
type(vector_field), optional, intent(in):: vector_x0
!! whether to check convergence (optional legacy argument to be passed straight on)
logical, optional, intent(in):: checkconvergence

  PetscReal norm
  PetscErrorCode ierr
  KSPConvergedReason reason
  PetscLogDouble flops1, flops2
  logical print_norms, timing
  real time1, time2
  
  timing=( debug_level()>=2 )
  
  print_norms=have_option(trim(solver_option_path)//'/diagnostics/print_norms')
  if (print_norms) then
     call VecNorm(b, NORM_2, norm, ierr)
     ewrite(2, *) '2-norm of RHS:', norm
     call VecNorm(b, NORM_INFINITY, norm, ierr)
     ewrite(2, *) 'inf-norm of RHS:', norm
     if (startfromzero) then
        call VecNorm(y, NORM_2, norm, ierr)
        ewrite(2, *) '2-norm of initial guess:', norm
        call VecNorm(y, NORM_INFINITY, norm, ierr)
        ewrite(2, *) 'inf-norm of initial guess:', norm
     end if
     call MatNorm(A, NORM_FROBENIUS, norm, ierr)
     ewrite(2, *) 'Frobenius norm of matrix:', norm
     call MatNorm(A, NORM_INFINITY, norm, ierr)
     ewrite(2, *) 'inf-norm of matrix:', norm
  end if
  
  if (timing) then
    call cpu_time(time1)
    call PetscGetFlops(flops1, ierr)
  end if
  
  ewrite(1, *) 'Entering solver.'

  call KSPSolve(ksp, b, y, ierr)
  call KSPGetConvergedReason(ksp, reason, ierr)
  call KSPGetIterationNumber(ksp, iterations, ierr)

  ewrite(1, *) 'Out of solver.'

  if (timing) then
    call cpu_time(time2)
    call PetscGetFlops(flops2, ierr)
    ewrite(2,*) 'CPU time spent in solver:', time2-time1
    ewrite(2,*) 'MFlops counted by Petsc:', (flops2-flops1)/1e6
    ewrite(2,*) 'MFlops/sec:', (flops2-flops1)/((time2-time1)*1e6)
  end if
  
  ! Check convergence and give warning+matrixdump if needed.
  ! This needs to be done before we copy back the result as
  ! x still contains the initial guess to be used in the matrixdump.
  call ConvergenceCheck(reason, iterations, name, solver_option_path, &
    startfromzero, A, b, petsc_numbering, &
    x0=x0, vector_x0=vector_x0, &
    checkconvergence=checkconvergence)

  ewrite(2, "(A, ' PETSc reason of convergence: ', I0)") trim(name), reason
  ewrite(2, "(A, ' PETSc n/o iterations: ', I0)") trim(name), iterations
    
  if (print_norms) then
     call VecNorm(y, NORM_2, norm, ierr)
     ewrite(2, *) '2-norm of solution:', norm
     call VecNorm(y, NORM_INFINITY, norm, ierr)
     ewrite(2, *) 'inf-norm of solution:', norm
  end if
  
end subroutine petsc_solve_core

subroutine petsc_solve_destroy(y, A, b, ksp, petsc_numbering, exact, &
     error_filename)
Vec, intent(inout):: y
Mat, intent(inout):: A
Vec, intent(inout):: b
KSP, intent(inout):: ksp
type(petsc_numbering_type), intent(inout):: petsc_numbering
type(scalar_field), intent(in), optional :: exact
character(len=*), optional, intent(in) :: error_filename

  PC pc
  PCType pctype
  integer ierr
  
  call VecDestroy(y, ierr)
  call MatDestroy(A, ierr)
  call VecDestroy(b, ierr)
  call KSPGetPC(ksp, pc, ierr)
  call PCGetType(pc, pctype, ierr)
  if (pctype==PCMG) then
    call DestroyMultigrid(pc)
  end if
  call KSPDestroy(ksp, ierr)

  !destroy callback (includes dumping output)
  if(present(exact)) then
     call petsc_solve_callback_destroy(error_filename)
  end if
  ! we could reuse this, but for the moment we don't:
  call deallocate(petsc_numbering)
  
end subroutine petsc_solve_destroy

subroutine petsc_solve_destroy_petsc_csr(y, b, ksp)
Vec, intent(inout):: y
Vec, intent(inout):: b
KSP, intent(inout):: ksp

  PC pc
  PCType pctype
  integer ierr
  
  call VecDestroy(y, ierr)
  call VecDestroy(b, ierr)
  call KSPGetPC(ksp, pc, ierr)
  call PCGetType(pc, pctype, ierr)
  if (pctype==PCMG) then
    call DestroyMultigrid(pc)
  end if
  call KSPDestroy(ksp, ierr)
  
end subroutine petsc_solve_destroy_petsc_csr

subroutine ConvergenceCheck(reason, iterations, name, solver_option_path, &
  startfromzero, A, b, petsc_numbering, x0, vector_x0, checkconvergence)
  !!< Checks reason of convergence. If negative (not converged)
  !!< writes out a scary warning and dumps matrix (if first time), 
  !!< and if reason<0 but reason/=-3
  !!< (i.e. not converged due to other reasons than reaching max_its)
  !!< it sets sig_int to .true. causing the run to halt and dump
  !!< at the end of the time step.
  integer, intent(in):: reason, iterations
  !! name of the thing we're solving for, used in log output:
  character(len=*), intent(in):: name
  !! for new options path to solver options
  character(len=*), intent(in):: solver_option_path  
  ! Arguments needed in the matrixdump:
  logical, intent(in):: startfromzero
  Mat, intent(in):: A
  Vec, intent(in):: b
  type(petsc_numbering_type), intent(in):: petsc_numbering
  ! initial guess to be written in matrixdump (if startfromzero==.false.)
  real, optional, dimension(:), intent(in):: x0
  type(vector_field), optional, intent(in):: vector_x0
  !! if present and .false. do not check, otherwise do check
  logical, optional, intent(in):: checkconvergence

  ! did we dump before? :
  logical, save:: matrixdumped=.false.
  Vec y0
  PetscErrorCode ierr
  character(len=30) reasons(10)
  real spin_up_time, current_time
  
  reasons(1)  = "Undefined"
  reasons(2)  = "KSP_DIVERGED_NULL"
  reasons(3)  = "KSP_DIVERGED_ITS"
  reasons(4)  = "KSP_DIVERGED_DTOL"
  reasons(5)  = "KSP_DIVERGED_BREAKDOWN"
  reasons(6)  = "KSP_DIVERGED_BREAKDOWN_BICG"
  reasons(7)  = "KSP_DIVERGED_NONSYMMETRIC"
  reasons(8)  = "KSP_DIVERGED_INDEFINITE_PC"
  reasons(9)  = "KSP_DIVERGED_NAN"
  reasons(10) = "KSP_DIVERGED_INDEFINITE_MAT"
  
  if (reason<=0) then
     if (present(checkconvergence)) then
        ! checkconvergence==.false. in iterative solver calls that will
        ! not always convergence within the allowed n/o iterations
        if (.not. checkconvergence .and. reason==-3) return
     end if
     ! write reason+iterations to STDERR so we never miss it:
     ewrite(-1,*) 'WARNING: Failed to converge.'
     ewrite(-1,*) "PETSc did not converge for matrix solve of: " // trim(name)
     if((reason>=-10) .and. (reason<=-1)) then
        ewrite(-1,*) 'Reason for non-convergence: ', reasons(-reason)
     else
        ewrite(-1,*) 'Reason for non-convergence is undefined: ', reason
     endif
     ewrite(-1,*) 'Number of iterations: ', iterations
     
     if (have_option(trim(solver_option_path)//'/ignore_all_solver_failures')) then
        ewrite(0,*) 'Specified ignore_all_solver_failures, therefore continuing'
     elseif (reason/=-3 .or. have_option(trim(solver_option_path)//'/never_ignore_solver_failures')) then
        ewrite(-1,*) "Sending signal to dump and finish"
        ! Setting SIGINT in Signal_Vars module will cause dump and crash
        sig_int=.true.
     elseif (have_option(trim(solver_option_path)//'/allow_non_convergence_during_spinup')) then
        call get_option(trim(solver_option_path)//'/allow_non_convergence_during_spin_up/spin_up_time', spin_up_time)
        call get_option('/timestepping/current_time', current_time)
        ewrite(2,*) 'current time:', current_time
         ewrite(2,*) 'spin up time:', spin_up_time
        if (current_time<spin_up_time) then
           ewrite(0,*) "Non-convergence during spin up, therefore continuing"
        else
           ewrite(-1,*) "Non-convergence after spin up period,", &
                        " therefore sending signal to dump and finish"
           sig_int=.true.
        end if
     elseif (have_option(trim(solver_option_path)//'/no_matrixdump')) then
        ! this is a special case internal option for petsc_readnsolve, where we don't want it to
        ! rewrite the matrixdump
        ! pretend we've already dumped before:
        matrixdumped=.true.
     else
        ! legacy case if none if the 3 choices are set:
        ewrite(-1,*) "Sending signal to dump and finish"
        sig_int=.true.
     endif
     
     if (.not. matrixdumped) then
        y0=PetscNumberingCreateVec(petsc_numbering)
        if (startfromzero) then
           call VecZeroEntries(y0, ierr)
        else if (present(x0)) then
           call array2petsc(x0, petsc_numbering, y0)
        else if (present(vector_x0)) then
           call field2petsc(vector_x0, petsc_numbering, y0)
        else
           ewrite(0,*) 'Initial guess not provided in ConvergenceCheck'
           ewrite(0,*) 'This is a bug!!!'
        end if
      
        call DumpMatrixEquation('matrixdump', y0, A, b)
        call VecDestroy(y0, ierr)
        matrixdumped=.true.
     end if
     
  end if
  
end subroutine ConvergenceCheck
#endif

#ifdef HAVE_PETSC
  subroutine SetupKSP(ksp, mat, pmat, solver_option_path, parallel, &       
       startfromzero_in, &
       prolongator,surface_node_list, matrix_csr, exact, &
       internal_smoothing_option)
  !!< Creates the KSP solver context and calls
  !!< setup_ksp_from_options
    KSP, intent(out) :: ksp
    Mat, intent(in) :: mat, pmat
    ! path to solver block (including '/solver')
    character(len=*), intent(in):: solver_option_path
    ! need to know this for setup:
    logical, intent(in):: parallel
    ! if true overrides what is set in the options:
    logical, optional, intent(in):: startfromzero_in
    ! prolongator to be used at the first level of 'mg'
    type(csr_matrix), optional, intent(in) :: prolongator
    ! exact solution for diagnostic purposes:
    type(scalar_field), optional, intent(in) :: exact
    ! additional information for multigrid smoother:
    integer, dimension(:), optional, intent(in) :: surface_node_list
    type(csr_matrix), optional, intent(in) :: matrix_csr
    integer, optional, intent(in) :: internal_smoothing_option
    
    PetscErrorCode ierr
    
    if (parallel) then
       call KSPCreate(MPI_COMM_WORLD, ksp, ierr)
    else
       call KSPCreate(MPI_COMM_SELF, ksp, ierr)
    end if
    call KSPSetOperators(ksp, mat, pmat, DIFFERENT_NONZERO_PATTERN, ierr)
    
    call setup_ksp_from_options(ksp, mat, pmat, solver_option_path, &
      startfromzero_in=startfromzero_in, &
      prolongator=prolongator, &
      surface_node_list=surface_node_list, &
      matrix_csr=matrix_csr, &
      exact=exact, &
      internal_smoothing_option=internal_smoothing_option)
      
  end subroutine SetupKSP
    
  recursive subroutine setup_ksp_from_options(ksp, mat, pmat, solver_option_path, &
      startfromzero_in, &
      prolongator, surface_node_list, matrix_csr, exact, &
      internal_smoothing_option)
  !!< Sets options for the given ksp according to the options
  !!< in the options tree.
    KSP, intent(out) :: ksp
    ! PETSc mat and pmat used to solve
    Mat, intent(in):: mat, pmat
    ! path to solver block (including '/solver')
    character(len=*), intent(in):: solver_option_path
    ! if true overrides what is set in the options:
    logical, optional, intent(in):: startfromzero_in
    ! prolongator to be used at the first level of 'mg'
    type(csr_matrix), optional, intent(in) :: prolongator
    ! exact solution for diagnostic purposes:
    type(scalar_field), optional, intent(in) :: exact
    ! additional information for multigrid smoother:
    integer, dimension(:), optional, intent(in) :: surface_node_list
    type(csr_matrix), optional, intent(in) :: matrix_csr
    integer, optional, intent(in) :: internal_smoothing_option
    
    ! hack to satisfy interface for MatNullSpaceCreate
    ! only works as the array won't actually be used
    PetscObject, dimension(1:0):: PETSC_NULL_OBJECT_ARRAY
    KSPType ksptype
    PC pc
    PetscReal rtol, atol, dtol
    PetscInt max_its
    PetscErrorCode ierr
    MatNullSpace sp ! Nullspace object
    
    logical startfromzero
    
    ewrite(1,*) "Inside setup_ksp_from_options"
    
    call KSPGetPC(ksp, pc, ierr)
    
    ! set ksptype and pctype, and external "PETSC_OPTIONS"
    ! =========================================================
    call get_option(trim(solver_option_path)//'/iterative_method[0]/name', &
      ksptype)
    call setup_pc_from_options(pc, pmat, &
       trim(solver_option_path)//'/preconditioner[0]', &
       prolongator=prolongator, surface_node_list=surface_node_list, &
       matrix_csr=matrix_csr, internal_smoothing_option=internal_smoothing_option)
    
    ! set type first, so ksptype specific PETSC_OPTIONS are picked up
    call KSPSetType(ksp, ksptype, ierr)
    ! set options that may have been supplied via the
    ! PETSC_OPTIONS env. variable
    ! this is done at first, so it will be overwritten by anything
    ! explicitly set via flml options
    call KSPSetFromOptions(ksp, ierr)
    ! set ksptype again to force the flml choice
    call KSPSetType(ksp, ksptype, ierr)
    ewrite(2, *) 'ksp_type:', trim(ksptype)

    ! set max. iterations and tolerances:
    ! =======================================
    call get_option(trim(solver_option_path)//'/relative_error', rtol, &
      default=real(0.0, kind = kind(rtol)))
    call get_option(trim(solver_option_path)//'/absolute_error', atol, &
      default=epsilon(atol))
    ! note that if neither is set the solve will never converge
    ! needs checking

    ! this may end up in the schema:
    dtol=PETSC_DEFAULT_DOUBLE_PRECISION
    ! maximum n/o iterations is required, so no default:
    call get_option(trim(solver_option_path)//'/max_iterations', max_its)
    
    ! set this choice as default (maybe overridden by PETSc options below)
    call KSPSetTolerances(ksp, rtol, atol, dtol, max_its, ierr)
    
    if (have_option(trim(solver_option_path)//'/start_from_zero') &
      .or. present_and_true(startfromzero_in)) then
      call KSPSetInitialGuessNonzero(ksp, PETSC_FALSE, ierr)
      startfromzero=.true.
    else
      call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
      startfromzero=.false.
    end if

    ! If requested, remove nullspace from residual field. This is similar
    ! to imposing a reference pressure node, but initial testing suggests that
    ! it leads to improved rates of convergence.    
    if (have_option(trim(solver_option_path)//'/remove_null_space')) then
       ewrite(2,*) 'Adding null-space removal options to KSP'
       call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL_OBJECT_ARRAY,sp,ierr)
       call KSPSetNullSpace(ksp,sp,ierr)
       call MatNullSpaceDestroy(sp,ierr)
    end if

    ! Inquire about settings as they may have changed by PETSc options:
    call KSPGetTolerances(ksp, rtol, atol, dtol, max_its, ierr)
    
    ewrite(2, *) 'ksp_max_it, ksp_atol, ksp_rtol, ksp_dtol: ', &
      max_its, atol, rtol, dtol
    ewrite(2, *) 'startfromzero:', startfromzero
    
    ! Set up the monitors:
    if (have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/preconditioned_residual')) then
        call KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL_INTEGER, &
           PETSC_NULL_FUNCTION, ierr)
    end if
    if (have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/true_residual')) then
        call KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, PETSC_NULL_INTEGER, &
           PETSC_NULL_FUNCTION, ierr)
    end if
    if (have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/preconditioned_residual_graph')) then
        call KSPMonitorSet(ksp, KSPMonitorLG, PETSC_NULL_INTEGER, &
           PETSC_NULL_FUNCTION, ierr)
    end if

    if(present(exact)) then
       call get_option(trim(solver_option_path)//'/max_iterations', max_its)
       call Petsc_solve_callback_setup(exact,max_its)
       call KSPMonitorSet(ksp,MyKSPMonitor,PETSC_NULL_INTEGER,          &
            &                     PETSC_NULL_FUNCTION,ierr)
    end if

  end subroutine setup_ksp_from_options
    
  recursive subroutine setup_pc_from_options(pc, pmat, option_path, &
    prolongator, surface_node_list, matrix_csr, internal_smoothing_option, &
    is_subpc)
  PC, intent(inout):: pc
  Mat, intent(in):: pmat
  character(len=*), intent(in):: option_path
  ! additional information for multigrid smoother:
  type(csr_matrix), optional, intent(in) :: prolongator
  integer, dimension(:), optional, intent(in) :: surface_node_list
  type(csr_matrix), optional, intent(in) :: matrix_csr
  integer, optional, intent(in) :: internal_smoothing_option
  ! if present and true, don't setup sor and eisenstat as subpc (again)
  logical, optional, intent(in) :: is_subpc
    
    KSP:: subksp
    PC:: subpc
    PCType:: pctype, hypretype
    PetscErrorCode:: ierr
    
    call get_option(trim(option_path)//'/name', pctype)
    if (pctype==PCMG) then
      call SetupMultigrid(pc, pmat, ierr, &
            external_prolongator=prolongator, &
            surface_node_list=surface_node_list, &
            matrix_csr=matrix_csr, &
            internal_smoothing_option=internal_smoothing_option)
      if (ierr/=0) then
         if (IsParallel()) then
           ! we give up as SOR is probably not good enough either
           ! for big paralel problems
           ewrite(-1,*) 'Set up of mg preconditioner failed for:'
           ewrite(-1,*) trim(option_path)
           FLAbort("MG failed: try another preconditioner or improve partioning")
         end if
         ewrite(0,*) 'Set up of mg preconditioner failed for:'
         ewrite(0,*) trim(option_path)
         ewrite(0,*) 'choosing "sor" instead for now'
         pctype=PCSOR
         call PCSetType(pc, pctype, ierr)
      end if
      
      ! set options that may have been supplied via the
      ! PETSC_OPTIONS env. variable for the preconditioner
      call PCSetFromOptions(pc, ierr)

    else if (pctype=='hypre') then
#ifdef HAVE_HYPRE
      call PCSetType(pc, pctype, ierr)
      call get_option(trim(option_path)//'/hypre_type[0]/name', &
        hypretype)
      call PCHYPRESetType(pc, hypretype, ierr)
#else
      ewrite(0,*) 'In solver option:', option_path
      FLAbort("The fluidity binary is built without hypre support!")
#endif
    else if (pctype==PCKSP) then
      
       ! this replaces the preconditioner by a complete solve
       ! (based on the pmat matrix)
       call PCSetType(pc, pctype, ierr)
       
       ! set the options for the ksp of this complete solve
       call PCKSPGetKSP(pc, subksp, ierr)
       ewrite(1,*) "Going into setup_ksp_from_options again to set the options&
          &for the complete ksp solve of the preconditioner"
       call setup_ksp_from_options(subksp, pmat, pmat, &
         trim(option_path)//'/solver')
       ewrite(1,*) "Returned from setup_ksp_from_options for the preconditioner solve&
          &, now setting options for the outer solve"
      
    else if (pctype==PCASM .or. pctype==PCBJACOBI) then
      
      call PCSetType(pc, pctype, ierr)
      ! need to call this before the subpc can be retrieved:
      call PCSetup(pc, ierr)
      
      if (pctype==PCBJACOBI) then
        call PCBJACOBIGetSubKSP(pc, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, subksp, ierr)
      else
        call PCASMGetSubKSP(pc, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, subksp, ierr)
      end if

      call KSPGetPC(subksp, subpc, ierr)
      ! recursively call to setup the subpc
      ewrite(2,*) "Going into setup_pc_from_options for the subpc within the local domain."
      call setup_pc_from_options(subpc, pmat, &
         trim(option_path)//'/preconditioner[0]', &
         prolongator=prolongator, surface_node_list=surface_node_list, &
         matrix_csr=matrix_csr, internal_smoothing_option=internal_smoothing_option, &
         is_subpc=.true.)
      ewrite(2,*) "Finished setting up subpc."      
      
    else if (IsParallel() .and. (pctype==PCSOR .or. &
      pctype==PCEISENSTAT) .and. .not. present_and_true(is_subpc)) then
        
       ! in parallel set sor and eisenstat up in combination with pcbjacobi
       ewrite(2,*) "In parallel sor and eisenstat are setup as bjacobi with&
          & sor/eisenstat as subpc in the local domain."
       call PCSetType(pc, PCBJACOBI, ierr)
       ! need to call this before the subpc can be retrieved:
       call PCSetup(pc, ierr)
       call PCBJACOBIGetSubKSP(pc, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, subksp, ierr)
       call KSPGetPC(subksp, subpc, ierr)
       call PCSetType(subpc, pctype, ierr)
       
    else
       
       ! this doesn't work for hypre
       call PCSetType(pc, pctype, ierr)
       ! set options that may have been supplied via the
       ! PETSC_OPTIONS env. variable for the preconditioner
       call PCSetFromOptions(pc, ierr)
       ! set pctype again to enforce flml choice
       call PCSetType(pc, pctype, ierr)

    end if
    
    ewrite(2, *) 'pc_type: ', trim(pctype)
    if (pctype=='hypre') then
      ewrite(2,*) 'pc_hypre_type:', trim(hypretype)
    end if
    
  end subroutine setup_pc_from_options
    
  subroutine ewrite_ksp_options(ksp)
    KSP, intent(in):: ksp
    
    PC:: pc
    KSPType:: ksptype
    PCType:: pctype
    PetscReal:: rtol, atol, dtol
    PetscInt:: maxits
    PetscTruth:: flag
    PetscErrorCode:: ierr
    
    ewrite(2, *) 'Using solver options from cache:'
    
    call KSPGetType(ksp, ksptype, ierr)
    ewrite(2, *) 'ksp_type: ', trim(ksptype)
    
    call KSPGetPC(ksp, pc, ierr)
    call PCGetType(pc, pctype, ierr)
    ewrite(2, *) 'pc_type: ', trim(pctype)
    
    call KSPGetTolerances(ksp, rtol, atol, dtol, maxits, ierr)
    ewrite(2, *) 'ksp_max_it, ksp_atol, ksp_rtol, ksp_dtol: ', &
      maxits, atol, rtol, dtol
    
    call KSPGetInitialGuessNonzero(ksp, flag, ierr)
    ewrite(2, *) 'startfromzero:', .not. flag
    
  end subroutine ewrite_ksp_options
#endif

subroutine set_solver_options_with_path(field_option_path, &
    ksptype, pctype, atol, rtol, max_its, &
    start_from_zero, petsc_options)
  character(len=*), intent(in):: field_option_path
  character(len=*), optional, intent(in):: ksptype
  character(len=*), optional, intent(in):: pctype
  real, optional, intent(in):: rtol, atol
  integer, optional, intent(in):: max_its
  logical, optional, intent(in):: start_from_zero
  character(len=*), optional, intent(in):: petsc_options
  
  character(len=OPTION_PATH_LEN):: option_path
  integer:: stat
  
  ! set the various options if supplied
  ! otherwise set a sensible default
  if (have_option(trim(field_option_path)//'/solver')) then
    option_path=trim(field_option_path)//'/solver'
  else if (have_option(trim(field_option_path)//'/prognostic/solver')) then
    option_path=trim(field_option_path)//'/prognostic/solver'
  else
    option_path=trim(field_option_path)//'/solver'
    call add_option(option_path, stat=stat)
  end if
  
  if (present(ksptype)) then
     call add_option(trim(option_path)//'/iterative_method::'//trim(ksptype), stat=stat)
  else
     call add_option(trim(option_path)//'/iterative_method::'//trim(KSPGMRES), stat=stat)
  endif
  
  if (present(pctype)) then
     call add_option(trim(option_path)//'/preconditioner::'//trim(pctype), stat=stat)
  else
     call add_option(trim(option_path)//'/preconditioner::'//trim(PCSOR), stat=stat)
  endif
  
  if (present(rtol)) then
     call set_option(trim(option_path)//'/relative_error', rtol, stat=stat)
  else
     call set_option(trim(option_path)//'/relative_error', 1.0e-7, stat=stat)
  end if
  
  if (present(atol)) then
     call set_option(trim(option_path)//'/absolute_error', atol, stat=stat)
  end if
  
  if (present(max_its)) then
     call set_option(trim(option_path)//'/max_iterations', max_its, stat=stat)
  else
     call set_option(trim(option_path)//'/max_iterations', 10000, stat=stat)
  end if
  
  if (present(start_from_zero)) then
     if (start_from_zero) then
        call add_option(trim(option_path)//'/start_from_zero', stat=stat)
     end if
  end if
  
  if (present(petsc_options)) then
     call set_option(trim(option_path)//'/petsc_options', petsc_options, stat=stat)
  end if
  
end subroutine set_solver_options_with_path

subroutine set_solver_options_scalar(field, &
    ksptype, pctype, atol, rtol, max_its, &
    start_from_zero, petsc_options)
  type(scalar_field), intent(inout):: field
  character(len=*), optional, intent(in):: ksptype
  character(len=*), optional, intent(in):: pctype
  real, optional, intent(in):: rtol, atol
  integer, optional, intent(in):: max_its
  logical, optional, intent(in):: start_from_zero
  character(len=*), optional, intent(in):: petsc_options
  
  integer:: stat

  if (field%option_path=="") then
     if (field%name=="") then
        FLAbort("In set_solver_options: if no option_path is supplied a field name is required.")
     end if
     call add_option("/solver_options/", stat=stat)
     field%option_path="/solver_options/scalar_field::"//trim(field%name)
     call add_option(field%option_path, stat=stat)  
  end if
  
  call set_solver_options_with_path(field%option_path, &
      ksptype=ksptype, pctype=pctype, atol=atol, rtol=rtol, max_its=max_its, &
      start_from_zero=start_from_zero, petsc_options=petsc_options)

end subroutine set_solver_options_scalar

subroutine set_solver_options_vector(field, &
    ksptype, pctype, atol, rtol, max_its, &
    start_from_zero, petsc_options)
  type(vector_field), intent(inout):: field
  character(len=*), optional, intent(in):: ksptype
  character(len=*), optional, intent(in):: pctype
  real, optional, intent(in):: rtol, atol
  integer, optional, intent(in):: max_its
  logical, optional, intent(in):: start_from_zero
  character(len=*), optional, intent(in):: petsc_options
  
  integer:: stat
  
  if (field%option_path=="") then
     if (field%name=="") then
        FLAbort("In set_solver_options: if no option_path is supplied a field name is required.")
     end if
     call add_option("/solver_options/", stat=stat)
     field%option_path="/solver_options/vector_field::"//trim(field%name)
     call add_option(field%option_path, stat=stat)
     
  end if
  
  call set_solver_options_with_path(field%option_path, &
      ksptype=ksptype, pctype=pctype, atol=atol, rtol=rtol, max_its=max_its, &
      start_from_zero=start_from_zero, petsc_options=petsc_options)

end subroutine set_solver_options_vector

subroutine set_solver_options_tensor(field, &
    ksptype, pctype, atol, rtol, max_its, &
    start_from_zero, petsc_options)
  type(tensor_field), intent(inout):: field
  character(len=*), optional, intent(in):: ksptype
  character(len=*), optional, intent(in):: pctype
  real, optional, intent(in):: rtol, atol
  integer, optional, intent(in):: max_its
  logical, optional, intent(in):: start_from_zero
  character(len=*), optional, intent(in):: petsc_options
  
  integer:: stat
  
  if (field%option_path=="") then
     if (field%name=="") then
        FLAbort("In set_solver_options: if no option_path is supplied a field name is required.")
     end if
     call add_option("/solver_options/", stat=stat)
     field%option_path="/solver_options/vector_field::"//trim(field%name)
     call add_option(field%option_path, stat=stat)
  end if
  
  call set_solver_options_with_path(field%option_path, &
      ksptype=ksptype, pctype=pctype, atol=atol, rtol=rtol, max_its=max_its, &
      start_from_zero=start_from_zero, petsc_options=petsc_options)

end subroutine set_solver_options_tensor

#ifdef HAVE_PETSC
subroutine petsc_solve_callback_setup(exact,max_its)
  !! sets up call back function to be called each iteration of a solve
  type(scalar_field), intent(in) :: exact
  integer, intent(in) :: max_its
  
  integer :: ierr, i
  integer, allocatable, dimension(:) :: numbering

  call VecCreateSeq(MPI_COMM_SELF, node_count(exact), &
       petsc_monitor_exact, ierr)

  allocate( numbering( node_count(exact) ) )

  do i = 1, size(numbering)
     numbering(i) = i-1
  end do

  call VecDuplicate(petsc_monitor_exact, petsc_monitor_x, ierr)

#ifdef DOUBLEP
  call VecSetValues(petsc_monitor_exact, size(exact%val), &
       numbering, exact%val, INSERT_VALUES, ierr)
#else
  call VecSetValues(petsc_monitor_exact, size(exact%val), &
       numbering, real(exact%val, kind = PetscScalar_kind), INSERT_VALUES, ierr)
#endif
       
  call VecAssemblyBegin(petsc_monitor_exact, ierr)
  call VecAssemblyEnd(petsc_monitor_exact, ierr)
  
  allocate( petsc_monitor_error(max_its+1) )
  allocate( petsc_monitor_flops(max_its+1) )
  petsc_monitor_error = 0.0
  petsc_monitor_flops = 0.0
  petsc_monitor_iteration=0

end subroutine Petsc_solve_callback_setup
  
subroutine petsc_solve_callback_destroy(error_filename)
!! Destroys everything asscociated with the call back function
  character(len=*), optional, intent(in) :: error_filename
  !
  integer :: ierr, i
  integer :: error_unit
  character(len = 100) :: format0

  if(present(error_filename)) then
     !dumping out errors and flops
     error_unit=free_unit()
     open(unit=error_unit, file=trim(error_filename), action="write")    
     format0="(i0,a," &
          & // real_format(padding = 1) // ",a," &
          & // real_format(padding = 1) //")" 
     do i = 1, petsc_monitor_iteration
        write(error_unit, format0) i , " ", &
             & petsc_monitor_error(i), " ", petsc_monitor_flops(i)
     end do
     close(error_unit)
  else
     ewrite(1,*) 'ERROR CALCULATIONS'
     ewrite(1,*) 'iteration, error'
     do i = 1, petsc_monitor_iteration
        ewrite(1,*) i, petsc_monitor_error(i), petsc_monitor_flops(i)
     end do
  end if

  call VecDestroy(petsc_monitor_exact, ierr)
  call VecDestroy(petsc_monitor_x, ierr)
  deallocate( petsc_monitor_error )
  deallocate( petsc_monitor_flops )

end subroutine petsc_solve_callback_destroy
  
subroutine MyKSPMonitor(ksp,n,rnorm,dummy,ierr)
!! The monitor function that gets called each iteration of petsc_solve
!! (if petsc_solve_callback_setup is called)
  PetscInt, intent(in) :: n,dummy
  KSP, intent(in) :: ksp
  PetscErrorCode, intent(out) :: ierr
  
  PetscScalar :: rnorm
  PetscLogDouble :: flops
  Vec:: dummy_vec

  ! Stop warnings.
  ierr = dummy  

  !  Build the solution vector  
  call VecZeroEntries(petsc_monitor_x,ierr)
  ! don't pass PETSC_NULL_OBJECT instead of dummy_vec, as petsc
  ! will clobber it (bug in fortran interface)
  call KSPBuildSolution(ksp,petsc_monitor_x, dummy_vec, ierr)

  ! Compare with exact solution
  call VecAXPY(petsc_monitor_x, real(-1.0, kind = PetscScalar_kind), petsc_monitor_exact, ierr)
  call VecNorm(petsc_monitor_x, NORM_INFINITY, rnorm, ierr)
  call PetscGetFlops(flops,ierr)

  petsc_monitor_error(n+1) = rnorm
  petsc_monitor_iteration = max(petsc_monitor_iteration,n+1)
  petsc_monitor_flops(n+1) = flops
  
  ierr=0
  
end subroutine MyKSPMonitor
#endif

end module solvers
