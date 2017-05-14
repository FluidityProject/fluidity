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
!    amcgsoftware@imperial.ac.uk
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
  use Global_Parameters
  use futils, only: present_and_true, int2str, free_unit, real_format
  use elements
  use spud
  use parallel_tools
#ifdef HAVE_PETSC_MODULES
  use petsc
#endif
  use Sparse_Tools
  use Fields
  use profiler
  use Petsc_tools
  use Signal_Vars
  use Multigrid
  use sparse_tools_petsc
  use sparse_matrices_fields
  use vtk_interfaces
  use halos
  use MeshDiagnostics
  implicit none
  ! Module to provide explicit interfaces to matrix solvers.

#include "petsc_legacy.h"

  ! stuff used in the PETSc monitor (see petsc_solve_callback_setup() below)
  integer :: petsc_monitor_iteration = 0
  Vec :: petsc_monitor_x
  !
  ! if .true. the code will compare with the provided exact answer, and
  ! give the error convergence each iteration:
  logical, save:: petsc_monitor_has_exact=.false.
  ! this requires the following:
  Vec :: petsc_monitor_exact
  real, dimension(:), pointer :: petsc_monitor_error => null()
  PetscLogDouble, dimension(:), pointer :: petsc_monitor_flops => null()
  type(scalar_field), save:: petsc_monitor_exact_sfield
  type(vector_field), save:: petsc_monitor_exact_vfield
  character(len=FIELD_NAME_LEN), save:: petsc_monitor_error_filename=""
  !
  ! if .true. a vtu will be written for each iteration
  logical, save:: petsc_monitor_iteration_vtus=.false.
  ! this requires the following:
  type(petsc_numbering_type), save:: petsc_monitor_numbering
  type(vector_field), target, save:: petsc_monitor_positions
  type(scalar_field), dimension(3), save:: petsc_monitor_sfields
  type(vector_field), dimension(3), save:: petsc_monitor_vfields
  character(len=FIELD_NAME_LEN), save:: petsc_monitor_vtu_name
  integer, save:: petsc_monitor_vtu_series=0
  
private

public petsc_solve, set_solver_options, &
   complete_solver_option_path, petsc_solve_needs_positions

! meant for unit-testing solver code only:
public petsc_solve_core, petsc_solve_destroy, &
  petsc_solve_copy_vectors_from_scalar_fields, &
  setup_ksp_from_options, create_ksp_from_options, petsc_solve_monitor_exact, &
  petsc_solve_monitor_iteration_vtus, attach_null_space_from_options, &
  petsc_solve_setup

interface petsc_solve
   module procedure petsc_solve_scalar, petsc_solve_vector, &
     petsc_solve_scalar_multiple, &
     petsc_solve_vector_components, &
     petsc_solve_tensor_components, &
     petsc_solve_scalar_petsc_csr, petsc_solve_vector_petsc_csr
end interface
  
interface set_solver_options
    module procedure set_solver_options_with_path, &
      set_solver_options_scalar, set_solver_options_vector, set_solver_options_tensor
end interface set_solver_options
  
interface petsc_solve_monitor_exact
    module procedure petsc_solve_monitor_exact_scalar
end interface petsc_solve_monitor_exact

contains

subroutine petsc_solve_scalar(x, matrix, rhs, option_path, &
  preconditioner_matrix, prolongators, surface_node_list, &
  internal_smoothing_option, iterations_taken)
  !!< Solve a linear system the nice way.
  type(scalar_field), intent(inout) :: x
  type(scalar_field), intent(in) :: rhs
  type(csr_matrix), intent(in) :: matrix
  character(len=*), optional, intent(in) :: option_path
  !! 2 experimental arguments to improve preconditioning with extra outside information
  !! provide approximation the matrix (only to be used in combination with pctype='KSP')
  type(csr_matrix), optional, intent(in) :: preconditioner_matrix
  !! prolongators to be used at the first levels of 'mg'
  type(petsc_csr_matrix), dimension(:), optional, intent(in) :: prolongators
  !! surface_node_list for internal smoothing
  integer, dimension(:), optional, intent(in) :: surface_node_list
  !! internal smoothing option
  integer, intent(in), optional :: internal_smoothing_option
  !! the number of petsc iterations taken
  integer, intent(out), optional :: iterations_taken
  
  KSP ksp
  Mat A
  Vec y, b

  character(len=OPTION_PATH_LEN):: solver_option_path
  type(petsc_numbering_type) petsc_numbering
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
  
  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
        solver_option_path, lstartfromzero, &
        matrix=matrix, &
        sfield=x, &
        option_path=option_path, &
        preconditioner_matrix=preconditioner_matrix, &
        prolongators=prolongators, surface_node_list=surface_node_list, &
        internal_smoothing_option=internal_smoothing_option)
 
  ! copy array into PETSc vecs
  call petsc_solve_copy_vectors_from_scalar_fields(y, b, x, &
       & matrix, rhs, petsc_numbering, lstartfromzero)
     
  ! the solve and convergence check
  call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
        solver_option_path, lstartfromzero, literations, &
        sfield=x, x0=x%val)
  
  ! set the optional variable passed out of this procedure 
  ! for the number of petsc iterations taken
  if (present(iterations_taken)) iterations_taken = literations
        
  ! Copy back the result using the petsc numbering:
  call petsc2field(y, petsc_numbering, x, rhs)
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy(y, A, b, ksp, petsc_numbering, &
       & solver_option_path)
  
end subroutine petsc_solve_scalar

subroutine petsc_solve_scalar_multiple(x, matrix, rhs, option_path)
  !!< Solves multiple scalar fields with the same matrix.
  !!< Need to specify an option_path as there's no default
  type(scalar_field), dimension(:), intent(inout) :: x
  type(scalar_field), dimension(:), intent(in) :: rhs
  type(csr_matrix), intent(in) :: matrix
  character(len=*), optional, intent(in) :: option_path

  KSP ksp
  Mat A
  Vec y, b

  type(petsc_numbering_type) petsc_numbering
  character(len=OPTION_PATH_LEN) solver_option_path
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
        solver_option_path, lstartfromzero, &
        matrix=matrix, sfield=x(1), &
        option_path=option_path)
  
  do i=1, size(x)
 
    ! copy array into PETSc vecs
    call petsc_solve_copy_vectors_from_scalar_fields(y, b, &
         x(i), matrix, rhs(i), &
         petsc_numbering, lstartfromzero)
     
    ! the solve and convergence check
    call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
        solver_option_path, lstartfromzero, literations, &
        sfield=x(i), x0=x(i)%val)
        
    ! Copy back the result using the petsc numbering:
    call petsc2field(y, petsc_numbering, x(i), rhs(i))
    
  end do
    
  ewrite(1,*) 'Finished solving all scalar fields'
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy(y, A, b, ksp, petsc_numbering, &
       & solver_option_path)
  
end subroutine petsc_solve_scalar_multiple

subroutine petsc_solve_vector(x, matrix, rhs, option_path, deallocate_matrix)
  !!< Solve a linear system the nice way. Options for this
  !!< come via the options mechanism. 
  type(vector_field), intent(inout) :: x
  type(vector_field), intent(in) :: rhs
  type(block_csr_matrix), intent(inout) :: matrix
  character(len=*), optional, intent(in) :: option_path
  !! deallocate the matrix after it's been copied
  logical, intent(in), optional :: deallocate_matrix

  KSP ksp
  Mat A
  Vec y, b

  type(petsc_numbering_type) petsc_numbering
  character(len=OPTION_PATH_LEN) solver_option_path
  integer literations
  logical lstartfromzero
  
  type(csr_matrix) :: matrixblock
  type(scalar_field) :: rhsblock, xblock
  integer :: i
  
  assert(x%dim==rhs%dim)
  assert(size(x%val(1,:))==size(rhs%val(1,:)))
  assert(size(x%val(1,:))==block_size(matrix,2))
  assert(size(rhs%val(1,:))==block_size(matrix,1))
  assert(x%dim==blocks(matrix,2))
  assert(rhs%dim==blocks(matrix,1))
  
  if(matrix%diagonal) then
    assert(blocks(matrix,1)==blocks(matrix,2))
    ! only want to solve using the diagonal blocks
    do i = 1, blocks(matrix,1)
      matrixblock=block(matrix,i,i)
      rhsblock = extract_scalar_field(rhs, i)
      xblock = extract_scalar_field(x, i)
  
      ! setup PETSc object and petsc_numbering from options and 
      call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
            solver_option_path, lstartfromzero, &
            matrix=matrixblock, &
            vfield=x, &
            option_path=option_path)
    
      ! copy array into PETSc vecs
      call petsc_solve_copy_vectors_from_scalar_fields(y, b, xblock, &
          & matrixblock, rhsblock, petsc_numbering, lstartfromzero)
      
      if(present_and_true(deallocate_matrix).and.(i==blocks(matrix,1))) then
        call deallocate(matrix)
      end if
        
      ! the solve and convergence check
      call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
            solver_option_path, lstartfromzero, literations, &
            vfield=x, x0=xblock%val)
            
      ! Copy back the result using the petsc numbering:
      call petsc2field(y, petsc_numbering, xblock, rhsblock)
      
      ! destroy all PETSc objects and the petsc_numbering
      call petsc_solve_destroy(y, A, b, ksp, petsc_numbering, &
         & solver_option_path)          
    end do
    
  else
  
    ! setup PETSc object and petsc_numbering from options and 
    call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
          solver_option_path, lstartfromzero, &
          block_matrix=matrix, &
          vfield=x, &
          option_path=option_path)
          
    if(present_and_true(deallocate_matrix)) then
      call deallocate(matrix)
    end if
  
    ! copy array into PETSc vecs
    call petsc_solve_copy_vectors_from_vector_fields(y, b, x, rhs, petsc_numbering, lstartfromzero)
      
    ! the solve and convergence check
    call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
            solver_option_path, lstartfromzero, literations, &
            vfield=x, vector_x0=x)
          
    ! Copy back the result using the petsc numbering:
    call petsc2field(y, petsc_numbering, x)
    
    ! destroy all PETSc objects and the petsc_numbering
    call petsc_solve_destroy(y, A, b, ksp, petsc_numbering, &
       & solver_option_path)
  end if
  
end subroutine petsc_solve_vector
  
subroutine petsc_solve_vector_components(x, matrix, rhs, option_path)
  !!< Solve a linear system the nice way. Options for this
  !!< come via the options mechanism. This version solves a linear system
  !!< for each of the components of rhs each time with the same matrix.
  type(vector_field), intent(inout) :: x
  type(vector_field), intent(in) :: rhs
  type(csr_matrix), intent(in) :: matrix
  character(len=*), optional, intent(in) :: option_path

  KSP ksp
  Mat A
  Vec y, b

  type(scalar_field) x_component, rhs_component
  type(petsc_numbering_type) petsc_numbering
  character(len=OPTION_PATH_LEN) solver_option_path, option_path_in
  integer literations, i
  logical lstartfromzero
  
  assert(x%dim==rhs%dim)
  assert(size(x%val(1,:))==size(rhs%val(1,:)))
  assert(size(x%val(1,:))==size(matrix,2))
  assert(size(rhs%val(1,:))==size(matrix,1))
  
  ! option_path_in may still point to field 
  ! (so we have to add "/prognostic/solver" below)
  if (present(option_path)) then
    option_path_in=option_path
  else
    option_path_in=x%option_path
  end if
  
  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
        solver_option_path, lstartfromzero, &
        matrix=matrix, &
        vfield=x, &
        option_path=option_path)
 
  ewrite(1,*) 'Solving for multiple components of a vector field'
  
  do i=1, x%dim
    
    ewrite(1, *) 'Now solving for component: ', i
     x_component=extract_scalar_field(x, i)
     rhs_component=extract_scalar_field(rhs, i)
     ! copy array into PETSc vecs
     call petsc_solve_copy_vectors_from_scalar_fields(y, b, x_component, matrix, rhs_component, petsc_numbering, lstartfromzero)
     
     ! the solve and convergence check
     call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
          solver_option_path, lstartfromzero, literations, &
          vfield=x, x0=x_component%val)
        
     ! Copy back the result using the petsc numbering:
     call petsc2field(y, petsc_numbering, x_component, rhs_component)
     
  end do
  
  ewrite(1,*) 'Finished solving all components.'
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy(y, A, b, ksp, petsc_numbering, &
       & solver_option_path)
  
end subroutine petsc_solve_vector_components

subroutine petsc_solve_scalar_petsc_csr(x, matrix, rhs, option_path, &
  prolongators, surface_node_list)
  !!< Solve a linear system the nice way. Options for this
  !!< come via the options mechanism. 
  type(scalar_field), intent(inout) :: x
  type(scalar_field), intent(in) :: rhs
  type(petsc_csr_matrix), intent(inout) :: matrix
  character(len=*), optional, intent(in) :: option_path
  !! prolongators to be used at the first levels of 'mg'
  type(petsc_csr_matrix), dimension(:), optional, intent(in) :: prolongators
  !! surface_node_list for internal smoothing
  integer, dimension(:), optional, intent(in) :: surface_node_list

  Vec y, b

  character(len=OPTION_PATH_LEN) solver_option_path
  integer literations
  logical lstartfromzero
  
  assert(size(x%val)==size(rhs%val))
  assert(size(x%val)==size(matrix,2))
  assert(size(rhs%val)==size(matrix,1))
  
  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup_petsc_csr(y, b, &
        solver_option_path, lstartfromzero, &
        matrix, &
        sfield=x, &
        option_path=option_path, &
        prolongators=prolongators, surface_node_list=surface_node_list)
        
  ! copy array into PETSc vecs
  call petsc_solve_copy_vectors_from_scalar_fields(y, b, x, rhs=rhs, &
     petsc_numbering=matrix%row_numbering, startfromzero=lstartfromzero)
    
  ! the solve and convergence check
  call petsc_solve_core(y, matrix%M, b, matrix%ksp, matrix%row_numbering, &
          solver_option_path, lstartfromzero, literations, &
          sfield=x, x0=x%val)
        
  ! Copy back the result using the petsc numbering:
  call petsc2field(y, matrix%column_numbering, x)
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy_petsc_csr(y, b, solver_option_path)
  
end subroutine petsc_solve_scalar_petsc_csr

subroutine petsc_solve_vector_petsc_csr(x, matrix, rhs, option_path, &
  prolongators, positions, rotation_matrix)
  !!< Solve a linear system the nice way. Options for this
  !!< come via the options mechanism. 
  type(vector_field), intent(inout) :: x
  type(vector_field), intent(in) :: rhs
  type(petsc_csr_matrix), intent(inout) :: matrix
  character(len=*), optional, intent(in) :: option_path
  !! prolongators to be used at the first levels of 'mg'
  type(petsc_csr_matrix), dimension(:), optional, intent(in) :: prolongators
  !! positions field is only used with remove_null_space/ or multigrid_near_null_space/ with rotational components
  type(vector_field), intent(in), optional :: positions
  !! with rotated bcs: matrix to transform from x,y,z aligned vectors to boundary aligned
  Mat, intent(in), optional:: rotation_matrix


  KSP ksp
  Vec y, b

  character(len=OPTION_PATH_LEN) solver_option_path
  integer literations
  logical lstartfromzero
  
  assert(x%dim==rhs%dim)
  assert(size(x%val(1,:))==size(rhs%val(1,:)))
  assert(size(x%val(1,:))==block_size(matrix,2))
  assert(size(rhs%val(1,:))==block_size(matrix,1))
  assert(x%dim==blocks(matrix,2))
  assert(rhs%dim==blocks(matrix,1))
  
  ! setup PETSc object and petsc_numbering from options and 
  call petsc_solve_setup_petsc_csr(y, b, &
        solver_option_path, lstartfromzero, &
        matrix, vfield=x, option_path=option_path, &
        prolongators=prolongators, &
        positions=positions, rotation_matrix=rotation_matrix)
        
  ! copy array into PETSc vecs
  call petsc_solve_copy_vectors_from_vector_fields(y, b, x, rhs, &
     matrix%row_numbering, lstartfromzero)
    
  ! the solve and convergence check
  call petsc_solve_core(y, matrix%M, b, matrix%ksp, matrix%row_numbering, &
          solver_option_path, lstartfromzero, literations, &
          vfield=x, vector_x0=x)
        
  ! Copy back the result using the petsc numbering:
  call petsc2field(y, matrix%column_numbering, x)
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy_petsc_csr(y, b, solver_option_path)
  
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

  KSP ksp
  Mat A
  Vec y, b

  type(scalar_field) x_component, rhs_component
  type(petsc_numbering_type) petsc_numbering
  character(len=OPTION_PATH_LEN) solver_option_path, option_path_in
  integer literations, i, j, startj
  logical lstartfromzero
  
  assert(all(x%dim==rhs%dim))
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
        solver_option_path, lstartfromzero, &
        matrix=matrix, &
        tfield=x, &
        option_path=option_path_in)
 
  ewrite(1,*) 'Solving for multiple components of a tensor field'
  
  startj=1
  do i=1, x%dim(1)
     
     if (present(symmetric)) then
       if (symmetric) then
         ! only computes with rhs(i,j) where j>=i
         startj=i
       end if
     end if
     
     do j=startj, x%dim(2)
       
        ewrite(1, *) 'Now solving for component: ', i, j
       
        x_component=extract_scalar_field(x, i, j)
        rhs_component=extract_scalar_field(rhs, i, j)
        ! copy array into PETSc vecs
        call petsc_solve_copy_vectors_from_scalar_fields(y, b, x_component, matrix, rhs_component, petsc_numbering, lstartfromzero)
         
        ! the solve and convergence check
        call petsc_solve_core(y, A, b, ksp, petsc_numbering, &
              solver_option_path, lstartfromzero, literations, &
              tfield=x, x0=x_component%val)
            
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
       do i=1, x%dim(1)
          do j=i, x%dim(2)
             x%val(j,i,:)=x%val(i,j,:)
          end do
       end do
     end if
  end if
  
  ! destroy all PETSc objects and the petsc_numbering
  call petsc_solve_destroy(y, A, b, ksp, petsc_numbering, solver_option_path)
  
end subroutine petsc_solve_tensor_components
  
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

subroutine petsc_solve_setup(y, A, b, ksp, petsc_numbering, &
  solver_option_path, startfromzero, &
  matrix, block_matrix, sfield, vfield, tfield, &
  option_path, startfromzero_in, &
  preconditioner_matrix, prolongators, surface_node_list, &
  internal_smoothing_option, positions)
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
!! returns the option path to solver/ block for new options, otherwise ""
character(len=*), intent(out):: solver_option_path
!! whether to start with zero initial guess
logical, intent(out):: startfromzero

!! Stuff that goes in:
!! 
!! provide either a matrix or block_matrix to be solved
type(csr_matrix), target, optional, intent(in):: matrix
type(block_csr_matrix), target, optional, intent(in):: block_matrix
!! provide either a scalar field or vector field to be solved for
type(scalar_field), optional, intent(in):: sfield
type(vector_field), optional, intent(in):: vfield
type(tensor_field), optional, intent(in):: tfield
  
!! if provided overrides sfield%option_path
character(len=*), optional, intent(in):: option_path
!! whether to start with zero initial guess (as passed in)
logical, optional, intent(in):: startfromzero_in
!! provide approximation the matrix (only to be used in combination with pctype='KSP')
type(csr_matrix), optional, intent(in) :: preconditioner_matrix
!! prolongators to be used at the first level of 'mg'
type(petsc_csr_matrix), dimension(:), optional, intent(in) :: prolongators
!! Stuff needed for internal smoother
integer, dimension(:), optional, intent(in) :: surface_node_list
integer, optional, intent(in) :: internal_smoothing_option
!! positions field is only used with remove_null_space/ or multigrid_near_null_space/ with rotational components
type(vector_field), intent(in), optional :: positions

  logical, dimension(:), pointer:: inactive_mask
  integer, dimension(:), allocatable:: ghost_nodes
  Mat:: pmat
  ! one of the PETSc supplied orderings see
  ! http://www-unix.mcs.anl.gov/petsc/petsc-as/snapshots/petsc-current/docs/manualpages/MatOrderings/MatGetOrdering.html
  MatOrderingType:: ordering_type
  logical:: use_reordering
  real time1, time2
  integer ierr
  logical:: parallel, timing, have_cache
  type(halo_type), pointer ::  halo
  integer i, j
  character(len=FIELD_NAME_LEN) :: name
  KSP, pointer:: ksp_pointer

  ! Initialise profiler
  if(present(sfield)) then
     call profiler_tic(sfield, "petsc_setup")
     name = sfield%name
  else if(present(vfield)) then
     call profiler_tic(vfield, "petsc_setup")
     name = vfield%name
  else if(present(tfield)) then
     call profiler_tic(tfield, "petsc_setup")
     name = tfield%name
  else
     FLAbort("petsc_solve_setup should be called with sfield, vfield or tfield")
  end if

  timing=(debug_level()>=2)
  if (timing) then
    call cpu_time(time1)
  end if
  
  
  if (present(option_path)) then
    solver_option_path=complete_solver_option_path(option_path)
  else if (present(sfield)) then
    solver_option_path=complete_solver_option_path(sfield%option_path)
  else if (present(vfield)) then
    solver_option_path=complete_solver_option_path(vfield%option_path)
  else if (present(tfield)) then
    solver_option_path=complete_solver_option_path(tfield%option_path)
  else
    FLAbort("Need to provide either sfield, vfield or tfield to petsc_solve_setup.")
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
    if (associated(block_matrix%ksp)) then
      ksp=block_matrix%ksp
    end if
  end if
  
  if (ksp/=PETSC_NULL_OBJECT) then
    ! oh goody, we've been left something useful!
    call KSPGetOperators(ksp, A, Pmat, ierr)
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
 
  ! Note the explicitly-described options rcm, 1wd and natural are now not
  ! listed explicitly in the schema (but can still be used by adding the
  ! appropriate string in the solver reordering node).
  call PetscOptionsGetString(PETSC_NULL_OBJECT, "", "-ordering_type", ordering_type, use_reordering, ierr)
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
  else

    if (present(preconditioner_matrix)) then
      ewrite(2,*)  'Using provided preconditioner matrix'
      pmat=csr2petsc(preconditioner_matrix, petsc_numbering)
    else
      pmat=A
    end if

    ewrite(2, *) 'Using solver options defined at: ', trim(solver_option_path)
    call attach_null_space_from_options(A, solver_option_path, pmat=pmat, &
      positions=positions, petsc_numbering=petsc_numbering)

    call create_ksp_from_options(ksp, A, pmat, solver_option_path, parallel, &
      petsc_numbering, &
      startfromzero_in=startfromzero_in, &
      prolongators=prolongators, surface_node_list=surface_node_list, &
      matrix_csr=matrix, &
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
    ewrite(2,*) trim(name)// " CPU time spent in PETSc setup: ", time2-time1
  end if

  if(present(sfield)) then
     call profiler_toc(sfield, "petsc_setup")
  else if(present(vfield)) then
     call profiler_toc(vfield, "petsc_setup")
  else if(present(tfield)) then
     call profiler_toc(tfield, "petsc_setup")
  end if
  
end subroutine petsc_solve_setup
  
subroutine petsc_solve_setup_petsc_csr(y, b, &
  solver_option_path, startfromzero, &
  matrix, sfield, vfield, tfield, &
  option_path, startfromzero_in, &
  prolongators,surface_node_list, &
  positions, rotation_matrix)
!!< sets up things needed to call petsc_solve_core
!! Stuff that comes out:
!!
!! PETSc solution vector
Vec, intent(out):: y
!! PETSc rhs vector
Vec, intent(out):: b
!! returns the option path to solver/ block for new options, otherwise ""
character(len=*), intent(out):: solver_option_path
!! whether to start with zero initial guess
logical, intent(out):: startfromzero

!! Stuff that goes in:
!! 
!! provide either a matrix or block_matrix to be solved
type(petsc_csr_matrix), intent(inout):: matrix
!! provide either a scalar field or vector field to be solved for
type(scalar_field), optional, intent(in):: sfield
type(vector_field), optional, intent(in):: vfield
type(tensor_field), optional, intent(in):: tfield

!! overrides sfield%option_path or vfield%option_path
character(len=*), intent(in), optional:: option_path
!! whether to start with zero initial guess (as passed in)
logical, optional, intent(in):: startfromzero_in

!! additional info for "mg" preconditioner:
!! prolongators to be used at the first levels of 'mg'
type(petsc_csr_matrix), dimension(:), optional, intent(in) :: prolongators
!! Stuff needed for internal smoother
integer, dimension(:), optional, intent(in) :: surface_node_list
!! positions field is only used with remove_null_space/ or multigrid_near_null_space/ with rotational components
type(vector_field), intent(in), optional :: positions
!! with rotated bcs: matrix to transform from x,y,z aligned vectors to boundary aligned
Mat, intent(in), optional:: rotation_matrix

  real time1, time2
  integer ierr
  logical parallel, timing
  character(len=FIELD_NAME_LEN) :: name

  ! Initialise profiler
  if(present(sfield)) then
     call profiler_tic(sfield, "petsc_setup")
     name = sfield%name
  else if(present(vfield)) then
     call profiler_tic(vfield, "petsc_setup")
     name = vfield%name
  else if(present(tfield)) then
     call profiler_tic(tfield, "petsc_setup")
     name = tfield%name
  else
     FLAbort("petsc_solve_setup should be called with sfield, vfield or tfield")
  end if

  timing=(debug_level()>=2)
  if (timing) then
    call cpu_time(time1)
  end if

  if (present(option_path)) then
    solver_option_path=complete_solver_option_path(option_path)
  else if (present(sfield)) then
    solver_option_path=complete_solver_option_path(sfield%option_path)
  else if (present(vfield)) then
    solver_option_path=complete_solver_option_path(vfield%option_path)
  else if (present(tfield)) then
    solver_option_path=complete_solver_option_path(tfield%option_path)
  else
    FLAbort("Need to provide either sfield, vfield or tfield to petsc_solve_setup.")
  end if
  
  call assemble(matrix)

  call attach_null_space_from_options(matrix%M, solver_option_path, &
    positions=positions, rotation_matrix=rotation_matrix, petsc_numbering=matrix%column_numbering)
  
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
  if (matrix%ksp==PETSC_NULL_OBJECT) then
  
    call create_ksp_from_options(matrix%ksp, matrix%M, matrix%M, solver_option_path, parallel, &
        matrix%column_numbering, &
        startfromzero_in=startfromzero_in, &
       prolongators=prolongators, surface_node_list=surface_node_list)
  else
    ewrite(2, *) "Reusing ksp from a previous solve"
    call setup_ksp_from_options(matrix%ksp, matrix%M, matrix%M, solver_option_path, &
        matrix%column_numbering, &
        startfromzero_in=startfromzero_in, &
       prolongators=prolongators, surface_node_list=surface_node_list)
  end if
  
  b=PetscNumberingCreateVec(matrix%column_numbering)
  call VecDuplicate(b, y, ierr)

  if (timing) then
    call cpu_time(time2)
    ewrite(2,*) trim(name)// " CPU time spent in PETSc setup: ", time2-time1
  end if

  if(present(sfield)) then
     call profiler_toc(sfield, "petsc_setup")
  else if(present(vfield)) then
     call profiler_toc(vfield, "petsc_setup")
  else if(present(tfield)) then
     call profiler_toc(tfield, "petsc_setup")
  end if
      
end subroutine petsc_solve_setup_petsc_csr

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
  
  call profiler_tic(x, "field2petsc")
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
  call profiler_toc(x, "field2petsc")
  
end subroutine petsc_solve_copy_vectors_from_scalar_fields

subroutine petsc_solve_copy_vectors_from_vector_fields(y, b,  x, rhs,  petsc_numbering, startfromzero)
Vec, intent(inout):: y, b
type(vector_field), intent(in):: x, rhs
type(petsc_numbering_type), intent(in):: petsc_numbering
logical, intent(in):: startfromzero

  call profiler_tic(x, "field2petsc")
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
  call profiler_toc(x, "field2petsc")
  
end subroutine petsc_solve_copy_vectors_from_vector_fields

subroutine petsc_solve_core(y, A, b, ksp, petsc_numbering, &
  solver_option_path, startfromzero, &
  iterations, &
  sfield, vfield, tfield, &
  x0, vector_x0, checkconvergence, nomatrixdump)
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
!! for new options option path to solver/ block
character(len=*), intent(in):: solver_option_path
!! whether to start with zero initial guess
logical, intent(in):: startfromzero
!! returns number of performed iterations
integer, intent(out):: iterations
!! provide either a scalar field or vector field to be solved for
!! (only for logging/timing purposes)
type(scalar_field), optional, intent(in):: sfield
type(vector_field), optional, intent(in):: vfield
type(tensor_field), optional, intent(in):: tfield
!! initial guess (written out in matrixdump after failed solve)
real, dimension(:), optional, intent(in):: x0
!! initial guess (written out in matrixdump after failed solve) in vector_field form
type(vector_field), optional, intent(in):: vector_x0
!! whether to check convergence (optional legacy argument to be passed straight on)
logical, optional, intent(in):: checkconvergence
!! logical to prevent dump of matrix equation (full projection solve eqn cannot be dumped):
logical, optional, intent(in):: nomatrixdump

  PetscReal norm
  PetscErrorCode ierr
  KSPConvergedReason reason
  MatNullSpace nullsp
  PetscLogDouble flops1, flops2
  Mat mat, pmat
  character(len=FIELD_NAME_LEN):: name
  logical print_norms, timing
  real time1, time2

  ! Initialise profiler
  if(present(sfield)) then
    name=sfield%name
    call profiler_tic(sfield, "solve")
  else if(present(vfield)) then
    name=vfield%name
    call profiler_tic(vfield, "solve")
  else if(present(tfield)) then
    name=tfield%name
    call profiler_tic(tfield, "solve")
  else
    FLAbort("petsc_solve_core should be called with sfield, vfield or tfield")
  end if
  
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

  ! if a null space is defined for the petsc matrix, make sure it's projected out of the rhs
  call KSPGetOperators(ksp, mat, pmat, ierr)
  call MatGetNullSpace(mat, nullsp, ierr)
  if (ierr==0  .and. nullsp/=PETSC_NULL_OBJECT) then
    ewrite(2,*) "Projecting nullspace from RHS"
    call MatNullSpaceRemove(nullsp, b, PETSC_NULL_OBJECT, ierr)
  end if

  call KSPSolve(ksp, b, y, ierr)
  call KSPGetConvergedReason(ksp, reason, ierr)
  call KSPGetIterationNumber(ksp, iterations, ierr)

  ewrite(1, *) 'Out of solver.'

  if (timing) then
    call cpu_time(time2)
    call PetscGetFlops(flops2, ierr)
    ewrite(2,*) trim(name)// ' CPU time spent in solver: ',time2-time1
    ewrite(2,*) trim(name)// ' MFlops counted by Petsc: ',(flops2-flops1)/1e6
    ewrite(2,*) trim(name)// ' MFlops/sec: ',(flops2-flops1)/((time2-time1)*1e6)
  end if
  
  if(have_option(trim(solver_option_path)//'/diagnostics/dump_matrix')) then
    if(present_and_true(nomatrixdump)) then
      ewrite(0,*) 'Requested to dump matrix on solve that is hard coded not to'
      ewrite(0,*) 'dump matrices.  Therefore ignoring the dump_matrix option request.'
    else
      call dump_matrix_option(solver_option_path, startfromzero, A, b, &
                              petsc_numbering, &
                              x0=x0, vector_x0=vector_x0)
    end if
  end if
  
  ! Check convergence and give warning+matrixdump if needed.
  ! This needs to be done before we copy back the result as
  ! x still contains the initial guess to be used in the matrixdump.
  call ConvergenceCheck(reason, iterations, name, solver_option_path, &
       startfromzero, A, b, petsc_numbering, &
       x0=x0, vector_x0=vector_x0, &
       checkconvergence=checkconvergence,nomatrixdump=nomatrixdump)

  ewrite(2, "(A, ' PETSc reason of convergence: ', I0)") trim(name), reason
  ewrite(2, "(A, ' PETSc n/o iterations: ', I0)") trim(name), iterations
    
  if (print_norms) then
     call VecNorm(y, NORM_2, norm, ierr)
     ewrite(2, *) '2-norm of solution:', norm
     call VecNorm(y, NORM_INFINITY, norm, ierr)
     ewrite(2, *) 'inf-norm of solution:', norm
  end if

  if(present(sfield)) then
    call profiler_toc(sfield, "solve")
  else if(present(vfield)) then
    call profiler_toc(vfield, "solve")
  else if(present(tfield)) then
    call profiler_toc(tfield, "solve")
  end if
  
end subroutine petsc_solve_core

subroutine petsc_solve_destroy(y, A, b, ksp, petsc_numbering, &
  solver_option_path)
Vec, intent(inout):: y
Mat, intent(inout):: A
Vec, intent(inout):: b
KSP, intent(inout):: ksp
type(petsc_numbering_type), intent(inout):: petsc_numbering
character(len=*), intent(in):: solver_option_path

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

  ! destroy everything associated with the monitors
  if(have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/true_error') .or. &
       have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/iteration_vtus')) then
     ! note we have to check the option itself and not the logicals
     ! as we may be in an inner solver, where only the outer solve 
     ! has the monitor set
     call petsc_monitor_destroy()
     petsc_monitor_has_exact=.false.
     petsc_monitor_iteration_vtus=.false.
  end if
  ! we could reuse this, but for the moment we don't:
  call deallocate(petsc_numbering)
  
end subroutine petsc_solve_destroy

subroutine petsc_solve_destroy_petsc_csr(y, b, solver_option_path)
Vec, intent(inout):: y
Vec, intent(inout):: b
character(len=*), intent(in):: solver_option_path

  integer ierr
  
  call VecDestroy(y, ierr)
  call VecDestroy(b, ierr)
  
  ! destroy everything associated with the monitors
  if(have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/true_error') .or. &
       have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/iteration_vtus')) then
     ! note we have to check the option itself and not the logicals
     ! as we may be in an inner solver, where only the outer solve 
     ! has the monitor set
     call petsc_monitor_destroy()
     petsc_monitor_has_exact=.false.
     petsc_monitor_iteration_vtus=.false.
  end if
  
end subroutine petsc_solve_destroy_petsc_csr

subroutine ConvergenceCheck(reason, iterations, name, solver_option_path, &
  startfromzero, A, b, petsc_numbering, x0, vector_x0, checkconvergence, nomatrixdump)
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
  !! if present do not dump matrix equation:
  logical, optional, intent(in):: nomatrixdump

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
     if(present_and_true(nomatrixdump)) matrixdumped = .true.    
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
  
subroutine dump_matrix_option(solver_option_path, startfromzero, A, b, &
                              petsc_numbering, &
                              x0, vector_x0)
                              
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
  
  character(len=FIELD_NAME_LEN) :: filename
  PetscErrorCode ierr
  Vec y0
  !! integer index for the optional dumping of the matrix
  integer, save :: dump_matrix_index=0

  call get_option(trim(solver_option_path)//'/diagnostics/dump_matrix/filename', filename)
  
  dump_matrix_index = dump_matrix_index + 1
  
  y0=PetscNumberingCreateVec(petsc_numbering)
  if (startfromzero) then
     call VecZeroEntries(y0, ierr)
  else if (present(x0)) then
     call array2petsc(x0, petsc_numbering, y0)
  else if (present(vector_x0)) then
     call field2petsc(vector_x0, petsc_numbering, y0)
  else
     ewrite(0,*) 'Initial guess not provided in dump_matrix_option'
     ewrite(0,*) 'This is a bug!!!'
  end if

  call DumpMatrixEquation(trim(filename)//"_"//int2str(dump_matrix_index), y0, A, b)
  call VecDestroy(y0, ierr)
 
end subroutine dump_matrix_option

subroutine create_ksp_from_options(ksp, mat, pmat, solver_option_path, parallel, &
       petsc_numbering, &
       startfromzero_in, &
       prolongators, surface_node_list, matrix_csr, &
       internal_smoothing_option)
  !!< Creates the KSP solver context and calls
  !!< setup_ksp_from_options
    KSP, intent(out) :: ksp
    Mat, intent(in) :: mat, pmat
    ! path to solver block (including '/solver')
    character(len=*), intent(in):: solver_option_path
    ! need to know this for setup:
    logical, intent(in):: parallel
    ! used in the monitors:
    type(petsc_numbering_type), intent(in):: petsc_numbering
    ! if true overrides what is set in the options:
    logical, optional, intent(in):: startfromzero_in
    ! prolongators to be used at the first levels of 'mg'
    type(petsc_csr_matrix), dimension(:), optional, intent(in) :: prolongators
    ! additional information for multigrid smoother:
    integer, dimension(:), optional, intent(in) :: surface_node_list
    type(csr_matrix), optional, intent(in) :: matrix_csr
    integer, optional, intent(in) :: internal_smoothing_option
    
    PetscErrorCode ierr
    
    if (parallel) then
       call KSPCreate(MPI_COMM_FEMTOOLS, ksp, ierr)
    else
       call KSPCreate(MPI_COMM_SELF, ksp, ierr)
    end if
    call KSPSetOperators(ksp, mat, pmat, ierr)
    
    call setup_ksp_from_options(ksp, mat, pmat, solver_option_path, &
      petsc_numbering=petsc_numbering, &
      startfromzero_in=startfromzero_in, &
      prolongators=prolongators, &
      surface_node_list=surface_node_list, &
      matrix_csr=matrix_csr, &
      internal_smoothing_option=internal_smoothing_option)
      
  end subroutine create_ksp_from_options
    
  recursive subroutine setup_ksp_from_options(ksp, mat, pmat, solver_option_path, &
      petsc_numbering, startfromzero_in, prolongators, surface_node_list, matrix_csr, &
      internal_smoothing_option)
  !!< Sets options for the given ksp according to the options
  !!< in the options tree.
    KSP, intent(out) :: ksp
    ! PETSc mat and pmat used to solve
    Mat, intent(in):: mat, pmat
    ! path to solver block (including '/solver')
    character(len=*), intent(in):: solver_option_path
    ! used in the monitors and fieldsplit
    type(petsc_numbering_type), optional, intent(in):: petsc_numbering
    ! if true overrides what is set in the options:
    logical, optional, intent(in):: startfromzero_in
    ! prolongators to be used at the first levels of 'mg'
    type(petsc_csr_matrix), dimension(:), optional, intent(in) :: prolongators
    ! additional information for multigrid smoother:
    integer, dimension(:), optional, intent(in) :: surface_node_list
    type(csr_matrix), optional, intent(in) :: matrix_csr
    integer, optional, intent(in) :: internal_smoothing_option
    
#if PETSC_VERSION_MINOR<6 || (PETSC_VERSION_MINOR==6 && PETSC_VERSION_SUBMINOR==0)
    MatNullSpace nullsp
#endif
    KSPType ksptype
    PC pc
    PetscReal rtol, atol, dtol
    PetscInt max_its, lrestart
    PetscErrorCode ierr
    PetscObject vf
    
    logical startfromzero, remove_null_space
    
    ewrite(1,*) "Inside setup_ksp_from_options"
    
    
    ! first set pc options
    ! =========================================================
    call KSPGetPC(ksp, pc, ierr)
    call setup_pc_from_options(pc, pmat, &
       trim(solver_option_path)//'/preconditioner[0]', &
       petsc_numbering=petsc_numbering, &
       prolongators=prolongators, surface_node_list=surface_node_list, &
       matrix_csr=matrix_csr, internal_smoothing_option=internal_smoothing_option)
    
    ! then ksp type
    ! =========================================================
    call get_option(trim(solver_option_path)//'/iterative_method[0]/name', &
      ksptype)
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

    if(trim(ksptype) == 'gmres') then
       call get_option(trim(solver_option_path)//&
            '/iterative_method::gmres/restart', lrestart, default=-1)
       if (lrestart >= 0) then
          call KSPGMRESSetRestart(ksp, lrestart, ierr)
          ewrite(2, *) 'restart:', lrestart
       end if
    end if

    ! set max. iterations and tolerances:
    ! =======================================
    call get_option(trim(solver_option_path)//'/relative_error', rtol)
    call get_option(trim(solver_option_path)//'/absolute_error', atol, &
      default=real(0.0, kind = kind(rtol)))
    ! note that if neither is set the solve will never converge
    ! needs checking

    ! this may end up in the schema:
    dtol=PETSC_DEFAULT_REAL
    ! maximum n/o iterations is required, so no default:
    call get_option(trim(solver_option_path)//'/max_iterations', max_its)
    
    ! set this choice as default (maybe overridden by PETSc options below)
    call KSPSetTolerances(ksp, rtol, atol, dtol, max_its, ierr)
    
    if (have_option(trim(solver_option_path)//'/start_from_zero') &
      .or. present_and_true(startfromzero_in) .or. ksptype==KSPPREONLY) then
      call KSPSetInitialGuessNonzero(ksp, PETSC_FALSE, ierr)
      startfromzero=.true.
    else
      call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
      startfromzero=.false.
    end if


    ! Inquire about settings as they may have changed by PETSc options:
    call KSPGetTolerances(ksp, rtol, atol, dtol, max_its, ierr)
    
    ewrite(2, *) 'ksp_max_it, ksp_atol, ksp_rtol, ksp_dtol: ', &
      max_its, atol, rtol, dtol
    ewrite(2, *) 'startfromzero:', startfromzero
    
    ! cancel all existing monitors (if reusing the same ksp)
    call KSPMonitorCancel(ksp, ierr)

    ! Set up the monitors:
    if (have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/preconditioned_residual')) then
        call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, &
           PETSC_VIEWER_DEFAULT,vf,ierr)
        call KSPMonitorSet(ksp, KSPMonitorDefault, vf, &
           PetscViewerAndFormatDestroy, ierr)
    end if
    if (have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/true_residual')) then
        call PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, &
           PETSC_VIEWER_DEFAULT,vf,ierr)
        call KSPMonitorSet(ksp, KSPMonitorTrueResidualNorm, vf, &
           PetscViewerAndFormatDestroy, ierr)
    end if

    if (have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/true_error') &
       .and. .not. petsc_monitor_has_exact) then
       ewrite(-1,*) "Solver option diagnostics/monitors/true_error set but "
       ewrite(0,*) "petsc_solve_monitor_exact() not called. This probably means"
       ewrite(0,*) "this version of petsc_solve() doesn't support the monitor."
       FLExit("petsc_solve_monitor_exact() not called")
    end if
    if (have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/iteration_vtus') &
       .and. .not. petsc_monitor_iteration_vtus) then
       ewrite(0,*) "Solver option diagnostics/monitors/iteration_vtus set but "
       ewrite(0,*) "petsc_solve_monitor_iteration_vtus() not called. This probably means"
       ewrite(0,*) "this version of petsc_solve() doesn't support the monitor."
       FLExit("petsc_solve_monitor_iteration_vtus() not called")
    end if
    if(have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/true_error') .or. &
       have_option(trim(solver_option_path)// &
       '/diagnostics/monitors/iteration_vtus')) then
       if (.not. present(petsc_numbering)) then
         FLAbort("Need petsc_numbering for monitor")
       end if
       call petsc_monitor_setup(petsc_numbering, max_its)
       call KSPMonitorSet(ksp,MyKSPMonitor,PETSC_NULL_OBJECT, &
            &                     PETSC_NULL_FUNCTION,ierr)
    end if

#if PETSC_VERSION_MINOR<6
    if (mat/=pmat) then
      ! This is to make things consistent with the situation in petsc>=3.6:
      ! there the nullspace is picked up directly from the Mat during KSPSolve.
      ! The routine KSPSetNullSpace no longer exists. Previously
      ! however, the nullspace was picked up from *pmat* inside KSPSetOperators
      ! At this point KSPSetOperators, has already been called, so if mat has 
      ! a nullspace we want it to be set as the nullspace of the KSP
      call MatGetNullSpace(mat, nullsp, ierr)
      if (ierr==0  .and. nullsp/=PETSC_NULL_OBJECT) then
        call KSPSetNullSpace(ksp, nullsp, ierr)
      else
        call MatGetNullSpace(pmat, nullsp, ierr)
        if (ierr==0  .and. nullsp/=PETSC_NULL_OBJECT) then
          FLAbort("Preconditioner matrix has nullspace whereas the matrix itself doesn't")
          ! This is a problem because the nullspace on the preconditioner matrix is now
          ! attached to the ksp already. Not sure how to remove it again; Can I just call
          ! KSPSetNullSpace with PETSC_NULL_OBJECT? I don't think this combination
          ! does anything useful anyway, so let's just error. You can try it out with 
          ! PETSc>=3.6 which should do the right thing.
        end if
      end if
    end if
#elif PETSC_VERSION_MINOR==6 && PETSC_VERSION_SUBMINOR==0
    ! this is 3.6.0 case where KSPSetNullSpace no longer exists, but the nullspace picked up
    ! in the krylov iteration is from *pmat* not mat (as it is in 3.6.1 and later)
    if (mat/=pmat) then
      call MatGetNullSpace(mat, nullsp, ierr)
      if (ierr==0  .and. nullsp/=PETSC_NULL_OBJECT) then
        ewrite(0,*) "Matrix and preconditioner matrix are different. For this case nullspaces"
        ewrite(0,*) "and petsc 3.6.0 are not supported. Please upgrade to petsc 3.6.1 or higher"
        FLExit("Cannot use petsc 3.6.0 with nullspaces when mat/=pmat")
      end if
    end if
#endif

  end subroutine setup_ksp_from_options

  recursive subroutine attach_null_space_from_options(mat, solver_option_path, pmat, &
      positions, rotation_matrix, petsc_numbering)
    !!< attach nullspace and multigrid near-nullspace 
    !!< if specified in solver options
    ! Petsc mat to attach nullspace to
    Mat, intent(inout):: mat
    ! path to solver block (including '/solver')
    character(len=*), intent(in):: solver_option_path
    ! the pmat (only required if different from mat)
    Mat, intent(inout), optional:: pmat
    ! positions field is only used for nullspaces with rotational components
    type(vector_field), intent(in), optional :: positions
    ! with rotated bcs: matrix to transform from x,y,z aligned vectors to boundary aligned
    Mat, intent(in), optional:: rotation_matrix
    type(petsc_numbering_type), optional, intent(in):: petsc_numbering

    MatNullSpace :: null_space
    PetscErrorCode :: ierr
    logical :: different_pmat

    if (present(pmat)) then
      different_pmat = mat/=pmat
    else
      different_pmat = .false.
    end if


    if(have_option(trim(solver_option_path)//"/multigrid_near_null_space")) then

       ! Check that we are using the gamg preconditioner:
       if(.not.(have_option(trim(solver_option_path)//"/preconditioner::gamg"))) then
          FLExit("multigrid_near_null_space removal only valid when using gamg preconditioner")
       end if

       if (.not. present(petsc_numbering)) then
          FLAbort("Need petsc_numbering for multigrid near null space")
       end if
       null_space = create_null_space_from_options_vector(mat, trim(solver_option_path)//"/multigrid_near_null_space", &
          petsc_numbering, positions=positions, rotation_matrix=rotation_matrix)
       if (different_pmat) then
         ! nns is only used in the preconditioner, so let's attach it to pmat only
         call MatSetNearNullSpace(pmat, null_space, ierr)
       else
         call MatSetNearNullSpace(mat, null_space, ierr)
       end if
       call MatNullSpaceDestroy(null_space, ierr)
    end if

    if (have_option(trim(solver_option_path)//'/remove_null_space')) then
       if (.not. present(petsc_numbering)) then
          FLAbort("Need petsc_numbering for null space removal")
       end if
       ewrite(2,*) "Attaching nullspace to matrix"
       if (size(petsc_numbering%gnn2unn,2)==1) then
          null_space = create_null_space_from_options_scalar(mat, trim(solver_option_path)//"/remove_null_space")
       else
          null_space = create_null_space_from_options_vector(mat, trim(solver_option_path)//"/remove_null_space", &
             petsc_numbering, positions=positions, rotation_matrix=rotation_matrix)
       end if
       call MatSetNullSpace(mat, null_space, ierr)
       call MatNullSpaceDestroy(null_space, ierr)
       if (have_option(trim(solver_option_path)//"/preconditioner::ksp/remove_null_space")) then
          if (.not. present(pmat)) then
             ewrite(-1,*) "For solver options specified at:", trim(solver_option_path)
             ewrite(-1,*) "Nullspace options were found under preconditioner::ksp but no"
             ewrite(-1,*) "separate preconditioner matrix is available."
             FLExit("Cannot use nullspace options under precondtioner::ksp for this solve")
          end if
          call attach_null_space_from_options(pmat, trim(solver_option_path)//"/preconditioner::ksp", &
             positions=positions, rotation_matrix=rotation_matrix, petsc_numbering=petsc_numbering)
       end if
    end if

    ! Check for nullspace options under the ksp preconditioner (near nullspaces aren't currently allowed in the schema)
    if (have_option(trim(solver_option_path)//"/preconditioner::ksp/solver/remove_null_space")) then
       if (.not. different_pmat) then
          ewrite(-1,*) "For solver options specified at:", trim(solver_option_path)
          ewrite(-1,*) "Nullspace options were found under preconditioner::ksp but no"
          ewrite(-1,*) "separate preconditioner matrix is available."
          FLExit("Cannot use nullspace options under precondtioner::ksp for this solve")
       end if
       call attach_null_space_from_options(pmat, trim(solver_option_path)//"/preconditioner::ksp", &
          positions=positions, rotation_matrix=rotation_matrix, petsc_numbering=petsc_numbering)
    end if
    
  end subroutine attach_null_space_from_options
    
  recursive subroutine setup_pc_from_options(pc, pmat, option_path, &
    petsc_numbering, prolongators, surface_node_list, matrix_csr, &
    internal_smoothing_option, is_subpc)
  PC, intent(inout):: pc
  Mat, intent(in):: pmat
  character(len=*), intent(in):: option_path
  ! needed for fieldsplit and to be pass down to pcksp:
  type(petsc_numbering_type), optional, intent(in):: petsc_numbering
  ! additional information for multigrid smoother:
  type(petsc_csr_matrix), dimension(:), optional, intent(in) :: prolongators
  integer, dimension(:), optional, intent(in) :: surface_node_list
  type(csr_matrix), optional, intent(in) :: matrix_csr
  integer, optional, intent(in) :: internal_smoothing_option
  ! if present and true, don't setup sor and eisenstat as subpc (again)
  logical, optional, intent(in) :: is_subpc
    
    KSP:: subksp
    PC:: subpc
    MatNullSpace:: nullsp
    PCType:: pctype, hypretype
    MatSolverPackage:: matsolverpackage
    PetscErrorCode:: ierr
    
    call get_option(trim(option_path)//'/name', pctype)

    if (pctype==PCMG) then
      call SetupMultigrid(pc, pmat, ierr, &
            external_prolongators=prolongators, &
            surface_node_list=surface_node_list, &
            matrix_csr=matrix_csr, &
            internal_smoothing_option=internal_smoothing_option)
      if (ierr/=0) then
         if (IsParallel()) then
           ! we give up as SOR is probably not good enough either
           ! for big paralel problems
           ewrite(-1,*) 'Set up of mg preconditioner failed for:'
           ewrite(-1,*) trim(option_path)
           FLExit("MG failed: try another preconditioner or improve partitioning")
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
      FLExit("The fluidity binary is built without hypre support!")
#endif

    else if (pctype==PCKSP) then
      
       ! this replaces the preconditioner by a complete solve
       ! (based on the pmat matrix)
       call PCSetType(pc, pctype, ierr)
       
       ! set the options for the ksp of this complete solve
       call PCKSPGetKSP(pc, subksp, ierr)
       ewrite(1,*) "Going into setup_ksp_from_options again to set the options "//&
          &"for the complete ksp solve of the preconditioner"
       call KSPSetOperators(subksp, pmat, pmat, ierr)
       call setup_ksp_from_options(subksp, pmat, pmat, &
         trim(option_path)//'/solver', petsc_numbering=petsc_numbering)
       ewrite(1,*) "Returned from setup_ksp_from_options for the preconditioner solve, "//&
          &"now setting options for the outer solve"
      
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
         petsc_numbering=petsc_numbering, &
         prolongators=prolongators, surface_node_list=surface_node_list, &
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

    else if (pctype==PCFIELDSPLIT) then

       if (.not. present(petsc_numbering)) then
         FLAbort("Need to pass down petsc numbering to set up fieldsplit")
       end if

       call setup_fieldsplit_preconditioner(pc, option_path, &
            petsc_numbering=petsc_numbering)

    else
      
       ! this doesn't work for hypre
       call PCSetType(pc, pctype, ierr)
       ! set options that may have been supplied via the
       ! PETSC_OPTIONS env. variable for the preconditioner
       call PCSetFromOptions(pc, ierr)
       ! set pctype again to enforce flml choice
       call PCSetType(pc, pctype, ierr)

       if (pctype==PCLU) then
          call get_option(trim(option_path)//'/factorization_package/name', matsolverpackage)
          call PCFactorSetMatSolverPackage(pc, matsolverpackage, ierr)
       end if

      if (pctype==PCGAMG) then
        ! we think this is a more useful default - the default value of 0.0
        ! causes spurious "unsymmetric" failures as well
        call PCGAMGSetThreshold(pc, 0.01, ierr)
        ! this was the old default:
        call PCGAMGSetCoarseEqLim(pc, 800, ierr)
        ! PC setup seems to be required so that the Coarse Eq Lim option is used.
        call PCSetup(pc,ierr)

        call MatGetNullSpace(pmat, nullsp, ierr)
        if (ierr==0  .and. nullsp/=PETSC_NULL_OBJECT) then
          ! if the preconditioner matrix has a nullspace, this may still be present
          ! at the coarsest level (the constant null vector always will be, the rotational
          ! are as well if a near-null-space is provided). In this case the default of 
          ! using a direct solver at the coarsest level causes issues. Instead we use
          ! a fixed number of SOR iterations
          call PCMGGETCoarseSolve(pc, subksp, ierr)
          call KSPSetType(subksp, KSPPREONLY, ierr)
          call KSPGetPC(subksp, subpc, ierr)
          call PCSetType(subpc, PCSOR, ierr)
          call KSPSetTolerances(subksp, 1e-50, 1e-50, 1e50, 10, ierr)
        end if
      end if
      
    end if

    ewrite(2, *) 'pc_type: ', trim(pctype)
    if (pctype=='hypre') then
      ewrite(2,*) 'pc_hypre_type:', trim(hypretype)
    end if
    
  end subroutine setup_pc_from_options

  recursive subroutine setup_fieldsplit_preconditioner(pc, option_path, &
    petsc_numbering)
  PC, intent(inout):: pc
  character(len=*), intent(in):: option_path
  type(petsc_numbering_type), intent(in):: petsc_numbering

    character(len=128):: fieldsplit_type
    KSP, dimension(size(petsc_numbering%gnn2unn,2)):: subksps
    Mat :: mat, pmat
    MatNullSpace :: null_space
    IS:: index_set
    PetscErrorCode:: ierr
    integer:: i, n

    call PCSetType(pc, "fieldsplit", ierr)

    call PCFieldSplitGetSubKSP(pc, n, subksps, ierr)
    if (n==0) then
      ! first time this pc set to type fieldplit: it's the first time we set it up,
      ! or it was previously set to a different type - in this case, PCSetType will
      ! have called PCCreate_FieldSplit which will have set n/o splits to zero
      do i=1, size(subksps)
        index_set = petsc_numbering_create_is(petsc_numbering, dim=i)
        call PCFieldSplitSetIS(pc, PETSC_NULL_CHARACTER, index_set, ierr)
        call ISDestroy(index_set, ierr)
      end do

    elseif (n/=size(subksps)) then

      ! if this pc is reused (and we've previously already set it up with fieldsplit)
      ! we need to check the n/o fieldsplits is the same

      FLAbort("PC being reused with different number of fieldsplits")

    end if

    call get_option(trim(option_path)//"/fieldsplit_type/name", &
      fieldsplit_type, ierr)
    select case (fieldsplit_type)
    case ("multiplicative")
      call pcfieldsplitsettype(pc, PC_COMPOSITE_MULTIPLICATIVE, ierr)
    case ("additive")
      call pcfieldsplitsettype(pc, PC_COMPOSITE_ADDITIVE, ierr)
    case ("symmetric_multiplicative")
      call pcfieldsplitsettype(pc, PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE, ierr)
    case default
      FLAbort("Unknown fieldsplit_type")
    end select

    call pcfieldsplitgetsubksp(pc, n, subksps, ierr)
    
    assert(n==size(subksps))

    do i=1, size(subksps)

      call KSPGetOperators(subksps(i), mat, pmat, ierr)
      if (have_option(trim(option_path)//"/remove_null_space")) then
        null_space = create_null_space_from_options_scalar(mat, trim(option_path)//"/remove_null_space")
        call MatSetNullSpace(mat, null_space, ierr)
        call MatNullSpaceDestroy(null_space, ierr)
      end if
      call setup_ksp_from_options(subksps(i), mat, pmat, option_path)

    end do

  end subroutine setup_fieldsplit_preconditioner
    
  subroutine ewrite_ksp_options(ksp)
    KSP, intent(in):: ksp
    
    PC:: pc
    KSPType:: ksptype
    PCType:: pctype
    PetscReal:: rtol, atol, dtol
    PetscInt:: maxits
    PetscBool:: flag
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

subroutine petsc_monitor_setup(petsc_numbering, max_its)
  ! sets up the petsc monitors "exact" or "iteration_vtus"
  type(petsc_numbering_type), intent(in):: petsc_numbering
  integer, intent(in) :: max_its
  
  type(mesh_type), pointer:: mesh
  integer :: ierr, ncomponents

  petsc_monitor_x=PetscNumberingCreateVec(petsc_numbering)
  petsc_monitor_numbering=petsc_numbering
  ncomponents=size(petsc_numbering%gnn2unn,2)
       
  if (petsc_monitor_has_exact) then
    
    call VecDuplicate(petsc_monitor_x, petsc_monitor_exact, ierr)
    
    if (ncomponents==1) then
      call field2petsc(petsc_monitor_exact_sfield, petsc_numbering, petsc_monitor_exact)
    else
      call field2petsc(petsc_monitor_exact_vfield, petsc_numbering, petsc_monitor_exact)
    end if
    
    allocate( petsc_monitor_error(max_its+1) )
    allocate( petsc_monitor_flops(max_its+1) )
    petsc_monitor_error = 0.0
    petsc_monitor_flops = 0.0
    petsc_monitor_iteration=0
    
  end if
  
  if (petsc_monitor_iteration_vtus) then
    mesh => petsc_monitor_positions%mesh
    if (ncomponents==1) then
      call allocate(petsc_monitor_sfields(1), mesh, name="X")
      call allocate(petsc_monitor_sfields(2), mesh, name="Residual")
      call allocate(petsc_monitor_sfields(3), mesh, name="PreconditionedResidual")
    else
      call allocate(petsc_monitor_vfields(1), ncomponents, mesh, name="X")
      call allocate(petsc_monitor_vfields(2), ncomponents, mesh, name="Residual")
      call allocate(petsc_monitor_vfields(3), ncomponents, mesh, name="PreconditionedResidual")
    end if
    petsc_monitor_vtu_name="petsc_monitor_"//int2str(petsc_monitor_vtu_series)
    petsc_monitor_vtu_series=petsc_monitor_vtu_series+1
  end if

end subroutine petsc_monitor_setup
  
subroutine petsc_solve_monitor_exact_scalar(exact, error_filename)
!! To be called before petsc_solve. Registers the exact solution field
!! to which the approximate solutions are compared each iteration.
type(scalar_field), intent(in):: exact
! if present write to this filename, otherwise writes to stdout
character(len=*), optional, intent(in):: error_filename
  
  petsc_monitor_exact_sfield=exact
  call incref(petsc_monitor_exact_sfield)
  if (present(error_filename)) then
    petsc_monitor_error_filename=error_filename
  else
    petsc_monitor_error_filename=""
  end if
  petsc_monitor_has_exact=.true.
  
end subroutine petsc_solve_monitor_exact_scalar
  
subroutine petsc_solve_monitor_iteration_vtus(positions)
!! To be called before petsc_solve. Registers the position field to be
!! used in the vtus written out by the "iteration_vtus" monitor. Needs
!! to be the exact same mesh as the solution field.
type(vector_field), intent(in):: positions
  
  petsc_monitor_positions=positions
  call incref(petsc_monitor_positions)
  petsc_monitor_iteration_vtus=.true.
  
end subroutine petsc_solve_monitor_iteration_vtus

subroutine petsc_monitor_destroy()
  ! Destroys everything asscociated with the petsc monitors
  integer :: ierr, i
  integer :: error_unit
  character(len = 100) :: format0

  call VecDestroy(petsc_monitor_x, ierr)
  
  if (petsc_monitor_has_exact) then

    if(petsc_monitor_error_filename/='') then
       !dumping out errors and flops
       error_unit=free_unit()
       open(unit=error_unit, file=trim(petsc_monitor_error_filename), action="write")    
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
       ewrite(1,*) 'iteration, error, flops'
       do i = 1, petsc_monitor_iteration
          ewrite(1,*) i, petsc_monitor_error(i), petsc_monitor_flops(i)
       end do
    end if
    
    if (size(petsc_monitor_numbering%gnn2unn,2)==1) then
      call deallocate(petsc_monitor_exact_sfield)
    else
      call deallocate(petsc_monitor_exact_vfield)
    end if

    call VecDestroy(petsc_monitor_exact, ierr)
    deallocate( petsc_monitor_error )
    deallocate( petsc_monitor_flops )
  end if
  
  if (petsc_monitor_iteration_vtus) then
    if (size(petsc_monitor_numbering%gnn2unn,2)==1) then
      do i=1, 3
        call deallocate(petsc_monitor_sfields(i))
      end do
    else
      do i=1, 3
        call deallocate(petsc_monitor_vfields(i))
      end do
    end if
    call deallocate(petsc_monitor_positions)
  end if

end subroutine petsc_monitor_destroy
  
subroutine MyKSPMonitor(ksp,n,rnorm,dummy,ierr)
!! The monitor function that gets called each iteration of petsc_solve
!! (if petsc_solve_callback_setup is called)
  PetscInt, intent(in) :: n,dummy
  KSP, intent(in) :: ksp
  PetscErrorCode, intent(out) :: ierr
  
  PetscScalar :: rnorm
  MatNullSpace :: nullsp
  PetscLogDouble :: flops
  Mat:: Amat, Pmat
  PC:: pc
  Vec:: dummy_vec, r, rhs

  !  Build the solution vector  
  call VecZeroEntries(petsc_monitor_x,ierr)
  ! don't pass PETSC_NULL_OBJECT instead of dummy_vec, as petsc
  ! will clobber it (bug in fortran interface)
  call KSPBuildSolution(ksp,petsc_monitor_x, dummy_vec, ierr)

  if (petsc_monitor_has_exact) then
    ! Compare with exact solution
    call VecAXPY(petsc_monitor_x, real(-1.0, kind = PetscScalar_kind), petsc_monitor_exact, ierr)
    call VecNorm(petsc_monitor_x, NORM_INFINITY, rnorm, ierr)
    call PetscGetFlops(flops,ierr)

    petsc_monitor_error(n+1) = rnorm
    petsc_monitor_iteration = max(petsc_monitor_iteration,n+1)
    petsc_monitor_flops(n+1) = flops
  end if
  
  if (petsc_monitor_iteration_vtus) then
    ! store the solution
    if (size(petsc_monitor_numbering%gnn2unn,2)==1) then
      call petsc2field(petsc_monitor_x, petsc_monitor_numbering, petsc_monitor_sfields(1))
    else
      call petsc2field(petsc_monitor_x, petsc_monitor_numbering, petsc_monitor_vfields(1))
    end if

    ! then (re)compute the (true) residual
    call KSPGetRhs(ksp, rhs, ierr)
    call KSPGetOperators(ksp, Amat, Pmat, ierr)
    call VecDuplicate(petsc_monitor_x, r, ierr)
    call MatMult(Amat, petsc_monitor_x, r, ierr)
    call VecAXPY(r, real(-1.0, kind = PetscScalar_kind), rhs, ierr)
    if (size(petsc_monitor_numbering%gnn2unn,2)==1) then
      call petsc2field(r, petsc_monitor_numbering, petsc_monitor_sfields(2))
    else
      call petsc2field(r, petsc_monitor_numbering, petsc_monitor_vfields(2))
    end if

    ! now (re)compute the preconditioned residual - which is what we usually look at for convergence
    call VecCopy(r, petsc_monitor_x, ierr)
    call KSPGetPC(ksp, pc, ierr)
    call PCApply(pc, petsc_monitor_x, r, ierr)
    ! within petsc the nullspace is removed directly after pcapply (see KSP_PCApply)
    call MatGetNullSpace(Pmat, nullsp, ierr)
    if (ierr==0  .and. nullsp/=PETSC_NULL_OBJECT) then
      call MatNullSpaceRemove(nullsp, r, PETSC_NULL_OBJECT, ierr)
    end if
    if (size(petsc_monitor_numbering%gnn2unn,2)==1) then
      call petsc2field(r, petsc_monitor_numbering, petsc_monitor_sfields(3))
    else
      call petsc2field(r, petsc_monitor_numbering, petsc_monitor_vfields(3))
    end if

    if (size(petsc_monitor_numbering%gnn2unn,2)==1) then
      call vtk_write_fields(petsc_monitor_vtu_name, index=n, &
        model=petsc_monitor_positions%mesh, position=petsc_monitor_positions, &
        sfields=petsc_monitor_sfields)
    else
      call vtk_write_fields(petsc_monitor_vtu_name, index=n, &
        model=petsc_monitor_positions%mesh, position=petsc_monitor_positions, &
        vfields=petsc_monitor_vfields)
    end if
    call VecDestroy(r, ierr)
  end if
  
  ierr=0
  
end subroutine MyKSPMonitor

function create_null_space_from_options_scalar(mat, null_space_option_path) &
    result (null_space)

   Mat, intent(in):: mat
   !! the option path to remove_null_space
   character(len=*), intent(in):: null_space_option_path

   ! hack to satisfy interface for MatNullSpaceCreate
   ! only works as the array won't actually be used
   PetscObject, dimension(1:0) :: PETSC_NULL_OBJECT_ARRAY
   MatNullSpace :: null_space
   PetscErrorCode :: ierr
   PetscBool :: isnull

   call MatNullSpaceCreate(MPI_COMM_FEMTOOLS, PETSC_TRUE, 0, PETSC_NULL_OBJECT_ARRAY, null_space, ierr)

   if(have_option(trim(null_space_option_path)//'/test_null_space')) then
     call MatNullSpaceTest(null_space, mat, isnull, ierr)
     ewrite(1,*) "For nullspace "//trim(null_space_option_path)//":"
     if (isnull) then
       ewrite(1,*) "PETSc's MatNullSpaceTest agrees that this is a null space"
     else
       ewrite(1,*) "PETSc's MatNullSpaceTest does not think this is a null space"
     end if
   end if

end function create_null_space_from_options_scalar

function create_null_space_from_options_vector(mat, null_space_option_path, &
      petsc_numbering, positions, rotation_matrix) result (null_space)
   Mat, intent(in):: mat
   !! the option path to remove_null_space or multigrid_near_space
   character(len=*), intent(in):: null_space_option_path
   type(petsc_numbering_type), intent(in):: petsc_numbering 
   ! positions field is only used with remove_null_space/ with rotational components
   type(vector_field), intent(in), optional :: positions
   ! with rotated bcs: matrix to transform from x,y,z aligned vectors to boundary aligned
   Mat, intent(in), optional:: rotation_matrix
   MatNullSpace :: null_space

   Vec, allocatable, dimension(:) :: null_space_array
   Vec :: aux_vec, swap
   PetscReal, dimension(:), allocatable :: dots
   PetscReal :: norm
   PetscErrorCode :: ierr
   PetscBool :: isnull

   integer :: i, nnulls, nnodes, comp, dim, universal_nodes
   logical, dimension(3) :: rot_mask
   logical, dimension(size(petsc_numbering%gnn2unn,2)) :: mask
   real, dimension(:,:), allocatable :: null_vector
   type(vector_field) :: null_vector_field
   type(vector_field), allocatable, dimension(:) :: vtk_vector_fields(:)
   integer, save :: vtk_index = 0

   integer, dimension(5), parameter:: permutations=(/ 1,2,3,1,2 /)

   nnodes=size(petsc_numbering%gnn2unn,1)
   dim=size(petsc_numbering%gnn2unn,2)

   if(have_option(trim(null_space_option_path)//'/specify_components')) then
     ! count how many null spaces we want
     mask = .false.
     mask(1) = have_option(trim(null_space_option_path)//&
                            &"/specify_components/x_component")
     if (have_option(trim(null_space_option_path)//&
                            &"/specify_components/y_component")) then
       if(dim<2) then
         FLExit("Requested the removal of a y component null space on a less than 2d vector.")
       end if
       mask(2) = .true.
     end if
     if (have_option(trim(null_space_option_path)//&
                            &"/specify_components/z_component")) then
       if(dim<3) then
         FLExit("Requested the removal of a z component null space on a less than 3d vector.")
       end if
       mask(3) = .true.
     end if
     if(.not. any(mask)) then
       FLExit("Requested null space removal on specific components but have not specified which components.")
     end if
   else if(have_option(trim(null_space_option_path)//'/all_components')) then
     mask = .true.
   else if(have_option(trim(null_space_option_path)//'/no_components')) then
     mask = .false.
   end if

   if(have_option(trim(null_space_option_path)//'/specify_rotations')) then
     rot_mask = .false.
     rot_mask(3)=have_option(trim(null_space_option_path)//&
                            &"/specify_rotations/xy_rotation")
     if (have_option(trim(null_space_option_path)//&
                            &"/specify_rotations/xz_rotation")) then

       if(dim<3) then
         FLExit("Requested the removal of xz rotation on a less than 3d vector.")
       end if
       rot_mask(2) = .true.
     end if
     if (have_option(trim(null_space_option_path)//&
                            &"/specify_rotations/yz_rotation")) then

       if(dim<3) then
         FLExit("Requested the removal of yz rotation on a less than 3d vector.")
       end if
       rot_mask(1) = .true.
     end if
     if(.not. any(rot_mask)) then
       FLExit("Requested null space removal on specific rotations but have not specified which rotations.")
     end if
   else if(have_option(trim(null_space_option_path)//'/all_rotations')) then
     rot_mask = .false.
     if (dim==3) then
       rot_mask = .true.
     else if (dim==2) then
       rot_mask(3) = .true.
     end if
   else
     rot_mask = .false.
   end if

   if (.not. (any(mask) .or. any(rot_mask)) ) then
      FLExit("You must remove either a component or a rotation.")
   end if

   if(any(rot_mask) .and. dim<2) then
     FLExit("Requested the removal of rotational component for a less than 2d vector.")
   end if 

   nnulls=count(mask)+count(rot_mask)
   ! allocate the array of null spaces
   allocate(null_space_array(1:nnulls))
   allocate(null_vector(nnodes,dim))
   universal_nodes=petsc_numbering%universal_length/dim

   ewrite(2,*) "Setting up array of "//int2str(nnulls)//" null spaces."
   
   ! now loop back over the components building up the null spaces we want
   i = 0
   do comp = 1, dim
     if (mask(comp)) then
       i = i + 1
       null_vector = 0.0
       ! ensure the translations are orthonormal:
       null_vector(:,comp)=1.0/sqrt(real(universal_nodes))
       null_space_array(i)=PetscNumberingCreateVec(petsc_numbering)
       call array2petsc(reshape(null_vector,(/nnodes*dim/)), petsc_numbering, null_space_array(i))
     end if
   end do

   deallocate(null_vector)

   if (any(rot_mask)) then

     if (.not. present(positions)) then
       ! when providing the option to remove rotational modes to the user, positions should have been passed in
       FLAbort("In create_null_space_from_options, need positions field")
     end if

     call allocate(null_vector_field, positions%dim, positions%mesh)

     do comp = 1, 3
       if (rot_mask(comp)) then
         i = i + 1
         if (dim==3) then
           ! for dim==2: comp==3 and both components will be set already
           call zero(null_vector_field, comp)
         end if
         call set(null_vector_field, permutations(comp+1), &
           extract_scalar_field(positions, permutations(comp+2)))
         call scale(null_vector_field, -1.0, dim=permutations(comp+1))
         call set(null_vector_field, permutations(comp+2), &
           extract_scalar_field(positions, permutations(comp+1)))
         null_space_array(i)=PetscNumberingCreateVec(petsc_numbering)
         call field2petsc(null_vector_field, petsc_numbering, null_space_array(i))
       end if
     end do
     call deallocate(null_vector_field)

   end if

   assert(i==nnulls)

   if (present(rotation_matrix) .and. nnulls>0) then
     call VecDuplicate(null_space_array(1), aux_vec, ierr)
     do i=1, nnulls
       ! rotate the null vector and store it in aux_vec
       call MatMultTranspose(rotation_matrix, null_space_array(i), aux_vec, ierr)
       ! swap the unrotated null_space_array(i) with aux_vec
       ! so that we store the rotated one in null_space_array(i) 
       ! and can use the unrotated as aux_vec in the next iteration
       swap = null_space_array(i)
       null_space_array(i) = aux_vec
       aux_vec = swap
     end do
     call VecDestroy(aux_vec, ierr)
   end if

   ! finally we need to ensure that the nullspace vectors are orthonormal
   if (any(rot_mask)) then
     ! but only the rotational ones, as the translations are orthonormal already
     allocate(dots(1:nnulls))
     do i=count(mask)+1, nnulls
       ! take the dot product with all previous vectors:
       call VecMDot(null_space_array(i), i-1, null_space_array(1:i-1), dots(1:i-1), ierr)
       dots = -dots
       ! then subtract their components
       call VecMAXPY(null_space_array(i), i-1, dots(1:i-1), null_space_array(1:i-1), ierr)
       call VecNormalize(null_space_array(i), norm, ierr)
     end do
     deallocate(dots)
   end if

   call MatNullSpaceCreate(MPI_COMM_FEMTOOLS, PETSC_FALSE, nnulls, &
       null_space_array, null_space, ierr)

   if(have_option(trim(null_space_option_path)//'/test_null_space')) then
     call MatNullSpaceTest(null_space, mat, isnull, ierr)
     ewrite(1,*) "For nullspace "//trim(null_space_option_path)//":"
     if (isnull) then
       ewrite(1,*) "PETSc's MatNullSpaceTest agrees that this is a null space"
     else
       ewrite(1,*) "PETSc's MatNullSpaceTest does not think this is a null space"
     end if
   end if

   if(have_option(trim(null_space_option_path)//'/write_null_space')) then
     allocate(vtk_vector_fields(1:nnulls))
     do i=1, nnulls
       call allocate(vtk_vector_fields(i), positions%dim, positions%mesh, name="NullVector"//int2str(i))
       call petsc2field(null_space_array(i), petsc_numbering, vtk_vector_fields(i))
     end do
     vtk_index = vtk_index + 1
     call vtk_write_fields("null_space", index=vtk_index, &
       model=positions%mesh, position=positions, &
       vfields=vtk_vector_fields)
     do i=1, nnulls
       call deallocate(vtk_vector_fields(i))
     end do
     deallocate(vtk_vector_fields)
   end if

   ! get rid of our Vec references
   do i=1, nnulls
     call VecDestroy(null_space_array(i), ierr)
   end do
   deallocate(null_space_array)


end function create_null_space_from_options_vector

function petsc_solve_needs_positions(solver_option_path)
  !!< Auxillary function to tell us if we need to pass in a positions field to petsc_solve
  !!< Currently only for vector solves with remove_null_space or multigrid_near_null_space
  !!< with specify_rotations or all_rotations
  character(len=*), intent(in):: solver_option_path
  logical:: petsc_solve_needs_positions

  petsc_solve_needs_positions = &
    have_option(trim(solver_option_path)//'/remove_null_space/specify_rotations') .or. &
    have_option(trim(solver_option_path)//'/remove_null_space/all_rotations') .or. &
    have_option(trim(solver_option_path)//'/multigrid_near_null_space/specify_rotations') .or. &
    have_option(trim(solver_option_path)//'/multigrid_near_null_space/all_rotations')

end function petsc_solve_needs_positions

end module solvers
