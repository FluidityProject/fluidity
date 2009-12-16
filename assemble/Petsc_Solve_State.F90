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

module petsc_solve_state_module
!!< This module provides an extension of the petsc_solve interface
!!< where state can supplied as an extra argument. This allows the use
!!< of extra geometric information pulled from state in the solver.
!!< Currently this is used for the "mg" preconditioner with 
!!< vertical_lumping option (with and without internal smoothing).
!!< This is put in a separate module as the way this information is stored
!!< in state is fluidity specific and should therefore not be dealt with
!!< in femtools/.
use spud
use sparse_tools
use fields
use solvers
use sparse_tools_petsc
use state_module
use global_parameters, only: OPTION_PATH_LEN
! modules from assemble:
use free_surface_module
implicit none

interface petsc_solve
   module procedure petsc_solve_scalar_state, petsc_solve_scalar_state_petsc_csr
end interface
  
private
public petsc_solve, petsc_solve_needs_state

contains

  subroutine petsc_solve_scalar_state(x, matrix, rhs, state, &
    option_path)
    !!< Solve a linear system the nice way.
    !!< This version uses state to pull geometric information from
    !!< if required for the specified options.
    type(scalar_field), intent(inout) :: x
    type(scalar_field), intent(in) :: rhs
    type(csr_matrix), intent(in) :: matrix
    type(state_type), intent(in):: state
    !! override x%option_path if provided:
    character(len=*), optional, intent(in):: option_path
    
    call petsc_solve_scalar_state_wrapper(x, matrix_csr=matrix, &
      rhs=rhs, state=state, option_path=option_path)
    
  end subroutine petsc_solve_scalar_state
  
  subroutine petsc_solve_scalar_state_petsc_csr(x, matrix, rhs, state, &
    option_path)
    !!< Solve a linear system the nice way.
    !!< This version uses state to pull geometric information from
    !!< if required for the specified options.
    type(scalar_field), intent(inout) :: x
    type(scalar_field), intent(in) :: rhs
    type(petsc_csr_matrix), intent(inout) :: matrix
    type(state_type), intent(in):: state
    !! override x%option_path if provided:
    character(len=*), optional, intent(in):: option_path
    
    call petsc_solve_scalar_state_wrapper(x, matrix_petsc_csr=matrix, &
      rhs=rhs, state=state, option_path=option_path)
    
  end subroutine petsc_solve_scalar_state_petsc_csr
  
  subroutine petsc_solve_scalar_state_wrapper(x, &
    matrix_csr, matrix_petsc_csr, rhs, state, &
    option_path)
    ! Wrapper for petsc_solve_scalar_state that takes either
    ! csr_matrix or petsc_csr_matrix
    type(scalar_field), intent(inout) :: x
    type(scalar_field), intent(in) :: rhs
    type(csr_matrix), optional, intent(in) :: matrix_csr
    type(petsc_csr_matrix), optional, intent(inout) :: matrix_petsc_csr
    type(state_type), intent(in):: state
    ! override x%option_path if provided:
    character(len=*), optional, intent(in):: option_path
      
    integer, dimension(:), pointer:: surface_nodes
    type(csr_matrix) prolongator
    type(scalar_field), pointer:: exact
    character(len=OPTION_PATH_LEN):: solver_option_path
    character(len=FIELD_NAME_LEN):: exact_field_name
    integer stat
      
    if (present(option_path)) then
       solver_option_path=complete_solver_option_path(option_path)
    else
       solver_option_path=complete_solver_option_path(x%option_path)
    end if
    
    call get_option( &
        trim(solver_option_path)//'/diagnostics/monitors/true_error/exact_solution_field', &
        exact_field_name, stat=stat)
    if (stat==0) then        
       exact => extract_scalar_field(state, exact_field_name)       
    else      
       nullify(exact)    
    end if
    
    if (have_option(trim(solver_option_path)//'/preconditioner::mg/vertical_lumping')) then
      
      prolongator=vertical_prolongator_from_free_surface(state, x%mesh)
      
      if (have_option(trim(solver_option_path)//'/preconditioner::mg/vertical_lumping/internal_smoother')) then
        
        surface_nodes => free_surface_nodes(state, x%mesh)
        
        call petsc_solve_exact(x, matrix_csr=matrix_csr, &
           matrix_petsc_csr=matrix_petsc_csr, rhs=rhs, exact=exact, &
           prolongator=prolongator, &
           surface_node_list=surface_nodes, option_path=option_path)
           
        deallocate(surface_nodes)
        
      else
      
        call petsc_solve_exact(x, matrix_csr=matrix_csr, &
           matrix_petsc_csr=matrix_petsc_csr, rhs=rhs, exact=exact, &
           prolongator=prolongator, option_path=option_path)
        
      end if
      
      call deallocate(prolongator)
    
    else
    
        call petsc_solve_exact(x, matrix_csr=matrix_csr, &
           matrix_petsc_csr=matrix_petsc_csr, rhs=rhs, exact=exact, &
           option_path=option_path)
      
    end if

  end subroutine petsc_solve_scalar_state_wrapper
    
  subroutine petsc_solve_exact(x, matrix_csr, matrix_petsc_csr, &
    rhs, exact, prolongator, surface_node_list, option_path)
    ! This wrapper passes on exact only if associated
    type(scalar_field), intent(inout) :: x
    type(scalar_field), intent(in) :: rhs
    type(csr_matrix), optional, intent(in) :: matrix_csr
    type(petsc_csr_matrix), optional, intent(inout) :: matrix_petsc_csr
    type(scalar_field), pointer :: exact
    ! optional arguments passed straight on:
    type(csr_matrix), optional, intent(in):: prolongator
    integer, dimension(:), optional, intent(in):: surface_node_list
    character(len=*), optional, intent(in):: option_path
      
    if (present(matrix_csr)) then
      if (associated(exact)) then
        
        call petsc_solve(x, matrix_csr, rhs, exact=exact, prolongator=prolongator, &
          surface_node_list=surface_node_list, option_path=option_path)
          
      else
      
        call petsc_solve(x, matrix_csr, rhs, prolongator=prolongator, &
          surface_node_list=surface_node_list, option_path=option_path)
          
      end if
    else if (present(matrix_petsc_csr)) then
      if (associated(exact)) then
        
        call petsc_solve(x, matrix_petsc_csr, rhs, exact=exact, prolongator=prolongator, &
          surface_node_list=surface_node_list, option_path=option_path)
          
      else
      
        call petsc_solve(x, matrix_petsc_csr, rhs, prolongator=prolongator, &
          surface_node_list=surface_node_list, option_path=option_path)
          
      end if
    else
      FLAbort("Either matrix_csr or matrix_petsc_csr required in petsc_solve_exact")
    end if
      
  end subroutine petsc_solve_exact
  
  logical function petsc_solve_needs_state(option_path)
  character(len=*), intent(in):: option_path
  
     character(len=OPTION_PATH_LEN) solver_option_path
     
     solver_option_path=complete_solver_option_path(option_path)
     
     petsc_solve_needs_state=have_option( &
        trim(solver_option_path)//'/preconditioner::mg/vertical_lumping') &
      .or. have_option( &
        trim(solver_option_path)//'/diagnostics/monitors/true_error')
  
  end function petsc_solve_needs_state
  
end module petsc_solve_state_module