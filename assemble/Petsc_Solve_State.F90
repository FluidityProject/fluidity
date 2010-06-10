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
use field_options
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
    
    integer, dimension(:), pointer:: surface_nodes
    type(csr_matrix):: prolongator
    character(len=OPTION_PATH_LEN):: solver_option_path
    
    type(mesh_type), pointer:: linear_mesh
    
    call petsc_solve_state_setup(solver_option_path, state, x, &
      option_path=option_path)
    
    if (have_option(trim(solver_option_path)//'/preconditioner::mg/vertical_lumping')) then
      
      prolongator=vertical_prolongator_from_free_surface(state, x%mesh)
      
      if (have_option(trim(solver_option_path)//'/preconditioner::mg/vertical_lumping/internal_smoother')) then
        
        surface_nodes => free_surface_nodes(state, x%mesh)
        
        call petsc_solve(x, matrix, rhs, &
           prolongator=prolongator, &
           surface_node_list=surface_nodes, option_path=option_path)
           
        deallocate(surface_nodes)
        
      else
      
        call petsc_solve(x, matrix, rhs, &
           prolongator=prolongator, option_path=option_path)
        
      end if
      
      call deallocate(prolongator)
    
    else if (have_option(trim(solver_option_path)//'/preconditioner::mg/higher_order_lumping')) then
    
      call find_linear_parent_mesh(state, x%mesh, linear_mesh)
      prolongator=higher_order_prolongator(linear_mesh, x%mesh)
    
      call petsc_solve(x, matrix, rhs, option_path=option_path, prolongator=prolongator)
      
      call deallocate(prolongator)
      
    else
    
      call petsc_solve(x, matrix, rhs, option_path=option_path)
      
    end if
    
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
    
    type(mesh_type), pointer:: linear_mesh
    integer, dimension(:), pointer:: surface_nodes
    type(csr_matrix):: prolongator
    character(len=OPTION_PATH_LEN):: solver_option_path
    
    call petsc_solve_state_setup(solver_option_path, state, x, &
      option_path=option_path)
    
    if (have_option(trim(solver_option_path)//'/preconditioner::mg/vertical_lumping')) then
      
      prolongator=vertical_prolongator_from_free_surface(state, x%mesh)
      
      if (have_option(trim(solver_option_path)//'/preconditioner::mg/vertical_lumping/internal_smoother')) then
        
        surface_nodes => free_surface_nodes(state, x%mesh)
        
        call petsc_solve(x, matrix, rhs, &
           prolongator=prolongator, &
           surface_node_list=surface_nodes, option_path=option_path)
           
        deallocate(surface_nodes)
        
      else
      
        call petsc_solve(x, matrix, rhs, &
           prolongator=prolongator, option_path=option_path)
        
      end if
      
      call deallocate(prolongator)
    
    else if (have_option(trim(solver_option_path)//'/preconditioner::mg/higher_order_lumping')) then
    
      call find_linear_parent_mesh(state, x%mesh, linear_mesh)
      prolongator=higher_order_prolongator(linear_mesh, x%mesh)
    
      call petsc_solve(x, matrix, rhs, option_path=option_path, prolongator=prolongator)
      
      call deallocate(prolongator)
      
    else
    
      call petsc_solve(x, matrix, rhs, option_path=option_path)
      
    end if
    
  end subroutine petsc_solve_scalar_state_petsc_csr
    
  subroutine petsc_solve_state_setup(solver_option_path, state, x, option_path)
    ! sets up monitors and returns solver_option_path
    character(len=*), intent(out):: solver_option_path
    type(state_type), intent(in):: state
    type(scalar_field), intent(in):: x
    character(len=*), intent(in), optional:: option_path
    
    type(vector_field):: positions
    type(scalar_field), pointer:: exact
    character(len=FIELD_NAME_LEN):: exact_field_name
    integer:: stat
    
    if (present(option_path)) then
       solver_option_path=complete_solver_option_path(option_path)
    else
       solver_option_path=complete_solver_option_path(x%option_path)
    end if
    
    call get_option(trim(solver_option_path)// &
        '/diagnostics/monitors/true_error/exact_solution_field', &
        exact_field_name, stat=stat)
    if (stat==0) then
       exact => extract_scalar_field(state, exact_field_name)
       call petsc_solve_monitor_exact(exact)
    end if
    
    if (have_option(trim(solver_option_path)// &
        '/diagnostics/monitors/iteration_vtus')) then
       positions=get_nodal_coordinate_field(state, x%mesh)
       ! creates its own reference that's cleaned up in petsc_solve:
       call petsc_solve_monitor_iteration_vtus(positions)
       ! so we're free to get rid of ours
       call deallocate(positions)
    end if

  end subroutine petsc_solve_state_setup
    
  logical function petsc_solve_needs_state(option_path)
  ! function used in petsc_readnsolve to work out whether it needs
  ! to read state and call the above petsc_solve_state, or can just
  ! go for the simple petsc_solve instead
  character(len=*), intent(in):: option_path
  
     character(len=OPTION_PATH_LEN) solver_option_path
     
     solver_option_path=complete_solver_option_path(option_path)
     
     petsc_solve_needs_state=have_option( &
        trim(solver_option_path)//'/preconditioner::mg/vertical_lumping') &
      .or. have_option( &
        trim(solver_option_path)//'/preconditioner::mg/higher_order_lumping') &
      .or. have_option( &
        trim(solver_option_path)//'/diagnostics/monitors/true_error') &
      .or. have_option( &
        trim(solver_option_path)//'/diagnostics/monitors/iteration_vtus')
  
  end function petsc_solve_needs_state
  
  
  function higher_order_prolongator(p1_mesh, pn_mesh) result (P)
  ! Creates the linear operator that extrapolates p1 fields to higher
  ! order pn meshes. This can be used as the first stage prolongator
  ! in the "mg" multigrid preconditioner
    type(mesh_type), intent(in):: p1_mesh, pn_mesh
    type(csr_matrix):: P
    
    type(csr_sparsity):: sparsity
    real, dimension(:), pointer:: Prow
    integer, dimension(:), pointer:: p1_nodes, pn_nodes, Prow_m
    integer:: rows, columns, entries, nloc
    integer:: i, j, k, node, ele
    
    rows=nowned_nodes(pn_mesh)
    columns=node_count(p1_mesh)
    ! need to adjust this for n>2
    nloc=ele_loc(p1_mesh,1)
    entries=rows*ele_loc(p1_mesh,1)
    
    call allocate(sparsity, rows, columns, entries, name="HigherOrderProlongatorSparsity")
    j=1
    do i=1, size(sparsity%findrm)
      sparsity%findrm(i)=j
      j=j+nloc
    end do
    if (associated(p1_mesh%halos)) then
      sparsity%column_halo => p1_mesh%halos(1)
      call incref(sparsity%column_halo)
    end if
      
    call allocate(P, sparsity, name="HigherOrderProlongator")
    call deallocate(sparsity)
    
    do ele=1, ele_count(pn_mesh)
      pn_nodes => ele_nodes(pn_mesh, ele)
      p1_nodes => ele_nodes(p1_mesh, ele)
      do j=1, size(pn_nodes)
        node=pn_nodes(j)
        if (node_owned(pn_mesh, node)) then
          Prow_m => row_m_ptr(P, node)
          Prow_m = p1_nodes
          Prow => row_val_ptr(P, node)        
          do k=1, size(p1_nodes)
            Prow(k)=eval_shape(p1_mesh%shape, k, &
                    local_coords(j, pn_mesh%shape))
          end do
        end if
      end do
    end do
    
  end function higher_order_prolongator
  
end module petsc_solve_state_module