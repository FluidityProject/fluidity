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
use fldebug
use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
use sparse_tools
use parallel_fields
use fields
use sparse_tools_petsc
use state_module
use solvers
use field_options
! modules from assemble:
use free_surface_module
implicit none

interface petsc_solve
   module procedure petsc_solve_scalar_state, &
       petsc_solve_scalar_state_petsc_csr, &
       petsc_solve_vector_state_petsc_csr
end interface
  
private
public petsc_solve, petsc_solve_needs_state, &
  petsc_solve_state_setup

contains

  subroutine petsc_solve_scalar_state(x, matrix, rhs, state, &
    option_path, iterations_taken)
    !!< Solve a linear system the nice way.
    !!< This version uses state to pull geometric information from
    !!< if required for the specified options.
    type(scalar_field), intent(inout) :: x
    type(scalar_field), intent(in) :: rhs
    type(csr_matrix), intent(in) :: matrix
    type(state_type), intent(in):: state
    !! override x%option_path if provided:
    character(len=*), optional, intent(in):: option_path
    !! the number of petsc iterations taken
    integer, intent(out), optional :: iterations_taken
    
    integer, dimension(:), pointer:: surface_nodes
    type(petsc_csr_matrix), dimension(:), pointer:: prolongators
    character(len=OPTION_PATH_LEN):: solver_option_path
    integer:: i
    
    call petsc_solve_state_setup(solver_option_path, prolongators, surface_nodes, &
      state, x%mesh, 1, x%option_path, has_solver_cache(matrix), option_path=option_path)

    if (associated(prolongators)) then
    
      if (associated(surface_nodes)) then
        
        call petsc_solve(x, matrix, rhs, &
           prolongators=prolongators, &
           surface_node_list=surface_nodes, option_path=option_path, &
           iterations_taken = iterations_taken)
           
        deallocate(surface_nodes)
        
      else
      
        call petsc_solve(x, matrix, rhs, &
           prolongators=prolongators, option_path=option_path, &
           iterations_taken = iterations_taken)
        
      end if
      
      do i=1, size(prolongators)
        call deallocate(prolongators(i))
      end do
      deallocate(prolongators)
    
    else
    
      call petsc_solve(x, matrix, rhs, option_path=option_path, &
                       iterations_taken = iterations_taken)
      
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
    
    integer, dimension(:), pointer:: surface_nodes
    type(petsc_csr_matrix), dimension(:), pointer:: prolongators
    character(len=OPTION_PATH_LEN):: solver_option_path
    integer:: i
    
    ! no solver cache for petsc_csr_matrices at the mo'
    call petsc_solve_state_setup(solver_option_path, prolongators, surface_nodes, &
      state, x%mesh, 1, x%option_path, .false., option_path=option_path)
    
    if (associated(prolongators)) then
    
      if (associated(surface_nodes)) then
        
        call petsc_solve(x, matrix, rhs, &
           prolongators=prolongators, &
           surface_node_list=surface_nodes, option_path=option_path)
           
        deallocate(surface_nodes)
        
      else
      
        call petsc_solve(x, matrix, rhs, &
           prolongators=prolongators, option_path=option_path)
        
      end if
      
      do i=1, size(prolongators)
        call deallocate(prolongators(i))
      end do
      deallocate(prolongators)
    
    else
    
      call petsc_solve(x, matrix, rhs, option_path=option_path)
      
    end if

  end subroutine petsc_solve_scalar_state_petsc_csr
  
  subroutine petsc_solve_vector_state_petsc_csr(x, matrix, rhs, state, &
    option_path)
    !!< Solve a linear system the nice way.
    !!< This version uses state to pull geometric information from
    !!< if required for the specified options.
    type(vector_field), intent(inout) :: x
    type(vector_field), intent(in) :: rhs
    type(petsc_csr_matrix), intent(inout) :: matrix
    type(state_type), intent(in):: state
    !! override x%option_path if provided:
    character(len=*), optional, intent(in):: option_path
    
    type(vector_field), pointer :: mesh_positions
    integer, dimension(:), pointer:: surface_nodes
    type(petsc_csr_matrix), dimension(:), pointer:: prolongators
    character(len=OPTION_PATH_LEN):: solver_option_path
    type(petsc_csr_matrix), pointer:: rotation_matrix
    integer:: i, rotation_stat
    
    ! no solver cache for petsc_csr_matrices at the mo'
    call petsc_solve_state_setup(solver_option_path, prolongators, surface_nodes, &
      state, x%mesh, x%dim, x%option_path, .false., option_path=option_path, &
      mesh_positions=mesh_positions)

    rotation_matrix => extract_petsc_csr_matrix(state, "RotationMatrix", stat=rotation_stat)
    if (rotation_stat==0 .and. associated(prolongators)) then
      FLExit("Rotated boundary conditions do not work with mg prolongators in the velocity solve")
    end if
    
    if (associated(prolongators) .and. associated(mesh_positions)) then
      call petsc_solve(x, matrix, rhs, &
           prolongators=prolongators, option_path=option_path, positions=mesh_positions)
    else if (associated(prolongators)) then
      call petsc_solve(x, matrix, rhs, &
           prolongators=prolongators, option_path=option_path)
    else if (associated(mesh_positions) .and. rotation_stat==0) then
      call petsc_solve(x, matrix, rhs, option_path=option_path, positions=mesh_positions, &
        rotation_matrix=rotation_matrix%M)
    else if (associated(mesh_positions)) then
      call petsc_solve(x, matrix, rhs, option_path=option_path, positions=mesh_positions)
    else if (rotation_stat==0) then
      call petsc_solve(x, matrix, rhs, option_path=option_path, rotation_matrix=rotation_matrix%M)
    else
      call petsc_solve(x, matrix, rhs, option_path=option_path)
    end if

    if (associated(mesh_positions)) then
      call deallocate(mesh_positions)
      deallocate(mesh_positions)
    end if

    if (associated(prolongators)) then
      do i=1, size(prolongators)
        call deallocate(prolongators(i))
      end do
      deallocate(prolongators)
    end if

  end subroutine petsc_solve_vector_state_petsc_csr
    
  subroutine petsc_solve_state_setup(solver_option_path, prolongators, surface_nodes, &
    state, mesh, field_dim, field_option_path, matrix_has_solver_cache, option_path, &
    mesh_positions)
    ! sets up monitors and returns solver_option_path,
    ! and prolongators and surface_nodes to be used in "mg" preconditioner
    character(len=*), intent(out):: solver_option_path
    ! if associated on return, this array of prolongators should be passed into petsc_solve
    type(petsc_csr_matrix), dimension(:), pointer:: prolongators
    ! if associated on return, this array of surface_nodes should be passed into petsc_solve
    integer, dimension(:), pointer:: surface_nodes

    type(state_type), intent(in):: state
    type(mesh_type), intent(in):: mesh ! mesh we're solving on
    integer, intent(in):: field_dim ! dimension of the field
    ! option_path of the provided field:
    character(len=*), intent(in):: field_option_path
    logical, intent(in):: matrix_has_solver_cache
    ! optional option_path that may be provided to override field option_path
    character(len=*), intent(in), optional:: option_path
    ! if associated on return, this mesh_positions field should be passed into petsc_solve
    ! currently only for vector solves
    type(vector_field), pointer, optional:: mesh_positions
    
    type(vector_field):: positions
    type(scalar_field), pointer:: exact
    type(mesh_type), pointer:: linear_mesh
    character(len=FIELD_NAME_LEN):: exact_field_name
    logical:: vertical_lumping, higher_order_lumping
    integer:: stat, no_prolongators
    
    if (present(option_path)) then
       solver_option_path=complete_solver_option_path(option_path)
    else
       solver_option_path=complete_solver_option_path(field_option_path)
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
       positions=get_nodal_coordinate_field(state, mesh)
       ! creates its own reference that's cleaned up in petsc_solve:
       call petsc_solve_monitor_iteration_vtus(positions)
       ! so we're free to get rid of ours
       call deallocate(positions)
    end if    
    
    nullify(prolongators)
    nullify(surface_nodes)
    
    higher_order_lumping = have_option(trim(solver_option_path)//'/preconditioner::mg/higher_order_lumping')
    vertical_lumping = have_option(trim(solver_option_path)//'/preconditioner::mg/vertical_lumping')
    no_prolongators = count( (/ higher_order_lumping, vertical_lumping /) )
    
    ! if the solver context has been cached from last time, we don't
    ! need to recreate the prolongation operators for "mg"
    if (no_prolongators>0 .and. .not. matrix_has_solver_cache) then
      allocate( prolongators(1:no_prolongators) )
      if (higher_order_lumping) then
        call find_linear_parent_mesh(state, mesh, linear_mesh)
        prolongators(1) = higher_order_prolongator(linear_mesh, mesh, field_dim)
      end if
      if (vertical_lumping) then
        if (field_dim>1) then
          FLExit("Cannot use vertical_lumping for vector fields")
        end if
        if (higher_order_lumping) then
          prolongators(2) = vertical_prolongator_from_free_surface(state, linear_mesh)
        else
          prolongators(1) = vertical_prolongator_from_free_surface(state, mesh)
        end if
        if (have_option(trim(solver_option_path)//'/preconditioner::mg/vertical_lumping/internal_smoother')) then
          surface_nodes => free_surface_nodes(state, mesh)
        end if
      end if
    end if

    if (petsc_solve_needs_positions(solver_option_path)) then
      if (.not. present(mesh_positions)) then
        ! currently this option only exists for vector solves, if it occurs in other places
        ! mesh_positions should be passed down
        FLAbort("mesh_positions should have been present")
      end if
      allocate(mesh_positions)
      ! get the positions of the nodes - for periodic this gives the aliased positions, is that right? who knows...
      mesh_positions = get_nodal_coordinate_field(state, mesh)
    else if (present(mesh_positions)) then
      nullify(mesh_positions)
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
        trim(solver_option_path)//'/diagnostics/monitors/iteration_vtus') &
      .or. petsc_solve_needs_positions(solver_option_path)
  
  end function petsc_solve_needs_state
  
  function higher_order_prolongator(p1_mesh, pn_mesh, ncomponents) result (P)
  ! Creates the linear operator that extrapolates p1 fields to higher
  ! order pn meshes. This can be used as the first stage prolongator
  ! in the "mg" multigrid preconditioner
    type(mesh_type), intent(in):: p1_mesh, pn_mesh
    integer, intent(in):: ncomponents
    type(petsc_csr_matrix):: P
    
    logical, dimension(:), allocatable:: nodes_visited
    integer, dimension(:), allocatable:: onnz, dnnz
    integer, dimension(:), pointer:: p1_nodes, pn_nodes
    integer:: rows, columns
    integer:: i, j, k, node, ele
    real, dimension(p1_mesh%shape%loc):: N
    
    rows=nowned_nodes(pn_mesh)
    columns=node_count(p1_mesh)
    allocate(dnnz(1:rows*ncomponents), onnz(1:rows*ncomponents))
    
    dnnz=0
    onnz=0
    do ele=1, ele_count(pn_mesh)
      pn_nodes => ele_nodes(pn_mesh, ele)
      p1_nodes => ele_nodes(p1_mesh, ele)
      do j=1, size(pn_nodes)
        node=pn_nodes(j)
        if (node_owned(pn_mesh, node)) then
          do k=1, size(p1_nodes)
            if (node_owned(p1_mesh, p1_nodes(k))) then
              dnnz(node)=dnnz(node)+1
            else
              onnz(node)=onnz(node)+1
            end if
          end do
        end if
      end do
    end do
    ! copy over to other components, if any
    do i=2, ncomponents
      dnnz( (i-1)*rows+1:i*rows )=dnnz(1:rows)
      onnz( (i-1)*rows+1:i*rows )=onnz(1:rows)
    end do
      
    call allocate(P, rows, columns, dnnz, onnz, (/ ncomponents, ncomponents /), name="HigherOrderProlongator")
    if (associated(P%column_halo)) then
      allocate(P%column_halo)
      P%column_halo = p1_mesh%halos(1)
      call incref(P%column_halo)
    end if
    call zero(P)
    
    allocate(nodes_visited(1:node_count(pn_mesh)))
    nodes_visited=.false.
    
    do ele=1, ele_count(pn_mesh)
      pn_nodes => ele_nodes(pn_mesh, ele)
      p1_nodes => ele_nodes(p1_mesh, ele)
      do j=1, size(pn_nodes)
        node=pn_nodes(j)
        if (node_owned(pn_mesh, node) .and. .not. nodes_visited(node)) then
          N = eval_shape(p1_mesh%shape, local_coords(j, pn_mesh%shape))
          do i=1, ncomponents
            do k=1, size(p1_nodes)
              call addto(P, i, i, node, p1_nodes(k), N(k))
            end do
          end do
          nodes_visited(node)=.true.
        end if
      end do
    end do
    
    call assemble(P)
    
  end function higher_order_prolongator
  
end module petsc_solve_state_module
