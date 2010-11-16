!    Copyright (C) 2006 Imperial College London and others.
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
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
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
!
module mmpde_module
!
  use elements
  use solvers
  use global_numbering
  use shape_functions
  use fields
  use element_numbering
  use field_derivatives
  use state_module
  use FEtools
  use sparsity_patterns
  use global_parameters
  use FLDebug
  use sparse_tools
  use vtk_interfaces
  use unittest_tools
  use spud
!
implicit none
!
private
!
public :: mmpde
!
contains
!
! A mesh movement algorithm
!
  subroutine mmpde(state)
    implicit none
    type(state_type), intent(in) :: state
    type(scalar_field), pointer :: monitor
    type(vector_field), pointer :: positions, u
    type(element_type) x_shape
    real, dimension(:), allocatable :: detwei,mass_lump,   xorig,yorig,zorig
    real, dimension(:,:), allocatable :: mass_matrix
    real, dimension(:,:,:), allocatable :: dn_t, Diffusivity_gi
    type(csr_sparsity) :: sparsity
    type(csr_matrix) :: M,A_x,A_y,A_z
    type(scalar_field) :: masslump,RHS_x,RHS_y,RHS_z,x_old,x_new,x_original,y_old,y_new,y_original,z_old,z_new,z_original
    type(tensor_field) :: Diffusivity

    integer :: ele,node,dim,meshdim,loc,nloc,snloc,ngi,sngi,number_of_iterations
    integer :: mesh_number,time_level
    real :: current_time, dt, ds, theta
    character(len=256) :: field_name
    integer, save :: dump = 0
       
    ewrite(3,*)'In mmpde'
    field_name = "Monitor"
!
    positions => extract_vector_field(state, "Coordinate")      
    u         => extract_vector_field(state, "Velocity")      
    monitor   => extract_scalar_field(state, "Monitor")
!
    meshdim   = mesh_dim(positions)      
    x_shape   = ele_shape(positions, 1)
    nloc      = ele_loc(positions, 1)
    snloc     = face_loc(positions, 1)
    ngi       = ele_ngi(positions, 1)
    sngi      = face_ngi(positions, 1)
!      
    sparsity = make_sparsity(monitor%mesh, monitor%mesh, name='Sparsity')
    call allocate(A_x, sparsity, name = 'LHS_Matrix_x')    
    call zero(A_x)
    call allocate(A_y, sparsity, name = 'LHS_Matrix_x')    
    call zero(A_y)
    call allocate(A_z, sparsity, name = 'LHS_Matrix_x')    
    call zero(A_z)
    call allocate(RHS_x, monitor%mesh, name='RHS_x')
    call zero(RHS_x)
    call allocate(RHS_y, monitor%mesh, name='RHS_y')
    call zero(RHS_y)
    call allocate(RHS_z, monitor%mesh, name='RHS_y')
    call zero(RHS_z)
    call allocate(M, sparsity, name='MassMatrix')
    call zero(M)
    call allocate(masslump,monitor%mesh)
    call zero(masslump)



    call get_option('/timestepping/current_time', current_time)
    call get_option('/timestepping/timestep', dt)
    ds = 0.0001  ! time step for our pseudo time stepping routine

    call get_option(trim(monitor%option_path)//"/prognostic/temporal_discretisation/theta", theta)

! At the moment the Diffusivity is the monitor function (or one over)
    Diffusivity = extract_tensor_field(state, trim(field_name)//"Diffusivity")

    allocate( Diffusivity_gi(meshdim,meshdim,ngi) )

    x_original = extract_scalar_field(positions,1) !original x locations
    y_original = extract_scalar_field(positions,2) !original y locations
    z_original = extract_scalar_field(positions,3)   
      
    call allocate(x_old, positions%mesh, "x_old")
    call allocate(x_new, positions%mesh, "x_new")
    call allocate(y_old, positions%mesh, "y_old")
    call allocate(y_new, positions%mesh, "y_new")
    call allocate(z_old, positions%mesh, "z_old")
    call allocate(z_new, positions%mesh, "z_new")
    
    x_old%val(:) = x_original%val(:)
    x_new%val(:) = x_original%val(:)
    y_old%val(:) = y_original%val(:)
    y_new%val(:) = y_original%val(:)
    z_old%val(:) = z_original%val(:)
    z_new%val(:) = z_original%val(:)

    allocate(xorig(node_count(positions)))
    allocate(yorig(node_count(positions)))
    allocate(zorig(node_count(positions)))
    
    do node = 1,node_count(positions)
       xorig(node) = positions%val(1,node) 
       yorig(node) = positions%val(2,node)
       zorig(node) = positions%val(3,node)
    end do

    call vtk_write_fields(trim("Mesh"), 0, positions, positions%mesh, sfields=(/x_new,y_new,z_new/), tfields=(/Diffusivity/) )    
    
    mesh_number = 0
    
    call vtk_write_fields(trim("Work"), mesh_number, positions, positions%mesh, sfields=(/x_new,y_new,z_new/) )    

   
    do time_level = 1,4
      element_loop: do ele=1, element_count(positions)      
        Diffusivity_gi = ele_val_at_quad(Diffusivity, ele)       
        call assemble_element_contribution_mmpde(A_x, RHS_x, positions, x_old, monitor,&
                                               ele, ds, theta, Diffusivity_gi)
        call assemble_element_contribution_mmpde(A_y, RHS_y, positions, y_old, monitor,&
                                               ele, ds, theta, Diffusivity_gi)
        call assemble_element_contribution_mmpde(A_z, RHS_z, positions, z_old, monitor,&
                                               ele, ds, theta, Diffusivity_gi)
      end do element_loop
! Add some boundary conditions which at the moment just fix nodes on border of square/box
      call mmpde_bc(A_x,RHS_x,positions,x_original)
      call mmpde_bc(A_y,RHS_y,positions,y_original)
      call mmpde_bc(A_z,RHS_z,positions,z_original)
      
! Solve linear systems using options taken from monitor field
      call zero(x_new)
      call zero(y_new)
      call zero(z_new)
      call petsc_solve(x_new, A_x, RHS_x, monitor%option_path)
      call petsc_solve(y_new, A_y, RHS_y, monitor%option_path)
      call petsc_solve(z_new, A_z, RHS_z, monitor%option_path)

      x_old%val(:) = x_new%val(:)
      y_old%val(:) = y_new%val(:)
      z_old%val(:) = z_new%val(:)
          
      mesh_number = mesh_number + 1
      call vtk_write_fields(trim("Work"), mesh_number, positions, positions%mesh, sfields=(/x_new,y_new,z_new,monitor/))            
    enddo


! Update the locations
    do node = 1,node_count(positions)
       positions%val(1,node) = x_new%val(node)
       positions%val(2,node) = y_new%val(node)
       positions%val(3,node) = z_new%val(node)
    end do
! Output the new mesh
    dump = dump+1
    call vtk_write_fields(trim("Mesh"), dump, positions, positions%mesh, sfields=(/x_new,y_new,z_new/), tfields=(/Diffusivity/) )    

! Now put them back for next time step - as my projection currently assumes base is comp coord.
    do node = 1,node_count(positions)
       positions%val(1,node) = xorig(node)
       positions%val(2,node) = yorig(node)
       positions%val(3,node) = zorig(node)
    end do

    deallocate(xorig)
    deallocate(yorig)
    deallocate(zorig)
            
!    stop      
  end subroutine mmpde
!
!
  subroutine assemble_element_contribution_mmpde(A, RHS, positions, xyz, monitor,&
                                                 ele, dt, theta, Diffusivity_gi)
    type(csr_matrix),   intent(inout) :: A
    type(scalar_field), intent(inout) :: RHS
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: xyz, monitor
    integer, intent(in) :: ele
    real, intent(in) :: dt, theta

    ! Locations of nodes.
    real, dimension(positions%dim,ele_loc(positions,ele)) :: X_ele
    real, dimension(ele_loc(positions,ele)) :: xyz_ele
    ! Locations of quadrature points
    real, dimension(positions%dim,ele_ngi(positions,ele)) :: X_quad
    ! Derivatives of shape function:
    real, dimension(ele_loc(monitor,ele), ele_ngi(monitor,ele), positions%dim) :: dshape_monitor
    ! Coordinate transform * quadrature weights.
    real, dimension(ele_ngi(positions,ele)) :: detwei    
    ! Node numbers of monitor element.
    integer, dimension(:), pointer :: ele_monitor
    ! Shape functions.
    type(element_type), pointer :: shape_monitor, shape_X
    ! Local discretisation matrix 
    real, dimension(ele_loc(monitor, ele), ele_loc(monitor, ele)) :: monitor_mat, mass_matrix
    ! Local right hand side.
    real, dimension(ele_loc(monitor, ele)) :: lrhs
    real, dimension(positions%dim,positions%dim,ele_ngi(positions,ele)) :: Diffusivity_gi
    integer :: loc
!
!
    ele_monitor   => ele_nodes(monitor, ele)
    shape_monitor => ele_shape(monitor, ele)
    shape_X       => ele_shape(positions, ele)

    ! Locations of local vertices.
    X_ele = ele_val(positions, ele)

    xyz_ele = ele_val(xyz, ele)

    ! Locations of quadrature points.
    X_quad = ele_val_at_quad(positions, ele)

    ! Transform derivatives and weights into physical space.
    call transform_to_physical(X_ele, shape_X, shape_X, dm_t=dshape_monitor, detwei=detwei)

    mass_matrix = shape_shape(shape_X, shape_X, detwei)    

    ! Local assembly:
    monitor_mat = dshape_tensor_dshape(dshape_monitor, Diffusivity_gi, dshape_monitor, detwei)

    ! Global assembly:
    call addto(A, ele_monitor, ele_monitor, mass_matrix + dt*theta*monitor_mat)

    call addto(RHS, ele_monitor, matmul(mass_matrix - dt*(1.0-theta)*monitor_mat, xyz_ele) )

  end subroutine assemble_element_contribution_mmpde
!
!
  subroutine mmpde_bc(A,RHS,positions,x_original)
    type(csr_matrix),   intent(inout) :: A
    type(scalar_field), intent(inout) :: RHS
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: x_original
    integer :: meshdim,nodes,nod
    real :: big,small
!    
    meshdim   = mesh_dim(positions)      
    nodes     = node_count(positions)
    big = 1.0e+40
    small = 1.0e-6
!    
! Add some BCs specific to a simple box/square domain
    do nod=1,nodes
       if( (abs(x_original%val(nod)-0.0)<small).or.&
           (abs(x_original%val(nod)-1.0)<small) ) then
         call addto_diag(A, nod, big)
         call addto(RHS, nod, big*x_original%val(nod))
       endif
    end do  
 !       
  end subroutine mmpde_bc
!
!
  end module mmpde_module
!
!
!
