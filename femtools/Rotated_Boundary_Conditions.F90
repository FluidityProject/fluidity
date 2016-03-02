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

module rotated_boundary_conditions

use spud
use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN
use fields
use sparse_tools_petsc
use state_module
use halos
use boundary_conditions

implicit none

contains
  
  function have_rotated_bcs(u)

    type(vector_field), intent(in):: u
    logical:: have_rotated_bcs

    character(len=FIELD_NAME_LEN) :: bctype
    character(len=OPTION_PATH_LEN):: bc_option_path
    integer, dimension(:), pointer:: surface_node_list
    integer:: i

    do i=1, get_boundary_condition_count(u)
       call get_boundary_condition(u, i, type=bctype, &
            surface_node_list=surface_node_list, &
            option_path=bc_option_path)
       if (bctype=="dirichlet" .and. &
            have_option(trim(bc_option_path)//"/type[0]/align_bc_with_surface")) then
          have_rotated_bcs=.true.
          return
       end if
    end do

    have_rotated_bcs=.false.

  end function have_rotated_bcs
  
  subroutine create_rotation_matrix(rotation_m, u)

    type(petsc_csr_matrix), intent(out):: rotation_m
    type(vector_field), intent(in):: u

    type(halo_type), pointer:: halo
    type(vector_field), pointer:: normal, tangent1, tangent2
    character(len=FIELD_NAME_LEN):: bctype
    character(len=OPTION_PATH_LEN):: bc_option_path
    integer, dimension(:), pointer:: surface_node_list
    real, dimension(u%dim, u%dim):: local_rotation
    integer, dimension(:), allocatable:: dnnz, onnz
    integer:: i, j, node, nodes, mynodes
    logical:: parallel

    ewrite(1,*) "Inside create_rotation_matrix"

    nodes=node_count(u)
    if (associated(u%mesh%halos)) then
       halo => u%mesh%halos(1)
       mynodes=halo_nowned_nodes(halo)
    else
       nullify(halo)
       mynodes=nodes
    end if
    parallel=IsParallel()

    allocate(dnnz(1:mynodes*u%dim), onnz(1:mynodes*u%dim))
    onnz=0
    ! default is just a 1.0 on the diagonal (no rotation)
    dnnz=1

    do i=1, get_boundary_condition_count(u)
       call get_boundary_condition(u, i, type=bctype, &
            surface_node_list=surface_node_list, &
            option_path=bc_option_path)
       if (bctype=="dirichlet" .and. &
            have_option(trim(bc_option_path)//"/type[0]/align_bc_with_surface")) then

          do j=1, size(surface_node_list)
             node=surface_node_list(j)
             if (parallel) then
                if (node>mynodes) cycle
             endif
             if (any(dnnz(node:node+(u%dim-1)*mynodes:mynodes)>1)) then
               FLExit("Two rotated boundary condition specifications for the same node.")
             end if
             dnnz( node:node+(u%dim-1)*mynodes:mynodes)=u%dim
          end do

       end if
    end do

    call allocate(rotation_m, nodes, nodes, &
         dnnz, onnz, (/ u%dim, u%dim /), "RotationMatrix", halo=halo)

    ! put a 1.0 on the diagonal for non-rotated nodes
    do i=1, mynodes
       ! skip rotated nodes
       if (dnnz(i)/=1) cycle

       do j=1, u%dim
          call addto(rotation_m, j, j, i, i, 1.0)
       end do
    end do

    ! insert the local rotation matrix as a diagonal block for the rotated nodes
    do i=1, get_boundary_condition_count(u)
       call get_boundary_condition(u, i, type=bctype, &
            surface_node_list=surface_node_list, &
            option_path=bc_option_path)

       if (bctype=="dirichlet" .and. &
            have_option(trim(bc_option_path)//"/type[0]/align_bc_with_surface")) then

          normal   => extract_surface_field(u, i, "normal")
          tangent1 => extract_surface_field(u, i, "tangent1")
          tangent2 => extract_surface_field(u, i, "tangent2")

          do j=1, size(surface_node_list)
             node=surface_node_list(j)
             if (node > mynodes) cycle
             local_rotation(:,1)=node_val(normal, j)
             local_rotation(:,2)=node_val(tangent1, j)
             if (u%dim>2) then
                local_rotation(:,3)=node_val(tangent2, j)
             end if

             call addto(rotation_m, node, node, local_rotation)
          end do

       end if
    end do

    call assemble(rotation_m)

  end subroutine create_rotation_matrix
    
  subroutine rotate_momentum_equation(big_m, rhs, u, state, dg)

    type(petsc_csr_matrix), intent(inout):: big_m
    type(vector_field), intent(inout):: rhs
    type(vector_field), intent(inout):: u
    type(state_type), intent(inout):: state
    logical, intent(in) :: dg

    type(petsc_csr_matrix), pointer:: rotation_m
    type(petsc_csr_matrix):: rotated_big_m
    type(vector_field):: result
    integer:: stat

    ewrite(1,*) "Inside rotate_momentum_equation"

    rotation_m => extract_petsc_csr_matrix(state, "RotationMatrix", stat=stat)

    if (stat/=0) then
       allocate(rotation_m)
       call create_rotation_matrix(rotation_m, u)
       call insert(state, rotation_m, "RotationMatrix")
    end if

    !call assemble(big_m)
    !call dump_matrix("bigm", big_m)

    ! rotate big_m:
    call ptap(rotated_big_m, big_m, rotation_m)

    !call dump_matrix("rotated_bigm", rotated_big_m)

    ! rotate rhs:
    ! need to have separate copy of the field, because of intent(out) and intent(in)
    ! of mult_T call, as result%val points at the same space as rhs%val, this directly
    ! puts the result in rhs as well 
    result=rhs 
    call mult_T(result, rotation_m, rhs)
    if (dg) then
      ! We have just poluted the halo rows of the rhs. This is incorrect
      ! in the dg case due to the non-local assembly system employed.
      call zero_non_owned(rhs)
    end if
    ! rotate u:
    if (dg) then
      call zero_non_owned(u)
    end if
    result=u ! same story
    call mult_T(result, rotation_m, u)

    ! throw out unrotated big_m and replace with rotated:
    call deallocate(big_m)
    big_m=rotated_big_m

    if (stat/=0) then
      call deallocate(rotation_m)
      deallocate(rotation_m)
    end if

  end subroutine rotate_momentum_equation
    
  subroutine rotate_ct_m(ct_m, u)

    type(block_csr_matrix), intent(inout):: ct_m
    type(vector_field), intent(in):: u

    type(vector_field), pointer:: normal, tangent1, tangent2
    character(len=FIELD_NAME_LEN):: bctype
    character(len=OPTION_PATH_LEN):: bc_option_path
    integer, dimension(:), pointer:: surface_node_list, rowcol
    integer, dimension(:), allocatable:: node2rotated_node
    real, dimension(u%dim, u%dim):: local_rotation
    real, dimension(u%dim):: ct_xyz, ct_rot
    real, dimension(:), pointer:: rowval
    integer:: bc, i, j, k, rotated_node

    ewrite(1,*) "Inside rotate_ct_m"

    assert( all(blocks(ct_m) == (/ 1, u%dim /)) )

    allocate( node2rotated_node(1:node_count(u)) )

    do bc=1, get_boundary_condition_count(u)
       call get_boundary_condition(u, bc, type=bctype, &
            surface_node_list=surface_node_list, &
            option_path=bc_option_path)
       if (bctype=="dirichlet" .and. &
            have_option(trim(bc_option_path)//"/type[0]/align_bc_with_surface")) then

          normal => extract_surface_field(u, bc, "normal")
          tangent1 => extract_surface_field(u, bc, "tangent1")
          tangent2 => extract_surface_field(u, bc, "tangent2")

          node2rotated_node=0
          node2rotated_node(surface_node_list)=(/ (j, j=1, size(surface_node_list)) /)

          do i=1, size(ct_m, 1)
             rowcol => row_m_ptr(ct_m, i)
             do j=1, size(rowcol)
                rotated_node=node2rotated_node(rowcol(j))
                if (rotated_node/=0) then
                   ! construct local rotation matrix
                   local_rotation(1,:)=node_val(normal, rotated_node)
                   local_rotation(2,:)=node_val(tangent1, rotated_node)
                   if (u%dim>2) then
                      local_rotation(3,:)=node_val(tangent2, rotated_node)
                   end if

                   ! look up ct_m values of row i, column rowcol(j) in xyz orientation
                   do k=1, blocks(ct_m,2)
                      rowval => row_val_ptr(ct_m, 1, k, i)
                      ct_xyz(k)=rowval(j)
                   end do
                   ! rotate to normal, tangent1, tangent2 orientation
                   ct_rot=matmul( local_rotation, ct_xyz)
                   ! put back in the matrix
                   do k=1, blocks(ct_m,2)
                      rowval => row_val_ptr(ct_m, 1, k, i)
                      rowval(j)=ct_rot(k)
                   end do
                end if
             end do
          end do
       end if
    end do

    deallocate(node2rotated_node)

  end subroutine rotate_ct_m
  
  subroutine rotate_velocity(vfield, state)

    type(vector_field), intent(inout):: vfield
    type(state_type), intent(inout):: state
    
    type(vector_field), pointer:: u
    type(vector_field):: result
    type(petsc_csr_matrix), pointer:: rotation_m
    integer:: stat
    
    rotation_m => extract_petsc_csr_matrix(state, "RotationMatrix", stat=stat)
    if (stat/=0) then
       allocate(rotation_m)
       ! the vector field we are rotating might not have the bcs attached to it:
       u => extract_vector_field(state, "Velocity")
       call create_rotation_matrix(rotation_m, u)
       call insert(state, rotation_m, "RotationMatrix")
    end if
    
    result=vfield ! see note in rotate_momentum_equation
    call mult_T(result, rotation_m, vfield)

    if (stat/=0) then
      call deallocate(rotation_m)
      deallocate(rotation_m)
    end if

  end subroutine rotate_velocity
  
  subroutine rotate_velocity_back(vfield, state)

    type(vector_field), intent(inout):: vfield
    type(state_type), intent(inout):: state
    
    type(vector_field), pointer:: u
    type(vector_field):: result
    type(petsc_csr_matrix), pointer:: rotation_m
    integer:: stat
    
    rotation_m => extract_petsc_csr_matrix(state, "RotationMatrix", stat=stat)
    if (stat/=0) then
       allocate(rotation_m)
       ! the vector field we are rotating might not have the bcs attached to it:
       u => extract_vector_field(state, "Velocity")
       call create_rotation_matrix(rotation_m, u)
       call insert(state, rotation_m, "RotationMatrix")
    end if
    
    result=vfield  ! see note in rotate_momentum_equation
    call mult(result, rotation_m, vfield)

    if (stat/=0) then
      call deallocate(rotation_m)
      deallocate(rotation_m)
    end if

  end subroutine rotate_velocity_back
  
end module rotated_boundary_conditions
