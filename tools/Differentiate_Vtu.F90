!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineeringp
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

subroutine differentiate_vtu(input_filename, input_filename_len, output_filename, output_filename_len, input_fieldname, input_fieldname_len)

! these 5 need to be on top and in this order, so as not to confuse silly old intel compiler 
  use quadrature
  use elements
  use sparse_tools
  use fields
  use state_module
!
  use field_derivatives
  use fldebug
  use state_module
  use vtk_interfaces
  
  implicit none
  
  integer, intent(in) :: input_filename_len
  integer, intent(in) :: output_filename_len
  integer, intent(in) :: input_fieldname_len
  character(len = input_filename_len), intent(in) :: input_filename
  character(len = output_filename_len), intent(in) :: output_filename
  character(len = input_fieldname_len), intent(in) :: input_fieldname
  
  integer :: dim, i, j, nfields
  logical :: allocated
  type(mesh_type), pointer :: mesh
  type(scalar_field) :: masslump
  type(scalar_field), dimension(:), allocatable :: s_fields
  type(scalar_field), pointer :: s_field
  type(state_type) :: collapsed_state, state
  type(vector_field) :: field_grad
  type(vector_field), dimension(:), allocatable :: field_grads
  type(vector_field), pointer :: positions
 
  ewrite(1, *) "In differentiate_vtu"
  
  call vtk_read_state(trim(input_filename), state)
  
  positions => extract_vector_field(state, "Coordinate")
  dim = positions%dim
  mesh => extract_mesh(state, "Mesh")
  
  if(len_trim(input_fieldname) == 0) then
    call collapse_fields_in_state(state, collapsed_state)
    nfields = scalar_field_count(collapsed_state)
    allocate(s_fields(nfields))
    do i = 1, scalar_field_count(collapsed_state)
      s_fields(i) = extract_scalar_field(collapsed_state, i)
    end do
    allocate(field_grads(nfields))
    do i = 1, nfields
      s_field => extract_scalar_field(collapsed_state, i)
      call allocate(field_grads(i), dim, mesh, trim(s_field%name) // "Gradient")
    end do
    call deallocate(collapsed_state)
    
    select case(continuity(positions))
      case(-1)
        do i = 1, ele_count(mesh)
          call solve_grad_ele(i, positions, s_fields, field_grads)
        end do
      case(0)
        call allocate(masslump, mesh, "LumpedMass")
        call zero(masslump)
        do i = 1, nfields
          call zero(field_grads(i))
        end do
        
        do i = 1, ele_count(mesh)
          call assemble_grad_ele(i, positions, s_fields, masslump, field_grads)
        end do
        
        do i = 1, nfields
          do j = 1, dim
            field_grads(i)%val(j)%ptr = field_grads(i)%val(j)%ptr / masslump%val
          end do
        end do
        call deallocate(masslump)
      case default
        ewrite(-1, *) "For continuity ", continuity(positions)
        FLAbort("Unrecognised continuity")
    end select
    
    do i = 1, nfields
      call insert(state, field_grads(i), field_grads(i)%name)
      call deallocate(field_grads(i))
    end do
    deallocate(field_grads)
    deallocate(s_fields)
  else
    s_field => extract_scalar_field(state, trim(input_fieldname), allocated = allocated)
  
    call allocate(field_grad, dim, mesh, trim(s_field%name) // "Gradient")
    call grad(s_field, positions, field_grad)
    call insert(state, field_grad, field_grad%name) 
    call deallocate(field_grad)   
    
    if(allocated) deallocate(s_field)
  end if
  
  call vtk_write_state(output_filename, state = (/state/))
  call deallocate(state)
  
  call print_references(0)
  
  ewrite(1, *) "Exiting differentate_vtu"
 
contains

  subroutine solve_grad_ele(ele, positions, fields, field_grads)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), dimension(:), intent(in) :: fields
    type(vector_field), dimension(size(fields)), intent(inout) :: field_grads
    
    integer :: i
    integer, dimension(:), pointer :: nodes
    real, dimension(ele_ngi(positions, ele)) :: detwei
    real, dimension(ele_loc(positions, ele), ele_loc(positions, ele)) :: little_mass
    real, dimension(ele_loc(positions, ele), size(fields) * positions%dim) :: little_rhs
    real, dimension(ele_loc(positions, ele), ele_ngi(positions, ele), positions%dim) :: dshape
    type(element_type), pointer :: shape
    
    shape => ele_shape(positions, ele)
    call transform_to_physical(positions, ele, shape, &
      & detwei = detwei, dshape = dshape)
      
    little_mass = shape_shape(shape, shape, detwei)
    do i = 1, size(fields)
      little_rhs(:, (i - 1) * positions%dim + 1:i * positions%dim) = transpose(shape_vector_rhs(shape, transpose(ele_grad_at_quad(fields(i), ele, dshape)), detwei))
    end do
    call solve(little_mass, little_rhs)
      
    nodes => ele_nodes(positions, ele)
    do i = 1, size(field_grads)
      call set(field_grads(i), nodes, little_rhs(:, (i - 1) * positions%dim + 1:i * positions%dim))
    end do
    
  end subroutine solve_grad_ele

  subroutine assemble_grad_ele(ele, positions, fields, masslump, rhs)
    integer, intent(in) :: ele
    type(vector_field), intent(in) :: positions
    type(scalar_field), dimension(:), intent(in) :: fields
    type(scalar_field), intent(inout) :: masslump
    type(vector_field), dimension(size(fields)), intent(inout) :: rhs
    
    integer :: i
    integer, dimension(:), pointer :: nodes
    real, dimension(ele_ngi(positions, ele)) :: detwei
    real, dimension(ele_loc(positions, ele), ele_ngi(positions, ele), positions%dim) :: dshape
    type(element_type), pointer :: shape
    
    shape => ele_shape(positions, ele)
    call transform_to_physical(positions, ele, shape, &
      & detwei = detwei, dshape = dshape)
      
    nodes => ele_nodes(positions, ele)
    call addto(masslump, nodes, sum(shape_shape(shape, shape, detwei), 2))
    do i = 1, size(rhs)
      call addto(rhs(i), nodes, shape_vector_rhs(shape, transpose(ele_grad_at_quad(fields(i), ele, dshape)), detwei))
    end do
      
  end subroutine assemble_grad_ele
       
end subroutine differentiate_vtu
