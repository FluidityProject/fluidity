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
!    C.Pain@Imperial.ac.uk
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

module field_copies_diagnostics

  use diagnostic_source_fields
  use field_options
  use fields
  use fldebug
  use solvers
  use sparse_tools
  use sparsity_patterns_meshes
  use state_module
  use transform_elements

  implicit none
  
  private
  
  public :: calculate_scalar_copy, calculate_vector_copy, calculate_tensor_copy
  public :: calculate_scalar_galerkin_projection, calculate_vector_galerkin_projection
  public :: extract_scalar_component_copy
contains

  subroutine calculate_scalar_copy(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(state, s_field)
    
    call remap_field(source_field, s_field)
    
  end subroutine calculate_scalar_copy
  
  subroutine calculate_vector_copy(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    
    type(vector_field), pointer :: source_field
    
    source_field => vector_source_field(state, v_field)
    
    call remap_field(source_field, v_field)
    
  end subroutine calculate_vector_copy
  
  subroutine calculate_tensor_copy(state, t_field)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(inout) :: t_field
    
    type(tensor_field), pointer :: source_field
    
    source_field => tensor_source_field(state, t_field)
    
    call remap_field(source_field, t_field)
    
  end subroutine calculate_tensor_copy

  subroutine extract_scalar_component_copy(state, s_field, dim)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    integer, intent(in) :: dim

    type(vector_field), pointer :: source_field
    type(scalar_field) :: component

    source_field => vector_source_field(state, s_field)
    component = extract_scalar_field(source_field, dim)

    call remap_field(component, s_field)
    
  end subroutine extract_scalar_component_copy
  
  subroutine calculate_scalar_galerkin_projection(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: positions
    
    source_field => scalar_source_field(state, s_field)
    positions => extract_vector_field(state, "Coordinate")
        
    select case(continuity(s_field))
      case(0)
        call gp_continuous
      case(-1)
        call gp_discontinuous
      case default
        ewrite(-1, *) "For mesh continuity", continuity(s_field)
        FLAbort("Unrecognised mesh continuity")
    end select
    
  contains
  
    subroutine gp_continuous
      type(csr_matrix) :: mass
      type(csr_sparsity), pointer :: sparsity
      type(scalar_field) :: rhs
      
      integer :: i
      
      sparsity => get_csr_sparsity_firstorder(state, s_field%mesh, s_field%mesh)
      call allocate(mass, sparsity, name = "MassMatrix")
      call allocate(rhs, s_field%mesh, "RHS")
      
      call zero(mass)
      call zero(rhs)
      do i = 1, ele_count(s_field)
        call assemble_gp_ele(i, s_field, positions, source_field, mass, rhs)
      end do
      
      call petsc_solve(s_field, mass, rhs, &
        & option_path = trim(complete_field_path(s_field%option_path)) // "/algorithm")
      
      call deallocate(mass)
      call deallocate(rhs)
      
    end subroutine gp_continuous
    
    subroutine assemble_gp_ele(ele, s_field, positions, source_field, mass, rhs)
      integer, intent(in) :: ele
      type(scalar_field), intent(in) :: s_field
      type(vector_field), intent(in) :: positions
      type(scalar_field), intent(in) :: source_field
      type(csr_matrix), intent(inout) :: mass
      type(scalar_field), intent(inout) :: rhs
      
      integer, dimension(:), pointer :: nodes
      real, dimension(ele_ngi(s_field, ele)) :: detwei
      type(element_type), pointer :: shape
      
      shape => ele_shape(s_field, ele)
      
      call transform_to_physical(positions, ele, detwei = detwei)
      
      nodes => ele_nodes(s_field, ele)
      
      call addto(mass, nodes, nodes, shape_shape(shape, shape, detwei))
      call addto(rhs, nodes, shape_rhs(shape, detwei * ele_val_at_quad(source_field, ele)))
    
    end subroutine assemble_gp_ele
    
    subroutine gp_discontinuous
      integer :: i
      
      do i = 1, ele_count(s_field)
        call solve_gp_ele(i, s_field, positions, source_field)
      end do
      
    end subroutine gp_discontinuous
    
    subroutine solve_gp_ele(ele, s_field, positions, source_field)
      integer, intent(in) :: ele
      type(scalar_field), intent(inout) :: s_field
      type(vector_field), intent(in) :: positions
      type(scalar_field), intent(in) :: source_field
      
      integer, dimension(:), pointer :: nodes
      real, dimension(ele_loc(s_field, ele)) :: little_rhs
      real, dimension(ele_loc(s_field, ele), ele_loc(s_field, ele)) :: little_mass
      real, dimension(ele_ngi(s_field, ele)) :: detwei
      type(element_type), pointer :: shape
      
      shape => ele_shape(s_field, ele)
      
      call transform_to_physical(positions, ele, detwei = detwei)
            
      little_mass = shape_shape(shape, shape, detwei)
      little_rhs = shape_rhs(shape, detwei * ele_val_at_quad(source_field, ele))      
      call solve(little_mass, little_rhs)
      
      nodes => ele_nodes(s_field, ele)
      
      call set(s_field, nodes, little_rhs)
    
    end subroutine solve_gp_ele
    
  end subroutine calculate_scalar_galerkin_projection
  
  subroutine calculate_vector_galerkin_projection(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    
    type(vector_field), pointer :: source_field
    type(vector_field), pointer :: positions
    
    source_field => vector_source_field(state, v_field)
    positions => extract_vector_field(state, "Coordinate")
        
    select case(continuity(v_field))
      case(0)
        call gp_continuous
      case(-1)
        call gp_discontinuous
      case default
        ewrite(-1, *) "For mesh continuity", continuity(v_field)
        FLAbort("Unrecognised mesh continuity")
    end select
    
  contains
  
    subroutine gp_continuous
      type(csr_matrix) :: mass
      type(csr_sparsity), pointer :: sparsity
      type(vector_field) :: rhs
      
      integer :: i
      
      sparsity => get_csr_sparsity_firstorder(state, v_field%mesh, v_field%mesh)
      call allocate(mass, sparsity, name = "MassMatrix")
      call allocate(rhs, v_field%dim, v_field%mesh, "RHS")
      
      call zero(mass)
      call zero(rhs)
      do i = 1, ele_count(v_field)
        call assemble_gp_ele(i, v_field, positions, source_field, mass, rhs)
      end do
      
      call petsc_solve(v_field, mass, rhs, &
        & option_path = trim(complete_field_path(v_field%option_path)) // "/algorithm")
      
      call deallocate(mass)
      call deallocate(rhs)
      
    end subroutine gp_continuous
    
    subroutine assemble_gp_ele(ele, v_field, positions, source_field, mass, rhs)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: v_field
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: source_field
      type(csr_matrix), intent(inout) :: mass
      type(vector_field), intent(inout) :: rhs
      
      integer, dimension(:), pointer :: nodes
      real, dimension(ele_ngi(v_field, ele)) :: detwei
      type(element_type), pointer :: shape
      
      shape => ele_shape(v_field, ele)
      
      call transform_to_physical(positions, ele, detwei = detwei)
      
      nodes => ele_nodes(v_field, ele)
      
      call addto(mass, nodes, nodes, shape_shape(shape, shape, detwei))
      call addto(rhs, nodes, shape_vector_rhs(shape, ele_val_at_quad(source_field, ele), detwei))
    
    end subroutine assemble_gp_ele
    
    subroutine gp_discontinuous
      integer :: i
      
      do i = 1, ele_count(v_field)
        call solve_gp_ele(i, v_field, positions, source_field)
      end do
      
    end subroutine gp_discontinuous
    
    subroutine solve_gp_ele(ele, v_field, positions, source_field)
      integer, intent(in) :: ele
      type(vector_field), intent(inout) :: v_field
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: source_field
      
      integer, dimension(:), pointer :: nodes
      real, dimension(ele_loc(v_field, ele), source_field%dim) :: little_rhs
      real, dimension(ele_loc(v_field, ele), ele_loc(v_field, ele)) :: little_mass
      real, dimension(ele_ngi(v_field, ele)) :: detwei
      type(element_type), pointer :: shape
      
      shape => ele_shape(v_field, ele)
      
      call transform_to_physical(positions, ele, detwei = detwei)
            
      little_mass = shape_shape(shape, shape, detwei)
      little_rhs = transpose(shape_vector_rhs(shape, ele_val_at_quad(source_field, ele), detwei)  )
      call solve(little_mass, little_rhs)
      
      nodes => ele_nodes(v_field, ele)
      
      call set(v_field, nodes, transpose(little_rhs))
    
    end subroutine solve_gp_ele
    
  end subroutine calculate_vector_galerkin_projection

end module field_copies_diagnostics
