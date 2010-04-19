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

module differential_operator_diagnostics

  use diagnostic_source_fields
  use divergence_matrix_cg
  use field_derivatives
  use field_options
  use fields
  use fldebug
  use geostrophic_pressure
  use global_parameters, only : OPTION_PATH_LEN
  use solvers
  use sparse_matrices_fields
  use sparsity_patterns_meshes
  use spud
  use state_fields_module
  use state_module

  implicit none
  
  private
  
  public :: calculate_grad, calculate_div, calculate_curl, calculate_curl_2d, &
    & calculate_finite_element_divergence, &
    & calculate_finite_element_divergence_transpose, &
    & calculate_scalar_advection, calculate_scalar_laplacian, &
    & calculate_vector_advection, calculate_vector_laplacian
  
contains

  subroutine calculate_grad(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: positions
    
    source_field => scalar_source_field(state, v_field)
    
    positions => extract_vector_field(state, "Coordinate")
    
    call grad(source_field, positions, v_field)
    
  end subroutine calculate_grad

  subroutine calculate_div(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    type(vector_field), pointer :: positions, source_field
    
    source_field => vector_source_field(state, s_field)
    
    positions => extract_vector_field(state, "Coordinate")
    
    call div(source_field, positions, s_field)
    
  end subroutine calculate_div
  
  subroutine calculate_finite_element_divergence(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    
    character(len = OPTION_PATH_LEN) :: path
    type(block_csr_matrix) :: CT_m
    type(csr_sparsity), pointer :: divergence_sparsity
    type(csr_matrix), pointer :: mass
    type(scalar_field) :: ctfield, ct_rhs
    type(scalar_field), pointer :: masslump
    type(vector_field), pointer :: positions, source_field
    
    source_field => vector_source_field(state, s_field)
    
    positions => extract_vector_field(state, "Coordinate")

    call allocate(ctfield, s_field%mesh, name="CTField")

    divergence_sparsity => get_csr_sparsity_firstorder(state, s_field%mesh, source_field%mesh)
    call allocate(CT_m, divergence_sparsity, (/1, source_field%dim/), name="DivergenceMatrix" )
    call allocate(ct_rhs, s_field%mesh, name="CTRHS")

    path = trim(complete_field_path(s_field%option_path)) // "/algorithm"
    call assemble_divergence_matrix_cg(CT_m, state, ct_rhs=ct_rhs, &
                                       test_mesh=s_field%mesh, field=source_field, &
                                       option_path=path)

    call mult(ctfield, CT_m, source_field)
    call addto(ctfield, ct_rhs, -1.0)
    call deallocate(ct_rhs)
    
    if(have_option(trim(path) // "/lump_mass")) then
      masslump => get_lumped_mass(state, s_field%mesh)
      s_field%val = ctfield%val / masslump%val
    else
      mass => get_mass_matrix(state, s_field%mesh)  
      call petsc_solve(s_field, mass, ctfield, option_path = path)
    end if

    call deallocate(CT_m)
    call deallocate(ctfield)
    
  end subroutine calculate_finite_element_divergence
  
  subroutine calculate_finite_element_divergence_transpose(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    
    character(len = OPTION_PATH_LEN) :: path
    type(cmc_matrices) :: matrices
    type(scalar_field), pointer :: source_field
    
    source_field => scalar_source_field(state, v_field)
    
    path = trim(complete_field_path(v_field%option_path)) // "/algorithm"
    call allocate(matrices, state, v_field, source_field, option_path = path, add_cmc = .false.)
    call compute_conservative(matrices, v_field, source_field)
    call deallocate(matrices)
    
  end subroutine calculate_finite_element_divergence_transpose

  subroutine calculate_curl_2d(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    
    character(len = OPTION_PATH_LEN) :: path
    integer :: i
    type(csr_matrix), pointer :: mass
    type(scalar_field), pointer :: masslump
    type(scalar_field) :: rhs
    type(vector_field), pointer :: positions, source_field
    
    ewrite(1, *) "In calculate_curl_2d"
        
    source_field => vector_source_field(state, s_field)
    ewrite(2, *) "Calculating curl of field: " // trim(source_field%name)
    ewrite(2, *) "On mesh: " // trim(source_field%mesh%name)
    ewrite(2, *) "Diagnostic field: " // trim(s_field%name)
    ewrite(2, *) "On mesh: " // trim(s_field%mesh%name)
        
    if(source_field%dim /= 2) then
      FLAbort("Can only calculate 2D curl in 2D")
    end if
    
    positions => extract_vector_field(state, "Coordinate")    
    path = trim(complete_field_path(s_field%option_path)) // "/algorithm"
    if(have_option(trim(path) // "/lump_mass")) then
      masslump => get_lumped_mass(state, s_field%mesh)
      call zero(s_field)
      do i = 1, ele_count(rhs)
        call assemble_curl_ele(i, positions, source_field, s_field)
      end do
      s_field%val = s_field%val / masslump%val
    else
      select case(continuity(s_field))
        case(0)
          if(.not. have_option(trim(path) // "/solver")) then
            FLExit("For continuous curl, must supply solver options when not lumping mass")
          end if
          mass => get_mass_matrix(state, s_field%mesh)
          call allocate(rhs, s_field%mesh, "CurlRHS")
          call zero(rhs)
          do i = 1, ele_count(rhs)
            call assemble_curl_ele(i, positions, source_field, rhs)
          end do
          call petsc_solve(s_field, mass, rhs, option_path = path)
          call deallocate(rhs)
        case(-1)
          do i = 1, ele_count(s_field)
            call solve_curl_ele(i, positions, source_field, s_field)
          end do
        case default
          ewrite(-1, *) "For mesh continuity: ", continuity(s_field)
          FLAbort("Unrecognised mesh continuity")
      end select
    end if
    
    ewrite_minmax(s_field%val)
    
    ewrite(1, *) "Exiting calculate_curl_2d"
    
  contains
  
    subroutine assemble_curl_ele(ele, positions, source, rhs)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: source
      type(scalar_field), intent(inout) :: rhs
      
      real, dimension(ele_ngi(positions, ele)) :: detwei
      real, dimension(ele_loc(source, ele), ele_ngi(source, ele), positions%dim) :: dshape
      
      call transform_to_physical(positions, ele, ele_shape(source, ele), &
        & dshape = dshape, detwei = detwei)
      
      call addto(rhs, ele_nodes(rhs, ele), &
        & shape_rhs(ele_shape(rhs, ele), &
        & ele_2d_curl_at_quad(source, ele, dshape) * detwei))
    
    end subroutine assemble_curl_ele
    
    subroutine solve_curl_ele(ele, positions, source, curl)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: source
      type(scalar_field), intent(inout) :: curl
      
      real, dimension(ele_ngi(positions, ele)) :: detwei
      real, dimension(ele_loc(curl, ele)) :: little_rhs
      real, dimension(ele_loc(curl, ele), ele_loc(curl, ele)) :: little_mass
      real, dimension(ele_loc(source, ele), ele_ngi(source, ele), positions%dim) :: dshape
      type(element_type), pointer :: shape
      
      call transform_to_physical(positions, ele, ele_shape(source, ele), &
        & dshape = dshape, detwei = detwei)
      
      shape => ele_shape(curl, ele)
      
      little_mass = shape_shape(shape, shape, detwei)
      little_rhs = shape_rhs(shape, ele_2d_curl_at_quad(source, ele, dshape) * detwei)
        
      call solve(little_mass, little_rhs)
      
      call set(curl, ele_nodes(curl, ele), little_rhs)
      
    end subroutine solve_curl_ele
    
  end subroutine calculate_curl_2d

  subroutine calculate_curl(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    
    character(len = OPTION_PATH_LEN) :: path
    integer :: i
    type(csr_matrix), pointer :: mass
    type(scalar_field), pointer :: masslump
    type(vector_field) :: rhs
    type(vector_field), pointer :: positions, source_field
    
    ewrite(1, *) "In calculate_curl"
    
    source_field => vector_source_field(state, v_field)
    ewrite(2, *) "Calculating curl of field: " // trim(source_field%name)
    ewrite(2, *) "On mesh: " // trim(source_field%mesh%name)
    ewrite(2, *) "Diagnostic field: " // trim(v_field%name)
    ewrite(2, *) "On mesh: " // trim(v_field%mesh%name)
    
    if(source_field%dim /= 3) then
      FLAbort("Can only calculate curl in 3D")
    end if
    
    positions => extract_vector_field(state, "Coordinate")    
    path = trim(complete_field_path(v_field%option_path)) // "/algorithm"
    if(have_option(trim(path) // "/lump_mass")) then
      masslump => get_lumped_mass(state, v_field%mesh)
      call zero(v_field)
      do i = 1, ele_count(rhs)
        call assemble_curl_ele(i, positions, source_field, v_field)
      end do
      do i = 1, v_field%dim
        v_field%val(i)%ptr = v_field%val(i)%ptr / masslump%val
      end do
    else
      select case(continuity(v_field))
        case(0)
          if(.not. have_option(trim(path) // "/solver")) then
            FLExit("For continuous curl, must supply solver options when not lumping mass")
          end if
          mass => get_mass_matrix(state, v_field%mesh)
          call allocate(rhs, v_field%dim, v_field%mesh, "CurlRHS")
          call zero(rhs)
          do i = 1, ele_count(rhs)
            call assemble_curl_ele(i, positions, source_field, rhs)
          end do
          call petsc_solve(v_field, mass, rhs, option_path = path)
          call deallocate(rhs)
        case(-1)
          do i = 1, ele_count(v_field)
            call solve_curl_ele(i, positions, source_field, v_field)
          end do
        case default
          ewrite(-1, *) "For mesh continuity: ", continuity(v_field)
          FLAbort("Unrecognised mesh continuity")
      end select
    end if
    
    do i = 1, v_field%dim
      ewrite_minmax(v_field%val(i)%ptr)
    end do
    
    ewrite(1, *) "Exiting calculate_curl"
    
  contains
  
    subroutine assemble_curl_ele(ele, positions, source, rhs)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: source
      type(vector_field), intent(inout) :: rhs
      
      real, dimension(ele_ngi(positions, ele)) :: detwei
      real, dimension(ele_loc(source, ele), ele_ngi(source, ele), positions%dim) :: dshape
      
      call transform_to_physical(positions, ele, ele_shape(source, ele), &
        & dshape = dshape, detwei = detwei)
      
      call addto(rhs, ele_nodes(rhs, ele), &
        & shape_vector_rhs(ele_shape(rhs, ele), &
        & ele_curl_at_quad(source, ele, dshape), detwei))
    
    end subroutine assemble_curl_ele
    
    subroutine solve_curl_ele(ele, positions, source, curl)
      integer, intent(in) :: ele
      type(vector_field), intent(in) :: positions
      type(vector_field), intent(in) :: source
      type(vector_field), intent(inout) :: curl
      
      real, dimension(ele_ngi(positions, ele)) :: detwei
      real, dimension(ele_loc(curl, ele), curl%dim) :: little_rhs
      real, dimension(ele_loc(curl, ele), ele_loc(curl, ele)) :: little_mass
      real, dimension(ele_loc(source, ele), ele_ngi(source, ele), positions%dim) :: dshape
      type(element_type), pointer :: shape
      
      call transform_to_physical(positions, ele, ele_shape(source, ele), &
        & dshape = dshape, detwei = detwei)
      
      shape => ele_shape(curl, ele)
      
      little_mass = shape_shape(shape, shape, detwei)
      little_rhs = transpose(shape_vector_rhs(shape, &
        & ele_curl_at_quad(source, ele, dshape), detwei))
        
      call solve(little_mass, little_rhs)
      
      call set(curl, ele_nodes(curl, ele), transpose(little_rhs))
      
    end subroutine solve_curl_ele
    
  end subroutine calculate_curl

  subroutine calculate_scalar_advection(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    integer :: stat
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: positions, velocity
    
    source_field => scalar_source_field(state, s_field)
    
    velocity => extract_vector_field(state, "Velocity", stat = stat)
    if(stat /= 0) then
      ewrite(0, *) "For field " // trim(s_field%name)
      ewrite(0, *) "Warning: Calculating advection with no velocity"
      call zero(s_field)
      return
    end if
    
    positions => extract_vector_field(state, "Coordinate")
    
    call u_dot_nabla(velocity, source_field, positions, s_field)
    
  end subroutine calculate_scalar_advection
  
  subroutine calculate_vector_advection(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    
    integer :: stat
    type(vector_field), pointer :: positions, source_field, velocity
    
    source_field => vector_source_field(state, v_field)
    
    velocity => extract_vector_field(state, "Velocity", stat = stat)
    if(stat /= 0) then
      ewrite(0, *) "For field " // trim(v_field%name)
      ewrite(0, *) "Warning: Calculating advection with no velocity"
      call zero(v_field)
      return
    end if
    
    positions => extract_vector_field(state, "Coordinate")
    
    call u_dot_nabla(velocity, source_field, positions, v_field)
    
  end subroutine calculate_vector_advection
  
  subroutine calculate_vector_laplacian(state, v_field)
    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    
    integer :: i
    type(scalar_field), pointer :: masslump
    type(vector_field), pointer :: positions, source_field
      
    ewrite(1, *) "In calculate_vector_laplacian"
    ewrite(2, *) "Computing laplacian for field " // trim(v_field%name)
      
    source_field => vector_source_field(state, v_field)
    assert(ele_count(source_field) == ele_count(v_field))
    assert(source_field%dim == v_field%dim)
    do i = 1, source_field%dim
      ewrite_minmax(source_field%val(i)%ptr)
    end do
    
    positions => extract_vector_field(state, "Coordinate")
    assert(ele_count(positions) == ele_count(v_field))
    
    call zero(v_field)
    do i = 1, ele_count(v_field)
      call assemble_vector_laplacian_ele(i, v_field, positions, source_field)
    end do
    do i = 1, v_field%dim
      ewrite_minmax(v_field%val(i)%ptr)
    end do
    
    masslump => get_lumped_mass(state, v_field%mesh)
    
    do i = 1, v_field%dim
      v_field%val(i)%ptr = v_field%val(i)%ptr / masslump%val
    end do
    do i = 1, v_field%dim
      ewrite_minmax(v_field%val(i)%ptr)
    end do
    
    ewrite(1, *) "Exiting calculate_vector_laplacian"
        
  end subroutine calculate_vector_laplacian
  
  subroutine assemble_vector_laplacian_ele(ele, v_field, positions, source_field)
    integer, intent(in) :: ele
    type(vector_field), intent(inout) :: v_field
    type(vector_field), intent(in) :: positions
    type(vector_field), intent(in) :: source_field
    
    integer :: i, j
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(v_field, ele)) :: detwei
    real, dimension(source_field%dim, ele_ngi(source_field, ele)) :: grad_gi
    real, dimension(ele_loc(source_field, ele), ele_ngi(source_field, ele), source_field%dim) :: dn_t
    
    assert(ele_ngi(positions, ele) == ele_ngi(v_field, ele))
    assert(ele_ngi(source_field, ele) == ele_ngi(v_field, ele))
        
    call transform_to_physical(positions, ele, ele_shape(source_field, ele), &
      & dshape = dn_t, detwei = detwei)
      
    element_nodes => ele_nodes(v_field, ele)
      
    do i = 1, source_field%dim
      do j = 1, source_field%dim
        grad_gi(j, :) = matmul(ele_val(source_field, i, ele), dn_t(:, :, j))
      end do
      
      call addto(v_field, i, element_nodes, -dshape_dot_vector_rhs(dn_t, grad_gi, detwei))
    end do
    
  end subroutine assemble_vector_laplacian_ele
  
  subroutine calculate_scalar_laplacian(state, s_field)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    
    integer :: i
    type(scalar_field), pointer :: masslump, source_field
    type(vector_field), pointer :: positions
      
    ewrite(1, *) "In calculate_scalar_laplacian"
    ewrite(2, *) "Computing laplacian for field " // trim(s_field%name)
      
    source_field => scalar_source_field(state, s_field)
    assert(ele_count(source_field) == ele_count(s_field))
    assert(mesh_dim(source_field) == mesh_dim(s_field))
    ewrite_minmax(source_field%val)
    
    positions => extract_vector_field(state, "Coordinate")
    assert(ele_count(positions) == ele_count(s_field))
    
    call zero(s_field)
    do i = 1, ele_count(s_field)
      call assemble_scalar_laplacian_ele(i, s_field, positions, source_field)
    end do
    ewrite_minmax(s_field%val)
    
    masslump => get_lumped_mass(state, s_field%mesh)
    
    s_field%val = s_field%val / masslump%val
    ewrite_minmax(s_field%val)
    
    ewrite(1, *) "Exiting calculate_scalar_laplacian"
        
  end subroutine calculate_scalar_laplacian
  
  subroutine assemble_scalar_laplacian_ele(ele, s_field, positions, source_field)
    integer, intent(in) :: ele
    type(scalar_field), intent(inout) :: s_field
    type(vector_field), intent(in) :: positions
    type(scalar_field), intent(in) :: source_field
    
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(s_field, ele)) :: detwei
    real, dimension(ele_loc(source_field, ele), ele_ngi(source_field, ele), mesh_dim(source_field)) :: dn_t
    
    assert(ele_ngi(positions, ele) == ele_ngi(s_field, ele))
    assert(ele_ngi(source_field, ele) == ele_ngi(s_field, ele))
        
    call transform_to_physical(positions, ele, ele_shape(source_field, ele), &
      & dshape = dn_t, detwei = detwei)
      
    element_nodes => ele_nodes(s_field, ele)
      
    call addto(s_field, element_nodes, -dshape_dot_vector_rhs(dn_t, transpose(ele_grad_at_quad(source_field, ele, dn_t)), detwei))
    
  end subroutine assemble_scalar_laplacian_ele

end module differential_operator_diagnostics
