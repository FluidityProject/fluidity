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

module field_copies_diagnostics

  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use spud
  use fldebug
  use sparse_tools
  use transform_elements
  use fields
  use state_module
  use field_options
  use diagnostic_source_fields
  use sparse_matrices_fields
  use solvers
  use smoothing_module
  use sparsity_patterns_meshes
  use state_fields_module

  implicit none
  
  private
  
  public :: calculate_scalar_copy, calculate_vector_copy, calculate_tensor_copy
  public :: calculate_extract_scalar_component
  public :: calculate_scalar_galerkin_projection, calculate_vector_galerkin_projection
  public :: calculate_helmholtz_smoothed_scalar, calculate_helmholtz_smoothed_vector, calculate_helmholtz_smoothed_tensor
  public :: calculate_lumped_mass_smoothed_scalar, calculate_lumped_mass_smoothed_vector
  public :: calculate_lumped_mass_smoothed_tensor
  public :: calculate_helmholtz_anisotropic_smoothed_scalar, calculate_helmholtz_anisotropic_smoothed_vector
  public :: calculate_helmholtz_anisotropic_smoothed_tensor

contains

  subroutine calculate_scalar_copy(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    type(scalar_field), pointer :: source_field
    integer :: stat
    
    source_field => scalar_source_field(state, s_field)
    
    call remap_field(source_field, s_field, stat)
    if(stat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
      if(.not.have_option(trim(complete_field_path(s_field%option_path))//"/algorithm/allow_discontinuous_continuous_remap")) then
        FLExit("In the scalar_copy diagnostic algorithm: remapping from a discontinuous mesh to a continuous mesh isn't allowed.")
      end if
    else if(stat==REMAP_ERR_UNPERIODIC_PERIODIC) then
      if(.not.have_option(trim(complete_field_path(s_field%option_path))//"/algorithm/allow_unperiodic_periodic_remap")) then
        FLExit("In the scalar_copy diagnostic algorithm: remapping from an unperiodic to a periodic mesh isn't allowed.")
      end if
    else if(stat==REMAP_ERR_HIGHER_LOWER_CONTINUOUS) then
      if(.not.have_option(trim(complete_field_path(s_field%option_path))//"/algorithm/allow_higher_lower_continuous_remap")) then
        FLExit("In the scalar_copy diagnostic algorithm: remapping from a higher order continuous mesh to a lower order continuous mesh isn't allowed.")
      end if
    else if(stat==REMAP_ERR_BUBBLE_LAGRANGE) then
      if(.not.have_option(trim(complete_field_path(s_field%option_path))//"/algorithm/allow_bubble_lagrange_remap")) then
        FLExit("In the scalar_copy diagnostic algorithm: remapping from a bubble mesh to a lagrange mesh isn't allowed.")
      end if
    end if

    
  end subroutine calculate_scalar_copy
  
  subroutine calculate_vector_copy(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    
    type(vector_field), pointer :: source_field
    integer :: stat
    
    source_field => vector_source_field(state, v_field)
    
    call remap_field(source_field, v_field, stat)
    if(stat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
      if(.not.have_option(trim(complete_field_path(v_field%option_path))//"/algorithm/allow_discontinuous_continuous_remap")) then
        FLExit("In the vector_copy diagnostic algorithm: remapping from a discontinuous mesh to a continuous mesh isn't allowed.")
      end if
    else if(stat==REMAP_ERR_UNPERIODIC_PERIODIC) then
      if(.not.have_option(trim(complete_field_path(v_field%option_path))//"/algorithm/allow_unperiodic_periodic_remap")) then
        FLExit("In the vector_copy diagnostic algorithm: remapping from an unperiodic to a periodic mesh isn't allowed.")
      end if
    else if(stat==REMAP_ERR_HIGHER_LOWER_CONTINUOUS) then
      if(.not.have_option(trim(complete_field_path(v_field%option_path))//"/algorithm/allow_higher_lower_continuous_remap")) then
        FLExit("In the vector_copy diagnostic algorithm: remapping from a higher order continuous mesh to a lower order continuous mesh isn't allowed.")
      end if
    else if(stat==REMAP_ERR_BUBBLE_LAGRANGE) then
      if(.not.have_option(trim(complete_field_path(v_field%option_path))//"/algorithm/allow_bubble_lagrange_remap")) then
        FLExit("In the vector_copy diagnostic algorithm: remapping from a bubble mesh to a lagrange mesh isn't allowed.")
      end if
    end if
    
  end subroutine calculate_vector_copy
  
  subroutine calculate_tensor_copy(state, t_field)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(inout) :: t_field
    
    type(tensor_field), pointer :: source_field
    integer :: stat
    
    source_field => tensor_source_field(state, t_field)
    
    call remap_field(source_field, t_field, stat)
    if(stat==REMAP_ERR_DISCONTINUOUS_CONTINUOUS) then
      if(.not.have_option(trim(complete_field_path(t_field%option_path))//"/algorithm/allow_discontinuous_continuous_remap")) then
        FLExit("In the tensor_copy diagnostic algorithm: remapping from a discontinuous mesh to a continuous mesh isn't allowed.")
      end if
    else if(stat==REMAP_ERR_UNPERIODIC_PERIODIC) then
      if(.not.have_option(trim(complete_field_path(t_field%option_path))//"/algorithm/allow_unperiodic_periodic_remap")) then
        FLExit("In the tensor_copy diagnostic algorithm: remapping from an unperiodic to a periodic mesh isn't allowed.")
      end if
    else if(stat==REMAP_ERR_HIGHER_LOWER_CONTINUOUS) then
      if(.not.have_option(trim(complete_field_path(t_field%option_path))//"/algorithm/allow_higher_lower_continuous_remap")) then
        FLExit("In the tensor_copy diagnostic algorithm: remapping from a higher order continuous mesh to a lower order continuous mesh isn't allowed.")
      end if
    else if(stat==REMAP_ERR_BUBBLE_LAGRANGE) then
      if(.not.have_option(trim(complete_field_path(t_field%option_path))//"/algorithm/allow_bubble_lagrange_remap")) then
        FLExit("In the tensor_copy diagnostic algorithm: remapping from a bubble mesh to a lagrange mesh isn't allowed.")
      end if
    end if
    
  end subroutine calculate_tensor_copy

  subroutine calculate_extract_scalar_component(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field

    logical :: allocated
    type(scalar_field), pointer :: source_field

    source_field => scalar_source_field(state, s_field, allocated = allocated)   
    call remap_field(source_field, s_field)
    if(allocated) deallocate(source_field)
    
  end subroutine calculate_extract_scalar_component
  
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
  
  subroutine calculate_helmholtz_smoothed_scalar(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    character(len = OPTION_PATH_LEN) :: path
    logical :: allocated
    real :: alpha
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: positions
    
    ewrite(1, *) "In calculate_helmholtz_smoothed_scalar"
    
    source_field => scalar_source_field(state, s_field, allocated = allocated)
    positions => extract_vector_field(state, "Coordinate")
    
    path = trim(complete_field_path(s_field%option_path)) // "/algorithm"
    call get_option(trim(path) // "/smoothing_scale_factor", alpha)
    ewrite(2, *) "alpha = ", alpha
    
    ewrite_minmax(source_field)
    call smooth_scalar(source_field, positions, s_field, alpha, path)
    ewrite_minmax(s_field)
    
    if(allocated) deallocate(source_field)
    
    ewrite(1, *) "Exiting calculate_helmholtz_smoothed_scalar"
    
  end subroutine calculate_helmholtz_smoothed_scalar

  subroutine calculate_helmholtz_smoothed_vector(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    
    character(len = OPTION_PATH_LEN) :: path
    real :: alpha
    type(vector_field), pointer :: source_field, positions
    
    ewrite(1, *) "In calculate_helmholtz_smoothed_vector"
    
    source_field => vector_source_field(state, v_field)
    positions => extract_vector_field(state, "Coordinate")
    
    path = trim(complete_field_path(v_field%option_path)) // "/algorithm"
    call get_option(trim(path) // "/smoothing_scale_factor", alpha)
    ewrite(2, *) "alpha = ", alpha
    
    call smooth_vector(source_field, positions, v_field, alpha, path)
    ewrite_minmax(v_field)
   
    ewrite(1, *) "Exiting calculate_helmholtz_smoothed_vector"
    
  end subroutine calculate_helmholtz_smoothed_vector

  subroutine calculate_helmholtz_smoothed_tensor(state, t_field)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(inout) :: t_field
    
    character(len = OPTION_PATH_LEN) :: path
    logical :: allocated
    real :: alpha
    type(vector_field), pointer :: positions
    type(tensor_field), pointer :: source_field

    ewrite(1, *) "In calculate_helmholtz_smoothed_tensor"

    positions => extract_vector_field(state, "Coordinate")
    source_field => tensor_source_field(state, t_field)

    path = trim(complete_field_path(t_field%option_path)) // "/algorithm"
    call get_option(trim(path) // "/smoothing_scale_factor", alpha)
    ewrite(2, *) "alpha = ", alpha
    call smooth_tensor(source_field, positions, t_field, alpha, path)

    ewrite_minmax(source_field)
    ewrite_minmax(t_field)
    
    ewrite(1, *) "Exiting calculate_helmholtz_smoothed_tensor"
    
  end subroutine calculate_helmholtz_smoothed_tensor

  subroutine calculate_helmholtz_anisotropic_smoothed_scalar(state, s_field)
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout) :: s_field
    
    character(len = OPTION_PATH_LEN) :: path
    logical :: allocated
    real :: alpha
    type(scalar_field), pointer :: source_field
    type(vector_field), pointer :: positions

    ewrite(1, *) "In calculate_helmholtz_anisotropic_smoothed_scalar"

    positions => extract_vector_field(state, "Coordinate")
    source_field => scalar_source_field(state, s_field, allocated = allocated)

    path = trim(complete_field_path(s_field%option_path)) // "/algorithm"
    call get_option(trim(path) // "/smoothing_scale_factor", alpha)
    ewrite(2, *) "alpha = ", alpha
    call anisotropic_smooth_scalar(source_field, positions, s_field, alpha, path)

    ewrite_minmax(source_field)
    ewrite_minmax(s_field)

    if(allocated) deallocate(source_field)
    
    ewrite(1, *) "Exiting calculate_helmholtz_anisotropic_smoothed_scalar"
    
  end subroutine calculate_helmholtz_anisotropic_smoothed_scalar

  subroutine calculate_helmholtz_anisotropic_smoothed_vector(state, v_field)
    type(state_type), intent(in) :: state
    type(vector_field), intent(inout) :: v_field
    
    character(len = OPTION_PATH_LEN) :: path
    real :: alpha
    type(vector_field), pointer :: positions, source_field

    ewrite(1, *) "In calculate_helmholtz_anisotropic_smoothed_vector"

    positions       => extract_vector_field(state, "Coordinate")
    source_field => vector_source_field(state, v_field)

    path = trim(complete_field_path(v_field%option_path)) // "/algorithm"
    call get_option(trim(path) // "/smoothing_scale_factor", alpha)
    ewrite(2, *) "alpha = ", alpha
    call anisotropic_smooth_vector(source_field, positions, v_field, alpha, path)

    ewrite_minmax(source_field)
    ewrite_minmax(v_field)
    
    ewrite(1, *) "Exiting calculate_helmholtz_anisotropic_smoothed_vector"
    
  end subroutine calculate_helmholtz_anisotropic_smoothed_vector

  subroutine calculate_helmholtz_anisotropic_smoothed_tensor(state, t_field)
    type(state_type), intent(in) :: state
    type(tensor_field), intent(inout) :: t_field
    
    character(len = OPTION_PATH_LEN) :: path
    real :: alpha
    type(vector_field), pointer :: positions
    type(tensor_field), pointer :: source_field

    ewrite(1, *) "In calculate_helmholtz_anisotropic_smoothed_tensor"

    positions       => extract_vector_field(state, "Coordinate")
    source_field => tensor_source_field(state, t_field)

    path = trim(complete_field_path(t_field%option_path)) // "/algorithm"
    call get_option(trim(path) // "/smoothing_scale_factor", alpha)
    ewrite(2, *) "alpha = ", alpha
    call anisotropic_smooth_tensor(source_field, positions, t_field, alpha, path)

    ewrite_minmax(source_field)
    ewrite_minmax(t_field)
    
    ewrite(1, *) "Exiting calculate_helmholtz_anisotropic_smoothed_tensor"
    
  end subroutine calculate_helmholtz_anisotropic_smoothed_tensor

  subroutine calculate_lumped_mass_smoothed_scalar(state, s_field)

    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: s_field
    type(scalar_field), pointer :: source_field, lumpedmass
    type(scalar_field) :: inverse_lumpedmass
    type(csr_matrix), pointer :: mass
    logical :: allocated

    ewrite(1, *) "In calculate_lumped_mass_smoothed_scalar"
    
    source_field => scalar_source_field(state, s_field, allocated = allocated)

    ! Apply smoothing filter
    call allocate(inverse_lumpedmass, source_field%mesh, "InverseLumpedMass")
    mass => get_mass_matrix(state, source_field%mesh)
    lumpedmass => get_lumped_mass(state, source_field%mesh)
    call invert(lumpedmass, inverse_lumpedmass)
    call mult( s_field, mass, source_field)
    call scale(s_field, inverse_lumpedmass) ! the averaging operator is [inv(ML)*M*]
    call deallocate(inverse_lumpedmass)
    if(allocated) deallocate(source_field)

    ewrite(1, *) "Exiting calculate_lumped_mass_smoothed_scalar"

  end subroutine calculate_lumped_mass_smoothed_scalar

  subroutine calculate_lumped_mass_smoothed_vector(state, v_field)

    type(state_type), intent(inout) :: state
    type(vector_field), intent(inout) :: v_field
    type(vector_field), pointer :: source_field
    type(scalar_field), pointer :: lumpedmass
    type(scalar_field) :: inverse_lumpedmass
    type(csr_matrix), pointer :: mass
    
    ewrite(1, *) "In calculate_lumped_mass_smoothed_vector"
    
    source_field => vector_source_field(state, v_field)

    ! Apply smoothing filter
    call allocate(inverse_lumpedmass, source_field%mesh, "InverseLumpedMass")
    mass => get_mass_matrix(state, source_field%mesh)
    lumpedmass => get_lumped_mass(state, source_field%mesh)
    call invert(lumpedmass, inverse_lumpedmass)
    call mult( v_field, mass, source_field)
    call scale(v_field, inverse_lumpedmass) ! the averaging operator is [inv(ML)*M*]
    call deallocate(inverse_lumpedmass)

    ewrite(1, *) "Exiting calculate_lumped_mass_smoothed_vector"

  end subroutine calculate_lumped_mass_smoothed_vector

  subroutine calculate_lumped_mass_smoothed_tensor(state, t_field)

    type(state_type), intent(inout) :: state
    type(tensor_field), intent(inout) :: t_field
    type(tensor_field), pointer :: source_field
    type(scalar_field), pointer :: lumpedmass
    type(scalar_field) :: inverse_lumpedmass
    type(csr_matrix), pointer :: mass
    
    ewrite(1, *) "In calculate_lumped_mass_smoothed_tensor"
    
    source_field => tensor_source_field(state, t_field)

    ! Apply smoothing filter
    call allocate(inverse_lumpedmass, source_field%mesh, "InverseLumpedMass")
    mass => get_mass_matrix(state, source_field%mesh)
    lumpedmass => get_lumped_mass(state, source_field%mesh)
    call invert(lumpedmass, inverse_lumpedmass)
    ! IS IT POSSIBLE TO MULTIPLY CSR_MATRIX BY TENSOR FIELD?
    ! SEE Sparse_Matrices_Fields/csr_mult_vector_vector
    !call mult( t_field, mass, source_field)
    call scale(t_field, inverse_lumpedmass) ! the averaging operator is [inv(ML)*M*]
    call deallocate(inverse_lumpedmass)

    ewrite(1, *) "Exiting calculate_lumped_mass_smoothed_tensor"

  end subroutine calculate_lumped_mass_smoothed_tensor

end module field_copies_diagnostics
