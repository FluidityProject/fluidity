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

module advection_diffusion_fv
  !!< This module contains the finite volume form of the advection
  !!< -diffusion equation for scalars.
  use quadrature
  use elements
  use sparse_tools
  use fields
  
  use fetools
  use state_module
  use shape_functions
  use transform_elements
  use fldebug
  use petsc_solve_state_module
  use boundary_conditions
  use boundary_conditions_from_options
  use spud
  use field_options
  use sparsity_patterns_meshes
  use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN
  use profiler
  
  implicit none

  private
  public solve_advection_diffusion_fv, advection_diffusion_fv_check_options

  ! Mass term?
  logical :: have_mass
  ! Advection?
  logical :: have_advection
  ! Source?
  logical :: have_source
  ! Absorption?
  logical :: have_absorption
  ! Diffusivity?
  logical :: have_diffusivity
  ! Isotropic diffusivity?
  logical :: isotropic_diffusivity
  ! Is the mesh moving?
  logical :: move_mesh
  
  ! timestepping parameters
  real :: theta, dt, dt_theta

contains

  subroutine solve_advection_diffusion_fv(field_name, state)
    !!< Construct and solve the advection-diffusion equation for the given
    !!< field unsing element centred finite volumes.
    
    !! Name of the field to be solved for.
    character(len=*), intent(in) :: field_name
    !! Collection of fields defining system state.
    type(state_type), intent(inout) :: state

    !! Tracer to be solved for.
    type(scalar_field), pointer :: T, T_old
    !! Change in T over one timestep.
    type(scalar_field) :: delta_T
    !! System matrix.
    type(csr_matrix) :: matrix
    !! Right hand side vector.
    type(scalar_field) :: rhs
    !! Sparsity of advection_diffusion matrix
    type(csr_sparsity), pointer :: sparsity
    
    ewrite(1,*) 'In solve_advection_diffusion_fv'

    t=>extract_scalar_field(state, field_name)
    if((continuity(T)>=0).or.(element_degree(T,1)/=0)) then
      FLExit("FV advection-diffusion requires a discontinuous mesh.")
    end if
    
    t_old=>extract_scalar_field(state, "Old"//field_name)

    ! Reset T to value at the beginning of the timestep.
    call set(t, t_old)

    sparsity => get_csr_sparsity_firstorder(state, t%mesh, t%mesh)
    
    call allocate(matrix, sparsity, name = trim(field_name)//"Matrix")
    call allocate(rhs, t%mesh, name = trim(field_name)//"RHS")
    
    call allocate(delta_t, t%mesh, "Delta"//trim(field_name))
    call zero(delta_t)
    
    call get_option("/timestepping/timestep", dt)
    
    call assemble_advection_diffusion_fv(t, matrix, rhs, state)
    
    call petsc_solve(delta_t, matrix, rhs, state, option_path = trim(t%option_path))
    
    ewrite_minmax(delta_t%val)
    
    call addto(t, delta_t, dt)
    
    ewrite_minmax(t%val)
    
    call deallocate(matrix)
    call deallocate(rhs)
    call deallocate(delta_t)

    ewrite(1,*) 'Exiting solve_advection_diffusion_fv'

  end subroutine solve_advection_diffusion_fv
  
  subroutine assemble_advection_diffusion_fv(t, matrix, rhs, state)
  
    type(scalar_field), intent(inout) :: t
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    type(state_type), intent(in) :: state
    
    type(vector_field), pointer :: coordinate, &
                                   old_coordinate, new_coordinate, &
                                   grid_velocity
    type(scalar_field), pointer :: source, absorption
    type(tensor_field), pointer :: diffusivity
    
    integer :: i, j, ele, stat
    
    ewrite(1,*) "In assemble_advection_diffusion_fv"
    
    coordinate => extract_vector_field(state, "Coordinate")
    assert(coordinate%dim == mesh_dim(t))
    assert(ele_count(coordinate) == ele_count(t))
    
    ! Source
    source => extract_scalar_field(state, trim(t%name)//"Source", stat = stat)
    have_source = stat == 0
    if(have_source) then
      assert(mesh_dim(source) == mesh_dim(t))
      assert(ele_count(source) == ele_count(t))
    
      ewrite_minmax(source%val)
    else
      ewrite(2,*) 'No source'
    end if

    ! Absorption
    absorption => extract_scalar_field(state, trim(t%name) // "Absorption", stat = stat)
    have_absorption = stat == 0
    if(have_absorption) then
      assert(mesh_dim(absorption) == mesh_dim(t))
      assert(ele_count(absorption) == ele_count(t))
    
      ewrite_minmax(absorption%val)
    else
      ewrite(2, *) "No absorption"
    end if

    ! Diffusivity
    diffusivity => extract_tensor_field(state, trim(t%name) // "Diffusivity", stat = stat)
    have_diffusivity = stat == 0
    if(have_diffusivity) then
      assert(diffusivity%dim == mesh_dim(t))
      assert(ele_count(diffusivity) == ele_count(t))
      
      isotropic_diffusivity = option_count(complete_field_path(diffusivity%option_path)) &
        & == option_count(trim(complete_field_path(diffusivity%option_path)) // "/value/isotropic")
        
      if(isotropic_diffusivity) then
        ewrite(2, *) "Isotropic diffusivity"
        assert(diffusivity%dim > 0)
        ewrite_minmax(diffusivity%val(1, 1, :))
      else
        do i = 1, diffusivity%dim
          do j = 1, diffusivity%dim
            ewrite_minmax(diffusivity%val(i, j, :))
          end do
        end do
      end if
    else
      isotropic_diffusivity = .false.
      ewrite(2, *) "No diffusivity"
    end if
    
    call get_option(trim(t%option_path) // "/prognostic/temporal_discretisation/theta", theta)
    assert(theta >= 0.0 .and. theta <= 1.0)
    ewrite(2, *) "Theta = ", theta
    dt_theta = dt*theta

    have_advection = .not. have_option(trim(t%option_path) // "/prognostic/spatial_discretisation/finite_volume/advection_terms/exclude_advection_terms")
    if(have_advection) then
      FLExit("Including advection not currently supported with FV")
      ewrite(2, *) "Including advection"
    else
      ewrite(2, *) "Excluding advection"
    end if
    
    have_mass = .not. have_option(trim(t%option_path) // "/prognostic/spatial_discretisation/finite_volume/mass_terms/exclude_mass_terms")
    if(have_mass) then
      ewrite(2, *) "Including mass"
    else
      ewrite(2, *) "Excluding mass"
    end if
    
    ! are we moving the mesh?
    move_mesh = (have_option("/mesh_adaptivity/mesh_movement") .and. have_mass)
    if(move_mesh) then
      FLExit("Moving the mesh not currently supported with FV")
      ewrite(2,*) "Moving the mesh"
      old_coordinate => extract_vector_field(state, "OldCoordinate")
      new_coordinate => extract_vector_field(state, "IteratedCoordinate")
      
      ! Grid velocity
      grid_velocity => extract_vector_field(state, "GridVelocity")
      assert(grid_velocity%dim == mesh_dim(t))
      assert(ele_count(grid_velocity) == ele_count(t))
      
      ewrite(2, *) "Grid velocity:"    
      do i = 1, grid_velocity%dim
        ewrite_minmax(grid_velocity%val(i)%ptr)
      end do
    else
      ewrite(2,*) "Not moving the mesh"
    end if
    
    call zero(matrix)
    call zero(rhs)
    
    do ele = 1, ele_count(t)
      call assemble_advection_diffusion_element_fv(ele, t, matrix, rhs, &
                                                   coordinate, &
                                                   source, absorption, diffusivity)
    end do

    ewrite(2, *) "Applying strong Dirichlet boundary conditions"
    call apply_dirichlet_conditions(matrix, rhs, t, dt)
    
    ewrite_minmax(rhs%val)

    ewrite(1,*) "Exiting assemble_advection_diffusion_fv"

  end subroutine assemble_advection_diffusion_fv
  
  subroutine assemble_advection_diffusion_element_fv(ele, t, matrix, rhs, &
                                                   coordinate, &
                                                   source, absorption, diffusivity)
                                                   
    integer, intent(in) :: ele
    type(scalar_field), intent(in) :: t
    type(csr_matrix), intent(inout) :: matrix
    type(scalar_field), intent(inout) :: rhs
    type(vector_field), intent(in) :: coordinate
    type(scalar_field), intent(in) :: source
    type(scalar_field), intent(in) :: absorption
    type(tensor_field), intent(in) :: diffusivity
    
    integer, dimension(:), pointer :: element_nodes
    real, dimension(ele_ngi(t, ele)) :: detwei
    type(element_type), pointer :: t_shape

    ! What we will be adding to the matrix and RHS - assemble these as we
    ! go, so that we only do the calculations we really need
    real, dimension(ele_loc(t, ele)) :: rhs_addto
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)) :: matrix_addto

    assert(element_degree(t,ele)==0)
    t_shape => ele_shape(t, ele)

    if(have_mass.or.have_source.or.have_absorption) then
      call transform_to_physical(coordinate, ele, detwei=detwei)
    end if
    
    ! Mass
    if(have_mass) call add_mass_element_fv(ele, t_shape, t, detwei, matrix_addto)
    
    ! Absorption
    if(have_absorption) call add_absorption_element_fv(ele, t_shape, t, absorption, detwei, matrix_addto, rhs_addto)
    
    ! Source
    if(have_source) call add_source_element_fv(ele, t_shape, t, source, detwei, rhs_addto)
     
    ! Diffusivity
!     if(have_diffusivity) call add_diffusivity_element_cg(ele, t, diffusivity, dt_t, detwei, matrix_addto, rhs_addto)
    
    element_nodes => ele_nodes(t, ele)
    call addto(matrix, element_nodes, element_nodes, matrix_addto)
    call addto(rhs, element_nodes, rhs_addto)

  end subroutine assemble_advection_diffusion_element_fv

  subroutine add_mass_element_fv(ele, t_shape, t, detwei, matrix_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: t_shape
    type(scalar_field), intent(in) :: t
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)), intent(inout) :: matrix_addto
    
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)) :: mass_matrix
    
    assert(have_mass)
    
    mass_matrix = shape_shape(t_shape, t_shape, detwei)
    
    matrix_addto = matrix_addto + mass_matrix
  
  end subroutine add_mass_element_fv

  subroutine add_absorption_element_fv(ele, t_shape, t, absorption, detwei, matrix_addto, rhs_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: t_shape
    type(scalar_field), intent(in) :: t
    type(scalar_field), intent(in) :: absorption
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)), intent(inout) :: matrix_addto
    real, dimension(ele_loc(t, ele)), intent(inout) :: rhs_addto
    
    real, dimension(ele_loc(t, ele), ele_loc(t, ele)) ::  absorption_mat
    
    assert(have_absorption)
    
    absorption_mat = shape_shape(t_shape, t_shape, detwei * ele_val_at_quad(absorption, ele))
    
    if(abs(dt_theta) > epsilon(0.0)) matrix_addto = matrix_addto + dt_theta * absorption_mat
    
    rhs_addto = rhs_addto - matmul(absorption_mat, ele_val(t, ele))
    
  end subroutine add_absorption_element_fv

  subroutine add_source_element_fv(ele, t_shape, t, source, detwei, rhs_addto)
    integer, intent(in) :: ele
    type(element_type), intent(in) :: t_shape
    type(scalar_field), intent(in) :: t
    type(scalar_field), intent(in) :: source
    real, dimension(ele_ngi(t, ele)), intent(in) :: detwei
    real, dimension(ele_loc(t, ele)), intent(inout) :: rhs_addto
   
    assert(have_source)
   
    rhs_addto = rhs_addto + shape_rhs(t_shape, detwei * ele_val_at_quad(source, ele))
    
  end subroutine add_source_element_fv
  
  subroutine advection_diffusion_fv_check_options
  
    character(len = FIELD_NAME_LEN) :: field_name, state_name
    character(len = OPTION_PATH_LEN) :: path
    integer :: i, j, stat
    real :: beta, l_theta
    
    if(option_count("/material_phase/scalar_field/prognostic/spatial_discretisation/finite_volume") == 0) then
      ! Nothing to check
      return
    end if
        
    ewrite(2, *) "Checking FV advection-diffusion options"
    
    do i = 0, option_count("/material_phase") - 1
      path = "/material_phase[" // int2str(i) // "]"
      call get_option(trim(path) // "/name", state_name)
      
      do j = 0, option_count(trim(path) // "/scalar_field") - 1
        path = "/material_phase[" // int2str(i) // "]/scalar_field[" // int2str(j) // "]"
        call get_option(trim(path) // "/name", field_name)
        
        if(field_name /= "Pressure") then
        
          path = trim(path) // "/prognostic"
          
          if(have_option(trim(path) // "/spatial_discretisation/finite_volume").and.&
             have_option(trim(path) // "/equation[0]")) then       
             
            call field_error(state_name, field_name, &
                             "Finite volume spatial_discretisation under development")
            
            call get_option(trim(path) // "/spatial_discretisation/conservative_advection", beta, stat)
            if(stat == SPUD_NO_ERROR) then
              if(beta < 0.0 .or. beta > 1.0) then
              
                call field_error(state_name, field_name, &
                  & "Conservative advection factor (beta) must be >= 0.0 and <= 1.0")
              end if
            else
              call field_error(state_name, field_name, &
                & "Conservative advection factor (beta) required")
            end if
            
            call get_option(trim(path) // "/temporal_discretisation/theta", l_theta, stat)
            if(stat == SPUD_NO_ERROR) then
              if(l_theta < 0. .or. l_theta > 1.0) then
                call field_error(state_name, field_name, &
                  &"Implicitness factor (theta) must be >= 0.0 and <= 1.0")
              end if
            else
              call field_error(state_name, field_name, &
                & "Implicitness factor (theta) required")
            end if
            if(have_option(trim(path) // "/spatial_discretisation/finite_volume/mass_terms/exclude_mass_terms") .and. &
              & abs(l_theta - 1.0) > epsilon(0.0)) then
              call field_warning(state_name, field_name, &
                & "Implicitness factor (theta) should = 1.0 when excluding mass")
            end if
  
            if(have_option(trim(path) // "/spatial_discretisation/finite_volume/advection_terms/exclude_advection_terms")) then
              if(have_option(trim(path) // "/scalar_field::SinkingVelocity")) then
                call field_warning(state_name, field_name, &
                  & "SinkingVelocity set, but advection terms have been excluded - SinkingVelocity will have no effect")
              end if
            end if
  
            if(option_count(trim(path) // "/boundary_conditions/type::neumann") > 0 &
              & .and. .not. (have_option(trim(path) // "/tensor_field::Diffusivity") &
              & .or. have_option(trim(path) // "/subgridscale_parameterisation::k-epsilon") &
              & .or. have_option(trim(path) // "/subgridscale_parameterisation::GLS"))) then
                call field_warning(state_name, field_name, &
                & "Neumann boundary condition set, but have no diffusivity - boundary condition will not be applied")
            end if
          end if
        end if
      end do
    end do
    
    ewrite(2, *) "Finished checking CG advection-diffusion options"

    contains
  
    subroutine field_warning(state_name, field_name, msg)
      character(len = *), intent(in) :: state_name
      character(len = *), intent(in) :: field_name
      character(len = *), intent(in) :: msg
      
      ewrite(0, *) "Warning: For field " // trim(field_name) // " in state " // trim(state_name)
      ewrite(0, *) trim(msg)
    
    end subroutine field_warning
  
    subroutine field_error(state_name, field_name, msg)
      character(len = *), intent(in) :: state_name
      character(len = *), intent(in) :: field_name
      character(len = *), intent(in) :: msg
      
      ewrite(-1, *) "For field " // trim(field_name) // " in state " // trim(state_name)
      FLExit(trim(msg))
    
    end subroutine field_error
  
  end subroutine advection_diffusion_fv_check_options


end module advection_diffusion_fv
