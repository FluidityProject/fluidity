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
module state_matrices_module
  !!< Module containing general tools for discretising Finite Element problems.

  use fields
  use state_module
  use global_parameters, only: FIELD_NAME_LEN
  use sparsity_patterns_meshes
  use divergence_matrix_cv, only: assemble_divergence_matrix_cv
  use divergence_matrix_cg, only: assemble_divergence_matrix_cg
  use gradient_matrix_cg, only: assemble_gradient_matrix_cg
  use eventcounter
  use field_options
  implicit none

  interface get_divergence_matrix_cv
    module procedure get_divergence_matrix_cv_single_state, get_divergence_matrix_cv_multiple_states
  end interface get_divergence_matrix_cv
  
  interface get_pressure_poisson_matrix
    module procedure get_pressure_poisson_matrix_single_state, get_pressure_poisson_matrix_multiple_states
  end interface get_pressure_poisson_matrix
  
  interface get_pressure_stabilisation_matrix
    module procedure get_pressure_stabilisation_matrix_single_state, &
                  get_pressure_stabilisation_matrix_multiple_states
  end interface get_pressure_stabilisation_matrix
  
  interface get_velocity_divergence_matrix
    module procedure get_velocity_divergence_matrix_single_state, get_velocity_divergence_matrix_multiple_states
  end interface get_velocity_divergence_matrix
  
  private
  public :: get_divergence_matrix_cv, get_pressure_poisson_matrix, &
            get_pressure_stabilisation_matrix, get_velocity_divergence_matrix, &
            get_hydrostatic_pressure_cg_matrix, get_vertical_balance_pressure_matrix

contains

  function get_divergence_matrix_cv_single_state(state, test_mesh, field, &
                                                 div_rhs, exclude_boundaries) result(div)
    !!< extracts the cv divergence matrix from state or creates it if it doesn't find it
    type(block_csr_matrix), pointer :: div
    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout) :: test_mesh
    type(vector_field), intent(inout) :: field
    type(scalar_field), intent(inout), optional :: div_rhs
    logical, intent(in), optional :: exclude_boundaries
    
    type(state_type), dimension(1) :: states
  
    states=(/state/)
    div=>get_divergence_matrix_cv(states, test_mesh, field, div_rhs, exclude_boundaries)
    state = states(1)
  
  end function get_divergence_matrix_cv_single_state

  function get_divergence_matrix_cv_multiple_states(states, test_mesh, field, &
                                                 div_rhs, exclude_boundaries) result(div)
    !!< extracts the cv divergence matrix from state or creates it if it doesn't find it
    type(block_csr_matrix), pointer :: div
    type(state_type), dimension(:), intent(inout) :: states
    type(mesh_type), intent(inout) :: test_mesh
    type(vector_field), intent(inout) :: field
    type(scalar_field), intent(inout), optional :: div_rhs
    logical, intent(in), optional :: exclude_boundaries
    
    integer :: stat
    character(len=FIELD_NAME_LEN) :: name
    type(block_csr_matrix) :: temp_div
    type(csr_sparsity), pointer :: temp_div_sparsity
    
    integer, save :: last_mesh_movement = -1
    
    if(present_and_true(exclude_boundaries)) then
      name = trim(test_mesh%name)//trim(field%mesh%name)//"CVDivergenceMatrixNoBoundaries"
    else
      name = trim(test_mesh%name)//trim(field%mesh%name)//"CVDivergenceMatrix"
    end if
    
    div => extract_block_csr_matrix(states, trim(name), stat)
    
    if(stat/=0) then
      ! couldn't find the matrix in state so we need to allocate and assemble it
      ! if div_rhs is present then we can assemble it here too (unless you're 
      ! excluding boundaries of course, in which case its strange you've passed in div_rhs!)
      temp_div_sparsity=>get_csr_sparsity_firstorder(states, test_mesh, field%mesh)
      call allocate(temp_div, temp_div_sparsity, (/1, field%dim/), name=trim(name))

      call assemble_divergence_matrix_cv(temp_div, states(1), ct_rhs=div_rhs, &
                                         test_mesh=test_mesh, field=field, &
                                         exclude_boundaries=exclude_boundaries)
      
      call insert(states, temp_div, trim(name))
      call deallocate(temp_div)
      
      div => extract_block_csr_matrix(states, trim(name))
      
      last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
    else if (eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement) then
      ! we found the matrix in state but the mesh has moved so we need to reassemble it
      call assemble_divergence_matrix_cv(div, states(1), ct_rhs=div_rhs, &
                                         test_mesh=test_mesh, field=field, &
                                         exclude_boundaries=exclude_boundaries)

      ! and record the mesh movement index during which we've just reassembled the matrix
      last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
    else if (present(div_rhs)) then
      ! found the div matrix but div_rhs always needs to be updated so call the assembly
      ! but tell it not to reassemble the div matrix.
      call assemble_divergence_matrix_cv(div, states(1), ct_rhs=div_rhs, &
                                         test_mesh=test_mesh, field=field, &
                                         get_ct=.false., exclude_boundaries=exclude_boundaries)
      
      ! not updating the matrix so no need to increment the mesh movement index
    end if
  
  end function get_divergence_matrix_cv_multiple_states
  
  function get_pressure_poisson_matrix_single_state(state, get_cmc) result(cmc_m)
    !!< extracts the cmc matrix from state, 
    !!< if it fails to find it it returns get_cmc=.true. to indicate that it needs assembling
    type(csr_matrix), pointer :: cmc_m
    type(state_type), intent(inout) :: state
    logical, intent(inout), optional :: get_cmc
  
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    cmc_m => get_pressure_poisson_matrix(states, get_cmc=get_cmc)
    state = states(1)
  
  end function get_pressure_poisson_matrix_single_state
  
  function get_pressure_poisson_matrix_multiple_states(states, get_cmc) result(cmc_m)
    !!< extracts the cmc matrix from states, 
    !!< if it fails to find it it returns get_cmc=.true. to indicate that it needs assembling
    type(csr_matrix), pointer :: cmc_m
    type(state_type), dimension(:), intent(inout) :: states
    logical, intent(inout), optional :: get_cmc
    
    integer :: stat
    type(mesh_type), pointer :: p_mesh, u_mesh
    type(csr_sparsity), pointer :: cmc_sparsity
    type(csr_matrix) :: temp_cmc_m
    
    integer, save :: last_mesh_movement = -1
    
    
    if(present(get_cmc)) get_cmc = .false.
    cmc_m => extract_csr_matrix(states, "PressurePoissonMatrix", stat)
    
    if(stat/=0) then
      if(present(get_cmc)) get_cmc = .true.
    
      p_mesh => extract_pressure_mesh(states)
      u_mesh => extract_velocity_mesh(states)
      
      cmc_sparsity => get_csr_sparsity_secondorder(states, p_mesh, u_mesh)
      
      call allocate(temp_cmc_m, cmc_sparsity, name="PressurePoissonMatrix")
      call insert(states, temp_cmc_m, name="PressurePoissonMatrix")
      call deallocate(temp_cmc_m)
      
      cmc_m => extract_csr_matrix(states, "PressurePoissonMatrix")
    else
      ! We found it in state so the only thing that will affect if we need to assemble it
      ! is whether the mesh has moved since we last called this subroutine.
      ! The actual assembly doesn't take place here though so we don't need to do anything
      ! else yet.
      if(present(get_cmc)) get_cmc = (eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement)
    end if
    
    last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
  
  end function get_pressure_poisson_matrix_multiple_states
  
  function get_pressure_stabilisation_matrix_single_state(state) result(kmk_m)
    !!< extracts the kmk matrix from state, 
    type(csr_matrix), pointer :: kmk_m
    type(state_type), intent(inout) :: state
  
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    kmk_m => get_pressure_stabilisation_matrix(states)
    state = states(1)
  
  end function get_pressure_stabilisation_matrix_single_state
  
  function get_pressure_stabilisation_matrix_multiple_states(states) result(kmk_m)
    !!< extracts the kmk matrix from states, 
    type(csr_matrix), pointer :: kmk_m
    type(state_type), dimension(:), intent(inout) :: states
    
    integer :: stat
    type(mesh_type), pointer :: p_mesh, u_mesh
    type(csr_sparsity), pointer :: cmc_sparsity
    type(csr_matrix) :: temp_cmc_m
    
    kmk_m => extract_csr_matrix(states, "PressureStabilisationMatrix", stat)
    
    if(stat/=0) then
      p_mesh => extract_pressure_mesh(states)
      u_mesh => extract_velocity_mesh(states)
      
      cmc_sparsity => get_csr_sparsity_secondorder(states, p_mesh, u_mesh)
      
      call allocate(temp_cmc_m, cmc_sparsity, name="PressureStabilisationMatrix")
      call insert(states, temp_cmc_m, name="PressureStabilisationMatrix")
      call deallocate(temp_cmc_m)
      
      kmk_m => extract_csr_matrix(states, "PressureStabilisationMatrix")
    end if
  
  end function get_pressure_stabilisation_matrix_multiple_states
  
  function get_velocity_divergence_matrix_single_state(state, get_ct) result(ct_m)
    !!< extracts the ct matrix from state, 
    !!< if it fails to find it it returns get_ct=.true. to indicate that it needs assembling
    type(block_csr_matrix), pointer :: ct_m
    type(state_type), intent(inout) :: state
    logical, intent(inout), optional :: get_ct
  
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    ct_m => get_velocity_divergence_matrix(states, get_ct=get_ct)
    state = states(1)
  
  end function get_velocity_divergence_matrix_single_state

  function get_velocity_divergence_matrix_multiple_states(states, get_ct) result(ct_m)
    !!< extracts the ct matrix from states, 
    !!< if it fails to find it it returns get_ct=.true. to indicate that it needs assembling
    type(block_csr_matrix), pointer :: ct_m
    type(state_type), dimension(:), intent(inout) :: states
    logical, intent(inout), optional :: get_ct
    
    integer :: stat, i
    type(mesh_type), pointer :: p_mesh, u_mesh
    type(vector_field), pointer :: velocity
    type(csr_sparsity), pointer :: ct_sparsity
    type(block_csr_matrix) :: temp_ct_m
    
    integer, save :: last_mesh_movement = -1
    
    
    if(present(get_ct)) get_ct = .false.
    ct_m => extract_block_csr_matrix(states, "VelocityDivergenceMatrix", stat)
    
    if(stat/=0) then
      if(present(get_ct)) get_ct = .true.
    
      p_mesh => extract_pressure_mesh(states)
      u_mesh => extract_velocity_mesh(states)
      do i = 1, size(states)
        velocity => extract_vector_field(states(i), "Velocity", stat)
        if(stat==0) exit
      end do
      
      ct_sparsity => get_csr_sparsity_firstorder(states, p_mesh, u_mesh)
      
      call allocate(temp_ct_m, ct_sparsity, blocks=(/1,velocity%dim/), name="VelocityDivergenceMatrix")
      call insert(states, temp_ct_m, name="VelocityDivergenceMatrix")
      call deallocate(temp_ct_m)
      
      ct_m => extract_block_csr_matrix(states, "VelocityDivergenceMatrix")
    else
      ! just check if we need to reassemble the matrix anyway
      if(present(get_ct)) get_ct = (eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement)
    end if
    
    ! record the last time this subroutine was called relative to movements of the mesh
    last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
  
  end function get_velocity_divergence_matrix_multiple_states
  
  function get_hydrostatic_pressure_cg_matrix(state, assemble_matrix) result(matrix)
    !!< extracts the continuous hydrostatic pressure matrix from state, 
    !!< if it fails to find it it returns assemble_matrix=.true. to indicate that it needs assembling
    type(csr_matrix), pointer :: matrix
    type(state_type), intent(inout) :: state
    logical, intent(inout), optional :: assemble_matrix
    
    integer :: stat
    type(scalar_field), pointer :: hp
    type(csr_sparsity), pointer :: matrix_sparsity
    type(csr_matrix) :: temp_matrix
    
    integer, save :: last_mesh_movement = -1
    
    if(present(assemble_matrix)) assemble_matrix = .false.
    matrix => extract_csr_matrix(state, "HydrostaticPressureCGMatrix", stat)
    
    if(stat/=0) then
      if(present(assemble_matrix)) assemble_matrix = .true.
    
      hp => extract_scalar_field(state, "HydrostaticPressure")
      
      matrix_sparsity => get_csr_sparsity_firstorder(state, hp%mesh, hp%mesh)
      
      call allocate(temp_matrix, matrix_sparsity, name="HydrostaticPressureCGMatrix")
      call insert(state, temp_matrix, name="HydrostaticPressureCGMatrix")
      call deallocate(temp_matrix)
      
      matrix => extract_csr_matrix(state, "HydrostaticPressureCGMatrix")
    else
      ! We found it in state so the only thing that will affect if we need to assemble it
      ! is whether the mesh has moved since we last called this subroutine.
      ! The actual assembly doesn't take place here though so we don't need to do anything
      ! else yet.
      if(present(assemble_matrix)) assemble_matrix = (eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement)
    end if
    
    last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
  
  end function get_hydrostatic_pressure_cg_matrix

  function get_vertical_balance_pressure_matrix(state, assemble_matrix) result(matrix)
    !!< extracts the vertical balance pressure matrix from state, 
    !!< if it fails to find it it returns assemble_matrix=.true. to indicate that it needs assembling
    type(csr_matrix), pointer :: matrix
    type(state_type), intent(inout) :: state
    logical, intent(inout), optional :: assemble_matrix
    
    integer :: stat
    type(scalar_field), pointer :: vbp
    type(csr_sparsity), pointer :: matrix_sparsity
    type(csr_matrix) :: temp_matrix
    
    integer, save :: last_mesh_movement = -1
    
    if(present(assemble_matrix)) assemble_matrix = .false.
    matrix => extract_csr_matrix(state, "VerticalBalancePressureMatrix", stat)
    
    if(stat/=0) then
      if(present(assemble_matrix)) assemble_matrix = .true.
    
      vbp => extract_scalar_field(state, "VerticalBalancePressure")
      
      matrix_sparsity => get_csr_sparsity_firstorder(state, vbp%mesh, vbp%mesh)
      
      call allocate(temp_matrix, matrix_sparsity, name="VerticalBalancePressureMatrix")
      call insert(state, temp_matrix, name="VerticalBalancePressureMatrix")
      call deallocate(temp_matrix)
      
      matrix => extract_csr_matrix(state, "VerticalBalancePressureMatrix")
    else
      ! We found it in state so the only thing that will affect if we need to assemble it
      ! is whether the mesh has moved since we last called this subroutine.
      ! The actual assembly doesn't take place here though so we don't need to do anything
      ! else yet.
      if(present(assemble_matrix)) assemble_matrix = (eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement)
    end if
    
    last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
  
  end function get_vertical_balance_pressure_matrix

end module state_matrices_module
