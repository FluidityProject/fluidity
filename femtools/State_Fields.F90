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
module state_fields_module
  !!< Module containing general tools for discretising Finite Element problems.

  use fields
  use state_module
  use fefields
  use global_parameters, only: FIELD_NAME_LEN
  use sparsity_patterns_meshes
  use dgtools, only: get_dg_inverse_mass_matrix
  use eventcounter
  implicit none

  interface get_lumped_mass
    module procedure get_lumped_mass_single_state, get_lumped_mass_multiple_states
  end interface get_lumped_mass
  
  interface get_mass_matrix
    module procedure get_mass_matrix_single_state, get_mass_matrix_multiple_states
  end interface get_mass_matrix
  
  interface get_dg_inverse_mass
    module procedure get_dg_inverse_mass_single_state, get_dg_inverse_mass_multiple_states
  end interface get_dg_inverse_mass
  
  interface get_lumped_mass_on_submesh
    module procedure get_lumped_mass_on_submesh_single_state, get_lumped_mass_on_submesh_multiple_states
  end interface get_lumped_mass_on_submesh
  
  private
  public :: get_lumped_mass, get_lumped_mass_on_submesh, get_mass_matrix, get_dg_inverse_mass

contains

  function get_lumped_mass_single_state(state, mesh) result(lumped_mass)
    !!< extracts the lumped mass from states or creates it if it doesn't find it
    type(scalar_field), pointer :: lumped_mass
    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout) :: mesh
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    lumped_mass => get_lumped_mass(states, mesh)
    state = states(1)
  
  end function get_lumped_mass_single_state

  function get_lumped_mass_multiple_states(states, mesh) result(lumped_mass)
    !!< extracts the lumped mass from states or creates it if it doesn't find it
    type(scalar_field), pointer :: lumped_mass
    type(state_type), dimension(:), intent(inout) :: states
    type(mesh_type), intent(inout) :: mesh
    
    integer :: stat
    character(len=FIELD_NAME_LEN) :: name
    type(scalar_field) :: temp_lumped_mass
    type(vector_field), pointer :: positions
    integer, save :: last_mesh_movement = -1
    
    name = trim(mesh%name)//"LumpedMass"
    
    lumped_mass => extract_scalar_field(states, trim(name), stat)
    
    if((stat/=0).or.(eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement)) then
    
      positions => extract_vector_field(states(1), "Coordinate")
      call allocate(temp_lumped_mass, mesh, name=trim(name))
      call compute_lumped_mass(positions, temp_lumped_mass)
      call insert(states, temp_lumped_mass, trim(name))
      call deallocate(temp_lumped_mass)
      
      lumped_mass => extract_scalar_field(states, trim(name))
      last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
    end if
  
  end function get_lumped_mass_multiple_states

  function get_mass_matrix_single_state(state, mesh) result(mass)
    !!< extracts the mass from states or creates it if it doesn't find it
    type(csr_matrix), pointer :: mass
    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout) :: mesh
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    mass => get_mass_matrix(states, mesh)
    state = states(1)
  
  end function get_mass_matrix_single_state

  function get_mass_matrix_multiple_states(states, mesh) result(mass)
    !!< extracts the mass from states or creates it if it doesn't find it
    type(csr_matrix), pointer :: mass
    type(state_type), dimension(:), intent(inout) :: states
    type(mesh_type), intent(inout) :: mesh
    
    integer :: stat
    character(len=FIELD_NAME_LEN) :: name
    type(csr_matrix) :: temp_mass
    type(csr_sparsity), pointer :: temp_mass_sparsity
    type(vector_field), pointer :: positions
    
    integer, save :: last_mesh_movement = -1
    
    name = trim(mesh%name)//"MassMatrix"
    
    mass => extract_csr_matrix(states, trim(name), stat)
    
    if((stat/=0).or.(eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement)) then
      positions => extract_vector_field(states(1), "Coordinate")
      
      temp_mass_sparsity => get_csr_sparsity_firstorder(states, mesh, mesh)
      call allocate(temp_mass, temp_mass_sparsity, name=trim(name))
      call compute_mass(positions, mesh, temp_mass)
      call insert(states, temp_mass, trim(name))
      call deallocate(temp_mass)
      
      mass => extract_csr_matrix(states, trim(name))
      last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
    end if
  
  end function get_mass_matrix_multiple_states

  function get_dg_inverse_mass_single_state(state, mesh) result(inverse_mass)
    !!< extracts the dg inverse mass from state or creates it if it doesn't find it
    type(csr_matrix), pointer :: inverse_mass
    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout) :: mesh
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    inverse_mass => get_dg_inverse_mass(states, mesh)
    state = states(1)
  
  end function get_dg_inverse_mass_single_state

  function get_dg_inverse_mass_multiple_states(states, mesh) result(inverse_mass)
    !!< extracts the dg inverse mass from states or creates it if it doesn't find it
    type(csr_matrix), pointer :: inverse_mass
    type(state_type), dimension(:), intent(inout) :: states
    type(mesh_type), intent(inout) :: mesh
    
    integer :: stat
    character(len=FIELD_NAME_LEN) :: name
    type(csr_matrix) :: temp_inverse_mass
    type(vector_field), pointer :: positions
    
    integer, save :: last_mesh_movement = -1
    
    name = trim(mesh%name)//"DGInverseMassMatrix"
    
    inverse_mass => extract_csr_matrix(states, trim(name), stat)
    
    if((stat/=0).or.(eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement)) then
      positions => extract_vector_field(states(1), "Coordinate")
      
      call get_dg_inverse_mass_matrix(temp_inverse_mass, mesh, positions)
      call insert(states, temp_inverse_mass, trim(name))
      call deallocate(temp_inverse_mass)
      
      inverse_mass => extract_csr_matrix(states, trim(name))
      last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
    end if
  
  end function get_dg_inverse_mass_multiple_states

  function get_lumped_mass_on_submesh_single_state(state, mesh) result(lumped_mass)
    !!< extracts the lumped mass from states or creates it if it doesn't find it
    type(scalar_field), pointer :: lumped_mass
    type(state_type), intent(inout) :: state
    type(mesh_type), intent(inout) :: mesh
    
    type(state_type), dimension(1) :: states
    
    states = (/state/)
    lumped_mass => get_lumped_mass_on_submesh(states, mesh)
    state = states(1)
  
  end function get_lumped_mass_on_submesh_single_state

  function get_lumped_mass_on_submesh_multiple_states(states, mesh) result(lumped_mass)
    !!< extracts the lumped mass from states or creates it if it doesn't find it
    type(scalar_field), pointer :: lumped_mass
    type(state_type), dimension(:), intent(inout) :: states
    type(mesh_type), intent(inout) :: mesh
    
    integer :: stat
    character(len=FIELD_NAME_LEN) :: name
    type(scalar_field) :: temp_lumped_mass
    
    integer, save :: last_mesh_movement = -1
    
    name = trim(mesh%name)//"SubMeshLumpedMass"
    
    lumped_mass => extract_scalar_field(states, trim(name), stat)
    
    if((stat/=0).or.(eventcount(EVENT_MESH_MOVEMENT)/=last_mesh_movement)) then
      call allocate(temp_lumped_mass, mesh, name=trim(name))
      call compute_lumped_mass_on_submesh(states(1), temp_lumped_mass)
      call insert(states, temp_lumped_mass, trim(name))
      call deallocate(temp_lumped_mass)
      
      lumped_mass => extract_scalar_field(states, trim(name))
      last_mesh_movement = eventcount(EVENT_MESH_MOVEMENT)
    end if
  
  end function get_lumped_mass_on_submesh_multiple_states

end module state_fields_module
