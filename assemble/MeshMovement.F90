#include "fdebug.h"

module meshmovement

  use global_parameters
  use fldebug
  use vector_tools, only: solve, mat_diag_mat, eigendecomposition_symmetric
  use element_numbering
  use elements
  use shape_functions
  use spud
  use sparse_tools
  use fields_base
  use global_numbering
  use eventcounter
  use fetools
  use unittest_tools
  use fields
  use state_module
  use vtk_interfaces
  use sparse_matrices_fields
  use solvers
  use fefields
  use field_derivatives
  use sparsity_patterns
  use sparsity_patterns_meshes

  implicit none
  integer,save :: MeshCount=0

  interface

     subroutine set_debug_level(level)
       implicit none
       integer, intent(in) :: level
     end subroutine set_debug_level

     subroutine reset_debug_level
     end subroutine reset_debug_level
  end interface
  
  private
  
  public :: move_mesh_imposed_velocity, move_mesh_pseudo_lagrangian

contains

  subroutine move_mesh_imposed_velocity(states)
  type(state_type), dimension(:), intent(inout) :: states
  
    type(vector_field), pointer :: coordinate, old_coordinate, new_coordinate
    type(vector_field), pointer :: velocity
    type(vector_field), pointer :: grid_velocity
    
    integer :: i, stat
    real :: itheta, dt
    logical :: found_velocity
    
    if(.not.have_option("/mesh_adaptivity/mesh_movement/imposed_grid_velocity")) return
    call IncrementEventCounter(EVENT_MESH_MOVEMENT)
    
    ewrite(1,*) 'Entering move_mesh_imposed_velocity'
    
    grid_velocity => extract_vector_field(states(1), "GridVelocity")
    
    coordinate => extract_vector_field(states(1), "Coordinate")
    old_coordinate => extract_vector_field(states(1), "OldCoordinate")
    new_coordinate => extract_vector_field(states(1), "IteratedCoordinate")
    
    call get_option("/timestepping/timestep", dt)
    
    found_velocity = .false.
    do i = 1, size(states)
      velocity => extract_vector_field(states(i), "Velocity", stat)
      if(stat==0 .and. .not. velocity%aliased) then
        call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation", itheta, stat)
        if(found_velocity.and.(stat==0)) then
          FLExit("Only one prognostic velocity allowed with imposed mesh movement.")
        else
          found_velocity = (stat==0)
        end if
      end if
    end do
    if(.not.found_velocity) then
      itheta = 0.5
    end if
    
    call set(new_coordinate, old_coordinate)
    call addto(new_coordinate, grid_velocity, scale=dt)
    
    call set(coordinate, new_coordinate, old_coordinate, itheta)
  
  end subroutine move_mesh_imposed_velocity

  subroutine move_mesh_pseudo_lagrangian(states)
  type(state_type), dimension(:), intent(inout) :: states
  
    type(vector_field), pointer :: coordinate, old_coordinate, new_coordinate
    type(vector_field), pointer :: velocity
    type(vector_field), pointer :: grid_velocity
    
    integer :: i, stat
    real :: itheta, dt
    logical :: found_velocity
    
    character(len=FIELD_NAME_LEN) :: state_name
    
    if(.not.have_option("/mesh_adaptivity/mesh_movement/pseudo_lagrangian")) return
    call IncrementEventCounter(EVENT_MESH_MOVEMENT)
    
    ewrite(1,*) 'Entering move_mesh_pseudo_lagrangian'
    
    grid_velocity => extract_vector_field(states(1), "GridVelocity")
    
    call get_option("/mesh_adaptivity/mesh_movement/pseudo_lagrangian/velocity_material_phase/material_phase_name", &
                    state_name, stat=stat)
    if(stat==0) then
      i = get_state_index(states, trim(state_name))
      velocity => extract_vector_field(states(i), "Velocity")
    else
      velocity => extract_vector_field(states(1), "Velocity")
    end if
    
    call set(grid_velocity, velocity)
    
    coordinate => extract_vector_field(states(1), "Coordinate")
    old_coordinate => extract_vector_field(states(1), "OldCoordinate")
    new_coordinate => extract_vector_field(states(1), "IteratedCoordinate")
    
    call get_option("/timestepping/timestep", dt)
    
    found_velocity = .false.
    do i = 1, size(states)
      velocity => extract_vector_field(states(i), "Velocity", stat)
      if(stat==0) then
        call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/relaxation", itheta, stat)
        if(found_velocity.and.(stat==0)) then
          FLExit("Only one prognostic velocity allowed with pseudo lagrangian mesh movement.")
        else
          found_velocity = (stat==0)
        end if
      end if
    end do
    if(.not.found_velocity) then
      itheta = 0.5
    end if
    
    call set(new_coordinate, old_coordinate)
    call addto(new_coordinate, grid_velocity, scale=dt)
    
    call set(coordinate, new_coordinate, old_coordinate, itheta)
  
  end subroutine move_mesh_pseudo_lagrangian

end module meshmovement


