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

module solid_update
  use fldebug
  use state_module
  use fields
  use spud
  use eventcounter
  implicit none
  
  private
  public :: displacement_update, solid_coordinate_update

contains

  subroutine displacement_update(state)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer :: displacement
    type(vector_field), pointer :: velocity
    real :: dt

    displacement => extract_vector_field(state, "Displacement")
    velocity => extract_vector_field(state, "Velocity")
    call get_option("/timestepping/timestep", dt)

    call addto(displacement, velocity, dt)

  end subroutine displacement_update

  subroutine solid_coordinate_update(state)

    type(state_type), intent(in) :: state
    type(vector_field), pointer :: coordinate
    type(vector_field), pointer :: gridvelocity
    real :: dt

    call IncrementEventCounter(EVENT_MESH_MOVEMENT)

    coordinate => extract_vector_field(state, "Coordinate")
    gridvelocity => extract_vector_field(state, "GridVelocity")
    call get_option("/timestepping/timestep", dt)

    call addto(coordinate, gridvelocity, dt)

  end subroutine solid_coordinate_update
end module solid_update


