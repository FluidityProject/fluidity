!    Copyright (C) 2008 Imperial College London and others.
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

module lagrangian_biology
  use fldebug
  use futils
  use spud
  use state_module
  use fields
  use global_parameters, only: PYTHON_FUNC_LEN
  use parallel_tools
  use pickers_inquire
  use detector_data_types
  use detector_tools
  use detector_move_lagrangian

implicit none

  private

  public :: initialise_lagrangian_biology, lagrangian_biology_cleanup, &
            calculate_lagrangian_biology

  type(detector_linked_list), dimension(:), allocatable, target, save :: agent_arrays

contains

  subroutine initialise_lagrangian_biology(state)
    type(state_type), intent(in) :: state

    type(detector_type), pointer :: agent
    type(vector_field), pointer :: xfield
    type(element_type), pointer :: shape
    character(len=PYTHON_FUNC_LEN) :: func
    character(len = 254) :: buffer
    real, allocatable, dimension(:,:) :: coords
    real:: current_time
    integer :: i, j, dim, n_agents, n_agent_arrays

    if (.not.have_option("/ocean_biology/lagrangian_ensemble")) return

    ewrite(1,*) "In initialise_lagrangian_biology"

    n_agent_arrays = option_count("/ocean_biology/lagrangian_ensemble/agents/agent_array")
    allocate(agent_arrays(n_agent_arrays))

    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)
    xfield=>extract_vector_field(state, "Coordinate")
    shape=>ele_shape(xfield,1)

    ewrite(2,*) "Found", n_agent_arrays, "agent arrays"

    do i = 1, n_agent_arrays
       write(buffer, "(a,i0,a)") "/ocean_biology/lagrangian_ensemble/agents/agent_array[",i-1,"]"
       call get_option(trim(buffer)//"/number_of_agents", n_agents)

       call get_option(trim(buffer)//"/initial_position", func)
       allocate(coords(dim,n_agents))
       call set_detector_coords_from_python(coords, n_agents, func, current_time)

       ! Create agent and insert into list
       do j = 1, n_agents
          allocate(agent)
          call insert(agent_arrays(i), agent)

          allocate(agent%position(dim))
          agent%position=coords(:,j)

          allocate(agent%local_coords(local_coord_count(shape))) 

          agent%id_number=j
          agent%name=int2str(j)
          agent%local=.true.
       end do
       deallocate(coords)

       ! Determine elements for current agent list
       call search_for_detectors(agent_arrays(i), xfield)

       ! Delete all agents whose element we don't own
       if (isparallel()) then
          agent => agent_arrays(i)%firstnode
          do while (associated(agent))
             if (element_owner(xfield%mesh,agent%element)/=getprocno() .or. agent%element<0) then 
                call delete(agent_arrays(i),agent)
             else
                agent%local=.true.
                agent%initial_owner=getprocno()            
                agent => agent%next
             end if
          end do
       end if

       ! Get options for lagrangian detector movement
       call read_detector_move_options(agent_arrays(i),"/ocean_biology/lagrangian_ensemble/agents")
    end do 

  end subroutine initialise_lagrangian_biology

  subroutine lagrangian_biology_cleanup()
    integer :: i, n_agent_arrays

    n_agent_arrays = option_count("/ocean_biology/lagrangian_ensemble/agents/agent_array")
    do i = 1, n_agent_arrays
       call delete_all(agent_arrays(i))
    end do
    deallocate(agent_arrays)

  end subroutine lagrangian_biology_cleanup

  subroutine calculate_lagrangian_biology(state, dt, timestep)
    type(state_type), intent(in) :: state
    real, intent(in) :: dt
    integer, intent(in) :: timestep

    integer :: i, n_agent_arrays

    n_agent_arrays = option_count("/ocean_biology/lagrangian_ensemble/agents/agent_array")
    do i = 1, n_agent_arrays
       ! Move lagrangian detectors
       if (check_any_lagrangian(agent_arrays(i))) then
          call move_lagrangian_detectors(state, agent_arrays(i), dt, timestep)
       end if
    end do

  end subroutine calculate_lagrangian_biology

end module lagrangian_biology
