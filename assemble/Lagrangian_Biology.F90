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
  use mpi_interfaces
  use pickers_inquire
  use detector_data_types
  use detector_tools
  use detector_move_lagrangian
  use detector_parallel
  use diagnostic_variables, only: initialise_constant_diagnostics, field_tag, write_detectors

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
    character(len = 254) :: buffer, schema_buffer, agent_array_name
    real, allocatable, dimension(:,:) :: coords
    real:: current_time
    integer :: i, j, dim, n_agents, n_agent_arrays, column, ierror, det_type

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
       write(schema_buffer, "(a,i0,a)") "/ocean_biology/lagrangian_ensemble/agents/agent_array[",i-1,"]"
       call get_option(trim(schema_buffer)//"/number_of_agents", n_agents)
       call get_option(trim(schema_buffer)//"/name", agent_array_name)

       call get_option(trim(schema_buffer)//"/initial_position", func)
       allocate(coords(dim,n_agents))
       call set_detector_coords_from_python(coords, n_agents, func, current_time)

       ! Register the agent array, so Zoltan/Adaptivity will not forget about it
       call register_detector_list(agent_arrays(i))

       ! Get options for lagrangian detector movement
       call read_detector_move_options(agent_arrays(i),"/ocean_biology/lagrangian_ensemble/agents")

       ! Collect other meta-information
       if (have_option(trim(schema_buffer)//"/binary_output")) then
          agent_arrays(i)%binary_output=.true.
       else
          agent_arrays(i)%binary_output=.false.
       end if
       if (have_option(trim(schema_buffer)//"/lagrangian")) then
          det_type=LAGRANGIAN_DETECTOR
       else
          det_type=STATIC_DETECTOR
       end if
       if (have_option(trim(schema_buffer)//"/exclude_from_advection")) then
          agent_arrays(i)%move_parameters%do_velocity_advect=.false.
       else
          agent_arrays(i)%move_parameters%do_velocity_advect=.true.
       end if
       if (have_option(trim(schema_buffer)//"/random_walk")) then
          agent_arrays(i)%move_parameters%do_random_walk=.true.
          call get_option(trim(schema_buffer)//"/random_walk", agent_arrays(i)%move_parameters%rw_pycode)
       else
          agent_arrays(i)%move_parameters%do_random_walk=.false.
       end if
       agent_arrays(i)%total_num_det=n_agents


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
          agent%type=det_type
       end do
       deallocate(coords)

       ! Determine elements for current agent list
       call search_for_detectors(agent_arrays(i), xfield)

       ! Create simple position-only agent I/O header 
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if (getprocno() == 1) then
          agent_arrays(i)%output_unit=free_unit()
          open(unit=agent_arrays(i)%output_unit, file=trim(agent_array_name)//'.detectors', action="write")
          write(agent_arrays(i)%output_unit, '(a)') "<header>"

          call initialise_constant_diagnostics(agent_arrays(i)%output_unit, binary_format = agent_arrays(i)%binary_output)

          ! Initial columns are elapsed time and dt.
          column=1
          buffer=field_tag(name="ElapsedTime", column=column, statistic="value")
          write(agent_arrays(i)%output_unit, '(a)') trim(buffer)
          column=column+1
          buffer=field_tag(name="dt", column=column, statistic="value")
          write(agent_arrays(i)%output_unit, '(a)') trim(buffer)

          ! Next columns contain the positions of all the detectors.
          agent => agent_arrays(i)%firstnode
          positionloop: do j=1, agent_arrays(i)%length
             buffer=field_tag(name=agent%name, column=column+1, statistic="position",components=size(agent%position))
             write(agent_arrays(i)%output_unit, '(a)') trim(buffer)
             column=column+size(agent%position)
             agent => agent%next
          end do positionloop

          write(agent_arrays(i)%output_unit, '(a)') "</header>"
          flush(agent_arrays(i)%output_unit)
          close(agent_arrays(i)%output_unit)
       end if

       ! bit of hack to delete any existing .detectors.dat file
       ! if we don't delete the existing .detectors.dat would simply be opened for random access and 
       ! gradually overwritten, mixing detector output from the current with that of a previous run
       call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(agent_array_name)//'.detectors.dat', MPI_MODE_CREATE + MPI_MODE_RDWR + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, agent_arrays(i)%mpi_fh, ierror)
       call MPI_FILE_CLOSE(agent_arrays(i)%mpi_fh, ierror)
    
       call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(agent_array_name)//'.detectors.dat', MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, agent_arrays(i)%mpi_fh, ierror)
       assert(ierror == MPI_SUCCESS)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! Delete all agents whose element we don't own
       if (isparallel()) then
          agent => agent_arrays(i)%firstnode
          do while (associated(agent))
             if (element_owner(xfield%mesh,agent%element)/=getprocno() .or. agent%element<0) then 
                call delete(agent_arrays(i),agent)
             else
                agent%local=.true. 
                agent => agent%next
             end if
          end do
       end if

       ewrite(2,*) "Found", agent_arrays(i)%length, "local agents in array", i
    end do 

  end subroutine initialise_lagrangian_biology

  subroutine lagrangian_biology_cleanup()
    integer :: i, ierror

    do i = 1, size(agent_arrays)
       call delete_all(agent_arrays(i))
       if (agent_arrays(i)%mpi_fh/=0) then
          call MPI_FILE_CLOSE(agent_arrays(i)%mpi_fh, ierror) 
          if(ierror /= MPI_SUCCESS) then
             ewrite(0,*) "Warning: failed to close .detector file open with mpi_file_open"
          end if
       end if
    end do
    deallocate(agent_arrays)

  end subroutine lagrangian_biology_cleanup

  subroutine calculate_lagrangian_biology(state, time, dt, timestep)
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time, dt
    integer, intent(in) :: timestep

    integer :: i

    do i = 1, size(agent_arrays)
       ! Move lagrangian detectors
       if (check_any_lagrangian(agent_arrays(i))) then
          call move_lagrangian_detectors(state, agent_arrays(i), dt, timestep)
       end if

       call write_detectors(state, agent_arrays(i), time, dt)
    end do

  end subroutine calculate_lagrangian_biology

end module lagrangian_biology
