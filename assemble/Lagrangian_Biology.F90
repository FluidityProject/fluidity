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
  use global_parameters, only: PYTHON_FUNC_LEN, OPTION_PATH_LEN
  use parallel_tools
  use mpi_interfaces
  use pickers_inquire
  use detector_data_types
  use detector_tools
  use detector_move_lagrangian
  use detector_parallel
  use diagnostic_variables, only: initialise_constant_diagnostics, field_tag, &
                                  write_detectors, create_single_detector
  use python_state

implicit none

  private

  public :: initialise_lagrangian_biology, lagrangian_biology_cleanup, &
            calculate_lagrangian_biology

  type(detector_linked_list), dimension(:), allocatable, target, save :: agent_arrays

  integer, parameter :: BIOFIELD_NONE=0, BIOFIELD_DIAG=1, BIOFIELD_UPTAKE=2, BIOFIELD_RELEASE=3

contains

  subroutine initialise_lagrangian_biology(state)
    type(state_type), intent(in) :: state

    type(detector_type), pointer :: agent
    type(vector_field), pointer :: xfield
    type(element_type), pointer :: shape
    character(len=PYTHON_FUNC_LEN) :: func
    character(len=OPTION_PATH_LEN) :: buffer, fg_buffer, stage_buffer, biovar_buffer
    character(len=FIELD_NAME_LEN) :: field_name, fg_name, stage_name
    real, allocatable, dimension(:,:) :: coords
    real, allocatable, dimension(:) :: initial_values
    real:: current_time
    integer :: i, j, fg, dim, n_fgroups, n_agents, n_agent_arrays, n_fg_arrays, column, &
               ierror, det_type, random_seed, biovar_index, biovar_total, biovar_internal, &
               biovar_uptake, biovar_release, index

    if (.not.have_option("/ocean_biology/lagrangian_ensemble")) return

    ewrite(1,*) "In initialise_lagrangian_biology"

    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)
    xfield=>extract_vector_field(state, "Coordinate")
    shape=>ele_shape(xfield,1)

    ! Determine how many arrays we need across all functional groups
    n_agent_arrays = 0
    n_fgroups = option_count("/ocean_biology/lagrangian_ensemble/functional_group")
    do fg=1, n_fgroups
       write(fg_buffer, "(a,i0,a)") "/ocean_biology/lagrangian_ensemble/functional_group[",i-1,"]"

       n_agent_arrays = n_agent_arrays + option_count(trim(fg_buffer)//"/stage_array")
    end do

    allocate(agent_arrays(n_agent_arrays))
    ewrite(2,*) "Found a total of ", n_agent_arrays, " agent arrays in ", n_fgroups, " functional groups"

    ! We create an agent array for each stage of each Functional Group
    do fg=1, n_fgroups
       write(fg_buffer, "(a,i0,a)") "/ocean_biology/lagrangian_ensemble/functional_group[",i-1,"]"
       call get_option(trim(fg_buffer)//"/name", fg_name)
       n_fg_arrays = option_count(trim(fg_buffer)//"/stage_array")

       do i = 1, n_fg_arrays
          write(stage_buffer, "(a,i0,a)") trim(fg_buffer)//"/stage_array[",i-1,"]"
          call get_option(trim(stage_buffer)//"/number_of_agents", n_agents)
          call get_option(trim(stage_buffer)//"/name", stage_name)
          agent_arrays(i)%name = trim(fg_name)//trim(stage_name)
          agent_arrays(i)%fg_id = fg

          ! Register the agent array, so Zoltan/Adaptivity will not forget about it
          call register_detector_list(agent_arrays(i))
          agent_arrays(i)%total_num_det=n_agents

          ! Get options for lagrangian detector movement
          call read_detector_move_options(agent_arrays(i),trim(fg_buffer))

          ! Get options for Random Walk
          if (have_option(trim(stage_buffer)//"/random_walk")) then
             agent_arrays(i)%move_parameters%do_random_walk=.true.
             call get_option(trim(stage_buffer)//"/random_walk/python", agent_arrays(i)%move_parameters%rw_pycode)
             call get_option(trim(stage_buffer)//"/random_walk/random_seed", random_seed)
          
             ! Initialise random number generator
             call python_run_string("numpy.random.seed("//trim(int2str(random_seed))//")")
          else
             agent_arrays(i)%move_parameters%do_random_walk=.false.
          end if
          agent_arrays(i)%total_num_det=n_agents

          ! Collect other meta-information
          if (have_option(trim(stage_buffer)//"/binary_output")) then
             agent_arrays(i)%binary_output=.true.
          else
             agent_arrays(i)%binary_output=.false.
          end if
          if (have_option(trim(stage_buffer)//"/lagrangian")) then
             det_type=LAGRANGIAN_DETECTOR
          else
             det_type=STATIC_DETECTOR
          end if
          if (have_option(trim(stage_buffer)//"/exclude_from_advection")) then
             agent_arrays(i)%move_parameters%do_velocity_advect=.false.
          else
             agent_arrays(i)%move_parameters%do_velocity_advect=.true.
          end if

          ! Create agent and insert into list
          call get_option(trim(stage_buffer)//"/initial_position", func)
          allocate(coords(dim,n_agents))
          call set_detector_coords_from_python(coords, n_agents, func, current_time)

          do j = 1, n_agents
             call create_single_detector(agent_arrays(i), xfield, coords(:,j), j, det_type, trim(int2str(j)))
          end do
          deallocate(coords)

          ewrite(2,*) "Found", agent_arrays(i)%length, "local agents in array", i

          ! Biology parameters
          biovar_internal = option_count(trim(fg_buffer)//"/variable") 
          biovar_uptake = option_count(trim(fg_buffer)//"/uptake_variable")
          biovar_release = option_count(trim(fg_buffer)//"/release_variable")
          biovar_total = biovar_internal + biovar_uptake + biovar_release
          if (biovar_total > 0) then

             ! Allocate variable arrays
             allocate(agent_arrays(i)%biovar_list(biovar_total))
             allocate(agent_arrays(i)%biofield_type(biovar_total))
             allocate(agent_arrays(i)%biofield_list(biovar_total))
             allocate(agent_arrays(i)%chemfield_list(biovar_total))
             allocate(initial_values(biovar_total))

             ! Add internal state variables
             index = 1
             do j=1, biovar_internal
                write(biovar_buffer, "(a,i0,a)") trim(fg_buffer)//"/variable[",j-1,"]"
                call get_option(trim(biovar_buffer)//"/name", agent_arrays(i)%biovar_list(index))
                call get_option(trim(biovar_buffer)//"/initial_value", initial_values(index))

                ! Record according diagnostic field
                if (have_option(trim(biovar_buffer)//"/scalar_field")) then
                   agent_arrays(i)%biofield_type(index) = BIOFIELD_DIAG
                   call get_option(trim(biovar_buffer)//"/scalar_field/name", field_name)
                   if (have_option(trim(biovar_buffer)//"/scalar_field/per_stage")) then
                      agent_arrays(i)%biofield_list(index) = trim(field_name)//trim(stage_name)
                   else
                      agent_arrays(i)%biofield_list(index) = trim(field_name)
                   end if
                else
                   agent_arrays(i)%biofield_type(index) = BIOFIELD_NONE
                end if
                index = index+1
             end do

             ! Add uptake variables
             do j=1, biovar_uptake
                write(biovar_buffer, "(a,i0,a)") trim(fg_buffer)//"/uptake_variable[",j-1,"]"
                call get_option(trim(biovar_buffer)//"/name", agent_arrays(i)%biovar_list(index))
                initial_values(index)=0.0
                agent_arrays(i)%biofield_type(index) = BIOFIELD_UPTAKE
                agent_arrays(i)%biofield_list(index) = trim(agent_arrays(i)%biovar_list(index))//"Request"
                call get_option(trim(biovar_buffer)//"/chemical_field/name", agent_arrays(i)%chemfield_list(index))

                index = index+1
             end do

             ! Add release variables
             do j=1, biovar_release
                write(biovar_buffer, "(a,i0,a)") trim(fg_buffer)//"/release_variable[",j-1,"]"
                call get_option(trim(biovar_buffer)//"/name", agent_arrays(i)%biovar_list(index))
                initial_values(index)=0.0
                agent_arrays(i)%biofield_type(index) = BIOFIELD_RELEASE
                agent_arrays(i)%biofield_list(index) = trim(agent_arrays(i)%biovar_list(index))//"Release"
                call get_option(trim(biovar_buffer)//"/chemical_field/name", agent_arrays(i)%chemfield_list(index))

                index = index+1
             end do

             ! Store the python update code
             if (have_option(trim(stage_buffer)//"/biology_update")) then
                call get_option(trim(stage_buffer)//"/biology_update", agent_arrays(i)%biovar_pycode)
                agent_arrays(i)%has_biology=.true.
             end if

             ! Initialise agent variables
             agent => agent_arrays(i)%first
             do while (associated(agent))
                allocate(agent%biology(biovar_total))
                do j=1, biovar_total
                   agent%biology(j) = initial_values(j)
                end do
                agent => agent%next
             end do
             deallocate(initial_values)
          end if

          ! Create simple position-only agent I/O header 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (getprocno() == 1) then
             agent_arrays(i)%output_unit=free_unit()
             open(unit=agent_arrays(i)%output_unit, file=trim(agent_arrays(i)%name)//'.detectors', action="write")
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
             positionloop: do j=1, n_agents
                buffer=field_tag(name=trim(int2str(j)), column=column+1, statistic="position",components=dim)
                write(agent_arrays(i)%output_unit, '(a)') trim(buffer)
                column=column+dim
             end do positionloop

             write(agent_arrays(i)%output_unit, '(a)') "</header>"
             flush(agent_arrays(i)%output_unit)
             close(agent_arrays(i)%output_unit)
          end if

          ! bit of hack to delete any existing .detectors.dat file
          ! if we don't delete the existing .detectors.dat would simply be opened for random access and 
          ! gradually overwritten, mixing detector output from the current with that of a previous run
          call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(agent_arrays(i)%name)//'.detectors.dat', MPI_MODE_CREATE + MPI_MODE_RDWR + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, agent_arrays(i)%mpi_fh, ierror)
          call MPI_FILE_CLOSE(agent_arrays(i)%mpi_fh, ierror)
    
          call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(agent_arrays(i)%name)//'.detectors.dat', MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, agent_arrays(i)%mpi_fh, ierror)
          assert(ierror == MPI_SUCCESS)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       end do 

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

  end subroutine lagrangian_biology_cleanup

  subroutine calculate_lagrangian_biology(state, time, dt, timestep)
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time, dt
    integer, intent(in) :: timestep

    type(detector_type), pointer :: agent
    type(scalar_field), pointer :: sfield, request_field, chemical_field, absorption_field, &
                                   depletion_field, release_field, source_field
    type(vector_field), pointer :: xfield
    integer :: i, j, n, current_fg

    ewrite(1,*) "In calculate_lagrangian_biology"

    xfield=>extract_vector_field(state(1), "Coordinate")

    do i = 1, size(agent_arrays)
       ! Move lagrangian detectors
       if (check_any_lagrangian(agent_arrays(i))) then
          call move_lagrangian_detectors(state, agent_arrays(i), dt, timestep)
       end if

       ! Prepare python state if the Random Walk didn't do so
       if (.not. agent_arrays(i)%move_parameters%do_random_walk) then
          call python_reset()
          call python_add_state(state(1))
       end if


       if (agent_arrays(i)%has_biology) then
          ! Compile python function to set bio-variable, and store in the global dictionary
          call python_run_detector_string(trim(agent_arrays(i)%biovar_pycode), trim(agent_arrays(i)%name), trim("biology_update"))

          ewrite(2,*) "Updating biology agents..."

          ! Update agent biology
          agent=>agent_arrays(i)%first
          do while (associated(agent))
             call python_calc_agent_biology(agent, xfield, dt, trim(agent_arrays(i)%name), trim("biology_update"))
             agent=>agent%next
          end do

          ! Reset the diagnostic fields associated with the agent array
          do j=1, size(agent_arrays(i)%biovar_list)
             if (agent_arrays(i)%biofield_type(j) /= BIOFIELD_NONE) then
                sfield=>extract_scalar_field(state(1), trim(agent_arrays(i)%biofield_list(j)))
                call zero(sfield)
             end if
          end do
       end if
    end do

    ! Calculate diagnostic fields from agents variables
    do i = 1, size(agent_arrays)
       if (agent_arrays(i)%has_biology) then
          do j=1, size(agent_arrays(i)%biovar_list)
             if (agent_arrays(i)%biofield_type(j) /= BIOFIELD_NONE) then
                sfield=>extract_scalar_field(state(1), trim(agent_arrays(i)%biofield_list(j)))
                call set_diagnostic_field_from_agents(agent_arrays(i), j, sfield)
             end if
          end do
       end if
    end do

    current_fg = 1
    do i = 1, size(agent_arrays)
       ! We only want to do this for each functional group
       if (agent_arrays(i)%has_biology.and.agent_arrays(i)%fg_id>=current_fg) then
          do j=1, size(agent_arrays(i)%biovar_list)
             ! Handle chemical uptake
             if (agent_arrays(i)%biofield_type(j) == BIOFIELD_UPTAKE) then
                request_field=>extract_scalar_field(state(1), trim(agent_arrays(i)%biofield_list(j)))
                chemical_field=>extract_scalar_field(state(1), trim(agent_arrays(i)%chemfield_list(j)))
                absorption_field=>extract_scalar_field(state(1), trim(agent_arrays(i)%chemfield_list(j))//"Absorption")
                depletion_field=>extract_scalar_field(state(1), trim(agent_arrays(i)%biovar_list(j))//"Depletion")
                call zero(depletion_field)

                do n=1, node_count(request_field)
                   if (node_val(request_field,n) > node_val(chemical_field,n)) then
                      ! Scale back the request
                      call set(depletion_field, n, node_val(chemical_field,n) / node_val(request_field,n))
                   else
                      call set(depletion_field, n, 1.0)
                   end if

                   call set(absorption_field, n, node_val(request_field,n) * node_val(depletion_field,n))
                end do
             end if

             ! Handle chemical release
             if (agent_arrays(i)%biofield_type(j) == BIOFIELD_RELEASE) then
                release_field=>extract_scalar_field(state(1), trim(agent_arrays(i)%biofield_list(j)))
                source_field=>extract_scalar_field(state(1), trim(agent_arrays(i)%chemfield_list(j))//"Source")
                do n=1, node_count(release_field)
                   call set(source_field, n, node_val(release_field,n))
                end do
             end if
          end do

          current_fg = current_fg + 1
       end if
    end do

    ! Output agent positions
    do i = 1, size(agent_arrays)
       call write_detectors(state, agent_arrays(i), time, dt)
    end do

  end subroutine calculate_lagrangian_biology

  subroutine set_diagnostic_field_from_agents(agent_list, biovar_id, sfield)
    type(detector_linked_list), intent(inout) :: agent_list
    integer, intent(in) :: biovar_id
    type(scalar_field), intent(inout) :: sfield

    type(detector_type), pointer :: agent
    type(element_type), pointer :: shape
    integer, dimension(ele_loc(sfield, agent_list%first%element)) :: element_nodes
    real :: value, scaled_value
    integer :: i

    ewrite(2,*) "In set_diagnostic_field_from_agents"

    agent => agent_list%first
    do while (associated(agent))
       value = agent%biology(biovar_id)
       shape => ele_shape(sfield, agent%element)
       element_nodes = ele_nodes(sfield, agent%element)
       do i=1, ele_loc(sfield, agent%element)
          scaled_value = value * eval_shape(shape, i, agent%local_coords)
          call addto(sfield, element_nodes(i), scaled_value)
       end do

       agent => agent%next
    end do

  end subroutine set_diagnostic_field_from_agents

end module lagrangian_biology
