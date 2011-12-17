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
  use detector_python
  use diagnostic_variables, only: initialise_constant_diagnostics, field_tag, &
                                  write_detectors, create_single_detector
  use python_state
  use diagnostic_fields
  use Profiler

implicit none

  private

  public :: initialise_lagrangian_biology, lagrangian_biology_cleanup, &
            update_lagrangian_biology, calculate_agent_diagnostics

  type(detector_linked_list), dimension(:), allocatable, target, save :: agent_arrays

  character(len=FIELD_NAME_LEN), dimension(:), pointer :: uptake_field_names, release_field_names

  integer, parameter :: BIOFIELD_NONE=0, BIOFIELD_DIAG=1, BIOFIELD_UPTAKE=2, BIOFIELD_RELEASE=3
  integer, parameter :: BIOVAR_STAGE=1, BIOVAR_SIZE=2

contains

  subroutine initialise_lagrangian_biology(state)
    type(state_type), intent(inout) :: state

    type(detector_type), pointer :: agent
    type(vector_field), pointer :: xfield
    type(element_type), pointer :: shape
    character(len=PYTHON_FUNC_LEN) :: func
    character(len=OPTION_PATH_LEN) :: buffer, fg_buffer, stage_buffer, biovar_buffer, &
                                      env_field_buffer
    character(len=FIELD_NAME_LEN) :: field_name, fg_name, stage_name, biovar_name
    real, allocatable, dimension(:,:) :: coords
    real:: current_time
    integer :: i, j, fg, dim, n_fgroups, n_agents, n_agent_arrays, n_fg_arrays, column, &
               ierror, det_type, biovar, array, biovar_total, biovar_state, rnd_seed_int, &
               biovar_chemical, biovar_uptake, biovar_release, rnd_dim, n_env_fields, chemvar_index
               
    integer, dimension(:), allocatable :: rnd_seed

    if (.not.have_option("/embedded_models/lagrangian_ensemble_biology")) return

    ewrite(1,*) "In initialise_lagrangian_biology"

    ! Initialise random number generator
    call random_seed(size=rnd_dim)
    allocate(rnd_seed(rnd_dim))
    call get_option("/embedded_models/lagrangian_ensemble_biology/random_seed", rnd_seed_int)
    do i=1, rnd_dim
       rnd_seed(i) = rnd_seed_int + i
    end do
    call random_seed(put=rnd_seed(1:rnd_dim))
    call python_run_string("numpy.random.seed("//trim(int2str(rnd_seed_int))//")")
    deallocate(rnd_seed)

    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)
    xfield=>extract_vector_field(state, "Coordinate")
    shape=>ele_shape(xfield,1)

    ! Determine how many arrays we need across all functional groups
    n_agent_arrays = 0
    n_fgroups = option_count("/embedded_models/lagrangian_ensemble_biology/functional_group")
    do fg=1, n_fgroups
       write(fg_buffer, "(a,i0,a)") "/embedded_models/lagrangian_ensemble_biology/functional_group[",fg-1,"]"
       n_agent_arrays = n_agent_arrays + option_count(trim(fg_buffer)//"/stages/stage")
    end do

    allocate(agent_arrays(n_agent_arrays))
    ewrite(2,*) "Found a total of ", n_agent_arrays, " agent arrays in ", n_fgroups, " functional groups"

    ! Create a persistent dictionary of FG variable and environment name mappings
    call python_run_string("persistent['fg_var_names'] = dict()")
    call python_run_string("persistent['fg_env_names'] = dict()")

    ! We create an agent array for each stage of each Functional Group
    array = 1
    do fg=1, n_fgroups
       write(fg_buffer, "(a,i0,a)") "/embedded_models/lagrangian_ensemble_biology/functional_group[",fg-1,"]"
       call get_option(trim(fg_buffer)//"/name", fg_name)
       n_fg_arrays = option_count(trim(fg_buffer)//"/stages/stage")

       do i = 1, n_fg_arrays
          write(stage_buffer, "(a,i0,a)") trim(fg_buffer)//"/stages/stage[",i-1,"]"
          call get_option(trim(stage_buffer)//"/number_of_agents", n_agents)
          call get_option(trim(stage_buffer)//"/name", stage_name)
          agent_arrays(array)%name = trim(fg_name)//trim(stage_name)
          agent_arrays(array)%fg_name = trim(fg_name)
          agent_arrays(array)%stage_name = trim(stage_name)
          agent_arrays(array)%fg_id = fg
          call get_option(trim(stage_buffer)//"/id", agent_arrays(array)%stage_id)

          ! Register the agent array, so Zoltan/Adaptivity will not forget about it
          call register_detector_list(agent_arrays(array))
          agent_arrays(array)%total_num_det=n_agents

          ! Get options for lagrangian detector movement
          call read_detector_move_options(agent_arrays(array),trim(fg_buffer))

          ! Get options for Random Walk
          call read_random_walk_options(agent_arrays(array), trim(stage_buffer))

          ! Set other meta-information
          agent_arrays(array)%binary_output=.true.
          det_type=LAGRANGIAN_DETECTOR
          if (have_option(trim(stage_buffer)//"/debug/exclude_from_advection")) then
             agent_arrays(array)%do_velocity_advect=.false.
          else
             agent_arrays(array)%do_velocity_advect=.true.
          end if

          ! Create agent and insert into list
          call get_option(trim(stage_buffer)//"/initial_position", func)
          allocate(coords(dim,n_agents))
          call set_detector_coords_from_python(coords, n_agents, func, current_time)

          do j = 1, n_agents
             call create_single_detector(agent_arrays(array), xfield, coords(:,j), j, det_type, trim(int2str(j)))
          end do
          deallocate(coords)

          ! Set agent%list_id, because we need it to detect stage changes
          agent=>agent_arrays(array)%first
          do while (associated(agent))
             agent%list_id=agent_arrays(array)%id
             agent=>agent%next
          end do

          ewrite(2,*) "Found", agent_arrays(array)%length, "local agents in array", i

          !!! Biology setup !!!

          ! First determine number of biovars
          biovar_state = option_count(trim(fg_buffer)//"/variables/state_variable") 
          biovar_chemical = option_count(trim(fg_buffer)//"/variables/chemical_variable")
          biovar_uptake = 0
          biovar_release = 0
          do j=1, biovar_chemical
             write(biovar_buffer, "(a,i0,a)") trim(fg_buffer)//"/variables/chemical_variable[",j-1,"]"
             if (have_option(trim(biovar_buffer)//"/uptake")) then
                biovar_uptake = biovar_uptake + 1
             end if
             if (have_option(trim(biovar_buffer)//"/release")) then
                biovar_release = biovar_release + 1
             end if
          end do
          biovar_total = biovar_state + biovar_chemical + biovar_uptake + biovar_release

          if (biovar_total > 0) then

             ! Allocate variable arrays
             allocate(agent_arrays(array)%biovars(biovar_total))

             ! Add internal state variables
             biovar = 1
             do j=1, biovar_state
                write(biovar_buffer, "(a,i0,a)") trim(fg_buffer)//"/variables/state_variable[",j-1,"]"
                call get_option(trim(biovar_buffer)//"/name", biovar_name)
                agent_arrays(array)%biovars(biovar)%name=trim(biovar_name)

                ! Record according diagnostic field
                if (have_option(trim(biovar_buffer)//"/scalar_field")) then
                   agent_arrays(array)%biovars(biovar)%field_type = BIOFIELD_DIAG
                   call get_option(trim(biovar_buffer)//"/scalar_field/name", field_name)                   
                   agent_arrays(array)%biovars(biovar)%field_name = trim(fg_name)//trim(field_name)//trim(biovar_name)
                   if (have_option(trim(biovar_buffer)//"/scalar_field/stage_aggregate")) then
                      agent_arrays(array)%biovars(biovar)%stage_aggregate=.true.
                   else
                      agent_arrays(array)%biovars(biovar)%stage_aggregate=.false.
                   end if
                else
                   agent_arrays(array)%biovars(biovar)%field_type = BIOFIELD_NONE
                end if
                biovar = biovar+1
             end do

             ! Add chemical variables
             do j=1, biovar_chemical
                write(biovar_buffer, "(a,i0,a)") trim(fg_buffer)//"/variables/chemical_variable[",j-1,"]"
                call get_option(trim(biovar_buffer)//"/name", biovar_name)

                ! Chemical pool variable and according diagnostic fields
                agent_arrays(array)%biovars(biovar)%name = trim(biovar_name)
                if (have_option(trim(biovar_buffer)//"/scalar_field")) then
                   agent_arrays(array)%biovars(biovar)%field_type = BIOFIELD_DIAG
                   call get_option(trim(biovar_buffer)//"/scalar_field/name", field_name)
                   agent_arrays(array)%biovars(biovar)%field_name = trim(fg_name)//trim(field_name)//trim(biovar_name)
                   if (have_option(trim(biovar_buffer)//"/scalar_field/stage_aggregate")) then
                      agent_arrays(array)%biovars(biovar)%stage_aggregate=.true.
                   else
                      agent_arrays(array)%biovars(biovar)%stage_aggregate=.false.
                   end if
                else
                   agent_arrays(array)%biovars(biovar)%field_type = BIOFIELD_NONE
                end if
                chemvar_index = biovar
                biovar = biovar+1

                ! Chemical uptake
                if (have_option(trim(biovar_buffer)//"/uptake")) then
                   if (.not.have_option(trim(biovar_buffer)//"/chemical_field")) then
                      FLExit("No chemical field specified for "//trim(biovar_name)//" uptake in functional group "//trim(fg_name))
                   end if

                   agent_arrays(array)%biovars(biovar)%name = trim(biovar_name)//"Uptake"
                   agent_arrays(array)%biovars(biovar)%field_type = BIOFIELD_UPTAKE
                   agent_arrays(array)%biovars(biovar)%field_name = trim(fg_name)//"Request"//trim(biovar_name)
                   call get_option(trim(biovar_buffer)//"/chemical_field/name", agent_arrays(array)%biovars(biovar)%chemfield)
                   call insert_global_uptake_field(agent_arrays(array)%biovars(biovar)%chemfield)

                   if (have_option(trim(biovar_buffer)//"/uptake/scalar_field::Request/stage_aggregate")) then
                      agent_arrays(array)%biovars(biovar)%stage_aggregate=.true.
                   else
                      agent_arrays(array)%biovars(biovar)%stage_aggregate=.false.
                   end if
                   agent_arrays(array)%biovars(biovar)%pool_index = chemvar_index
                   biovar = biovar+1
                end if

                ! Chemical release
                if (have_option(trim(biovar_buffer)//"/release")) then
                   if (.not.have_option(trim(biovar_buffer)//"/chemical_field")) then
                      FLExit("No chemical field specified for "//trim(biovar_name)//" release in functional group "//trim(fg_name))
                   end if

                   agent_arrays(array)%biovars(biovar)%name = trim(biovar_name)//"Release"
                   agent_arrays(array)%biovars(biovar)%field_type = BIOFIELD_RELEASE
                   agent_arrays(array)%biovars(biovar)%field_name = trim(fg_name)//"Release"//trim(biovar_name)
                   call get_option(trim(biovar_buffer)//"/chemical_field/name", agent_arrays(array)%biovars(biovar)%chemfield)
                   call insert_global_release_field(agent_arrays(array)%biovars(biovar)%chemfield)

                   if (have_option(trim(biovar_buffer)//"/release/scalar_field::Release/stage_aggregate")) then
                      agent_arrays(array)%biovars(biovar)%stage_aggregate=.true.
                   else
                      agent_arrays(array)%biovars(biovar)%stage_aggregate=.false.
                   end if
                   agent_arrays(array)%biovars(biovar)%pool_index = chemvar_index
                   biovar = biovar+1
                end if

             end do

             ! Store the python update code
             if (have_option(trim(stage_buffer)//"/biology_update")) then
                call get_option(trim(stage_buffer)//"/biology_update", agent_arrays(array)%biovar_pycode)
                agent_arrays(array)%has_biology=.true.
             end if

             ! Now we know the variable names for the FG, so we populate the mapping dict
             call python_run_string("persistent['fg_var_names']['"//trim(agent_arrays(array)%name)//"'] = []")
             do j=1, biovar_total
                call python_run_string("persistent['fg_var_names']['"//trim(agent_arrays(array)%name)//"'].append('"//trim(agent_arrays(array)%biovars(j)%name)//"')")
             end do

             ! Initialise agent variables
             call get_option(trim(stage_buffer)//"/initial_state", func)
             agent => agent_arrays(array)%first
             do while (associated(agent))
                allocate(agent%biology(biovar_total))
                call python_init_agent_biology(agent, agent_arrays(array), trim(func))
                agent => agent%next
             end do

             ! Record which environment fields to evaluate and pass to the update function
             n_env_fields = option_count(trim(fg_buffer)//"/environment/field")
             call python_run_string("persistent['fg_env_names']['"//trim(agent_arrays(array)%name)//"'] = []")
             allocate(agent_arrays(array)%env_field_name(n_env_fields))
             do j=1, n_env_fields
                write(env_field_buffer, "(a,i0,a)") trim(fg_buffer)//"/environment/field[",j-1,"]"
                call get_option(trim(env_field_buffer)//"/name", agent_arrays(array)%env_field_name(j))
                call python_run_string("persistent['fg_env_names']['"//trim(agent_arrays(array)%name)//"'].append('"//trim(agent_arrays(array)%env_field_name(j))//"')")
             end do

             ! Add Particle Management options
             if (have_option(trim(stage_buffer)//"/particle_management")) then
                agent_arrays(array)%do_particle_management = .true.
                call get_option(trim(stage_buffer)//"/particle_management/minimum", agent_arrays(array)%pm_min)
                call get_option(trim(stage_buffer)//"/particle_management/maximum", agent_arrays(array)%pm_max)
             end if
          end if

          ! Create simple position-only agent I/O header 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (getprocno() == 1) then
             agent_arrays(array)%output_unit=free_unit()
             open(unit=agent_arrays(array)%output_unit, file=trim(agent_arrays(array)%name)//'.detectors', action="write")
             write(agent_arrays(array)%output_unit, '(a)') "<header>"

             call initialise_constant_diagnostics(agent_arrays(array)%output_unit, binary_format = agent_arrays(array)%binary_output)

             ! Initial columns are elapsed time and dt.
             column=1
             buffer=field_tag(name="ElapsedTime", column=column, statistic="value")
             write(agent_arrays(array)%output_unit, '(a)') trim(buffer)
             column=column+1
             buffer=field_tag(name="dt", column=column, statistic="value")
             write(agent_arrays(array)%output_unit, '(a)') trim(buffer)

             ! Next columns contain the positions of all the detectors.
             positionloop: do j=1, n_agents
                buffer=field_tag(name=trim(int2str(j)), column=column+1, statistic="position",components=dim)
                write(agent_arrays(array)%output_unit, '(a)') trim(buffer)
                column=column+dim
             end do positionloop

             write(agent_arrays(array)%output_unit, '(a)') "</header>"
             flush(agent_arrays(array)%output_unit)
             close(agent_arrays(array)%output_unit)
          end if

          ! bit of hack to delete any existing .detectors.dat file
          ! if we don't delete the existing .detectors.dat would simply be opened for random access and 
          ! gradually overwritten, mixing detector output from the current with that of a previous run
          call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(agent_arrays(array)%name)//'.detectors.dat', MPI_MODE_CREATE + MPI_MODE_RDWR + MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, agent_arrays(array)%mpi_fh, ierror)
          call MPI_FILE_CLOSE(agent_arrays(array)%mpi_fh, ierror)    
          call MPI_FILE_OPEN(MPI_COMM_FEMTOOLS, trim(agent_arrays(array)%name)//'.detectors.dat', MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, agent_arrays(array)%mpi_fh, ierror)
          assert(ierror == MPI_SUCCESS)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          array = array + 1
       end do 

    end do

  end subroutine initialise_lagrangian_biology

  subroutine lagrangian_biology_cleanup(state)
    type(state_type), intent(inout) :: state

    integer :: i, ierror
    type(mesh_type) :: biology_mesh

    do i = 1, size(agent_arrays)
       call delete_all(agent_arrays(i))
       if (agent_arrays(i)%mpi_fh/=0) then
          call MPI_FILE_CLOSE(agent_arrays(i)%mpi_fh, ierror) 
          if(ierror /= MPI_SUCCESS) then
             ewrite(0,*) "Warning: failed to close .detector file open with mpi_file_open"
          end if
       end if
    end do

    biology_mesh = extract_mesh(state, "BiologyMesh")
    call deallocate(biology_mesh)

    deallocate(uptake_field_names)
    deallocate(release_field_names)

  end subroutine lagrangian_biology_cleanup

  subroutine update_lagrangian_biology(state, time, dt, timestep)
    type(state_type), dimension(:), intent(inout) :: state
    real, intent(in) :: time, dt
    integer, intent(in) :: timestep

    type(detector_type), pointer :: agent
    type(detector_linked_list) :: stage_change_list
    type(vector_field), pointer :: xfield
    integer :: i, j

    ewrite(1,*) "In calculate_lagrangian_biology"
    call profiler_tic("/update_lagrangian_biology")

    xfield=>extract_vector_field(state(1), "Coordinate")

    ! Prepare python-state
    call python_reset()
    call python_add_state(state(1))

    do i = 1, size(agent_arrays)
       ! Move lagrangian detectors
       if (check_any_lagrangian(agent_arrays(i))) then
          call move_lagrangian_detectors(state, agent_arrays(i), dt, timestep)
       end if

       if (agent_arrays(i)%has_biology) then
          ! Compile python function to set bio-variable, and store in the global dictionary
          call python_run_detector_string(trim(agent_arrays(i)%biovar_pycode), trim(agent_arrays(i)%name), trim("biology_update"))

          ewrite(2,*) "Updating biology agents..."

          ! Update agent biology
          agent=>agent_arrays(i)%first
          do while (associated(agent))
             call python_calc_agent_biology(agent, agent_arrays(i), xfield, state(1), dt, trim(agent_arrays(i)%name), trim("biology_update"))

             ! Check for stage change
             if (agent%biology(1) /= agent_arrays(i)%stage_id) then
                call move(agent, agent_arrays(i), stage_change_list)
             end if
             agent=>agent%next
          end do

       end if
    end do

    ewrite(2,*) "Handling stage changes..."

    ! Handle stage changes within FG
    agent=>stage_change_list%first
    do while (associated(agent))
       do j=1, size(agent_arrays)
          if (agent_arrays(j)%fg_id == agent_arrays(agent%list_id)%fg_id) then
             if (agent_arrays(j)%stage_id==agent%biology(1)) then
                call move(agent, stage_change_list, agent_arrays(j))
             end if
          end if
       end do
       agent=>agent%next
    end do

    call calculate_agent_diagnostics(state(1))

    call chemical_exchange(state(1))

    ewrite(2,*) "Particle Management..."

    ! Particle Management
    do i = 1, size(agent_arrays)
       if (agent_arrays(i)%do_particle_management) then
          call particle_management(agent_arrays(i), xfield)
       end if
    end do

    ! Output agent positions
    do i = 1, size(agent_arrays)
       call write_detectors(state, agent_arrays(i), time, dt)
    end do

    call profiler_toc("/update_lagrangian_biology")

  end subroutine update_lagrangian_biology

  subroutine calculate_agent_diagnostics(state)
    type(state_type), intent(inout) :: state

    type(scalar_field), pointer :: diagfield_stage, diagfield_agg, agent_count_field
    type(detector_type), pointer :: agent
    integer :: i, j, ele

    ewrite(1,*) "In calculate_agent_diagnostics"
    call profiler_tic("/calculate_agent_diagnostics")

    ! First we derive the per-stage diagnostic quantities from the agent data
    do i = 1, size(agent_arrays)
       if (agent_arrays(i)%has_biology) then
          call derive_per_stage_diagnostics(agent_arrays(i), state)
       end if
    end do

    ! Then we reset all aggregated diagnostics
    do i = 1, size(agent_arrays)
       if (agent_arrays(i)%has_biology) then
          do j=1, size(agent_arrays(i)%biovars)
             ! Reset stage-aggregated diagnostics
             if (agent_arrays(i)%biovars(j)%stage_aggregate.and.agent_arrays(i)%biovars(j)%field_type /= BIOFIELD_NONE) then
                diagfield_agg=>extract_scalar_field(state, trim(agent_arrays(i)%biovars(j)%field_name))
                call zero(diagfield_agg)
             end if

             ! Reset request global aggregate fields
             if (agent_arrays(i)%biovars(j)%field_type == BIOFIELD_UPTAKE) then
                diagfield_agg=>extract_scalar_field(state, trim(agent_arrays(i)%biovars(j)%chemfield)//"Request")
                call zero(diagfield_agg)
             end if

             ! Reset release global aggregate fields
             if (agent_arrays(i)%biovars(j)%field_type == BIOFIELD_RELEASE) then
                diagfield_agg=>extract_scalar_field(state, trim(agent_arrays(i)%biovars(j)%chemfield)//"Release")
                call zero(diagfield_agg)
             end if
          end do
       end if
    end do

    ! Now we derive all aggregated quantities
    do i = 1, size(agent_arrays)
       if (agent_arrays(i)%has_biology) then
          do j=1, size(agent_arrays(i)%biovars)

             ! Aggregate stage-aggregated diagnostic fields
             if (agent_arrays(i)%biovars(j)%stage_aggregate.and.agent_arrays(i)%biovars(j)%field_type /= BIOFIELD_NONE) then
                diagfield_stage=>extract_scalar_field(state, trim(agent_arrays(i)%biovars(j)%field_name)//trim(agent_arrays(i)%stage_name))
                diagfield_agg=>extract_scalar_field(state, trim(agent_arrays(i)%biovars(j)%field_name))

                do ele=1, ele_count(diagfield_agg)
                   call addto(diagfield_agg, ele, node_val(diagfield_stage, ele))
                end do
             end if

             ! Aggregate chemical uptake request
             if (agent_arrays(i)%biovars(j)%field_type == BIOFIELD_UPTAKE) then
                diagfield_stage=>extract_scalar_field(state, trim(agent_arrays(i)%biovars(j)%field_name)//trim(agent_arrays(i)%stage_name))
                diagfield_agg=>extract_scalar_field(state, trim(agent_arrays(i)%biovars(j)%chemfield)//"Request")

                do ele=1, ele_count(diagfield_agg)
                   call addto(diagfield_agg, ele, node_val(diagfield_stage, ele))
                end do
             end if

             ! Aggregate chemical release quantity
             if (agent_arrays(i)%biovars(j)%field_type == BIOFIELD_RELEASE) then
                diagfield_stage=>extract_scalar_field(state, trim(agent_arrays(i)%biovars(j)%field_name)//trim(agent_arrays(i)%stage_name))
                diagfield_agg=>extract_scalar_field(state, trim(agent_arrays(i)%biovars(j)%chemfield)//"Release")

                do ele=1, ele_count(diagfield_agg)
                   call addto(diagfield_agg, ele, node_val(diagfield_stage, ele))
                end do
             end if
          end do
       end if
    end do

    call profiler_toc("/calculate_agent_diagnostics")
    ewrite(2,*) "Exiting calculate_agent_diagnostics"

  end subroutine calculate_agent_diagnostics

  subroutine derive_per_stage_diagnostics(agent_list, state)
    ! Set per-stage diagnostic fields from agent variables, 
    ! including agent counts and chemical request/release fields
    type(detector_linked_list), intent(inout) :: agent_list
    type(state_type), intent(inout) :: state

    type(vector_field), pointer :: xfield
    type(scalar_field), dimension(size(agent_list%biovars)) :: diagfields
    type(scalar_field), pointer :: agent_count_field
    type(detector_type), pointer :: agent
    integer :: i
    real :: ele_volume

    call profiler_tic(trim(agent_list%name)//"::derive_per_stage_diagnostics")

    xfield=>extract_vector_field(state, "Coordinate")

    ! Pull and reset agent count field
    agent_count_field=>extract_scalar_field(state, trim(agent_list%fg_name)//"Agents"//trim(agent_list%stage_name))
    call zero(agent_count_field)

    ! Pull and reset all per-stage-array diagnostic fields
    do i=1, size(agent_list%biovars)
       if (agent_list%biovars(i)%field_type /= BIOFIELD_NONE) then
          diagfields(i)=extract_scalar_field(state, trim(agent_list%biovars(i)%field_name)//trim(agent_list%stage_name))
          call zero(diagfields(i))
       end if
    end do

    ! Loop over all agents in this list and determine element volume
    agent => agent_list%first
    do while (associated(agent))

       ! Increase agent density field
       ele_volume = element_volume(xfield, agent%element)
       call addto(agent_count_field, agent%element, 1.0/ele_volume)

       ! Add diagnostic quantities to field for all variables
       do i=1, size(agent_list%biovars)
          if (agent_list%biovars(i)%field_type /= BIOFIELD_NONE) then
             if (i == BIOVAR_SIZE) then
                ! Don't multiply size by size
                call addto(diagfields(i), agent%element, agent%biology(i))
             else
                ! Any other quantity gets scales by size(plankters per agent) and element volume
                call addto(diagfields(i), agent%element, agent%biology(i)*agent%biology(BIOVAR_SIZE))
             end if
          end if
       end do

       agent => agent%next
    end do

    call profiler_toc(trim(agent_list%name)//"::derive_per_stage_diagnostics")

  end subroutine derive_per_stage_diagnostics

  subroutine chemical_exchange(state)
    ! Handle uptake and release of chemicals
    type(state_type), intent(inout) :: state

    type(scalar_field), pointer :: request_field, chemfield, depletion_field
    type(vector_field), pointer :: xfield
    type(detector_type), pointer :: agent
    integer :: i, j, n, ele, poolvar
    integer, dimension(:), pointer :: element_nodes
    real :: chemval, chem_integral, request, depletion, ingested_amount

    call profiler_tic("/chemical_exchange")

    xfield=>extract_vector_field(state, "Coordinate")

    ! Modify quantities on the chemical fields
    do i = 1, size(uptake_field_names)

       request_field=>extract_scalar_field(state, trim(uptake_field_names(i))//"Request")
       chemfield=>extract_scalar_field(state, trim(uptake_field_names(i)))
       depletion_field=>extract_scalar_field(state, trim(uptake_field_names(i))//"Depletion")

       ! Loop over all elements in chemical fields
       do ele=1,ele_count(chemfield)
          chem_integral = integral_element(chemfield, xfield, ele)
          request = node_val(request_field, ele)

          ! Derive depletion factor
          if (request > chem_integral) then 
             depletion = chem_integral / request
          else
             depletion = 1.0
          end if
          call set(depletion_field, ele, depletion)

          ! Set new chemical concentration
          element_nodes=>ele_nodes(chemfield, ele)
          do n=1, size(element_nodes)
             chemval = node_val(chemfield, element_nodes(n))
             ! Avoid div-by-zero error
             if (chem_integral /= 0.0) then 
                call set(chemfield, element_nodes(n), chemval * (1.0 - (request * depletion / chem_integral) ))
             else
                call set(chemfield, element_nodes(n), 0.0)
             end if
          end do
       end do
    end do

    do i = 1, size(agent_arrays)
       if (agent_arrays(i)%has_biology) then
          do j=1, size(agent_arrays(i)%biovars)
             ! Adjust agent pools
             if (agent_arrays(i)%biovars(j)%field_type == BIOFIELD_UPTAKE) then

                depletion_field=>extract_scalar_field(state, trim(agent_arrays(i)%biovars(j)%chemfield)//"Depletion")
                poolvar = agent_arrays(i)%biovars(j)%pool_index

                agent => agent_arrays(i)%first
                do while (associated(agent))
                   ingested_amount = agent%biology(j) * node_val(depletion_field, agent%element) / agent%biology(BIOVAR_SIZE)
                   agent%biology(poolvar) = agent%biology(poolvar) + ingested_amount

                   agent => agent%next
                end do
             end if
          end do

          ! We need to re-do stage-aggregation as well...
          call derive_per_stage_diagnostics(agent_arrays(i), state)
       end if
    end do

    call profiler_toc("/chemical_exchange")

  end subroutine chemical_exchange

  subroutine particle_management(agent_list, xfield)
    type(detector_linked_list), intent(inout) :: agent_list
    type(vector_field), pointer, intent(inout) :: xfield

    type(detector_type), pointer :: agent
    type(detector_type), pointer :: new_agent
    type(element_type), pointer :: shape
    integer :: i

    shape=>ele_shape(xfield,1)

    agent => agent_list%first
    do while (associated(agent))
       if (agent%biology(2)>agent_list%pm_max) then
          ewrite(2,*) "ml805 agent over PM limit detected. Splitting..."

          agent_list%total_num_det=agent_list%total_num_det + 1

          ! Allocate and insert detector
          new_agent=>null()
          call allocate(new_agent, agent)
          call insert(new_agent,agent_list)

          ! Populate new agent
          new_agent%id_number=agent_list%total_num_det
          new_agent%name=trim(int2str(agent_list%total_num_det))
          new_agent%position=agent%position
          new_agent%element=agent%element
          new_agent%local_coords=agent%local_coords
          new_agent%type=agent%type
          new_agent%list_id=agent%list_id

          ! Populate new agent's biology
          allocate(new_agent%biology(size(agent%biology)))
          do i=1, size(agent%biology)
             new_agent%biology(i)=agent%biology(i)
          end do

          agent%biology(2)=agent%biology(2) / 2.0
          new_agent%biology(2)=new_agent%biology(2) / 2.0
       end if

       if (agent%biology(2)<agent_list%pm_min) then
          ewrite(2,*) "ml805 agent under PM limit detected"
       end if

       agent => agent%next
    end do

  end subroutine particle_management

  subroutine insert_global_uptake_field(field_name)
    character(len=FIELD_NAME_LEN), intent(in) :: field_name

    character(len=FIELD_NAME_LEN), dimension(:), pointer :: tmp_list
    integer :: i, old_size

    if (.not. associated(uptake_field_names)) then
       allocate(uptake_field_names(1))
       uptake_field_names(1) = trim(field_name)
    end if

    do i=1, size(uptake_field_names)
      if (trim(uptake_field_names(i)) == trim(field_name)) then
         return
      end if
    end do

    old_size=size(uptake_field_names)
    tmp_list=>uptake_field_names
    allocate(uptake_field_names(old_size+1))
    uptake_field_names(1:old_size)=tmp_list
    uptake_field_names(old_size+1)=trim(field_name)
    deallocate(tmp_list)

  end subroutine insert_global_uptake_field

  subroutine insert_global_release_field(field_name)
    character(len=FIELD_NAME_LEN), intent(in) :: field_name

    character(len=FIELD_NAME_LEN), dimension(:), pointer :: tmp_list
    integer :: i, old_size

    if (.not. associated(release_field_names)) then
       allocate(release_field_names(1))
       release_field_names(1) = trim(field_name)
    end if

    do i=1, size(release_field_names)
      if (trim(release_field_names(i)) == trim(field_name)) then
         return
      end if
    end do

    old_size=size(release_field_names)
    tmp_list=>release_field_names
    allocate(release_field_names(old_size+1))
    release_field_names(1:old_size)=tmp_list
    release_field_names(old_size+1)=trim(field_name)
    deallocate(tmp_list)

  end subroutine insert_global_release_field

end module lagrangian_biology
