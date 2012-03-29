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
  use python_state
  use diagnostic_fields
  use Profiler
  use integer_hash_table_module
  use LERM

implicit none

  private

  public :: initialise_lagrangian_biology_metamodel, initialise_lagrangian_biology_agents, &
            lagrangian_biology_cleanup, update_lagrangian_biology, calculate_agent_diagnostics, &
            get_num_functional_groups, get_functional_group
  public :: BIOFIELD_NONE, BIOFIELD_DIAG, BIOFIELD_UPTAKE, BIOFIELD_RELEASE, BIOFIELD_INGESTED

  type(detector_linked_list), dimension(:), allocatable, target, save :: agent_arrays
  type(functional_group), dimension(:), allocatable, target, save :: functional_groups

  character(len=FIELD_NAME_LEN), dimension(:), pointer :: uptake_field_names, release_field_names

  integer, parameter :: BIOFIELD_NONE=0, BIOFIELD_DIAG=1, BIOFIELD_UPTAKE=2, BIOFIELD_RELEASE=3, BIOFIELD_INGESTED=4
  integer, parameter :: BIOVAR_STAGE=1, BIOVAR_SIZE=2

contains

  function get_num_functional_groups()
    integer :: get_num_functional_groups
    get_num_functional_groups = size(functional_groups)
  end function get_num_functional_groups

  function get_functional_group(i) result(fg)
    integer, intent(in) :: i
    type(functional_group), pointer :: fg

    fg => functional_groups(i)
  end function get_functional_group

  subroutine initialise_lagrangian_biology_metamodel()
    character(len=OPTION_PATH_LEN) :: fg_buffer, stage_buffer, env_field_buffer
    character(len=FIELD_NAME_LEN) :: stage_name
    integer, dimension(:), allocatable :: rnd_seed
    integer :: i, rnd_seed_int, rnd_dim
    integer :: j, array, fg, stage, n_agent_arrays, n_fgroups, n_stages, n_env_fields

    if (.not.have_option("/embedded_models/lagrangian_ensemble_biology")) return

    ewrite(1,*) "Lagrangian biology: Initialising metamodel..."

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
    ewrite(2,*) "Lagrangian biology: Initialised RNG"

    ! Determine how many arrays we need across all functional groups
    n_agent_arrays = 0
    n_fgroups = option_count("/embedded_models/lagrangian_ensemble_biology/functional_group")
    allocate(functional_groups(n_fgroups))
    ewrite(2,*) "Lagrangian biology: Allocated ", n_fgroups, " functional groups"

    do fg=1, n_fgroups
       write(fg_buffer, "(a,i0,a)") "/embedded_models/lagrangian_ensemble_biology/functional_group[",fg-1,"]"
       n_agent_arrays = n_agent_arrays + option_count(trim(fg_buffer)//"/stages/stage")
    end do

    allocate(agent_arrays(n_agent_arrays))
    ewrite(2,*) "Lagrangian biology: Allocated ", n_agent_arrays, " agent arrays"

    ! Create a persistent dictionary of FG variable and environment name mappings
    call python_run_string("persistent['fg_var_names'] = dict()")
    call python_run_string("persistent['fg_env_names'] = dict()")

    ! We create an agent array for each stage of each Functional Group
    array = 1
    do fg=1, n_fgroups
       write(fg_buffer, "(a,i0,a)") "/embedded_models/lagrangian_ensemble_biology/functional_group[",fg-1,"]"
       n_stages = option_count(trim(fg_buffer)//"/stages/stage")

       ! Get biology meta-data
       if (have_option(trim(fg_buffer)//"/variables")) then
          call read_functional_group(functional_groups(fg), trim(fg_buffer))

       ! Get biology meta-data for LERM diatoms
       elseif (have_option(trim(fg_buffer)//"/variables_lerm")) then
          call LERM_get_diatom_fgroup(functional_groups(fg), trim(fg_buffer))

          allocate(uptake_field_names(3))
          uptake_field_names(1)="DissolvedAmmonium"
          uptake_field_names(2)="DissolvedNitrate"
          uptake_field_names(3)="DissolvedSilicate"

          allocate(release_field_names(2))
          release_field_names(1)="DissolvedAmmonium"
          release_field_names(2)="DissolvedSilicate"

       else
          FLExit("No variables defined for functional group under: "//trim(fg_buffer))
       end if

       do stage = 1, n_stages
          write(stage_buffer, "(a,i0,a)") trim(fg_buffer)//"/stages/stage[",stage-1,"]"
          agent_arrays(array)%stage_options = trim(stage_buffer)
          call get_option(trim(stage_buffer)//"/name", stage_name)
          agent_arrays(array)%fgroup => functional_groups(fg)

          agent_arrays(array)%name = trim(functional_groups(fg)%name)//trim(stage_name)
          agent_arrays(array)%stage_name = trim(stage_name)
          call get_option(trim(stage_buffer)//"/id", agent_arrays(array)%stage_id)
          agent_arrays(array)%id=array

          ! Register the agent array, so Zoltan/Adaptivity will not forget about it
          call register_detector_list(agent_arrays(array))

          ! Get options for lagrangian detector movement
          call read_detector_move_options(agent_arrays(array), trim(stage_buffer)//"/movement")

          ! Get options for Random Walk
          call read_random_walk_options(agent_arrays(array), trim(stage_buffer)//"/movement")

          ! Now we know the variable names for the FG, so we record them in Python
          call python_run_string("persistent['fg_var_names']['"//trim(agent_arrays(array)%name)//"'] = []")
          do j=1, size(functional_groups(fg)%variables)
             call python_run_string("persistent['fg_var_names']['"//trim(agent_arrays(array)%name)//"'].append('"//trim(functional_groups(fg)%variables(j)%name)//"')")
          end do

          ! Record which environment fields to evaluate and pass to the update function
          if (have_option(trim(fg_buffer)//"/environment")) then
             n_env_fields = option_count(trim(fg_buffer)//"/environment/field")
             call python_run_string("persistent['fg_env_names']['"//trim(agent_arrays(array)%name)//"'] = []")
             allocate(agent_arrays(array)%env_field_name(n_env_fields))
             do j=1, n_env_fields
                write(env_field_buffer, "(a,i0,a)") trim(fg_buffer)//"/environment/field[",j-1,"]"
                call get_option(trim(env_field_buffer)//"/name", agent_arrays(array)%env_field_name(j))
                call python_run_string("persistent['fg_env_names']['"//trim(agent_arrays(array)%name)//"'].append('"//trim(agent_arrays(array)%env_field_name(j))//"')")
             end do

          elseif (have_option(trim(fg_buffer)//"/environment_lerm")) then
             call LERM_get_diatom_env_fields(agent_arrays(array)%env_field_name)
          else
             FLExit("No environment fields defined for functional group "//trim(agent_arrays(array)%fgroup%name))
          end if

          ! Store the python update code
          if (have_option(trim(stage_buffer)//"/biology/python")) then
             call get_option(trim(stage_buffer)//"/biology/python", agent_arrays(array)%biovar_pycode)
          end if

          array = array + 1
       end do 
    end do

    ewrite(1,*) "Lagrangian biology: Initialised metamodel"

  end subroutine initialise_lagrangian_biology_metamodel

  subroutine initialise_lagrangian_biology_agents(state)
    type(state_type), dimension(:), intent(inout) :: state

    type(detector_type), pointer :: agent
    type(vector_field), pointer :: xfield
    character(len=PYTHON_FUNC_LEN) :: func
    character(len=OPTION_PATH_LEN) :: fg_buffer, stage_buffer
    character(len=FIELD_NAME_LEN) :: field_name, fg_name, stage_name
    real, allocatable, dimension(:,:) :: coords
    real:: current_time
    integer :: j, dim, n_fgroups, n_agents, n_fg_arrays, &
               ierror, det_type, array

    if (.not.have_option("/embedded_models/lagrangian_ensemble_biology")) return

    ewrite(1,*) "Lagrangian biology: Initialising agents..."

    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)
    xfield=>extract_vector_field(state(1), "Coordinate")

    ! Initialise the shared ID counter to generate agent IDs
    call init_id_counter()
    ewrite(2,*) "Lagrangian biology: Initialised ID counter"

    ! We create an agent array for each stage of each Functional Group
    do array=1, size(agent_arrays)
       stage_buffer = trim(agent_arrays(array)%stage_options)

       call get_option(trim(stage_buffer)//"/initial_state/number_of_agents", n_agents)
       agent_arrays(array)%total_num_det=n_agents

       ! Set other meta-information
       agent_arrays(array)%binary_output=.true.
       det_type=LAGRANGIAN_DETECTOR

       ! Create agent and insert into list
       call get_option(trim(stage_buffer)//"/initial_state/position", func)
       allocate(coords(dim,n_agents))
       call set_detector_coords_from_python(coords, n_agents, func, current_time)

       do j = 1, n_agents
          call create_single_detector(agent_arrays(array), xfield, coords(:,j), det_type)
       end do
       deallocate(coords)

       ! Set agent%list_id, because we need it to detect stage changes
       agent=>agent_arrays(array)%first
       do while (associated(agent))
          agent%list_id=agent_arrays(array)%id
          agent=>agent%next
       end do

       ewrite(2,*) "Lagrangian biology: Initialised ", agent_arrays(array)%length, "agents locally for ", trim(agent_arrays(array)%name)

       ! Initialise agent biology variables
       if (associated(agent_arrays(array)%fgroup)) then
          if (have_option(trim(stage_buffer)//"/initial_state/biology")) then
             call get_option(trim(stage_buffer)//"/initial_state/biology", func)
             agent => agent_arrays(array)%first
             do while (associated(agent))
                allocate(agent%biology(size(agent_arrays(array)%fgroup%variables)))
                call python_init_agent_biology(agent, agent_arrays(array), trim(func))
                agent => agent%next
             end do
          elseif (have_option(trim(stage_buffer)//"/initial_state/biology_lerm_living_diatom")) then
             agent => agent_arrays(array)%first
             do while (associated(agent))
                allocate(agent%biology(size(agent_arrays(array)%fgroup%variables)))
                call LERM_initialise_living_diatom(agent%biology)
                agent => agent%next
             end do
          end if
       end if

       call write_detector_header(state, agent_arrays(array))

    end do

    ewrite(1,*) "Lagrangian biology: Initialised agents"

  end subroutine initialise_lagrangian_biology_agents

  subroutine read_functional_group(fgroup, fg_path) 
    type(functional_group), intent(out) :: fgroup
    character(len=*), intent(in) :: fg_path

    character(len=OPTION_PATH_LEN) :: var_buffer, stage_buffer
    character(len=FIELD_NAME_LEN) :: biovar_name, field_name
    integer :: i, vars_state, vars_chem, vars_uptake, vars_release, &
               vars_total, var_index, chemvar_index, stage_count

    call get_option(trim(fg_path)//"/name", fgroup%name)

    ! Record all associated stage names
    stage_count = option_count(trim(fg_path)//"/stages/stage")
    if (stage_count>0) then
       allocate(fgroup%stage_names%ptr(stage_count))
       do i=1, stage_count
          write(stage_buffer, "(a,i0,a)") trim(fg_path)//"/stages/stage[",i-1,"]"
          call get_option(trim(stage_buffer)//"/name", fgroup%stage_names%ptr(i))
       end do
    end if

    ! Record the agent count field path
    if (have_option(trim(fg_path)//"/scalar_field::Agents")) then
       fgroup%agents_field_path = trim(fg_path)//"/scalar_field::Agents"
    else
       FLExit("No Agents field specified for functional group "//trim(fgroup%name))
    end if

    ! Determine number of variables
    vars_state = option_count(trim(fg_path)//"/variables/state_variable") 
    vars_chem = option_count(trim(fg_path)//"/variables/chemical_variable")
    vars_uptake = 0
    vars_release = 0
    do i=1, vars_chem
       write(var_buffer, "(a,i0,a)") trim(fg_path)//"/variables/chemical_variable[",i-1,"]"
       if (have_option(trim(var_buffer)//"/uptake")) then
          vars_uptake = vars_uptake + 2
       end if
       if (have_option(trim(var_buffer)//"/release")) then
          vars_release = vars_release + 1
       end if
    end do
    vars_total = vars_state + vars_chem + vars_uptake + vars_release

    if (vars_total > 0) then

       ! Allocate variable arrays
       allocate(fgroup%variables(vars_total))

       ! Add internal state variables
       var_index = 1
       do i=1, vars_state
          write(var_buffer, "(a,i0,a)") trim(fg_path)//"/variables/state_variable[",i-1,"]"
          call get_option(trim(var_buffer)//"/name", biovar_name)
          fgroup%variables(var_index)%name=trim(biovar_name)
          if (have_option(trim(var_buffer)//"/include_in_io")) then
             fgroup%variables(var_index)%write_to_file=.true.
          end if

          ! Record according diagnostic field
          if (have_option(trim(var_buffer)//"/scalar_field")) then
             fgroup%variables(var_index)%field_type = BIOFIELD_DIAG
             call get_option(trim(var_buffer)//"/scalar_field/name", field_name)                   
             fgroup%variables(var_index)%field_name = trim(fgroup%name)//trim(field_name)//trim(biovar_name)
             fgroup%variables(var_index)%field_path = trim(var_buffer)//"/scalar_field"
             if (have_option(trim(var_buffer)//"/scalar_field/stage_aggregate")) then
                fgroup%variables(var_index)%stage_aggregate=.true.
             else
                fgroup%variables(var_index)%stage_aggregate=.false.
             end if
          else
             fgroup%variables(var_index)%field_type = BIOFIELD_NONE
          end if
          var_index = var_index+1
       end do

       ! Add chemical variables
       do i=1, vars_chem
          write(var_buffer, "(a,i0,a)") trim(fg_path)//"/variables/chemical_variable[",i-1,"]"
          call get_option(trim(var_buffer)//"/name", biovar_name)
          if (have_option(trim(var_buffer)//"/include_in_io")) then
             fgroup%variables(var_index)%write_to_file=.true.
          end if

          ! Chemical pool variable and according diagnostic fields
          fgroup%variables(var_index)%name = trim(biovar_name)
          if (have_option(trim(var_buffer)//"/scalar_field")) then
             fgroup%variables(var_index)%field_type = BIOFIELD_DIAG
             call get_option(trim(var_buffer)//"/scalar_field/name", field_name)
             fgroup%variables(var_index)%field_name = trim(fgroup%name)//trim(field_name)//trim(biovar_name)
             fgroup%variables(var_index)%field_path = trim(var_buffer)//"/scalar_field"
             if (have_option(trim(var_buffer)//"/scalar_field/stage_aggregate")) then
                fgroup%variables(var_index)%stage_aggregate=.true.
             else
                fgroup%variables(var_index)%stage_aggregate=.false.
             end if
          else
             fgroup%variables(var_index)%field_type = BIOFIELD_NONE
          end if
          chemvar_index = var_index
          var_index = var_index+1

          ! Chemical uptake
          ! Note check for existence of Request and Depletion field, and record Depletion field
          if (have_option(trim(var_buffer)//"/uptake")) then
             if (.not.have_option(trim(var_buffer)//"/chemical_field")) then
                FLExit("No chemical field specified for "//trim(biovar_name)//" uptake in functional group "//trim(fgroup%name))
             end if

             fgroup%variables(var_index)%name = trim(biovar_name)//"Uptake"
             fgroup%variables(var_index)%field_type = BIOFIELD_UPTAKE
             fgroup%variables(var_index)%field_name = trim(fgroup%name)//"Request"//trim(biovar_name)
             fgroup%variables(var_index)%field_path = trim(var_buffer)//"/uptake/scalar_field::Request"
             fgroup%variables(var_index)%depletion_field_path = trim(var_buffer)//"/uptake/scalar_field::Depletion"
             call get_option(trim(var_buffer)//"/chemical_field/name", fgroup%variables(var_index)%chemfield)
             call insert_global_uptake_field(fgroup%variables(var_index)%chemfield)

             if (have_option(trim(var_buffer)//"/uptake/scalar_field::Request/stage_aggregate")) then
                fgroup%variables(var_index)%stage_aggregate=.true.
             else
                fgroup%variables(var_index)%stage_aggregate=.false.
             end if
             fgroup%variables(var_index)%pool_index = chemvar_index
             var_index = var_index+1

             ! Create a 'ChemIngested' variable
             fgroup%variables(var_index)%name = trim(biovar_name)//"Ingested"
             fgroup%variables(var_index)%field_type = BIOFIELD_INGESTED
             fgroup%variables(var_index)%field_name = trim(fgroup%name)//"Ingested"//trim(biovar_name)
             fgroup%variables(var_index)%field_path = trim(var_buffer)//"/uptake/scalar_field::Ingested"
             call get_option(trim(var_buffer)//"/chemical_field/name", fgroup%variables(var_index)%chemfield)

             if (have_option(trim(var_buffer)//"/uptake/scalar_field::Ingested/stage_aggregate")) then
                fgroup%variables(var_index)%stage_aggregate=.true.
             else
                fgroup%variables(var_index)%stage_aggregate=.false.
             end if
             fgroup%variables(var_index)%pool_index = chemvar_index
             fgroup%variables(var_index)%request_index = var_index - 1
             var_index = var_index+1
          end if

          ! Chemical release
          if (have_option(trim(var_buffer)//"/release")) then
             if (.not.have_option(trim(var_buffer)//"/chemical_field")) then
                FLExit("No chemical field specified for "//trim(biovar_name)//" release in functional group "//trim(fgroup%name))
             end if

             fgroup%variables(var_index)%name = trim(biovar_name)//"Release"
             fgroup%variables(var_index)%field_type = BIOFIELD_RELEASE
             fgroup%variables(var_index)%field_name = trim(fgroup%name)//"Release"//trim(biovar_name)
             fgroup%variables(var_index)%field_path = trim(var_buffer)//"/release/scalar_field::Release"
             call get_option(trim(var_buffer)//"/chemical_field/name", fgroup%variables(var_index)%chemfield)
             call insert_global_release_field(fgroup%variables(var_index)%chemfield)

             if (have_option(trim(var_buffer)//"/release/scalar_field::Release/stage_aggregate")) then
                fgroup%variables(var_index)%stage_aggregate=.true.
             else
                fgroup%variables(var_index)%stage_aggregate=.false.
             end if
             fgroup%variables(var_index)%pool_index = chemvar_index
             var_index = var_index+1
          end if

       end do
    end if

  end subroutine read_functional_group

  subroutine lagrangian_biology_cleanup(state)
    type(state_type), intent(inout) :: state

    integer :: i, ierror
    type(mesh_type) :: biology_mesh

    do i = 1, size(agent_arrays)
       call delete_all(agent_arrays(i))
       if (agent_arrays(i)%output_unit/=0) then
          call MPI_FILE_CLOSE(agent_arrays(i)%output_unit, ierror) 
          if(ierror /= MPI_SUCCESS) then
             ewrite(0,*) "Warning: failed to close .detector file open with mpi_file_open"
          end if
       end if
    end do

    biology_mesh = extract_mesh(state, "BiologyMesh")
    call deallocate(biology_mesh)

    if (associated(uptake_field_names)) then
       deallocate(uptake_field_names)
    end if
    if (associated(release_field_names)) then
       deallocate(release_field_names)
    end if

    call delete_id_counter()

  end subroutine lagrangian_biology_cleanup

  subroutine update_lagrangian_biology(state, time, dt, timestep)
    type(state_type), dimension(:), intent(inout) :: state
    real, intent(in) :: time, dt
    integer, intent(in) :: timestep

    type(detector_type), pointer :: agent, agent_to_move
    type(detector_linked_list) :: stage_change_list
    type(vector_field), pointer :: xfield
    integer :: i, j, pm_period
    logical :: python_update, lerm_living_update, lerm_dead_update

    ewrite(1,*) "Lagrangian biology: Updating agents..."
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

       if (have_option(trim(agent_arrays(i)%stage_options)//"/biology")) then

          ewrite(2,*) "Lagrangian biology: Updating ", trim(agent_arrays(i)%name)

          if (have_option(trim(agent_arrays(i)%stage_options)//"/biology/python")) then
             python_update=.true.
          else
             python_update=.false.
          end if
          if (have_option(trim(agent_arrays(i)%stage_options)//"/biology/lerm_living_diatom")) then
             lerm_living_update=.true.
          else
             lerm_living_update=.false.
          end if
          if (have_option(trim(agent_arrays(i)%stage_options)//"/biology/lerm_dead_diatom")) then
             lerm_dead_update=.true.
          else
             lerm_dead_update=.false.
          end if

          if (python_update) then
             ! Compile python function to set bio-variable, and store in the global dictionary
             call python_run_detector_string(trim(agent_arrays(i)%biovar_pycode), trim(agent_arrays(i)%name), trim("biology_update"))
          end if

          ! Update agent biology
          agent=>agent_arrays(i)%first
          do while (associated(agent))

             if (python_update) then
                call python_calc_agent_biology(agent, agent_arrays(i), state(1), dt, trim(agent_arrays(i)%name), trim("biology_update"))

             elseif (lerm_living_update) then
                call LERM_update_living_diatom(agent, agent_arrays(i), state(1), dt)
             elseif (lerm_dead_update) then
                call LERM_update_dead_diatom(agent, agent_arrays(i), state(1), dt)
             end if

             ! Check for stage change
             if (agent%biology(BIOVAR_STAGE) /= agent_arrays(i)%stage_id) then
                agent_to_move=>agent
                agent=>agent%next
                call move(agent_to_move, agent_arrays(i), stage_change_list)
             else
                agent=>agent%next
             end if
          end do

       end if
    end do

    ewrite(2,*) "Handling stage changes across all agent lists"

    ! Handle stage changes within FG
    agent=>stage_change_list%first
    stage_change_loop: do while (associated(agent))
       do j=1, size(agent_arrays)
          if (trim(agent_arrays(j)%fgroup%name) == trim(agent_arrays(agent%list_id)%fgroup%name)) then
             if (agent_arrays(j)%stage_id==agent%biology(BIOVAR_STAGE)) then
                agent_to_move=>agent
                agent=>agent%next
                call move(agent_to_move, stage_change_list, agent_arrays(j))
                cycle stage_change_loop
             end if
          end if
       end do
    end do stage_change_loop

    ! Derive the full set af eulerian diagnostic fields
    ! This includes per-stage and stage-aggregated, 
    ! as well as request/release quantities 
    call calculate_agent_diagnostics(state(1))

    ! Execute chemical uptake/release
    call chemical_uptake(state(1))
    call chemical_release(state(1))

    ! Particle Management
    do i = 1, size(agent_arrays)
       if (have_option(trim(agent_arrays(i)%stage_options)//"/particle_management")) then
          call get_option(trim(agent_arrays(i)%stage_options)//"/particle_management/period_in_timesteps", pm_period)
          if (timestep == 1 .or. mod(timestep, pm_period) == 0) then
             call particle_management(state(1), agent_arrays(i))
          end if
       end if
    end do

    ! Re-derive the required agent diagnostic variables
    call calculate_agent_diagnostics(state(1))

    ! Output agent positions
    do i = 1, size(agent_arrays)
       call write_detectors(state, agent_arrays(i), time, dt, timestep)
    end do

    call profiler_toc("/update_lagrangian_biology")

  end subroutine update_lagrangian_biology

  subroutine calculate_agent_diagnostics(state)
    type(state_type), intent(inout) :: state

    type(scalar_field), pointer :: diagfield_stage, diagfield_agg, agent_count_field
    type(detector_type), pointer :: agent
    integer :: i, j, ele

    ewrite(1,*) "In calculate_agent_diagnostics"

    ! First we derive the per-stage diagnostic quantities from the agent data
    do i = 1, size(agent_arrays)
       if (have_option(trim(agent_arrays(i)%stage_options)//"/biology")) then
          call derive_per_stage_diagnostics(agent_arrays(i), state)
       end if
    end do

    ! Then we reset all aggregated diagnostics
    do i = 1, size(agent_arrays)
       call profiler_tic(trim(agent_arrays(i)%name)//"::diagnostics")

       if (have_option(trim(agent_arrays(i)%stage_options)//"/biology")) then
          do j=1, size(agent_arrays(i)%fgroup%variables)
             ! Reset stage-aggregated diagnostics
             if (agent_arrays(i)%fgroup%variables(j)%stage_aggregate.and.agent_arrays(i)%fgroup%variables(j)%field_type /= BIOFIELD_NONE) then
                diagfield_agg=>extract_scalar_field(state, trim(agent_arrays(i)%fgroup%variables(j)%field_name))
                call zero(diagfield_agg)
             end if
          end do
       end if

       call profiler_toc(trim(agent_arrays(i)%name)//"::diagnostics")
    end do

    ! Reset global request fields
    if (associated(uptake_field_names)) then
       do i = 1, size(uptake_field_names)
          diagfield_agg=>extract_scalar_field(state, trim(uptake_field_names(i))//"Request")
          call zero(diagfield_agg)
       end do
    end if

    ! Reset global release fields
    if (associated(release_field_names)) then
       do i = 1, size(release_field_names)
          diagfield_agg=>extract_scalar_field(state, trim(release_field_names(i))//"Release")
          call zero(diagfield_agg)
       end do
    end if

    ! Now we derive all aggregated quantities
    do i = 1, size(agent_arrays)
       call profiler_tic(trim(agent_arrays(i)%name)//"::diagnostics")

       if (have_option(trim(agent_arrays(i)%stage_options)//"/biology")) then
          do j=1, size(agent_arrays(i)%fgroup%variables)

             ! Aggregate stage-aggregated diagnostic fields
             if (agent_arrays(i)%fgroup%variables(j)%stage_aggregate.and.agent_arrays(i)%fgroup%variables(j)%field_type /= BIOFIELD_NONE) then
                diagfield_stage=>extract_scalar_field(state, trim(agent_arrays(i)%fgroup%variables(j)%field_name)//trim(agent_arrays(i)%stage_name))
                diagfield_agg=>extract_scalar_field(state, trim(agent_arrays(i)%fgroup%variables(j)%field_name))

                do ele=1, ele_count(diagfield_agg)
                   call addto(diagfield_agg, ele_nodes(diagfield_stage, ele), ele_val(diagfield_stage, ele))
                end do
             end if

             ! Aggregate chemical uptake request
             if (agent_arrays(i)%fgroup%variables(j)%field_type == BIOFIELD_UPTAKE) then
                diagfield_stage=>extract_scalar_field(state, trim(agent_arrays(i)%fgroup%variables(j)%field_name)//trim(agent_arrays(i)%stage_name))
                diagfield_agg=>extract_scalar_field(state, trim(agent_arrays(i)%fgroup%variables(j)%chemfield)//"Request")

                do ele=1, ele_count(diagfield_agg)
                   call addto(diagfield_agg, ele_nodes(diagfield_stage, ele), ele_val(diagfield_stage, ele))
                end do
             end if

             ! Aggregate chemical release quantity
             if (agent_arrays(i)%fgroup%variables(j)%field_type == BIOFIELD_RELEASE) then
                diagfield_stage=>extract_scalar_field(state, trim(agent_arrays(i)%fgroup%variables(j)%field_name)//trim(agent_arrays(i)%stage_name))
                diagfield_agg=>extract_scalar_field(state, trim(agent_arrays(i)%fgroup%variables(j)%chemfield)//"Release")

                do ele=1, ele_count(diagfield_agg)
                   call addto(diagfield_agg, ele_nodes(diagfield_stage, ele), ele_val(diagfield_stage, ele))
                end do
             end if
          end do
       end if

       call profiler_toc(trim(agent_arrays(i)%name)//"::diagnostics")
    end do

    ewrite(2,*) "Exiting calculate_agent_diagnostics"

  end subroutine calculate_agent_diagnostics

  subroutine derive_per_stage_diagnostics(agent_list, state)
    ! Set per-stage diagnostic fields from agent variables, 
    ! including agent counts and chemical request/release fields
    type(detector_linked_list), intent(inout) :: agent_list
    type(state_type), intent(inout) :: state

    type(vector_field), pointer :: xfield
    type(scalar_field), dimension(size(agent_list%fgroup%variables)) :: diagfields
    type(scalar_field), pointer :: agent_count_field
    type(detector_type), pointer :: agent
    integer :: i
    real :: ele_volume, release_amount

    call profiler_tic(trim(agent_list%name)//"::diagnostics")

    xfield=>extract_vector_field(state, "Coordinate")

    ! Pull and reset agent density field
    agent_count_field=>extract_scalar_field(state, trim(agent_list%fgroup%name)//"Agents"//trim(agent_list%stage_name))
    call zero(agent_count_field)

    ! Pull and reset all per-stage-array diagnostic fields
    do i=1, size(agent_list%fgroup%variables)
       if (agent_list%fgroup%variables(i)%field_type /= BIOFIELD_NONE) then
          diagfields(i)=extract_scalar_field(state, trim(agent_list%fgroup%variables(i)%field_name)//trim(agent_list%stage_name))
          call zero(diagfields(i))
       end if
    end do

    agent => agent_list%first
    do while (associated(agent))
       ele_volume = element_volume(xfield, agent%element)

       ! Increase agent density field
       call addto(agent_count_field, agent%element, 1.0/ele_volume)

       ! Add diagnostic quantities to the field for all variables
       do i=1, size(agent_list%fgroup%variables)

          ! All diagnostic agent quantities get divided by element volume
          ! and multiplied by the number of individuals an agent represents
          if (agent_list%fgroup%variables(i)%field_type == BIOFIELD_DIAG) then

             ! Agent size (number of individuals represented) does not get multiplied by itself
             if (i == BIOVAR_SIZE) then
                call addto(diagfields(i), agent%element, agent%biology(i) / ele_volume)
             else
                call addto(diagfields(i), agent%element, agent%biology(i)*agent%biology(BIOVAR_SIZE) / ele_volume)
             end if

          ! Uptake request and release are total amounts, so don't divide by element volume
          elseif (agent_list%fgroup%variables(i)%field_type == BIOFIELD_UPTAKE) then 
                call addto(diagfields(i), agent%element, agent%biology(i)*agent%biology(BIOVAR_SIZE))
          elseif (agent_list%fgroup%variables(i)%field_type == BIOFIELD_RELEASE) then
                ! Make sure we don't release more than we have
                !release_amount = min(agent%biology(i), agent%biology(agent_list%fgroup%variables(i)%pool_index))
                call addto(diagfields(i), agent%element, agent%biology(i)*agent%biology(BIOVAR_SIZE))
          end if
       end do

       agent => agent%next
    end do

    call profiler_toc(trim(agent_list%name)//"::diagnostics")

  end subroutine derive_per_stage_diagnostics

  subroutine chemical_uptake(state)
    ! Handle uptake and release of chemicals
    type(state_type), intent(inout) :: state

    type(scalar_field), pointer :: request_field, chemfield, depletion_field
    type(vector_field), pointer :: xfield
    type(detector_type), pointer :: agent
    integer :: i, j, n, ele
    integer, dimension(:), pointer :: element_nodes
    real :: chemval, chemval_new, chem_integral
    real, dimension(1) :: request, depletion

    ! Exit if there are no uptake fields
    if (.not.associated(uptake_field_names)) return

    call profiler_tic("/chemical_exchange")

    ewrite(2,*) "In chemical_uptake"

    xfield=>extract_vector_field(state, "Coordinate")

    ! Modify quantities on the chemical fields
    do i = 1, size(uptake_field_names)

       request_field=>extract_scalar_field(state, trim(uptake_field_names(i))//"Request")
       chemfield=>extract_scalar_field(state, trim(uptake_field_names(i)))
       depletion_field=>extract_scalar_field(state, trim(uptake_field_names(i))//"Depletion")

       ! Loop over all elements in chemical fields
       do ele=1,ele_count(chemfield)
          chem_integral = integral_element(chemfield, xfield, ele)
          request = ele_val(request_field, ele)

          ! Derive depletion factor
          if (request(1) > chem_integral .and. request(1) > 0.0) then 
             depletion(1) = chem_integral / request(1)
          else
             depletion(1) = 1.0
          end if
          call set(depletion_field, ele, depletion(1))

          ! Set new chemical concentration
          element_nodes=>ele_nodes(chemfield, ele)
          do n=1, size(element_nodes)
             chemval = node_val(chemfield, element_nodes(n))
             ! Avoid div-by-zero errors, and make sure not to write 0.0 into any field,
             ! since it will cause NaNs in the solvers. Instead use tiny()
             if (chem_integral > 0.0) then 
                ! C := C (1 - S/C_bar )
                chemval_new = chemval * (1.0 - (request(1) * depletion(1) / chem_integral))
                if (chemval_new > 0.0) then
                   call set(chemfield, element_nodes(n), chemval_new)
                else
                   call set(chemfield, element_nodes(n), tiny(1.0))
                end if
             else
                call set(chemfield, element_nodes(n), tiny(1.0))
             end if
          end do
       end do
    end do

    ! Apply depletion factor to uptake variables
    do i = 1, size(agent_arrays)
       if (have_option(trim(agent_arrays(i)%stage_options)//"/biology")) then
          do j=1, size(agent_arrays(i)%fgroup%variables)

             if (agent_arrays(i)%fgroup%variables(j)%field_type == BIOFIELD_INGESTED) then
                depletion_field=>extract_scalar_field(state, trim(agent_arrays(i)%fgroup%variables(j)%chemfield)//"Depletion")

                agent => agent_arrays(i)%first
                do while (associated(agent))
                   depletion = ele_val(depletion_field, agent%element)
                   agent%biology(j) = depletion(1) * agent%biology(agent_arrays(i)%fgroup%variables(j)%request_index)

                   agent => agent%next
                end do
             end if
          end do
       end if
    end do

    call profiler_toc("/chemical_exchange")

  end subroutine chemical_uptake

  subroutine chemical_release(state)
    ! Handle chemical release
    type(state_type), intent(inout) :: state

    type(scalar_field), pointer :: release_field, chemfield
    type(vector_field), pointer :: xfield
    type(detector_type), pointer :: agent
    integer :: i, j, n, ele, poolvar
    integer, dimension(:), pointer :: element_nodes
    real :: chemval, chemval_new, chem_integral, ele_volume
    real, dimension(1) :: release

    ! Exit if there are no release fields
    if (.not.associated(release_field_names)) return

    call profiler_tic("/chemical_exchange")

    ewrite(2,*) "In chemical_release"

    xfield=>extract_vector_field(state, "Coordinate")

    ! Modify quantities on the chemical fields
    do i = 1, size(release_field_names)

       release_field=>extract_scalar_field(state, trim(release_field_names(i))//"Release")
       chemfield=>extract_scalar_field(state, trim(release_field_names(i)))

       ! Loop over all elements in chemical fields
       do ele=1,ele_count(chemfield)
          chem_integral = integral_element(chemfield, xfield, ele)
          release = ele_val(release_field, ele)

          ! Set new chemical concentration
          element_nodes=>ele_nodes(chemfield, ele)
          do n=1, size(element_nodes)
             chemval = node_val(chemfield, element_nodes(n))
             ! Avoid div-by-zero error, and make sure not to write 0.0 into any field,
             ! since it will cause NaNs in the solvers. Instead use tiny()
             if (chem_integral > 0.0) then 
                ! C := C (1 + S/C_bar )
                chemval_new = chemval * (1.0 + (release(1) / chem_integral))
                if (chemval_new > 0.0) then
                   call set(chemfield, element_nodes(n), chemval_new)
                else
                   call set(chemfield, element_nodes(n), tiny(1.0))
                end if
             else
                ! If there is zero chemical in this element 
                ! we turn the total quantity into an evenly distributed concentration
                ele_volume = element_volume(xfield, ele)
                if (release(1) > 0.0) then
                   call set(chemfield, element_nodes(n), release(1) / ele_volume)
                else
                   call set(chemfield, element_nodes(n), tiny(1.0))
                end if 
             end if
          end do
       end do
    end do

    call profiler_toc("/chemical_exchange")

  end subroutine chemical_release

  subroutine particle_management(state, agent_list)
    type(state_type), intent(inout) :: state
    type(detector_linked_list), intent(inout) :: agent_list

    type(detector_type), pointer :: agent, agent_to_move, merge_target, agent_to_merge
    type(scalar_field), pointer :: agent_density_field
    type(vector_field), pointer :: xfield
    character(len=PYTHON_FUNC_LEN) :: pm_pycode
    type(ilist) :: elements_split_list, elements_merge_list
    type(integer_hash_table) :: split_index_mapping, merge_index_mapping, element_minima, element_maxima
    type(detector_linked_list), dimension(:), allocatable :: split_lists, merge_lists
    integer, dimension(:), allocatable :: elements_split, elements_merge, xbiggest_ids, xsmallest_ids
    real, dimension(:), allocatable :: xbiggest_sizes, xsmallest_sizes
    integer :: i, ele, index, splits_wanted, merges_wanted, found_agent, element_minimum, element_maximum
    real :: agent_density, ele_volume, minimum_density, maximum_density

    call profiler_tic(trim(agent_list%name)//"::particle_management")

    ewrite(2,*) "Performing Particle Management for list: ", trim(agent_list%name)

    xfield=>extract_vector_field(state, "Coordinate")
    agent_density_field=>extract_scalar_field(state, trim(agent_list%fgroup%name)//"Agents"//trim(agent_list%stage_name))

    ! Initialise the Python functions to determine PM limits
    if (have_option(trim(agent_list%stage_options)//"/particle_management/minimum/python")) then
       call get_option(trim(agent_list%stage_options)//"/particle_management/minimum/python", pm_pycode)
       call python_run_detector_string(trim(pm_pycode), trim(agent_list%name), trim("pm_minimum"))
    end if
    if (have_option(trim(agent_list%stage_options)//"/particle_management/maximum/python")) then
       call get_option(trim(agent_list%stage_options)//"/particle_management/maximum/python", pm_pycode)
       call python_run_detector_string(trim(pm_pycode), trim(agent_list%name), trim("pm_maximum"))
    end if

    ! First establish in which elements we need to split/merge
    ! We skip elements with no agents...
    ewrite(2,*) "Particle Management: Establishing elements to split/merge in"
    call allocate(element_minima)
    call allocate(element_maxima)
    do ele=1, ele_count(agent_density_field)
       agent_density = node_val(agent_density_field, ele)
       ele_volume = element_volume(xfield, ele)

       ! Here we either use the global minima/maxima from the options,
       ! or run Python code to determine the minimum/maximum for each element
       if (have_option(trim(agent_list%stage_options)//"/particle_management/minimum/python")) then
          call python_get_element_limit(ele, xfield, trim(agent_list%name), trim("pm_minimum"), minimum_density)
       else
          call get_option(trim(agent_list%stage_options)//"/particle_management/minimum/constant", minimum_density)
       end if

       if (have_option(trim(agent_list%stage_options)//"/particle_management/maximum/python")) then
          call python_get_element_limit(ele, xfield, trim(agent_list%name), trim("pm_maximum"), maximum_density)
       else
          call get_option(trim(agent_list%stage_options)//"/particle_management/maximum/constant", maximum_density)
       end if

       ! Store element numbers and minima/maxima for split and merge
       if (agent_density < minimum_density .and. agent_density > 0) then
          call insert(elements_split_list, ele)
          call insert(element_minima, ele, nint(minimum_density * ele_volume))
       end if

       if (agent_density > maximum_density .and. agent_density > 0) then
          call insert(elements_merge_list, ele)
          call insert(element_maxima, ele, nint(maximum_density * ele_volume))
       end if
    end do

    ! Convert linked lists of elements to split/merge to vectors
    allocate(elements_split(elements_split_list%length))
    elements_split = list2vector(elements_split_list)
    allocate(elements_merge(elements_merge_list%length))
    elements_merge = list2vector(elements_merge_list)

    ! Create a hashtable from element -> index into split_lists
    call allocate(split_index_mapping)
    do i=1, size(elements_split)
       call insert(split_index_mapping, elements_split(i), i)
    end do

    ! Create a hashtable from element -> index into merge_lists
    call allocate(merge_index_mapping)
    do i=1, size(elements_merge)
       call insert(merge_index_mapping, elements_merge(i), i)
    end do

    ! Put agents into temporary agent lists for each element to split/merge
    ewrite(2,*) "Particle Management: Moving agents to split/merge into temporary lists"
    allocate(split_lists(elements_split_list%length))
    allocate(merge_lists(elements_merge_list%length))
    agent=>agent_list%first
    do while(associated(agent))
       ! 
       if (has_key(split_index_mapping, agent%element)) then
          index = fetch(split_index_mapping,agent%element)
          agent_to_move=>agent
          agent=>agent%next
          call move(agent_to_move, agent_list, split_lists(index))
       ! 
       elseif (has_key(merge_index_mapping, agent%element)) then
          index = fetch(merge_index_mapping,agent%element)
          agent_to_move=>agent
          agent=>agent%next
          call move(agent_to_move, agent_list, merge_lists(index))

       else
          agent=>agent%next
       end if
    end do 

    ! Perform splits over the temporary agent arrays
    ! Each array now represents all agents in one element
    ! the according element number is stored in elements_split(i)
    ewrite(2,*) "Particle Management: Splitting agents per element"
    do i=1, size(split_lists)
       element_minimum = fetch(element_minima, elements_split(i))

       !!!!!!!!!!!!!!!!!!!!!
       ! Original VEW split algorithm:
       do while(split_lists(i)%length < element_minimum)
          ! 1) Find x, the number of splits we want
          splits_wanted = element_minimum - split_lists(i)%length
          splits_wanted = min(splits_wanted, split_lists(i)%length)

          ! 2) Get the x biggest agents
          allocate(xbiggest_ids(splits_wanted))
          allocate(xbiggest_sizes(splits_wanted))
          xbiggest_sizes = 0.0

          agent=>split_lists(i)%first
          do while(associated(agent))
             if (minval(xbiggest_sizes) < agent%biology(BIOVAR_SIZE)) then
                index = minval(minloc(xbiggest_sizes))
                xbiggest_ids(index) = agent%id_number 
                xbiggest_sizes(index) = agent%biology(BIOVAR_SIZE)
             end if

             agent=>agent%next
          end do

          ! 3) Split each of the x agents...
          agent=>split_lists(i)%first
          do while(associated(agent))
             if (find(xbiggest_ids, agent%id_number) > 0) then
                call pm_split(agent, split_lists(i))
             end if
             agent=>agent%next
          end do

          ! ...and recurse
          deallocate(xbiggest_ids)
          deallocate(xbiggest_sizes)
       end do
       !!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!
       ! Simple algorithm: always split first agent in the list
       !do while(agent_density < element_minimum)
       !   call pm_split(split_lists(i)%first, split_lists(i))
       !   agent_density = split_lists(i)%length / ele_volume
       !end do
       !!!!!!!!!!!!!!!!!!!!!

    end do

    ! Copy all agents back to the original list
    do i=1, size(split_lists)
       call move_all(split_lists(i), agent_list)
    end do

    ! Deallocate everythin split related
    deallocate(split_lists)
    deallocate(elements_split)
    call flush_list(elements_split_list)
    call deallocate(split_index_mapping)
    call deallocate(element_minima)

    ! Perform merges over the temporary agent arrays
    ewrite(2,*) "Particle Management: Merging agents per element"
    do i=1, size(merge_lists)
       element_maximum = fetch(element_maxima, elements_merge(i))

       !!!!!!!!!!!!!!!!!!!!!
       ! Original VEW merge algorithm:
       do while(merge_lists(i)%length > element_maximum)
          ! 1) Find x, the number of merges we want
          merges_wanted = merge_lists(i)%length - element_maximum
          ! We need to make sure we don't have more merge requests than agent pairs
          merges_wanted = floor(min(real(merges_wanted), real(merge_lists(i)%length) / 2.0))

          ! 2) Get the 2*x smallest agents (x pairs)
          allocate(xsmallest_ids(2*merges_wanted))
          allocate(xsmallest_sizes(2*merges_wanted))
          xsmallest_sizes = huge(1.0)

          agent=>merge_lists(i)%first
          do while(associated(agent))
             if (maxval(xsmallest_sizes) > agent%biology(BIOVAR_SIZE)) then
                index = minval(maxloc(xsmallest_sizes))
                xsmallest_ids(index) = agent%id_number 
                xsmallest_sizes(index) = agent%biology(BIOVAR_SIZE)
             end if

             agent=>agent%next
          end do

          ! 3) Merge pairwise...
          ! Note: We don't sort, but select merge pairs as we find them
          merge_target=>null()
          agent=>merge_lists(i)%first
          do while(associated(agent))
             found_agent = find(xsmallest_ids, agent%id_number)
             ! First hit is the target...
             if (found_agent > 0 .and. .not.associated(merge_target)) then
                merge_target=>agent
                agent=>agent%next
             ! ... second hit; we merge and reset the target
             elseif (found_agent > 0 .and. associated(merge_target)) then
                agent_to_merge=>agent
                agent=>agent%next
                call pm_merge(xfield, merge_target, agent_to_merge, merge_lists(i))
                merge_target=>null()
             else
                agent=>agent%next
             end if
          end do

          ! ...and recurse
          deallocate(xsmallest_ids)
          deallocate(xsmallest_sizes)
       end do
       !!!!!!!!!!!!!!!!!!!!!

       !!!!!!!!!!!!!!!!!!!!!
       ! Simple algorithm: always merge the first two agents in the list
       !do while(agent_density > agent_list%pm_max)
       !   call pm_merge(merge_lists(i)%first, merge_lists(i)%first%next, merge_lists(i))
       !   agent_density = merge_lists(i)%length / ele_volume
       !end do
       !!!!!!!!!!!!!!!!!!!!!

    end do

    ! Copy all agents back to the original list
    do i=1, size(merge_lists)
       call move_all(merge_lists(i), agent_list)
    end do

    ! Deallocate eveything merge related
    deallocate(merge_lists)
    deallocate(elements_merge)
    call flush_list(elements_merge_list)
    call deallocate(merge_index_mapping)
    call deallocate(element_maxima)

    call profiler_toc(trim(agent_list%name)//"::particle_management")

  contains

    function find(array, val) result(loc)
      !!< Find the first instance of val in array.
      integer, intent(in), dimension(:) :: array
      integer, intent(in) :: val
      integer :: i, loc

      loc = -1
      do i=1,size(array)
        if (array(i) == val) then
          loc = i
          return
        end if
      end do
    end function find

  end subroutine particle_management

  subroutine pm_split(agent, agent_list)
    type(detector_type), pointer, intent(inout) :: agent
    type(detector_linked_list), intent(inout) :: agent_list

    type(detector_type), pointer :: new_agent
    integer :: i

    ! Allocate and insert agent
    new_agent=>null()
    call allocate(new_agent, agent)
    call insert(new_agent,agent_list)

    ! Populate new agent
    call get_next_detector_id(new_agent%id_number)
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

    ! Distribute biomass
    agent%biology(BIOVAR_SIZE)=agent%biology(BIOVAR_SIZE) / 2.0
    new_agent%biology(BIOVAR_SIZE)=new_agent%biology(BIOVAR_SIZE) / 2.0

    ! Update the list
    agent_list%total_num_det=agent_list%total_num_det + 1

  end subroutine pm_split

  subroutine pm_merge(xfield, agent1, agent2, agent_list)
    type(vector_field), pointer, intent(in) :: xfield
    type(detector_type), pointer, intent(inout) :: agent1, agent2
    type(detector_linked_list), intent(inout) :: agent_list

    integer :: i

    ! Both agents are in the same element, 
    ! so we weight-average each physical coordinate
    ! and re-set the local coordinates
    do i=1, size(agent1%position)
       agent1%position(i)=wtavg(agent1%position(i),agent1%biology(BIOVAR_SIZE),agent2%position(i),agent2%biology(BIOVAR_SIZE))
    end do
    agent1%local_coords=local_coords(xfield,agent1%element,agent1%position)

    ! Weight-average the biology variables
    if (size(agent1%biology)>2) then
       do i=3, size(agent1%biology)
          agent1%biology(i) = wtavg(agent1%biology(i),agent1%biology(BIOVAR_SIZE),agent2%biology(i),agent2%biology(BIOVAR_SIZE))
       end do
    end if

    ! Add the sizes
    agent1%biology(BIOVAR_SIZE) = agent1%biology(BIOVAR_SIZE) + agent2%biology(BIOVAR_SIZE)

    call delete(agent2, agent_list)
    agent_list%total_num_det=agent_list%total_num_det - 1

  contains

    function wtavg(v1, w1, v2, w2) result(val)
      real, intent(in) :: v1, w1, v2, w2
      real :: val

      val = ((v1*w1)+(v2*w2))/(w1+w2)
    end function wtavg

  end subroutine pm_merge

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
