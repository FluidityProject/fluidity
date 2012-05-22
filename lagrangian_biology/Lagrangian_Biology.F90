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
  use lagrangian_biology_pm
  use lebiology_python

implicit none

  private

  public :: initialise_lagrangian_biology_metamodel, initialise_lagrangian_biology_agents, &
            lagrangian_biology_cleanup, update_lagrangian_biology, stage_aggregate, &
            get_num_functional_groups, get_functional_group
  public :: BIOFIELD_NONE, BIOFIELD_DIAG, BIOFIELD_UPTAKE, BIOFIELD_RELEASE, BIOFIELD_INGESTED, &
            BIOFIELD_FOOD_REQUEST, BIOFIELD_FOOD_INGEST

  type(functional_group), dimension(:), allocatable, target, save :: functional_groups

  character(len=FIELD_NAME_LEN), dimension(:), pointer :: uptake_field_names, release_field_names

  integer, parameter :: BIOFIELD_NONE=0, BIOFIELD_DIAG=1, BIOFIELD_UPTAKE=2, BIOFIELD_RELEASE=3, BIOFIELD_INGESTED=4, &
                        BIOFIELD_FOOD_REQUEST=5, BIOFIELD_FOOD_INGEST=6
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
    type(functional_group), pointer :: fgroup
    type(detector_linked_list), pointer :: agent_array
    integer, dimension(:), allocatable :: rnd_seed
    integer :: i, rnd_seed_int, rnd_dim
    integer :: j, fg, stage, n_fgroups, n_stages, n_env_fields, array

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

    n_fgroups = option_count("/embedded_models/lagrangian_ensemble_biology/functional_group")
    allocate(functional_groups(n_fgroups))
    ewrite(2,*) "Lagrangian biology: Allocated ", n_fgroups, " functional groups"

    ! Create a persistent dictionary of FG variable and environment name mappings
    call python_run_string("persistent['fg_var_names'] = dict()")
    call python_run_string("persistent['fg_env_names'] = dict()")

    ! Initialise Python module 'lebiology'
    call lebiology_init_module()

    ! We create an agent array for each stage of each Functional Group
    array=1
    do fg=1, n_fgroups
       fgroup => functional_groups(fg)
       write(fg_buffer, "(a,i0,a)") "/embedded_models/lagrangian_ensemble_biology/functional_group[",fg-1,"]"

       ! Get biology meta-data
       if (have_option(trim(fg_buffer)//"/variables")) then
          call read_functional_group(fgroup, trim(fg_buffer))

          ! Add variable names to the Python module
          call lebiology_add_variables(fgroup)
          call lebiology_add_foods(fgroup)
       else
          FLExit("No variables defined for functional group under: "//trim(fg_buffer))
       end if

       ! Record the environment fields to sample before agent update
       if (have_option(trim(fg_buffer)//"/environment")) then
          n_env_fields = option_count(trim(fg_buffer)//"/environment/field")
          allocate(fgroup%envfield_names(n_env_fields))
          allocate(fgroup%envfield_integrate(n_env_fields))
          do j=1, n_env_fields
             write(env_field_buffer, "(a,i0,a)") trim(fg_buffer)//"/environment/field[",j-1,"]"
             call get_option(trim(env_field_buffer)//"/name", fgroup%envfield_names(j))

             if (have_option(trim(env_field_buffer)//"/integrate_along_path")) then
                fgroup%envfield_integrate(j) = .true.
             else
                fgroup%envfield_integrate(j) = .false.
             end if
          end do

          ! Add field names to the Python module
          call lebiology_add_envfields(fgroup)
       else
          FLExit("No environment fields defined for FG::"//trim(fgroup%name))
       end if

       ! Get initialisation options
       if (have_option(trim(fg_buffer)//"/initial_state")) then
          fgroup%init_options = trim(fg_buffer)//"/initial_state"
       else
          FLExit("No initialisation options found for FG::"//trim(fgroup%name))
       end if

       n_stages = option_count(trim(fg_buffer)//"/stages/stage")
       allocate(fgroup%agent_arrays(n_stages))
       ewrite(2,*) "Lagrangian biology: Allocated ", n_stages, " agent arrays for FG::", trim(fgroup%name)
       do stage = 1, n_stages
          agent_array => fgroup%agent_arrays(stage)
          write(stage_buffer, "(a,i0,a)") trim(fg_buffer)//"/stages/stage[",stage-1,"]"
          agent_array%stage_options = trim(stage_buffer)
          call get_option(trim(stage_buffer)//"/name", stage_name)
          agent_array%fgroup => fgroup

          agent_array%name = trim(functional_groups(fg)%name)//trim(stage_name)
          agent_array%stage_name = trim(stage_name)
          agent_array%id=array
          call lebiology_set_stage_id(fgroup, stage_name, float(array))

          ! Register the agent array, so Zoltan/Adaptivity will not forget about it
          call register_detector_list(agent_array)

          ! Get options for lagrangian detector movement
          call read_detector_move_options(agent_array, trim(fg_buffer)//"/movement")

          ! Get options for Random Walk
          call read_random_walk_options(agent_array, trim(fg_buffer)//"/movement")

          ! Now we know the variable names for the FG, so we record them in Python
          call python_run_string("persistent['fg_var_names']['"//trim(agent_array%name)//"'] = []")
          do j=1, size(functional_groups(fg)%variables)
             call python_run_string("persistent['fg_var_names']['"//trim(agent_array%name)//"'].append('"//trim(functional_groups(fg)%variables(j)%name)//"')")
          end do

          ! Store the python update code
          if (have_option(trim(stage_buffer)//"/biology/python")) then
             call get_option(trim(stage_buffer)//"/biology/python", agent_array%biovar_pycode)
          end if

          if (have_option(trim(stage_buffer)//"/io_period")) then
             call get_option(trim(stage_buffer)//"/io_period", agent_array%output_period)
          else
             agent_array%output_period = 1.0
          end if

          array = array + 1
       end do 
    end do

    ewrite(1,*) "Lagrangian biology: Initialised metamodel"

  end subroutine initialise_lagrangian_biology_metamodel

  subroutine initialise_lagrangian_biology_agents(state)
    type(state_type), dimension(:), intent(inout) :: state

    type(functional_group), pointer :: fgroup
    type(detector_linked_list), pointer :: agent_array
    type(detector_linked_list) :: init_array
    type(detector_type), pointer :: agent
    type(vector_field), pointer :: xfield
    character(len=PYTHON_FUNC_LEN) :: pos_func, bio_func
    character(len=OPTION_PATH_LEN) :: stage_buffer
    real, allocatable, dimension(:,:) :: coords
    real:: current_time
    integer :: ag, fg, stage, dim, n_agents

    if (.not.have_option("/embedded_models/lagrangian_ensemble_biology")) return

    ewrite(1,*) "Lagrangian biology: Initialising agents..."

    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)
    xfield=>extract_vector_field(state(1), "Coordinate")

    ! Initialise the shared ID counter to generate agent IDs
    call init_id_counter()
    ewrite(2,*) "Lagrangian biology: Initialised ID counter"

    ! We create an agent array for each stage of each Functional Group
    do fg=1, get_num_functional_groups()
       fgroup => get_functional_group(fg)

       ! Create agents and insert into list
       call get_option(trim(fgroup%init_options)//"/number_of_agents", n_agents)
       call get_option(trim(fgroup%init_options)//"/position", pos_func)
       allocate(coords(dim,n_agents))
       call set_detector_coords_from_python(coords, n_agents, pos_func, current_time)

       do ag=1, n_agents
          call create_single_detector(init_array, xfield, coords(:,ag), LAGRANGIAN_DETECTOR)
       end do
       deallocate(coords)

       ! Initialise agent biology variables
       if (allocated(fgroup%variables)) then
          if (have_option(trim(fgroup%init_options)//"/biology")) then
             call get_option(trim(fgroup%init_options)//"/biology", bio_func)
             call lebiology_prepare_pyfunc(fgroup, "Agent_Init", bio_func)

             agent => init_array%first
             do while (associated(agent))
                allocate(agent%biology(size(fgroup%variables)))
                call lebiology_initialise_agent(fgroup, "Agent_Init", agent)

                ! Allocate request/ingest variety vectors
                if ( size(fgroup%food_sets) > 0) then
                   allocate(agent%food_requests( size(fgroup%food_sets(1)%varieties) ))
                   agent%food_requests = 0.0
                   allocate(agent%food_ingests( size(fgroup%food_sets(1)%varieties) ))
                   agent%food_ingests = 0.0
                end if
                agent => agent%next
             end do
          end if
       end if

       ewrite(2,*) "Lagrangian biology: Initialised ", init_array%length, " agents for FG::", trim(fgroup%name)

       ! After initialising agents in a temporary array 
       ! we now move them according to their stage
       call distribute_by_stage(fgroup, init_array)

       do stage=1, size(fgroup%agent_arrays)
          agent_array => fgroup%agent_arrays(stage)
          stage_buffer = trim(agent_array%stage_options)

          ! Only allow binary output
          agent_array%binary_output=.true.

          ! Set agent%list_id, because we need it to detect stage changes
          agent=>agent_array%first
          do while (associated(agent))
             agent%list_id=agent_array%id
             agent=>agent%next
          end do

          call write_detector_header(state, agent_array)

          ! Derive diagnostics from initial state...
          call derive_primary_diagnostics(state(1), agent_array)

       end do  ! Stage

       ! ... and aggregate them
       call aggregate_diagnostics_by_stage(state(1), fgroup)

    end do  ! FGroup

    ewrite(1,*) "Lagrangian biology: Initialised agents"

  end subroutine initialise_lagrangian_biology_agents

  subroutine read_functional_group(fgroup, fg_path) 
    type(functional_group), pointer, intent(inout) :: fgroup
    character(len=*), intent(in) :: fg_path

    character(len=OPTION_PATH_LEN) :: var_buffer, stage_buffer, food_buffer, food_type_buffer
    character(len=FIELD_NAME_LEN) :: biovar_name, field_name, food_name, vname
    type(ilist) :: motion_variables, history_variables, food_chem_inds, food_target_inds, uptake_vars, release_vars
    type(functional_group), pointer :: target_fg
    type(food_variety), pointer :: fvariety
    integer :: i, f, t, v, vars_state, vars_hist, vars_chem, vars_ingest, vars_uptake, vars_release, &
               vars_total, var_index, chemvar_index, stage_count, &
               n_food_sets, n_food_types, n_fvarieties, hist, history_depth

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
    vars_hist = 0
    do i=1, vars_state
       write(var_buffer, "(a,i0,a)") trim(fg_path)//"/variables/state_variable[",i-1,"]"
       if (have_option(trim(var_buffer)//"/history")) then
          call get_option(trim(var_buffer)//"/history", history_depth)
          vars_hist = vars_hist + history_depth - 1
       end if
    end do

    vars_ingest = 0
    vars_uptake = 0
    vars_release = 0
    do i=1, vars_chem
       write(var_buffer, "(a,i0,a)") trim(fg_path)//"/variables/chemical_variable[",i-1,"]"
       if (have_option(trim(var_buffer)//"/scalar_field::Ingested")) then
          vars_ingest = vars_ingest + 1
       end if
       if (have_option(trim(var_buffer)//"/uptake")) then
          vars_uptake = vars_uptake + 1
       end if
       if (have_option(trim(var_buffer)//"/release")) then
          vars_release = vars_release + 1
       end if
    end do
    vars_total = vars_state + vars_hist + vars_chem + vars_ingest + vars_uptake + vars_release

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

          ! Record indices of motion variables
          if (have_option(trim(var_buffer)//"/include_in_motion")) then
             call insert(motion_variables, var_index)
          end if

          ! Record according diagnostic field
          if (have_option(trim(var_buffer)//"/scalar_field")) then
             fgroup%variables(var_index)%field_type = BIOFIELD_DIAG
             call get_option(trim(var_buffer)//"/scalar_field/name", field_name)                   
             fgroup%variables(var_index)%field_name = trim(fgroup%name)//trim(field_name)//trim(biovar_name)
             fgroup%variables(var_index)%field_path = trim(var_buffer)//"/scalar_field"
          else
             fgroup%variables(var_index)%field_type = BIOFIELD_NONE
          end if
          var_index = var_index+1

          if (have_option(trim(var_buffer)//"/history")) then
             call get_option(trim(var_buffer)//"/history", history_depth)
             do hist = 2, history_depth
                fgroup%variables(var_index)%name=trim(biovar_name)//"_"//int2str(hist - 1)
                fgroup%variables(var_index)%field_type = BIOFIELD_NONE
                fgroup%variables(var_index)%pool_index = var_index - 1
                call insert(history_variables, var_index)
                var_index = var_index+1
             end do
          end if
       end do

       allocate(fgroup%history_var_inds(history_variables%length))
       fgroup%history_var_inds = list2vector(history_variables)

       allocate(fgroup%motion_var_inds(motion_variables%length))
       fgroup%motion_var_inds = list2vector(motion_variables)

       ! Add chemical variables
       do i=1, vars_chem
          write(var_buffer, "(a,i0,a)") trim(fg_path)//"/variables/chemical_variable[",i-1,"]"
          call get_option(trim(var_buffer)//"/name", biovar_name)
          if (have_option(trim(var_buffer)//"/include_in_io")) then
             fgroup%variables(var_index)%write_to_file = .true.
          end if

          ! Chemical pool variable and according diagnostic fields
          fgroup%variables(var_index)%name = trim(biovar_name)
          if (have_option(trim(var_buffer)//"/scalar_field::Particulate")) then
             fgroup%variables(var_index)%field_type = BIOFIELD_DIAG
             call get_option(trim(var_buffer)//"/scalar_field/name", field_name)
             fgroup%variables(var_index)%field_name = trim(fgroup%name)//trim(field_name)//trim(biovar_name)
             fgroup%variables(var_index)%field_path = trim(var_buffer)//"/scalar_field"
          else
             fgroup%variables(var_index)%field_type = BIOFIELD_NONE
          end if
          chemvar_index = var_index
          var_index = var_index+1

          ! Create a 'ChemIngested' variable
          if (have_option(trim(var_buffer)//"/scalar_field::Ingested")) then
             fgroup%variables(var_index)%name = trim(biovar_name)//"Ingested"
             fgroup%variables(var_index)%field_type = BIOFIELD_INGESTED
             fgroup%variables(var_index)%field_name = trim(fgroup%name)//"Ingested"//trim(biovar_name)
             fgroup%variables(var_index)%field_path = trim(var_buffer)//"/scalar_field::Ingested"

             fgroup%variables(var_index)%pool_index = chemvar_index
             !fgroup%variables(var_index)%request_index = var_index - 1
             fgroup%variables(chemvar_index)%ingest_index = var_index

             if (have_option(trim(var_buffer)//"/scalar_field::Ingested/include_in_io")) then
                fgroup%variables(var_index)%write_to_file = .true.
             end if

             var_index = var_index+1
          end if

          ! Chemical uptake
          ! Note check for existence of Request and Depletion field, and record Depletion field
          if (have_option(trim(var_buffer)//"/uptake")) then
             fgroup%variables(var_index)%name = trim(biovar_name)//"Uptake"
             fgroup%variables(var_index)%field_type = BIOFIELD_UPTAKE
             fgroup%variables(var_index)%field_name = trim(fgroup%name)//"Request"//trim(biovar_name)
             fgroup%variables(var_index)%field_path = trim(var_buffer)//"/uptake/scalar_field::Request"
             fgroup%variables(var_index)%depletion_field_path = trim(var_buffer)//"/uptake/scalar_field::Depletion"

             call get_option(trim(var_buffer)//"/uptake/source_field/name", fgroup%variables(var_index)%chemfield)
             fgroup%variables(var_index)%i_chemfield = insert_global_uptake_field(fgroup%variables(var_index)%chemfield)

             if (have_option(trim(var_buffer)//"/uptake/integrate_along_path")) then
                fgroup%variables(var_index)%path_integration = .true.
             end if

             if (have_option(trim(var_buffer)//"/uptake/include_in_io")) then
                fgroup%variables(var_index)%write_to_file = .true.
             end if

             fgroup%variables(var_index)%pool_index = chemvar_index
             call insert(uptake_vars, var_index)
             var_index = var_index+1
          end if

          ! Chemical release
          if (have_option(trim(var_buffer)//"/release")) then
             fgroup%variables(var_index)%name = trim(biovar_name)//"Release"
             fgroup%variables(var_index)%field_type = BIOFIELD_RELEASE
             fgroup%variables(var_index)%field_name = trim(fgroup%name)//"Release"//trim(biovar_name)
             fgroup%variables(var_index)%field_path = trim(var_buffer)//"/release/scalar_field::Release"
             call get_option(trim(var_buffer)//"/release/target_field/name", fgroup%variables(var_index)%chemfield)
             fgroup%variables(var_index)%i_chemfield = insert_global_release_field(fgroup%variables(var_index)%chemfield)

             if (have_option(trim(var_buffer)//"/release/integrate_along_path")) then
                fgroup%variables(var_index)%path_integration = .true.
             end if

             if (have_option(trim(var_buffer)//"/release/include_in_io")) then
                fgroup%variables(var_index)%write_to_file = .true.
             end if

             fgroup%variables(var_index)%pool_index = chemvar_index
             call insert(release_vars, var_index)
             var_index = var_index+1
          end if
       end do

       allocate(fgroup%ivars_uptake(uptake_vars%length))
       fgroup%ivars_uptake = list2vector(uptake_vars)

       allocate(fgroup%ivars_release(release_vars%length))
       fgroup%ivars_release = list2vector(release_vars)

    end if ! vars > 0

    n_food_sets = option_count(trim(fg_path)//"/food_set")
    n_fvarieties = 0
    do i=1, n_food_sets
       write(food_buffer, "(a,i0,a)") trim(fg_path)//"/food_set[",i-1,"]"
       n_food_types = option_count(trim(food_buffer)//"/food_type")
       n_fvarieties = n_fvarieties + n_food_types
    end do

    allocate(fgroup%food_sets(n_food_sets))
    do i=1, n_food_sets
       write(food_buffer, "(a,i0,a)") trim(fg_path)//"/food_set[",i-1,"]"

       call get_option(trim(food_buffer)//"/name", food_name)
       fgroup%food_sets(i)%name = trim(food_name)
       call get_option(trim(food_buffer)//"/functional_group", fgroup%food_sets(i)%target_fgroup)

       n_food_types = option_count(trim(food_buffer)//"/food_type")
       allocate(fgroup%food_sets(i)%varieties(n_food_types))
       do f=1, n_food_types
          fvariety => fgroup%food_sets(i)%varieties(f)
          write(food_type_buffer, "(a,i0,a)") trim(food_buffer)//"/food_type[",f-1,"]"
          call get_option(trim(food_type_buffer)//"/name", vname)
          fvariety%name = trim(vname)

          fvariety%vrequest%name = trim(food_name)//"Request"//trim(vname)
          fvariety%vrequest%field_type = BIOFIELD_FOOD_REQUEST
          fvariety%vrequest%field_name = trim(fgroup%name)//trim(food_name)//"Request"//trim(vname)
          fvariety%vrequest%field_path = trim(food_buffer)//"/scalar_field::Request"
          fvariety%vrequest%depletion_field_path = trim(food_buffer)//"/scalar_field::Depletion"
          if (have_option(trim(food_buffer)//"/scalar_field::Request/stage_diagnostic")) then
             fvariety%vrequest%stage_diagnostic = .true.
          end if
          if (have_option(trim(food_buffer)//"/scalar_field::Request/include_in_io")) then
             fvariety%vrequest%write_to_file = .true.
          end if

          fvariety%vingest%name = trim(food_name)//"IngestedCells"//trim(vname)
          fvariety%vingest%field_type = BIOFIELD_FOOD_INGEST
          fvariety%vingest%field_name = trim(fgroup%name)//trim(food_name)//"IngestedCells"//trim(vname)
          fvariety%vingest%field_path = trim(food_buffer)//"/scalar_field::IngestedCells"
          if (have_option(trim(food_buffer)//"/scalar_field::IngestedCells/stage_diagnostic")) then
             fvariety%vingest%stage_diagnostic = .true.
          end if
          if (have_option(trim(food_buffer)//"/scalar_field::IngestedCells/include_in_io")) then
             fvariety%vingest%write_to_file = .true.
          end if
       end do

       ! Now we collect the indices of the chemical pools to ingest
       ! This should probably happen after all FGs are created, but it works for now...

       ! First find our target FG
       target_fg => null()
       do v=1, get_num_functional_groups()
          if (trim(functional_groups(v)%name) == trim(fgroup%food_sets(i)%target_fgroup)) then
             target_fg => functional_groups(v)
          end if
       end do
       if (.not. associated(target_fg)) then
          ewrite(-1,*) "Target for food set ", trim(fgroup%name), "::", trim(fgroup%food_sets(i)%name), " not found!"
          FLExit("Unknown functional group "//trim(fgroup%food_sets(i)%target_fgroup))
       end if

       do v=1, var_index
          if (fgroup%variables(v)%field_type==BIOFIELD_INGESTED) then
             ! We have an Ingested field, check if target FG has the according Pool
             vname = fgroup%variables( fgroup%variables(v)%pool_index )%name
             do t=1, size(target_fg%variables)
                if (trim(vname)==trim(target_fg%variables(t)%name)) then
                   call insert(food_chem_inds, fgroup%variables(v)%pool_index)
                end if
             end do
          end if
       end do

       ! Now store our index pairs
       allocate(fgroup%food_sets(i)%ingest_chem_inds(food_chem_inds%length))
       fgroup%food_sets(i)%ingest_chem_inds = list2vector(food_chem_inds)

       ! And finally we go through all target stages and 
       ! associate a pointer with the according agent list
       do t=1, size(fgroup%food_sets(i)%varieties)
          do v=1, size(target_fg%agent_arrays)
             if (trim(target_fg%agent_arrays(v)%stage_name) == trim(fgroup%food_sets(i)%varieties(t)%name) ) then
                 fgroup%food_sets(i)%varieties(t)%target_list%ptr => target_fg%agent_arrays(v)
                 fgroup%food_sets(i)%varieties(t)%conc_field = trim(target_fg%variables(BIOVAR_SIZE)%field_name)//trim(fgroup%food_sets(i)%varieties(t)%name)
             end if
          end do
       end do
    end do

  end subroutine read_functional_group

  subroutine lagrangian_biology_cleanup(state)
    type(state_type), intent(inout) :: state

    type(functional_group), pointer :: fgroup
    type(detector_linked_list), pointer :: agent_array
    type(mesh_type) :: biology_mesh
    integer :: fg, stage, ierror

    do fg=1, get_num_functional_groups()
       fgroup => get_functional_group(fg)
       do stage=1, size(fgroup%agent_arrays)
          agent_array=>fgroup%agent_arrays(stage)

          call delete_all(agent_array)
          if (agent_array%output_unit/=0) then
             call MPI_FILE_CLOSE(agent_array%output_unit, ierror) 
             if(ierror /= MPI_SUCCESS) then
                ewrite(0,*) "Warning: failed to close .detector file open with mpi_file_open"
             end if
          end if
       end do
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

    type(functional_group), pointer :: fgroup
    type(detector_linked_list), pointer :: agent_array, new_agent_list
    type(detector_linked_list) :: stage_change_list
    type(detector_type), pointer :: agent, agent_to_move
    type(vector_field), pointer :: xfield
    type(scalar_field_pointer), dimension(:), pointer :: env_fields, food_fields
    character(len=FIELD_NAME_LEN) :: foodname
    type(food_set) :: fset
    integer :: i, j, f, v, env, pm_period, hvar, hvar_ind, hvar_src_ind, fg, stage

    ewrite(1,*) "Lagrangian biology: Updating agents..."
    call profiler_tic("/update_lagrangian_biology")

    xfield=>extract_vector_field(state(1), "Coordinate")

    ! Prepare python-state
    ! Note: Currently disabled for performance reasons. 
    ! Actually only needed for Python RW with field sampling

    !call profiler_tic("/update_lagrangian_biology::python_reload")
    !call python_reset()
    !call python_add_state(state(1))
    !call profiler_toc("/update_lagrangian_biology::python_reload")

    do fg=1, get_num_functional_groups()
       fgroup => get_functional_group(fg)

       ! Aggregate food concentrations
       if (size(fgroup%food_sets) > 0) then
          !call aggregate_food_diagnostics(state(1), fgroup)

          ! Note: Assume there is only one FoodSet for now...
          fset = fgroup%food_sets(1)
          foodname = trim(fset%name)
          allocate(food_fields(size(fset%varieties)))
          do f=1, size(fset%varieties)
             food_fields(f)%ptr => extract_scalar_field(state(1), trim(fset%varieties(f)%conc_field))
          end do
       else
          allocate(food_fields(0))
       end if

       if (allocated(fgroup%envfield_names)) then
          allocate(env_fields(size(fgroup%envfield_names)))
          do env=1, size(fgroup%envfield_names)
             env_fields(env)%ptr => extract_scalar_field(state(1), fgroup%envfield_names(env))

             ! If we want to integrate along the agent path we need to check that field is P0
             if (fgroup%envfield_integrate(env) .and. &
                 element_degree(env_fields(env)%ptr, 1) /= 0) then
                FLExit("Path integration only available for P0 environment fields")
             end if
          end do
       end if

       do stage=1, size(fgroup%agent_arrays)
          agent_array=>fgroup%agent_arrays(stage)

          ! Advance history variables
          if (have_option(trim(agent_array%stage_options)//"/biology")) then
             agent=>agent_array%first
             do while (associated(agent))
                ! Loop over variables backwards
                do hvar=size(fgroup%history_var_inds), 1, -1
                   hvar_ind = fgroup%history_var_inds(hvar)
                   hvar_src_ind = fgroup%variables(hvar_ind)%pool_index
                   agent%biology(hvar_ind) = agent%biology(hvar_src_ind)
                end do
                agent=>agent%next
             end do
          end if

          ! Move lagrangian detectors
          if (check_any_lagrangian(agent_array)) then
             call move_lagrangian_detectors(state, agent_array, dt, timestep)
          end if

          if (have_option(trim(agent_array%stage_options)//"/biology") .and. agent_array%length > 0 ) then
             ewrite(2,*) "Lagrangian biology: Updating ", trim(agent_array%name)

             if (have_option(trim(agent_array%stage_options)//"/biology/python")) then
                ! Compile python function to set bio-variable, and store in the global dictionary
                call lebiology_prepare_pyfunc(fgroup, trim(agent_array%stage_name)//"_Update", agent_array%biovar_pycode)
             end if

             ! Update agent biology
             agent=>agent_array%first
             do while (associated(agent))

                ! Reset Request variables
                do v=1, size(fgroup%variables)
                   if ( fgroup%variables(v)%field_type==BIOFIELD_UPTAKE .or. &
                        fgroup%variables(v)%field_type==BIOFIELD_RELEASE) then
                      agent%biology(v) = 0.0
                   end if
                end do

                if (size(fgroup%food_sets) > 0) then
                   do v=1, size(agent%food_requests)
                      agent%food_requests(v) = 0.0
                   end do
                end if

                ! Update agent via the Python module
                call lebiology_update_agent(fgroup, trim(agent_array%stage_name)//"_Update", &
                          trim(foodname), agent, xfield, env_fields, food_fields, dt)

                ! Reset Ingested variables
                do v=1, size(fgroup%variables)
                   if ( fgroup%variables(v)%field_type==BIOFIELD_INGESTED) then
                      agent%biology(v) = 0.0
                   end if
                end do

                if (size(fgroup%food_sets) > 0) then
                   do v=1, size(agent%food_ingests)
                      agent%food_ingests(v) = 0.0
                   end do
                end if

                ! Check for stage change
                if (nint(agent%biology(BIOVAR_STAGE)) /= agent_array%id) then
                   agent_to_move=>agent
                   agent=>agent%next
                   call move(agent_to_move, agent_array, stage_change_list)
                else
                   agent=>agent%next
                end if
             end do

          end if  ! have_biology

       end do  ! stages

       if (associated(env_fields)) then
          deallocate(env_fields)
       end if

       ! Deal with newly created agents
       new_agent_list => get_new_agent_list()
       agent => new_agent_list%first
       do while (associated(agent))
          ! Let the picker determine parametric coordinates
          ! Note: Picker will only be run locally, since it does MPI comms
          ! The assumption is that an agent created on this proc
          ! will exist on the local partition
          call picker_inquire(xfield,agent%position,agent%element,local_coord=agent%local_coords,global=.false.)
          if (agent%element <= 0) then
             ewrite(-1,*) "Agent added outside the computational domain!"
             ewrite(-1,*) "Position:", agent%position, " , element:", agent%element, ", local coordinates:", agent%local_coords
             FLAbort("Error establishing parametric coordinates for new agent")
          end if

          if (size(fgroup%food_sets) > 0) then
             if (.not.allocated(agent%food_requests)) then
                allocate(agent%food_requests(size(fgroup%food_sets(1)%varieties)))
             end if
             if (.not.allocated(agent%food_ingests)) then
                allocate(agent%food_ingests(size(fgroup%food_sets(1)%varieties)))
             end if
          end if

          agent=>agent%next
       end do
       call distribute_by_stage(fgroup, new_agent_list)

       ! Handle stage changes within FG
       call distribute_by_stage(fgroup, stage_change_list)

    end do  ! FGroup

    ! Execute chemical uptake/release
    call profiler_tic("/update_lagrangian_biology::chemical_exchange")
    call aggregate_chemical_diagnostics(state(1))
    call chemical_uptake(state(1))
    call chemical_release(state(1))
    call profiler_toc("/update_lagrangian_biology::chemical_exchange")

    call profiler_tic("/update_lagrangian_biology::ingestion")
    call ingestion(state(1))
    call profiler_toc("/update_lagrangian_biology::ingestion")

    ! Particle Management
    do fg=1, get_num_functional_groups()
       fgroup => get_functional_group(fg)
       do stage=1, size(fgroup%agent_arrays)
          agent_array => fgroup%agent_arrays(stage)

          ! First remove all agent with zero biomass to avoid div-by-zero errors
          call pm_strip_insignificant(agent_array)

          ! Re-derive agent density field
          call derive_primary_diagnostics(state(1), agent_array)

          if (have_option(trim(agent_array%stage_options)//"/particle_management") .and. &
              agent_array%length > 0) then
             call get_option(trim(agent_array%stage_options)//"/particle_management/period_in_timesteps", pm_period)
             if (timestep == 1 .or. mod(timestep, pm_period) == 0) then
                call particle_management(state(1), agent_array)
             end if
          end if

          ! Output agent positions after re-sampling
          if (timestep == 1 .or. mod(time, agent_array%output_period) == 0) then
             call write_detectors(state, agent_array, time, dt, timestep)
          end if

          ! Re-derive the required agent diagnostic variables...
          call derive_primary_diagnostics(state(1), agent_array)
       end do
       ! ... and aggregate them
       call aggregate_diagnostics_by_stage(state(1), fgroup)
    end do

    call profiler_toc("/update_lagrangian_biology")

  end subroutine update_lagrangian_biology

  function stage_aggregate(var)
    type(le_variable), intent(in) :: var
    logical :: stage_aggregate

    if (have_option(trim(var%field_path)//"/stage_aggregate")) then
       stage_aggregate = .true.
    else
       stage_aggregate = .false.
    end if
  end function

  subroutine distribute_by_stage(fgroup, agent_list)
    type(functional_group), pointer, intent(inout):: fgroup
    type(detector_linked_list), intent(inout):: agent_list

    type(detector_type), pointer :: agent, agent_to_move
    integer :: j

    if (agent_list%length <= 0) then
       return
    end if

    ! Handle stage changes within FG
    ewrite(2,*) "Lagrangian biology: Distributing ", agent_list%length," agents by stage for FG::", fgroup%name
    agent=>agent_list%first
    agent_loop: do while (associated(agent))
       do j=1, size(fgroup%agent_arrays)
          if (fgroup%agent_arrays(j)%id==nint(agent%biology(BIOVAR_STAGE))) then
             agent_to_move=>agent
             agent=>agent%next
             call move(agent_to_move, agent_list, fgroup%agent_arrays(j))
             cycle agent_loop
          end if
       end do
       ewrite(-1,*) "Lagrangian biology: Target stage ID ", agent%biology(BIOVAR_STAGE), "not defined"
       FLExit("Lagrangian biology: Target stage not found")
    end do agent_loop

  end subroutine distribute_by_stage

  subroutine derive_primary_diagnostics(state, agent_list)
    ! Set per-stage diagnostic fields from agent variables, 
    ! including agent counts and chemical request/release fields
    type(state_type), intent(inout) :: state
    type(detector_linked_list), intent(inout) :: agent_list

    type(vector_field), pointer :: xfield
    type(scalar_field_pointer), dimension(size(agent_list%fgroup%variables)) :: diagfields
    type(scalar_field), pointer :: agent_count_field
    type(detector_type), pointer :: agent
    type(food_variety), pointer :: fvariety
    integer :: i, ele
    real :: ele_volume, release_amount

    ewrite(2,*) "Lagrangian biology: Deriving primary diagnostic fields for ", (agent_list%name)

    call profiler_tic(trim(agent_list%name)//"::primary_diagnostics")

    xfield=>extract_vector_field(state, "Coordinate")

    ! Pull and reset agent density field
    agent_count_field=>extract_scalar_field(state, trim(agent_list%fgroup%name)//"Agents"//trim(agent_list%stage_name))
    call zero(agent_count_field)

    ! Pull and reset all per-stage-array diagnostic fields
    do i=1, size(agent_list%fgroup%variables)
       if (agent_list%fgroup%variables(i)%field_type == BIOFIELD_DIAG .or. &
                  agent_list%fgroup%variables(i)%field_type == BIOFIELD_INGESTED ) then

          diagfields(i)%ptr => extract_scalar_field(state, trim(agent_list%fgroup%variables(i)%field_name)//trim(agent_list%stage_name))
          call zero(diagfields(i)%ptr)
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
                call addto(diagfields(i)%ptr, agent%element, agent%biology(i) / ele_volume)
             else
                call addto(diagfields(i)%ptr, agent%element, agent%biology(i)*agent%biology(BIOVAR_SIZE) / ele_volume)
             end if

          ! Ingestion, uptake and release requests are total amounts, so don't divide by element volume
          elseif (agent_list%fgroup%variables(i)%field_type == BIOFIELD_INGESTED ) then

             if (agent_list%fgroup%variables(i)%path_integration) then
                call integrate_along_path(diagfields(i)%ptr, xfield, agent, agent%biology(i))
             else
                call addto(diagfields(i)%ptr, agent%element, agent%biology(i)*agent%biology(BIOVAR_SIZE) / ele_volume)
             end if
          end if
       end do

       agent => agent%next
    end do

    call profiler_toc(trim(agent_list%name)//"::primary_diagnostics")

  end subroutine derive_primary_diagnostics

  subroutine aggregate_diagnostics_by_stage(state, fgroup)
    type(state_type), intent(inout) :: state
    type(functional_group), pointer, intent(inout):: fgroup

    type(scalar_field), pointer :: diagfield_agg, diagfield_stage
    type(le_variable), pointer :: var
    type(food_variety), pointer :: fvariety
    integer :: v, s

    ewrite(2,*) "Lagrangian biology: Aggregating stage-diagnostics"

    ! Now we derive all aggregated quantities
    do v=1, size(fgroup%variables)
       var => fgroup%variables(v)

       ! Aggregate stage-aggregated diagnostic fields
       if (stage_aggregate(var) .and. var%field_type /= BIOFIELD_NONE) then
          diagfield_agg=>extract_scalar_field(state, trim(var%field_name))
          call zero(diagfield_agg)

          do s=1, size(fgroup%stage_names%ptr)
             diagfield_stage=>extract_scalar_field(state, trim(var%field_name)//trim(fgroup%stage_names%ptr(s)))
             call addto(diagfield_agg, diagfield_stage)
          end do
       end if
    end do
  end subroutine aggregate_diagnostics_by_stage

  subroutine aggregate_chemical_diagnostics(state)
    type(state_type), intent(inout) :: state

    type(vector_field), pointer :: xfield
    type(functional_group), pointer :: fgroup
    type(detector_type), pointer :: agent
    type(scalar_field_pointer), dimension(:), allocatable :: uptake_diagfields, release_diagfields
    type(scalar_field_pointer), dimension(size(uptake_field_names)) :: uptake_fields
    type(scalar_field_pointer), dimension(size(release_field_names)) :: release_fields
    type(le_variable), pointer :: var
    real :: ele_volume
    integer :: f, v, fg, stage, ivar

    ewrite(2,*) "Lagrangian biology: Aggregating chemical request/release fields"

    xfield=>extract_vector_field(state, "Coordinate")

    ! Pull and reset global fields
    do f=1, size(uptake_field_names)
       uptake_fields(f)%ptr => extract_scalar_field(state, trim(uptake_field_names(f))//"Request")
       call zero(uptake_fields(f)%ptr)
    end do

    do f=1, size(release_field_names)
       release_fields(f)%ptr => extract_scalar_field(state, trim(release_field_names(f))//"Release")
       call zero(release_fields(f)%ptr)
    end do

    do fg=1, get_num_functional_groups()
       fgroup => get_functional_group(fg)

       ! Pull and reset FG-level diagnostic fields
       allocate(uptake_diagfields(size(fgroup%ivars_uptake)))
       do v=1, size(fgroup%ivars_uptake)
          var => fgroup%variables( fgroup%ivars_uptake(v) )
          uptake_diagfields(v)%ptr => extract_scalar_field(state, var%field_name)
          call zero(uptake_diagfields(v)%ptr)
       end do

       allocate(release_diagfields(size(fgroup%ivars_release)))
       do v=1, size(fgroup%ivars_release)
          var => fgroup%variables( fgroup%ivars_release(v) )
          release_diagfields(v)%ptr => extract_scalar_field(state, var%field_name)
          call zero(release_diagfields(v)%ptr)
       end do

       ! Derive FG-level diagnostic fields
       do stage=1, size(fgroup%agent_arrays)
          agent => fgroup%agent_arrays(stage)%first
          do while (associated(agent))

             do v=1, size(fgroup%ivars_uptake)
                ivar = fgroup%ivars_uptake(v)
                if (fgroup%variables(ivar)%path_integration) then
                   call integrate_along_path(uptake_diagfields(v)%ptr, xfield, agent, agent%biology(ivar) )
                else
                   ele_volume = element_volume(xfield, agent%element)
                   call addto(uptake_diagfields(v)%ptr, agent%element, agent%biology(ivar)*agent%biology(BIOVAR_SIZE) / ele_volume)
                end if
             end do

             do v=1, size(fgroup%ivars_release)
                ivar = fgroup%ivars_release(v)
                if (fgroup%variables(ivar)%path_integration) then
                   call integrate_along_path(release_diagfields(v)%ptr, xfield, agent, agent%biology(ivar) )
                else
                   ele_volume = element_volume(xfield, agent%element)
                   call addto(release_diagfields(v)%ptr, agent%element, agent%biology(ivar)*agent%biology(BIOVAR_SIZE) / ele_volume)
                end if
             end do

             agent => agent%next
          end do
       end do

       ! Aggregate global requests from FG-level diagnostic fields
       do v=1, size(fgroup%ivars_uptake)
          var => fgroup%variables( fgroup%ivars_uptake(v) )
          call addto(uptake_fields(var%i_chemfield)%ptr, uptake_diagfields(v)%ptr)
       end do
       deallocate(uptake_diagfields)

       do v=1, size(fgroup%ivars_release)
          var => fgroup%variables( fgroup%ivars_release(v) )
          call addto(release_fields(var%i_chemfield)%ptr, release_diagfields(v)%ptr)
       end do
       deallocate(release_diagfields)
    end do

  end subroutine aggregate_chemical_diagnostics

  subroutine chemical_uptake(state)
    ! Handle uptake and release of chemicals
    type(state_type), intent(inout) :: state

    type(functional_group), pointer :: fgroup
    type(detector_linked_list), pointer :: agent_array
    type(scalar_field), pointer :: request_field, chemfield, depletion_field
    type(vector_field), pointer :: xfield
    type(detector_type), pointer :: agent
    integer :: i, j, n, ele, ingest_ind, fg, stage
    integer, dimension(:), pointer :: element_nodes
    real :: chemval, chemval_new, chem_integral, path_total, ele_volume, request
    real, dimension(1) :: depletion

    ! Exit if there are no uptake fields
    if (.not.associated(uptake_field_names)) return

    ewrite(2,*) "In chemical_uptake"

    xfield=>extract_vector_field(state, "Coordinate")

    ! Modify quantities on the chemical fields
    do i = 1, size(uptake_field_names)

       request_field=>extract_scalar_field(state, trim(uptake_field_names(i))//"Request")
       chemfield=>extract_scalar_field(state, trim(uptake_field_names(i)))
       depletion_field=>extract_scalar_field(state, trim(uptake_field_names(i))//"Depletion")

       ! Loop over all elements in chemical fields
       do ele=1,ele_count(chemfield)
          ele_volume = element_volume(xfield, ele)

          chem_integral = integral_element(chemfield, xfield, ele)
          request = integral_element(request_field, xfield, ele)

          ! Derive depletion factor
          if (request > chem_integral .and. request > 0.0) then 
             depletion(1) = chem_integral / request
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
                chemval_new = chemval * (1.0 - (request * depletion(1) / chem_integral))
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
    do fg=1, get_num_functional_groups()
       fgroup => get_functional_group(fg)
       do stage=1, size(fgroup%agent_arrays)
          agent_array => fgroup%agent_arrays(stage)

          if (have_option(trim(agent_array%stage_options)//"/biology")) then
             do j=1, size(fgroup%variables)

                if (fgroup%variables(j)%field_type == BIOFIELD_UPTAKE) then
                   depletion_field=>extract_scalar_field(state, trim(fgroup%variables(j)%chemfield)//"Depletion")

                   agent => agent_array%first
                   do while (associated(agent))
                      ingest_ind = fgroup%variables( fgroup%variables(j)%pool_index )%ingest_index

                      if (fgroup%variables(j)%path_integration) then
                         path_total = sum(agent%ele_dist)
                         do ele=1, size(agent%ele_path)
                            depletion = ele_val(depletion_field, agent%ele_path(ele))
                            agent%biology(ingest_ind) = agent%biology(ingest_ind) + depletion(1) * agent%biology(j) * agent%ele_dist(ele) / path_total
                         end do
                      else
                         depletion = ele_val(depletion_field, agent%element)
                         agent%biology(ingest_ind) = agent%biology(ingest_ind) + depletion(1) * agent%biology(j)
                      end if

                      agent => agent%next
                   end do
                end if
             end do
          end if
       end do
    end do

  end subroutine chemical_uptake

  subroutine chemical_release(state)
    ! Handle chemical release
    type(state_type), intent(inout) :: state

    type(scalar_field), pointer :: release_field, chemfield
    type(vector_field), pointer :: xfield
    type(detector_type), pointer :: agent
    integer :: i, j, n, ele, poolvar
    integer, dimension(:), pointer :: element_nodes
    real :: chemval, chemval_new, chem_integral, ele_volume, release

    ! Exit if there are no release fields
    if (.not.associated(release_field_names)) return

    ewrite(2,*) "In chemical_release"

    xfield=>extract_vector_field(state, "Coordinate")

    ! Modify quantities on the chemical fields
    do i = 1, size(release_field_names)

       release_field=>extract_scalar_field(state, trim(release_field_names(i))//"Release")
       chemfield=>extract_scalar_field(state, trim(release_field_names(i)))

       ! Loop over all elements in chemical fields
       do ele=1,ele_count(chemfield)
          chem_integral = integral_element(chemfield, xfield, ele)
          release = integral_element(release_field, xfield, ele)

          ! Set new chemical concentration
          element_nodes=>ele_nodes(chemfield, ele)
          do n=1, size(element_nodes)
             chemval = node_val(chemfield, element_nodes(n))
             ! Avoid div-by-zero error, and make sure not to write 0.0 into any field,
             ! since it will cause NaNs in the solvers. Instead use tiny()
             if (chem_integral > 0.0) then 
                ! C := C (1 + S/C_bar )
                chemval_new = chemval * (1.0 + (release / chem_integral))
                if (chemval_new > 0.0) then
                   call set(chemfield, element_nodes(n), chemval_new)
                else
                   call set(chemfield, element_nodes(n), tiny(1.0))
                end if
             else
                ! If there is zero chemical in this element 
                ! we turn the total quantity into an evenly distributed concentration
                ele_volume = element_volume(xfield, ele)
                if (release > 0.0) then
                   call set(chemfield, element_nodes(n), release / ele_volume)
                else
                   call set(chemfield, element_nodes(n), tiny(1.0))
                end if 
             end if
          end do
       end do
    end do

  end subroutine chemical_release

  subroutine aggregate_food_diagnostics(state, fgroup)
    type(state_type), intent(inout) :: state
    type(functional_group), intent(inout) :: fgroup

    type(scalar_field), pointer :: concentration
    type(scalar_field), pointer :: source_field
    type(scalar_field), pointer :: chem_field
    type(food_set) :: food
    type(le_variable) :: chem_var
    integer :: c, f, s

    ewrite(2,*) "Lagrangian biology: Aggregating food sets for FG::", fgroup%name

    do f=1, size(fgroup%food_sets)
       food = fgroup%food_sets(f)

       !concentration => extract_scalar_field(state, trim(food%conc_field_name))
       call zero(concentration)

       do s=1, size(food%varieties)
          source_field => extract_scalar_field(state, trim(food%varieties(s)%conc_field))
          call addto(concentration, source_field)
       end do

       ! Add the chemical concentrations across all target stages
       do c=1, size(food%ingest_chem_inds)
          chem_var = fgroup%variables( food%ingest_chem_inds(c) )
          chem_field => extract_scalar_field(state, trim(fgroup%name)//trim(food%name)//trim(chem_var%name) )
          call zero(chem_field)

          do s=1, size(food%varieties)
             source_field => extract_scalar_field(state, trim(food%target_fgroup)//"Particulate"//trim(chem_var%name)//trim(food%varieties(s)%name) )
             call addto(chem_field, source_field)
          end do
       end do
    end do

  end subroutine aggregate_food_diagnostics

  subroutine ingestion(state)
    type(state_type), intent(inout) :: state

    type(functional_group), pointer :: fgroup, target_fg
    type(detector_linked_list), pointer :: agent_array
    type(scalar_field), pointer :: conc_field, request_field, depletion_field
    type(scalar_field_pointer), dimension(:), allocatable :: prey_chem_fields
    type(vector_field), pointer :: xfield
    type(food_set), pointer :: fset 
    type(food_variety), pointer :: fvariety
    type(le_variable), pointer :: chempool_var
    type(detector_type), pointer :: agent
    real :: conc, request, chem_conc, prop, old_size
    real, dimension(1) :: depletion
    integer :: i, c, fs, fg, fv, t, ele, ingest_ind, stage

    ewrite(2,*) "Lagrangian_biology: Handling ingestion"

    xfield=>extract_vector_field(state, "Coordinate")

    do fg=1, get_num_functional_groups()
       fgroup => get_functional_group(fg)
       if (allocated(fgroup%food_sets)) then

          ! Derive the required Request fields
          if (size(fgroup%food_sets) > 0) then
             call ingestion_derive_requests(state, fgroup)
          end if

          do fs=1, size(fgroup%food_sets)
             fset => fgroup%food_sets(fs)
             do fv=1, size(fset%varieties)
                fvariety => fset%varieties(fv)
                target_fg => fvariety%target_list%ptr%fgroup

                conc_field => extract_scalar_field(state, fvariety%conc_field)
                request_field => extract_scalar_field(state, fvariety%vrequest%field_name)
                depletion_field => extract_scalar_field(state, trim(fvariety%vrequest%field_name)//"Depletion")       

                ! Loop over all elements in the source concentration field
                do ele=1,ele_count(conc_field)
                   conc = integral_element(conc_field, xfield, ele)
                   request = integral_element(request_field, xfield, ele)

                   ! Derive depletion factor
                   if (request > conc .and. request > 0.0) then 
                      depletion(1) = conc / request
                   else
                      depletion(1) = 1.0
                   end if
                   call set(depletion_field, ele, depletion(1))
                end do

                ! Loop over predator agents to set Ingest variables
                allocate(prey_chem_fields(size(fset%ingest_chem_inds)))
                do c=1, size(fset%ingest_chem_inds)
                   chempool_var => fgroup%variables( fset%ingest_chem_inds(c) )
                   prey_chem_fields(c)%ptr => extract_scalar_field(state, trim(fset%target_fgroup)//"Particulate"//trim(chempool_var%name)//trim(fvariety%name) )
                end do

                do stage=1, size(fgroup%agent_arrays)
                   agent_array => fgroup%agent_arrays(stage)

                   agent => agent_array%first
                   do while (associated(agent))
                      ! Set 'IngestedCells' variable
                      depletion = ele_val(depletion_field, agent%element)
                      agent%food_ingests(fv) = depletion(1) * agent%food_requests(fv)

                      ! Set ChemIngested pools
                      conc = integral_element(conc_field, xfield, agent%element)
                      do c=1, size(fset%ingest_chem_inds)
                         ingest_ind = fgroup%variables( fset%ingest_chem_inds(c) )%ingest_index
                         chem_conc = integral_element(prey_chem_fields(c)%ptr, xfield, agent%element)
                         if (conc > 0.0) then
                            agent%biology(ingest_ind) = agent%biology(ingest_ind) + ( chem_conc * (agent%food_ingests(fv) / conc) )
                         end if
                      end do

                      agent => agent%next
                   end do
                end do
                deallocate(prey_chem_fields)

                ! Loop over prey agents to adjust ensemble size
                agent => fvariety%target_list%ptr%first
                do while (associated(agent))
                   request = integral_element(request_field, xfield, agent%element)
                   if (request > 0.0) then     
                      conc = integral_element(conc_field, xfield, agent%element)
                      if (conc > 0.0) then                 
                         depletion = ele_val(depletion_field, agent%element)
                         old_size = agent%biology(BIOVAR_SIZE)

                         ! Proportion of food ingested
                         prop = old_size / conc
                         agent%biology(BIOVAR_SIZE) = old_size - (prop * request * depletion(1))
                      end if
                   end if

                   agent => agent%next
                end do

             end do  ! Food Variety
          end do  ! FoodSet

          ! Derive the resulting Ingested fields
          if (size(fgroup%food_sets) > 0) then
             call ingestion_set_ingests(state, fgroup)
          end if

       end if
    end do  ! FGroup

  end subroutine ingestion

  subroutine ingestion_derive_requests(state, fgroup)
    type(state_type), intent(inout) :: state
    type(functional_group), pointer, intent(inout) :: fgroup

    type(vector_field), pointer :: xfield
    type(scalar_field_pointer), dimension(size(fgroup%food_sets(1)%varieties)) :: frequest_fields
    type(food_variety), pointer :: variety
    type(detector_type), pointer :: agent
    real :: ele_volume
    integer :: i, v, vsize

    xfield=>extract_vector_field(state, "Coordinate")

    ! Pull requests fields and zero them
    vsize = size(fgroup%food_sets(1)%varieties)
    do v=1, vsize
       variety => fgroup%food_sets(1)%varieties(v)

       frequest_fields(v)%ptr => extract_scalar_field(state, trim(variety%vrequest%field_name))
       call zero(frequest_fields(v)%ptr)
    end do

    do i=1, size(fgroup%agent_arrays)       
       agent => fgroup%agent_arrays(i)%first
       do while (associated(agent))
          do v=1, vsize
             variety => fgroup%food_sets(1)%varieties(v)

             if (variety%vrequest%path_integration) then
                call integrate_along_path(frequest_fields(v)%ptr, xfield, agent, agent%food_requests(v))
             else
                ele_volume = element_volume(xfield, agent%element)
                call addto(frequest_fields(v)%ptr, agent%element, agent%food_requests(v)*agent%biology(BIOVAR_SIZE) / ele_volume)
             end if
          end do
          agent => agent%next
       end do
    end do
  end subroutine ingestion_derive_requests

  subroutine ingestion_set_ingests(state, fgroup)
    type(state_type), intent(inout) :: state
    type(functional_group), pointer, intent(inout) :: fgroup

    type(vector_field), pointer :: xfield
    type(scalar_field_pointer), dimension(size(fgroup%food_sets(1)%varieties)) :: fingest_fields
    type(food_variety), pointer :: variety
    type(detector_type), pointer :: agent
    real :: ele_volume
    integer :: i, v, vsize

    xfield=>extract_vector_field(state, "Coordinate")

    ! Pull requests fields and zero them
    vsize = size(fgroup%food_sets(1)%varieties)
    do v=1, vsize
       variety => fgroup%food_sets(1)%varieties(v)

       fingest_fields(v)%ptr => extract_scalar_field(state, trim(variety%vingest%field_name))
       call zero(fingest_fields(v)%ptr)
    end do

    do i=1, size(fgroup%agent_arrays)       
       agent => fgroup%agent_arrays(i)%first
       do while (associated(agent))
          do v=1, vsize
             variety => fgroup%food_sets(1)%varieties(v)

             if (variety%vingest%path_integration) then
                call integrate_along_path(fingest_fields(v)%ptr, xfield, agent, agent%food_ingests(v))
             else
                ele_volume = element_volume(xfield, agent%element)
                call addto(fingest_fields(v)%ptr, agent%element, agent%food_ingests(v)*agent%biology(BIOVAR_SIZE) / ele_volume)
             end if
          end do
          agent => agent%next
       end do
    end do
  end subroutine ingestion_set_ingests

  subroutine integrate_along_path(diagfield, xfield, agent, agent_variable)
    type(scalar_field), pointer, intent(inout) :: diagfield
    type(vector_field), pointer, intent(inout) :: xfield
    type(detector_type), pointer, intent(inout) :: agent
    real, intent(in) :: agent_variable

    real :: path_total, quantity, ele_path_volume
    integer :: ele

    ! Integrate along the path of the agent
    if (allocated(agent%ele_path)) then                   
       path_total = sum(agent%ele_dist)
       quantity = agent_variable * agent%biology(BIOVAR_SIZE)

       do ele=1, size(agent%ele_path)
          ele_path_volume = element_volume(xfield, agent%ele_path(ele))
          call addto(diagfield, agent%ele_path(ele), (agent%ele_dist(ele) / path_total) * (quantity / ele_path_volume))
       end do
    else
       FLAbort("Cannot integrate along path without a recorded path!")
    end if
  end subroutine integrate_along_path

  function insert_global_uptake_field(field_name) result (index)
    character(len=FIELD_NAME_LEN), intent(in) :: field_name

    character(len=FIELD_NAME_LEN), dimension(:), pointer :: tmp_list
    integer :: i, old_size, index

    if (.not. associated(uptake_field_names)) then
       allocate(uptake_field_names(1))
       uptake_field_names(1) = trim(field_name)
       index = 1
       return
    end if

    do i=1, size(uptake_field_names)
      if (trim(uptake_field_names(i)) == trim(field_name)) then
         index = i
         return
      end if
    end do

    old_size=size(uptake_field_names)
    tmp_list=>uptake_field_names
    allocate(uptake_field_names(old_size+1))
    uptake_field_names(1:old_size)=tmp_list
    uptake_field_names(old_size+1)=trim(field_name)
    deallocate(tmp_list)

    index = old_size+1
  end function insert_global_uptake_field

  function insert_global_release_field(field_name) result (index)
    character(len=FIELD_NAME_LEN), intent(in) :: field_name

    character(len=FIELD_NAME_LEN), dimension(:), pointer :: tmp_list
    integer :: i, old_size, index

    if (.not. associated(release_field_names)) then
       allocate(release_field_names(1))
       release_field_names(1) = trim(field_name)
       index = 1
       return
    end if

    do i=1, size(release_field_names)
      if (trim(release_field_names(i)) == trim(field_name)) then
         index = i
         return
      end if
    end do

    old_size=size(release_field_names)
    tmp_list=>release_field_names
    allocate(release_field_names(old_size+1))
    release_field_names(1:old_size)=tmp_list
    release_field_names(old_size+1)=trim(field_name)
    deallocate(tmp_list)

    index = old_size+1
  end function insert_global_release_field

end module lagrangian_biology
