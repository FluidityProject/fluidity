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
  use diagnostic_fields
  use Profiler

implicit none

  private

  public :: initialise_lagrangian_biology, lagrangian_biology_cleanup, &
            update_lagrangian_biology, calculate_agent_diagnostics

  type(detector_linked_list), dimension(:), allocatable, target, save :: agent_arrays

  integer, parameter :: BIOFIELD_NONE=0, BIOFIELD_DIAG=1, BIOFIELD_UPTAKE=2, BIOFIELD_RELEASE=3

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
               ierror, det_type, biovar_index, biovar_total, biovar_state, &
               biovar_chemical, biovar_uptake, biovar_release, index, rnd_dim, n_env_fields
    integer, dimension(1) :: rnd_seed

    if (.not.have_option("/embedded_models/lagrangian_ensemble_biology")) return

    ewrite(1,*) "In initialise_lagrangian_biology"

    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)
    xfield=>extract_vector_field(state, "Coordinate")
    shape=>ele_shape(xfield,1)

    ! Determine how many arrays we need across all functional groups
    n_agent_arrays = 0
    n_fgroups = option_count("/embedded_models/lagrangian_ensemble_biology/functional_group")
    do fg=1, n_fgroups
       write(fg_buffer, "(a,i0,a)") "/embedded_models/lagrangian_ensemble_biology/functional_group[",fg-1,"]"

       n_agent_arrays = n_agent_arrays + option_count(trim(fg_buffer)//"/stage_array")
    end do

    allocate(agent_arrays(n_agent_arrays))
    ewrite(2,*) "Found a total of ", n_agent_arrays, " agent arrays in ", n_fgroups, " functional groups"

    ! Create a persistent dictionary of FG variable name mappings
    call python_run_string("persistent['fg_var_names'] = dict()")
    call python_run_string("persistent['fg_env_names'] = dict()")

    ! We create an agent array for each stage of each Functional Group
    do fg=1, n_fgroups
       write(fg_buffer, "(a,i0,a)") "/embedded_models/lagrangian_ensemble_biology/functional_group[",fg-1,"]"
       call get_option(trim(fg_buffer)//"/name", fg_name)
       n_fg_arrays = option_count(trim(fg_buffer)//"/stage_array")

       do i = 1, n_fg_arrays
          write(stage_buffer, "(a,i0,a)") trim(fg_buffer)//"/stage_array[",i-1,"]"
          call get_option(trim(stage_buffer)//"/number_of_agents", n_agents)
          call get_option(trim(stage_buffer)//"/name", stage_name)
          agent_arrays(i)%name = trim(fg_name)//trim(stage_name)
          agent_arrays(i)%fg_name = trim(fg_name)
          agent_arrays(i)%stage_name = trim(stage_name)
          agent_arrays(i)%fg_id = fg
          call get_option(trim(stage_buffer)//"/id", agent_arrays(i)%stage_id)

          ! Register the agent array, so Zoltan/Adaptivity will not forget about it
          call register_detector_list(agent_arrays(i))
          agent_arrays(i)%total_num_det=n_agents

          ! Get options for lagrangian detector movement
          call read_detector_move_options(agent_arrays(i),trim(fg_buffer))

          ! Get options for Random Walk
          if (have_option(trim(stage_buffer)//"/random_walk")) then
             agent_arrays(i)%move_parameters%do_random_walk=.true.
             call get_option("/embedded_models/lagrangian_ensemble_biology/random_seed", rnd_seed(1))

             ! Initialise random number generator in Python
             call python_run_string("numpy.random.seed("//trim(int2str(rnd_seed(1)))//")")
             if (have_option(trim(stage_buffer)//"/random_walk/python")) then 
                call get_option(trim(stage_buffer)//"/random_walk/python", agent_arrays(i)%move_parameters%rw_pycode)
             end if

             if (have_option(trim(stage_buffer)//"/random_walk/diffusive_random_walk")) then 
                agent_arrays(i)%move_parameters%use_internal_rw=.true.
                call get_option(trim(stage_buffer)//"/random_walk/diffusive_random_walk/diffusivity_field", &
                       agent_arrays(i)%move_parameters%diffusivity_field)
                call get_option(trim(stage_buffer)//"/random_walk/diffusive_random_walk/diffusivity_gradient", &
                       agent_arrays(i)%move_parameters%diffusivity_grad)
                ! Initialise random number generator
                rnd_dim=1
                call random_seed(size=rnd_dim)
                call randoM_seed(put=rnd_seed(1:rnd_dim))
             end if

          else
             agent_arrays(i)%move_parameters%do_random_walk=.false.
          end if

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
          if (have_option(trim(stage_buffer)//"/debug/exclude_from_advection")) then
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

          ! Set agent%list_id, because we need it to detect stage changes
          agent=>agent_arrays(i)%first
          do while (associated(agent))
             agent%list_id=agent_arrays(i)%id
             agent=>agent%next
          end do

          ewrite(2,*) "Found", agent_arrays(i)%length, "local agents in array", i

          !!! Biology setup !!!

          ! First determine number of biovars
          biovar_state = option_count(trim(fg_buffer)//"/state_variable") 
          biovar_chemical = option_count(trim(fg_buffer)//"/chemical_variable")
          biovar_uptake = 0
          biovar_release = 0
          do j=1, biovar_chemical
             write(biovar_buffer, "(a,i0,a)") trim(fg_buffer)//"/chemical_variable[",j-1,"]"
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
             allocate(agent_arrays(i)%biovar_name(biovar_total))
             allocate(agent_arrays(i)%biofield_type(biovar_total))
             allocate(agent_arrays(i)%biofield_dg_name(biovar_total))
             allocate(agent_arrays(i)%biofield_cg_name(biovar_total))
             allocate(agent_arrays(i)%chemfield_name(biovar_total))

             ! Add internal state variables
             index = 1
             do j=1, biovar_state
                write(biovar_buffer, "(a,i0,a)") trim(fg_buffer)//"/state_variable[",j-1,"]"
                call get_option(trim(biovar_buffer)//"/name", agent_arrays(i)%biovar_name(index))

                ! Record according diagnostic field
                if (have_option(trim(biovar_buffer)//"/scalar_field")) then
                   agent_arrays(i)%biofield_type(index) = BIOFIELD_DIAG
                   call get_option(trim(biovar_buffer)//"/scalar_field/name", field_name)
                   if (have_option(trim(biovar_buffer)//"/scalar_field/per_stage")) then
                      agent_arrays(i)%biofield_dg_name(index) = trim(fg_name)//trim(field_name)//trim(stage_name)
                   else
                      agent_arrays(i)%biofield_dg_name(index) = trim(fg_name)//trim(field_name)
                   end if
                else
                   agent_arrays(i)%biofield_type(index) = BIOFIELD_NONE
                end if
                index = index+1
             end do

             ! Add chemical variables
             do j=1, biovar_chemical
                write(biovar_buffer, "(a,i0,a)") trim(fg_buffer)//"/chemical_variable[",j-1,"]"
                call get_option(trim(biovar_buffer)//"/name", biovar_name)

                ! Chemical pool variable and according diagnostic fields
                agent_arrays(i)%biovar_name(index) = trim(biovar_name)//"Pool"
                if (have_option(trim(biovar_buffer)//"scalar_field")) then
                   agent_arrays(i)%biofield_type(index) = BIOFIELD_DIAG
                   call get_option(trim(biovar_buffer)//"/scalar_field/name", field_name)
                   if (have_option(trim(biovar_buffer)//"/scalar_field/per_stage")) then
                      agent_arrays(i)%biofield_dg_name(index) = trim(fg_name)//trim(field_name)//trim(biovar_name)//trim(stage_name)
                   else
                      agent_arrays(i)%biofield_dg_name(index) = trim(fg_name)//trim(field_name)//trim(biovar_name)
                   end if
                else
                   agent_arrays(i)%biofield_type(index) = BIOFIELD_NONE
                end if
                index = index+1

                ! Chemical uptake
                if (have_option(trim(biovar_buffer)//"/uptake")) then
                   if (.not.have_option(trim(biovar_buffer)//"/chemical_field")) then
                      FLExit("No chemical field specified for "//trim(biovar_name)//" uptake in functional group "//trim(fg_name))
                   end if
                   agent_arrays(i)%biovar_name(index) = trim(biovar_name)//"Ing"
                   agent_arrays(i)%biofield_type(index) = BIOFIELD_UPTAKE
                   agent_arrays(i)%biofield_dg_name(index) = trim(fg_name)//"DGRequest"//trim(biovar_name)
                   agent_arrays(i)%biofield_cg_name(index) = trim(fg_name)//"CGRequest"//trim(biovar_name)
                   call get_option(trim(biovar_buffer)//"/chemical_field/name", agent_arrays(i)%chemfield_name(index))
                   index = index+1
                end if

                ! Chemical release
                if (have_option(trim(biovar_buffer)//"/release")) then
                   if (.not.have_option(trim(biovar_buffer)//"/chemical_field")) then
                      FLExit("No chemical field specified for "//trim(biovar_name)//" release in functional group "//trim(fg_name))
                   end if
                   agent_arrays(i)%biovar_name(index) = trim(biovar_name)//"Rel"
                   agent_arrays(i)%biofield_type(index) = BIOFIELD_RELEASE
                   agent_arrays(i)%biofield_dg_name(index) = trim(fg_name)//"DGRelease"//trim(biovar_name)
                   agent_arrays(i)%biofield_cg_name(index) = trim(fg_name)//"CGRelease"//trim(biovar_name)
                   call get_option(trim(biovar_buffer)//"/chemical_field/name", agent_arrays(i)%chemfield_name(index))
                   index = index+1
                end if

             end do

             ! Store the python update code
             if (have_option(trim(stage_buffer)//"/biology_update")) then
                call get_option(trim(stage_buffer)//"/biology_update", agent_arrays(i)%biovar_pycode)
                agent_arrays(i)%has_biology=.true.
             end if

             ! Now we know the variable names for the FG, so we populate the mapping dict
             call python_run_string("persistent['fg_var_names']['"//trim(agent_arrays(i)%name)//"'] = []")
             do j=1, biovar_total
                call python_run_string("persistent['fg_var_names']['"//trim(agent_arrays(i)%name)//"'].append('"//trim(agent_arrays(i)%biovar_name(j))//"')")
             end do

             ! Initialise agent variables
             agent => agent_arrays(i)%first
             do while (associated(agent))
                allocate(agent%biology(biovar_total))
                call get_option(trim(stage_buffer)//"/initial_state/values", agent%biology)
                agent => agent%next
             end do

             ! Record which environment fields to evaluate and pass to the update function
             n_env_fields = option_count(trim(fg_buffer)//"/environment_field")
             call python_run_string("persistent['fg_env_names']['"//trim(agent_arrays(i)%name)//"'] = []")
             allocate(agent_arrays(i)%env_field_name(n_env_fields))
             do j=1, n_env_fields
                write(env_field_buffer, "(a,i0,a)") trim(fg_buffer)//"/environment_field[",j-1,"]"
                call get_option(trim(env_field_buffer)//"/name", agent_arrays(i)%env_field_name(j))
                call python_run_string("persistent['fg_env_names']['"//trim(agent_arrays(i)%name)//"'].append('"//trim(agent_arrays(i)%env_field_name(j))//"')")
             end do
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

    ! Output agent positions
    do i = 1, size(agent_arrays)
       call write_detectors(state, agent_arrays(i), time, dt)
    end do

    call profiler_toc("/update_lagrangian_biology")

  end subroutine update_lagrangian_biology

  subroutine calculate_agent_diagnostics(state)
    type(state_type), intent(inout) :: state

    type(scalar_field), pointer :: sfield, request_field_dg, request_field_cg, chemical_field, absorption_field, &
                                   depletion_field, release_field_dg, release_field_cg, source_field
    type(vector_field), pointer :: xfield
    integer :: i, j, n, current_fg

    ewrite(1,*) "In calculate_agent_diagnostics"
    call profiler_tic("/calculate_agent_diagnostics")

    xfield=>extract_vector_field(state, "Coordinate")

    do i = 1, size(agent_arrays)
       if (agent_arrays(i)%has_biology) then

          sfield=>extract_scalar_field(state, trim(agent_arrays(i)%fg_name)//"Agents"//trim(agent_arrays(i)%stage_name))
          call zero(sfield)

          ! Reset the diagnostic fields associated with the agent array
          do j=1, size(agent_arrays(i)%biovar_name)
             if (agent_arrays(i)%biofield_type(j) /= BIOFIELD_NONE) then
                sfield=>extract_scalar_field(state, trim(agent_arrays(i)%biofield_dg_name(j)))
                call zero(sfield)
             end if
             if (agent_arrays(i)%biofield_type(j) == BIOFIELD_UPTAKE) then
                sfield=>extract_scalar_field(state, trim(agent_arrays(i)%chemfield_name(j))//"Absorption")
                call zero(sfield)
             end if
             if (agent_arrays(i)%biofield_type(j) == BIOFIELD_RELEASE) then
                sfield=>extract_scalar_field(state, trim(agent_arrays(i)%chemfield_name(j))//"Source")
                call zero(sfield)
             end if
          end do
       end if
    end do

    do i = 1, size(agent_arrays)
       if (agent_arrays(i)%has_biology) then
          sfield=>extract_scalar_field(state, trim(agent_arrays(i)%fg_name)//"Agents"//trim(agent_arrays(i)%stage_name))
          call addto_agent_count(agent_arrays(i), sfield)

          do j=1, size(agent_arrays(i)%biovar_name)
             if (agent_arrays(i)%biofield_type(j) /= BIOFIELD_NONE) then
                sfield=>extract_scalar_field(state, trim(agent_arrays(i)%biofield_dg_name(j)))
                call addto_diagnostic_field_from_agents(agent_arrays(i), j, sfield, xfield)
             end if
          end do
       end if
    end do

    current_fg = 1
    do i = 1, size(agent_arrays)
       ! We only want to do this for each functional group
       if (agent_arrays(i)%has_biology.and.agent_arrays(i)%fg_id>=current_fg) then
          do j=1, size(agent_arrays(i)%biovar_name)
             ! Handle chemical uptake
             if (agent_arrays(i)%biofield_type(j) == BIOFIELD_UPTAKE) then

                request_field_dg=>extract_scalar_field(state, trim(agent_arrays(i)%biofield_dg_name(j)))
                request_field_cg=>extract_scalar_field(state, trim(agent_arrays(i)%biofield_cg_name(j)))
                chemical_field=>extract_scalar_field(state, trim(agent_arrays(i)%chemfield_name(j)))
                absorption_field=>extract_scalar_field(state, trim(agent_arrays(i)%chemfield_name(j))//"Absorption")
                depletion_field=>extract_scalar_field(state, trim(agent_arrays(i)%chemfield_name(j))//"Depletion")

                ewrite(2,*) "Galerkin projecting request field: ", trim(agent_arrays(i)%biofield_dg_name(j))
                call profiler_tic("/lagrangian_biology_galerkin_projection")
                call calculate_galerkin_projection(state, request_field_cg)
                call profiler_toc("/lagrangian_biology_galerkin_projection")

                do n=1, node_count(request_field_cg)
                   if (node_val(request_field_cg,n) > node_val(chemical_field,n)) then
                      ! Scale back the request
                      call set(depletion_field, n, node_val(chemical_field,n) / node_val(request_field_cg,n))
                   else
                      call set(depletion_field, n, 1.0)
                   end if

                   call set(absorption_field, n, node_val(request_field_cg,n) * node_val(depletion_field,n))
                end do
             end if

             ! Handle chemical release
             if (agent_arrays(i)%biofield_type(j) == BIOFIELD_RELEASE) then
                release_field_dg=>extract_scalar_field(state, trim(agent_arrays(i)%biofield_dg_name(j)))
                release_field_cg=>extract_scalar_field(state, trim(agent_arrays(i)%biofield_cg_name(j)))
                source_field=>extract_scalar_field(state, trim(agent_arrays(i)%chemfield_name(j))//"Source")

                ewrite(2,*) "Galerkin projecting release field: ", trim(agent_arrays(i)%biofield_dg_name(j))
                call profiler_tic("/lagrangian_biology_galerkin_projection")
                call calculate_galerkin_projection(state, release_field_cg)
                call profiler_toc("/lagrangian_biology_galerkin_projection")

                do n=1, node_count(release_field_cg)
                   call set(source_field, n, node_val(release_field_cg,n))
                end do
             end if
          end do

          current_fg = current_fg + 1
       end if
    end do

    ewrite(2,*) "Exiting calculate_agent_diagnostics"
    call profiler_toc("/calculate_agent_diagnostics")

  end subroutine calculate_agent_diagnostics

  subroutine addto_diagnostic_field_from_agents(agent_list, biovar_id, sfield, xfield)
    type(detector_linked_list), intent(inout) :: agent_list
    integer, intent(in) :: biovar_id
    type(scalar_field), pointer, intent(inout) :: sfield
    type(vector_field), pointer, intent(inout) :: xfield

    type(detector_type), pointer :: agent
    type(element_type), pointer :: shape
    integer, dimension(ele_loc(sfield, agent_list%first%element)) :: element_nodes
    real :: value, scaled_value, ele_volume, shape_val
    integer :: i, j

    ewrite(2,*) "In addto_diagnostic_field_from_agents: ", sfield%name

    agent => agent_list%first
    do while (associated(agent))
       value = agent%biology(biovar_id)
       ele_volume = element_volume(xfield, agent%element)
       element_nodes = ele_nodes(sfield, agent%element)
       scaled_value = value / ele_volume

       do i=1, ele_loc(sfield, agent%element)
          call addto(sfield, element_nodes(i), scaled_value)
       end do

       agent => agent%next
    end do

  end subroutine addto_diagnostic_field_from_agents

  subroutine addto_agent_count(agent_list, sfield)
    type(detector_linked_list), intent(inout) :: agent_list
    type(scalar_field), pointer, intent(inout) :: sfield

    type(detector_type), pointer :: agent
    integer, dimension(ele_loc(sfield, agent_list%first%element)) :: element_nodes
    integer :: i

    agent => agent_list%first
    do while (associated(agent))
       element_nodes = ele_nodes(sfield, agent%element)
       do i=1, ele_loc(sfield, agent%element)
          call addto(sfield, element_nodes(i), 1.0)
       end do

       agent => agent%next
    end do

  end subroutine addto_agent_count

end module lagrangian_biology
