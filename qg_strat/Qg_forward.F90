#include "fdebug.h"
module qg_forward
  use populate_state_module
  use initialise_fields_module
  use fields_manipulation
  use state_module
  use fields
  use FLDebug
  use form_metric_field
  use global_parameters, only: OPTION_PATH_LEN
  use metric_assemble
  use adapt_state_module 
  use spud
  use sparse_tools
  use sparsity_patterns
  use boundary_conditions
  use advection_diffusion_cg
  use vtk_interfaces
  use global_parameters, only:FIELD_NAME_LEN, global_debug_level
  use pv_inversion
  use timeloop_utilities
  implicit none

CONTAINS

  subroutine qg_forward_run()
    !!< Integrate the forward model for the quasi-geostrophic equations 
    type(state_type), pointer, dimension(:) :: state
    type(state_type), pointer, dimension(:) :: bc_states
    real :: t
    real :: dt
    real :: dump_period
    real :: tdump
    real :: adapt_period
    real :: tadapt
    integer :: dump_count=1
    !
    type(scalar_field), pointer :: PV, PV_Old, streamfunction, time
    type(scalar_field), pointer :: buoyancy, buoyancy_old
    type(vector_field), pointer :: X, beta
    type(csr_sparsity), pointer :: pv_sparsity
    type(csr_sparsity), pointer :: buoyancy_sparsity
    type(mesh_type), pointer :: mesh
    type(tensor_field) :: metric
    integer :: nstages, stage
    integer :: i, n_buoyancy_bcs
    integer, dimension(2) :: shape_option

    character(len=FIELD_NAME_LEN) :: bc_name
    character(len=FIELD_NAME_LEN) :: filename=""

    call get_option("/simulation_name",filename)

    ! Get timestepping options
    call get_option("/timestepping/current_time", t)
    call get_option("/timestepping/timestep", dt)
    call get_option("/io/dump_period", dump_period)

    call populate_state(state)

    call allocate_and_insert_auxilliary_qg_fields(state(1))

    ! Setup adaptivity if used
    if(have_option("/mesh_adaptivity")) then
       mesh => extract_mesh(state(1), "CoordinateMesh")
       call allocate(metric, mesh, "ErrorMetric")
    end if

    if(have_option("/mesh_adaptivity/hr_adaptivity/adapt_at_first_timestep")) &
         then
       ewrite(1,*) "adapting mesh at first timestep"
       call assemble_metric(state, metric)
       call adapt_state(state, metric)
       call allocate_and_insert_auxilliary_qg_fields(state(1))
    end if

    PV => extract_scalar_field(state(1), 'PotentialVorticity')
    PV_Old => extract_scalar_field(state(1), 'OldPotentialVorticity')
    streamfunction => extract_scalar_field(state(1),'Streamfunction')
    X => extract_vector_field(state(1), 'Coordinate')
    
    ! Solve for streamfunction if it is not prescribed
    n_buoyancy_bcs=0
    if(have_option("/material_phase[0]/scalar_field::Streamfunction&
         /prognostic")) then
       ! set up bc_states if prognostic streamfunction has buoyancy
       ! boundary conditions
       n_buoyancy_bcs=option_count("/material_phase["//int2str(0)//&
            "]/scalar_field::Streamfunction/prognostic/boundary_conditions/&
            type::buoyancy")
       if(n_buoyancy_bcs>=1) then
          allocate(bc_states(1:n_buoyancy_bcs))
          do i=1, n_buoyancy_bcs
             call nullify(bc_states(i))
          end do
          call construct_2d_state_for_boundary_conditions(state(1), bc_states)
       end if
       ! Can now solve for streamfunction
       call solve_streamfunction_qg(state(1))
    end if

    ! Always calculate velocity from streamfunction
    call streamfunction2velocity(state(1))

    ! We have now allocated, inserted and calculated all fields so can
    ! write out state at t=0
    call vtk_write_state(trim(filename), index=0, state=(/state/))
    if(n_buoyancy_bcs>=1) then
       call vtk_write_state(trim(filename)//"_bcs", index=0, state=(/bc_states/))
    end if

    tdump = 0.

    timestep_loop: do
       if(simulation_completed(t)) exit

       t = t + dt
       !update time in options dictionary
       call set_option("/timestepping/current_time", t)
       !update time field in state (for writing to vtus)
       time => extract_scalar_field(state(1), 'Time')
       call set(time, t)

       tdump = tdump + dt
       tadapt = tadapt + dt

       ! for time dependent prescribed fields
       if(have_option("/material_phase[0]/scalar_field::Streamfunction&
         /prescribed")) then
          call set_prescribed_field_values(state)
          call streamfunction2velocity(state(1))
       end if

       ! Adapt mesh if required
       if(have_option("/mesh_adaptivity").and.tadapt>=adapt_period) then
          ewrite(1,*) "adapting mesh"
          mesh => extract_mesh(state(1), "CoordinateMesh")
          call allocate(metric, mesh, "ErrorMetric")
          call assemble_metric(state, metric)
          call adapt_state(state, metric)
          ! Reallocate and insert stuff into state:
          call allocate_and_insert_auxilliary_qg_fields(state(1))
          tadapt=0.
          ! Re-extract fields and recalculate velocity from streamfunction
          PV => extract_scalar_field(state(1), 'PotentialVorticity')
          PV_Old => extract_scalar_field(state(1), 'OldPotentialVorticity')
          streamfunction => extract_scalar_field(state(1),'Streamfunction')
          X => extract_vector_field(state(1), 'Coordinate')
          call streamfunction2velocity(state(1))
       end if

       call set(PV_old, PV)

       pv_sparsity => extract_csr_sparsity(state(1), &
            'PotentialVorticitySparsity')
       call solve_field_equation_cg('PotentialVorticity', &
            state(1), dt, velocity_name='GeostrophicVelocity')
       PV => extract_scalar_field(state(1), "PotentialVorticity")

       ! If streamfunction is prognostic, solve for it and 
       ! recalculate velocity
       if(have_option("/material_phase[0]/scalar_field::&
            Streamfunction/prognostic")) then
          ! If we are advecting buoyancy on the top and bottom
          ! surfaces, we must solve for it first:
          if(n_buoyancy_bcs>=1) then
             call update_2d_state_for_boundary_conditions(state(1), bc_states)
             do i=0, n_buoyancy_bcs-1
                buoyancy=>extract_scalar_field(bc_states(i+1), "Buoyancy")
                buoyancy_old=>extract_scalar_field(bc_states(i+1), &
                     "OldBuoyancy")
                call set(buoyancy_old, buoyancy)
                buoyancy_sparsity=>extract_csr_sparsity(bc_states(i+1), &
                     "BuoyancySparsity")
                call set(buoyancy_old, buoyancy)
                call solve_field_equation_cg('Buoyancy', bc_states(i+1), &
                     dt, velocity_name='GeostrophicVelocity')
                call get_option("/material_phase["//int2str(0)//&
                     "]/scalar_field::Streamfunction/prognostic/&
                     boundary_conditions["//int2str(i)//"]/name",&
                     bc_name)
                call insert_surface_field(streamfunction, trim(bc_name), &
                     buoyancy)
             end do
          end if
          call solve_streamfunction_qg(state(1))
          call streamfunction2velocity(state(1))
       end if

       if(tdump>=dump_period) then 
          ewrite(1,*) 'dump', dump_count
          tdump = 0.

          call vtk_write_state(trim(filename), dump_count, state=state)
          if(n_buoyancy_bcs>=1) then
             call vtk_write_state(trim(filename)//"_bcs", dump_count, state=bc_states)
          end if
          dump_count = dump_count + 1
       end if

    end do timestep_loop

    do i=1, size(state)
       call deallocate(state(i))
    end do

    if(n_buoyancy_bcs>=1) then
       do i=1, size(bc_states)
          call deallocate(bc_states(i))
       end do
    end if

    ewrite(0,*) "Printing references at the end of Qg_forward "
    ewrite(0,*) "- there shouldn't be any..."
    call print_references(0)

  end subroutine qg_forward_run

  subroutine allocate_and_insert_auxilliary_qg_fields(state)
    !!< This subroutine allocates and inserts into state 
    !!< fields that are unique to qg_strat

    type(state_type), intent(inout) :: state

    type(scalar_field), pointer :: PV, streamfunction
    type(scalar_field) :: OldPV, Time
    type(vector_field), pointer :: X
    type(vector_field) :: V, beta
    type(csr_sparsity) :: pv_sparsity
    type(mesh_type) :: V_mesh

    integer :: stage

    ewrite(1,*) "In allocate_and_insert_auxilliary_qg_fields"

    PV => extract_scalar_field(state, 'PotentialVorticity')
    streamfunction => extract_scalar_field(state, 'Streamfunction')
    X => extract_vector_field(state, 'Coordinate')

    ! allocate, zero and insert OldPotentialVorticity
    call allocate(OldPV, PV%mesh, 'OldPotentialVorticity')
    call zero(OldPV)
    call insert(state, OldPV, 'OldPotentialVorticity')
    call deallocate(OldPV)

    ! construct mesh for velocity field (it's DG because it is the skew
    ! gradient of a streamfunction)
    V_mesh=make_mesh(streamfunction%mesh, streamfunction%mesh%shape, continuity=-1, &
         name='VelocityMesh')
    call allocate(V, mesh_dim(X), V_mesh, 'GeostrophicVelocity')
    call zero(V)
    call insert(state, V, 'GeostrophicVelocity')
    call deallocate(V_mesh)
    call deallocate(V)

    call allocate(beta, mesh_dim(X), X%mesh, 'Beta')
    call zero(beta)
    if(have_option("/physical_parameters/coriolis/beta_plane")) then
       call initialise_field(beta, "/physical_parameters/coriolis/beta_plane/beta/vector_field/prescribed/value/", X)
    end if
    call insert(state, beta, 'Beta')
    call deallocate(beta)

    pv_sparsity=make_sparsity(PV%mesh, PV%mesh,  &
         "PotentialVorticitySparsity")
    call insert(state, pv_sparsity, "PotentialVorticitySparsity")
    call deallocate(pv_sparsity)    

    ! Ensure that time is output in vtu files.
    call allocate(Time, X%mesh, 'Time', field_type=FIELD_TYPE_CONSTANT)
    call get_option('/timestepping/current_time', current_time)
    call set(Time, current_time)
    Time%option_path = ""
    call insert(state, Time, 'Time')
    call deallocate(Time)


  end subroutine allocate_and_insert_auxilliary_qg_fields

  subroutine construct_2d_state_for_boundary_conditions(state, bc_states)

    type(state_type), intent(in) :: state
    type(state_type), dimension(:), intent(inout) ::bc_states

    type(mesh_type), pointer :: surface_mesh
    type(mesh_type) :: surface_v_mesh
    type(scalar_field), pointer :: streamfunction
    type(scalar_field) :: buoyancy, OldBuoyancy
    type(vector_field), pointer :: position, velocity
    type(vector_field) :: tmp_surface_position, surface_position
    type(vector_field) :: tmp_surface_velocity, surface_velocity
    type(tensor_field) :: diffusivity
    type(csr_sparsity), allocatable, dimension(:) :: buoyancy_sparsity

    integer :: i, n_bcs, n_buoyancy_bcs
    integer, pointer, dimension(:) :: surface_element_list
    character(len=FIELD_NAME_LEN) :: bc_path, bc_name

    ewrite(1,*) "In construct_2d_state_for_boundary_conditions"

    ! Extract required fields from state
    streamfunction=>extract_scalar_field(state, 'Streamfunction') 
    position=>extract_vector_field(state, 'Coordinate')
    velocity=>extract_vector_field(state, 'GeostrophicVelocity')

    ! Get number of boundary conditions
    n_bcs=option_count("/material_phase["//int2str(0)//&
         "]/scalar_field::Streamfunction/prognostic/boundary_conditions")

    ! Get number of buoyancy boundary conditions
    n_buoyancy_bcs=option_count("/material_phase["//int2str(0)//&
         "]/scalar_field::Streamfunction/prognostic/&
         boundary_conditions/type::buoyancy")

    ! Allocate sparsity array
    allocate(buoyancy_sparsity(1:n_buoyancy_bcs))

    ! Loop over boundary conditions
    do i=0,n_bcs-1

       bc_path="/material_phase["//int2str(0)//&
            "]/scalar_field::Streamfunction/prognostic/boundary_conditions["&
            //int2str(i)//"]"

       if(have_option(trim(bc_path)//"/type::buoyancy")) then

          ! get boundary condition name
          call get_option(trim(bc_path)//"/name", bc_name)
          bc_states(i+1)%name=trim(bc_name)
          ewrite(2,*) "boundary condition name is: ", trim(bc_name)

          ! get surface_element_list and surface mesh
          call get_boundary_condition(streamfunction, bc_name, &
               surface_mesh=surface_mesh, &
               surface_element_list=surface_element_list)

          ! put surface mesh into bc_state(i)
          call insert(bc_states(i+1), surface_mesh, "surface_mesh")

          ! get a positions field on this surface
          call allocate(tmp_surface_position, position%dim, surface_mesh, &
               "tmp_surface_position")
          call remap_field_to_surface(position, tmp_surface_position, &
               surface_element_list)
          call allocate(surface_position, position%dim-1, surface_mesh, &
               'Coordinate')
          surface_position%val(1)%ptr=tmp_surface_position%val(1)%ptr
          surface_position%val(2)%ptr=tmp_surface_position%val(2)%ptr
          call insert(bc_states(i+1), surface_position, 'Coordinate')
          call deallocate(tmp_surface_position)
          call deallocate(surface_position)

          ! get a velocity field on this surface
          surface_v_mesh=make_mesh(surface_mesh, surface_mesh%shape, &
               continuity=-1, name='VelocityMesh')
          call allocate(tmp_surface_velocity, velocity%dim, surface_v_mesh, &
               "tmp_surface_velocity")
          call remap_field_to_surface(velocity, tmp_surface_velocity, &
               surface_element_list)
          call allocate(surface_velocity, velocity%dim-1, surface_v_mesh, &
               'GeostrophicVelocity')
          surface_velocity%val(1)%ptr=tmp_surface_velocity%val(1)%ptr
          surface_velocity%val(2)%ptr=tmp_surface_velocity%val(2)%ptr
          call insert(bc_states(i+1), surface_velocity, 'GeostrophicVelocity')
          call deallocate(tmp_surface_velocity)
          call deallocate(surface_velocity)
          call deallocate(surface_v_mesh)

          ! get buoyancy field on this surface and insert into state
          buoyancy=extract_surface_field(streamfunction, bc_name, "value")
          buoyancy%name="Buoyancy"
          buoyancy%option_path=trim(bc_path)//"/type::buoyancy/scalar_field"
          call insert(bc_states(i+1), buoyancy, "Buoyancy")

          ! need an OldBuoyancy field to solve advection diffusion equation
          call allocate(OldBuoyancy, surface_mesh, "OldBuoyancy")
          call zero(OldBuoyancy)
          call insert(bc_states(i+1), OldBuoyancy, "OldBuoyancy")
          call deallocate(OldBuoyancy)

          ! get buoyancy sparsity and insert into state
          buoyancy_sparsity(i+1)=make_sparsity(surface_mesh, surface_mesh, &
               "BuoyancySparsity")
          call insert(bc_states(i+1), buoyancy_sparsity(i+1), &
               "BuoyancySparsity")
          call deallocate(buoyancy_sparsity(i+1))

          ! get the associated diffusivity field
          call allocate(diffusivity, surface_mesh, &
               name='BuoyancyDiffusivity')
          diffusivity%option_path=trim(bc_path)//"/type::buoyancy/&
               scalar_field/prognostic/&
               tensor_field::Diffusivity/prescribed/value["//int2str(0)//"]"
          call initialise_field(diffusivity, diffusivity%option_path, &
               surface_position)
          call insert(bc_states, diffusivity, 'BuoyancyDiffusivity')
          call deallocate(diffusivity)

       end if

    end do

    deallocate(buoyancy_sparsity)

  end subroutine construct_2d_state_for_boundary_conditions

  subroutine update_2d_state_for_boundary_conditions(state, bc_states)

    type(state_type), intent(in) :: state
    type(state_type), dimension(:), intent(inout) ::bc_states

    type(scalar_field), pointer :: streamfunction
    type(vector_field), pointer :: velocity, surface_velocity
    type(vector_field) :: tmp_surface_velocity
    integer, pointer, dimension(:) :: surface_element_list

    integer :: i
    character(len=FIELD_NAME_LEN) :: bc_path

    ewrite(1,*) "In update_2d_state_for_boundary_conditions"

    streamfunction=>extract_scalar_field(state, 'Streamfunction')
    velocity=>extract_vector_field(state, 'GeostrophicVelocity')

    ! Loop over boundary conditions
    do i=0,size(bc_states)-1

       bc_path="/material_phase["//int2str(0)//&
            "]/scalar_field::Streamfunction/prognostic/boundary_conditions["&
            //int2str(i)//"]"

       if(have_option(trim(bc_path)//"/type::buoyancy")) then

          ! extract old velocity field
          surface_velocity=>extract_vector_field(bc_states(i+1),&
               'GeostrophicVelocity')

          ! get surface_element_list
          call get_boundary_condition(streamfunction, i+1, &
               surface_element_list=surface_element_list)

          ! update velocity field
          call allocate(tmp_surface_velocity, velocity%dim,&
               surface_velocity%mesh, 'tmp_surface_velocity')
          call remap_field_to_surface(velocity, tmp_surface_velocity,&
               surface_element_list)
          surface_velocity%val(1)%ptr=tmp_surface_velocity%val(1)%ptr
          surface_velocity%val(2)%ptr=tmp_surface_velocity%val(2)%ptr
          call deallocate(tmp_surface_velocity)

       end if

    end do

  end subroutine update_2d_state_for_boundary_conditions

end module qg_forward
