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
  use SPUD
  use sparse_tools
  use sparsity_patterns
  use boundary_conditions
  use advection_diffusion_cg
  use vtk_interfaces
  use global_parameters, only:FIELD_NAME_LEN, global_debug_level
  use PV_INVERSION
  implicit none

CONTAINS

  subroutine qg_forward_run()
    !!< Integrate the forward model for the quasi-geostrophic equations 
    type(state_type), pointer, dimension(:) :: state
    type(state_type), pointer, dimension(:) :: bc_states
    real :: t
    real :: dt
    real :: tmax
    real :: dump_period
    real :: tdump
    real :: adapt_period
    real :: tadapt
    integer :: dump_count=1
    !
    type(scalar_field), pointer :: PV, PV_Old, streamfunction
    type(scalar_field), pointer :: buoyancy, buoyancy_old
    type(vector_field), pointer :: X
    type(scalar_field), dimension(:), allocatable :: PV_RK_Stages
    type(csr_sparsity), pointer :: pv_sparsity
    type(csr_sparsity), pointer :: buoyancy_sparsity
    type(mesh_type), pointer :: mesh
    type(tensor_field) :: metric
    !RK stage coefficients
    real, dimension(:), allocatable :: Rk_a, RK_b
    integer :: nstages, stage
    integer :: i, n_bcs
    integer, dimension(2) :: shape_option

    character(len=FIELD_NAME_LEN) :: bc_name
    character(len=FIELD_NAME_LEN) :: filename=""

    call get_option("/simulation_name",filename)

    ! Get timestepping options
    call get_option("/timestepping/current_time", t)
    call get_option("/timestepping/finish_time", tmax)
    call get_option("/timestepping/timestep", dt)
    call get_option("/io/dump_period", dump_period)

    call populate_state(state)

    ! RK stages
    shape_option = option_shape("/timestepping/rungekutta/stages")
    nstages = shape_option(1)
    allocate( rk_a(nstages), rk_b(nstages), PV_RK_Stages(nstages) )
    call get_option("/timestepping/rungekutta/stages",rk_a)
    call get_option("/timestepping/rungekutta/weights",rk_b)

    call allocate_and_insert_auxilliary_qg_fields(state(1), PV_RK_Stages)

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
       call allocate_and_insert_auxilliary_qg_fields(state(1), PV_RK_Stages)
    end if

    PV => extract_scalar_field(state(1), 'PotentialVorticity')
    PV_Old => extract_scalar_field(state(1), 'OldPotentialVorticity')
    streamfunction => extract_scalar_field(state(1),'Streamfunction')
    X => extract_vector_field(state(1), 'Coordinate')
    
    ! Solve for streamfunction if it is not prescribed
    if(have_option("/material_phase[0]/scalar_field::Streamfunction&
         /prognostic")) then
       ! set up bc_states if prognostic streamfunction has buoyancy
       ! boundary conditions
       n_bcs=option_count("/material_phase["//int2str(0)//&
            "]/scalar_field::Streamfunction/prognostic/boundary_conditions/&
            type::buoyancy")
       if(n_bcs>=1) then
          allocate(bc_states(1:n_bcs))
          do i=1, n_bcs
             call nullify(bc_states(i))
          end do
          call construct_2d_state_for_boundary_conditions(state(1), bc_states)
       end if
       ! Can now solve for streamfunction
       call solve_streamfunction_qg(state(1), PV)
    end if

    ! Always calculate velocity from streamfunction
    call streamfunction2velocity(state(1))

    ! We have now allocated, inserted and calculated all fields so can
    ! write out state at t=0
    call vtk_write_state(trim(filename), 0, state=state)
    if(n_bcs>=1) then
       call vtk_write_state(trim(filename)//"_bcs", 0, state=bc_states)
    end if

    tdump = 0.

    timestep_loop: do
       if(t>tmax) exit

       t = t + dt 
       tdump = tdump + dt
       tadapt = tadapt + dt

       ! Adapt mesh if required
       if(have_option("/mesh_adaptivity").and.tadapt>=adapt_period) then
          ewrite(1,*) "adapting mesh"
          mesh => extract_mesh(state(1), "CoordinateMesh")
          call allocate(metric, mesh, "ErrorMetric")
          call assemble_metric(state, metric)
          call adapt_state(state, metric)
          ! Reallocate and insert stuff into state:
          call allocate_and_insert_auxilliary_qg_fields(state(1), PV_RK_Stages)
          tadapt=0.
          ! Re-extract fields and recalculate velocity from streamfunction
          PV => extract_scalar_field(state(1), 'PotentialVorticity')
          PV_Old => extract_scalar_field(state(1), 'OldPotentialVorticity')
          streamfunction => extract_scalar_field(state(1),'Streamfunction')
          X => extract_vector_field(state(1), 'Coordinate')
          call streamfunction2velocity(state(1))
       end if

       call set(PV_old, PV)

       rkstages: do stage = 1, nstages

          pv_sparsity => extract_csr_sparsity(state(1), &
               'PotentialVorticitySparsity')
          call solve_advection_diffusion_cg('PotentialVorticity', &
          state(1), pv_sparsity, rk_a(stage))
          PV => extract_scalar_field(state(1), "PotentialVorticity")
          call set(PV_RK_Stages(stage), PV)

          ! If streamfunction is prognostic, solve for it and 
          ! recalculate velocity
          if(have_option("/material_phase[0]/scalar_field::&
               Streamfunction/prognostic")) then
             ! If we are advecting buoyancy on the top and bottom
             ! surfaces, we must solve for it first:
             call update_2d_state_for_boundary_conditions(state(1), bc_states)
             do i=0, n_bcs-1
                buoyancy=>extract_scalar_field(bc_states(i+1), "Buoyancy")
                buoyancy_old=>extract_scalar_field(bc_states(i+1), &
                     "OldBuoyancy")
                call set(buoyancy_old, buoyancy)
                buoyancy_sparsity=>extract_csr_sparsity(bc_states(i+1), &
                     "BuoyancySparsity")
                call set(buoyancy_old, buoyancy)
                call solve_advection_diffusion_cg('Buoyancy', bc_states(i+1), &
                  buoyancy_sparsity)
                call get_option("/material_phase["//int2str(0)//&
                     "]/scalar_field::Streamfunction/prognostic/&
                     boundary_conditions["//int2str(i)//"]/name",&
                     bc_name)
                call insert_surface_field(streamfunction, trim(bc_name), &
                     buoyancy)
             end do
             call solve_streamfunction_qg(state(1), PV)
             call streamfunction2velocity(state(1))
          end if

       end do rkstages

       call zero(PV)
       rkstage_step: do stage = 1, nstages
          call addto(PV, PV_RK_Stages(stage), rk_b(stage))
       end do rkstage_step

       if(tdump>=dump_period) then 
          ewrite(1,*) 'dump', dump_count
          tdump = 0.
          
          call vtk_write_state(trim(filename), dump_count, state=state)
          if(n_bcs>=1) then
             call vtk_write_state(trim(filename)//"_bcs", dump_count, state=bc_states)
          end if
          dump_count = dump_count + 1
       end if

    end do timestep_loop

    do i=1, size(state)
       call deallocate(state(i))
    end do

    do i=1, size(bc_states)
       call deallocate(bc_states(i))
    end do

    ewrite(0,*) "Printing references at the end of Qg_forward "
    ewrite(0,*) "- there shouldn't be any..."
    call print_references(0)

  end subroutine qg_forward_run

  subroutine allocate_and_insert_auxilliary_qg_fields(state, PV_RK_Stages)
    !!< This subroutine allocates and inserts into state 
    !!< fields that are unique to qg_strat

    type(state_type), intent(inout) :: state

    type(scalar_field), pointer :: PV
    type(scalar_field) :: OldPV
    type(scalar_field), dimension(:) :: PV_RK_Stages
    type(vector_field), pointer :: X
    type(vector_field) :: V
    type(csr_sparsity) :: pv_sparsity
    type(mesh_type) :: V_mesh

    integer :: stage

    ewrite(1,*) "In allocate_and_insert_auxilliary_qg_fields"

    PV => extract_scalar_field(state, 'PotentialVorticity')
    X => extract_vector_field(state, 'Coordinate')

    ! allocate, zero and insert OldPotentialVorticity
    call allocate_scalar_field(OldPV, PV%mesh, 'OldPotentialVorticity')
    call zero(OldPV)
    call insert(state, OldPV, 'OldPotentialVorticity')
    call deallocate(OldPV)

    ! construct mesh for velocity field (it's DG because it is the skew
    ! gradient of a streamfunction)
    V_mesh=make_mesh(X%mesh, PV%mesh%shape, continuity=-1, &
         name='VelocityMesh')
    call allocate_vector_field(V, mesh_dim(PV), V_mesh, 'NonlinearVelocity')
    call zero(V)
    call insert(state, V, 'NonlinearVelocity')
    call deallocate(V_mesh)
    call deallocate(V)

    pv_sparsity=make_sparsity(PV%mesh, PV%mesh,  &
         "PotentialVorticitySparsity")
    call insert(state, pv_sparsity, "PotentialVorticitySparsity")
    call deallocate(pv_sparsity)    

    do stage = 1, size(PV_RK_Stages)
       call allocate_scalar_field( PV_RK_Stages(stage), PV%mesh, &
            'PV_RK_Stages'//int2str(stage) )
       call zero(PV_RK_Stages(stage))
       call insert(state, PV_RK_Stages(stage), 'PV_RK_Stages'//int2str(stage))
       call deallocate(PV_RK_Stages(stage))
    end do

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
    velocity=>extract_vector_field(state, 'NonlinearVelocity')

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
               'NonlinearVelocity')
          surface_velocity%val(1)%ptr=tmp_surface_velocity%val(1)%ptr
          surface_velocity%val(2)%ptr=tmp_surface_velocity%val(2)%ptr
          call insert(bc_states(i+1), surface_velocity, 'NonlinearVelocity')
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
    velocity=>extract_vector_field(state, 'NonlinearVelocity')

    ! Loop over boundary conditions
    do i=0,size(bc_states)-1

       bc_path="/material_phase["//int2str(0)//&
            "]/scalar_field::Streamfunction/prognostic/boundary_conditions["&
            //int2str(i)//"]"

       if(have_option(trim(bc_path)//"/type::buoyancy")) then

          ! extract old velocity field
          surface_velocity=>extract_vector_field(bc_states(i+1),&
               'NonlinearVelocity')

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
