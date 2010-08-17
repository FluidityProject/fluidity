!    Copyright (C) 2007 Imperial College London and others.
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
module boundary_conditions_from_options
use fldebug
use quadrature
use elements
use fields
use field_options
use state_module
use sparse_tools_petsc
use boundary_conditions
use initialise_fields_module
use global_parameters, only: OPTION_PATH_LEN, PYTHON_FUNC_LEN, pi, current_debug_level
use tidal_module
use SampleNetCDF
use coordinates
use synthetic_bc
use spud
use vector_tools
use vtk_interfaces
use pickers_inquire

implicit none

  private
  public populate_boundary_conditions, set_boundary_conditions_values, &
       apply_dirichlet_conditions_inverse_mass, impose_reference_pressure_node, &
       impose_reference_velocity_node
  public :: populate_scalar_boundary_conditions, &
    & populate_vector_boundary_conditions

  interface apply_dirichlet_conditions_inverse_mass
     module procedure apply_dirichlet_conditions_inverse_mass_vector, &
                      apply_dirichlet_conditions_inverse_mass_vector_lumped
  end interface apply_dirichlet_conditions_inverse_mass

contains

  subroutine populate_boundary_conditions(states)
    ! Populate the boundary conditions of all fields
    ! This is called as part of populate_state but also
    ! after an adapt.
    type(state_type), dimension(:), intent(in):: states

    ! these must be pointers as bc's should be added to the original field
    type(scalar_field), pointer:: sfield
    type(vector_field), pointer:: vfield

    type(vector_field), pointer:: position
    character(len=OPTION_PATH_LEN) phase_path, field_path
    integer, dimension(:), allocatable:: surface_ids
    integer shape_option(2)
    integer p, f, nphases, nfields

    ewrite(1,*) "In populate_boundary_conditions"

    nphases = option_count('/material_phase')  
    do p = 0, nphases-1

       phase_path = '/material_phase['//int2str(p)//']'

       position => extract_vector_field(states(p+1), "Coordinate")

       ! Scalar fields:

       nfields = scalar_field_count(states(p+1))
       do f = 1, nfields
          sfield => extract_scalar_field(states(p+1),f)
          field_path=sfield%option_path

          if (.not. have_option(trim(field_path)//'/prognostic')) cycle

          ! only prognostic fields from here:
          call populate_scalar_boundary_conditions(sfield, &
               trim(field_path)//'/prognostic/boundary_conditions', position)

       end do

       ! Vector fields:

       nfields = vector_field_count(states(p+1))
       do f = 1, nfields
          vfield => extract_vector_field(states(p+1), f)
          field_path=vfield%option_path
          
          if (.not. have_option(trim(field_path)//'/prognostic')) cycle

          ! only prognostic fields from here:
          call populate_vector_boundary_conditions(states(p+1),vfield, &
               trim(field_path)//'/prognostic/boundary_conditions', position)

       end do

    end do

    ! special case 'boundary conditions', include:
    ! - ocean boundaries
    ! - ocean forcing
    ! - GLS stable boundaries
    if (have_option('/geometry/ocean_boundaries')) then
       ! NOTE: has to be a pointer, as bcs should be added to original field
       sfield => extract_scalar_field(states(1), "DistanceToTop")
       ! Get vector of surface ids
       shape_option=option_shape('/geometry/ocean_boundaries/top_surface_ids')
       allocate(surface_ids(1:shape_option(1)))
       call get_option('/geometry/ocean_boundaries/top_surface_ids', surface_ids)
       ! Add boundary condition that marks the top of the domain
       call add_boundary_condition(sfield, "top", "surface", surface_ids)
       deallocate(surface_ids)

       ! NOTE: has to be a pointer, as bcs should be added to original field
       sfield => extract_scalar_field(states(1), "DistanceToBottom")
       ! Get vector of surface ids
       shape_option=option_shape('/geometry/ocean_boundaries/bottom_surface_ids')
       allocate(surface_ids(1:shape_option(1)))
       call get_option('/geometry/ocean_boundaries/bottom_surface_ids', surface_ids)
       ! Add boundary condition that marks the top of the domain
       call add_boundary_condition(sfield, "bottom", "surface", surface_ids)
       deallocate(surface_ids)
       
       call ocean_boundaries_stats(states(1))

       if (have_option('/ocean_forcing')) then
           call populate_ocean_forcing_boundary_conditions(states(1))
       end if

       if (have_option('/material_phase[0]/subgridscale_parameterisations/GLS/calculate_boundaries')) then
           call populate_gls_boundary_conditions(states(1))
       end if

    end if
    if (have_option('/turbine_model')) then
       call populate_flux_turbine_boundary_conditions(states(1))
    end if

  end subroutine populate_boundary_conditions

  subroutine populate_scalar_boundary_conditions(field, bc_path, position)
    ! Populate the boundary conditions of one scalar field
    ! needs to be a pointer:
    type(scalar_field), pointer:: field
    character(len=*), intent(in):: bc_path
    type(vector_field), intent(in):: position

    type(mesh_type), pointer:: surface_mesh
    type(scalar_field) surface_field
    type(vector_field) bc_position
    integer, dimension(:), pointer:: surface_element_list
    character(len=OPTION_PATH_LEN) bc_path_i
    character(len=FIELD_NAME_LEN) bc_name, bc_type
    integer, dimension(:), allocatable:: surface_ids
    integer i, nbcs, shape_option(2)

    ! Get number of boundary conditions
    nbcs=option_count(trim(bc_path))

    ! Loop over boundary conditions
    boundary_conditions: do i=0, nbcs-1

       bc_path_i=trim(bc_path)//"["//int2str(i)//"]"

       ! Get vector of surface ids
       shape_option=option_shape(trim(bc_path_i)//"/surface_ids")
       allocate(surface_ids(1:shape_option(1)))
       call get_option(trim(bc_path_i)//"/surface_ids", surface_ids)

       ! Get name of boundary
       call get_option(trim(bc_path_i)//"/name", bc_name)

       ! Get type of boundary condition
       call get_option(trim(bc_path_i)//"/type[0]/name", bc_type)

       if(have_option(trim(bc_path_i)//"/type[0]/apply_weakly")) then
         bc_type = "weak"//trim(bc_type)
       end if

       ! Add boundary condition
       call add_boundary_condition(field, trim(bc_name), trim(bc_type), &
            surface_ids, option_path=bc_path_i)
       deallocate(surface_ids)

       ! mesh of only the part of the surface where this b.c. applies
       call get_boundary_condition(field, i+1, surface_mesh=surface_mesh, &
          surface_element_list=surface_element_list)
       call allocate(bc_position, position%dim, surface_mesh)
       call remap_field_to_surface(position, bc_position, surface_element_list)

       ! Dirichlet and Neumann boundary conditions require one input
       ! while a Robin boundary condition requires two. This input can
       ! be constant or set from a generic or python function.
       select case(trim(bc_type))

       case("dirichlet", "neumann", "weakdirichlet", &
            "buoyancy")

          call allocate(surface_field, surface_mesh, name="value")
          call insert_surface_field(field, i+1, surface_field)
          call deallocate(surface_field)

       case("robin")

          call allocate(surface_field, surface_mesh, name="order_zero_coefficient")
          call insert_surface_field(field, i+1, surface_field)
          call deallocate(surface_field)

          call allocate(surface_field, surface_mesh, name="order_one_coefficient")
          call insert_surface_field(field, i+1, surface_field)
          call deallocate(surface_field)
          
       case("zero_flux")
       
          ! nothing to be done here
       
       case( "k_epsilon" )
       
           if(.not. have_option &
           ("/material_phase[0]/subgridscale_parameterisations/k-epsilon/") ) then
               FLAbort("Incorrect boundary condition type for field")
           end if

       case default

          ! This really shouldn't happen
          FLAbort("Incorrect boundary condition type for field")

       end select
       
       call deallocate(bc_position)

    end do boundary_conditions

  end subroutine populate_scalar_boundary_conditions

  subroutine populate_vector_boundary_conditions(state, field, bc_path, position)
    ! Populate the boundary conditions of one vector field
    ! needs to be a pointer:
    type(state_type), intent(in) :: state
    type(vector_field), pointer:: field
    character(len=*), intent(in):: bc_path
    type(vector_field), intent(in):: position

    ! possible vector components for vector b.c.s
    ! either carteisan aligned or aligned with the surface
    character(len=20), parameter, dimension(3) :: &
         cartesian_aligned_components=(/ &
           "x_component", &
           "y_component", &
           "z_component" /), &
         surface_aligned_components=(/ &
           "normal_component   ", & 
           "tangent_component_1", &
           "tangent_component_2" /)

    character(len=20), dimension(3) :: aligned_components

    type(mesh_type), pointer:: mesh, surface_mesh
    type(vector_field) surface_field, surface_field2, bc_position
    type(vector_field):: normal, tangent_1, tangent_2
    type(scalar_field) :: scalar_surface_field, scalar_surface_field2
    character(len=OPTION_PATH_LEN) bc_path_i, bc_type_path, bc_component_path
    character(len=FIELD_NAME_LEN) bc_name, bc_type
    logical applies(3), have_sem_bc, debugging_mode
    integer, dimension(:), allocatable:: surface_ids
    integer, dimension(:), pointer:: surface_element_list, surface_node_list
    integer i, j, nbcs, shape_option(2)

    nbcs=option_count(trim(bc_path))

    boundary_conditions: do i=0, nbcs-1
       bc_path_i=trim(bc_path)//"["//int2str(i)//"]"
       shape_option=option_shape(trim(bc_path_i)//"/surface_ids")
       allocate(surface_ids(1:shape_option(1)))
       call get_option(trim(bc_path_i)//"/surface_ids", surface_ids)
       call get_option(trim(bc_path_i)//"/name", bc_name)

       call get_option(trim(bc_path_i)//"/type[0]/name", bc_type)

       if(have_option(trim(bc_path_i)//"/type[0]/apply_weakly")) then
         bc_type = "weak"//trim(bc_type)
       end if

       select case(trim(bc_type))
       case("dirichlet", "neumann", "weakdirichlet")

          if(have_option(trim(bc_path_i)//"/type[0]/align_bc_with_cartesian")) then
             aligned_components=cartesian_aligned_components
             bc_type_path=trim(bc_path_i)//"/type[0]/align_bc_with_cartesian"
          else
             aligned_components=surface_aligned_components             
             bc_type_path=trim(bc_path_i)//"/type[0]/align_bc_with_surface"
          end if

          have_sem_bc=.false.
          do j=1,3
             bc_component_path=trim(bc_type_path)//"/"//aligned_components(j)
             applies(j)=have_option(trim(bc_component_path))
             ! check for SEM bc:
             bc_component_path=trim(bc_component_path)//'/synthetic_eddy_method'
             have_sem_bc=have_sem_bc .or. (applies(j) .and. have_option(bc_component_path))
          end do
          call add_sem_bc(have_sem_bc)
          
          call add_boundary_condition(field, trim(bc_name), trim(bc_type),&
               & surface_ids, applies=applies, option_path=bc_path_i)
          deallocate(surface_ids)
          
          call get_boundary_condition(field, i+1, surface_mesh=surface_mesh)
          call allocate(surface_field, field%dim, surface_mesh, name="value")
          
          call insert_surface_field(field, i+1, surface_field)
          call deallocate(surface_field)
          
          if (have_sem_bc) then
             call allocate(surface_field, field%dim, surface_mesh, name="TurbulenceLengthscale")
             call insert_surface_field(field, i+1, surface_field)
             call deallocate(surface_field)

             call allocate(surface_field, field%dim, surface_mesh, name="MeanProfile")
             call insert_surface_field(field, i+1, surface_field)
             call deallocate(surface_field)

             call allocate(surface_field, field%dim, surface_mesh, name="ReStressesProfile")
             call insert_surface_field(field, i+1, surface_field)
             call deallocate(surface_field)
          end if

       case("robin")

          if(have_option(trim(bc_path_i)//"/type[0]/align_bc_with_cartesian")) then
             aligned_components=cartesian_aligned_components
             bc_type_path=trim(bc_path_i)//"/type[0]/align_bc_with_cartesian"
          else
             aligned_components=surface_aligned_components             
             bc_type_path=trim(bc_path_i)//"/type[0]/align_bc_with_surface"
          end if

          do j=1,3
             bc_component_path=trim(bc_type_path)//"/"//aligned_components(j)
             applies(j)=have_option(trim(bc_component_path))
          end do

          call add_boundary_condition(field, trim(bc_name), trim(bc_type),&
               & surface_ids, applies=applies, option_path=bc_path_i)
          deallocate(surface_ids)

          call get_boundary_condition(field, i+1, surface_mesh=surface_mesh)

          call allocate(surface_field, field%dim, surface_mesh, name="order_zero_coefficient")
          call allocate(surface_field2, field%dim, surface_mesh, name="order_one_coefficient")
          call insert_surface_field(field, i+1, surface_field)
          call insert_surface_field(field, i+1, surface_field2)
          call deallocate(surface_field)
          call deallocate(surface_field2)

       case("drag")

          call add_boundary_condition(field, trim(bc_name), trim(bc_type), &
               & surface_ids, option_path=bc_path_i)
          deallocate(surface_ids)

          call get_boundary_condition(field, i+1, surface_mesh=surface_mesh)
          call allocate(scalar_surface_field, surface_mesh, name="DragCoefficient")
          call insert_surface_field(field, i+1, scalar_surface_field)
          call deallocate(scalar_surface_field)

       case ("wind_forcing")

          call add_boundary_condition(field, trim(bc_name), trim(bc_type), &
               & surface_ids, option_path=bc_path_i)
          deallocate(surface_ids)

          call get_boundary_condition(field, i+1, surface_mesh=surface_mesh)
          call allocate(surface_field, field%dim-1, surface_mesh, name="WindSurfaceField")
          call insert_surface_field(field, i+1, surface_field)
          call deallocate(surface_field)
          bc_path_i=trim(bc_path_i)//"/type[0]/wind_velocity"
          if (have_option(bc_path_i)) then
             call allocate(scalar_surface_field, surface_mesh, name="WindDragCoefficient")
             call insert_surface_field(field, i+1, scalar_surface_field)
             call deallocate(scalar_surface_field)
          end if
          
       case ("free_surface", "no_normal_flow")

          ! these are marked as applying in the 1st direction only
          ! so they could potentially be combined with rotated bcs
          ! applying in the tangential directions only
          call add_boundary_condition(field, trim(bc_name), trim(bc_type), &
               & surface_ids, option_path=bc_path_i, &
               & applies=(/ .true., .false., .false. /) )
          deallocate(surface_ids)

          if (trim(bc_type)=="free_surface") then
             bc_path_i=trim(bc_path_i)//"/type[0]/wetting_drying"
             if(have_option(trim(bc_path_i))) then
                  ! Wetting and drying needs an auxiliary field on the pressure mesh
                  mesh => extract_pressure_mesh(state)
                  call get_boundary_condition(field, i+1, surface_element_list=surface_element_list)
                  allocate(surface_mesh)
                  call create_surface_mesh(surface_mesh, surface_node_list, mesh, surface_element_list, "PressureSurfaceMesh")
                  call allocate(scalar_surface_field, surface_mesh, name="WettingDryingAlpha")
                  call allocate(scalar_surface_field2, surface_mesh, name="WettingDryingOldAlpha")
                  call insert_surface_field(field, i+1, scalar_surface_field)
                  call insert_surface_field(field, i+1, scalar_surface_field2)
                  call deallocate(scalar_surface_field)
                  call deallocate(scalar_surface_field2)
                  deallocate(surface_mesh)
             end if
          end if

       case ("near_wall_treatment", "log_law_of_wall")
         
          ! these are marked as applying in the 2nd (and 3d if present) only
          ! so they could potentially be combined with a no_normal_flow
          ! or a rotatted bc applying in the normal direction only
          call add_boundary_condition(field, trim(bc_name), trim(bc_type), &
               & surface_ids, option_path=bc_path_i, &
               & applies=(/ .false., .true., .true. /) )
          deallocate(surface_ids)

       case ("outflow")

          ! dummy bc for outflow planes
          call add_boundary_condition(field, trim(bc_name), trim(bc_type), &
               & surface_ids, option_path=bc_path_i, &
               & applies=(/ .true., .true., .true. /) )
          deallocate(surface_ids) 
 
       case default
          FLAbort("Incorrect boundary condition type for field")
       end select

       ! now check for user-specified normal/tangent vectors (rotation matrix)
       select case (bc_type)
       case ("dirichlet", "neumann", "robin", "weakdirichlet", "near_wall_treatment", "log_law_of_wall")
          ! this is the same for all 3 b.c. types

          bc_type_path=trim(bc_path_i)//"/type[0]/align_bc_with_surface"

          ! map the coordinate field onto this mesh
          call get_boundary_condition(field, i+1, surface_mesh=surface_mesh, &
               surface_element_list=surface_element_list)
          call allocate(bc_position, position%dim, surface_mesh)
          call remap_field_to_surface(position, bc_position, surface_element_list)

          if (have_option(bc_type_path) .or. bc_type=="near_wall_treatment" &
              .or. bc_type=="log_law_of_wall") then

             call allocate(normal, field%dim, surface_mesh, name="normal")
             bc_component_path=trim(bc_type_path)//"/normal_direction"
             if (have_option(bc_component_path)) then
                call initialise_field(normal, bc_component_path, bc_position)
             else
                call zero(normal)
             end if
             call insert_surface_field(field, i+1, normal)

             call allocate(tangent_1, field%dim, surface_mesh, name="tangent1")
             bc_component_path=trim(bc_type_path)//"/tangent_direction_1"
             if (have_option(bc_component_path)) then
                call initialise_field(tangent_1, bc_component_path, bc_position)
             else
                call zero(tangent_1)
             end if
             call insert_surface_field(field, i+1, tangent_1)

             call allocate(tangent_2, field%dim, surface_mesh, name="tangent2")
             bc_component_path=trim(bc_type_path)//"/tangent_direction_2"
             if (have_option(bc_component_path)) then
                call initialise_field(tangent_2, bc_component_path, bc_position)
             else
                call zero(tangent_2)
             end if
             call insert_surface_field(field, i+1, tangent_2)

             debugging_mode=have_option(trim(bc_type_path)//"/debugging_mode")
            
             ! calculate the normal, tangent_1 and tangent_2 on every boundary node
             call initialise_rotated_bcs(surface_element_list, &
                  position, debugging_mode, normal, tangent_1, tangent_2)
             call deallocate(normal)
             call deallocate(tangent_1)
             call deallocate(tangent_2)
             
          end if
          call deallocate(bc_position)

       case default
          ! nothing to do for other bcs
       end select

    end do boundary_conditions

  end subroutine populate_vector_boundary_conditions

  subroutine set_boundary_conditions_values(states, shift_time)
    !!< Set the values of the boundary conditions of all fields
    !!< This is called each time step.
    type(state_type), dimension(:), intent(in):: states
    !! if present and true the time level at which the bcs are evaluated
    !! is shifted according to:
    !! "dirichlet": current_time+dt
    !! all others:  current_time+theta*dt
    !! Otherwise (no shift_time) current_time is used, which is how this
    !! routine should be called for initialisation, so that for instance
    !! the fields can be overwritten with the right initial bc values.
    logical, optional, intent(in):: shift_time

    type(scalar_field), pointer:: sfield
    type(vector_field), pointer:: vfield

    type(vector_field), pointer:: position
    character(len=OPTION_PATH_LEN) phase_path, field_path
    integer p, f, nphases, nfields

    ewrite(1,*) "In set_boundary_conditions"

    nphases = size(states)
    do p = 0, nphases-1

       phase_path = '/material_phase['//int2str(p)//']'

       position => extract_vector_field(states(p+1), "Coordinate")

       ! Scalar fields:

       nfields = scalar_field_count(states(p+1))
       do f = 1, nfields
          sfield => extract_scalar_field(states(p+1),f)
          field_path=sfield%option_path
          if (.not. have_option(trim(field_path)//'/prognostic')) cycle

          ! only prognostic fields from here:
          call set_scalar_boundary_conditions_values(states(p+1), sfield, &
               trim(field_path)//'/prognostic/boundary_conditions', &
               position, shift_time=shift_time)

       end do

       ! Vector fields:

       nfields = vector_field_count(states(p+1))
       do f = 1, nfields
          vfield => extract_vector_field(states(p+1), f)
          field_path=vfield%option_path
          if (.not. have_option(trim(field_path)//'/prognostic')) cycle

          ! only prognostic fields from here:
          call set_vector_boundary_conditions_values(vfield, &
               trim(field_path)//'/prognostic/boundary_conditions', &
               position, shift_time=shift_time)

       end do

    end do

    if (have_option('/ocean_forcing')) then
        call set_ocean_forcings_boundary_conditions(states(1))
    end if

    if (have_option('/turbine_model')) then
        call set_dirichlet_turbine_boundary_conditions(states(1))
        call set_flux_turbine_boundary_conditions(states(1))
    end if

  end subroutine set_boundary_conditions_values

  subroutine set_scalar_boundary_conditions_values(state, field, bc_path, position, shift_time)
    ! Set the boundary condition values of one scalar field
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout):: field
    character(len=*), intent(in):: bc_path
    type(vector_field), intent(in):: position
    ! see above in set_boundary_conditions:
    logical, optional, intent(in):: shift_time

    type(mesh_type), pointer:: surface_mesh
    type(scalar_field), pointer:: surface_field
    type(vector_field) bc_position, temp_position
    character(len=OPTION_PATH_LEN) bc_path_i, bc_type_path
    character(len=FIELD_NAME_LEN) bc_name, bc_type
    real:: time, theta, dt
    integer, dimension(:), pointer:: surface_element_list
    integer i, nbcs

    ! Get number of boundary conditions
    nbcs=option_count(trim(bc_path))

    ! Loop over boundary conditions
    boundary_conditions: do i=0, nbcs-1

       bc_path_i=trim(bc_path)//"["//int2str(i)//"]"

       ! Get name of boundary (or set a default if none present)
       call get_option(trim(bc_path_i)//"/name", bc_name)

       ! Get type of boundary condition
       call get_option(trim(bc_path_i)//"/type[0]/name", bc_type)

       if(have_option(trim(bc_path_i)//"/type[0]/apply_weakly")) then
         bc_type = "weak"//trim(bc_type)
       end if

       ! mesh of only the part of the surface where this b.c. applies
       call get_boundary_condition(field, i+1, surface_mesh=surface_mesh, &
            surface_element_list=surface_element_list)
       ! map the coordinate field onto this mesh
       call allocate(bc_position, position%dim, surface_mesh, "BCPositions")

       if((surface_mesh%shape%degree==0).and.(bc_type=="dirichlet")) then
         ! if the boundary condition is on a 0th degree mesh and is of type strong dirichlet
         ! then the positions used to calculate the bc should be body element centred not
         ! surface element centred
         call allocate(temp_position, position%dim, field%mesh, "TemporaryPositions")
         ! first remap to body element centred positions
         call remap_field(position, temp_position)
         ! then remap these to the surface
         call remap_field_to_surface(temp_position, bc_position, surface_element_list)
         call deallocate(temp_position)

       else
         ! in all other cases the positions are remapped to the actual surface
         call remap_field_to_surface(position, bc_position, surface_element_list)
       end if

       ! Dirichlet and Neumann boundary conditions require one input
       ! while a Robin boundary condition requires two. This input can
       ! be constant or set from a generic or python function.
       select case(trim(bc_type))

       case("dirichlet", "neumann", "weakdirichlet")

          bc_type_path=trim(bc_path_i)//"/type[0]"

          surface_field => extract_surface_field(field, bc_name, "value")
          
          ! work out time level at which to evaluate:
          call get_option("/timestepping/current_time", time)
          if (present_and_true(shift_time)) then
            call get_option("/timestepping/timestep", dt)
            if (bc_type=="dirichlet") then
              time=time+dt
            else
              call get_option( trim(field%option_path)// &
                "/prognostic/temporal_discretisation/theta", theta, default=0.5)
              time=time+theta*dt
            end if
          end if

          ! Tidal: set free surface height at the boundary.
          if (have_option(trim(bc_type_path)//"/from_file")) then
            ! Special case for tidal harmonic boundary conditions
            call set_tidal_bc_value(surface_field, bc_position, trim(bc_type_path), field%name)
          else if (have_option(trim(bc_type_path)//"/NEMO_data")) then
            call set_nemo_bc_value(state, surface_field, bc_position, trim(bc_type_path), field%name, surface_element_list)
          else
            call initialise_field(surface_field, bc_type_path, bc_position, &
              time=time)
          end if

       case("robin")

          bc_type_path=trim(bc_path_i)//"/type[0]/order_zero_coefficient"
          surface_field => extract_surface_field(field, bc_name, name="order_zero_coefficient")
          call initialise_field(surface_field, bc_type_path, bc_position)

          bc_type_path=trim(bc_path_i)//"/type[0]/order_one_coefficient"
          surface_field => extract_surface_field(field, bc_name, name="order_one_coefficient")
          call initialise_field(surface_field, bc_type_path, bc_position)
          
       case( "buoyancy")

          bc_type_path=trim(bc_path_i)//"/type::buoyancy/scalar_field/prognostic/initial_condition"
          surface_field => extract_surface_field(field, bc_name, "value")
          call initialise_field(surface_field, bc_type_path, bc_position)
          
       case( "zero_flux" )
       
          ! nothing to be done here

       case( "k_epsilon" )
       
           if(.not. have_option &
           ("/material_phase[0]/subgridscale_parameterisations/k-epsilon/") ) then
               FLAbort("Incorrect boundary condition type for field")
           end if

       case default

          ! This really shouldn't happen
          FLAbort("Incorrect boundary condition type for field")

       end select

       call deallocate(bc_position)

    end do boundary_conditions

  end subroutine set_scalar_boundary_conditions_values

  subroutine set_vector_boundary_conditions_values(field, bc_path, position, &
    shift_time)
    ! Set the boundary condition values of one vector field
    type(vector_field), intent(inout):: field
    character(len=*), intent(in):: bc_path
    type(vector_field), intent(in):: position
    ! see above in set_boundary_conditions:
    logical, optional, intent(in):: shift_time

    ! possible vector components for vector b.c.s
    ! either cartesian aligned or aligned with the surface
    character(len=20), parameter, dimension(3) :: &
         cartesian_aligned_components=(/ &
           "x_component", &
           "y_component", &
           "z_component" /), &
         surface_aligned_components=(/ &
           "normal_component   ", &
           "tangent_component_1", &
           "tangent_component_2" /)

    character(len=20), dimension(3) :: aligned_components

    ! for sem
    logical have_sem_bc
    integer ns, nots

    type(mesh_type), pointer:: surface_mesh
    type(scalar_field) surface_field_component
    type(scalar_field), pointer:: scalar_surface_field
    type(vector_field), pointer:: surface_field, surface_field11
    type(vector_field), pointer:: surface_field2, surface_field21, surface_field22
    type(vector_field) bc_position, temp_position
    character(len=OPTION_PATH_LEN) bc_path_i, bc_type_path, bc_component_path
    character(len=FIELD_NAME_LEN) bc_name, bc_type
    logical applies(3)
    real:: time, theta, dt
    integer, dimension(:), pointer:: surface_element_list
    integer i, j, k, nbcs

    ns=1
    nbcs=option_count(trim(bc_path))
   
    boundary_conditions: do i=0, nbcs-1
       bc_path_i=trim(bc_path)//"["//int2str(i)//"]"
       call get_option(trim(bc_path_i)//"/name", bc_name)

       bc_path_i=trim(bc_path_i)//'/type[0]'
       call get_option(trim(bc_path_i)//"/name", bc_type)

       if(have_option(trim(bc_path_i)//"/apply_weakly")) then
         bc_type = "weak"//trim(bc_type)
       end if
       
       ! work out time level at which to evaluate:
       call get_option("/timestepping/current_time", time)
       if (present_and_true(shift_time)) then
         call get_option("/timestepping/timestep", dt)
         if (bc_type=="dirichlet") then
           time=time+dt
         else
           call get_option( trim(field%option_path)// &
             "/prognostic/temporal_discretisation/theta", theta, default=0.5)
           time=time+theta*dt
         end if
       end if

       select case(trim(bc_type))
       case("dirichlet", "neumann", "weakdirichlet")

          if(have_option(trim(bc_path_i)//"/align_bc_with_cartesian")) then
             aligned_components=cartesian_aligned_components
             bc_type_path=trim(bc_path_i)//"/align_bc_with_cartesian"
          else
             aligned_components=surface_aligned_components             
             bc_type_path=trim(bc_path_i)//"/align_bc_with_surface"
          end if

          have_sem_bc=.false.
          do j=1,3
             bc_component_path=trim(bc_type_path)//"/"//aligned_components(j)
             applies(j)=have_option(trim(bc_component_path))
             ! check for SEM bc:
             bc_component_path=trim(bc_component_path)//'/synthetic_eddy_method'
             have_sem_bc=have_sem_bc .or. (applies(j) .and. have_option(bc_component_path))
          end do

          call get_boundary_condition(field, i+1, surface_mesh=surface_mesh, &
               surface_element_list=surface_element_list)
          surface_field => extract_surface_field(field, bc_name, name="value")
          ! map the coordinate field onto this mesh
          call allocate(bc_position, position%dim, surface_mesh)
          
          if((surface_mesh%shape%degree==0).and.(bc_type=="dirichlet")) then
            ! if the boundary condition is on a 0th degree mesh and is of type strong dirichlet
            ! then the positions used to calculate the bc should be body element centred not
            ! surface element centred
            call allocate(temp_position, position%dim, field%mesh, "TemporaryPositions")
            ! first remap to body element centred positions
            call remap_field(position, temp_position)
            ! then remap these to the surface
            call remap_field_to_surface(temp_position, bc_position, surface_element_list)
            call deallocate(temp_position)
    
         else
            ! in all other cases the positions are remapped to the actual surface
            call remap_field_to_surface(position, bc_position, surface_element_list)
         end if
         
         if (have_sem_bc) then
            surface_field11 => extract_surface_field(field, bc_name, name="MeanProfile")
            do k=1,3
               surface_field_component=extract_scalar_field(surface_field11, k)
               bc_component_path=trim(bc_type_path)//"/"//trim(aligned_components(k))//"/synthetic_eddy_method/mean_profile"
               call initialise_field(surface_field_component, trim(bc_component_path), bc_position, &
                  time=time)
            enddo
            
            surface_field22 => extract_surface_field(field, bc_name, name="TurbulenceLengthscale")
            do k=1,3
               surface_field_component=extract_scalar_field(surface_field22, k)
               bc_component_path=trim(bc_type_path)//"/"//trim(aligned_components(k))//"/synthetic_eddy_method/turbulence_lengthscale"
               call initialise_field(surface_field_component, trim(bc_component_path), bc_position, &
                  time=time)
            enddo
            
            surface_field21 => extract_surface_field(field, bc_name, name="ReStressesProfile")
            do k=1,3
               surface_field_component=extract_scalar_field(surface_field21, k)
               bc_component_path=trim(bc_type_path)//"/"//trim(aligned_components(k))//"/synthetic_eddy_method/Re_stresses_profile"
               call initialise_field(surface_field_component, trim(bc_component_path), bc_position, &
                 time=time)
            enddo
           
            bc_component_path=trim(bc_type_path)//"/"//trim(aligned_components(1))//"/synthetic_eddy_method/number_of_eddies"
            call get_option(trim(bc_component_path),nots)

            ! allocate memory for eddies...
            call initialise_sem_memory(ns,nots)
            ! calculate the boundary condition...
            call synthetic_eddy_method(surface_field, surface_field11, surface_field21, surface_field22,  &
                 bc_position, bc_component_path, ns)
            ns=ns+1
         
         else
            
            do j=1,3
               if (applies(j)) then
                  if (j>surface_field%dim) then
                     FLExit("Too many dimensions in boundary condition")
                  end if
                  
                  bc_component_path=trim(bc_type_path)//"/"//aligned_components(j)
                  surface_field_component=extract_scalar_field(surface_field, j)
                  
                  call initialise_field(surface_field_component, bc_component_path, bc_position, &
                       time=time)
                  
               end if
            end do
         end if

         call deallocate(bc_position)
         
      case("robin")

          if(have_option(trim(bc_path_i)//"/align_bc_with_cartesian")) then
             aligned_components=cartesian_aligned_components
             bc_type_path=trim(bc_path_i)//"/align_bc_with_cartesian"
          else
             aligned_components=surface_aligned_components             
             bc_type_path=trim(bc_path_i)//"/align_bc_with_surface"
          end if

          do j=1,3
             bc_component_path=trim(bc_type_path)//"/"//aligned_components(j)
             applies(j)=have_option(trim(bc_component_path))
          end do

          call get_boundary_condition(field, i+1, surface_mesh=surface_mesh, &
               surface_element_list=surface_element_list)
          ! map the coordinate field onto this mesh
          call allocate(bc_position, position%dim, surface_mesh)
          call remap_field_to_surface(position, bc_position, surface_element_list)

          surface_field => extract_surface_field(field, bc_name, name="order_zero_coeffcient")
          surface_field2 => extract_surface_field(field, bc_name, name="order_one_coeffcient")
          do j=1,3
             if (j>surface_field%dim) then
                FLExit("Too many dimensions in boundary condition")
             end if
             bc_component_path=trim(bc_type_path)//"/"//aligned_components(j)//"/order_zero_coefficient"
             surface_field_component=extract_scalar_field(surface_field, j)
             call initialise_field(surface_field_component, bc_component_path, bc_position, &
               time=time)

             bc_component_path=trim(bc_type_path)//"/"//aligned_components(j)//"/order_one_coefficient"
             surface_field_component=extract_scalar_field(surface_field2, j)
             call initialise_field(surface_field_component, bc_component_path, bc_position, &
               time=time)
          end do
          call deallocate(bc_position)

       case("drag")

          call get_boundary_condition(field, i+1, surface_mesh=surface_mesh, &
               surface_element_list=surface_element_list)
          scalar_surface_field => extract_scalar_surface_field(field, bc_name, name="DragCoefficient")
          ! map the coordinate field onto this mesh
          call allocate(bc_position, position%dim, surface_mesh)
          call remap_field_to_surface(position, bc_position, surface_element_list)

          call initialise_field(scalar_surface_field, bc_path_i, bc_position, &
            time=time)
          call deallocate(bc_position)

       case("wind_forcing")

          call get_boundary_condition(field, i+1, surface_mesh=surface_mesh, &
               surface_element_list=surface_element_list)
          surface_field => extract_surface_field(field, bc_name, name="WindSurfaceField")
          ! map the coordinate field onto this mesh
          call allocate(bc_position, position%dim, surface_mesh)
          call remap_field_to_surface(position, bc_position, surface_element_list)

          if (have_option(trim(bc_path_i)//"/wind_stress")) then
             bc_type_path=trim(bc_path_i)//"/wind_stress"
             call initialise_field(surface_field, bc_type_path, bc_position, &
               time=time)
          else if (have_option(trim(bc_path_i)//"/wind_velocity")) then
             bc_type_path=trim(bc_path_i)//"/wind_velocity"

             scalar_surface_field => extract_scalar_surface_field(field, bc_name, name="WindDragCoefficient")
             call initialise_field(scalar_surface_field, &
                  trim(bc_type_path)//"/wind_drag_coefficient", bc_position, &
                  time=time)

             call initialise_field(surface_field, &
                  trim(bc_type_path)//"/wind_velocity", bc_position, &
                  time=time)
          end if
          call deallocate(bc_position)

        case("free_surface")
            
           bc_path_i=trim(bc_path_i)//"/wetting_drying"
           if(have_option(trim(bc_path_i))) then
              scalar_surface_field => extract_scalar_surface_field(field, bc_name, name="WettingDryingAlpha")
              call zero(scalar_surface_field)
              scalar_surface_field => extract_scalar_surface_field(field, bc_name, name="WettingDryingOldAlpha")
              call zero(scalar_surface_field)
           end if

         case ("no_normal_flow", "near_wall_treatment", "log_law_of_wall", "outflow")

          ! nothing to be done (yet?)
          
       case default
          FLAbort("Incorrect boundary condition type for field")
       end select
       
    end do boundary_conditions
    
  end subroutine set_vector_boundary_conditions_values
  
  subroutine set_tidal_bc_value(surface_field, bc_position, bc_type_path, field_name)
    ! tidal_forcing - asc
    type(scalar_field), intent(inout):: surface_field
    type(vector_field), intent(in):: bc_position
    character(len=*), intent(in):: bc_type_path, field_name

    real :: current_time, frequency, amplitude_factor
    integer :: constituent_count, i, j, id, stat
    character(len=3) :: constituent_name
    character(len=4096) :: file_name, variable_name_amplitude, variable_name_phase
    real, dimension(:), allocatable::amplitude, phase
    real :: xyz(3), longitude, latitude
    real :: gravty

    allocate(amplitude(node_count(surface_field)), phase(node_count(surface_field)))
    
    if(have_option(trim(bc_type_path)//"/from_file/tidal")) then
       call set(surface_field, 0.0)
       call get_option("/timestepping/current_time", current_time)
       constituent_count = option_count(trim(bc_type_path)//"/from_file/tidal")

       do i=0, constituent_count-1
          amplitude = 0.0
          phase = 0.0
          call get_option(trim(bc_type_path)//"/from_file/tidal["//int2str(i)//"/amplitude_factor", amplitude_factor, stat=stat)
          if (stat/=0) then
             amplitude_factor = 1.0
          end if
          
          ! Taken from E.W. Schwiderski - Rev. Geophys. Space Phys. Vol. 18 No. 1 pp. 243--268, 1980
          call get_option(trim(bc_type_path)//"/from_file/tidal["//int2str(i)//"]/name", constituent_name, stat=stat)

          frequency = get_tidal_frequency(constituent_name)

          call get_option(trim(bc_type_path)//"/from_file/tidal["//int2str(i)//"]/file_name", file_name, stat=stat)
          call get_option(trim(bc_type_path)//"/from_file/tidal["//int2str(i)//"]/variable_name_amplitude", variable_name_amplitude, stat=stat)
          call get_option(trim(bc_type_path)//"/from_file/tidal["//int2str(i)//"]/variable_name_phase", variable_name_phase, stat=stat)

          call SampleNetCDF_Open(trim(file_name), id)
          call SampleNetCDF_SetVariable(id, trim(variable_name_amplitude))
          do j=1, node_count(bc_position)
             xyz = node_val(bc_position, j)
             call LongitudeLatitude(xyz, longitude, latitude)
             call SampleNetCDF_GetValue(id, longitude, latitude, amplitude(j))
          end do
          call SampleNetCDF_SetVariable(id, trim(variable_name_phase))
          do j=1, node_count(bc_position)
             xyz = node_val(bc_position, j)
             call LongitudeLatitude(xyz, longitude, latitude)
             call SampleNetCDF_GetValue(id, longitude, latitude, phase(j))
          end do
          
          call get_option('/physical_parameters/gravity/magnitude', gravty)
          do j=1, node_count(bc_position)
            if (field_name=="Pressure") then
              call addto(surface_field, j, &
                  gravty * amplitude(j) * cos( frequency * current_time - ( phase(j) * pi / 180 ) ))
            else
              call addto(surface_field, j, &
                  amplitude(j) * cos( frequency * current_time - ( phase(j) * pi / 180 ) ))
            end if
          end do
          call SampleNetCDF_Close(id)
       end do
    end if
  end subroutine set_tidal_bc_value

  subroutine set_nemo_bc_value(state, surface_field, bc_position, bc_type_path, field_name, surface_element_list)
    ! This subroutine sets the pressure at the boundary from NEMO data
    type(state_type), intent(in) :: state
    type(scalar_field), intent(inout):: surface_field
    type(vector_field), intent(in):: bc_position
    integer, dimension(:), intent(in) :: surface_element_list
    character(len=*), intent(in):: bc_type_path, field_name

    type(scalar_field) :: nemo_pressure_bc
    type(scalar_field), pointer :: NEMOpressure
    character(len=4096) :: data_field_name

    integer :: j
    real :: gravty, boundary_pressure
    
    if(have_option(trim(bc_type_path)//"/NEMO_data")) then
       call set(surface_field, 0.0)
          
       
       call get_option(trim(bc_type_path)//"/NEMO_data/field_name", data_field_name)
       call get_option('/physical_parameters/gravity/magnitude', gravty)

       NEMOpressure => extract_scalar_field(state, data_field_name)
       call allocate(nemo_pressure_bc, surface_field%mesh, name="NPBC")
       call remap_field_to_surface(NEMOpressure, nemo_pressure_bc, surface_element_list)

       do j=1, node_count(bc_position)
         boundary_pressure=node_val(nemo_pressure_bc,j)
         if (field_name=="Pressure") then
           call addto(surface_field, j, boundary_pressure)
         else
           call addto(surface_field, j, boundary_pressure/gravty)
         end if
       end do
    end if
  end subroutine set_nemo_bc_value

  subroutine apply_dirichlet_conditions_inverse_mass_vector(inverse_mass, field)
    !!< Zeroes the rows of dirichlet boundary conditions in 
    !!< the inverse mass matrix
    type(block_csr_matrix), intent(inout) :: inverse_mass
    type(vector_field), intent(in) :: field

    logical, dimension(field%dim):: applies
    character(len=FIELD_NAME_LEN):: bctype
    real, dimension(:), pointer:: vals
    integer, dimension(:), pointer:: cols
    integer, dimension(:), pointer:: surface_node_list
    logical, dimension(:), allocatable:: dirichlet_mask
    integer :: i,j,k

    allocate(dirichlet_mask(1:node_count(field)))
    
    bcloop: do i=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, i, type=bctype, &
          surface_node_list=surface_node_list, applies=applies)

       if (bctype/="dirichlet") cycle bcloop
       
       ! first zero rows:
       do j=1,size(surface_node_list)
          do k = 1, field%dim
            if(applies(k)) then
              call zero_row(inverse_mass, k, surface_node_list(j))
            end if
          end do
       end do
      
       ! then zero columns:
       dirichlet_mask=.false.
       dirichlet_mask(surface_node_list)=.true.
       do j=1, node_count(field)
         cols => row_m_ptr(inverse_mass, j)
         do k=1, field%dim
           if (applies(k)) then
             vals => row_val_ptr(inverse_mass, k, k, j)
             where (dirichlet_mask(cols))
               vals=0.0
             end where
           end if
         end do
       end do

    end do bcloop

  end subroutine apply_dirichlet_conditions_inverse_mass_vector

  subroutine apply_dirichlet_conditions_inverse_mass_vector_lumped(inverse_masslump, field)
    !!< Zeroes the rows of dirichlet boundary conditions in 
    !!< the lumped mass vector field
    type(vector_field), intent(inout) :: inverse_masslump
    type(vector_field), intent(in) :: field

    logical, dimension(field%dim):: applies
    character(len=FIELD_NAME_LEN):: bctype
    integer, dimension(:), pointer:: surface_node_list
    integer :: i,j,k

    bcloop: do i=1, get_boundary_condition_count(field)
       call get_boundary_condition(field, i, type=bctype, &
          surface_node_list=surface_node_list, applies=applies)

       if (bctype/="dirichlet") cycle bcloop
       
       do j=1,size(surface_node_list)
          do k = 1, field%dim
            if(applies(k)) then
              
              call set(inverse_masslump, k, surface_node_list(j), 0.0)

            end if
         end do
       end do

    end do bcloop

  end subroutine apply_dirichlet_conditions_inverse_mass_vector_lumped
  
  subroutine ocean_boundaries_stats(state)
  type(state_type), intent(in):: state
    
    if (current_debug_level>=1) then
      call single_ocean_boundary_stats(state, "DistanceToTop", "top")
      call single_ocean_boundary_stats(state, "DistanceToBottom", "bottom")
    end if
    
  end subroutine ocean_boundaries_stats

  subroutine single_ocean_boundary_stats(state, field_name, surface_name)
  type(state_type), intent(in):: state
  character(len=*), intent(in):: field_name, surface_name
  
    type(vector_field), pointer:: gravity_direction, positions
    type(scalar_field), pointer:: distance_field
    real, dimension(:,:), allocatable:: face_normal, grav_normal
    real, dimension(:), allocatable:: detwei_f
    real area, max_inn, min_inn, inn
    integer, dimension(:), pointer:: surface_element_list, surface_node_list
    integer i, j, ele, sele, sngi
          
    ewrite(1,*) "Stats for "//trim(surface_name)//" ocean boundary:"
    
    distance_field => extract_scalar_field(state, field_name)
    gravity_direction => extract_vector_field(state, "GravityDirection")
    positions => extract_vector_field(state, "Coordinate")
    call get_boundary_condition(distance_field, 1, &
      surface_element_list=surface_element_list, &
      surface_node_list=surface_node_list)
      
    ewrite(2,*) "Number of surface nodes:", size(surface_node_list)
    ewrite(2,*) "Number of surface elements:", size(surface_element_list)
    
    sngi = face_ngi(distance_field, 1)
    
    allocate( detwei_f(1:sngi), &
      face_normal(1:gravity_direction%dim, 1:sngi), &
      grav_normal(1:gravity_direction%dim, 1:sngi) )
    
    area=0.0
    min_inn=huge(1.0)
    max_inn=-min_inn
    do i=1, size(surface_element_list)
      
      sele=surface_element_list(i)
      ! 3D element behind the surface element:
      ele=face_ele(distance_field, sele)
      
      call transform_facet_to_physical(positions, sele, &
        detwei_f=detwei_f, normal=face_normal)
        
      area=area+sum(detwei_f)
      
      ! gravity normal at face gauss points
      grav_normal=face_val_at_quad(gravity_direction, sele)
      
      ! inner product of face_normal and grav_normal at gauss points:
      do j=1, sngi
        inn=dot_product(face_normal(:,j), grav_normal(:,j))
        
        min_inn=min(min_inn, inn)
        max_inn=max(max_inn, inn)
      end do
      
    end do
      
    ewrite(2, *) "Surface area:", area
    ewrite(2, *) "Maximum of inner product gravity and face normal", min_inn
    ewrite(2, *) "Minimum of inner product gravity and face normal", max_inn
    ewrite(2, *) ""
    
  end subroutine single_ocean_boundary_stats

 subroutine populate_ocean_forcing_boundary_conditions(state)
    type(state_type), intent(in) :: state
    
    integer, dimension(:), pointer :: surface_element_list, surface_nodes
    type(scalar_field), pointer :: distanceToTop
    type(scalar_field) :: scalar_surface_field
    type(vector_field) :: vector_surface_field
    type(scalar_field), pointer :: scalar_source_field
    type(vector_field), pointer :: vector_source_field
    type(mesh_type) :: ocean_mesh, input_mesh
    type(mesh_type), pointer :: surface_mesh
    character(len=FIELD_NAME_LEN) input_mesh_name
    integer, dimension(:), allocatable:: surface_ids
    integer, dimension(2) :: shape_option

    distanceToTop => extract_scalar_field(state, "DistanceToTop")
    call get_option('/ocean_forcing/mesh_choice/mesh/name', input_mesh_name)
    input_mesh = extract_mesh(state, input_mesh_name);
    call get_boundary_condition(distanceToTop, name='top', surface_element_list=surface_element_list)
    call create_surface_mesh(ocean_mesh, surface_nodes, input_mesh, surface_element_list, 'OceanSurface')

    ! Add a BC to the relavent field for each BC enabled
    shape_option=option_shape("geometry/ocean_boundaries/top_surface_ids")
    allocate(surface_ids(1:shape_option(1)))
    call get_option("geometry/ocean_boundaries/top_surface_ids", surface_ids)

    if (have_option('/ocean_forcing/surface_stress')) then
      vector_source_field => extract_vector_field(state, 'Velocity')
      call add_boundary_condition(vector_source_field, 'ocean_forcing_boundary_stress',&
                                  'wind_forcing', surface_ids)
      call get_boundary_condition(vector_source_field, 'ocean_forcing_boundary_stress', &
                                  surface_mesh=surface_mesh)
      call allocate(vector_surface_field, 2, surface_mesh, &
                    name="WindSurfaceField")
      call insert_surface_field(vector_source_field, 'ocean_forcing_boundary_stress', &
                                vector_surface_field)
      call deallocate(vector_surface_field)
    end if
    if (have_option('/ocean_forcing/temperature_flux')) then
      scalar_source_field => extract_scalar_field(state, 'Temperature')
      call add_boundary_condition(scalar_source_field, 'ocean_forcing_boundary_heat', &
                                  'neumann', surface_ids)
      call get_boundary_condition(scalar_source_field, 'ocean_forcing_boundary_heat', &
                                  surface_mesh=surface_mesh)
      call allocate(scalar_surface_field, surface_mesh, name="value");
      call insert_surface_field(scalar_source_field, 'ocean_forcing_boundary_heat', &
                                scalar_surface_field)
      call deallocate(scalar_surface_field)
    end if
    if (have_option('/ocean_forcing/salinity_flux')) then
      scalar_source_field => extract_scalar_field(state, 'Salinity')
      call add_boundary_condition(scalar_source_field, 'ocean_forcing_boundary_salinity',&
                                  'neumann', surface_ids)
      call get_boundary_condition(scalar_source_field, 'ocean_forcing_boundary_salinity', &
                                  surface_mesh=surface_mesh)
      call allocate(scalar_surface_field, surface_mesh, name="value")
      call insert_surface_field(scalar_source_field, 'ocean_forcing_boundary_salinity', &
                                scalar_surface_field)
      call deallocate(scalar_surface_field)
    end if
    if (have_option('/ocean_forcing/solar_radiation')) then
       scalar_source_field => extract_scalar_field(state, 'PhotosyntheticRadiation')
       call add_boundary_condition(scalar_source_field, 'ocean_forcing_boundary_solar', &
                                   'weakdirichlet', surface_ids)
       call get_boundary_condition(scalar_source_field, 'ocean_forcing_boundary_solar', &
                                  surface_mesh=surface_mesh)                            
       call allocate(scalar_surface_field, surface_mesh, name="value")
       call insert_surface_field(scalar_source_field, 'ocean_forcing_boundary_solar', &
                                scalar_surface_field)
       call deallocate(scalar_surface_field)
    end if

    call deallocate(ocean_mesh)
    deallocate(surface_ids)
    
  end subroutine populate_ocean_forcing_boundary_conditions

  subroutine set_ocean_forcings_boundary_conditions(state)
    type(state_type), intent(in) :: state

    type(mesh_type) :: ocean_mesh, input_mesh
    ! output from the get_fluxes call on the ocean_mesh
    type(scalar_field) :: salinity_flux, heat_flux, solar_flux
    type(vector_field) :: stress_flux
    ! the current state to be put on the ocean_mesh - input to the fluxes call
    type(scalar_field) :: temperature, salinity
    type(vector_field) :: velocity, position
    ! these are pointers to the fields in the state
    type(scalar_field), pointer :: p_temperature, p_salinity
    type(vector_field), pointer :: p_velocity, p_position

    ! some temporary storage arrays
    ! some of this need hiding inside get_fluxes
    real, dimension(3) :: temp_vector_3D, transformation
    real, dimension(2) :: temp_vector_2D
    real, dimension(:), allocatable :: temp, sal, X, Y, Z, Vx, Vy, Vz
    real, dimension(:), allocatable :: F_as, Q_as, Q_s, Tau_u, Tau_v
    real, dimension(:), allocatable :: d1, d2, d3, d4, d5
    integer :: NNodes, i
    integer, dimension(:), pointer :: surface_element_list, surface_nodes
    real :: current_time
    type(scalar_field), pointer :: scalar_source_field, sfield
    type(vector_field), pointer :: vector_source_field, vfield
    type(scalar_field), pointer :: scalar_surface
    type(vector_field), pointer :: vector_surface
    character(len=FIELD_NAME_LEN) input_mesh_name
    type(scalar_field), pointer :: distanceToTop
    logical*1 :: on_sphere  ! needs to be handed over to C, which has 1 bit booleans
    real, dimension(:), allocatable:: lat_long
    integer shape_option(2)
    integer stat

    distanceToTop => extract_scalar_field(state, "DistanceToTop")
    call get_option('/ocean_forcing/mesh_choice/mesh/name', input_mesh_name)
    input_mesh = extract_mesh(state, input_mesh_name);
    call get_boundary_condition(distanceToTop, name='top', surface_element_list=surface_element_list)
    call create_surface_mesh(ocean_mesh, surface_nodes, input_mesh, surface_element_list, 'OceanSurface')
    NNodes = node_count(ocean_mesh) 

    ! temp arrays for fluxes
    allocate(temp(NNodes))
    allocate(sal(NNodes))
    allocate(X(NNodes))
    allocate(Y(NNodes))
    allocate(Z(NNodes))
    allocate(Vx(NNodes))
    allocate(Vy(NNodes))
    allocate(Vz(NNodes))
    allocate(F_as(NNodes))
    allocate(Q_as(NNodes))
    allocate(Q_s(NNodes))
    allocate(Tau_u(NNodes))
    allocate(Tau_v(NNodes))
    allocate(d1(NNodes))
    allocate(d2(NNodes))
    allocate(d3(NNodes))
    allocate(d4(NNodes))
    allocate(d5(NNodes))


    ! allocate field on ocean mesh to store output of get_fluxes
    if (have_option('/ocean_forcing/surface_stress')) then
      call allocate(stress_flux, 2, ocean_mesh, name="stress_flux")
    end if
    if (have_option('/ocean_forcing/temperature_flux')) then
        call allocate(heat_flux, ocean_mesh, name="heat_flux")
    end if
    if (have_option('/ocean_forcing/salinity_flux')) then
        call allocate(salinity_flux, ocean_mesh, name="salinity_flux")
    end if
    if (have_option('/ocean_forcing/solar_radiation')) then
        call allocate(solar_flux, ocean_mesh, name="solar_flux")
    end if
    ! allocate space to store current state of parameters required
    ! to get the fluxes
    call allocate(temperature, ocean_mesh, name="temperature")
    call allocate(velocity, 3, ocean_mesh, name="velocity")
    call allocate(position, 3, ocean_mesh, name="position")
    call allocate(salinity, ocean_mesh, name="salinity")

    ! grab current state, this needs doing regardless of which BCs
    ! are applied
    p_temperature => extract_scalar_field(state, "Temperature")
    p_velocity => extract_vector_field(state, "Velocity")
    p_position => extract_vector_field(state, "Coordinate")

    
    ! remap modelled params onto the appropriate field in ocean_mesh
    call remap_field_to_surface(p_temperature, temperature, &
                                surface_element_list)
    call remap_field_to_surface(p_velocity, velocity, &
                                surface_element_list)
    call remap_field_to_surface(p_position, position, &
                                surface_element_list)

    ! check if we are transforming the coordinates to a specified lat/long
    if (have_option('/ocean_forcing/position')) then
        shape_option=option_shape('/ocean_forcing/position')
        if (shape_option(1) .ne. 2) then
            FLExit("Only specify a latitude and longitude under /ocean_forcing/positions, i.e. two numbers expected.")
        end if
        allocate(lat_long(1:shape_option(1)))
        call get_option('/ocean_forcing/position', lat_long)
        transformation(1) = lat_long(2) ! Longtitude
        transformation(2) = lat_long(1) ! latitude
        transformation(3) = 0.0
        call projections_spherical_cartesian(1, transformation(1), transformation(2), transformation(3))
    else
        transformation = 0.0
    end if

    ! we now have the modelled parameters, temperature, etc on the same
    ! mesh as the surface - we can now grab these, so...
    ! loop over surface mesh points, grabbing field values at each and
    ! shoving them unceremoniously into my temporary arrays
    do i=1,NNodes
      temp_vector_3D = node_val(position,i)
      X(i) = temp_vector_3D(1)+transformation(1) 
      Y(i) = temp_vector_3D(2)+transformation(2)
      Z(i) = temp_vector_3D(3)+transformation(3)
      temp(i) = node_val(temperature,i)
      temp_vector_3D = node_val(velocity,i)
      Vx(i) = temp_vector_3D(1) 
      Vy(i) = temp_vector_3D(2)
      Vz(i) = temp_vector_3D(3)
    end do

    ! finally, check if the single position option is on, if so, make all
    ! positions the same
    if (have_option('/ocean_forcing/position/single_location')) then
        do i=1, NNodes
            X(i) = transformation(1) 
            Y(i) = transformation(2)
            Z(i) = transformation(3)
        end do
    end if

    if (have_option('/ocean_forcing/salinity_flux')) then
        ! we only need to worry about salinity if the flux is on
        p_salinity => extract_scalar_field(state, "Salinity", stat)
        if (stat /= 0) then
            FLExit("If you switch on a salinity flux, you'd better have a Salinity field...")
        end if
        call remap_field_to_surface(p_salinity, salinity, &
                                surface_element_list)
        do i=1,NNodes
            sal(i) = node_val(salinity,i)
        end do
    else
        ! use a decent estimate for the surface salinity
        do i=1,NNodes
            sal(i) = 35.0
        end do
    end if

    call get_option("/timestepping/current_time", current_time)
    on_sphere=have_option('/geometry/spherical_earth')

    call get_era40_fluxes(current_time, X, Y, Z, temp, Vx, Vy, Vz, sal, &
                          F_as, Q_as, Tau_u, Tau_v, Q_s, &
                          NNodes,on_sphere,0)

    ! finally, we need to reverse-map the temporary fields on the ocean mesh
    ! to the actual fields in state
    ! using remap for now, but this assumes the same type of field
    ! otherwise need to project
    if (have_option('/ocean_forcing/surface_stress')) then
        do i=1,NNodes
            temp_vector_2D(1) = Tau_u(i)
            temp_vector_2D(2) = Tau_v(i)
            call set(stress_flux,i,temp_vector_2D)
        end do
        vector_source_field => extract_vector_field(state, 'Velocity')
        vector_surface => extract_surface_field(vector_source_field, 'ocean_forcing_boundary_stress',&
                                             "WindSurfaceField")
        call remap_field(stress_flux, vector_surface)
        if (have_option("/ocean_forcing/output_fluxes_diagnostics/vector_field::MomentumFlux")) then
            vector_source_field => extract_vector_field(state, 'MomentumFlux')
            ! copy the values onto the mesh using the global node id
            do i=1,size(surface_nodes)
                call set(vector_source_field,1,surface_nodes(i),node_val(vector_surface,1,i))
            end do
            do i=1,size(surface_nodes)
                call set(vector_source_field,2,surface_nodes(i),node_val(vector_surface,2,i))
            end do
            do i=1,size(surface_nodes)
                call set(vector_source_field,3,surface_nodes(i),0.0)
            end do
            vfield => extract_vector_field(state, 'OldMomentumFlux',stat)
           if (stat == 0) then
                call set(vfield,vector_source_field)
            end if
        end if
        call deallocate(stress_flux)                            
    end if
    if (have_option('/ocean_forcing/temperature_flux')) then
        do i=1,NNodes
           call set(heat_flux,i,Q_as(i))
        end do 
        scalar_source_field => extract_scalar_field(state, 'Temperature')
        scalar_surface => extract_surface_field(scalar_source_field, 'ocean_forcing_boundary_heat',&
                                             "value")
        call remap_field(heat_flux, scalar_surface)
        if (have_option("/ocean_forcing/output_fluxes_diagnostics/scalar_field::HeatFlux")) then
            scalar_source_field => extract_scalar_field(state, 'HeatFlux')
            ! copy the values onto the mesh using the global node id
            do i=1,size(surface_nodes)
                call set(scalar_source_field,surface_nodes(i),node_val(scalar_surface,i))
            end do   
            sfield => extract_scalar_field(state, 'OldHeatFlux',stat)
            if (stat == 0) then
                call set(sfield,scalar_source_field)
            end if
        end if
        call deallocate(heat_flux)
    end if
    if (have_option('/ocean_forcing/salinity_flux')) then
        do i=1,NNodes
            call set(salinity_flux,i,F_as(i))
        end do
        scalar_source_field => extract_scalar_field(state, 'Salinity')
        scalar_surface => extract_surface_field(scalar_source_field, 'ocean_forcing_boundary_salinity',&
                                             "value")
        call remap_field(salinity_flux, scalar_surface)
        if (have_option("/ocean_forcing/output_fluxes_diagnostics/scalar_field::SalinityFlux")) then
            scalar_source_field => extract_scalar_field(state, 'SalinityFlux')
            ! copy the values onto the mesh using the global node id
            do i=1,size(surface_nodes)
                call set(scalar_source_field,surface_nodes(i),node_val(scalar_surface,i))
            end do 
            sfield => extract_scalar_field(state, 'OldSalinityFlux',stat)
            if (stat == 0) then
                call set(sfield,scalar_source_field)
            end if
        end if
        call deallocate(salinity_flux)
    end if
    if (have_option('/ocean_forcing/solar_radiation')) then
        do i=1,NNodes
            call set(solar_flux,i,Q_s(i))
        end do
        scalar_source_field => extract_scalar_field(state, 'PhotosyntheticRadiation')
        scalar_surface => extract_surface_field(scalar_source_field, 'ocean_forcing_boundary_solar',&
                                             "value")
        call remap_field(solar_flux, scalar_surface)
        if (have_option("/ocean_forcing/output_fluxes_diagnostics/scalar_field::PhotosyntheticRadiationDownward")) then
            scalar_source_field => extract_scalar_field(state, 'PhotosyntheticRadiationDownward')
            ! copy the values onto the mesh using the global node id
            do i=1,size(surface_nodes)
                call set(scalar_source_field,surface_nodes(i),node_val(scalar_surface,i))
            end do 
            sfield => extract_scalar_field(state, 'OldPhotosyntheticRadiationDownward',stat)
            if (stat == 0) then
                call set(sfield,scalar_source_field)
            end if
        end if
        call deallocate(solar_flux)
    end if
    
    call deallocate(temperature)
    call deallocate(salinity)
    call deallocate(velocity)
    call deallocate(position) 
    call deallocate(ocean_mesh)
    deallocate(temp)
    deallocate(sal)
    deallocate(X)
    deallocate(Y)
    deallocate(Z)
    deallocate(Vx)
    deallocate(Vy)
    deallocate(Vz)
    deallocate(F_as)
    deallocate(Q_as)
    deallocate(Q_s)
    deallocate(Tau_u)
    deallocate(Tau_v)
    deallocate(d1)
    deallocate(d2)
    deallocate(d3)
    deallocate(d4)
    deallocate(d5)

  end subroutine set_ocean_forcings_boundary_conditions

  subroutine set_dirichlet_turbine_boundary_conditions(state)
    type(state_type), intent(in) :: state
    type(vector_field), pointer :: vel_field,  coord_field
    type(vector_field), pointer :: surface_field
    type(vector_field) :: bc_positions
    type(scalar_field) :: surface_field_normal
    type(scalar_field), pointer :: fs_field

    integer :: notur, i, j,ele, stat
    character(len=FIELD_NAME_LEN) :: turbine_path
    character(len=FIELD_NAME_LEN), dimension(2) ::  bc_name
    character(len=PYTHON_FUNC_LEN) :: func
    real, dimension(3) :: picker
    real,dimension(:),allocatable :: local_coord
    real, dimension(2) :: fs_val, surface_area
    real :: flux
    integer, dimension(:), pointer :: surface_element_list
    real :: time
    logical :: have_fs
    ewrite(1,*) "In dirichlet turbine model"

    have_fs=.False.
    ! We made sure in populate_boundary_conditions that FreeSurface exists.
    vel_field => extract_vector_field(state, "Velocity")
    fs_field => extract_scalar_field(state, "FreeSurface", stat)
    if(stat==0) have_fs = .True.
    coord_field => extract_vector_field(state, "Coordinate")
    allocate(local_coord(coord_field%dim+1))

    ! loop through turbines
    notur = option_count("/turbine_model/turbine")
    do i=0, notur-1
       turbine_path="/turbine_model/turbine["//int2str(i)//"]"
       if (.not. have_option(trim(turbine_path)//"/dirichlet")) cycle
       ! We need a FreeSurface field for the dirchlet turbine model.
       if (.not. have_fs) FLExit("Turbine error: No FreeSurface field found.")
       do j=1,2
         call get_option(trim(turbine_path)//"/dirichlet/boundary_condition_name_"//int2str(j)//"/name", bc_name(j))
       end do

       call get_option("/timestepping/current_time", time)

       ! Get free surface values at the user specified points
       do j=1,2
         call get_option(trim(turbine_path)//"/dirichlet/free_surface_point_"//int2str(j), picker(1:coord_field%dim))
         if (coord_field%dim==1) then
             call picker_inquire(coord_field, coordx=picker(1), ele=ele, local_coord=local_coord, global=.true.)
         elseif (coord_field%dim==2) then
             call picker_inquire(coord_field, coordx=picker(1), coordy=picker(2), ele=ele, local_coord=local_coord, global=.true.)
         elseif (coord_field%dim==3) then
             call picker_inquire(coord_field, picker(1), picker(2), picker(3), ele, local_coord, global=.true.)
         end if
         if (ele<0) then
             FLExit("Turbine error: The point defined in "//trim(turbine_path)//"/free_surface_point_"//int2str(j)//" is not located in a mesh element")
         end if
         fs_val(j) = eval_field_scalar(ele, fs_field, local_coord)
       end do
       ! Function head -> outflow
       call get_option(trim(turbine_path)//"/dirichlet/head_flux", func)
       call real_from_python(func, fs_val(1)-fs_val(2), flux)

       ! Get surface areas of the boundaries
       ! This could be done only once if no free surface mesh movement is performed.
       do j=1,2
           surface_field => extract_surface_field(vel_field, bc_name(j), name="value")
           call get_boundary_condition(vel_field, bc_name(j), surface_element_list=surface_element_list)
           ! Map the coord_field field to the surface mesh
           call allocate(bc_positions, coord_field%dim, surface_field%mesh)
           call remap_field_to_surface(coord_field, bc_positions, surface_element_list)
           surface_area(j)=get_surface_area(bc_positions)
           call deallocate(bc_positions)
       end do

       ! Overwrite the normal component of the dichilet boundaries
       do j=1,2
         surface_field => extract_surface_field(vel_field, bc_name(j), name="value")
         ! Extract the scalar normal component
         surface_field_normal = extract_scalar_field(surface_field, 1)
         ! speed = outflow/area
         ewrite(3,*) "Surface area of boundary ", trim(bc_name(j)), ": ",  surface_area(j)
         ewrite(3,*) "Setting normal flow speed at boundary ", trim(bc_name(j)), " to: ", flux/surface_area(j)
         call set(surface_field_normal, (-1)**(j+1)*flux/surface_area(j))
       end do
    end do
    ! Tidying up
    deallocate(local_coord)
    ewrite(3,*) "Out dirichlet turbine model."

  contains 
  ! Calculates the surface area by integrating the unit function over the domain defined by "positions"
  function get_surface_area(positions) result(surface_area)
       type(vector_field), intent(in) :: positions
       integer :: i
       real :: surface_area
       real, dimension(1:ele_ngi(positions%mesh,1)):: detwei

       surface_area=0.0
       do i=1, element_count(positions)
          call transform_to_physical(positions, i, detwei)
          surface_area=surface_area+sum(detwei)
       end do
   end function get_surface_area

  end subroutine set_dirichlet_turbine_boundary_conditions

  subroutine set_flux_turbine_boundary_conditions(state)
    type(state_type), intent(in) :: state
    type(vector_field), pointer :: vel_field,  coord_field, surface_field
    type(scalar_field) :: scalar_surface_field

    integer :: notur, i, j
    character(len=FIELD_NAME_LEN) :: turbine_path, turbine_type
    character(len=FIELD_NAME_LEN), dimension(2) ::  bc_name
    character(len=PYTHON_FUNC_LEN) :: func
    real :: flux, h
    logical :: active

    ewrite(1,*) "In flux turbine model"
    vel_field => extract_vector_field(state, "Velocity")
    coord_field => extract_vector_field(state, "Coordinate")

    ! loop through turbines
    notur = option_count("/turbine_model/turbine")
    do i=0, notur-1
       turbine_path="/turbine_model/turbine["//int2str(i)//"]"
       if (.not. have_option(trim(turbine_path)//"/flux")) cycle
       if (have_option(trim(turbine_path)//"/flux/penalty")) then
          turbine_type="penalty"
       else
          turbine_type="dg"
       end if
       do j=1,2
         call get_option(trim(turbine_path)//"/flux/boundary_condition_name_"//int2str(j)//"/name", bc_name(j))
       end do
       ! TODO: Set h to pressure jump or free surface jump?
       h=1.0
       call get_option(trim(turbine_path)//"/flux/"//trim(turbine_type)//"/factor", func)
       call real_from_python(func, h, flux) 

       ! Set the turbine dg fluxes
       do j=1,2
         surface_field => extract_surface_field(vel_field, trim(bc_name(j))//"_turbine", name="value")
         scalar_surface_field = extract_scalar_field(surface_field, 1)
         !scalar_surface_field => extract_scalar_surface_field(vel_field, trim(bc_name(j))//"_dgflux", name="value")
         call set(scalar_surface_field, flux)
       end do

       ! Decide if turbine is active or not
       if (have_option(trim(turbine_path)//"/flux/always_on")) then
           active=.true.
       elseif (have_option(trim(turbine_path)//"/flux/always_off")) then
           active=.false.
       end if

       ! De/Activate the Dirichlet conditions and turbine boundary condition by changing their "applied" flag
       if (active) then
         do j=1,2
            call set_boundary_condition_applies_flag(vel_field, trim(bc_name(j))//"_turbine", (/.true., .true., .true./))
            call set_boundary_condition_applies_flag(vel_field, trim(bc_name(j)), (/.false., .false., .false./))
         end do
       else
         do j=1,2
            call set_boundary_condition_applies_flag(vel_field, trim(bc_name(j))//"_turbine", (/.false., .false., .false./))
            call set_boundary_condition_applies_flag_from_options_path(vel_field, trim(bc_name(j))) 
         end do
       end if
    end do
    ewrite(3,*) "Out flux turbine model."


   contains 
   subroutine set_boundary_condition_applies_flag(field, name, applies)
     type(vector_field), intent(in):: field
     ! Name of the boundary condition
     character(len=*), intent(in):: name
     logical, dimension(3):: applies
     integer i
     do i=1, size(field%bc%boundary_condition)
       if (field%bc%boundary_condition(i)%name==name) then
         field%bc%boundary_condition(i)%applies=applies
         return
       end if
     end do
     FLExit("Internal Error: Field "//trim(field%name)//" does not have boundary condition of name "//trim(name))
   end subroutine set_boundary_condition_applies_flag

   ! Sets the "applies" flags according to the option path.
   subroutine set_boundary_condition_applies_flag_from_options_path(field, name)
     type(vector_field), intent(in):: field
     ! Name of the boundary condition
     character(len=*), intent(in):: name
     logical, dimension(3):: applies
     logical :: bc_found
     integer i, j
     character(len=OPTION_PATH_LEN) bc_path, bc_component_path
     character(len=OPTION_PATH_LEN), dimension(3) :: components

     do i=1, size(field%bc%boundary_condition)
       if (field%bc%boundary_condition(i)%name==name) then
         bc_path=field%bc%boundary_condition(i)%option_path
         bc_found=.true.
         exit
       end if
     end do
     if (.not. bc_found) FLExit("Internal Error: Field "//trim(field%name)//" does not have boundary condition of name "//trim(name))
     bc_path=trim(bc_path)//"/type"
     if(have_option(trim(bc_path)//"/align_bc_with_cartesian")) then
             bc_path=trim(bc_path)//"/align_bc_with_cartesian"
             components=(/ &
                "x_component", &
                "y_component", &
                "z_component" /)
     else
             bc_path=trim(bc_path)//"/align_bc_with_surface"
             components=(/ &
                "normal_component   ", &
                "tangent_component_1", &
                "tangent_component_2" /)
     end if
     do j=1,3
             bc_component_path=trim(bc_path)//"/"//trim(components(j))
             applies(j)=have_option(trim(bc_component_path))
     end do
     ! Finally set the "applies" flags
     field%bc%boundary_condition(i)%applies=applies
   end subroutine set_boundary_condition_applies_flag_from_options_path

  end subroutine set_flux_turbine_boundary_conditions


  subroutine initialise_rotated_bcs(surface_element_list, x, & 
    debugging_mode, normal, tangent_1, tangent_2)
    
    integer, dimension(:),intent(in):: surface_element_list
    ! vector fields on the surface mesh
    type(vector_field),intent(inout):: normal
    type(vector_field),intent(inout), optional::tangent_1, tangent_2
    ! positions on the entire mesh (may not be same order as surface mesh!!)
    type(vector_field),intent(in)   :: x
    logical, intent(in):: debugging_mode
      
    type(vector_field)              :: bc_position

    real, dimension(x%dim, face_ngi(x, 1)):: normal_bdy
    real, dimension(face_ngi(x, 1))       :: detwei_bdy

    integer                       :: i, bcnod
    integer                       :: sele
    integer                       :: t1_max
    real, dimension(x%dim)        :: t1, t2, t1_norm, n
    real                          :: proj1, det


    do i=1, size(surface_element_list)
       
       sele = surface_element_list(i)
       
       call transform_facet_to_physical(x, sele, &
            detwei_f=detwei_bdy, normal=normal_bdy)
       
       call addto(normal, ele_nodes(normal, i), &
         shape_vector_rhs(ele_shape(normal, i), normal_bdy, detwei_bdy))
    end do ! surface_element_list

    t1 = 0.0
    t2 = 0.0

    bcnod=normal%mesh%nodes
    do i=1, bcnod

       ! get node normal
       n=node_val(normal,i)

       ! normalise it
       n=n/sqrt(sum(n**2))

       call set(normal, i, n)
       if (x%dim>1) then
       
         assert(present(tangent_1))
           
         t1_max=minloc( abs(n), dim=1 )
         t1_norm=0.
         t1_norm(t1_max)=1.
       
         proj1=dot_product(n, t1_norm)
         t1= t1_norm - proj1 * n
       
         ! normalise it
         t1=t1/sqrt(sum(t1**2))

         call set( tangent_1, i, t1 )

         if (x%dim>2)then
          
            assert(present(tangent_2))
            
            t2 = cross_product(n, t1)
          
            call set( tangent_2, i, t2 )
          
         end if
         
       endif

       ! dump normals when debugging
       if (debugging_mode) then
        
          det = abs( &
                n(1) * (t1(2) * t2(3) - t2(2) * t1(3) ) + &
                t1(1) * (t2(2) *  n(3) -  n(2) * t2(3) ) + &
                t2(1) * ( n(2) * t1(3) - t1(2) *  n(3) ) )
          if (  abs( det - 1.) > 1.e-5) then
            call allocate(bc_position, normal%dim, normal%mesh, "BoundaryPosition")
            call remap_field_to_surface(x, bc_position, surface_element_list)
            call vtk_write_fields( "normals", 0, bc_position, bc_position%mesh, &
                  vfields=(/ normal, tangent_1, tangent_2/))
            call deallocate(bc_position)
            ewrite(-1,*) "rotation matrix determinant", det
            FLExit("rotation matrix is messed up, rotated bcs have exploded...")
          end if
       end if
       
    end do ! bcnod

    ! dump normals when debugging
    if (debugging_mode) then
      call allocate(bc_position, normal%dim, normal%mesh, "BoundaryPosition")
      call remap_field_to_surface(x, bc_position, surface_element_list)
      call vtk_write_fields( "normals", 0, bc_position, bc_position%mesh, &
                vfields=(/ normal, tangent_1, tangent_2/))
      call deallocate(bc_position)
    end if
    
  end subroutine initialise_rotated_bcs

  subroutine populate_gls_boundary_conditions(state)
    type(state_type), intent(in) :: state
    
    integer, dimension(:), pointer :: surface_element_list, surface_nodes
    type(scalar_field), pointer :: distanceToTop, distanceToBottom
    type(scalar_field) :: scalar_surface_field
    type(scalar_field), pointer :: tke, psi
    type(mesh_type) :: ocean_mesh, input_mesh
    type(mesh_type), pointer :: surface_mesh
    character(len=FIELD_NAME_LEN) :: input_mesh_name
    character(len=FIELD_NAME_LEN) :: bc_type
    integer, dimension(:), allocatable:: surface_ids
    integer, dimension(2) :: shape_option
    integer :: stat

    ewrite(1,*) "Initialising GLS stable boundaries"

    distanceToTop => extract_scalar_field(state, "DistanceToTop")
    distanceToBottom => extract_scalar_field(state, "DistanceToBottom")

    tke  => extract_scalar_field(state, "GLSTurbulentKineticEnergy",stat)
    if(stat/=0) FLAbort("Need GLSTurbulentKineticEnergy field")
    psi => extract_scalar_field(state, "GLSGenericSecondQuantity",stat)
    if(stat/=0) FLAbort("Need GLSGenericSecondQuantity field")
    
    call get_option(trim(tke%option_path)//'/prognostic/mesh/name', input_mesh_name)
    input_mesh = extract_mesh(state, input_mesh_name);
    
    ! Add a top BC to the two GLS fields
    call get_boundary_condition(distanceToTop, name='top', surface_element_list=surface_element_list)
    call create_surface_mesh(ocean_mesh, surface_nodes, input_mesh, surface_element_list, 'OceanSurface')

    shape_option=option_shape("geometry/ocean_boundaries/top_surface_ids")
    allocate(surface_ids(1:shape_option(1)))
    call get_option("geometry/ocean_boundaries/top_surface_ids", surface_ids)

    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/calculate_boundaries/", bc_type)
    
    call add_boundary_condition(tke, 'tke_top_boundary', bc_type, surface_ids)
    call get_boundary_condition(tke, 'tke_top_boundary', surface_mesh=surface_mesh)
    call allocate(scalar_surface_field, surface_mesh, name="value")
    call insert_surface_field(tke, 'tke_top_boundary', scalar_surface_field)
    call deallocate(scalar_surface_field)

    call add_boundary_condition(psi, 'psi_top_boundary', bc_type, surface_ids)
    call get_boundary_condition(psi, 'psi_top_boundary', surface_mesh=surface_mesh)
    call allocate(scalar_surface_field, surface_mesh, name="value")
    call insert_surface_field(psi, 'psi_top_boundary', scalar_surface_field)
    call deallocate(scalar_surface_field)
    deallocate(surface_ids)

    call deallocate(ocean_mesh)

    ! Bottom boundary on tke and psi
    call get_boundary_condition(distanceToBottom, name='bottom', surface_element_list=surface_element_list)
    call create_surface_mesh(ocean_mesh, surface_nodes, input_mesh, surface_element_list, 'OceanBottom')

    shape_option=option_shape("geometry/ocean_boundaries/bottom_surface_ids")
    allocate(surface_ids(1:shape_option(1)))
    call get_option("geometry/ocean_boundaries/bottom_surface_ids", surface_ids)

    call add_boundary_condition(tke, 'tke_bottom_boundary', bc_type, surface_ids)
    call get_boundary_condition(tke, 'tke_bottom_boundary', surface_mesh=surface_mesh)
    call allocate(scalar_surface_field, surface_mesh, name="value")
    call insert_surface_field(tke, 'tke_bottom_boundary', scalar_surface_field)
    call deallocate(scalar_surface_field)

    call add_boundary_condition(psi, 'psi_bottom_boundary', bc_type, surface_ids)
    call get_boundary_condition(psi, 'psi_bottom_boundary', surface_mesh=surface_mesh)
    call allocate(scalar_surface_field, surface_mesh, name="value")
    call insert_surface_field(psi, 'psi_bottom_boundary', scalar_surface_field)
    call deallocate(scalar_surface_field)
    call deallocate(ocean_mesh)
    deallocate(surface_ids)
    
  end subroutine populate_gls_boundary_conditions


  subroutine populate_flux_turbine_boundary_conditions(state)
    type(state_type), intent(in) :: state
    
    type(vector_field), pointer :: velocity
    character(len=OPTION_PATH_LEN) :: bc_path_i 
    character(len=FIELD_NAME_LEN) :: bc_name_ij, turbine_type
    integer :: i, j, nbtur
 
    ewrite(1,*) "Populate flux turbine boundaries"
    velocity => extract_vector_field(state, "Velocity")
    nbtur=option_count("/turbine_model/turbine")
  
    do i=0, nbtur-1
      bc_path_i="/turbine_model/turbine["//int2str(i)//"]"
      if (.not. have_option(trim(bc_path_i)//"/flux"))  cycle
      if (have_option(trim(bc_path_i)//"/flux/penalty")) then 
          turbine_type="penalty"
      else 
          turbine_type="dg" 
      end if
      do j=1,2
          call get_option(trim(bc_path_i)//"/flux/boundary_condition_name_"//int2str(j)//"/name", bc_name_ij)
          call insert_flux_turbine_boundary_condition(velocity, bc_name_ij, turbine_type)
      end do 
    end do

    contains 
    subroutine insert_flux_turbine_boundary_condition(field, bc_name, turbine_type)
       character(len=FIELD_NAME_LEN), intent(in) :: bc_name 
       type(vector_field), pointer :: field
       integer, dimension(:), allocatable:: surface_ids
       character(len=FIELD_NAME_LEN) :: bc_type, turbine_type ! turbine type is either "dg" or "penalty"
       character(len=OPTION_PATH_LEN) :: option_path
       integer, dimension(2) :: shape_option
       !type(scalar_field) :: scalar_surface_field
       type(vector_field) :: surface_field
       type(mesh_type), pointer :: surface_mesh

       bc_type="turbine_flux_"//trim(turbine_type)
       call get_boundary_condition(field, bc_name, option_path=option_path)
       shape_option=option_shape(trim(option_path)//"/surface_ids")
       allocate(surface_ids(1:shape_option(1)))
       call get_option(trim(option_path)//"/surface_ids", surface_ids)
       call add_boundary_condition(velocity, trim(bc_name)//"_turbine", bc_type, surface_ids)     
       deallocate(surface_ids)
 
       call get_boundary_condition(field, trim(bc_name)//"_turbine", surface_mesh=surface_mesh)
       call allocate(surface_field, field%dim, surface_mesh, name="value")
       call insert_surface_field(field, trim(bc_name)//"_turbine", surface_field)
       call deallocate(surface_field)
       !call allocate(scalar_surface_field, surface_mesh, name="value")
       !call insert_surface_field(field, trim(bc_name)//"_flux", scalar_surface_field)
       !call deallocate(scalar_surface_field)
    end subroutine insert_flux_turbine_boundary_condition
  end subroutine populate_flux_turbine_boundary_conditions

  
  subroutine impose_reference_pressure_node(cmc_m, rhs, option_path)
    !!< If there are only Neumann boundaries on P, it is necessary to pin
    !!< the value of the pressure at one point. As the rhs of the equation
    !!< needs to be zeroed for this node, you will have to call this for
    !!< both of the pressure equations.
    type(csr_matrix), intent(inout) :: cmc_m
    type(scalar_field), intent(inout):: rhs
    character(len=*), intent(in) :: option_path

    integer :: reference_node, stat, stat2

    call get_option(trim(complete_field_path(option_path, stat2))//&
        &"/reference_node", reference_node, &
        & stat=stat)
    if (stat==0) then
       ! all processors now have to call this routine, although only
       ! process 1 sets it
       ewrite(1,*) 'Imposing_reference_pressure_node'    
       call set_reference_node(cmc_m, reference_node, rhs)
    end if

  end subroutine impose_reference_pressure_node

  subroutine impose_reference_velocity_node(big_m, rhs, option_path)
    !!< If solving the Stokes equation and there 
    !!< are only Neumann boundaries on u, it is necessary to pin
    !!< the value of the velocity at one point.
    !!< This is currently done using a big spring (unlike for pressure).
    type(petsc_csr_matrix), intent(inout) :: big_m
    type(vector_field), intent(inout):: rhs
    character(len=*), intent(in) :: option_path

    integer :: reference_node, stat, stat2

    call get_option(trim(complete_field_path(option_path, stat2))//&
        &"/reference_node", reference_node, &
        & stat=stat)
    if (stat==0) then
       ! all processors now have to call this routine, although only
       ! process 1 sets it
       ewrite(1,*) 'Imposing_reference_velocity_node'    
       call set_reference_node(big_m, reference_node, rhs)
    end if

  end subroutine impose_reference_velocity_node
  
end module boundary_conditions_from_options
