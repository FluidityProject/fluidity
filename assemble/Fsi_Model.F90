!    Copyright (C) 2006 Imperial College London and others.
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
#define INLINE_MATMUL

module fsi_model
! these 5 need to be on top and in this order, 
! so as not to confuse silly old intel compiler 
  use quadrature
  use elements
  use sparse_tools
  use fields
  use state_module
!
  use vtk_interfaces
  use linked_lists
  use intersection_finder_module
  use tetrahedron_intersection_module
  use unify_meshes_module
  use unittest_tools
  use global_parameters, only: FIELD_NAME_LEN, OPTION_PATH_LEN, &
       PYTHON_FUNC_LEN, dt, timestep, current_time
  use spud
  use timeloop_utilities
  use fefields, only: compute_lumped_mass, project_field
  use parallel_tools
  use diagnostic_variables
  use qmesh_module
  use mesh_files
  use read_triangle
  use fields_manipulation
  use solvers
  use pickers_inquire
  use transform_elements
  use field_derivatives
  use FLDebug
  use supermesh_construction
  use futils
  use meshdiagnostics
  use sparsity_patterns
  use vector_tools
  use tensors
  use fetools
  use interpolation_module
  use adjacency_lists
  use sparse_matrices_fields
  use bound_field_module
  use halos
!  use diagnostic_fields
  use boundary_conditions
  use data_structures
  use edge_length_module
  use conservative_interpolation_module
  use state_fields_module
  use fsi_projections

  implicit none

  logical, save :: adapt_at_previous_dt
!  logical, save :: have_prescribed_solid_movement
  type(state_type), save :: global_fluid_state, global_solid_state

  private

  public :: fsi_modelling, fsi_model_compute_diagnostics, &
             fsi_model_nonlinear_iteration_converged, &
             fsi_model_register_diagnostic, &
             fsi_model_pre_adapt_cleanup, fsi_post_adapt_operations, &
             fsi_ibm_all_alpha_projections, &
             fsi_add_dalpha_solid_dt

  contains

    subroutine fsi_modelling(state, solid_states, its)
    !! Main routine for FSI, being called every picard iteration
        type(state_type), intent(inout) :: state
        type(state_type), intent(inout), dimension(:) :: solid_states
        integer, intent(in) :: its

        ewrite(2, *) "inside fsi_modelling"

        ! If this is not the first picard iteration, set solid specific fiels on the fluid mesh
        ! to be as the last iterated field:
        if (its > 1) then
            call fsi_set_from_iterated_fields(state)
        end if

        ! 1-WAY COUPLING (prescribed solid velocity)
        ! First check if prescribed solid movement is enabled,
        ! and if so, move the solid mesh and update the new
        ! solid volume fractions as well
        if (its == 1) then
            call fsi_apply_prescribed_solid_velocity(state, solid_states)
        end if
        ! Some tests for fsi_interface
        ! call set_fsi_interface(state, solid_states)
        ! call set_fsi_interface_correction(state)
        ! call set_fsi_interface_from_alpha(state)

        ! Set absorption term/sigma
        call compute_fluid_absorption(state)
        ! Always set source term:
        call compute_source_term(state)
        ! Set the fluid velocity (only) field:
        call set_fsi_fluidvelocity(state)

        ewrite(2, *) 'leaving fsi_modelling'

    end subroutine fsi_modelling

    !----------------------------------------------------------------------------

    subroutine fsi_set_from_iterated_fields(state)
    !! This subroutine sets the solid related fields to be as the 
    !! previous iterated fields
        type(state_type), intent(inout) :: state
        type(scalar_field), pointer :: salpha, iter_salpha
        ! type(scalar_field), pointer :: sinterface, iter_sinterface
        type(vector_field), pointer :: svel, iter_svel
        type(vector_field), pointer :: sforce, iter_sforce
        character(len=OPTION_PATH_LEN) :: mesh_name
        integer :: i, num_solid_mesh

        ewrite(2,*) "Inside fsi_set_from_iterated_fields"

        ! Get number of solid meshes:
        num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

        ! Loop over solid meshes, to set solid volume fraction, force and velocity based on iterated
        ! values:
        solid_mesh_loop: do i=0, num_solid_mesh-1

            ! Get mesh name:
            call get_option('/embedded_models/fsi_model/geometry/mesh['//int2str(i)//']/name', mesh_name)

            ! Get relevant fields:
            salpha => extract_scalar_field(state, trim(mesh_name)//'SolidConcentration')
            sforce => extract_vector_field(state, trim(mesh_name)//'SolidForce')
            svel => extract_vector_field(state, trim(mesh_name)//'SolidVelocity')
            ! And their iterated counterparts:
            iter_salpha => extract_scalar_field(state, 'Iterated'//trim(mesh_name)//'SolidConcentration')
            iter_sforce => extract_vector_field(state, 'Iterated'//trim(mesh_name)//'SolidForce')
            iter_svel => extract_vector_field(state, 'Iterated'//trim(mesh_name)//'SolidVelocity')

            ! Now set the fields to equal as the iterated fields (these fields do not change after the first picard iteration):
            call set(salpha, iter_salpha)
            call set(sforce, iter_sforce)
            call set(svel, iter_svel)

        end do solid_mesh_loop

        ! Now do the same for the global fields (on the fluid mesh of course):
        salpha => extract_scalar_field(state, 'SolidConcentration')
        ! sinterface => extract_scalar_field(state, 'SolidPhase') ! SolidPhase only exists as a global field!
        sforce => extract_vector_field(state, 'SolidForce')
        svel => extract_vector_field(state, 'SolidVelocity')
        iter_salpha => extract_scalar_field(state, 'IteratedSolidConcentration')
        ! iter_sinterface => extract_scalar_field(state, 'IteratedSolidPhase') ! SolidPhase only exists as a global field!
        iter_sforce => extract_vector_field(state, 'IteratedSolidForce')
        iter_svel => extract_vector_field(state, 'IteratedSolidVelocity')
        call set(salpha, iter_salpha)
        ! call set(sinterface, iter_sinterface)
        call set(sforce, iter_sforce)
        call set(svel, iter_svel)

        ewrite(2,*) "Leaving fsi_set_from_iterated_fields"

    end subroutine fsi_set_from_iterated_fields

    !----------------------------------------------------------------------------

    function fsi_recompute_alpha(its, solid_moved) result(recompute)
    !! Function to set a logical that determines whether alpha is recomputed or not 
    !! at the next timestep
        integer, intent(in) :: its
        logical, optional, intent(in) :: solid_moved
        logical :: recompute

        recompute = .false.

        if (its == 1) then
            if (.not. present(solid_moved)) then
!                solid_moved = .true.
                if (adapt_at_previous_dt) then
                    recompute = .true.
                else
                    recompute = .false.
                end if
            else
                if (solid_moved .or. adapt_at_previous_dt) then
                    recompute = .true.
                else
                    recompute = .false.
                end if
            end if
        else ! its /= 1
            recompute = .false.
        end if

    end function fsi_recompute_alpha

    !----------------------------------------------------------------------------

    subroutine fsi_apply_prescribed_solid_velocity(state, solid_states)
    !! Subroutine that loops over all solid meshes that have a prescribed
    !! velocity or movement, applies solid movement and 
    !! updates the solid volume fraction fields afterwards
        type(state_type), intent(inout) :: state
        type(state_type), intent(inout), dimension(:) :: solid_states

        type(scalar_field), pointer :: solid_alpha_mesh, solid_alpha_global
        type(vector_field), pointer :: solid_velocity_mesh, solid_velocity_global

        character(len=OPTION_PATH_LEN) :: mesh_name
        character(len=OPTION_PATH_LEN) :: fsi_path="/embedded_models/fsi_model/"
        integer :: i, num_solid_mesh

        logical :: solid_moved

        ewrite(2,*) "Inside fsi_apply_prescribed_solid_velocity"

        ! Set logical to false, as at this moment, no solid has moved yet:
        solid_moved = .false.

        ! Get number of solid meshes:
        num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

        ! Loop over number of solid meshes defined:
        solid_mesh_position_loop: do i=0, num_solid_mesh-1

            ! Get mesh name:
            call get_option('/embedded_models/fsi_model/geometry/mesh['//int2str(i)//']/name', mesh_name)

            ! Check if the current mesh has a prescribed velocity/movement:
            if (have_option(trim(fsi_path)//"solid_phase::"//trim(mesh_name)//"/vector_field::SolidVelocity/prescribed")) then
                ! Move the mesh:
                call fsi_move_solid_mesh(state, solid_states(i+1), mesh_name, solid_moved)

                if (solid_moved) then
                  ewrite(1,*) "moved solid mesh: ", trim(mesh_name)
                else
                  ewrite(1,*) "solid mesh "//trim(mesh_name)//" did not move this timestep, although it has a prescribed velocity/position/movement."
                end if
            end if

            ! Recompute solid volume fraction, if solid just moved:
            if (solid_moved) then
                ewrite(1,*) "update alpha of solid mesh: ", trim(mesh_name)
                call fsi_ibm_projections_mesh(state, solid_states, mesh_name, project_solid_velocity=.true.)
            end if

            solid_moved = .false.

        end do solid_mesh_position_loop

        ! Update the global solid volume fraction:
        solid_alpha_global => extract_scalar_field(state, 'SolidConcentration')
        call zero(solid_alpha_global)
        ! And the solid velocity field:
        solid_velocity_global => extract_vector_field(state, 'SolidVelocity')
        call zero(solid_velocity_global)

        ! Loop over solid meshes again, to update global solid volume fraction field:
        solid_mesh_addto_global_fields_loop: do i=0, num_solid_mesh-1

            ! Get mesh name:
            call get_option('/embedded_models/fsi_model/geometry/mesh['//int2str(i)//']/name', mesh_name)

            ! Get solid volume fraction field of this mesh:
            solid_alpha_mesh => extract_scalar_field(state, trim(mesh_name)//'SolidConcentration')
            ! And the solid velocity:
            solid_velocity_mesh => extract_vector_field(state, trim(mesh_name)//'SolidVelocity')

            ! Add mesh solid volume fraction to the global solid volume fraction:
            call addto(solid_alpha_global, solid_alpha_mesh)
            ! And for the solid velocity:
            call addto(solid_velocity_global, solid_velocity_mesh)

        end do solid_mesh_addto_global_fields_loop

        ! Make sure that the solid volume fraction is max 1.0,
        ! only necessary if more than one solid mesh is present
        if (num_solid_mesh .gt. 1) then
            call sum_union_solid_volume_fraction(state)
        end if

        ! Set iterated global fields:
        call fsi_set_iterated_global_fields(state)

        ewrite(2,*) "Leaving fsi_apply_prescribed_solid_velocity"

    end subroutine fsi_apply_prescribed_solid_velocity

    !----------------------------------------------------------------------------

    subroutine fsi_set_iterated_global_fields(state)
      type(state_type), intent(inout) :: state
      type(scalar_field), pointer :: alpha_global, iter_alpha_global
      type(vector_field), pointer :: svel_global, iter_svel_global

        ewrite(2,*) "In fsi_set_iterated_global_fields"

        ! Get global alpha fields from state:
        alpha_global => extract_scalar_field(state, "SolidConcentration")
        iter_alpha_global => extract_scalar_field(state, "IteratedSolidConcentration")
        call zero(iter_alpha_global)
        call set(iter_alpha_global, alpha_global)

        svel_global => extract_vector_field(state, "SolidVelocity")
        iter_svel_global => extract_vector_field(state, "IteratedSolidVelocity")
        call zero(iter_svel_global)
        call set(iter_svel_global, svel_global)

    end subroutine fsi_set_iterated_global_fields

    !----------------------------------------------------------------------------

    subroutine fsi_ibm_all_alpha_projections(states, solid_states)
    !! Subroutine that loops over all solid meshes and calls the corresponding 
    !! subroutines to obtain the solid volume fraction
        type(state_type), intent(inout), dimension(:) :: states
        type(state_type), intent(inout), dimension(:) :: solid_states
        type(scalar_field), pointer :: alpha_global
        type(scalar_field), pointer :: solid_alpha_mesh

        character(len=OPTION_PATH_LEN) :: mesh_name
        integer :: num_solid_mesh
        integer :: i

        ewrite(2,*) "inside fsi_ibm_all_alpha_projections"

        ! Get relevant fields:
        alpha_global => extract_scalar_field(states(1), "SolidConcentration")
        ! Since we are (re)computing the solid volume fraction field, set the 'old' field to zero:
        call zero(alpha_global)

        ! Get number of solid meshes:
        num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

        ! Loop over number of solid meshes defined:
        solid_mesh_loop: do i=0, num_solid_mesh-1

            ! Get mesh name:
            call get_option('/embedded_models/fsi_model/geometry/mesh['//int2str(i)//']/name', mesh_name)

            ! Now do the projection to obtain the alpha field of the current solid mesh:
            call fsi_ibm_projections_mesh(states(1), solid_states, mesh_name)

            ! And add the alpha of the current solid to the global alpha (of all solids):
            solid_alpha_mesh => extract_scalar_field(states(1), trim(mesh_name)//'SolidConcentration')
            call addto(alpha_global, solid_alpha_mesh)

        end do solid_mesh_loop

        ! Make sure that the solid volume fraction is max 1.0,
        ! only necessary if more than one solid mesh is present
        if (num_solid_mesh .gt. 1) then
            call sum_union_solid_volume_fraction(states(1))
        end if

        ! At this point, the new solid volume fraction of all given solid meshes has been (re)computed and stored in SolidConcentration in state!

        ewrite(2,*) "leaving fsi_ibm_all_alpha_projections"

    end subroutine fsi_ibm_all_alpha_projections

    !----------------------------------------------------------------------------

    subroutine fsi_ibm_projections_mesh(state, solid_states, mesh_name, project_solid_velocity)
        type(state_type), intent(inout) :: state
        type(state_type), intent(inout), dimension(:) :: solid_states
        character(len=OPTION_PATH_LEN), intent(in) :: mesh_name
        logical, intent(in), optional :: project_solid_velocity

        type(vector_field), pointer :: fluid_position, fluid_velocity
        type(scalar_field), pointer :: alpha_solid_fluidmesh
        !type(scalar_field), pointer :: alpha_solid_solidmesh
        type(vector_field) :: solid_position_mesh ! these live on the SOLID mesh
        type(vector_field), pointer :: solid_velocity_mesh ! these live on the SOLID mesh
        type(vector_field), pointer :: solid_velocity_fluidmesh ! these live on the FLUID mesh
        type(scalar_field) :: alpha_tmp ! these live on the FLUID mesh
        type(vector_field) :: solid_velocity_fluidmesh_tmp ! these live on the FLUID mesh

        character(len=OPTION_PATH_LEN) :: proj_path, proj_type='no_interpolation'
        character(len=OPTION_PATH_LEN) :: proj_solver_path

        ewrite(2,*) "inside fsi_ibm_projections_mesh"

        ! Assemble path to the projection settings in the schema:
        proj_path = '/embedded_models/fsi_model/solid_phase::'//trim(mesh_name)//'/fsi_projection'

        ! Get relevant fields:
        fluid_velocity => extract_vector_field(state, "Velocity")
        fluid_position => extract_vector_field(state, "Coordinate")

        ! Get coordinate field of the current solid (mesh):
        solid_position_mesh = extract_vector_field(solid_states, trim(mesh_name)//'SolidCoordinate')
        ! Also the solid volume fraction field (on fluid mesh):
        alpha_solid_fluidmesh => extract_scalar_field(state, trim(mesh_name)//'SolidConcentration')
        ! And the solid volume fraction field (on solid mesh):
!        alpha_solid_solidmesh => extract_scalar_field(solid_states, trim(mesh_name)//'SolidConcentration')
        ! And the solid velocity field on the fluid mesh:
        solid_velocity_fluidmesh => extract_vector_field(state, trim(mesh_name)//'SolidVelocity')

        ! Allocate temp solid volume fraction field for the projection:
        call allocate(alpha_tmp, alpha_solid_fluidmesh%mesh, 'TMPSolidConcentration')
        call zero(alpha_tmp)

        ! Distinguish between different projection methods:
        if (have_option(trim(proj_path)//'/galerkin_projection')) then
            proj_type = 'galerkin_projection'
        else if (have_option(trim(proj_path)//'/grandy_interpolation')) then
            proj_type = 'grandy_interpolation'
        else if (have_option(trim(proj_path)//'/consistent_interpolation')) then
            proj_type = 'consistent_interpolation'
        end if


        ! Galerkin projection via supermesh:
        select case (proj_type)
            case ('galerkin_projection')
                proj_solver_path = trim(proj_path)//'/'//trim(proj_type)//'/continuous'
                if (present_and_true(project_solid_velocity)) then
                    ! Get pointer to the solid velocity
                    solid_velocity_mesh => extract_vector_field(solid_states, trim(mesh_name)//"SolidVelocity")
                    ! Allocate temp solid velocity field (on fluid mesh) for the projection:
                    call allocate(solid_velocity_fluidmesh_tmp, solid_velocity_fluidmesh%dim, solid_velocity_fluidmesh%mesh, 'TMP'//trim(mesh_name)//'SolidVelocity')
                    call zero(solid_velocity_fluidmesh_tmp)
                    ! Computing alpha_s^f by projecting unity and the u_s^f by projection u_s^s from the solid mesh to the fluid mesh:
                    call fsi_one_way_galerkin_projection(fluid_velocity, fluid_position, solid_position_mesh, &
                              & alpha_tmp, proj_solver_path, &
                              & solid_velocity_mesh, solid_velocity_fluidmesh_tmp)
                else
                    ! Computing /alpha_s^f by projecting unity from the solid mesh to the fluid mesh:
                    call fsi_one_way_galerkin_projection(fluid_velocity, fluid_position, solid_position_mesh, &
                              & alpha_tmp, proj_solver_path)
                end if

            case ('grandy_interpolation')
                ! Grandy interpolation:
                call fsi_one_way_grandy_interpolation(fluid_position, solid_position_mesh, alpha_tmp)

            ! Add more interpolations here

            case default
                FLAbort("Unrecognised interpolation algorithm")
        end select


        ! After projection, set the solid volume fraction of the current solid:
        call set(alpha_solid_fluidmesh, alpha_tmp)
        ! And solid velocity (on fluid mesh), IFF it was projected (thus IFF d(s_solid)/dt /= 0 or in other words u_solid /= 0)
        if (present_and_true(project_solid_velocity)) then
            call set(solid_velocity_fluidmesh, solid_velocity_fluidmesh_tmp)
        else
            call zero(solid_velocity_fluidmesh)
        end if

        ! Deallocate temp arrays:
        call deallocate(alpha_tmp)
        if (present_and_true(project_solid_velocity)) then
            call deallocate(solid_velocity_fluidmesh_tmp)
        end if

        ewrite(2,*) "leaving fsi_ibm_projections_mesh"

    end subroutine fsi_ibm_projections_mesh

    !----------------------------------------------------------------------------

    subroutine sum_union_solid_volume_fraction(state)
    !! In case two solid meshes are in close neighbourhood, the global solid
    !! volume fraction field could be significantly > 1.0, which is unphysical, 
    !! but could results from the bounded, minimal diffusive Galerkin projection
    !! via supermesh. 
    !! This subroutine ensures that in such a case, the solid volume fraction is 
    !! set to 1.0 at all nodes >1.0. 
    !! As a result, we loose conservation, but the loss should be neglectably
    !! small, and furthermore is expected to be within a solid.

        type(state_type), intent(inout) :: state
        type(scalar_field), pointer :: alpha_global
        integer, dimension(:), pointer :: nodes
        integer :: i, ele

        ! Get global alpha field:
        alpha_global => extract_scalar_field(state, "SolidConcentration")

        ! Looping over ele, looping over nodes:
        do ele = 1, ele_count(alpha_global%mesh)
            nodes => ele_nodes(alpha_global%mesh, ele)
            do i = 1, size(nodes)
                if (node_val(alpha_global, nodes(i)) .gt. 1.0) then
                    call set(alpha_global, nodes(i), 1.0)
                end if
            end do
        end do

    end subroutine sum_union_solid_volume_fraction

    !----------------------------------------------------------------------------

    subroutine set_fsi_interface_from_alpha(state)
    !! Set the fsi interface based on alpha values:
        type(state_type), intent(inout) :: state

        type(scalar_field), pointer :: alpha
        type(scalar_field), pointer :: fsi_interface, fsi_interface_iter
        integer, dimension(:), pointer :: nodes
        integer :: i, j, ele

        alpha => extract_scalar_field(state, 'IteratedSolidConcentration')
        fsi_interface_iter => extract_scalar_field(state, 'IteratedSolidPhase')
        call zero(fsi_interface_iter)

        do ele=1,ele_count(alpha%mesh)
            nodes => ele_nodes(alpha%mesh, ele)
            do i = 1, size(nodes)
                if (node_val(alpha, nodes(i)) .gt. 0.0 .and. node_val(alpha, nodes(i)) .lt. 1) then
                    call set(fsi_interface_iter, nodes(i), 1.0)
                end if
            end do
        end do

        fsi_interface => extract_scalar_field(state, 'SolidPhase')
        call zero(fsi_interface)
        call set(fsi_interface, fsi_interface_iter)

    end subroutine set_fsi_interface_from_alpha

    !----------------------------------------------------------------------------

    subroutine set_fsi_interface_correction(state)
    !! Set the fsi interface based on alpha values:
        type(state_type), intent(inout) :: state

        type(scalar_field), pointer :: alpha
        type(scalar_field), pointer :: fsi_interface, fsi_interface_iter
        integer, dimension(:), pointer :: nodes
        integer :: i, j, ele

        alpha => extract_scalar_field(state, 'IteratedSolidConcentration')
        fsi_interface => extract_scalar_field(state, 'SolidPhase')
        fsi_interface_iter => extract_scalar_field(state, 'IteratedSolidPhase')
        call zero(fsi_interface_iter)

        do ele=1,ele_count(alpha%mesh)
            nodes => ele_nodes(alpha%mesh, ele)
            do i = 1, size(nodes)
                if (node_val(fsi_interface, nodes(i)) == 1.0) then
                    if (node_val(alpha, nodes(i)) .gt. 0.0 .and. node_val(alpha, nodes(i)) .lt. 1) then
                        call set(fsi_interface_iter, nodes(i), 1.0)
                    end if
                end if
            end do
        end do

        fsi_interface => extract_scalar_field(state, 'SolidPhase')
        call zero(fsi_interface)
        call set(fsi_interface, fsi_interface_iter)

    end subroutine set_fsi_interface_correction

    !----------------------------------------------------------------------------

    subroutine set_fsi_interface(state, solid_states)
    !! This subroutine sets a scalar field which is 1 only at the fluid elements
    !! which intersect with the solid boundary elements
        type(state_type), intent(inout) :: state
        type(state_type), intent(inout), dimension(:) :: solid_states

        type(vector_field), pointer :: fluid_position, solid_position
        
        type(scalar_field) :: fsi_interface_mesh
        type(scalar_field), pointer :: fsi_interface, fsi_interface_iter

        character(len=OPTION_PATH_LEN) :: mesh_path, mesh_name
        integer :: i, num_solid_mesh

        ! Get fluid coordinate mesh:
        fluid_position => extract_vector_field(state, 'Coordinate')
        ! Get the global interface scalar field (on the fluid mesh):
        fsi_interface => extract_scalar_field(state, 'SolidPhase')
        fsi_interface_iter => extract_scalar_field(state, 'IteratedSolidPhase')
        call zero(fsi_interface)
        call zero(fsi_interface_iter)

        num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

        ! Loop over number of solid meshes defined:
        solid_mesh_loop: do i=0, num_solid_mesh-1
            ! Get mesh name:
            mesh_path="/embedded_models/fsi_model/geometry/mesh["//int2str(i)//"]"
            call get_option(trim(mesh_path)//'/name', mesh_name)
            ! Solid coordinate mesh:
            solid_position => extract_vector_field(solid_states, trim(mesh_name)//'SolidCoordinate')

            ! Allocate local interface scalar field:
            call allocate(fsi_interface_mesh, fsi_interface%mesh, 'TMPLocalSolidPhase')

            ! Get the interface for this solid mesh:
            ewrite(2,*) "Get FSI interface for solid '"//trim(mesh_name)//"'."
            call fsi_get_interface(fsi_interface_mesh, fluid_position, solid_position)

            ! Add the local interface field to the global field:
            call addto(fsi_interface_iter, fsi_interface_mesh)

            ! Deallocate:
            call deallocate(fsi_interface_mesh)

        end do solid_mesh_loop

        call set(fsi_interface, fsi_interface_iter)

    end subroutine set_fsi_interface

    !----------------------------------------------------------------------------

    subroutine fsi_compute_intersection_map(state, solid_states, solid_position, map_SF)
    !! This subroutine computes a list of intersecting elements between a solid and fluid mesh
    !! The rtree intersection finder is used as it is more robust in parallel than the 
    !! advancing front intersection finder, e.g. it doesn't abort if no intersection was found,
    !! which is possible with the solid domain being immersed in a fluid domain
        type(state_type), intent(inout) :: state
        type(state_type), intent(inout), dimension(:) :: solid_states
        type(vector_field), pointer, intent(in) :: solid_position
        type(ilist), dimension(:), allocatable, intent(inout) :: map_SF
        type(ilist) :: map_SF_tmp

        type(vector_field), pointer :: fluid_position
        integer :: ele_F, ele_S
        integer :: num_intersections
        integer :: i, dim

        ewrite(2,*) "inside fsi_compute_intersection_map"

        ! Get fluid positions:
        fluid_position => extract_vector_field(state, "Coordinate")

        ! Set the dimension for the intersection finder:
        dim = mesh_dim(solid_position)
        call intersector_set_dimension(dim) 
        ! Use rtree intersection finder because it is more robust, e.g.
        ! it does not fail when running in parallel and no solid element is found in the 
        ! composition of the fluid mesh:
        ! Set input for rtree intersection finder:
        call rtree_intersection_finder_set_input(fluid_position)

        do ele_S=1,ele_count(solid_position)
            ! Via rtree, find intersection of solid element with fluid elements:
            call rtree_intersection_finder_find(solid_position, ele_S)
            ! Fetch output, the number of intersections
            call rtree_intersection_finder_query_output(num_intersections)
            ! Generate an ilist of elements in that intersect with ele_S:
            if (num_intersections .ne. 0) then
                do i=1, num_intersections
                    ! Get the donor (fluid) element which intersects with ele_S
                    call rtree_intersection_finder_get_output(ele_F, i)
                    call insert(map_SF_tmp, ele_F)
                end do
            end if
        end do

        ! Copy values over to map_SF
        map_SF = map_SF_tmp

        ewrite(2,*) "leaving fsi_compute_intersection_map"

    end subroutine fsi_compute_intersection_map

    !----------------------------------------------------------------------------

    subroutine fsi_remove_from_state(state)
    !! This subroutine removes additional fields from state, e.g. solidvolumefraction
    !! per solid mesh, solidforce per solid mesh etc, such that state only contains
    !! the fields stated in the flml, and thus the fluid mesh can be adapted

        type(state_type), intent(inout) :: state
        type(scalar_field), pointer :: alpha_solidmesh
        type(vector_field), pointer :: solidforce_mesh, solidvelocity_mesh
        character(len=OPTION_PATH_LEN) :: mesh_path, mesh_name
        integer :: i, num_solid_mesh

        ewrite(2,*) "inside fsi_remove_from_state"

        ! figure out if we want to print out diagnostics and initialise files
        ! check for mutiple solids and get translation coordinates
        num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

        ! Loop over number of solid meshes defined:
        solid_mesh_loop: do i=0, num_solid_mesh-1
           ! Get mesh name:
           mesh_path="/embedded_models/fsi_model/geometry/mesh["//int2str(i)//"]"
           call get_option(trim(mesh_path)//'/name', mesh_name)
           
           ! Change the option_path:
           alpha_solidmesh => extract_scalar_field(state, trim(mesh_name)//'SolidConcentration')
           alpha_solidmesh%option_path = 'testtmp/'//trim(alpha_solidmesh%option_path)
           solidforce_mesh => extract_vector_field(state, trim(mesh_name)//'SolidForce')
           solidforce_mesh%option_path = 'testtmp/'//trim(solidforce_mesh%option_path)
           solidvelocity_mesh => extract_vector_field(state, trim(mesh_name)//'SolidForce')
           solidvelocity_mesh%option_path = 'testtmp/'//trim(solidvelocity_mesh%option_path)

        end do solid_mesh_loop

        ewrite(2,*) "leaving fsi_remove_from_state"

    end subroutine fsi_remove_from_state

    !----------------------------------------------------------------------------

    subroutine fsi_insert_state(state)
    !! This subroutine inserts additional fields to state, e.g. solidvolumefraction
    !! per solid mesh, solidforce per solid mesh etc. 
    !! This should happen after adapting the mesh, such that the neccessary fields
    !! are allocated.

        type(state_type), intent(inout) :: state
        type(mesh_type), pointer :: solid_mesh
        type(vector_field), pointer :: solid_position, fluid_position, solid_force, solid_velocity
        type(scalar_field), pointer :: alpha_global
        type(scalar_field) :: alpha_solidmesh
        type(vector_field) :: solidforce_mesh, solidvelocity_mesh
        character(len=OPTION_PATH_LEN) :: mesh_path, mesh_name
        integer :: i, num_solid_mesh

        ewrite(2,*) "inside fsi_insert_state"

        if (adapt_at_previous_dt .and. .not.(timestep == 1)) then

          ! Get fields:
          fluid_position => extract_vector_field(state, 'Coordinate')
          alpha_global => extract_scalar_field(state, 'SolidConcentration')
          solid_force => extract_vector_field(state, "SolidForce")
          solid_velocity => extract_vector_field(state, "SolidVelocity")

          ! figure out if we want to print out diagnostics and initialise files
          ! check for mutiple solids and get translation coordinates
          num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

          ! Loop over number of solid meshes defined:
          solid_mesh_loop: do i=0, num_solid_mesh-1
             ! Get mesh name:
             mesh_path="/embedded_models/fsi_model/geometry/mesh["//int2str(i)//"]"
             call get_option(trim(mesh_path)//'/name', mesh_name)

             ! Extract solid mesh from original state, and 
             ! insert solid mesh into additional solid state:
             solid_mesh => extract_mesh(global_solid_state, trim(mesh_name))
             ! Add to extra state:
             call insert(state, solid_mesh, trim(solid_mesh%name))

             ! Solid positions:
             solid_position => extract_vector_field(global_solid_state, trim(solid_mesh%name)//'SolidCoordinate')
             ! Add to extra state:
             call insert(state, solid_position, trim(solid_position%name))

             ! SolidConcentration per solid mesh:
             call allocate(alpha_solidmesh, fluid_position%mesh, trim(solid_mesh%name)//"SolidConcentration")
             alpha_solidmesh%option_path = alpha_global%option_path
             call zero(alpha_solidmesh)
             call insert(state, alpha_solidmesh, trim(alpha_solidmesh%name))

             ! SolidForce per solid mesh:
             call allocate(solidforce_mesh, solid_force%dim, solid_force%mesh, trim(solid_mesh%name)//"SolidForce")
             solidforce_mesh%option_path = solid_force%option_path
             call zero(solidforce_mesh)
             call insert(state, solidforce_mesh, trim(solid_mesh%name)//"SolidForce")
             
             ! SolidVelocity per solid mesh:
             if (have_option("/embedded_models/fsi_model/solid_phase::"//trim(mesh_name)//"/vector_field::SolidVelocity/prescribed")) then
                call allocate(solidvelocity_mesh, solid_velocity%dim, solid_velocity%mesh, trim(solid_mesh%name)//"SolidVelocity")
                solidvelocity_mesh%option_path = solid_velocity%option_path
                call zero(solidforce_mesh)
                call insert(state, solidvelocity_mesh, trim(solid_mesh%name)//"SolidVelocity")
             end if

          end do solid_mesh_loop

        end if

        ewrite(2,*) "leaving fsi_insert_state"

    end subroutine fsi_insert_state

    !----------------------------------------------------------------------------

    subroutine fsi_assemble_fs_states(fluid_state, solid_state, solid_position, fluid_position)
    !! This subroutine assembles a solid and a fluid state, that can be used to add more fields to
    !! it. Currently it is not used, but their might be some useful need for it in the future
        type(state_type), dimension(:), allocatable, intent(inout) :: fluid_state, solid_state
        type(vector_field), pointer, intent(in) :: fluid_position, solid_position

        type(mesh_type), pointer :: solid_mesh, fluid_mesh

        ewrite(2,*) "inside fsi_assemble_fs_states"

        ! Get solid and fluid mesh and store them in the new states:
        solid_mesh => solid_position%mesh
        fluid_mesh => fluid_position%mesh

        ! Insert meshes fields into the new states:
        call insert(solid_state, solid_mesh, "Mesh")
        call insert(fluid_state, fluid_mesh, "Mesh")
        call insert(solid_state, solid_position, "Coordinate")
        call insert(fluid_state, fluid_position, "Coordinate")

        ewrite(2,*) "leaving fsi_assemble_fs_states"

    end subroutine fsi_assemble_fs_states

  !----------------------------------------------------------------------------

    subroutine compute_fluid_absorption(state)
    !! This subroutine computes the absorption of the fluid
    !! within a solid, also called the time relaxation of the fluid
        type(state_type), intent(inout) :: state
        type(vector_field), pointer :: absorption_iter, absorption
        type(scalar_field), pointer :: alpha_global
        integer :: i, j, ele
        real :: sigma, beta
        integer, dimension(:), pointer :: nodes

        ewrite(2, *) 'inside compute_fluid_absorption'

        alpha_global => extract_scalar_field(state, 'IteratedSolidConcentration')
        ! Get absorption velocity from state:
        absorption_iter => extract_vector_field(state, 'IteratedVelocityAbsorption')
        call zero(absorption_iter)

        do ele = 1, ele_count(alpha_global%mesh)
            nodes => ele_nodes(alpha_global%mesh, ele)
            do i = 1, size(nodes)
                sigma = node_val(alpha_global, nodes(i)) / dt
                do j = 1, mesh_dim(absorption_iter%mesh)
                    call set(absorption_iter, j, nodes(i), sigma)
                end do
            end do
        end do

        if (have_option('/embedded_models/fsi_model/beta')) then
            call get_option('/embedded_models/fsi_model/beta', beta)
            call scale(absorption_iter, beta)
        end if

        ! Now set the absorption velocity:
        ! Get absorption velocity from state:
        absorption => extract_vector_field(state, 'VelocityAbsorption')
        call set(absorption, absorption_iter)

        ewrite(2, *) "leaving compute_fluid_absorption"

    end subroutine compute_fluid_absorption

  !----------------------------------------------------------------------------

    subroutine fsi_move_solid_mesh(state, solid_state, mesh_name, solid_moved)
    !! Moving a solid mesh
      type(state_type), intent(in) :: state
      type(state_type), intent(inout) :: solid_state
      character(len=OPTION_PATH_LEN), intent(in) :: mesh_name
      logical, intent(inout) :: solid_moved

      type(vector_field), pointer :: solid_velocity_global
      type(vector_field), pointer :: solid_position_mesh, new_solid_position_mesh, solid_velocity_mesh
      type(vector_field) :: solid_python_return_field, solid_movement_mesh

      character(len=PYTHON_FUNC_LEN) :: func
      character(len=OPTION_PATH_LEN) :: fsi_path="/embedded_models/fsi_model/"
      integer :: d

      ewrite(2,*) "inside move_solid_mesh"

      ! 1st get coordinate field of this solid mesh, and its velocity field (which is on the fluid mesh):
      solid_position_mesh => extract_vector_field(solid_state, trim(mesh_name)//"SolidCoordinate")
      ! Get solid velocity field and reset it:
      solid_velocity_mesh => extract_vector_field(solid_state, trim(mesh_name)//'SolidVelocity')
      call zero(solid_velocity_mesh)
      ! and create new field in which the return field from the python interface is stored:
      call allocate(solid_python_return_field, solid_position_mesh%dim, solid_position_mesh%mesh, name=trim(mesh_name)//"SolidPythonReturnField")
      call zero(solid_python_return_field)
      ! Allocate field for the travelled distance of a solid:
      call allocate(solid_movement_mesh, solid_position_mesh%dim, solid_position_mesh%mesh, name=trim(mesh_name)//"SolidMovement")
      call zero(solid_movement_mesh)


      ! 2nd get the python function for this solid mesh:
      ! Here we have to differentiate between the user using a prescribed 'SolidVelocity' and 'SolidMovement':
      if (have_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_velocity')) then
        call get_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_velocity/python', func)
      else if (have_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_movement')) then
        call get_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_movement/python', func)
      end if
      ! 3rd Set the new coordinates from the python function:
      ! Here we have to differentiate between 'SolidVelocity' and 'SolidMovement' again:
      ! 'SolidVelocity' gets back the solid velocity,
      ! and 'SolidMovement' gets back the travelled distance within dt.
      ! This allows to easily set prescribed movement for different kinds of problems.
      ! Also, catch what was specified by the user to pass into the Python interface, either current_time or timestep (dt):
      ! First, for prescribed SolidVelocity do:
      if (have_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_velocity')) then
        if (have_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_velocity/time_variable_in_python/current_time')) then
          call set_from_python_function(solid_python_return_field, func, solid_position_mesh, current_time)
        else if (have_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_velocity/time_variable_in_python/current_timestep')) then
          call set_from_python_function(solid_python_return_field, func, solid_position_mesh, dt)
        else
          FLAbort('Neither current_time nor current_timestep was enabled in the schema. This must not happen')
        end if
      ! For prescribed SolidMovement, do:
      else if (have_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_movement')) then
        if (have_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_movement/time_variable_in_python/current_time')) then
          call set_from_python_function(solid_python_return_field, func, solid_position_mesh, current_time)
        else if (have_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_movement/time_variable_in_python/current_timestep')) then
          call set_from_python_function(solid_python_return_field, func, solid_position_mesh, dt)
        else
          FLAbort('Neither current_time nor current_timestep was enabled in the schema. This must not happen')
        end if
      end if

      ! 4th Check if the mesh has actually moved:
      ! If a prescribed SolidVelocity or SolidMovement was set for this solid mesh, 
      ! the solid did NOT move, if the return field from the python interface is zero,
      ! if a prescribed SolidVelocity or SolidMovement was set, the solid did NOT move, 
      ! if the difference of old position and new position is zero:
      if (maxval(solid_python_return_field%val(:,:)) == 0.0 .and. &
          minval(solid_python_return_field%val(:,:)) == 0.0) then
        solid_moved = .false.
      else
        solid_moved = .true.
      end if


      ! Set fields (new coordinates, and velocity), IFF the solid moved:
      if (solid_moved) then
        ! 4th Set the new coordinates and velocity of the solid:
        ! First for prescribed SolidVelocity:
        if (have_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_velocity')) then
            ! Set solid velocity:
            call set(solid_velocity_mesh, solid_python_return_field) ! Here solid_python_return_field is the solid velocity
            ! Compute the distance travelled within dt:
            call addto(solid_movement_mesh, solid_python_return_field, dt) ! Here solid_python_return_field is the solid velocity
            ! Set new coordinates:
            call addto(solid_position_mesh, solid_movement_mesh)
        ! Then for SolidMovement:
        else if (have_option(trim(fsi_path)//'solid_phase::'//trim(mesh_name)//'/vector_field::SolidVelocity/prescribed/python_movement')) then
            ! Set solid velocity:
            call addto(solid_velocity_mesh, solid_python_return_field, 1.0/dt) ! Here solid_python_return_field is the travelled distance
            ! Set new coordinates:
            call addto(solid_position_mesh, solid_python_return_field) ! Here solid_python_return_field is the travelled distance
        end if
      end if

      ! 5th Deallocate:
      call deallocate(solid_python_return_field)
      call deallocate(solid_movement_mesh)

      ewrite(2,*) 'end of move_solid_mesh'

    end subroutine fsi_move_solid_mesh

    !----------------------------------------------------------------------------

    subroutine compute_source_term(state)

      type(state_type), intent(inout) :: state

      type(vector_field), pointer :: source_term_iter, source_term
      type(vector_field), pointer :: fluid_velocity, solid_velocity, fluid_absorption
      type(scalar_field), pointer :: alpha, fsi_interface

      integer, dimension(:), pointer :: nodes
      integer :: i, j, ele
      real :: beta
      
      ewrite(2,*) "Inside compute_source_term"

      fluid_velocity => extract_vector_field(state, "IteratedVelocity")
      !fluid_absorption => extract_vector_field(state, "IteratedVelocityAbsorption")
      solid_velocity => extract_vector_field(state, "IteratedSolidVelocity")
      alpha => extract_scalar_field(state, "IteratedSolidConcentration")
      ! fsi_interface => extract_scalar_field(state, 'IteratedSolidPhase')

      ! Setting Iterated source term:
      source_term_iter => extract_vector_field(state, "IteratedVelocitySource")
      call zero(source_term_iter)

      ! beta:
      if (have_option('/embedded_models/fsi_model/beta')) then
          call get_option('/embedded_models/fsi_model/beta', beta)
      else
          beta = 1
      end if

!      ! (\rho_f \alpha_s / \Delta t) (\hat{u_f} or u) - F_s/ \Delta t
!      do ele = 1, ele_count(source_term_iter%mesh)
!          nodes => ele_nodes(source_term_iter%mesh, ele)
!          do i = 1, size(nodes)
!              ! Originally:
!              ! The if statement below causes instability for the source term on a DG mesh:
!!              if (node_val(alpha, nodes(i)) .gt. 0.0 .and. node_val(alpha, nodes(i)) .lt. 1.0) then
!              ! Source term should be at the FSI interface, so cap short before inside the solid volume
!              ! here and for now by comparing against the value of alpha, must be in interval [0, 0.999]:
!              if (node_val(alpha, nodes(i)) .gt. 0.0 .and. node_val(alpha, nodes(i)) .lt. 0.9999) then
!          !    if (node_val(alpha, nodes(i)) .gt. 0.0) then
!                  do j = 1, source_term_iter%dim
!!                      ! Originally used, and theorhetically correct:
!!                      ! Originally used this source term: u_s*sigma - sigma*(alpha*u_f + u_s), with sigma=1/dt
!!                      call set(source_term_iter, j, nodes(i), (node_val(solid_velocity,j,nodes(i)) / dt) - &
!!                        & (node_val(solid_velocity,j,nodes(i)) + &
!!                        & node_val(fluid_velocity,j,nodes(i)) ) / dt )
!                      call set(source_term_iter, j, nodes(i), ( (node_val(solid_velocity,j,nodes(i)) + &
!                        & node_val(fluid_velocity,j,nodes(i)) ) / dt ) - &
!                        & (node_val(solid_velocity,j,nodes(i)) / dt) )
!                      ! Below was used before, which is the same as above but with alpha:
!                      !  & node_val(fluid_velocity,j,nodes(i)) * node_val(alpha,nodes(i))) / dt )
!
!                  end do
!              end if
!          end do
!      end do

      do ele = 1, ele_count(source_term_iter%mesh)
          nodes => ele_nodes(source_term_iter%mesh, ele)
          do i = 1, size(nodes)
              if (node_val(alpha, nodes(i)) .gt. 0.0) then
                  do j = 1, source_term_iter%dim
                      call set(source_term_iter, j, nodes(i), ( (node_val(solid_velocity,j,nodes(i))*(beta/dt) ) - &
                        & node_val(fluid_velocity,j,nodes(i)) ) )
                  end do
              end if
          end do
      end do

      ! Setting source term:
      source_term => extract_vector_field(state, "VelocitySource")
      call set(source_term, source_term_iter)

      !end if

    ewrite(2,*) "Leaving compute_source_term"

    end subroutine compute_source_term

    !----------------------------------------------------------------------------

    subroutine set_fsi_fluidvelocity(state)
      type(state_type), intent(inout) :: state

      type(vector_field), pointer :: fsifluidvelocity_iter, fsifluidvelocity
      type(vector_field), pointer :: bulk_velocity, solid_velocity


      ewrite(2,*) "Inside set_fsi_fluidvelocity"

      fsifluidvelocity_iter => extract_vector_field(state, 'IteratedFSIFluidVelocity')
      fsifluidvelocity => extract_vector_field(state, 'FSIFluidVelocity')
      bulk_velocity => extract_vector_field(state, 'IteratedVelocity')
      solid_velocity => extract_vector_field(state, 'IteratedSolidVelocity')

      call zero(fsifluidvelocity)
      call zero(fsifluidvelocity_iter)

      call set(fsifluidvelocity_iter, bulk_velocity)
      call addto(fsifluidvelocity_iter, solid_velocity, -1.0)
      call set(fsifluidvelocity, fsifluidvelocity_iter)

      ewrite(2,*) "Leaving set_fsi_fluidvelocity"
    end subroutine set_fsi_fluidvelocity

    !----------------------------------------------------------------------------

    subroutine fsi_add_dalpha_solid_dt(state, ct_rhs)

      type(state_type), intent(in) :: state
      type(scalar_field), intent(inout) :: ct_rhs

      type(scalar_field), pointer :: alpha, old_alpha
      type(scalar_field) :: alpha_projected, old_alpha_projected
      type(scalar_field) :: alpha_pmesh, old_alpha_pmesh
      type(scalar_field), pointer :: p
      type(vector_field), pointer :: x
      type(element_type), pointer :: p_shape 
      type(element_type) :: test_function
      integer :: ele

      ewrite(2,*) "Inside fsi_add_dalpha_solid_dt"

      alpha => extract_scalar_field(state, "IteratedSolidConcentration")
      old_alpha => extract_scalar_field(state, "OldSolidConcentration")
      x => extract_vector_field(state, "Coordinate")
      p => extract_scalar_field(state, "IteratedPressure")

      ! Now project alpha to the cg coordinate mesh, then remap to the pressure mesh:
      ! First for alpha:
      call allocate(alpha_projected, x%mesh, 'SolidConcentrationCoordinateMesh')
      call allocate(alpha_pmesh, p%mesh, 'SolidConcentrationPressureMesh')
      call zero(alpha_projected)
      call zero(alpha_pmesh)
      call project_field(alpha, alpha_projected, x)
      call remap_field(alpha_projected, alpha_pmesh)
      ! And for old solid position (old_alpha):
      call allocate(old_alpha_projected, x%mesh, 'OldSolidConcentrationCoordinateMesh')
      call allocate(old_alpha_pmesh, p%mesh, 'OldSolidConcentrationPressureMesh')
      call zero(old_alpha_projected)
      call zero(old_alpha_pmesh)
      call project_field(old_alpha, old_alpha_projected, x)
      call remap_field(old_alpha_projected, old_alpha_pmesh)

      ewrite_minmax(alpha)
      ewrite_minmax(old_alpha)
      ewrite_minmax(alpha_pmesh)
      ewrite_minmax(old_alpha_pmesh)
      ewrite_minmax(p)
      ewrite_minmax(ct_rhs)

      do ele = 1, element_count(p)
         p_shape => ele_shape(p, ele)
         test_function = p_shape

         call add_ct_rhs_element_cg(ele, test_function, &
              p_shape, x, p, alpha_pmesh, old_alpha_pmesh, ct_rhs)
      end do

      ewrite_minmax(ct_rhs)

      call deallocate(alpha_projected)
      call deallocate(old_alpha_projected)
      call deallocate(alpha_pmesh)
      call deallocate(old_alpha_pmesh)

    contains

      subroutine add_ct_rhs_element_cg(ele, test_function, &
           p_shape, x, p, alpha, old_alpha, ct_rhs)

        type(scalar_field), intent(inout) :: ct_rhs
        integer, intent(in) :: ele
        type(element_type), intent(in) :: test_function, p_shape
        type(vector_field), intent(in) :: x
        type(scalar_field), intent(in) :: alpha, old_alpha, p

        real, dimension(ele_ngi(p, ele)) :: detwei
        real, dimension(ele_loc(p, ele), ele_ngi(p, ele), x%dim) :: dp_t
        integer, dimension(:), pointer :: p_ele
        real, dimension(ele_loc(p, ele)) :: mat
        real, dimension(ele_ngi(alpha, ele)) :: ds

        call transform_to_physical(x, ele, &
             p_shape, dshape=dp_t, detwei=detwei)

        ds = ele_val_at_quad(alpha, ele) - ele_val_at_quad(old_alpha, ele)
        mat = shape_rhs(test_function, detwei*ds/dt)

        p_ele => ele_nodes(p, ele)

        call addto(ct_rhs, p_ele, mat)

      end subroutine add_ct_rhs_element_cg

    end subroutine fsi_add_dalpha_solid_dt

    !----------------------------------------------------------------------------

    subroutine fsi_model_solid_force_computation(state, solid_force_diag, mesh_name)
    ! Computes the solid force field of the solid with name 'mesh_name'
      type(state_type), intent(inout) :: state
      real, dimension(:), intent(inout) :: solid_force_diag
      character(len=OPTION_PATH_LEN), optional ::  mesh_name

      type(vector_field), pointer :: fluid_velocity, fluid_coord, fluid_absorption
      type(vector_field), pointer :: solidforce, solidforce_global
      type(vector_field), pointer :: solid_velocity
      type(scalar_field), pointer :: alpha
      integer :: i, j, ele
      integer, dimension(:), pointer :: nodes
      real :: beta

      ewrite(2, *) "inside fsi_model_solid_force_computation"

      ! beta:
      if (have_option('/embedded_models/fsi_model/beta')) then
          call get_option('/embedded_models/fsi_model/beta', beta)
      else
          beta = 1
      end if

      ! Get relevant quantities from state
      fluid_velocity => extract_vector_field(state, "Velocity")
      fluid_absorption => extract_vector_field(state, "VelocityAbsorption")
      fluid_coord => extract_vector_field(state, "Coordinate")
      if(present(mesh_name)) then
          solidforce_global => extract_vector_field(state, "SolidForce")
          solidforce => extract_vector_field(state, trim(mesh_name)//"SolidForce")
          alpha => extract_scalar_field(state, trim(mesh_name)//"SolidConcentration")
      else
          solidforce => extract_vector_field(state, "SolidForce")
          alpha => extract_scalar_field(state, "SolidConcentration")
      end if
      ! solidvelocity doesn't have to be of that specific solid, but can be the global solidvelocity field,
      ! as for the computation of the solidforce (for a specific solid), only the elements where alpha>0.0
      ! are taken into account. Since alpha is solid-mesh specific, this works!
      solid_velocity => extract_vector_field(state, trim(mesh_name)//'SolidVelocity')

      ! 1-way coupling:
!      if (have_option('/embedded_models/fsi_model/one_way_coupling') &
!           & .and. (.not. (have_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed') .or. &
!           & have_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidMovement/prescribed') )) ) then
!
!          ! Looping over ele, looping over nodes:
!          call zero(solidforce)
!          do ele = 1, ele_count(solidforce%mesh)
!              nodes => ele_nodes(solidforce%mesh, ele)
!              do i = 1, size(nodes)
!                  if (node_val(alpha, nodes(i)) .gt. 0.0) then
!                      do j = 1, solidforce%dim
!                          call set(solidforce, j, nodes(i), (node_val(alpha, nodes(i)) / dt) * (node_val(fluid_velocity,j,nodes(i))) )
!                      end do
!                  end if
!              end do
!          end do
!
!          if (present(mesh_name)) then
!              ! Add the computed solid force to the field holding all the force-fields:
!              call addto(solidforce_global, solidforce)
!          end if !else solidforce points to the global field, e.g. only 1 solid mesh provided
!
!          ! computing the integral of the solidforce (x-component)
!          solid_force_diag = 0.0
!          solid_force_diag = field_integral(solidforce, fluid_coord)
!
      ! 1-way coupling with prescribed movement:
      !else if (have_option('/embedded_models/fsi_model/one_way_coupling') &
      !     & .and. (have_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed') .or. &
      !     & have_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidMovement/prescribed') ) ) then

          ! Looping over ele, looping over nodes:
          call zero(solidforce)
          do ele = 1, ele_count(solidforce%mesh)
              nodes => ele_nodes(solidforce%mesh, ele)
              do i = 1, size(nodes)
                  if (node_val(alpha, nodes(i)) .gt. 0.0) then
                      do j = 1, solidforce%dim
                          call set(solidforce, j, nodes(i), (beta/dt) * (node_val(fluid_velocity,j,nodes(i)) - &
                                   & (node_val(solid_velocity,j,nodes(i))) ) )
                      end do
                  end if
              end do
          end do

          if (present(mesh_name)) then
              ! Add the computed solid force to the field holding all the force-fields:
              call addto(solidforce_global, solidforce)
          end if !else solidforce points to the global field, e.g. only 1 solid mesh provided

          ! computing the integral of the solidforce
          solid_force_diag = 0.0
          solid_force_diag = field_integral(solidforce, fluid_coord)

      ewrite(2, *) "leaving fsi_model_solid_force_computation"

    ! The following subroutine is an alternative way of computing it and left here for possible future use:
    contains
      subroutine compute_solid_drag_force_ele(field, u, ele)
         type(vector_field), intent(inout) :: field
         type(vector_field), intent(in) :: u
         integer, intent(in) :: ele

         type(vector_field), pointer :: x, fluid_absorption, fluid_velocity
         real, dimension(u%dim, ele_ngi(u,ele)) :: velocity_gi, absorption_gi
         real, dimension(u%dim, ele_ngi(u,ele)) :: field_gi
         real, dimension(u%dim, ele_loc(u,ele),ele_loc(u,ele)) :: field_mat
         type(element_type), pointer :: velocity_shape
         ! Current element global node numbers.
         integer, dimension(:), pointer :: field_nodes
         real, dimension(ele_ngi(field,ele)) :: detwei
         integer :: dim

         ewrite(2,*) "inside compute_solid_drag_force_ele"

         x => extract_vector_field(state, "Coordinate")
         fluid_absorption => extract_vector_field(state, "VelocityAbsorption")
         fluid_velocity => extract_vector_field(state, "Velocity")

         field_nodes => ele_nodes(field, ele)
         velocity_shape => ele_shape(fluid_velocity, ele)

         call transform_to_physical(x, ele, detwei = detwei)

         absorption_gi = ele_val_at_quad(fluid_absorption, ele)
         velocity_gi = ele_val_at_quad(fluid_velocity, ele)
         ! Compute the solid force at gauss points:
         do dim = 1, u%dim
            field_gi(dim,:) = absorption_gi(dim,:) * velocity_gi(dim,:)
         end do

         ! To get the solid force, compute mass matrix, invert it and multiply by mass times force at gauss points
         do dim = 1, fluid_velocity%dim
            field_mat(dim,:,:) = matmul(inverse(shape_shape(velocity_shape, velocity_shape, detwei)), &
            shape_shape(velocity_shape, velocity_shape, detwei*field_gi(dim,:)))
         end do

         do dim = 1, field%dim
            field%val(dim, field_nodes) = sum(field_mat(dim,:,:), 2)
!            call addto(field, field_nodes, field_nodes, field_mat)
         end do
      end subroutine compute_solid_drag_force_ele

    end subroutine fsi_model_solid_force_computation

    !----------------------------------------------------------------------------

    subroutine fsi_model_solid_volume_computation(state, solid_volume_diag, mesh_name)
    !! Computes the solid volume fraction of the solid with name 'mesh_name'
        type(state_type), intent(in) :: state
        real, intent(inout) :: solid_volume_diag
        character(len=OPTION_PATH_LEN), optional ::  mesh_name
        type(vector_field), pointer :: fluid_coord
        type(scalar_field), pointer :: alpha

        ! Get relevant quantities from state
        if(present(mesh_name)) then
            alpha => extract_scalar_field(state, trim(mesh_name)//"SolidConcentration")
        else
            alpha => extract_scalar_field(state, "SolidConcentration")
        end if
        fluid_coord => extract_vector_field(state, "Coordinate")

        ! Compute the integral, thus the volume of immersed body:
        solid_volume_diag = field_integral(alpha, fluid_coord)

    end subroutine fsi_model_solid_volume_computation

    !----------------------------------------------------------------------------

    subroutine fsi_model_compute_diagnostics(state)
    !! Compute all sorts of FSI diagnostics
      type(state_type), intent(inout) :: state
      type(vector_field), pointer :: fluid_coord
      type(vector_field), pointer :: solidforce

      real, dimension(:), allocatable :: solid_force_diag, pre_solid_vel
      real :: solid_volume_diag
      integer :: num_solid_mesh, num_pre_solid_vel, num_pre_solid_pos
      character(len=OPTION_PATH_LEN) :: mesh_name
      integer :: i

      ewrite(2, *) "inside fsi_model_compute_diagnostics"

      fluid_coord => extract_vector_field(state, "Coordinate")
      ! Get solid force (global) field and set it to zero as it is being recomputed below:
      solidforce => extract_vector_field(state, "SolidForce")
      call zero(solidforce)

      ! For 1-way coupling:
      if (have_option('/embedded_models/fsi_model')) then

         ! Compute diagnostic variables, e.g. Force on solid and prescribed Solid Velocity
         if (.not. have_option('/embedded_models/fsi_model/stat/exclude_in_stat')) then 

            num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')
            allocate(solid_force_diag(fluid_coord%dim))

        ! Compute Force and volume of solids:
            solid_mesh_loop: do i=0, num_solid_mesh-1

               ! Get mesh name:
               call get_option('/embedded_models/fsi_model/geometry/mesh['//int2str(i)//']/name', mesh_name)

               ! Resetting force and volume:
               solid_force_diag = 0.0; solid_volume_diag = 0.0
               ! Compute volume of solid:
               call fsi_model_solid_volume_computation(state, solid_volume_diag, mesh_name=mesh_name)
               ! Compute Force of solid:
               call fsi_model_solid_force_computation(state, solid_force_diag, mesh_name=mesh_name)

               ! Set the force of the solid with name 'mesh_name'
               call set_diagnostic(name='ForceOnSolid_'//trim(mesh_name), statistic='Value', value=(/ solid_force_diag /))
               ! Set the volume of the solid with name 'mesh_name'
               call set_diagnostic(name='VolumeOfSolid_'//trim(mesh_name), statistic='Value', value=(/ solid_volume_diag /))

            end do solid_mesh_loop
            deallocate(solid_force_diag)

        ! Now compute diagnostics for the velocity of the solids:
!            num_pre_solid_vel = option_count('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh')
!            
!            ! Now loop over the number of prescribed solid velocities and set those diagnostics:
!            pre_solid_vel_loop: do i=0, num_pre_solid_vel-1
!            allocate(pre_solid_vel(fluid_coord%dim))
!
!               ! Get mesh name:
!               call get_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh['//int2str(i)//']/name', mesh_name)
!               ! Resetting force:
!               pre_solid_vel = 0.0
!               ! Get min/max values of x/y/z component of solid velocity
!               !/* still to be implemented */!
!
!               ! Set the force of the solid with name 'mesh_name'
!               call set_diagnostic(name='VelocityOfSolid_'//trim(mesh_name), statistic='Value', value=(/ pre_solid_vel, pre_solid_vel /))
!               
!               deallocate(pre_solid_vel)
!            end do pre_solid_vel_loop
!
!            ! Same if prescribed solid position is switched on:
!            num_pre_solid_vel = option_count('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidPosition/prescribed/mesh')
!            
!            ! Now loop over the number of prescribed solid velocities and set those diagnostics:
!            pre_solid_pos_loop: do i=0, num_pre_solid_pos-1
!            allocate(pre_solid_vel(fluid_coord%dim))
!
!               ! Get mesh name:
!               call get_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidPosition/prescribed/mesh['//int2str(i)//']/name', mesh_name)
!               ! Resetting force:
!               pre_solid_vel = 0.0
!               ! Get min/max values of x/y/z component of solid velocity
!               !/* still to be implemented */!
!
!               ! Set the force of the solid with name 'mesh_name'
!               call set_diagnostic(name='VelocityOfSolid_'//trim(mesh_name), statistic='Value', value=(/ pre_solid_vel, pre_solid_vel /))
!               
!               deallocate(pre_solid_vel)
!            end do pre_solid_pos_loop
!
         end if

      end if

      ewrite(2, *) "leaving fsi_model_compute_diagnostics"

    end subroutine fsi_model_compute_diagnostics

    !----------------------------------------------------------------------------

    subroutine fsi_model_nonlinear_iteration_converged(state)
      type(state_type), intent(inout) :: state
    !! If a tolerance for nonlinear iterations has been set, this subroutine is called
    !! once the tolerance has been reached, so that the FSI diagnostics are being computed
      if (do_adapt_mesh(current_time, timestep)) then
         if (have_option('/embedded_models/fsi_model')) then
             call fsi_model_compute_diagnostics(state)
         end if
         adapt_at_previous_dt = .true.
      end if

    end subroutine fsi_model_nonlinear_iteration_converged

    !----------------------------------------------------------------------------

    subroutine fsi_initialise(state)
    !! Initialise fields and meshes for FSI problems
      type(state_type), intent(inout) :: state
      type(vector_field), pointer :: fluid_position, solid_force, solid_velocity
      type(vector_field) :: solid_position, solidforce_mesh, solidvelocity_mesh
      type(scalar_field) :: alpha_solidmesh
      type(scalar_field), pointer :: alpha_global
      type(mesh_type) :: solid_mesh
      integer :: quad_degree
      character(len=OPTION_PATH_LEN) :: mesh_path
      character(len=OPTION_PATH_LEN) :: mesh_name, mesh_filename, meshformat
      integer :: i, num_solid_mesh

      ! pointer to vector field of coordinates of fluids mesh:
      fluid_position => extract_vector_field(state, "Coordinate")
      ! Get global solid volume fraction field:
      alpha_global => extract_scalar_field(state, "SolidConcentration")
      ! and the solid force field in state:
      solid_force => extract_vector_field(state, "SolidForce")
      solid_velocity => extract_vector_field(state, "SolidVelocity")
      
      ! figure out if we want to print out diagnostics and initialise files
      ! check for mutiple solids and get translation coordinates
      num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

      ! If we have multiple solids with the exact same geometry, we could just translate their positions via a python function:
      ! The code below might work, but routines are called that are not designed for this, 
      ! so the code below needs to be reviewed and changed if this functionality is wanted!
      !if (translate_solid) then
         !allocate(translation_coordinates(positions%dim, number_of_solids))
         !call get_option("/embedded_models/fsi_model/mesh[0]/python", python_function)
         !call set_detector_coords_from_python(translation_coordinates, &
         !     number_of_solids, python_function, current_time)
      !end if

      ! Loop over number of solid meshes defined:
      solid_mesh_loop: do i=0, num_solid_mesh-1
         ! Get mesh name:
         mesh_path="/embedded_models/fsi_model/geometry/mesh["//int2str(i)//"]"
         call get_option(trim(mesh_path)//'/name', mesh_name)
         call get_option(trim(mesh_path)//'/from_file/file_name', mesh_filename)
         call get_option(trim(mesh_path)//'/from_file/format/name', meshformat)
         call get_option('/embedded_models/fsi_model/geometry/quadrature/degree', quad_degree)
         ! Read in the serial mesh of solid body
         ! For now the FS coupling is parallel for the fluid domain, but serial for the solid domain:
!         solid_position = read_triangle_serial(mesh_filename, quad_degree=quad_degree)
!         solid_position = read_mesh_files(mesh_filename, quad_degree=quad_degree, solidmesh=.true., format=meshformat)
          solid_position = read_mesh_files(mesh_filename, quad_degree=quad_degree, format=meshformat, solid=1)



         ! Set mesh parameters
         solid_mesh = solid_position%mesh
         solid_mesh%name = trim(mesh_name)
         ! Now copy those back to the solid_position field:
         solid_position%mesh = solid_mesh
         solid_position%name = trim(solid_mesh%name)//'SolidCoordinate'

         ! Insert solid_mesh and solid_position into state:
         call insert(state, solid_mesh, trim(solid_mesh%name))
         call insert(state, solid_position, trim(solid_position%name))
         ! Add to extra state:
         call insert(global_solid_state, solid_mesh, trim(solid_mesh%name))
         call insert(global_solid_state, solid_position, trim(solid_position%name))

         ! Also set-up solidvolumefraction field for all solids:
         call allocate(alpha_solidmesh, fluid_position%mesh, trim(solid_mesh%name)//"SolidConcentration")
         alpha_solidmesh%option_path = alpha_global%option_path
         call zero(alpha_solidmesh)
         ! And insert it into state so that we have one solidconcentration field per solidmesh:
         call insert(state, alpha_solidmesh, trim(solid_mesh%name)//"SolidConcentration")
!         ! Add to extra state:
!         call insert(global_fluid_state, alpha_solidmesh, trim(alpha_solidmesh%name))
         call deallocate(alpha_solidmesh)

         ! And the solidforce:
         call allocate(solidforce_mesh, solid_force%dim, solid_force%mesh, trim(solid_mesh%name)//"SolidForce")
         solidforce_mesh%option_path = solid_force%option_path
         call zero(solidforce_mesh)
         call insert(state, solidforce_mesh, trim(solid_mesh%name)//"SolidForce")
!         ! Add to extra state:
!         call insert(global_fluid_state, solidforce_mesh, trim(solidforce_mesh%name))
         call deallocate(solidforce_mesh)
         
         !if (have_option("/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed")) then
            ! And the solidvelocity:
            call allocate(solidvelocity_mesh, solid_velocity%dim, solid_velocity%mesh, trim(solid_mesh%name)//"SolidVelocity")
            solidvelocity_mesh%option_path = solid_velocity%option_path
            call zero(solidvelocity_mesh)
            call insert(state, solidvelocity_mesh, trim(solid_mesh%name)//"SolidVelocity")
            call deallocate(solidvelocity_mesh)
         !end if

         ! Abort if dimensions of fluid and solid mesh don't add up
         !assert(fluid_position%dim == solid_position%dim)         
         
      end do solid_mesh_loop

      ! 1D set-ups not supported:
      assert(fluid_position%dim >= 2)

    end subroutine fsi_initialise

  !----------------------------------------------------------------------------

    subroutine fsi_model_register_diagnostic
    !! Registering FSI diagnostics
      integer :: ndimension
      integer :: num_solid_mesh, num_pre_solid_vel, num_pre_solid_pos
      character(len=OPTION_PATH_LEN) :: mesh_name
      integer :: i

      if (have_option('/embedded_models/fsi_model')) then

         ! Get dimension of problem (from fluid mesh):
         call get_option('/geometry/dimension', ndimension)
         num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

         ! Loop over number of solid meshes defined:
         solid_mesh_loop: do i=0, num_solid_mesh-1

            ! Get mesh name:
            call get_option('/embedded_models/fsi_model/geometry/mesh['//int2str(i)//']/name', mesh_name)

            ! register Force acting on Solid if stat is enabled:
            if (.not. have_option('/embedded_models/fsi_model/stat/exclude_in_stat')) then
               call register_diagnostic(dim=ndimension, name='ForceOnSolid_'//trim(mesh_name), statistic='Value')
               call register_diagnostic(dim=1, name='VolumeOfSolid_'//trim(mesh_name), statistic='Value')
            end if

         end do solid_mesh_loop

!         ! Get number of solid meshes that have prescribed velocity:
!         num_pre_solid_vel = option_count('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh')
!         ! Loop over number of solid meshes that have a prescribed velocity:
!         pre_solid_vel_loop: do i=0, num_pre_solid_vel-1
!
!            ! Store prescribed solid velocity (of each solid, thus of each solid mesh) in statfile
!            if (have_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh['//int2str(i)//']/stat/include_in_stat')) then
!               ! Get corresponding mesh name:
!               call get_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh['//int2str(i)//']/name', mesh_name)
!               ! Register diagnostic variable of corresponding mesh (with mesh_name):
!               call register_diagnostic(dim=ndimension, name='VelocityOfSolid_'//trim(mesh_name), statistic='Value')
!            end if
!
!         end do pre_solid_vel_loop
!
!         ! Same for Prescribed Solid Position:
!         ! Get number of solid meshes that have prescribed position:
!         num_pre_solid_pos = option_count('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidPosition/prescribed/mesh')
!         ! Loop over number of solid meshes that have a prescribed velocity:
!         pre_solid_pos_loop: do i=0, num_pre_solid_pos-1
!            if (have_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidPosition/prescribed/mesh['//int2str(i)//']/stat/include_in_stat')) then
!               ! Get corresponding mesh name:
!               call get_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidPosition/prescribed/mesh['//int2str(i)//']/name', mesh_name)
!               ! Register diagnostic variable of corresponding mesh (with mesh_name):
!               call register_diagnostic(dim=ndimension, name='VelocityOfSolid_'//trim(mesh_name), statistic='Value')
!            end if
!         end do pre_solid_pos_loop

      end if

    end subroutine fsi_model_register_diagnostic

    !----------------------------------------------------------------------------

    subroutine fsi_model_pre_adapt_cleanup(states, solid_states)
    !! Initialise fields and meshes for FSI problems
      type(state_type), dimension(:) :: states
      type(state_type), dimension(:) :: solid_states
      type(scalar_field), pointer :: alpha_solidmesh
      type(vector_field), pointer :: solidforce_mesh, solidvelocity_mesh
      character(len=OPTION_PATH_LEN) :: mesh_path, mesh_name
      integer :: i, j, num_solid_mesh, nstates

      ewrite(2,*) "inside fsi_model_pre_adapt_cleanup"

      ! get number of fluid-states and solid meshes:
      nstates=option_count("/material_phase")
      num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

      state_loop: do j=0, nstates-1

        ewrite(2,*) "state number: ", j

        ! Loop over number of solid meshes defined:
        solid_mesh_loop: do i=0, num_solid_mesh-1

          ! Get mesh name:
          mesh_path="/embedded_models/fsi_model/geometry/mesh["//int2str(i)//"]"
          call get_option(trim(mesh_path)//'/name', mesh_name)
          ewrite(2,*) "removing fields from solid: ", mesh_name

          ! not removing, but changing option paths:
          alpha_solidmesh => extract_scalar_field(states(j+1), trim(mesh_name)//'SolidConcentration')
          solidforce_mesh => extract_vector_field(states(j+1), trim(mesh_name)//'SolidForce')
          solidvelocity_mesh => extract_vector_field(states(j+1), trim(mesh_name)//'SolidVelocity')
          alpha_solidmesh%option_path = 'solidtmp/'//trim(alpha_solidmesh%option_path)
          solidforce_mesh%option_path = 'solidtmp/'//trim(solidforce_mesh%option_path)
          solidvelocity_mesh%option_path = 'solidtmp/'//trim(solidvelocity_mesh%option_path)

        end do solid_mesh_loop

      end do state_loop

      ewrite(2,*) "leaving fsi_model_pre_adapt_cleanup"

    end subroutine fsi_model_pre_adapt_cleanup

    !----------------------------------------------------------------------------

    subroutine fsi_post_adapt_operations(states, solid_states)
    !! Initialise fields and meshes for FSI problems
      type(state_type), dimension(:) :: states
      type(state_type), dimension(:) :: solid_states
      type(scalar_field), pointer :: alpha_global
      type(vector_field), pointer :: solid_velocity, solid_force
      type(scalar_field) :: alpha_solid_fluidmesh
      type(vector_field) :: solidforce_mesh, solidvelocity_mesh
      character(len=OPTION_PATH_LEN) :: mesh_path, state_path, mesh_name
      integer :: i, j, num_solid_mesh, nstates

      ewrite(2,*) "inside fsi_post_adapt_operations"

      ! get number of fluid-states and solid meshes:
      nstates=option_count("/material_phase")
      num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

      ! get some fields from state that we need:
      alpha_global => extract_scalar_field(states(1), "SolidConcentration")
      solid_force => extract_vector_field(states(1), "SolidForce")
      solid_velocity => extract_vector_field(states(1), "SolidVelocity")

      state_loop: do j=0, nstates-1

        ewrite(2,*) "state number: ", j
        state_path = '/material_phase['//int2str(j)//']/'

        ! Loop over number of solid meshes defined:
        solid_mesh_loop: do i=0, num_solid_mesh-1

          ! Get mesh name:
          mesh_path="/embedded_models/fsi_model/geometry/mesh["//int2str(i)//"]"
          call get_option(trim(mesh_path)//'/name', mesh_name)
          ewrite(2,*) "resetting option_paths of solid fields of solid mesh: ", mesh_name

          ! Setting up additional fields on the fluid mesh, e.g. 
          ! one solid volume fraction field per solid mesh on the fluid mesh:
          call allocate(alpha_solid_fluidmesh, alpha_global%mesh, trim(mesh_name)//"SolidConcentration")
          alpha_solid_fluidmesh%option_path = alpha_global%option_path
          call zero(alpha_solid_fluidmesh)
          call insert(states(j+1), alpha_solid_fluidmesh, trim(mesh_name)//"SolidConcentration")
          call deallocate(alpha_solid_fluidmesh)

          ! And the solidforce:
          call allocate(solidforce_mesh, solid_force%dim, solid_force%mesh, trim(mesh_name)//"SolidForce")
          solidforce_mesh%option_path = solid_force%option_path
          call zero(solidforce_mesh)
          call insert(states(j+1), solidforce_mesh, trim(mesh_name)//"SolidForce")
          call deallocate(solidforce_mesh)

          ! And the solidvelocity:
          call allocate(solidvelocity_mesh, solid_velocity%dim, solid_velocity%mesh, trim(mesh_name)//"SolidVelocity")
          solidvelocity_mesh%option_path = solid_velocity%option_path
          call zero(solidvelocity_mesh)
          call insert(states(j+1), solidvelocity_mesh, trim(mesh_name)//"SolidVelocity")
          call deallocate(solidvelocity_mesh)

        end do solid_mesh_loop

      end do state_loop

      ewrite(2,*) "leaving fsi_post_adapt_operations"

    end subroutine fsi_post_adapt_operations

end module fsi_model
