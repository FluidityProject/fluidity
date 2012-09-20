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
             fsi_model_register_diagnostic, fsi_insert_state, &
             fsi_model_pre_adapt_cleanup, fsi_post_adapt_operations

  contains

    subroutine fsi_modelling(state, solid_states, its, itinoi)
    !! Main routine for FSI, being called every picard iteration
        type(state_type), intent(inout) :: state
        type(state_type), intent(inout), dimension(:) :: solid_states
        integer, intent(in) :: its, itinoi
        type(scalar_field), pointer :: alpha_global
        type(scalar_field), pointer :: old_alpha_global
        logical :: recompute_alpha
        logical :: solid_moved
        !type(vector_field), pointer :: test

        ewrite(2, *) "inside fsi_modelling"

        if (timestep == 1 .and. its == 1) then
!            call fsi_initialise(state)
            ! Make sure /alpha_s^f is computed at first timestep:
            adapt_at_previous_dt = .true.
            recompute_alpha = .true.
            solid_moved = .true.
            ! At this stage, everything is initialised
        end if
        
            
        !test => extract_vector_field(state, "frontcylinderSolidCoordinate")


!        ! For future use:
!        ! 1-WAY COUPLING (prescribed solid velocity)
!        ! First check if prescribed solid movement is enabled,
!        ! and if so, move the solid mesh
!        if (have_option("/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed") .or. &
!            have_option("/embedded_models/fsi_model/one_way_coupling/vector_field::SolidPosition/prescribed") &
!            & .and. its == 1) then
!           call fsi_move_solid_mesh(state)
!           ! Check if solid actually moved:
!
!           ! If solid moved, store new coordinates in state:
!           solid_moved = .true.
!           ! ... And recalculate the volume fraction:
!           
!        else
!           solid_moved = .false.
!
!        end if

!        ! Check if alpha needs to be recomputed (at this timestep):
!        recompute_alpha = fsi_recompute_alpha(its, solid_moved=solid_moved)


        ! Do we need to compute the new alpha field(s)?
!        if (recompute_alpha) then
        if (adapt_at_previous_dt) then
            ! projection between solid/fluid mesh:
            call fsi_ibm_projections(state, solid_states)
        end if

        ! 1-WAY COUPLING (prescribed solid velocity)
        ! First check if prescribed solid movement is enabled,
        ! and if so, move the solid mesh and update the new
        ! solid volume fractions as well
        if (have_option("/embedded_models/fsi_model/one_way_coupling/vector_field::SolidPosition/prescribed") &
            & .and. its == 1) then
            call fsi_update_solid_position_volume_fraction(state, solid_states)
        end if


        ! Get alpha (of all solids) from state
        alpha_global => extract_scalar_field(state, "SolidConcentration")
        ! Get old alpha, so that we can overwrite it:
        old_alpha_global => extract_scalar_field(state, "OldSolidConcentration")
        call set(old_alpha_global, alpha_global)

        ! Set absorption term /sigma
        call compute_fluid_absorption(state)
        ! Set source term:
        if (have_option("/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed") .and. its == 1) then
            call compute_source_term(state)
        end if

        ! the variable name might be misleading when reading it as 
        ! holding the value at the end of the timestep, but when 
        ! reading it at the beginning of the next timestep, the name makes sense:
        if (do_adapt_mesh(current_time, timestep) .and. its==itinoi) then
            adapt_at_previous_dt = .true.
        else
            adapt_at_previous_dt = .false.
        end if

        ewrite(2, *) 'leaving fsi_modelling'

    end subroutine fsi_modelling

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

    subroutine fsi_update_solid_position_volume_fraction(state, solid_states)
    !! Subroutine that loops over all solid meshes that have a prescribed
    !! velocity or movement, applies solid movement and 
    !! updates the solid volume fraction fields afterwards
        type(state_type), intent(inout) :: state
        type(state_type), intent(inout), dimension(:) :: solid_states

        type(vector_field) :: solid_movement
        type(vector_field), pointer :: fluid_position
        type(vector_field), pointer :: solid_position, solid_velocity_global
        type(scalar_field), pointer :: solid_alpha_mesh, solid_alpha_global

        character(len=OPTION_PATH_LEN) :: mesh_name
        character(len=OPTION_PATH_LEN) :: fsi_path="/embedded_models/fsi_model/one_way_coupling/vector_field::"
        character(len=PYTHON_FUNC_LEN) :: func
        integer :: num_solid_mesh, num_pre_solid_pos, num_pre_solid_vel
        integer :: i

        logical :: solid_moved

        ewrite(2,*) "Inside fsi_update_solid_position_volume_fraction"

        ! Set logical to false, as at this moment, no solid has moved yet:
        solid_moved = .false.

        ! Set global solid velocity field to zero:
!        solid_velocity_global => extract_vector_field(state, "SolidVelocity")
!        call zero(solid_velocity_global)

        ! Get number of solid meshes:
        num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

        ! Loop over number of solid meshes defined:
        solid_mesh_position_loop: do i=0, num_solid_mesh-1

            ! Get mesh name:
            call get_option('/embedded_models/fsi_model/geometry/mesh['//int2str(i)//']/name', mesh_name)

            ! Check for user error:
            if (have_option(trim(fsi_path)//"SolidPosition/prescribed/mesh::"//trim(mesh_name))) then
            
                call fsi_move_solid_mesh(state, solid_states, mesh_name, solid_moved=solid_moved)
                ewrite(1,*) "moved solid mesh: ", trim(mesh_name)
                solid_moved = .true.
            else
                solid_moved = .false.
            end if

            ! Recompute solid volume fraction, if solid just moved:
            if (solid_moved) then
                ewrite(1,*) "update alpha of solid mesh: ", trim(mesh_name)
                call fsi_ibm_projections_mesh(state, solid_states, mesh_name)
            end if

        end do solid_mesh_position_loop


        ! Update the global solid volume fraction:
        solid_alpha_global => extract_scalar_field(state, 'SolidConcentration')
        call zero(solid_alpha_global)

        ! Loop over solid meshes again, to update global solid volume fraction field:
        solid_mesh_alpha_loop: do i=0, num_solid_mesh-1

            ! Get mesh name:
            call get_option('/embedded_models/fsi_model/geometry/mesh['//int2str(i)//']/name', mesh_name)
            
            ! Get solid volume fraction field of this mesh:
            solid_alpha_mesh => extract_scalar_field(state, trim(mesh_name)//'SolidConcentration')
            
            ! Add mesh solid volume fraction to the global solid volume fraction:
            call addto(solid_alpha_global, solid_alpha_mesh)
            
        end do solid_mesh_alpha_loop

        ! Make sure that the solid volume fraction is max 1.0,
        ! only necessary if more than one solid mesh is present
        if (num_solid_mesh .gt. 1) then
            call sum_union_solid_volume_fraction(state)
        end if

        ewrite(2,*) "Leaving fsi_update_solid_position_volume_fraction"


    end subroutine fsi_update_solid_position_volume_fraction

    !----------------------------------------------------------------------------

    subroutine fsi_ibm_projections(state, solid_states)
    !! Subroutine that loops over all solid meshes and calls the corresponding 
    !! subroutines to obtain the solid volume fraction
        type(state_type), intent(inout) :: state
        type(state_type), intent(inout), dimension(:) :: solid_states
        type(scalar_field), pointer :: alpha_global

        character(len=OPTION_PATH_LEN) :: mesh_name
        integer :: num_solid_mesh
        integer :: i

        ewrite(2,*) "inside fsi_ibm_projections"

        ! Get relevant 
        alpha_global => extract_scalar_field(state, "SolidConcentration")
        ! Since we are (re)computing the solid volume fraction field, set the 'old' field to zero:
        call zero(alpha_global)

        ! Get number of solid meshes:
        num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

        ! Loop over number of solid meshes defined:
        solid_mesh_loop: do i=0, num_solid_mesh-1

            ! Get mesh name:
            call get_option('/embedded_models/fsi_model/geometry/mesh['//int2str(i)//']/name', mesh_name)

            ! Now do the projection to obtain the alpha field of the current solid mesh:
            call fsi_ibm_projections_mesh(state, solid_states, mesh_name)

        end do solid_mesh_loop

        ! Make sure that the solid volume fraction is max 1.0,
        ! only necessary if more than one solid mesh is present
        if (num_solid_mesh .gt. 1) then
            call sum_union_solid_volume_fraction(state)
        end if

        ! At this point, the new solid volume fraction of all given solid meshes has been (re)computed and stored in SolidConcentration in state!

        ewrite(2,*) "leaving fsi_ibm_projections"

    end subroutine fsi_ibm_projections

    !----------------------------------------------------------------------------

    subroutine fsi_ibm_projections_mesh(state, solid_states, mesh_name)
        type(state_type), intent(inout) :: state
        type(state_type), intent(inout), dimension(:) :: solid_states
        character(len=OPTION_PATH_LEN), intent(in) :: mesh_name

        type(vector_field), pointer :: fluid_position, fluid_velocity
        type(scalar_field), pointer :: alpha_global, alpha_solidmesh
        type(vector_field) :: solid_position_mesh
        type(scalar_field) :: alpha_tmp

        ewrite(2,*) "inside fsi_ibm_projections_mesh"

        ! Get relevant fields:
        fluid_velocity => extract_vector_field(state, "Velocity")
        fluid_position => extract_vector_field(state, "Coordinate")
        alpha_global => extract_scalar_field(state, "SolidConcentration")

        ! Get coordinate field of the current solid (mesh):
        solid_position_mesh = extract_vector_field(state, trim(mesh_name)//'SolidCoordinate')
        ! Also the according solid volume fraction field:
        alpha_solidmesh => extract_scalar_field(state, trim(mesh_name)//'SolidConcentration')

        ! Allocate temp solid volume fraction field for the projection:
        call allocate(alpha_tmp, alpha_solidmesh%mesh, 'TMPSolidConcentration')
        call zero(alpha_tmp)


        ! Distinguish between different projection methods:

        ! Galerkin projection via supermesh:
        if (have_option("/embedded_models/fsi_model/one_way_coupling/inter_mesh_projection/galerkin_projection")) then

            ! Computing /alpha_s^f by projecting unity from the solid mesh to the fluid mesh:
            call fsi_one_way_galerkin_projection(fluid_velocity, fluid_position, solid_position_mesh, alpha_tmp)

        else if (have_option("/embedded_models/fsi_model/one_way_coupling/inter_mesh_projection/grandy_interpolation")) then
            ! Grandy interpolation:
            call fsi_one_way_grandy_interpolation(fluid_position, solid_position_mesh, alpha_tmp)

        !else
            ! Add more interpolations here

        end if


        ! After projection, set the solid volume fraction of the current solid:
        call set(alpha_solidmesh, alpha_tmp)
        ewrite_minmax(alpha_solidmesh)

        ! And add the alpha of the current solid to the global alpha (of all solids):
        call addto(alpha_global, alpha_solidmesh)

        ! Deallocate temp arrays:
        call deallocate(alpha_tmp)

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
          if (have_option("/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed")) then
             solid_velocity => extract_vector_field(state, "SolidVelocity")
          end if

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
             if (have_option("/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed")) then
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
        type(vector_field), pointer :: absorption
        type(scalar_field), pointer :: alpha_global
        integer :: i, j, ele
        real :: sigma, beta
        integer, dimension(:), pointer :: nodes

        ewrite(2, *) 'inside compute_fluid_absorption'
      
        !fluid_velocity => extract_vector_field(state, 'Velocity')
        alpha_global => extract_scalar_field(state, "SolidConcentration")
        ! Get absorption velocity from state:
        absorption => extract_vector_field(state, 'VelocityAbsorption')
        call zero(absorption)

        do ele = 1, ele_count(alpha_global%mesh)
            nodes => ele_nodes(alpha_global%mesh, ele)
            do i = 1, size(nodes)
                sigma = node_val(alpha_global, nodes(i)) / dt
                do j = 1, mesh_dim(absorption%mesh)
                    call set(absorption, j, nodes(i), sigma)
                end do
            end do
        end do
        
        if (have_option('/embedded_models/fsi_model/one_way_coupling/beta')) then
            call get_option('/embedded_models/fsi_model/one_way_coupling/beta', beta)
            call scale(absorption, beta)
        end if

        ewrite(2, *) "leaving compute_fluid_absorption"

    end subroutine compute_fluid_absorption

  !----------------------------------------------------------------------------

    subroutine fsi_move_solid_mesh(state, solid_states, mesh_name, solid_moved)
    !! Moving a solid mesh
      type(state_type), intent(in) :: state
      type(state_type), intent(inout), dimension(:) :: solid_states
      character(len=OPTION_PATH_LEN), intent(in) :: mesh_name
      logical, intent(inout), optional :: solid_moved

      type(vector_field), pointer :: solid_position_mesh, solid_velocity_global
      type(vector_field) :: solid_movement_mesh, solid_velocity_mesh

      character(len=PYTHON_FUNC_LEN) :: func      


      ewrite(2,*) "inside move_solid_mesh"

      ! 1st get coordinate field of this solid mesh, and its velocity field (which is on the fluid mesh):
      solid_position_mesh => extract_vector_field(state, trim(mesh_name)//"SolidCoordinate")
      ewrite(2,*) "after getting solid position field"
      ! and create new field in which the travelled distance is stored:
      call allocate(solid_movement_mesh, solid_position_mesh%dim, solid_position_mesh%mesh, name=trim(mesh_name)//"SolidMovement")
      call zero(solid_movement_mesh)
      ewrite(2,*) "after allocating and setting solid movement field"
      ! and create field for its velocity: (for now, we'll comment it out):
!      call allocate(solid_velocity_mesh, solid_position_mesh%dim, solid_position_mesh%mesh, name=trim(mesh_name)//"SolidVelocity")
!      call zero(solid_velocity_mesh)
!      ewrite(2,*) "after getting solid velocity field"

      ! 2nd get the python function for this solid mesh:
      call get_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidPosition/prescribed/mesh::'//trim(mesh_name)//'/python', func)
      ! 3rd Set the new coordinates from the python function:
      call set_from_python_function(solid_movement_mesh, func, solid_position_mesh, dt)
      ewrite(2,*) "after getting python function"
      
      ! 4th Check if the mesh has actually moved:
!      if (.not. solid_movement_mesh .eq. zero) then
!      
!      end if
!      ewrite(2,*) "max(solid_movement_mesh) = ", maxval(solid_movement_mesh(0))
!      ewrite(2,*) "min(solid_movement_mesh) = ", minval(solid_movement_mesh(0))

      ! 4th Set the new coordinates of the solid:
      call addto(solid_position_mesh, solid_movement_mesh)
      ewrite(2,*) "after adding solid movement to position field"

      ! 5th Set the solid velocity field of this solid, and addto global velocity field:
      ! comment it out for now, as global and local velocity fields live on different coordinate meshes:
!      call addto(solid_velocity_mesh, solid_movement_mesh, 1.0/dt)
!      ewrite(2,*) "after adding solid movement per dt  to velocity field"
!      solid_velocity_global => extract_vector_field(state, "SolidVelocity")
!      ewrite(2,*) "after getting pointer to global solid velocity field"
      ! solid velocity lives on solid coordinate mesh, 
      ! thus solid_velocity_mesh cannot be added to the global solid velocity field
      ! as this lives on the fluid coordinate mesh. 
      ! To get the solid velocity on the fluid mesh, we have to project it from the
      ! solid mesh. For now, we'll comment it out.
!      call addto(solid_velocity_global, solid_velocity_mesh)
!      ewrite(2,*) "after adding local solid velocity to global velocity field"

      ! Deallocate:
      call deallocate(solid_movement_mesh)
      ewrite(2,*) "after deallocating solid movement field"


      ewrite(2,*) 'end of move_solid_mesh'

    end subroutine fsi_move_solid_mesh

    !----------------------------------------------------------------------------

! For future use:
    subroutine compute_source_term(state)

      type(state_type), intent(inout) :: state

      type(vector_field), pointer :: source_term
      type(vector_field), pointer :: fluid_velocity, solid_velocity, fluid_absorption
      type(scalar_field), pointer :: alpha

      integer, dimension(:), pointer :: nodes
      integer :: i, j, ele

      if (have_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed')) then

         fluid_velocity => extract_vector_field(state, "Velocity")
         fluid_absorption => extract_vector_field(state, "VelocityAbsorption")
         solid_velocity => extract_vector_field(state, "SolidVelocity")
         alpha => extract_scalar_field(state, "SolidConcentration")

         source_term => extract_vector_field(state, "VelocitySource")
         call zero(source_term)

         ! (\rho_f \alpha_s / \Delta t) (\hat{u_f} or u) - F_s/ \Delta t
         do ele = 1, ele_count(source_term%mesh)
             nodes => ele_nodes(source_term%mesh, ele)
             do i = 1, size(nodes)
                 if (node_val(alpha, nodes(i)) .gt. 0.0) then
                     do j = 1, source_term%dim
                         call set(source_term, j, nodes(i), node_val(solid_velocity,j,nodes(i)) / dt  + (-1.0 * node_val(alpha, nodes(i)) / dt) * (node_val(fluid_velocity,j,nodes(i)) )  )
!                         call set(source_term, j, nodes(i), (node_val(alpha, nodes(i)) / dt) * (node_val(fluid_velocity,j,nodes(i)) - (node_val(solid_velocity,j,nodes(i))) ) )
                     end do
                 end if
             end do
         end do

      end if


!      if (have_pressure_gradient) then
!
!         source => extract_vector_field(state, "VelocitySource")
!         call zero(source)
!         call set(source, pressure_gradient)
!         ! remove the source from the solids
!         do i = 1, source%dim
!            source%val(i,:) = source%val(i,:) * (1. - alpha_sf%val)
!         end do
!
!      end if
!
    end subroutine compute_source_term

    !----------------------------------------------------------------------------
! For future use:
!    subroutine add_mass_source_absorption(ct_rhs, state)
!
!      type(state_type), intent(in) :: state
!      type(scalar_field), intent(inout) :: ct_rhs
!  
!      type(scalar_field), pointer :: solid, old_solid, p
!      type(vector_field), pointer :: x
!      type(element_type), pointer :: p_shape 
!      type(element_type) :: test_function
!      integer :: ele
!
!      ! if solving for the bulk
!      ! velocity nothing to be done
!
!      solid => extract_scalar_field(state, "SolidConcentration")
!      old_solid => extract_scalar_field(state, "OldSolidConcentration")
!      x => extract_vector_field(state, "Coordinate")
!      p => extract_scalar_field(state, "Pressure")
!
!      do ele = 1, element_count(p)
!         p_shape => ele_shape(p, ele)
!         test_function = p_shape
!
!         call add_ct_rhs_element_cg(ele, test_function, &
!              p_shape, x, p, solid, old_solid, ct_rhs)
!      end do
!
!    contains
!
!      subroutine add_ct_rhs_element_cg(ele, test_function, &
!           p_shape, x, p, solid, old_solid, ct_rhs)
!
!        type(scalar_field), intent(inout) :: ct_rhs
!        integer, intent(in) :: ele
!        type(element_type), intent(in) :: test_function, p_shape
!        type(vector_field), intent(in) :: x
!        type(scalar_field), intent(in) :: solid, old_solid, p
!
!        real, dimension(ele_ngi(p, ele)) :: detwei
!        real, dimension(ele_loc(p, ele), ele_ngi(p, ele), x%dim) :: dp_t
!        integer, dimension(:), pointer :: p_ele
!        real, dimension(ele_loc(p, ele)) :: mat
!        real, dimension(ele_ngi(solid, ele)) :: ds
!
!        call transform_to_physical(x, ele, &
!             p_shape, dshape=dp_t, detwei=detwei)
!
!        ds = ele_val_at_quad(solid, ele) - ele_val_at_quad(old_solid, ele)
!        mat = shape_rhs(test_function, detwei*ds/dt)
!
!        p_ele => ele_nodes(p, ele)
!
!        call addto(ct_rhs, p_ele, mat)
!
!      end subroutine add_ct_rhs_element_cg
!
!    end subroutine add_mass_source_absorption

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

      ewrite(2, *) "inside fsi_model_solid_force_computation"

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
      if (have_option("/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed")) then
         solid_velocity => extract_vector_field(state, "SolidVelocity")
      end if

      ! 1-way coupling:
      if (have_option('/embedded_models/fsi_model/one_way_coupling') &
           & .and. (.not. have_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed'))) then

          ! Looping over ele, looping over nodes:
          call zero(solidforce)
          do ele = 1, ele_count(solidforce%mesh)
              nodes => ele_nodes(solidforce%mesh, ele)
              do i = 1, size(nodes)
                  if (node_val(alpha, nodes(i)) .gt. 0.0) then
                      do j = 1, solidforce%dim
                          call set(solidforce, j, nodes(i), (node_val(alpha, nodes(i)) / dt) * (node_val(fluid_velocity,j,nodes(i))) )
                      end do
                  end if
              end do
          end do
          
          if (present(mesh_name)) then
              ! Add the computed solid force to the field holding all the force-fields:
              call addto(solidforce_global, solidforce)
          end if !else solidforce points to the global field, e.g. only 1 solid mesh provided
          
          ! computing the integral of the solidforce (x-component)
          solid_force_diag = 0.0
          solid_force_diag = field_integral(solidforce, fluid_coord)

      ! 1-way coupling with prescribed movement:
      else if (have_option('/embedded_models/fsi_model/one_way_coupling') &
           & .and. have_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed')) then

          ! Looping over ele, looping over nodes:
          call zero(solidforce)
          do ele = 1, ele_count(solidforce%mesh)
              nodes => ele_nodes(solidforce%mesh, ele)
              do i = 1, size(nodes)
                  if (node_val(alpha, nodes(i)) .gt. 0.0) then
                      do j = 1, solidforce%dim
                          call set(solidforce, j, nodes(i), (node_val(alpha, nodes(i)) / dt) * (node_val(fluid_velocity,j,nodes(i)) - (node_val(solid_velocity,j,nodes(i))) ) )
                      end do
                  end if
              end do
          end do

          if (present(mesh_name)) then
              ! Add the computed solid force to the field holding all the force-fields:
              call addto(solidforce_global, solidforce)
          end if !else solidforce points to the global field, e.g. only 1 solid mesh provided
          
          ! computing the integral of the solidforce (x-component)
          solid_force_diag = 0.0
          solid_force_diag = field_integral(solidforce, fluid_coord)

      ! 2-way coupling
      else if (have_option('/embedded_models/fsi_model/two_way_coupling')) then
         ! Do nothing at this stage
      end if

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
      integer :: num_solid_mesh, num_pre_solid_vel
      character(len=OPTION_PATH_LEN) :: mesh_name
      integer :: i

      ewrite(2, *) "inside fsi_model_compute_diagnostics"

      fluid_coord => extract_vector_field(state, "Coordinate")
      ! Get solid force (global) field and set it to zero as it is being recomputed below:
      solidforce => extract_vector_field(state, "SolidForce")
      call zero(solidforce)

      ! For 1-way coupling:
      if (have_option('/embedded_models/fsi_model/one_way_coupling')) then

         ! Compute diagnostic variables, e.g. Force on solid and prescribed Solid Velocity
         if (.not. have_option('/embedded_models/fsi_model/one_way_coupling/stat/exclude_in_stat')) then 

            ! Currently wrong, as the force over the whole fluid domain would be computed,
            ! thus it would be the sum of all forces on all solids,
            ! but it works when only one solid mesh is present and
            ! when computing the force on the solid mesh, which should be done anyway!!!
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
            num_pre_solid_vel = option_count('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh')
            
            ! Now loop over the number of prescribed solid velocities and set those diagnostics:
            pre_solid_vel_loop: do i=0, num_pre_solid_vel-1
            allocate(pre_solid_vel(fluid_coord%dim))

               ! Get mesh name:
               call get_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh['//int2str(i)//']/name', mesh_name)
               ! Resetting force:
               pre_solid_vel = 0.0
               ! Get min/max values of x/y/z component of solid velocity
               !/* still to be implemented */!

               ! Set the force of the solid with name 'mesh_name'
               call set_diagnostic(name='VelocityOfSolid_'//trim(mesh_name), statistic='Value', value=(/ pre_solid_vel /))
               
               deallocate(pre_solid_vel)
            end do pre_solid_vel_loop

         end if
      else if (have_option('/embedded_models/fsi_model/two_way_coupling')) then
         ! Do nothing at this stage
      end if
      
      ! If the mesh is about to be adapted at the end of this timestep,
      ! remove additional fields from state:
!      if (adapt_at_previous_dt) then
!         call fsi_remove_from_state(state)
!      end if

      ewrite(2, *) "leaving fsi_model_compute_diagnostics"

    end subroutine fsi_model_compute_diagnostics

    !----------------------------------------------------------------------------

    subroutine fsi_model_nonlinear_iteration_converged(state)
      type(state_type), intent(inout) :: state
    !! If a tolerance for nonlinear iterations has been set, this subroutine is called
    !! once the tolerance has been reached, so that the FSI diagnostics are being computed
      if (do_adapt_mesh(current_time, timestep)) then
         if (have_option('/embedded_models/fsi_model/one_way_coupling')) then
             call fsi_model_compute_diagnostics(state)
         else if (have_option('/embedded_models/fsi_model/two_way_coupling')) then
            ! Do nothing
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
      if (have_option("/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed")) then
         ! and the solid velocity field in state:
         solid_velocity => extract_vector_field(state, "SolidVelocity")
      end if
      
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
         
         if (have_option("/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed")) then
            ! And the solidvelocity:
            call allocate(solidvelocity_mesh, solid_velocity%dim, solid_velocity%mesh, trim(solid_mesh%name)//"SolidVelocity")
            solidvelocity_mesh%option_path = solid_velocity%option_path
            call zero(solidvelocity_mesh)
            call insert(state, solidvelocity_mesh, trim(solid_mesh%name)//"SolidVelocity")
            call deallocate(solidvelocity_mesh)
         end if

         ! Abort if dimensions of fluid and solid mesh don't add up
         assert(fluid_position%dim == solid_position%dim)         

      end do solid_mesh_loop

      ! 1D set-ups not supported:
      assert(fluid_position%dim >= 2)

    end subroutine fsi_initialise

  !----------------------------------------------------------------------------

    subroutine fsi_model_register_diagnostic
    !! Registering FSI diagnostics
      integer :: ndimension
      integer :: num_solid_mesh, num_pre_solid_vel
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
            if (.not. have_option('/embedded_models/fsi_model/one_way_coupling/stat/exclude_in_stat')) then
               call register_diagnostic(dim=ndimension, name='ForceOnSolid_'//trim(mesh_name), statistic='Value')
               call register_diagnostic(dim=1, name='VolumeOfSolid_'//trim(mesh_name), statistic='Value')
            end if

         end do solid_mesh_loop

         ! Get number of solid meshes that have prescribed velocity:
         num_pre_solid_vel = option_count('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh')
         ! Loop over number of solid meshes that have a prescribed velocity:
         pre_solid_vel_loop: do i=0, num_pre_solid_vel-1

            ! Store prescribed solid velocity (of each solid, thus of each solid mesh) in statfile
            if (have_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh['//int2str(i)//']/stat/include_in_stat')) then
               ! Get corresponding mesh name:
               call get_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh['//int2str(i)//']/name', mesh_name)
               ! Register diagnostic variable of corresponding mesh (with mesh_name):
               call register_diagnostic(dim=ndimension, name='VelocityOfSolid_'//trim(mesh_name), statistic='Value')
            end if

         end do pre_solid_vel_loop

      end if

    end subroutine fsi_model_register_diagnostic

    !----------------------------------------------------------------------------

    subroutine fsi_model_pre_adapt_cleanup(states, solid_states)
    !! Initialise fields and meshes for FSI problems
      type(state_type), dimension(:) :: states
      type(state_type), dimension(:) :: solid_states
      type(scalar_field), pointer :: alpha_solidmesh
      type(vector_field), pointer :: solid_position_mesh, solidforce_mesh, solidvelocity_mesh
      type(mesh_type), pointer :: solid_mesh
      character(len=OPTION_PATH_LEN) :: mesh_path, mesh_name
      integer :: i, j, num_solid_mesh, nstates
      integer :: stat

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

!          call remove_scalar_field(states(j+1), trim(mesh_name)//"SolidConcentration", stat)
!          ewrite(2,*) "stat = ", stat
!          call remove_vector_field(states(j+1), trim(mesh_name)//"SolidForce", stat)
!          ewrite(2,*) "stat = ", stat
!          call remove_vector_field(states(j+1), trim(mesh_name)//"SolidVelocity", stat)
!          ewrite(2,*) "stat = ", stat

          ! Get solid specific fields and mesh:
          solid_mesh => extract_mesh(states(j+1), trim(mesh_name)//"SolidMesh")
          solid_position_mesh => extract_vector_field(states(j+1), trim(mesh_name)//"SolidCoordinate")
          alpha_solidmesh => extract_scalar_field(states(j+1), trim(mesh_name)//'SolidConcentration')
          solidforce_mesh => extract_vector_field(states(j+1), trim(mesh_name)//'SolidForce')
          solidvelocity_mesh => extract_vector_field(states(j+1), trim(mesh_name)//'SolidVelocity')
          
          call remove_scalar_field(states(j+1), trim(mesh_name)//"SolidConcentration", stat)
          alpha_solidmesh => extract_scalar_field(states(j+1), trim(mesh_name)//'SolidConcentration')
          
          !test => extract_scalar_field(states(j+1), trim(mesh_name)//'TEST')
          !FLExit("end")

          ! not removing, but changing option paths:
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
      type(scalar_field), pointer :: alpha_solidmesh
      type(vector_field), pointer :: solid_position_mesh, solidforce_mesh, solidvelocity_mesh
      character(len=OPTION_PATH_LEN) :: mesh_path, state_path, mesh_name
      integer :: i, j, num_solid_mesh, nstates
      integer :: stat

      ewrite(2,*) "inside fsi_post_adapt_operations"

      ! get number of fluid-states and solid meshes:
      nstates=option_count("/material_phase")
      num_solid_mesh = option_count('/embedded_models/fsi_model/geometry/mesh')

      state_loop: do j=0, nstates-1

        ewrite(2,*) "state number: ", j
        state_path = '/material_phase['//int2str(j)//']/'

        ! Loop over number of solid meshes defined:
        solid_mesh_loop: do i=0, num_solid_mesh-1

          ! Get mesh name:
          mesh_path="/embedded_models/fsi_model/geometry/mesh["//int2str(i)//"]"
          call get_option(trim(mesh_path)//'/name', mesh_name)
          ewrite(2,*) "resetting option_paths of solid fields of solid mesh: ", mesh_name

          ! not removing, but changing option paths:
          !test => extract_scalar_field(states(j+1), trim(mesh_name)//'TEST')
          
          solid_position_mesh => extract_vector_field(states(j+1), trim(mesh_name)//"SolidCoordinate")
          alpha_solidmesh => extract_scalar_field(states(j+1), trim(mesh_name)//'SolidConcentration')
          alpha_solidmesh%option_path = trim(state_path)//'scalar_field::SolidConcentration'
          solidforce_mesh => extract_vector_field(states(j+1), trim(mesh_name)//'SolidForce')
          solidforce_mesh%option_path = trim(state_path)//'vector_field::SolidForce'
          solidvelocity_mesh => extract_vector_field(states(j+1), trim(mesh_name)//'SolidVelocity')
          solidvelocity_mesh%option_path = trim(state_path)//'vector_field::SolidVelocity'
        end do solid_mesh_loop
      
      end do state_loop

      ewrite(2,*) "leaving fsi_post_adapt_operations"

    end subroutine fsi_post_adapt_operations

end module fsi_model
