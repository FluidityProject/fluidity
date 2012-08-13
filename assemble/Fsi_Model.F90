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

  private
  
  public :: fsi_modelling, fsi_model_compute_diagnostics, &
            fsi_model_nonlinear_iteration_converged, &
            fsi_model_check_options, fsi_model_register_diagnostic

  contains

    subroutine fsi_modelling(state, its, itinoi)
    !! Main routine for FSI, being called every picard iteration
        type(state_type), intent(inout) :: state
        integer, intent(in) :: its, itinoi
        type(scalar_field), pointer :: alpha_global
        type(scalar_field), pointer :: old_alpha_global
        logical :: recompute_alpha
        logical :: solid_moved

        ewrite(2, *) "inside fsi_modelling"

        if (timestep == 1 .and. its == 1) then
            call fsi_initialise(state)
            ! Make sure /alpha_s^f is computed at first timestep:
            adapt_at_previous_dt = .true.
            recompute_alpha = .true.
            solid_moved = .true.
            ! At this stage, everything is initialised:
        end if

        ! Check if alpha needs to be recomputed (at this timestep):
        recompute_alpha = fsi_recompute_alpha(its, solid_moved=solid_moved)

        ! For future use:
        ! 1-WAY COUPLING (prescribed solid velocity)
        ! First check if prescribed solid movement is enabled,
        ! and if so, move the solid mesh
        if (have_option("/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed") &
            & .and. its == 1) then
           call fsi_move_solid_mesh(state)
           ! Check if solid actually moved:
        
           ! If solid moved, store new coordinates in state:
           solid_moved = .true.
           ! ... And recalculate the volume fraction:
           
        else
           solid_moved = .false.
        
        end if

        ! Do we need to compute the new alpha field(s)?
        if (recompute_alpha) then
            call fsi_ibm_projections(state)
        end if

        ! Get alpha (of all solids) from state
        alpha_global => extract_scalar_field(state, "SolidConcentration")
        ! Get old alpha, so that we can overwrite it:
        old_alpha_global => extract_scalar_field(state, "OldSolidConcentration")
        call set(old_alpha_global, alpha_global)

        ! Set absorption term /sigma
        call compute_fluid_absorption(state)

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

    subroutine fsi_ibm_projections(state)
    !! Subroutine that loops over all solid meshes and calls the corresponding 
    !! subroutines to obtain the solid volume fraction
        type(state_type), intent(inout) :: state
        
        type(vector_field), pointer :: solid_position_mesh
        type(scalar_field), pointer :: alpha_global, alpha_solidmesh
        
        type(scalar_field) :: alpha_tmp
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

            ! Get coordinate field of the current solid (mesh):
            solid_position_mesh => extract_vector_field(state, trim(mesh_name)//'SolidCoordinate')
            ! Also the according solid volume fraction field:
            alpha_solidmesh => extract_scalar_field(state, trim(mesh_name)//'SolidConcentration')

            ! allocate temp solid volume fraction field for the projection:
            call allocate(alpha_tmp, alpha_solidmesh%mesh, 'TMPSolidConcentration')
            call zero(alpha_tmp)

            ! Now do the projection to obtain new alpha field:
            call fsi_ibm_projections_mesh(state, solid_position_mesh, alpha_tmp)

            ! After projection, set the solid volume fraction of the current solid:
            call set(alpha_solidmesh, alpha_tmp)

            ! And add the alpha of the current solid to the global alpha (of all solids):
            call addto(alpha_global, alpha_solidmesh)

            ! Deallocate temp arrays:
            call deallocate(alpha_tmp)

        end do solid_mesh_loop
        
        ! At this point, the new solid volume fraction of all given solid meshes has been (re)computed and stored in SolidConcentration in state!

        ewrite(2,*) "leaving fsi_ibm_projections"
    
    end subroutine fsi_ibm_projections

    !----------------------------------------------------------------------------

    subroutine fsi_ibm_projections_mesh(state, solid_position, alpha_tmp)
        type(state_type), intent(inout) :: state
        type(vector_field), pointer, intent(inout) :: solid_position
        type(scalar_field), intent(inout) :: alpha_tmp
        type(vector_field), pointer :: fluid_position, fluid_velocity
        
        ewrite(2,*) "inside fsi_ibm_projections_mesh"
        
        ! Get relevant fields:
        fluid_velocity => extract_vector_field(state, "Velocity")
        fluid_position => extract_vector_field(state, "Coordinate")

        ! Distinguish between different projection methods:

        ! Galerkin projection via supermesh:
        if (have_option("/embedded_models/fsi_model/one_way_coupling/inter_mesh_projection/galerkin_projection")) then

            ! Computing /alpha_s^f by projecting unity from the solid mesh to the fluid mesh:
            call fsi_one_way_galerkin_projection(fluid_velocity, fluid_position, solid_position, alpha_tmp)

        else if (have_option("/embedded_models/fsi_model/one_way_coupling/inter_mesh_projection/grandy_interpolation")) then
            ! Grandy interpolation:
            call fsi_one_way_grandy_interpolation(fluid_position, solid_position, alpha_tmp)
        
        !else
            ! Add more interpolations here
        
        end if

        ewrite_minmax(alpha_tmp)

        ewrite(2,*) "leaving fsi_ibm_projections_mesh"
    
    end subroutine fsi_ibm_projections_mesh

    !----------------------------------------------------------------------------

    subroutine fsi_compute_intersection_map(state, solid_position, map_SF)
    !! This subroutine computes a list of intersecting elements between a solid and fluid mesh
    !! The rtree intersection finder is used as it is more robust in parallel than the 
    !! advancing front intersection finder, e.g. it doesn't abort if no intersection was found,
    !! which is possible with the solid domain being immersed in a fluid domain
        type(state_type), intent(inout) :: state
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
        real :: sigma
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

        ewrite(2, *) "leaving compute_fluid_absorption"

    end subroutine compute_fluid_absorption

  !----------------------------------------------------------------------------

    subroutine fsi_move_solid_mesh(state)
    !! Moving a solid mesh
      type(state_type), intent(in) :: state
      type(vector_field), pointer :: solid_position
      type(vector_field) :: solid_movement
      integer :: i, num_pre_solid_vel, num_solid_mesh
      character(len=PYTHON_FUNC_LEN) :: func
      character(len=OPTION_PATH_LEN) :: mesh_name

      ewrite(2,*) "inside move_solid_mesh"

      ! Get number of solid meshes that have a prescribed velocity:
      num_pre_solid_vel = option_count('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh')

         ! Loop over number of solid meshes that have a prescribed velocity and set their new coordinates:
         pre_solid_vel_loop: do i=0, num_pre_solid_vel-1

            ! 1st get mesh name:
            call get_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh['//int2str(i)//']/name', mesh_name)

            ! 2nd get coordinate field of this solid mesh:
            solid_position => extract_vector_field(state, trim(mesh_name)//"SolidCoordinate")

            ! 3nd get the python function for this solid mesh:
            call get_option('/embedded_models/fsi_model/one_way_coupling/vector_field::SolidVelocity/prescribed/mesh['//int2str(i)//']/python', func)

            ! 4 set the field based on the python function:
            call allocate(solid_movement, solid_position%dim, solid_position%mesh, name=trim(mesh_name)//"SolidMovement")
            call zero(solid_movement)

            call set_from_python_function(solid_movement, func, solid_position, current_time)
            
            ! 5 set solid position field:
            solid_position => extract_vector_field(state, trim(mesh_name)//"SolidCoordinate")
            call addto(solid_position, solid_movement, dt)

            ! 6 Deallocate dummy vector field:
            call deallocate(solid_movement)

         end do pre_solid_vel_loop

  !    call vtk_write_fields("solid_mesh_test", index=2, position=solid_position, model=solid_position%mesh, vfields=(/solid_movement/))

      ewrite(2,*) 'end of move_solid_mesh'

    end subroutine fsi_move_solid_mesh

    !----------------------------------------------------------------------------

! For future use:
!    subroutine set_source(state)
!
!      type(state_type), intent(inout) :: state
!
!      type(vector_field), pointer :: source, positions
!      type(vector_field), pointer :: velocity, absorption
!      type(scalar_field), pointer :: Tsource
!      integer :: i, particle
!      real :: x0, y0, z0, x, y, z, sigma
!
!      if (two_way_coupling) then
!
!         velocity => extract_vector_field(state, "Velocity")
!         absorption => extract_vector_field(state, "VelocityAbsorption")
!
!         source => extract_vector_field(state, "VelocitySource")
!         call zero(source)
!
!         ! (\rho_f \alpha_s / \Delta t) (\hat{u_f} or u) - F_s/ \Delta t
!         do i = 1, node_count(source)
!            call set(source, i, &
!                 node_val(absorption, i) * node_val(velocity, i) &
!                 - node_val(fl_pos_solid_force, i)/dt)
!         end do
!
!      end if
!
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
!    end subroutine set_source

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

      ! 1-way coupling:
      if (have_option('/embedded_models/fsi_model/one_way_coupling')) then

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
      type(vector_field), pointer :: fluid_position, solid_force
      type(vector_field) :: solid_position, solidforce_mesh
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
         solid_position = read_triangle_serial(mesh_filename, quad_degree=quad_degree)
!         solid_position = read_mesh_files(mesh_filename, quad_degree=quad_degree, solidmesh=.true., format=meshformat)

         ! Set mesh parameters
         solid_mesh = solid_position%mesh
         solid_mesh%name = trim(mesh_name)
         ! Now copy those back to the solid_position field:
         solid_position%mesh = solid_mesh
         solid_position%name = trim(solid_mesh%name)//'SolidCoordinate'

         ! Insert solid_mesh and solid_position into state:
         call insert(state, solid_mesh, trim(solid_mesh%name))
         call insert(state, solid_position, trim(solid_position%name))

         ! Also set-up solidvolumefraction field for all solids:
         call allocate(alpha_solidmesh, fluid_position%mesh, trim(solid_mesh%name)//"SolidConcentration")
         alpha_solidmesh%option_path = alpha_global%option_path
         call zero(alpha_solidmesh)
         ! And insert it into state so that we have one solidconcentration field per solidmesh:
         call insert(state, alpha_solidmesh, trim(solid_mesh%name)//"SolidConcentration")
         call deallocate(alpha_solidmesh)

         ! And the solidforce:
         call allocate(solidforce_mesh, solid_force%dim, solid_force%mesh, trim(solid_mesh%name)//"SolidForce")
         solidforce_mesh%option_path = solid_force%option_path
         call zero(solidforce_mesh)
         call insert(state, solidforce_mesh, trim(solid_mesh%name)//"SolidForce")
         call deallocate(solidforce_mesh)

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

    subroutine fsi_model_check_options
       integer :: ndim
       ! Get dimension:
       call get_option('/geometry/dimension', ndim)
       ! Check options for Implicit Solids:
       if (have_option('/embedded_models/fsi_model/one_way_coupling') .and. ndim==1) then 
          ewrite(-1,*) 'Error: The 1-way Fluid-Structure Interactions are not supported for 1D simulations via the FSI Model'
          FLExit('Use a 2D or 3D set-up when using a 1-way coupled simulation')
       end if

       if (.not. have_option('/material_phase[0]/vector_field::Velocity/prognostic/vector_field::Absorption/diagnostic')) then
          ewrite(-1,*) 'This models relies on an internal algorithm to compute the absorption. Please enable the vector_field "Absorption"'
          ewrite(-1,*) 'under vector_field "Veclotiy" in your options file'
       end if

       if (.not. have_option('/material_phase[0]/scalar_field::SolidConcentration/diagnostic')) then
          ewrite(-1,*) 'This models relies on an internal algorithm to compute the solid volume fraction on the fluid mesh.'
          ewrite(-1,*) 'Please enable the scalar_field "SolidConcentration" under "material_phase" in your options file'
       end if
       
       if (.not. have_option('/material_phase[0]/vector_field::SolidForce')) then
          ewrite(-1,*) 'This models relies on an internal algorithm to compute the force acting on the solid body.'
          ewrite(-1,*) 'Please enable the vector_field "SolidForce" under "material_phase" in your options file'
       end if

       if (have_option('/embedded_models/fsi_model/two_way_coupled')) then
          ewrite(-1,*) 'Error: The FSI Model does not support 2-way Fluid-Solid coupling'
          ewrite(-1,*) 'Use the Fluidity/FEMDEM approach for 2-way coupling instead'
          FLExit('2-way coupling of Fluids and Solids is not supported here. Use the Fluidity/FEMDEM approach for 2-way coupling instead (implicit_solids in the schema)')
       end if
    end subroutine fsi_model_check_options

  !----------------------------------------------------------------------------

end module fsi_model
