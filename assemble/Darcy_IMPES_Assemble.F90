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

module darcy_impes_assemble_module

   ! keep in this order, please:
   use quadrature
   use elements
   use sparse_tools
   use fields
   !
   use field_derivatives
   use cv_shape_functions
   use cv_faces
   use cvtools
   use cv_fields
   use cv_upwind_values
   use cv_face_values
   use cv_options
   use boundary_conditions
   use fefields, only: compute_cv_mass
   use petsc_solve_state_module
   use transform_elements, only: transform_cvsurf_to_physical, &
                                 transform_cvsurf_facet_to_physical
   use state_module
   use fldebug
   use spud
   use porous_media
   use parallel_tools
   use adaptive_timestepping
   use signal_vars, only : SIG_INT
   use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN

   implicit none

   private

   public :: darcy_impes_type, &
             darcy_impes_cv_options_type, &
             darcy_impes_copy_to_old, &
             darcy_impes_calculate_gradient_pressures, &
             darcy_impes_calculate_non_first_phase_pressures, &
             darcy_impes_calculate_phase_one_saturation_diagnostic, &
             darcy_impes_calculate_vel_mob_ff_and_cfl_fields, &
             darcy_impes_calculate_sum_saturation, &
             darcy_impes_calculate_densities, &
             darcy_impes_calculate_cflnumber_field_based_dt, &
             darcy_impes_initialise_cached_face_value, &
             darcy_impes_calculate_inverse_characteristic_length, &
             darcy_impes_calculate_relperm_den_first_face_values
   
   ! Options associated with the CV discretisation for the darcy impes solver
   type darcy_impes_cv_options_type
      ! what CV_FACEVALUE_ to use
      integer :: facevalue
      ! whether to limit the face value spatially
      logical :: limit_facevalue
      ! what CV_LIMITER_ to use
      integer :: limiter
      ! the slopes of the limiter being used (if Sweby)
      real, dimension(2) :: limiter_slopes
      ! what upwind scheme are we using (if any)
      integer :: upwind_scheme   
   end type darcy_impes_cv_options_type
   
   ! Options associated with explicit advection subcycling for CV scalar field solver
   type darcy_impes_subcycle_options_type
      logical :: have
      integer :: number
   end type darcy_impes_subcycle_options_type
   
   ! Options associated with adaptive time stepping
   type darcy_impes_adaptive_dt_options_type
      logical :: have
      real    :: requested_cfl
      real    :: min_dt
      real    :: max_dt
      real    :: increase_tolerance
      logical :: min_dt_terminate_if_reached
      logical :: at_first_dt
   end type darcy_impes_adaptive_dt_options_type
   
   ! Options associated with EoS for a phase
   type darcy_impes_eos_options_type
      logical :: have_fluids_linear
      real    :: fluids_linear_reference_density
   end type darcy_impes_eos_options_type
   
   ! Data associated with the cached CV face values
   type cached_face_value_type
      real, dimension(:,:,:,:), pointer :: relperm      ! nsub, ngi_vele, vele, nphase
      real, dimension(:,:,:,:), pointer :: relperm_bdy  ! nsub, ngi_sele, sele, nphase
      real, dimension(:,:,:),   pointer :: den          !       ngi_vele, vele, nphase
      real, dimension(:,:,:),   pointer :: den_bdy      !       ngi_sele, sele, nphase
      real, dimension(:,:),     pointer :: detwei       !       ngi_vele, vele
      real, dimension(:,:),     pointer :: detwei_bdy   !       ngi_sele, sele
      real, dimension(:,:,:),   pointer :: normal       ! ndim, ngi_vele, vele
      real, dimension(:,:,:),   pointer :: normal_bdy   ! ndim, ngi_sele, sele
      real, dimension(:,:,:,:), pointer :: p_dshape     ! nloc, ngi_vele, ndim, vele
!!!!! *** THIS IS NOT POSSIBLE YET ***
      real, dimension(:,:,:,:), pointer :: p_dshape_bdy ! nloc, ngi_sele, ndim, sele
      logical                           :: cached_detwei_normal
      logical                           :: cached_p_dshape
   end type cached_face_value_type
   
   type darcy_impes_type
      ! *** Pointers to fields from state that have array length of number of phases ***
      type(vector_field_pointer), dimension(:), pointer :: darcy_velocity
      type(vector_field_pointer), dimension(:), pointer :: old_darcy_velocity
      type(scalar_field_pointer), dimension(:), pointer :: mobility      
      type(scalar_field_pointer), dimension(:), pointer :: fractional_flow
      type(scalar_field_pointer), dimension(:), pointer :: saturation
      type(scalar_field_pointer), dimension(:), pointer :: old_saturation
      type(scalar_field_pointer), dimension(:), pointer :: saturation_source
      type(scalar_field_pointer), dimension(:), pointer :: old_saturation_source
      type(scalar_field_pointer), dimension(:), pointer :: relative_permeability
      type(scalar_field_pointer), dimension(:), pointer :: old_relative_permeability
      type(scalar_field_pointer), dimension(:), pointer :: viscosity
      type(scalar_field_pointer), dimension(:), pointer :: old_viscosity
      type(scalar_field_pointer), dimension(:), pointer :: cfl
      type(scalar_field_pointer), dimension(:), pointer :: old_cfl
      type(scalar_field_pointer), dimension(:), pointer :: pressure
      type(scalar_field_pointer), dimension(:), pointer :: old_pressure
      type(scalar_field_pointer), dimension(:), pointer :: capilliary_pressure
      type(scalar_field_pointer), dimension(:), pointer :: old_capilliary_pressure
      type(vector_field_pointer), dimension(:), pointer :: gradient_pressure
      type(vector_field_pointer), dimension(:), pointer :: old_gradient_pressure
      type(vector_field_pointer), dimension(:), pointer :: iterated_gradient_pressure      
      type(scalar_field_pointer), dimension(:), pointer :: density
      type(scalar_field_pointer), dimension(:), pointer :: old_density
      ! *** Pointers to fields from state that are NOT phase dependent ***
      type(mesh_type),    pointer :: pressure_mesh
      type(mesh_type),    pointer :: gradient_pressure_mesh
      type(mesh_type),    pointer :: velocity_mesh      
      type(mesh_type),    pointer :: elementwise_mesh
      type(scalar_field), pointer :: average_pressure
      type(scalar_field), pointer :: porosity
      type(scalar_field), pointer :: old_porosity
      type(scalar_field), pointer :: absolute_permeability
      type(scalar_field), pointer :: old_absolute_permeability
      type(vector_field), pointer :: positions
      type(vector_field), pointer :: total_darcy_velocity
      type(vector_field), pointer :: total_mobility
      type(scalar_field), pointer :: sum_saturation
      type(scalar_field), pointer :: div_total_darcy_velocity
      type(vector_field), pointer :: gravity_direction
      ! *** Pointer to the pressure mesh - pressure mesh sparsity, used for pressure matrix and finding CV upwind values ***
      type(csr_sparsity), pointer :: sparsity_pmesh_pmesh
      ! *** The gravity magnitude ***
      real :: gravity_magnitude
      ! *** The full gravity field - direction * magnitude, allocated here ***
      type(vector_field) :: gravity
      type(vector_field) :: old_gravity
      ! *** Fields allocated here used in assemble algorithm ***
      type(vector_field) :: positions_pressure_mesh
      type(csr_matrix)   :: pressure_matrix
      type(scalar_field) :: lhs
      type(scalar_field) :: rhs
      type(scalar_field) :: rhs_adv
      type(scalar_field) :: rhs_time
      type(scalar_field) :: inverse_cv_mass_velocity_mesh
      type(scalar_field) :: inverse_cv_mass_pressure_mesh
      type(scalar_field) :: cv_mass_pressure_mesh_with_porosity   
      type(scalar_field) :: cv_mass_pressure_mesh_with_old_porosity 
      type(scalar_field) :: cv_mass_pressure_mesh
      type(scalar_field_pointer), dimension(:), pointer :: old_relperm_subcycle
      type(scalar_field_pointer), dimension(:), pointer :: old_saturation_subcycle
      ! *** Data associated with v, pressure and saturation BC allocated here ***
      type(mesh_type)                          :: bc_surface_mesh
      type(scalar_field)                       :: v_bc_value
      integer,           dimension(:), pointer :: v_bc_flag
      type(scalar_field)                       :: pressure_bc_value
      integer,           dimension(:), pointer :: pressure_bc_flag
      type(scalar_field)                       :: saturation_bc_value
      integer,           dimension(:), pointer :: saturation_bc_flag
      type(scalar_field)                       :: inverse_characteristic_length
      real                                     :: weak_pressure_bc_coeff
      ! *** The number of phase and the CV surface quadrature degree to use ***
      integer :: number_phase, cv_surface_quaddegree   
      ! *** Data specifically associated with the CV discretisation allocated here ***
      type(darcy_impes_cv_options_type)       :: old_relperm_cv_options
      type(darcy_impes_cv_options_type)       :: density_cv_options
      type(cv_faces_type)                     :: cvfaces
      type(element_type)                      :: x_cvshape_full
      type(element_type)                      :: p_cvshape_full
      type(element_type)                      :: gradp_cvshape_full
      type(element_type)                      :: x_cvshape
      type(element_type)                      :: p_cvshape
      type(element_type)                      :: gradp_cvshape
      type(element_type)                      :: x_cvbdyshape
      type(element_type)                      :: p_cvbdyshape      
      type(element_type)                      :: gradp_cvbdyshape
      type(darcy_impes_subcycle_options_type) :: subcy_opt_sat
      type(csr_matrix)                        :: relperm_upwind
      type(csr_matrix)                        :: density_upwind
      ! *** Flag for whether the first phase saturation is diagnostic, else it is prognostic ***
      logical                                 :: phase_one_saturation_diagnostic      
      ! *** Flag for whether the first phase pressure is prognostic, else it is prescribed *** 
      logical :: first_phase_pressure_prognostic
      ! *** The cached face values for all phases ***
      type(cached_face_value_type) :: cached_face_value
      ! *** The advection subcycle timestep size ***
      real    :: dt_subcycle
      ! *** Time time step size, also stored here for convenience ***
      real :: dt
      ! *** The current time, also stored here for convenience ***
      real :: current_time 
      ! *** Non linear iteration, also stored here for convencience ***
      integer :: nonlinear_iter
      ! *** Max Non linear iteration, also stored here for convencience ***
      integer :: max_nonlinear_iter
      ! *** Max Non linear iteration first timestep after adapt, also stored here for convencience ***
      integer :: max_nonlinear_iter_first_timestep_after_adapt
      ! *** Max Non linear iteration for this timestep, also stored here for convencience ***
      integer :: max_nonlinear_iter_this_timestep
      ! *** The number of volume and surface elements, stored here for convenience ***
      integer :: number_vele
      integer :: number_sele
      ! *** Geometric dimension, also stored here for convenience ***
      integer :: ndim
      ! *** Options data associated with adaptive time stepping stored here ***
      type(darcy_impes_adaptive_dt_options_type) :: adaptive_dt_options
      ! *** Options data associated with each phase EoS stored here *** 
      type(darcy_impes_eos_options_type), dimension(:), pointer :: eos_options
      ! *** Pointer to main state array, for convenience ***
      type(state_type), dimension(:), pointer :: state
   end type darcy_impes_type
   
   contains

! ----------------------------------------------------------------------------
  
   subroutine darcy_impes_copy_to_old(di)
            
      !!< Copy to old darcy impes data not associated with di%state

      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p
      
      call set(di%old_gravity, di%gravity)
      
      call set(di%cv_mass_pressure_mesh_with_old_porosity, di%cv_mass_pressure_mesh_with_porosity)
      
   end subroutine darcy_impes_copy_to_old

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve_phase_pressures(di)
      
      !!< Assemble and solve the phase pressures. The first phase pressure 
      !!< is solved via a matrix equation (if prognostic) and the other
      !!< phases pressures are calculated from the first phase pressure
      !!< and their own capilliary pressures
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! Assemble and solve the first phase pressure if it is prognostic
      if (di%first_phase_pressure_prognostic) call darcy_impes_assemble_and_solve_first_phase_pressure(di)
      
      ! Calculate the non first phase pressure's
      call darcy_impes_calculate_non_first_phase_pressures(di)       
      
      ! Calculate the average pressure
      call darcy_impes_calculate_average_pressure(di)
      
   end subroutine darcy_impes_assemble_and_solve_phase_pressures

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve_first_phase_pressure(di)
      
      !!< Assemble and solve the first phase pressure. A source is included 
      !!< due to the capilliary pressures of non first phases as well as gravity 
      !!< and the rate of change of porosity and the individual phase sources.
      !!< Face values for relative permeability and density have already been evaluated.
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p, vele, sele, iloc, oloc, jloc, face, gi, ggi
      real    :: face_value, grad_cap_p_dot_n, g_dot_n
      real,    dimension(1)              :: absperm_ele, visc_ele
      real,    dimension(:,:),   pointer :: grav_ele
      real,    dimension(:,:),   pointer :: grad_cap_pressure_face_quad
      real,    dimension(:,:),   pointer :: x_ele
      real,    dimension(:,:,:), pointer :: p_dshape
      real,    dimension(:,:),   pointer :: normal
      real,    dimension(:),     pointer :: detwei
      real,    dimension(:),     pointer :: normgi
      logical, dimension(:),     pointer :: notvisited
      real,    dimension(:,:),   pointer :: p_mat_local
      real,    dimension(:),     pointer :: p_rhs_local
      real,    dimension(:,:),   pointer :: x_face_quad
      integer, dimension(:),     pointer :: x_pmesh_nodes
      integer, dimension(:),     pointer :: p_nodes      
      real,    dimension(1)              :: absperm_ele_bdy, visc_ele_bdy
      real,    dimension(:,:),   pointer :: grav_ele_bdy
!!!!! *** THIS IS NOT POSSIBLE YET ***
!!!!!      real,    dimension(:,:),   pointer :: grad_cap_pressure_face_quad_bdy
!!!!!      real,    dimension(:,:,:), pointer :: p_dshape_bdy
      real,    dimension(:),     pointer :: bc_sele_val
      real,    dimension(:),     pointer :: inv_char_len_ele_bdy
      real,    dimension(:,:),   pointer :: normal_bdy
      real,    dimension(:),     pointer :: detwei_bdy
      real,    dimension(:,:),   pointer :: x_ele_bdy
      real,    dimension(:),     pointer :: p_rhs_local_bdy
      real,    dimension(:,:),   pointer :: p_matrix_local_bdy
      integer, dimension(:),     pointer :: p_nodes_bdy
      integer, parameter :: PRESSURE_BC_TYPE_WEAKDIRICHLET   = 1, PRESSURE_BC_TYPE_DIRICHLET = 2 
      integer, parameter :: SATURATION_BC_TYPE_WEAKDIRICHLET = 1
      integer, parameter :: V_BC_TYPE_NORMAL_FLOW            = 1, V_BC_TYPE_NO_NORMAL_FLOW   = 2
            
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      if (.not. di%cached_face_value%cached_p_dshape) then
         allocate(p_dshape(ele_loc(di%pressure_mesh,1), di%p_cvshape%ngi, mesh_dim(di%pressure_mesh)))
      end if
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(normal(di%ndim,di%x_cvshape%ngi))
         allocate(detwei(di%x_cvshape%ngi))      
      end if
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(p_mat_local(ele_loc(di%pressure_mesh,1), ele_loc(di%pressure_mesh,1)))
      allocate(p_rhs_local(ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(grav_ele(di%ndim,1))
      allocate(grad_cap_pressure_face_quad(di%ndim, di%p_cvshape%ngi))

      allocate(grav_ele_bdy(di%ndim,1))
      allocate(bc_sele_val(face_loc(di%pressure_mesh,1)))
      allocate(inv_char_len_ele_bdy(face_loc(di%pressure_mesh,1)))
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(detwei_bdy(di%x_cvbdyshape%ngi))
         allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      end if
      allocate(p_rhs_local_bdy(face_loc(di%pressure_mesh,1)))
      allocate(p_matrix_local_bdy(face_loc(di%pressure_mesh,1),face_loc(di%pressure_mesh,1)))
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(p_nodes_bdy(face_loc(di%pressure_mesh,1)))
!!!!! *** THIS IS NOT POSSIBLE YET ***
!!!!!      allocate(grad_cap_pressure_face_quad_bdy(di%ndim, di%p_cvbdyshape_full%ngi))
!!!!!      if (.not. di%cached_face_value%cached_p_dshape) then
!!!!!         allocate(p_dshape_bdy(face_loc(di%pressure_mesh,1), di%p_cvbdyshape_full%ngi, mesh_dim(di%pressure_mesh)))
!!!!!      end if
      
      ewrite(1,*) 'Solve first phase Pressure'
      
      ! Initialise the matrix and rhs
      call zero(di%pressure_matrix)
      call zero(di%rhs)

      ewrite(1,*) 'Add phase sources to global continuity'
      
      src_phase_loop: do p = 1, di%number_phase
      
         call addto(di%rhs, di%saturation_source(p)%ptr)
      
      end do src_phase_loop
      
      ! Should this include the porosity ...?!?!
      call scale(di%rhs, di%cv_mass_pressure_mesh)
            
      ewrite(1,*) 'Add rate of change of porosity to global continuity equation'
      
      ! Add rate of change of porosity to rhs    
      call addto(di%rhs, di%cv_mass_pressure_mesh_with_old_porosity, scale = 1.0/di%dt)
                  
      call addto(di%rhs, di%cv_mass_pressure_mesh_with_porosity, scale = -1.0/di%dt)
      
      ! Assemble a contribution from each phase to form a global continuity equation to solve for first phase pressure
      phase_loop: do p = 1, di%number_phase
         
         ewrite(1,*) 'Assemble volume contribution to global continuity from phase ',p
            
         ! Loop volume elements assembling local contributions     
         vol_element_loop: do vele = 1, di%number_vele

            ! The node indices of the pressure mesh
            p_nodes => ele_nodes(di%pressure_mesh, vele)
            
            ! The gravity values for this element for each direction
            grav_ele = ele_val(di%gravity, vele) 

            ! get the viscosity value for this element for this phase
            visc_ele = ele_val(di%viscosity(p)%ptr, vele)

            ! get the absolute permeability value for this element
            absperm_ele = ele_val(di%absolute_permeability, vele)         
            
            ! obtain the transformed determinant*weight and normals
            if (di%cached_face_value%cached_detwei_normal) then
               
               detwei => di%cached_face_value%detwei(:,vele)
               
               normal => di%cached_face_value%normal(:,:,vele)
               
            else

               ! get the coordinate values for this element for each positions local node
               x_ele = ele_val(di%positions, vele)         

               ! get the coordinate values for this element for each quadrature point
               x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)

               ! The node indices of the positions projected to the pressure mesh
               x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)
            
               call transform_cvsurf_to_physical(x_ele, di%x_cvshape, detwei, normal, di%cvfaces)
            
            end if
            
            ! obtain the derivative of the pressure mesh shape function at the CV face quadrature points
            if (di%cached_face_value%cached_p_dshape) then
               
               p_dshape => di%cached_face_value%p_dshape(:,:,:,vele)
               
            else
            
               call transform_to_physical(di%positions, vele, x_shape = di%x_cvshape_full, &
                                          shape = di%p_cvshape_full, dshape = p_dshape)
            
            end if
            
            ! get the gradient capilliary pressure at the cv surface quadrature points for each direction 
            if (p > 1) then
               grad_cap_pressure_face_quad = &
              &darcy_impes_ele_grad_at_quad_scalar(di%capilliary_pressure(p)%ptr, vele, dn = p_dshape)
            end if
            
            ! Initialise array for the quadrature points of this 
            ! element for whether it has already been visited
            notvisited = .true.

            ! Initialise the local p matrix and rhs to assemble for this element
            p_mat_local = 0.0
            p_rhs_local = 0.0
            
            ! loop over local nodes within this element
            nodal_loop_i: do iloc = 1, di%pressure_mesh%shape%loc

               ! loop over CV faces internal to this element
               face_loop: do face = 1, di%cvfaces%faces

                  ! is this a face neighbouring iloc?
                  is_neigh: if(di%cvfaces%neiloc(iloc, face) /= 0) then

                     ! find the opposing local node across the CV face
                     oloc = di%cvfaces%neiloc(iloc, face)

                     ! loop over gauss points on face
                     quadrature_loop: do gi = 1, di%cvfaces%shape%ngi

                        ! global gauss pt index for this element
                        ggi = (face-1)*di%cvfaces%shape%ngi + gi

                        ! check if this quadrature point has already been visited
                        check_visited: if(notvisited(ggi)) then

                           notvisited(ggi) = .false.

                           if (di%cached_face_value%cached_detwei_normal) then

                              ! cached normal has correct orientation already
                              normgi = normal(:,ggi)

                           else

                              ! correct the orientation of the normal so it points away from iloc
                              normgi = orientate_cvsurf_normgi(node_val(di%positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                              &x_face_quad(:,ggi), normal(:,ggi))

                           end if

                           ! Form the face value = detwei * (relperm*absperm/visc)
                           face_value = detwei(ggi) * &
                                        sum(di%cached_face_value%relperm(:,ggi,vele,p)) * &
                                        absperm_ele(1) / &
                                        visc_ele(1) 
   
                           ! Form the local matrix given by - n_i . sum_{phase} ( relperm*absperm/visc ) dP_1/dx_j
                           do jloc = 1,di%pressure_mesh%shape%loc
                              
                              p_mat_local(iloc,jloc) = p_mat_local(iloc,jloc) - &
                                                       sum(p_dshape(jloc, ggi, :)*normgi, 1)*face_value

                              p_mat_local(oloc,jloc) = p_mat_local(oloc,jloc) - &
                                                       sum(p_dshape(jloc, ggi, :)*(-normgi), 1)*face_value

                           end do

                           ! Add gravity term to rhs = - n_i . sum_{phase} ( relperm*absperm/visc ) * den*grav

                           ! Find g dot n
                           g_dot_n = dot_product(grav_ele(:,1), normgi)

                           p_rhs_local(iloc)  = p_rhs_local(iloc) - &
                                                face_value * &
                                                di%cached_face_value%den(ggi,vele,p) * &
                                                g_dot_n

                           ! Add capilliary pressure term to rhs = n_i . sum_{phase} ( relperm*absperm/visc ) dP_c/dx_j
                           ! only for phase > 1
                           if (p > 1) then

                              ! Find grad_P_c dot n
                              grad_cap_p_dot_n  = dot_product(grad_cap_pressure_face_quad(:,ggi), normgi)

                              p_rhs_local(iloc) = p_rhs_local(iloc) + &
                                                  face_value * &
                                                  grad_cap_p_dot_n

                           end if

                        end if check_visited

                     end do quadrature_loop

                  end if is_neigh

               end do face_loop

            end do nodal_loop_i

            ! Add volume element contribution to global pressure matrix and rhs
            call addto(di%pressure_matrix, p_nodes, p_nodes, p_mat_local)
            
            call addto(di%rhs, p_nodes, p_rhs_local)

         end do vol_element_loop
         
      end do phase_loop
                  
      ! Add normal v BC integrals for each phase, else include if required 
      ! a weak pressure BC (which if not given is assumed zero).

      phase_loop_bc: do p = 1, di%number_phase

         ewrite(1,*) 'Assemble boundary contribution to global continuity from phase ',p
         
         ! Get the phase pressure BC - if no v and weak pressure then extra integrals are added
         ! For phase 1 also get strong pressure BC as this implies no need to add any surface integrals
         if (p == 1) then

            call get_entire_saturation_or_pressure_boundary_condition(di%pressure(p)%ptr, &
                                                                      (/"weakdirichlet", &
                                                                        "dirichlet    "/), &
                                                                      di%bc_surface_mesh, &
                                                                      di%pressure_bc_value, &
                                                                      di%pressure_bc_flag)
         
         else
         
            call get_entire_saturation_or_pressure_boundary_condition(di%pressure(p)%ptr, &
                                                                      (/"weakdirichlet"/), &
                                                                      di%bc_surface_mesh, &
                                                                      di%pressure_bc_value, &
                                                                      di%pressure_bc_flag)
         
         end if
         
         ! Get this phase v BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
         call get_v_boundary_condition(di%darcy_velocity(p)%ptr, &
                                       (/"normal_flow   ", &
                                         "no_normal_flow"/), &
                                       di%bc_surface_mesh, &
                                       di%v_bc_value, &
                                       di%v_bc_flag)
                  
         sele_loop: do sele = 1, di%number_sele
            
            ! A no_normal_flow BC adds nothing to the matrix and rhs so cycle sele loop
            if (di%v_bc_flag(sele) == V_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop
            
            ! If phase 1 and strong pressure BC then no need to add integrals so cycle sele loop
            if ((p == 1) .and. (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_DIRICHLET)) cycle sele_loop 
            
            p_nodes_bdy = face_global_nodes(di%pressure_mesh, sele)

            ! obtain the transformed determinant*weight and normals
            if (di%cached_face_value%cached_detwei_normal) then

               detwei_bdy => di%cached_face_value%detwei_bdy(:,sele)
               
               normal_bdy => di%cached_face_value%normal_bdy(:,:,sele)
            
            else

               x_ele     = ele_val(di%positions, face_ele(di%positions, sele))
               x_ele_bdy = face_val(di%positions, sele)
            
               call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)
            
            end if
            
            if (di%v_bc_flag(sele) == V_BC_TYPE_NORMAL_FLOW) then
                              
               bc_sele_val = ele_val(di%v_bc_value, sele)
                        
            else
            
               ! get the viscosity value for this sele for this phase
               visc_ele_bdy = face_val(di%viscosity(p)%ptr, sele)

               ! get the absolute permeability value for this sele
               absperm_ele_bdy = face_val(di%absolute_permeability, sele)         
               
               ! the inverse characteristic length used for pressure weak bc
               inv_char_len_ele_bdy = ele_val(di%inverse_characteristic_length, sele)
                              
               ! Get the pressure BC value if required
               if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then
               
                  bc_sele_val = ele_val(di%pressure_bc_value, sele)
               
               end if

               ! The gravity values for this element for each direction
               grav_ele_bdy = face_val(di%gravity, sele) 

!!!!! *** THIS IS NOT POSSIBLE YET ***
!!!!!               ! get the gradient capilliary pressure at the cv surface quadrature points for each direction 
!!!!!               if (p > 1) then
!!!!!                  grad_cap_pressure_face_quad_bdy = &
!!!!!                 &darcy_impes_face_grad_at_quad_scalar(di%capilliary_pressure(p)%ptr, sele, dn = p_dshape)
!!!!!               end if

!!!!!               ! obtain the derivative of the pressure mesh shape function at the CV face quadrature points
!!!!!               if (di%cached_face_value%cached_p_dshape) then

!!!!!                  p_dshape_bdy => di%cached_face_value%p_dshape_bdy(:,:,:,sele)

!!!!!               else

!!!!!                  call transform_to_physical(di%positions, sele, x_shape = di%x_cvbdyshape_full, &
!!!!!                                             shape = di%p_cvbdyshape_full, dshape = p_dshape_bdy)

!!!!!               end if
               
            end if 
            
            p_matrix_local_bdy = 0.0
            p_rhs_local_bdy    = 0.0

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi
                                                
                        if (di%v_bc_flag(sele) == V_BC_TYPE_NORMAL_FLOW) then
                           
                           ! If have normal_flow BC for this phase then include CV integral of value*number_subcycle
                                                      
                           p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - &
                                                   di%subcy_opt_sat%number * &
                                                   bc_sele_val(iloc) * &
                                                   detwei_bdy(ggi)
                        
                        else 
                                                                                 
                           ! Form the face value = detwei * (relperm*absperm/visc)
                           face_value = detwei_bdy(ggi) * &
                                        sum(di%cached_face_value%relperm_bdy(:,ggi,sele,p)) * &
                                        absperm_ele_bdy(1) / &
                                        visc_ele_bdy(1) 
                           
                           ! Add gravity term to rhs = - n_i . sum_{phase} ( relperm*absperm/visc ) * den*grav

                           ! Find g dot n
                           g_dot_n = dot_product(grav_ele_bdy(:,1), normal_bdy(:,ggi))
                           
                           p_rhs_local_bdy(iloc)  = p_rhs_local_bdy(iloc) - &
                                                    face_value * &
                                                    di%cached_face_value%den_bdy(ggi,sele,p) * &
                                                    g_dot_n

!!!!! *** THIS IS NOT POSSIBLE YET ***
!!!!!                           ! Add capilliary pressure term to rhs = n_i . sum_{phase} ( relperm*absperm/visc ) dP_c/dx_j
!!!!!                           ! only for phase > 1
!!!!!                           if (p > 1) then

!!!!!                              ! Find grad_P_c dot n
!!!!!                              grad_cap_p_dot_n  = dot_product(grad_cap_pressure_face_quad_bdy(:,ggi), normal_bdy(:,ggi))

!!!!!                              p_rhs_local(iloc) = p_rhs_local(iloc) + &
!!!!!                                                  face_value * &
!!!!!                                                  grad_cap_p_dot_n

!!!!!                           end if
                           
!!!!!                           ! Add implicit gradient of pressure phase 1 to matrix with coeff 
!!!!!                           ! - n_i . sum_{phase} ( relperm*absperm/visc ) 

                           
                           ! If have phase weak pressure BC then include integral of BC value in rhs
                           ! and an implicit surface mass matrix term
                           if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then
                                                            
                              do jloc = 1,di%pressure_mesh%faces%shape%loc 
                              
                                 p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - &
                                                         inv_char_len_ele_bdy(jloc) * &
                                                         di%p_cvbdyshape%n(jloc,ggi) * &
                                                         bc_sele_val(jloc) * &
                                                         sum(normal_bdy(:,ggi)) * &
                                                         face_value * &
                                                         di%weak_pressure_bc_coeff
                                                            
                              end do
                                                         
                              do jloc = 1,di%pressure_mesh%faces%shape%loc 

                                 p_matrix_local_bdy(iloc,jloc) = p_matrix_local_bdy(iloc,jloc) - &
                                                                 inv_char_len_ele_bdy(jloc) * &
                                                                 di%p_cvbdyshape%n(jloc,ggi) * &
                                                                 sum(normal_bdy(:,ggi)) * &
                                                                 face_value * &
                                                                 di%weak_pressure_bc_coeff

                              end do

                           end if 

                        end if
                        
                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop
            
            call addto(di%pressure_matrix, p_nodes_bdy, p_nodes_bdy, p_matrix_local_bdy)

            call addto(di%rhs, p_nodes_bdy, p_rhs_local_bdy)

         end do sele_loop      

      end do phase_loop_bc
            
      ! Apply any strong dirichlet BC's - that can only be applied to the first phase pressure
      call apply_dirichlet_conditions(di%pressure_matrix, di%rhs, di%pressure(1)%ptr)
      
      ! Solve the pressure
      call petsc_solve(di%pressure(1)%ptr, di%pressure_matrix, di%rhs, di%state(1))
      
      ! Set the strong BC nodes to the values to be consistent
      call set_dirichlet_consistent(di%pressure(1)%ptr) 
           
      ! deallocate local variables as required
      deallocate(x_ele)
      if (di%cached_face_value%cached_p_dshape) then
         nullify(p_dshape)      
      else
         deallocate(p_dshape)
      end if
      if (di%cached_face_value%cached_detwei_normal) then
         nullify(detwei)
         nullify(normal)      
      else
         deallocate(detwei)
         deallocate(normal)
      end if
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(p_mat_local)
      deallocate(p_rhs_local)
      deallocate(x_face_quad)
      deallocate(grav_ele)
      deallocate(grad_cap_pressure_face_quad)

      deallocate(grav_ele_bdy)
      deallocate(bc_sele_val)
      deallocate(inv_char_len_ele_bdy)
      if (di%cached_face_value%cached_detwei_normal) then
         nullify(detwei_bdy)
         nullify(normal_bdy)      
      else
         deallocate(detwei_bdy)
         deallocate(normal_bdy)
      end if
      deallocate(p_rhs_local_bdy)
      deallocate(p_matrix_local_bdy)
      deallocate(x_ele_bdy)      
      deallocate(p_nodes_bdy)
!!!!! *** THIS IS NOT POSSIBLE YET ***
!!!!!      deallocate(grad_cap_pressure_face_quad_bdy)
      
      ewrite(1,*) 'Finished solve first phase Pressure'
      
      ewrite_minmax(di%pressure(1)%ptr)
      
   end subroutine darcy_impes_assemble_and_solve_first_phase_pressure

! ----------------------------------------------------------------------------

   subroutine get_v_boundary_condition(v, &
                                       types, &
                                       bc_surface_mesh, &
                                       v_bc_value, &
                                       v_bc_flag)

      !!< Form the data associated with any v BC by returning 
      !!< full surface mesh arrays of a flag indicating  either 
      !!< no_normal_flow or normal_flow and the value.

      type(vector_field),               intent(in),   target :: v
      character(len=*),   dimension(:), intent(in)           :: types
      type(mesh_type),                  intent(in)           :: bc_surface_mesh
      type(scalar_field),               intent(inout)        :: v_bc_value
      integer,            dimension(:), intent(inout)        :: v_bc_flag

      ! Local variables
      character(len=FIELD_NAME_LEN)                       :: bctype
      type(scalar_field),                         pointer :: scalar_surface_field
      integer,                      dimension(:), pointer :: surface_element_list
      integer                                             :: i, j, k, sele

      ewrite(1,*) 'Get DarcyVelocity boundary condition data'

      ! Zero the normal_flow value surface field on whole boundary mesh  
      call zero(v_bc_value)

      ! Initialise flag for whether surface element has normal_flow BC
      v_bc_flag = 0

      ! Loop each BC object instance for the v
      ! (May have multiple normal_flow BC's applied to different surface id's)
      BC_loop: do i=1, get_boundary_condition_count(v)

         ! Get this BC info
         call get_boundary_condition(v, i, type = bctype, &
                                     surface_element_list = surface_element_list)

         ! check this is a normal_flow BC  or no_normal_flow (nothing else is permitted)
         if ((trim(bctype) /= 'normal_flow') .and. (trim(bctype) /= 'no_normal_flow')) then
            FLAbort('Have unknown BC type for a DarcyVelocity')
         end if

         ! Extract the scalar_surface_field for this BC for normal_flow
         if (trim(bctype) == 'normal_flow') then
            if (associated(v%bc%boundary_condition(i)%scalar_surface_fields)) then
               scalar_surface_field => v%bc%boundary_condition(i)%scalar_surface_fields(1)
            else
               FLAbort('Component scalar_surface_fields for DarcyVelocity BC type not associated')
            end if
         end if 
         
         ! Loop the surface elements associated with this BC instance
         ! and place the required BC value in whole boundary field
         ! for the normal_flow BC type
         BC_sele_loop: do k = 1, size(surface_element_list)

            ! Find the whole domain surface element number
            sele = surface_element_list(k)

            ! Check that there is only 1 BC applied per surface element
            if (v_bc_flag(sele) /= 0) then             
               FLExit('Cannot apply more than 1 BC to a surface element for DarcyVelocity')
            end if

            ! Set the sele flag to indicate the BC type
            do j = 1, size(types)
               if (trim(types(j)) == trim(bctype)) exit
            end do
            if (j > size(types)) then
               FLAbort('Cannot find no_normal_flow or normal_flow bctype for DarcyVelocity')
            end if
            
            v_bc_flag(sele) = j

            ! Set the normal_flow field values from this BC for its sele
            if (trim(bctype) == 'normal_flow') then
               call set(v_bc_value, &
                        ele_nodes(bc_surface_mesh, sele), &
                        ele_val(scalar_surface_field, k))
            end if
            
         end do BC_sele_loop

      end do BC_loop

      ewrite_minmax(v_bc_value)

      ewrite(1,*) 'Finished get DarcyVelocity boundary condition data'

   end subroutine get_v_boundary_condition
  
! ----------------------------------------------------------------------------

   subroutine get_entire_saturation_or_pressure_boundary_condition(sfield, &
                                                                   types, &
                                                                   bc_surface_mesh, &
                                                                   sfield_bc_value, &
                                                                   sfield_bc_flag)
      
      !!< Get the entire BC data for either the saturation or pressure field
      
      type(scalar_field),               intent(in),   target :: sfield
      character(len=*),   dimension(:), intent(in)           :: types
      type(mesh_type),                  intent(in)           :: bc_surface_mesh
      type(scalar_field),               intent(inout)        :: sfield_bc_value
      integer,            dimension(:), intent(inout)        :: sfield_bc_flag
      
      ! local variables
      character(len=FIELD_NAME_LEN)                       :: bctype
      type(scalar_field),                         pointer :: surface_field
      integer,                      dimension(:), pointer :: surface_element_list
      integer                                             :: i, j, k, sele

      ewrite(1,*) 'Get ',trim(sfield%name),' boundary condition data'
      
      ! zero the bc value field
      call zero(sfield_bc_value)
      
      ! Initialise the bc type list for each sele
      sfield_bc_flag = 0
      
      ! Loop each BC object for the sfield 
      ! (May have multiple BC's applied to different suface id's)
      BC_loop: do i=1, get_boundary_condition_count(sfield)
         
         ! Get this BC info
         call get_boundary_condition(sfield, i, type = bctype, &
                                     surface_element_list = surface_element_list)

         ! see if we're interested in this one, if not skip it
         do j = 1, size(types)
            if (trim(types(j)) == trim(bctype)) exit
         end do
         if (j > size(types)) cycle
         
         ! Extract the surface field for this BC
         if (associated(sfield%bc%boundary_condition(i)%surface_fields)) then
            ! extract 1st surface field
            surface_field => sfield%bc%boundary_condition(i)%surface_fields(1)
         else
            FLAbort('Component surface_fields for Saturation or Pressure BC type not associated')
         end if

         ! Loop the surface elements associated with this BC instance
         ! and place the required BC value in whole boundary field
         BC_sele_loop: do k = 1, size(surface_element_list)
            
            ! Find the whole domain surface element number
            sele = surface_element_list(k)
            
            ! Check that there is only 1 BC applied per surface element
            if (sfield_bc_flag(sele)/=0) then
               ewrite(0,*) 'Requested types:', types
               ewrite(0,*) 'Of these boundary condition types only one may be applied'
               ewrite(0,*) 'to each surface element.'
               ewrite(0,*) 'Surface element nr.:', sele
               ewrite(0,*) 'has types', types(sfield_bc_flag(sele)), bctype
               ewrite(0,*) 'on field: ', sfield%name
               FLAbort('Issue with BC data for Saturation or Pressure')
            end if
            
            ! Set the sele flag to indicate the BC type
            sfield_bc_flag(sele) = j
            
            ! Set the BC sfield values from this BC for its sele
            call set(sfield_bc_value, &
                     ele_nodes(bc_surface_mesh, sele), &
                     ele_val(surface_field, k))

         end do BC_sele_loop
         
      end do BC_loop      

      ewrite_minmax(sfield_bc_value)
      
      ewrite(1,*) 'Finished get ',trim(sfield%name),' boundary condition data'
      
   end subroutine get_entire_saturation_or_pressure_boundary_condition

! ----------------------------------------------------------------------------

   subroutine darcy_impes_calculate_non_first_phase_pressures(di)
      
      !!< Calcuate the non first phase pressures which = first_phase_pressure + capilliary pressure
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p

      ewrite(1,*) 'Calculate Non first phase Pressures'
      
      phase_loop: do p = 2, di%number_phase
      
         call set(di%pressure(p)%ptr, di%pressure(1)%ptr)
         
         call addto(di%pressure(p)%ptr, di%capilliary_pressure(p)%ptr)
         
         ewrite_minmax(di%pressure(p)%ptr) 
         
      end do phase_loop

      ewrite(1,*) 'Finished Calculate Non first phase Pressures'

   end subroutine darcy_impes_calculate_non_first_phase_pressures

! ----------------------------------------------------------------------------

   subroutine darcy_impes_calculate_average_pressure(di)
      
      !!< Calcuate the average phase pressure
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p

      ewrite(1,*) 'Calculate AveragePressure'
      
      call zero(di%average_pressure)
      
      phase_loop: do p = 1, di%number_phase
               
         call addto(di%average_pressure, di%pressure(p)%ptr)
                  
      end do phase_loop
      
      call scale(di%average_pressure, 1.0/real(di%number_phase))
      
      ewrite_minmax(di%average_pressure) 

      ewrite(1,*) 'Finished Calculate AveragePressure'

   end subroutine darcy_impes_calculate_average_pressure

! ----------------------------------------------------------------------------

   subroutine darcy_impes_calculate_gradient_pressures(di)
      
      !!< Calculate the gradient pressures for each phase
      
      ! This routine requires optimisation and would become 
      ! not required if ele_grad_at_quad is used for pressure 
      ! and a face_grad_at_quad is made available

      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p
      
      ewrite(1,*) 'Calculate GradientPressure for each phase'
      
      phase_loop: do p = 1, di%number_phase
      
         call grad(di%pressure(p)%ptr, &
                   di%positions, &
                   di%gradient_pressure(p)%ptr)

         ewrite_minmax(di%gradient_pressure(p)%ptr)
      
      end do phase_loop

      ewrite(1,*) 'Finished Calculate GradientPressure for each phase'
      
   end subroutine darcy_impes_calculate_gradient_pressures

! ----------------------------------------------------------------------------

   subroutine darcy_impes_calculate_divergence_total_darcy_velocity(di)
      
      !!< Calculate the divergence of the total divergence velocity, 
      !!< which may include terms summed over subcycles

      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: vele, p, iloc, oloc, jloc, face, gi, ggi, sele, dim
      real    :: darcy_vel_face_value_dot_n, v_over_relperm_dot_n, face_value
      real,    dimension(1)            :: absperm_ele, visc_ele
      real,    dimension(:,:), pointer :: grad_pressure_face_quad
      real,    dimension(:,:), pointer :: grav_ele
      real,    dimension(:,:), pointer :: x_ele
      real,    dimension(:,:), pointer :: normal
      real,    dimension(:),   pointer :: detwei
      real,    dimension(:),   pointer :: normgi
      logical, dimension(:),   pointer :: notvisited
      real,    dimension(:),   pointer :: div_tvel_rhs_local
      real,    dimension(:,:), pointer :: x_face_quad
      integer, dimension(:),   pointer :: x_pmesh_nodes
      integer, dimension(:),   pointer :: p_nodes      
      real,    dimension(1)            :: absperm_ele_bdy, visc_ele_bdy
      real,    dimension(:,:), pointer :: grav_ele_bdy
      real,    dimension(:,:), pointer :: grad_pressure_face_quad_bdy
      real,    dimension(:,:), pointer :: v_over_relperm_face_quad_bdy
      real,    dimension(:),   pointer :: bc_sele_val
      real,    dimension(:),   pointer :: pressure_ele_bdy
      real,    dimension(:),   pointer :: inv_char_len_ele_bdy
      real,    dimension(:,:), pointer :: normal_bdy
      real,    dimension(:),   pointer :: detwei_bdy
      real,    dimension(:,:), pointer :: x_ele_bdy
      real,    dimension(:),   pointer :: div_tvel_rhs_local_bdy
      integer, dimension(:),   pointer :: p_nodes_bdy
      integer, parameter :: V_BC_TYPE_NORMAL_FLOW          = 1, V_BC_TYPE_NO_NORMAL_FLOW = 2
      integer, parameter :: PRESSURE_BC_TYPE_WEAKDIRICHLET = 1
                        
      ewrite(1,*) 'Calculate the DivergenceTotalDarcyVelocity'
      
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(normal(di%ndim,di%x_cvshape%ngi))
         allocate(detwei(di%x_cvshape%ngi))      
      end if
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(div_tvel_rhs_local(ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(grav_ele(di%ndim,1))
      
      allocate(grav_ele_bdy(di%ndim,1))
      allocate(grad_pressure_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      allocate(v_over_relperm_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      allocate(bc_sele_val(face_loc(di%pressure_mesh,1)))
      allocate(pressure_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(inv_char_len_ele_bdy(face_loc(di%pressure_mesh,1)))
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(detwei_bdy(di%x_cvbdyshape%ngi))
         allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      end if
      allocate(div_tvel_rhs_local_bdy(face_loc(di%pressure_mesh,1)))      
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(p_nodes_bdy(face_loc(di%pressure_mesh,1)))

      call zero(di%div_total_darcy_velocity)

      phase_loop: do p = 1, di%number_phase
                       
         vele_loop: do vele = 1, di%number_vele

            ! The node indices of the pressure field
            p_nodes => ele_nodes(di%pressure_mesh, vele)
            
            ! The gravity values for this element for each direction
            grav_ele = ele_val(di%gravity, vele) 

            ! get the viscosity value for this element for this phase
            visc_ele = ele_val(di%viscosity(p)%ptr, vele)

            ! get the absolute permeability value for this element
            absperm_ele = ele_val(di%absolute_permeability, vele)         

            ! get the latest gradient pressure at the cv surface quadrature points for each direction for this phase
            grad_pressure_face_quad = ele_val_at_quad(di%gradient_pressure(p)%ptr, vele, di%gradp_cvshape)

            ! obtain the transformed determinant*weight and normals
            if (di%cached_face_value%cached_detwei_normal) then
               
               detwei => di%cached_face_value%detwei(:,vele)
               
               normal => di%cached_face_value%normal(:,:,vele)
               
            else

               ! get the coordinate values for this element for each positions local node
               x_ele = ele_val(di%positions, vele)         

               ! get the coordinate values for this element for each quadrature point
               x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)

               ! The node indices of the positions projected to the pressure mesh
               x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)
            
               call transform_cvsurf_to_physical(x_ele, di%x_cvshape, detwei, normal, di%cvfaces)
            
            end if

            ! Initialise array for the quadrature points of this 
            ! element for whether it has already been visited
            notvisited = .true.

            ! Initialise the local rhs to assemble for this element
            div_tvel_rhs_local = 0.0

            ! loop over local nodes within this element
            nodal_loop_i: do iloc = 1, di%pressure_mesh%shape%loc

               ! loop over CV faces internal to this element
               face_loop: do face = 1, di%cvfaces%faces

                  ! is this a face neighbouring iloc?
                  is_neigh: if(di%cvfaces%neiloc(iloc, face) /= 0) then

                     ! find the opposing local node across the CV face
                     oloc = di%cvfaces%neiloc(iloc, face)

                     ! loop over gauss points on face
                     quadrature_loop: do gi = 1, di%cvfaces%shape%ngi

                        ! global gauss pt index for this element
                        ggi = (face-1)*di%cvfaces%shape%ngi + gi

                        ! check if this quadrature point has already been visited
                        check_visited: if(notvisited(ggi)) then

                           notvisited(ggi) = .false.

                           if (di%cached_face_value%cached_detwei_normal) then

                              ! cached normal has correct orientation already
                              normgi = normal(:,ggi)

                           else

                              ! correct the orientation of the normal so it points away from iloc
                              normgi = orientate_cvsurf_normgi(node_val(di%positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                              &x_face_quad(:,ggi), normal(:,ggi))

                           end if

                           ! Form the face value = detwei * (relperm*absperm/visc)
                           face_value = detwei(ggi) * &
                                        sum(di%cached_face_value%relperm(:,ggi,vele,p)) * &
                                        absperm_ele(1) / &
                                        visc_ele(1) 

                           darcy_vel_face_value_dot_n = - dot_product( normgi, &
                          &face_value * (grad_pressure_face_quad(:,ggi) - di%cached_face_value%den(ggi,vele,p) * grav_ele(:,1)) )

                           div_tvel_rhs_local(iloc) = div_tvel_rhs_local(iloc) + &
                                                      darcy_vel_face_value_dot_n 

                           div_tvel_rhs_local(oloc) = div_tvel_rhs_local(oloc) - &
                                                      darcy_vel_face_value_dot_n 

                        end if check_visited

                     end do quadrature_loop

                  end if is_neigh

               end do face_loop

            end do nodal_loop_i

            call addto(di%div_total_darcy_velocity, p_nodes, div_tvel_rhs_local)

         end do vele_loop

         ! Get this phase v BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
         call get_v_boundary_condition(di%darcy_velocity(p)%ptr, &
                                       (/"normal_flow   ", &
                                         "no_normal_flow"/), &
                                       di%bc_surface_mesh, &
                                       di%v_bc_value, &
                                       di%v_bc_flag)
         
         ! Get the pressure BC - required if no v given and weak pressure dirichlet given for extra integrals
         call get_entire_saturation_or_pressure_boundary_condition(di%pressure(p)%ptr, &
                                                                   (/"weakdirichlet"/), &
                                                                   di%bc_surface_mesh, &
                                                                   di%pressure_bc_value, &
                                                                   di%pressure_bc_flag)
         
         sele_loop: do sele = 1, di%number_sele
            
            if (di%v_bc_flag(sele) == V_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop

            p_nodes_bdy = face_global_nodes(di%pressure_mesh, sele)

            ! obtain the transformed determinant*weight and normals
            if (di%cached_face_value%cached_detwei_normal) then

               detwei_bdy => di%cached_face_value%detwei_bdy(:,sele)
               
               normal_bdy => di%cached_face_value%normal_bdy(:,:,sele)
            
            else

               x_ele     = ele_val(di%positions, face_ele(di%positions, sele))
               x_ele_bdy = face_val(di%positions, sele)
            
               call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)
            
            end if
                     
            if (di%v_bc_flag(sele) == V_BC_TYPE_NORMAL_FLOW) then
              
               bc_sele_val = ele_val(di%v_bc_value, sele)
            
            else
            
               visc_ele_bdy         = face_val(di%viscosity(p)%ptr, sele)
               absperm_ele_bdy      = face_val(di%absolute_permeability, sele)
               inv_char_len_ele_bdy = ele_val(di%inverse_characteristic_length, sele)
               if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then
                   bc_sele_val = ele_val(di%pressure_bc_value, sele)     
               end if
               grad_pressure_face_quad_bdy = face_val_at_quad(di%gradient_pressure(p)%ptr, sele, di%gradp_cvbdyshape)         
               grav_ele_bdy                = face_val(di%gravity, sele) 
               pressure_ele_bdy            = face_val(di%pressure(p)%ptr, sele) 

               ! the latest DarcyVelocity/relperm at the quadrature points
               do dim = 1,di%ndim

                  v_over_relperm_face_quad_bdy(dim,:) =  &
- (absperm_ele_bdy(1) / visc_ele_bdy(1)) * &
(grad_pressure_face_quad_bdy(dim,:) - face_val_at_quad(di%density(p)%ptr, sele, di%p_cvbdyshape) * grav_ele_bdy(dim,1))   

               end do

            end if

            ! Initialise the local rhs to assemble for this element
            div_tvel_rhs_local_bdy = 0.0

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi
                                                
                        if (di%v_bc_flag(sele) == V_BC_TYPE_NORMAL_FLOW) then

                           div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) + &
                                                          di%subcy_opt_sat%number * &
                                                          bc_sele_val(iloc) * &
                                                          detwei_bdy(ggi)                   
                        
                        else
                         
                           ! NOTE this v_over_relperm_dot_n does not contain the relperm terms
                           v_over_relperm_dot_n = dot_product(v_over_relperm_face_quad_bdy(:,ggi), normal_bdy(:,ggi))
                           
                           ! Adding this term includes the gradient of pressure, gradient of
                           ! capilliary pressure and the gravity term. 
                           div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) + &
                                                          sum(di%cached_face_value%relperm_bdy(:,ggi,sele,p)) * &
                                                          v_over_relperm_dot_n * &
                                                          detwei_bdy(ggi)
                           
                           ! Add two extra terms associated with weak pressure BC's
                                                      
                           ! Add two extra terms associated with weak pressure BC's
                           if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then

                              ! Form the face value = detwei * (relperm*absperm/visc)
                              face_value = detwei_bdy(ggi) * &
                                           sum(di%cached_face_value%relperm_bdy(:,ggi,sele,p)) * &
                                           absperm_ele_bdy(1) / &
                                           visc_ele_bdy(1) 

                              do jloc = 1, di%pressure_mesh%faces%shape%loc 

                                 div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) - &
                                                                inv_char_len_ele_bdy(jloc) * &
                                                                di%p_cvbdyshape%n(jloc,ggi) * &
                                                                bc_sele_val(jloc) * &
                                                                sum(normal_bdy(:,ggi)) * &
                                                                face_value * &
                                                                di%weak_pressure_bc_coeff

                              end do

                              do jloc = 1, di%pressure_mesh%faces%shape%loc 

                                 div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) + &
                                                                inv_char_len_ele_bdy(jloc) * &
                                                                di%p_cvbdyshape%n(jloc,ggi) * &
                                                                pressure_ele_bdy(jloc) * &
                                                                sum(normal_bdy(:,ggi)) * &
                                                                face_value * &
                                                                di%weak_pressure_bc_coeff

                              end do

                           end if 
                        
                        end if
                        
                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop

            call addto(di%div_total_darcy_velocity, p_nodes_bdy, div_tvel_rhs_local_bdy)

         end do sele_loop
         
      end do phase_loop
      
      call scale(di%div_total_darcy_velocity, di%inverse_cv_mass_pressure_mesh)
      
      ewrite_minmax(di%div_total_darcy_velocity)
      
      ! deallocate local variables as required
      deallocate(x_ele)
      if (di%cached_face_value%cached_detwei_normal) then
         nullify(detwei)
         nullify(normal)      
      else
         deallocate(detwei)
         deallocate(normal)
      end if
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(div_tvel_rhs_local)      
      deallocate(x_face_quad)
      deallocate(grad_pressure_face_quad)
      deallocate(grav_ele)

      deallocate(grav_ele_bdy)
      deallocate(grad_pressure_face_quad_bdy)
      deallocate(v_over_relperm_face_quad_bdy)
      deallocate(bc_sele_val)
      deallocate(pressure_ele_bdy)
      deallocate(inv_char_len_ele_bdy)
      if (di%cached_face_value%cached_detwei_normal) then
         nullify(detwei_bdy)
         nullify(normal_bdy)      
      else
         deallocate(detwei_bdy)
         deallocate(normal_bdy)
      end if
      deallocate(div_tvel_rhs_local_bdy)
      deallocate(x_ele_bdy)
      deallocate(p_nodes_bdy)

      ewrite(1,*) 'Finished Calculate the DivergenceTotalDarcyVelocity'
            
   end subroutine darcy_impes_calculate_divergence_total_darcy_velocity

! ----------------------------------------------------------------------------
                  
   subroutine assemble_rhs_adv(di, &
                               p, &
                               isub, &
                               form_new_subcycle_relperm_face_values)

      !!< Assemble the rhs advection contribtion for saturation phase p
      
      type(darcy_impes_type), intent(inout) :: di
      integer,                intent(in)    :: p      
      integer,                intent(in)    :: isub    
      logical,                intent(in)    :: form_new_subcycle_relperm_face_values    
      logical,                intent(in)    :: isub    

      ! local variables
      logical :: inflow, determine_face_value
      integer :: vele, iloc, oloc, jloc, face, gi, ggi, sele, upwind_pos, dim
      real    :: income, face_value, v_over_relperm_dot_n, old_relperm_face_value
      real,    dimension(1)            :: absperm_ele, visc_ele
      real,    dimension(:,:), pointer :: grad_pressure_face_quad
      real,    dimension(:,:), pointer :: grav_ele
      real,    dimension(:,:), pointer :: v_over_relperm_face_quad
      real,    dimension(:),   pointer :: old_relperm_subcycle_ele
      real,    dimension(:,:), pointer :: x_ele
      real,    dimension(:,:), pointer :: normal
      real,    dimension(:),   pointer :: detwei
      real,    dimension(:),   pointer :: normgi
      logical, dimension(:),   pointer :: notvisited
      real,    dimension(:),   pointer :: s_rhs_local
      real,    dimension(:,:), pointer :: x_face_quad
      integer, dimension(:),   pointer :: x_pmesh_nodes
      integer, dimension(:),   pointer :: p_nodes      
      integer, dimension(:),   pointer :: upwind_nodes
      real,    dimension(1)            :: absperm_ele_bdy, visc_ele_bdy
      real,    dimension(:,:), pointer :: grav_ele_bdy
      real,    dimension(:,:), pointer :: grad_pressure_face_quad_bdy
      real,    dimension(:,:), pointer :: v_over_relperm_face_quad_bdy
      real,    dimension(:),   pointer :: bc_sele_val
      real,    dimension(:),   pointer :: pressure_ele_bdy
      real,    dimension(:),   pointer :: inv_char_len_ele_bdy
      real,    dimension(:),   pointer :: old_relperm_subcycle_ele_bdy
      real,    dimension(:,:), pointer :: normal_bdy
      real,    dimension(:),   pointer :: detwei_bdy
      real,    dimension(:,:), pointer :: x_ele_bdy
      real,    dimension(:),   pointer :: s_rhs_local_bdy
      integer, dimension(:),   pointer :: p_nodes_bdy
      integer, parameter :: SATURATION_BC_TYPE_WEAKDIRICHLET = 1, SATURATION_BC_TYPE_DIRICHLET    = 2
      integer, parameter :: V_OVER_S_BC_TYPE_NORMAL_FLOW     = 1, V_OVER_S_BC_TYPE_NO_NORMAL_FLOW = 2
      integer, parameter :: PRESSURE_BC_TYPE_WEAKDIRICHLET   = 1
      
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(normal(di%ndim,di%x_cvshape%ngi))
         allocate(detwei(di%x_cvshape%ngi))      
      end if
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(s_rhs_local(ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(old_relperm_subcycle_ele(ele_loc(di%pressure_mesh,1)))
      allocate(grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(grav_ele(di%ndim,1))
      allocate(v_over_relperm_face_quad(di%ndim,di%p_cvshape%ngi))

      allocate(grav_ele_bdy(di%ndim,1))
      allocate(grad_pressure_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      allocate(v_over_relperm_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      allocate(bc_sele_val(face_loc(di%pressure_mesh,1)))
      allocate(pressure_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(inv_char_len_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(old_relperm_subcycle_ele_bdy(face_loc(di%pressure_mesh,1)))      
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(detwei_bdy(di%x_cvbdyshape%ngi))
         allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      end if
      allocate(s_rhs_local_bdy(face_loc(di%pressure_mesh,1)))
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(s_nodes_bdy(face_loc(di%pressure_mesh,1)))      
      
      ! Inititalise rhs advection field
      call zero(di%rhs_adv)

      ! Decide if a face value needs to be determined for this subcycle
      determine_face_value = (form_new_subcycle_relperm_face_values .and. isub > 1)

      if (determine_face_value) then

         ! Determine the upwind relperm values if required for higher order CV face value
         if(need_upwind_values(di%relperm_cv_options)) then

           call find_upwind_values(di%state, &
                                   di%positions_pressure_mesh, &
                                   di%old_relperm_subcycle(p)%ptr, &
                                   di%old_relperm_upwind, &
                                   di%old_relperm_subcycle(p)%ptr, &
                                   di%old_relperm_upwind, &
                                   option_path = trim(di%relative_permeability(p)%ptr%option_path))

         else

           call zero(di%old_relperm_upwind)

         end if

      end if

      ! Initialise optimisation flag used in finding upwind value in high resolution schemes
      upwind_pos = 0

      ! Loop volume elements assembling local contributions    
      vol_element_loop: do vele = 1, di%number_vele

         ! The node indices of the pressure mesh (which saturation is associated with)
         p_nodes => ele_nodes(di%pressure_mesh, vele)

         ! The node indices of the positions projected to the pressure mesh (which saturation is associated with)
         x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)

         ! The gravity values for this element for each direction
         grav_ele = ele_val(di%gravity, vele) 

         ! get the viscosity value for this element
         visc_ele = ele_val(di%viscosity(p)%ptr, vele)

         ! get the absolute permeability value for this element
         absperm_ele = ele_val(di%absolute_permeability, vele)         

         ! get the latest gradient pressure at the cv surface quadrature points for each direction 
         grad_pressure_face_quad = ele_val_at_quad(di%gradient_pressure(p)%ptr, vele, di%gradp_cvshape)

         ! obtain the transformed determinant*weight and normals
         if (di%cached_face_value%cached_detwei_normal) then

            detwei => di%cached_face_value%detwei(:,vele)

            normal => di%cached_face_value%normal(:,:,vele)

         else

            ! get the coordinate values for this element for each positions local node
            x_ele = ele_val(di%positions, vele)         

            ! get the coordinate values for this element for each quadrature point
            x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)

            call transform_cvsurf_to_physical(x_ele, di%x_cvshape, detwei, normal, di%cvfaces)

         end if

         ! Get necessary element data for determining relperm CV face value
         if (determine_face_value) then

            ! get the old relperm ele values from start of subcycle
            old_relperm_subcycle_ele = ele_val(di%old_relperm_subcycle(p)%ptr, vele)

            ! Determine the node numbers to use to determine the upwind values
            if((di%relperm_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
               (di%relperm_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

               upwind_nodes => x_pmesh_nodes

            else

               upwind_nodes => p_nodes

            end if

            ! the latest DarcyVelocity/relperm at the quadrature points for start of subcycle
            ! determined from FE interpolation of each component, only used to determine upwind.
            do dim = 1, di%ndim

               v_over_relperm_face_quad(dim,:) = - (absperm_ele(1) / visc_ele(1)) * &
(grad_pressure_face_quad(dim,:) - ele_val_at_quad(di%density(p)%ptr, vele, di%p_cvshape) * grav_ele(dim,1))   

            end do

         end if

         ! Initialise array for the quadrature points of this 
         ! element for whether it has already been visited
         notvisited = .true.

         ! Initialise the local rhs to assemble for this element
         s_rhs_local = 0.0

         ! loop over local nodes within this element
         nodal_loop_i: do iloc = 1, di%pressure_mesh%shape%loc

            ! loop over CV faces internal to this element
            face_loop: do face = 1, di%cvfaces%faces

               ! is this a face neighbouring iloc?
               is_neigh: if(di%cvfaces%neiloc(iloc, face) /= 0) then

                  ! find the opposing local node across the CV face
                  oloc = di%cvfaces%neiloc(iloc, face)

                  ! loop over gauss points on face
                  quadrature_loop: do gi = 1, di%cvfaces%shape%ngi

                     ! global gauss pt index for this element
                     ggi = (face-1)*di%cvfaces%shape%ngi + gi

                     ! check if this quadrature point has already been visited
                     check_visited: if(notvisited(ggi)) then

                        notvisited(ggi) = .false.

                        if (di%cached_face_value%cached_detwei_normal) then

                           ! cached normal has correct orientation already
                           normgi = normal(:,ggi)

                        else

                           ! correct the orientation of the normal so it points away from iloc
                           normgi = orientate_cvsurf_normgi(node_val(di%positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                           &x_face_quad(:,ggi), normal(:,ggi))

                        end if

                        ! Determine the relperm CV face value
                        if (determine_face_value) then

                           ! determine if the flow is in or out of the face at this quadrature
                           ! with respect to the normal orientation using latest v_over_relperm
                           v_over_relperm_dot_n = dot_product(v_over_relperm_face_quad(:,ggi), normgi(:))

                           inflow = (v_over_relperm_dot_n<=0.0)

                           income = merge(1.0,0.0,inflow)

                           ! evaluate the nonlinear face value for relperm
                           call darcy_impes_evaluate_face_val(old_relperm_face_value, &
                                                              iloc, &
                                                              oloc, &
                                                              ggi, &
                                                              upwind_nodes, &
                                                              di%p_cvshape, &
                                                              old_relperm_subcycle_ele, &
                                                              di%old_relperm_upwind, &
                                                              inflow, &
                                                              income, &
                                                              di%relperm_cv_options, &
                                                              save_pos = upwind_pos)

                           ! Cache the relperm face value if required
                           if (cache_new_subcycle_relperm_face_values) then

                              di%cached_face_value%relperm(isub,ggi,vele,p) = old_relperm_face_value

                           end if

                        else 

                           old_relperm_face_value = di%cached_face_value%relperm(isub,ggi,vele,p)

                        end if

                        face_value = detwei(ggi) * old_relperm_face_value * absperm_ele(1) * &
dot_product((grad_pressure_face_quad(:,ggi) - di%cached_face_value%den(ggi,vele,p) * grav_ele(:,1)), normgi)/ visc_ele(1)

                        ! Form the local rhs for iloc and opposing oloc with normal vector sign change
                        s_rhs_local(iloc) = s_rhs_local(iloc) + face_value

                        s_rhs_local(oloc) = s_rhs_local(oloc) - face_value

                     end if check_visited

                  end do quadrature_loop

               end if is_neigh

            end do face_loop

         end do nodal_loop_i

         ! Add volume element contribution to global rhs advection field
         call addto(di%rhs_adv, p_nodes, s_rhs_local)

      end do vol_element_loop

      ! Add BC integrals

      sele_loop: do sele = 1, di%number_sele

         if (di%v_bc_flag(sele) == V_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop

         p_nodes_bdy = face_global_nodes(di%pressure_mesh, sele)

         ! obtain the transformed determinant*weight and normals
         if (di%cached_face_value%cached_detwei_normal) then

            detwei_bdy => di%cached_face_value%detwei_bdy(:,sele)

            normal_bdy => di%cached_face_value%normal_bdy(:,:,sele)

         else

            x_ele     = ele_val(di%positions, face_ele(di%positions, sele))
            x_ele_bdy = face_val(di%positions, sele)

            call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)

         end if

         if (di%v_bc_flag(sele) == V_BC_TYPE_NORMAL_FLOW) then

            bc_sele_val = ele_val(di%v_bc_value, sele)

         else

            visc_ele_bdy         = face_val(di%viscosity(p)%ptr, sele)
            absperm_ele_bdy      = face_val(di%absolute_permeability, sele)
            inv_char_len_ele_bdy = ele_val(di%inverse_characteristic_length, sele)
            if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then
                bc_sele_val = ele_val(di%pressure_bc_value, sele)     
            end if
            grad_pressure_face_quad_bdy = face_val_at_quad(di%gradient_pressure(p)%ptr, sele, di%gradp_cvbdyshape)         
            grav_ele_bdy                = face_val(di%gravity, sele) 
            pressure_ele_bdy            = face_val(di%pressure(p)%ptr, sele) 

            old_relperm_subcycle_ele_bdy = face_val(di%old_relperm_subcycle(p)%ptr, sele)

            ! the latest DarcyVelocity/relperm at the quadrature points for start of subcycle
            ! determined from FE interpolation of each component, 
            do dim = 1, di%ndim

               v_over_relperm_face_quad_bdy(dim,:) =  &
- (absperm_ele_bdy(1) / visc_ele_bdy(1)) * &
(grad_pressure_face_quad_bdy(dim,:) - face_val_at_quad(di%density(p)%ptr, sele, di%p_cvbdyshape) * grav_ele_bdy(dim,1))   

            end do

         end if               

         ! Initialise the local rhs to assemble for this element
         s_rhs_local_bdy = 0.0

         bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

            bc_face_loop: do face = 1, di%cvfaces%sfaces

               bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                  bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                     ggi = (face-1)*di%cvfaces%shape%ngi + gi

                     if (di%v_bc_flag(sele) == V_BC_TYPE_NORMAL_FLOW) then

                        s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - &
                                                bc_sele_val(iloc) * &
                                                detwei_bdy(ggi)

                     else                            

                        ! NOTE this v_over_relperm_dot_n does not contain the relperm terms
                        v_over_relperm_dot_n = dot_product(v_over_relperm_face_quad_bdy(:,ggi), normal_bdy(:,ggi))

                        if (determine_face_value) then

                           old_relperm_face_value = old_relperm_subcycle_ele_bdy(iloc)

                           ! Cache the saturation face value if required
                           if (cache_new_subcycle_sat_face_values) then

                              di%cached_face_value%relperm_bdy(isub,ggi,sele,p) = old_relperm_face_value

                           end if 

                        else 

                           old_relperm_face_value = di%cached_face_value%relperm_bdy(isub,ggi,sele,p)                              

                        end if

                        ! This term includes the gradient pressure, gradient capilliary pressure 
                        ! and the gravity term. 
                        s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - &
                                                old_relperm_face_value * &
                                                v_over_relperm_dot_n * &
                                                detwei_bdy(ggi)

                        ! Add two extra terms associated with weak pressure BC's                              
                        if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then

                           ! Form the face value = detwei * (relperm*absperm/visc)
                           face_value = detwei_bdy(ggi) * &
                                        old_relperm_face_value * &
                                        absperm_ele_bdy(1) / &
                                        visc_ele_bdy(1) 

                           if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then

                              do jloc = 1, di%pressure_mesh%faces%shape%loc 

                                 s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - &
                                                         inv_char_len_ele_bdy(jloc) * &
                                                         di%p_cvbdyshape%n(jloc,ggi) * &
                                                         bc_sele_val(jloc) * &
                                                         sum(normal_bdy(:,ggi)) * &
                                                         face_value * &
                                                         di%weak_pressure_bc_coeff

                              end do

                           end if 

                           do jloc = 1, di%pressure_mesh%faces%shape%loc 

                              s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) + &
                                                      inv_char_len_ele_bdy(jloc) * &
                                                      di%p_cvbdyshape%n(jloc,ggi) * &
                                                      pressure_ele_bdy(jloc) * &
                                                      sum(normal_bdy(:,ggi)) * &
                                                      face_value * &
                                                      di%weak_pressure_bc_coeff

                           end do

                        end if

                     end if

                  end do bc_quad_loop

               end if bc_neigh_if

            end do bc_face_loop

         end do bc_iloc_loop

         call addto(di%rhs_adv, p_nodes_bdy, s_rhs_local_bdy)

      end do sele_loop

      ! deallocate local variables as required
      deallocate(x_ele)
      if (di%cached_face_value%cached_detwei_normal) then
         nullify(detwei)
         nullify(normal)      
      else
         deallocate(detwei)
         deallocate(normal)
      end if
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(s_rhs_local)
      deallocate(x_face_quad)
      deallocate(old_relperm_subcycle_ele)
      deallocate(grad_pressure_face_quad)
      deallocate(grav_ele)
      deallocate(v_over_relperm_face_quad)

      deallocate(grav_ele_bdy)
      deallocate(grad_pressure_face_quad_bdy)
      deallocate(v_over_relperm_face_quad_bdy)
      deallocate(bc_sele_val)
      deallocate(pressure_ele_bdy)
      deallocate(inv_char_len_ele_bdy)
      deallocate(old_relperm_subcycle_ele_bdy)      
      if (di%cached_face_value%cached_detwei_normal) then
         nullify(detwei_bdy)
         nullify(normal_bdy)      
      else
         deallocate(detwei_bdy)
         deallocate(normal_bdy)
      end if
      deallocate(s_rhs_local_bdy)
      deallocate(x_ele_bdy)
      deallocate(s_nodes_bdy)

   end subroutine assemble_rhs_adv
   
! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_phase_one_saturation_diagnostic(di)
      
      !!< Calculate the first phase diagnostic saturation as 1.0 - the rest
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p
      
      ewrite(1,*) 'Calculate phase 1 diagnostic Saturation'
      
      call set(di%saturation(1)%ptr, 1.0)
      
      do p = 2, di%number_phase
         
         call addto(di%saturation(1)%ptr, di%saturation(p)%ptr, scale = -1.0)
         
      end do
      
      ewrite_minmax(di%saturation(1)%ptr)
       
   end subroutine darcy_impes_calculate_phase_one_saturation_diagnostic

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_vel_mob_ff_and_cfl_fields(di)
      
      !!< Calculate the various velocities, mobilities, fractional flows, and CFL fields
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      logical :: inflow
      integer :: vele, p, dim, iloc, oloc, face, gi, ggi, sele
      real    :: income, v_dot_n_face_value
      real,    dimension(1)              :: visc_ele, absperm_ele
      real,    dimension(:,:),   pointer :: grad_pressure_face_quad
      real,    dimension(:,:),   pointer :: grad_pressure_ele
      real,    dimension(:),     pointer :: v_local
      real,    dimension(:),     pointer :: m_local
      real,    dimension(:,:),   pointer :: grav_ele
      real,    dimension(:),     pointer :: den_ele
      real,    dimension(:),     pointer :: relperm_ele      
      real,    dimension(:,:),   pointer :: x_ele
      real,    dimension(:,:),   pointer :: normal
      real,    dimension(:),     pointer :: detwei
      real,    dimension(:),     pointer :: normgi
      logical, dimension(:),     pointer :: notvisited
      real,    dimension(:),     pointer :: cfl_rhs_local
      real,    dimension(:,:),   pointer :: x_face_quad
      integer, dimension(:),     pointer :: x_pmesh_nodes
      integer, dimension(:),     pointer :: p_nodes      
      real,    dimension(1)              :: visc_ele_bdy, absperm_ele_bdy
      real,    dimension(:,:),   pointer :: grav_ele_bdy
      real,    dimension(:,:),   pointer :: grad_pressure_face_quad_bdy
      real,    dimension(:,:),   pointer :: normal_bdy
      real,    dimension(:),     pointer :: detwei_bdy
      real,    dimension(:,:),   pointer :: x_ele_bdy
      real,    dimension(:),     pointer :: cfl_rhs_local_bdy
      integer, dimension(:),     pointer :: p_nodes_bdy
      real,    dimension(:),     pointer :: bc_sele_val
      integer, parameter :: V_BC_TYPE_NORMAL_FLOW = 1, V_BC_TYPE_NO_NORMAL_FLOW = 2
            
      ewrite(1,*) 'Calculate CFL, Velocity, Mobilities and FractionalFlow fields'

      ! allocate arrays used in CFL assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(normal(di%ndim,di%x_cvshape%ngi))
         allocate(detwei(di%x_cvshape%ngi))      
      end if
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(cfl_rhs_local(ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(grav_ele(di%ndim,1))
      
      allocate(grav_ele_bdy(di%ndim,1))
      allocate(grad_pressure_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(detwei_bdy(di%x_cvbdyshape%ngi))
         allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      end if
      allocate(cfl_rhs_local_bdy(face_loc(di%pressure_mesh,1)))
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(p_nodes_bdy(face_loc(di%pressure_mesh,1)))
      allocate(bc_sele_val(face_loc(di%pressure_mesh,1)))
      
      allocate(den_ele(ele_loc(di%pressure_mesh,1)))
      allocate(relperm_ele(ele_loc(di%pressure_mesh,1)))
      allocate(grad_pressure_ele(di%ndim,1))
      allocate(v_local(ele_loc(di%velocity_mesh,1)))
      allocate(m_local(ele_loc(di%velocity_mesh,1)))

      phase_loop: do p = 1, di%number_phase
         
         call zero(di%cfl(p)%ptr)
              
         vele_loop: do vele = 1, di%number_vele

            ! The node indices of the pressure field
            p_nodes => ele_nodes(di%pressure_mesh, vele)
            
            ! The gravity values for this element for each direction
            grav_ele = ele_val(di%gravity, vele) 

            ! get the viscosity value for this element for this phase
            visc_ele = ele_val(di%viscosity(p)%ptr, vele)

            ! get the absolute permeability value for this element
            absperm_ele = ele_val(di%absolute_permeability, vele)         

            ! get the latest gradient pressure at the cv surface quadrature points for each direction for this phase
            grad_pressure_face_quad = ele_val_at_quad(di%gradient_pressure(p)%ptr, vele, di%gradp_cvshape)

            ! obtain the transformed determinant*weight and normals
            if (di%cached_face_value%cached_detwei_normal) then
               
               detwei => di%cached_face_value%detwei(:,vele)
               
               normal => di%cached_face_value%normal(:,:,vele)
               
            else

               ! get the coordinate values for this element for each positions local node
               x_ele = ele_val(di%positions, vele)         

               ! get the coordinate values for this element for each quadrature point
               x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)

               ! The node indices of the positions projected to the pressure mesh
               x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)
            
               call transform_cvsurf_to_physical(x_ele, di%x_cvshape, detwei, normal, di%cvfaces)
            
            end if

            ! Initialise array for the quadrature points of this 
            ! element for whether it has already been visited
            notvisited = .true.

            ! Initialise the local rhs's to assemble for this element
            cfl_rhs_local = 0.0

            ! loop over local nodes within this element
            nodal_loop_i: do iloc = 1, di%pressure_mesh%shape%loc

               ! loop over CV faces internal to this element
               face_loop: do face = 1, di%cvfaces%faces

                  ! is this a face neighbouring iloc?
                  is_neigh: if(di%cvfaces%neiloc(iloc, face) /= 0) then

                     ! find the opposing local node across the CV face
                     oloc = di%cvfaces%neiloc(iloc, face)

                     ! loop over gauss points on face
                     quadrature_loop: do gi = 1, di%cvfaces%shape%ngi

                        ! global gauss pt index for this element
                        ggi = (face-1)*di%cvfaces%shape%ngi + gi

                        ! check if this quadrature point has already been visited
                        check_visited: if(notvisited(ggi)) then

                           notvisited(ggi) = .false.

                           if (di%cached_face_value%cached_detwei_normal) then

                              ! cached normal has correct orientation already
                              normgi = normal(:,ggi)

                           else

                              ! correct the orientation of the normal so it points away from iloc
                              normgi = orientate_cvsurf_normgi(node_val(di%positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                              &x_face_quad(:,ggi), normal(:,ggi))

                           end if

                           v_dot_n_face_value = &
di%cached_face_value%relperm(1,ggi,vele,p) * absperm_ele(1) * &
dot_product((grad_pressure_face_quad(:,ggi) - di%cached_face_value%den(ggi,vele,p) * grav_ele(:,1)), normgi)/ visc_ele(1)

                           inflow = (v_dot_n_face_value<=0.0)

                           income = merge(1.0,0.0,inflow)

                           cfl_rhs_local(iloc) = cfl_rhs_local(iloc) + &
                                                 abs(v_dot_n_face_value) * &
                                                 detwei(ggi) * &
                                                 (1.0 - income)

                           cfl_rhs_local(oloc) = cfl_rhs_local(oloc) + &
                                                 abs(v_dot_n_face_value) * &
                                                 detwei(ggi) * &
                                                 income

                        end if check_visited

                     end do quadrature_loop

                  end if is_neigh

               end do face_loop

            end do nodal_loop_i

            call addto(di%cfl(p)%ptr, p_nodes, cfl_rhs_local)

         end do vele_loop

         ! Get this phase v BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
         call get_v_boundary_condition(di%darcy_velocity(p)%ptr, &
                                       (/"normal_flow   ", &
                                         "no_normal_flow"/), &
                                       di%bc_surface_mesh, &
                                       di%v_bc_value, &
                                       di%v_bc_flag)
         
         sele_loop: do sele = 1, di%number_sele
            
            if (di%v_bc_flag(sele) == V_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop

            p_nodes_bdy = face_global_nodes(di%pressure_mesh, sele)

            ! obtain the transformed determinant*weight and normals
            if (di%cached_face_value%cached_detwei_normal) then

               detwei_bdy => di%cached_face_value%detwei_bdy(:,sele)
               
               normal_bdy => di%cached_face_value%normal_bdy(:,:,sele)
            
            else

               x_ele     = ele_val(di%positions, face_ele(di%positions, sele))
               x_ele_bdy = face_val(di%positions, sele)
            
               call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)
            
            end if
            
            if (di%v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NORMAL_FLOW) then
            
               bc_sele_val = ele_val(di%v_over_s_bc_value, sele)
            
            else
            
               visc_ele_bdy                = face_val(di%viscosity(p)%ptr, sele)
               absperm_ele_bdy             = face_val(di%absolute_permeability, sele)
               grad_pressure_face_quad_bdy = face_val_at_quad(di%gradient_pressure(p)%ptr, sele, di%gradp_cvbdyshape)         
               grav_ele_bdy                = face_val(di%gravity, sele) 
            
            end if

            ! Initialise the local rhs to assemble for this element
            cfl_rhs_local_bdy = 0.0

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi
                        
                        if (di%v_bc_flag(sele) == V_BC_TYPE_NORMAL_FLOW) then
                        
                           v_dot_n_face_value = bc_sele_val(iloc)
                                                      
                        else

                           v_dot_n_face_value = &
di%cached_face_value%relperm_bdy(1,ggi,sele,p) * absperm_ele_bdy(1) * &
dot_product((grad_pressure_face_quad_bdy(:,ggi) - di%cached_face_value%den_bdy(ggi,sele,p) * grav_ele_bdy(:,1)), normal_bdy)/ visc_ele_bdy(1)
                                                                              
                        end if
                        
                        if (v_dot_n_face_value > 0.0) then
                           income = 0.0
                        else
                           income = 1.0
                        end if
                        
                        cfl_rhs_local_bdy(iloc) = cfl_rhs_local_bdy(iloc) + &
                                                  abs(v_dot_n_face_value) * &
                                                  detwei_bdy(ggi) * &
                                                  (1.0 - income)                     

                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop

            call addto(di%cfl(p)%ptr, p_nodes_bdy, cfl_rhs_local_bdy)

         end do sele_loop
                  
         di%cfl(p)%ptr%val = di%cfl(p)%ptr%val * di%dt / di%cv_mass_pressure_mesh_with_porosity%val

         ewrite_minmax(di%cfl(p)%ptr)

      end do phase_loop      

      
      ewrite(1,*) 'Calculate DarcyVelocity and Mobility'
      
      ! Calculate the darcy_velocity field = - (absperm * relperm / visc) * ( grad_pressure - den * grav), 
      
      ! Calculate mobility = relperm / visc
      
      do p = 1, di%number_phase
                  
         do vele = 1, element_count(di%velocity_mesh)

            ! Find the element wise absolute_permeability
            absperm_ele = ele_val(di%absolute_permeability, vele)

            ! Find the element wise viscosity
            visc_ele = ele_val(di%viscosity(p)%ptr, vele)

            ! Find the element wise local values for relperm
            relperm_ele = ele_val(di%relative_permeability(p)%ptr, vele)

            ! Find the element wise local values for density
            den_ele = ele_val(di%density(p)%ptr, vele)
         
            ! The gravity values for this element for each direction
            grav_ele = ele_val(di%gravity, vele) 

            ! Find the element wise gradient pressure 
            grad_pressure_ele = ele_val(di%gradient_pressure(p)%ptr, vele)
                        
            do iloc = 1, di%velocity_mesh%shape%loc
               
               m_local(iloc) = relperm_ele(iloc) / visc_ele(1)
                              
            end do 
            
            call set(di%mobility(p)%ptr, &
                     ele_nodes(di%velocity_mesh, vele), &
                     m_local)
            
            do dim = 1,di%ndim

               do iloc = 1, di%velocity_mesh%shape%loc
               
                  v_local(iloc) = - (relperm_ele(iloc) * absperm_ele(1) / visc_ele(1)) * &
                                  & ( grad_pressure_ele(dim,1) - den_ele(iloc) * grav_ele(dim,1))
                                 
               end do

               call set(di%darcy_velocity(p)%ptr, &
                        dim, &
                        ele_nodes(di%velocity_mesh, vele), &
                        v_over_s_local)               
               
            end do

         end do
         
      end do 
                  
      do p = 1,di%number_phase
         ewrite_minmax(di%darcy_velocity(p)%ptr)
      end do
            
      ewrite(1,*) 'Calculate TotalDarcyVelocity and TotalMobility'
      
      call set(di%total_darcy_velocity, di%darcy_velocity(1)%ptr)
      call set(di%total_mobility, di%mobility(1)%ptr)
      
      do p = 2, di%number_phase
         
         call addto(di%total_darcy_velocity, di%darcy_velocity(p)%ptr)
         call addto(di%total_mobility, di%mobility(p)%ptr)
         
      end do

      ewrite_minmax(di%total_darcy_velocity)
            
      ! calculate the fractional flow for each phase
      do p = 1, di%number_phase
         
         ewrite(1,*) 'Calculate FractionalFlow for phase ',p
                  
         call set(di%fractional_flow(p)%ptr, di%total_mobility)
         
         call invert(di%fractional_flow(p)%ptr)
         
         call scale(di%fractional_flow(p)%ptr, di%mobility(p)%ptr)
                  
         ewrite_minmax(di%fractional_flow(p)%ptr)
      
      end do      
      
      ! deallocate local variables as required
      deallocate(x_ele)
      if (di%cached_face_value%cached_detwei_normal) then
         nullify(detwei)
         nullify(normal)      
      else
         deallocate(detwei)
         deallocate(normal)
      end if
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(cfl_rhs_local)
      deallocate(x_face_quad)
      deallocate(grad_pressure_face_quad)
      deallocate(grav_ele)
      
      deallocate(grav_ele_bdy)
      deallocate(grad_pressure_face_quad_bdy)
      if (di%cached_face_value%cached_detwei_normal) then
         nullify(detwei_bdy)
         nullify(normal_bdy)      
      else
         deallocate(detwei_bdy)
         deallocate(normal_bdy)
      end if
      deallocate(x_ele_bdy)
      deallocate(p_nodes_bdy)
      deallocate(cfl_rhs_local_bdy)
      deallocate(bc_sele_val)      

      deallocate(den_ele)
      deallocate(relperm_ele)
      deallocate(grad_pressure_ele)
      deallocate(v_local)
      deallocate(m_local)
             
   end subroutine darcy_impes_calculate_vel_mob_ff_and_cfl_fields

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_sum_saturation(di)
      
      !!< Calculate the sum of the saturation fields
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p
      
      ewrite(1,*) 'Calculate SumSaturation'
      
      call zero(di%sum_saturation)
      
      do p = 1, di%number_phase
         
         call addto(di%sum_saturation, di%saturation(p)%ptr)
         
      end do
      
      ewrite_minmax(di%sum_saturation)
       
   end subroutine darcy_impes_calculate_sum_saturation

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_densities(di)
      
      !!< Calculate the density field of each phase
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p
      
      ewrite(1,*) 'Calculate Density of each phase'
            
      do p = 1, di%number_phase
         
         if (di%eos_options(p)%have_fluids_linear) then
            
            ! Currently no variation of the reference density by any pertubation
            
            call set(di%density(p)%ptr, di%eos_options(p)%fluids_linear_reference_density)
         
         end if

         ewrite_minmax(di%density(p)%ptr)
         
      end do 

      ewrite(1,*) 'Finished Calculate Density of each phase'
       
   end subroutine darcy_impes_calculate_densities

! ----------------------------------------------------------------------------
    
   subroutine darcy_impes_calculate_cflnumber_field_based_dt(di)
      
      !!< Calculate the cfl number field based adaptive time step
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p
      real, dimension(di%number_phase) :: phase_dt
      
      ewrite(1,*) 'Calculate CFL number field based adaptive timestep'
      
      ! Find the adaptive dt for each phase then take the minimum
      
      phase_dt = di%dt
      
      phase_loop: do p = 1, di%number_phase
         
         phase_dt(p) = cflnumber_field_based_dt(di%cfl(p)%ptr, &
                                                di%dt, &
                                                di%adaptive_dt_options%requested_cfl, &
                                                di%adaptive_dt_options%min_dt, &
                                                di%adaptive_dt_options%max_dt, &
                                                di%adaptive_dt_options%increase_tolerance)         
         
      end do phase_loop
      
      di%dt = minval(phase_dt)

      if(di%adaptive_dt_options%min_dt_terminate_if_reached) then
         if(di%dt <= di%adaptive_dt_options%min_dt + spacing(di%adaptive_dt_options%min_dt)) then
            ewrite(0, *) "Minimum timestep reached - terminating"
            SIG_INT = .true.
         end if
      end if

      ewrite(1,*) 'Finished Calculate CFL number field based adaptive timestep'
       
   end subroutine darcy_impes_calculate_cflnumber_field_based_dt

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_initialise_cached_face_value(di)
   
      !!< Initialise the cached_face_value variable's by 
      !!< allocating and zeroing the values and setting the 
      !!< detwei, normal and p_dshape components as required.

      type(darcy_impes_type), intent(inout) :: di
      
      ! local variable
      integer :: vele, sele, iloc, face, gi, ggi
      real,    dimension(:),   allocatable :: normgi
      logical, dimension(:),   allocatable :: notvisited
      real,    dimension(:,:), allocatable :: x_face_quad
      integer, dimension(:),   pointer     :: x_pmesh_nodes

      ewrite(1,*) 'Initialise cached phase face value types'
            
      allocate(di%cached_face_value%relperm(di%p_cvshape%ngi, di%number_vele, di%number_phase))
      allocate(di%cached_face_value%relperm_bdy(di%p_cvbdyshape%ngi, di%number_sele, di%number_phase))
      
      di%cached_face_value%relperm     = 0.0
      di%cached_face_value%relperm_bdy = 0.0
      
      allocate(di%cached_face_value%den(di%p_cvshape%ngi, di%number_vele, di%number_phase))
      allocate(di%cached_face_value%den_bdy(di%p_cvbdyshape%ngi, di%number_sele, di%number_phase))
      
      di%cached_face_value%den     = 0.0
      di%cached_face_value%den_bdy = 0.0
      
      if (di%cached_face_value%cached_detwei_normal) then
      
         allocate(di%cached_face_value%detwei(di%p_cvshape%ngi, di%number_vele))
         allocate(di%cached_face_value%detwei_bdy(di%p_cvbdyshape%ngi, di%number_sele))
     
         di%cached_face_value%detwei     = 0.0
         di%cached_face_value%detwei_bdy = 0.0
      
         allocate(di%cached_face_value%normal(di%ndim, di%p_cvshape%ngi, di%number_vele))
         allocate(di%cached_face_value%normal_bdy(di%ndim, di%p_cvbdyshape%ngi, di%number_sele))

         di%cached_face_value%normal     = 0.0
         di%cached_face_value%normal_bdy = 0.0
      
      end if
      
      if (di%cached_face_value%cached_p_dshape) then
      
         allocate(di%cached_face_value%p_dshape(di%pressure_mesh%shape%loc, di%p_cvshape%ngi, di%ndim, di%number_vele))
!!!!! *** THIS IS NOT POSSIBLE YET ***
!!!!!         allocate(di%cached_face_value%p_dshape_bdy(di%pressure_mesh%faces%shape%loc, di%p_cvbdyshape_full%ngi, di%ndim, di%number_sele))

         di%cached_face_value%p_dshape     = 0.0
!!!!! *** THIS IS NOT POSSIBLE YET ***
!!!!!         di%cached_face_value%p_dshape_bdy = 0.0
      
      end if

      allocate(normgi(di%ndim))
      allocate(notvisited(di%p_cvshape%ngi))
      allocate(x_face_quad(di%ndim, di%p_cvshape%ngi))
            
      ! Determine the detwei and normal if required
      find_cached_detwei_normal: if (di%cached_face_value%cached_detwei_normal) then
      
         vele_loop_detwei_normal: do vele = 1, di%number_vele

            call transform_cvsurf_to_physical(ele_val(di%positions, vele), &
                                              di%x_cvshape, &
                                              di%cached_face_value%detwei(:,vele), &
                                              di%cached_face_value%normal(:,:,vele), &
                                              di%cvfaces)
            
            ! Correct the orientation of the normal

            ! get the coordinate values for this element for each quadrature point
            x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)

            ! The node indices of the positions projected to the pressure mesh
            x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)
            
            ! Initialise array for the quadrature points of this 
            ! element for whether it has already been visited
            notvisited = .true.
            
            ! loop over local nodes within this element
            nodal_loop_i: do iloc = 1, di%pressure_mesh%shape%loc

              ! loop over CV faces internal to this element
              face_loop: do face = 1, di%cvfaces%faces

                ! is this a face neighbouring iloc?
                is_neigh: if(di%cvfaces%neiloc(iloc, face) /= 0) then

                  ! loop over gauss points on face
                  quadrature_loop: do gi = 1, di%cvfaces%shape%ngi

                    ! global gauss pt index for this element
                    ggi = (face-1)*di%cvfaces%shape%ngi + gi

                    ! check if this quadrature point has already been visited
                    check_visited: if(notvisited(ggi)) then

                       notvisited(ggi) = .false.

                       ! correct the orientation of the normal so it points away from iloc
                       normgi = orientate_cvsurf_normgi(node_val(di%positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                        x_face_quad(:,ggi), &
                                                        di%cached_face_value%normal(:,ggi,vele))
                       
                       di%cached_face_value%normal(:,ggi,vele) = normgi
                       
                    end if check_visited

                  end do quadrature_loop

                end if is_neigh

              end do face_loop

            end do nodal_loop_i
            
         end do vele_loop_detwei_normal

         sele_loop_detwei_normal: do sele = 1, di%number_sele

            call transform_cvsurf_facet_to_physical(ele_val(di%positions, face_ele(di%positions, sele)), &
                                                    face_val(di%positions, sele), &
                                                    di%x_cvbdyshape, &
                                                    di%cached_face_value%normal_bdy(:,:,sele), &
                                                    di%cached_face_value%detwei_bdy(:,sele))

         end do sele_loop_detwei_normal
      
      end if find_cached_detwei_normal

      deallocate(normgi)
      deallocate(notvisited)
      deallocate(x_face_quad)
      
      ! Determine the p_dshape if required
      find_cached_p_dshape: if (di%cached_face_value%cached_p_dshape) then

         vele_loop_p_dshape: do vele = 1, di%number_vele
            
            call transform_to_physical(di%positions, &
                                       vele, &
                                       x_shape = di%x_cvshape_full, &
                                       shape   = di%p_cvshape_full, &
                                       dshape  = di%cached_face_value%p_dshape(:,:,:,vele))

         end do vele_loop_p_dshape
         
!!!!! *** THIS IS NOT POSSIBLE YET *** 
!!!!!         sele_loop_p_dshape: do sele = 1, di%number_sele
            
!!!!!            call transform_to_physical(di%positions, &
!!!!!                                       sele, &
!!!!!                                       x_shape = di%x_cvbdyshape_full, &
!!!!!                                       shape   = di%p_cvbdyshape_full, &
!!!!!                                       dshape  = di%cached_face_value%p_dshape_bdy(:,:,:,sele))

!!!!!         end do sele_loop_p_dshape                 
         
      end if find_cached_p_dshape
      
      ewrite(1,*) 'Finished Initialise cached face value types'
      
   end subroutine darcy_impes_initialise_cached_face_value

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_relperm_den_first_face_values(di)
   
      !!< Find the relperm and density CV first face values for all phase and cache.
            
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      logical :: inflow
      integer :: p, vele, sele, iloc, oloc, face, gi, ggi
      integer :: r_upwind_pos, d_upwind_pos, dim
      real    :: income, v_over_relperm_dot_n, old_relperm_face_value, den_face_value
      real,    dimension(1)              :: absperm_ele, visc_ele
      real,    dimension(:),     pointer :: old_relperm_ele
      real,    dimension(:),     pointer :: den_ele
      real,    dimension(:,:),   pointer :: grav_ele
      real,    dimension(:,:),   pointer :: grad_pressure_face_quad
      real,    dimension(:,:),   pointer :: v_over_relperm_face_quad
      real,    dimension(:,:),   pointer :: x_ele
      real,    dimension(:,:,:), pointer :: p_dshape
      real,    dimension(:,:),   pointer :: normal
      real,    dimension(:),     pointer :: detwei
      real,    dimension(:),     pointer :: normgi
      logical, dimension(:),     pointer :: notvisited
      real,    dimension(:,:),   pointer :: x_face_quad
      integer, dimension(:),     pointer :: x_pmesh_nodes
      integer, dimension(:),     pointer :: p_nodes      
      integer, dimension(:),     pointer :: r_upwind_nodes
      integer, dimension(:),     pointer :: d_upwind_nodes
      real,    dimension(:),     pointer :: den_ele_bdy
      real,    dimension(:),     pointer :: old_relperm_ele_bdy
      integer, parameter :: V_BC_TYPE_NORMAL_FLOW = 1, V_BC_TYPE_NO_NORMAL_FLOW = 2
            
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      if (.not. di%cached_face_value%cached_p_dshape) then
         allocate(p_dshape(ele_loc(di%pressure_mesh,1), di%x_cvshape%ngi, mesh_dim(di%pressure_mesh)))
      end if
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(normal(di%ndim,di%x_cvshape%ngi))
         allocate(detwei(di%x_cvshape%ngi))      
      end if
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(old_relperm_ele(ele_loc(di%pressure_mesh,1)))
      allocate(den_ele(ele_loc(di%pressure_mesh,1)))
      allocate(grav_ele(di%ndim,1))
      allocate(grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(v_over_relperm_face_quad(di%ndim,di%p_cvshape%ngi))

      allocate(den_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(old_relperm_ele_bdy(face_loc(di%pressure_mesh,1)))      
      
      ewrite(1,*) 'Calculate relperm and density and first face values'      

      ! Initialise the cached face values
      di%cached_face_value%relperm     = 0.0
      di%cached_face_value%relperm_bdy = 0.0
      di%cached_face_value%den         = 0.0
      di%cached_face_value%den_bdy     = 0.0

      phase_loop: do p = 1, di%number_phase
          
         ! Determine the upwind relperm values of all node pairs if required for higher order CV face value
         if(need_upwind_values(di%modrelperm_cv_options)) then

           ! NOTE only new old values are required but this procedure 
           !      expects latest and old. So pass in old twice - not optimal
           call find_upwind_values(di%state, &
                                   di%positions_pressure_mesh, &
                                   di%old_relative_permeability(p)%ptr, &
                                   di%old_relperm_upwind, &
                                   di%old_relative_permeability(p)%ptr, &
                                   di%old_relperm_upwind, &
                                   option_path = trim(di%relative_permeability(p)%ptr%option_path))

         else

           call zero(di%modrelperm_upwind)

         end if
 
         ! Determine the upwind density values of all node pairs if required for higher order CV face value
         if(need_upwind_values(di%modrelperm_cv_options)) then

           ! NOTE only new upwind values are required but this procedure 
           !      expects latest and old. So pass in new twice - not optimal
           call find_upwind_values(di%state, &
                                   di%positions_pressure_mesh, &
                                   di%density(p)%ptr, &
                                   di%density_upwind, &
                                   di%density(p)%ptr, &
                                   di%density_upwind, &
                                   option_path = trim(di%density(p)%ptr%option_path))

         else

           call zero(di%density_upwind)

         end if
                          
         ! Initialise optimisation flag used in finding upwind value in high resolution schemes
         r_upwind_pos = 0
         d_upwind_pos = 0
            
         vol_element_loop: do vele = 1, di%number_vele
            
            ! The node indices of the pressure mesh
            p_nodes => ele_nodes(di%pressure_mesh, vele)

            ! The node indices of the positions projected to the pressure mesh
            x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)
                        
            ! Determine the node numbers to use to determine the relperm upwind values
            if((di%relperm_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
               (di%relperm_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

               r_upwind_nodes => x_pmesh_nodes

            else

               r_upwind_nodes => p_nodes

            end if
            
            ! Determine the node numbers to use to determine the density upwind values
            if((di%density_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
               (di%density_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

               d_upwind_nodes => x_pmesh_nodes

            else

               d_upwind_nodes => p_nodes

            end if

            ! get the old_relperm ele values for this phase
            old_relperm_ele = ele_val(di%old__relative_permeability(p)%ptr, vele)

            ! get the density value for this element for this phase
            den_ele = ele_val(di%density(p)%ptr, vele)
            
            ! The gravity values for this element for each direction
            grav_ele = ele_val(di%gravity, vele) 

            ! get the viscosity value for this element for this phase
            visc_ele = ele_val(di%viscosity(p)%ptr, vele)

            ! get the absolute permeability value for this element
            absperm_ele = ele_val(di%absolute_permeability, vele)         
            
            ! obtain the transformed determinant*weight and normals
            if (di%cached_face_value%cached_detwei_normal) then
               
               detwei => di%cached_face_value%detwei(:,vele)
               
               normal => di%cached_face_value%normal(:,:,vele)
               
            else

               ! get the coordinate values for this element for each positions local node
               x_ele = ele_val(di%positions, vele)         

               ! get the coordinate values for this element for each quadrature point
               x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)
            
               call transform_cvsurf_to_physical(x_ele, di%x_cvshape, detwei, normal, di%cvfaces)
            
            end if
            
            ! obtain the derivative of the pressure mesh shape function at the CV face quadrature points
            if (di%cached_face_value%cached_p_dshape) then
               
               p_dshape => di%cached_face_value%p_dshape(:,:,:,vele)
               
            else
            
               call transform_to_physical(di%positions, vele, x_shape = di%x_cvshape_full, &
                                          shape = di%p_cvshape_full, dshape = p_dshape)
            
            end if
            
            ! get the gradient pressure at the cv surface quadrature points for each direction for this phase
            grad_pressure_face_quad = darcy_impes_ele_grad_at_quad_scalar(di%pressure(p)%ptr, vele, dn = p_dshape)
                        
            ! the latest DarcyVelocity/relperm at the quadrature points
            ! determined from FE interpolation of each component, only used to determine upwind.
            do dim = 1,di%ndim

               v_over_relperm_face_quad(dim,:) = - (absperm_ele(1) / visc_ele(1)) * &
(grad_pressure_face_quad(dim,:) - ele_val_at_quad(di%density(p)%ptr, vele, di%p_cvshape) * grav_ele(dim,1))

            end do
            
            ! Initialise array for the quadrature points of this 
            ! element for whether it has already been visited
            notvisited = .true.
            
            ! loop over local nodes within this element
            nodal_loop_i: do iloc = 1, di%pressure_mesh%shape%loc

              ! loop over CV faces internal to this element
              face_loop: do face = 1, di%cvfaces%faces

                ! is this a face neighbouring iloc?
                is_neigh: if(di%cvfaces%neiloc(iloc, face) /= 0) then

                  ! find the opposing local node across the CV face
                  oloc = di%cvfaces%neiloc(iloc, face)

                  ! loop over gauss points on face
                  quadrature_loop: do gi = 1, di%cvfaces%shape%ngi

                    ! global gauss pt index for this element
                    ggi = (face-1)*di%cvfaces%shape%ngi + gi

                    ! check if this quadrature point has already been visited
                    check_visited: if(notvisited(ggi)) then

                       notvisited(ggi) = .false.

                       if (di%cached_face_value%cached_detwei_normal) then
                          
                          ! cached normal has correct orientation already
                          normgi = normal(:,ggi)
                          
                       else
                       
                          ! correct the orientation of the normal so it points away from iloc
                          normgi = orientate_cvsurf_normgi(node_val(di%positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                          &x_face_quad(:,ggi), normal(:,ggi))
                       
                       end if
                       
                       ! determine if the flow is in or out of the face at this quadrature
                       ! with respect to the normal orientation using the latest v_over_relperm
                       v_over_relperm_dot_n = dot_product(v_over_relperm_face_quad(:,ggi), normgi)

                       inflow = (v_over_relperm_dot_n<=0.0)

                       income = merge(1.0,0.0,inflow)

                       ! Evaluate the face value for old_relperm 
                       call darcy_impes_evaluate_face_val(old_relperm_face_value, &
                                                          iloc, &
                                                          oloc, &
                                                          ggi, &
                                                          r_upwind_nodes, &
                                                          di%p_cvshape, &
                                                          old_relperm_ele, &
                                                          di%old_relperm_upwind, &
                                                          inflow, &
                                                          income, &
                                                          di%relperm_cv_options, &
                                                          save_pos = r_upwind_pos)
                       
                       ! Evaluate the density face value 
                       call darcy_impes_evaluate_face_val(den_face_value, &
                                                          iloc, &
                                                          oloc, &
                                                          ggi, &
                                                          d_upwind_nodes, &
                                                          di%p_cvshape, &
                                                          den_ele, &
                                                          di%density_upwind, &
                                                          inflow, &
                                                          income, &
                                                          di%density_cv_options, &
                                                          save_pos = d_upwind_pos)
                       
                       di%cached_face_value%modrelperm(1,ggi,vele,p) = old_relperm_face_value
                      
                       di%cached_face_value%den(ggi,vele,p) = den_face_value
                                           
                    end if check_visited

                  end do quadrature_loop

                end if is_neigh

              end do face_loop

            end do nodal_loop_i

         end do vol_element_loop
                  
         ! Get this phase v BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
         call get_v_boundary_condition(di%darcy_velocity(p)%ptr, &
                                       (/"normal_flow   ", &
                                         "no_normal_flow"/), &
                                       di%bc_surface_mesh, &
                                       di%v_over_s_bc_value, &
                                       di%v_over_s_bc_flag)
         
         sele_loop: do sele = 1, di%number_sele
            
            ! A no_normal_flow BC adds nothing to the matrix and rhs so cycle sele loop
            if (di%v_bc_flag(sele) == V_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop
            
            ! A normal_flow BC adds only the given value - no need for cached face values          
            if (di%v_bc_flag(sele) == V_BC_TYPE_NORMAL_FLOW) cycle sele_loop
            
            ! get the density sele values for this phase
            den_ele_bdy = face_val(di%density(p)%ptr, sele)
            
            ! get the relperm domain value
            old_relperm_ele_bdy = face_val(di%old_relative_permeability(p)%ptr, sele)
                        
            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi
                                                       
                           ! Evaluate the old_relperm and density face value via taking the CV value
                           ! - no other choice currently ...
                           
                           di%cached_face_value%relperm_bdy(1,ggi,sele,p) = old_relperm_ele_bdy(iloc)

                           di%cached_face_value%den_bdy(ggi,sele,p) = den_ele_bdy(iloc)
                           
                        end if
                        
                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop

         end do sele_loop      

      end do phase_loop
      
      ! deallocate local variables as required
      deallocate(x_ele)
      if (di%cached_face_value%cached_p_dshape) then
         nullify(p_dshape)      
      else
         deallocate(p_dshape)
      end if
      if (di%cached_face_value%cached_detwei_normal) then
         nullify(detwei)
         nullify(normal)      
      else
         deallocate(detwei)
         deallocate(normal)
      end if
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(x_face_quad)
      deallocate(old_relperm_ele)
      deallocate(den_ele)
      deallocate(grav_ele)
      deallocate(grad_pressure_face_quad)
      deallocate(v_over_relperm_face_quad)

      deallocate(den_ele_bdy)
      deallocate(old_relperm_ele_bdy)
      
      ewrite(1,*) 'Finished calculate relperm and density first face values'      
      
   end subroutine darcy_impes_calculate_relperm_den_first_face_values

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_inverse_characteristic_length(di)
   
      !!< Calculate the inverse_characteristic_length variable
      
      ! This is only set up for linear mesh
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variable
      integer :: vele, sele
      real    :: subcv_area, subcv_volume
      real, dimension(:), pointer :: len_vals

      ewrite(1,*) 'Calculate inverse characteristic length'
      
      allocate(len_vals(face_loc(di%pressure_mesh,1)))
      
      sele_loop: do sele = 1, di%number_sele
         
         ! Find corresponding volume element for this surface facet
         vele = face_ele(di%pressure_mesh, sele)
         
         ! Find the subcv volumes assuming linear mesh
         subcv_volume = element_volume(di%positions, vele) / real(ele_loc(di%pressure_mesh, vele))
         
         ! Find the subcv areas assuming linear mesh         
         subcv_area = surface_element_area(di%positions, sele) / real(face_loc(di%pressure_mesh, sele))
         
         ! Find characteristic length assuming linear mesh
         len_vals = subcv_area / subcv_volume
         
         ! Set the characteristic length for the local dg nodes
         call set(di%inverse_characteristic_length, &
                  ele_nodes(di%inverse_characteristic_length, sele), &
                  len_vals)
         
      end do sele_loop
      
      deallocate(len_vals)
      
      call invert(di%inverse_characteristic_length)
      
      ewrite_minmax(di%inverse_characteristic_length)
      
      ewrite(1,*) 'Finished Calculate inverse characteristic length'
      
   end subroutine darcy_impes_calculate_inverse_characteristic_length

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_evaluate_face_val()
   
      !!< 
            
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variable

                  
   end subroutine darcy_impes_evaluate_face_val

! ----------------------------------------------------------------------------

  function darcy_impes_ele_grad_at_quad_scalar(field, ele_number, dn) result (quad_grad)
    
    !!< Return the grad of field at the quadrature points of
    !!< ele_number. dn is the transformed element gradient.
    
    type(scalar_field),intent(in) :: field
    integer, intent(in) :: ele_number
    real, dimension(:,:,:), intent(in) :: dn
    real, dimension(mesh_dim(field), size(dn,2)) :: quad_grad
    
    ! local variables
    integer :: i
    
    ! sanity check
    assert(field%mesh%shape%loc == size(dn,1))
    assert(field%mesh%shape%dim == size(dn,3))
    
    do i=1, mesh_dim(field)
       quad_grad(i,:) = matmul(ele_val(field, ele_number),dn(:,:,i))
    end do
    
  end function darcy_impes_ele_grad_at_quad_scalar
   
! ----------------------------------------------------------------------------

end module darcy_impes_assemble_module
