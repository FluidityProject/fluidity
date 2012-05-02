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
             darcy_impes_assemble_and_solve, &
             darcy_impes_calculate_gradient_pressures, &
             darcy_impes_calculate_non_first_phase_pressures, &
             darcy_impes_calculate_phase_one_saturation_diagnostic, &
             darcy_impes_calculate_velocity_and_cfl_fields, &
             darcy_impes_calculate_sum_saturation, &
             darcy_impes_calculate_densities, &
             darcy_impes_calculate_cflnumber_field_based_dt, &
             darcy_impes_initialise_cached_phase_face_value, &
             darcy_impes_calculate_inverse_characteristic_length
   
   ! Options associated with explicit advection subcycling for CV scalar field solver
   type darcy_impes_subcycle_options_type
      logical :: have_max_cfl
      logical :: have_number
      integer :: number_advection_subcycle
      real    :: max_courant_per_advection_subcycle         
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
   
   ! Data associated with a cached CV face value
   type cached_face_value_type
      real,    dimension(:), pointer :: value
      integer, dimension(:), pointer :: ele_base
   end type cached_face_value_type
   
   type darcy_impes_type
      ! *** Pointers to fields from state that have array length of number of phases ***
      type(vector_field_pointer), dimension(:), pointer :: inter_velocity_porosity
      type(vector_field_pointer), dimension(:), pointer :: old_inter_velocity_porosity
      type(vector_field_pointer), dimension(:), pointer :: darcy_velocity   
      type(vector_field_pointer), dimension(:), pointer :: fractional_flow
      type(scalar_field_pointer), dimension(:), pointer :: saturation
      type(scalar_field_pointer), dimension(:), pointer :: old_saturation
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
      type(scalar_field_pointer), dimension(:), pointer :: density
      type(scalar_field_pointer), dimension(:), pointer :: old_density
      ! *** Pointers to fields from state that are NOT phase dependent ***
      type(mesh_type),    pointer :: pressure_mesh
      type(mesh_type),    pointer :: velocity_mesh      
      type(mesh_type),    pointer :: elementwise_mesh
      type(scalar_field), pointer :: average_pressure
      type(scalar_field), pointer :: porosity
      type(scalar_field), pointer :: old_porosity
      type(scalar_field), pointer :: absolute_permeability
      type(scalar_field), pointer :: old_absolute_permeability
      type(vector_field), pointer :: positions
      type(vector_field), pointer :: total_darcy_velocity   
      type(scalar_field), pointer :: sum_saturation
      type(scalar_field), pointer :: old_sum_saturation
      type(scalar_field), pointer :: div_total_darcy_velocity
      type(vector_field), pointer :: gravity_direction
      ! *** The gravity magnitude ***
      real :: gravity_magnitude
      ! *** The full gravity field - direction * magnitude, allocated here ***
      type(vector_field) :: gravity
      ! *** Fields allocated here used in assemble algorithm ***
      type(vector_field) :: positions_pressure_mesh
      type(csr_matrix)   :: pressure_matrix
      type(scalar_field) :: lhs
      type(scalar_field) :: rhs
      type(scalar_field) :: rhs_adv
      type(scalar_field) :: rhs_time
      type(scalar_field) :: inverse_cv_mass_cfl_mesh
      type(scalar_field) :: inverse_cv_mass_pressure_mesh
      type(scalar_field) :: cv_mass_pressure_mesh_with_porosity   
      type(scalar_field) :: cv_mass_pressure_mesh_with_old_porosity 
      type(scalar_field) :: cfl_subcycle
      type(scalar_field) :: old_sfield_subcycle   
      ! *** Data associated with vphi, pressure and saturation BC allocated here ***
      type(mesh_type)                          :: bc_surface_mesh
      type(scalar_field)                       :: vphi_bc_value
      integer,           dimension(:), pointer :: vphi_bc_flag
      type(scalar_field)                       :: pressure_bc_value
      integer,           dimension(:), pointer :: pressure_bc_flag
      type(scalar_field)                       :: saturation_bc_value
      integer,           dimension(:), pointer :: saturation_bc_flag
      type(scalar_field)                       :: inverse_characteristic_length
      real                                     :: weak_pressure_bc_coeff
      ! *** The pressure mesh - pressure mesh sparsity, used for pressure matrix and finding CV upwind values ***
      type(csr_sparsity), pointer :: sparsity_pmesh_pmesh
      ! *** The number of phase and the CV surface quadrature degree to use ***
      integer :: number_phase, quaddegree   
      ! *** Data specifically associated with the CV discretisation allocated here ***
      type(cv_options_type)                   :: saturation_cv_options
      type(cv_faces_type)                     :: cvfaces
      type(element_type)                      :: x_cvshape_full
      type(element_type)                      :: p_cvshape_full
      type(element_type)                      :: x_cvshape
      type(element_type)                      :: p_cvshape
      type(element_type)                      :: x_cvbdyshape
      type(element_type)                      :: p_cvbdyshape      
      logical                                 :: phase_one_saturation_diagnostic
      type(darcy_impes_subcycle_options_type) :: subcy_opt_sat
      type(csr_matrix)                        :: old_sfield_upwind
      ! *** Flag for whether the first phase pressure is prognostic, else it is prescribed *** 
      logical :: first_phase_pressure_prognostic
      ! *** Flag for whether the Porosity is diagnostic, else it is prescribed *** 
      logical :: porosity_is_diagnostic
      ! *** Flag for whether the AbsolutePermeability is diagnostic, else it is prescribed *** 
      logical :: absolute_permeability_is_diagnostic
      ! *** The cached phase face value at each quadrature point summed over subcycles if necessary for each phase ***
      type(cached_face_value_type), dimension(:), pointer :: cached_phase_face_value_domain
      type(cached_face_value_type), dimension(:), pointer :: cached_phase_face_value_boundary
      ! *** Time time step size, also stored here for convenience ***
      real :: dt
      ! *** The current time, also stored here for convenience ***
      real :: current_time 
      ! *** Non linear iteration, also stored here for convencience ***
      integer :: nonlinear_iter
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

   subroutine darcy_impes_assemble_and_solve(di)
      
      !!< Assemble and solve the Darcy equations using an CIMPESS algorithm
      !!< which is a modification of the IMPES to include consistent subcycling.
      
      type(darcy_impes_type), intent(inout) :: di
      
      ewrite(1,*) 'Start Darcy IMPES assemble and solve'
      
      ! Form the full gravity field
      call set(di%gravity, di%gravity_direction)
      call scale(di%gravity, di%gravity_magnitude)
      
      ! Copy to Old the CV mass on the pressure mesh with porosity included
      call set(di%cv_mass_pressure_mesh_with_old_porosity, di%cv_mass_pressure_mesh_with_porosity)
      
      ! Calculate the latest CV mass on the pressure mesh with porosity included
      call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_porosity, di%porosity)      
            
      ! Assemble and solve the phase pressures
      call darcy_impes_assemble_and_solve_phase_pressures(di)
            
      ! Calculate the gradient pressures 
      call darcy_impes_calculate_gradient_pressures(di)

      ! Calculate the divergence total darcy velocity
      call darcy_impes_calculate_divergence_total_darcy_velocity(di)
      
      ! Assemble and solve the phase saturations
      call darcy_impes_assemble_and_solve_phase_saturations(di)
      
      ! Calculate the sum of the saturations
      call darcy_impes_calculate_sum_saturation(di)
      
      ! Calculate the density field of each phase
      call darcy_impes_calculate_densities(di)
            
      ewrite(1,*) 'Finished Darcy IMPES assemble and solve'
           
   end subroutine darcy_impes_assemble_and_solve

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
      
      !!< Assemble and solve the first phase pressure. For the first nonlinear iteration
      !!< the phase face value is calculated here and cached to be used 
      !!< for the calculation of total darcy velocity divergence and volume fraction  
      !!< advection. Any further nonlinear iterations use the cached face value 
      !!< which may contain the summed subcycled face values. A source is included 
      !!< due to the capilliary pressures of non first phases as well as gravity.
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p, vele, sele, iloc, oloc, jloc, face, gi, ggi, f_ele_counter, f_ele_base, upwind_pos, dim
      real    :: income, face_value, vphi_dot_n, old_relperm_absperm_over_visc
      real    :: old_saturation_face_value, old_relperm_face_value
      logical :: inflow, determine_face_value
      real,    dimension(1)                  :: old_absperm_ele, old_visc_ele
      real,    dimension(:,:),   allocatable :: grad_pressure_ele
      real,    dimension(:,:),   allocatable :: vphi_face
      real,    dimension(:),     allocatable :: old_saturation_ele
      real,    dimension(:),     allocatable :: old_relperm_ele
      real,    dimension(:),     allocatable :: cfl_ele
      real,    dimension(:,:),   allocatable :: x_ele
      real,    dimension(:,:,:), allocatable :: p_dshape
      real,    dimension(:,:),   allocatable :: normal
      real,    dimension(:),     allocatable :: detwei
      real,    dimension(:),     allocatable :: normgi
      logical, dimension(:),     allocatable :: notvisited
      real,    dimension(:,:),   allocatable :: p_mat_local
      real,    dimension(:,:),   allocatable :: x_face_quad
      integer, dimension(:),     pointer     :: x_pmesh_nodes
      integer, dimension(:),     pointer     :: p_nodes      
      integer, dimension(:),     pointer     :: upwind_nodes
      real,    dimension(1)                  :: old_absperm_ele_bdy, old_visc_ele_bdy
      real,    dimension(:),     allocatable :: old_relperm_ele_bdy
      real,    dimension(:,:),   allocatable :: grad_pressure_ele_bdy
      real,    dimension(:),     allocatable :: bc_sele_val
      real,    dimension(:),     allocatable :: inv_char_len_ele_bdy
      real,    dimension(:),     allocatable :: old_saturation_ele_bdy
      real,    dimension(:),     allocatable :: ghost_old_saturation_ele_bdy
      real,    dimension(:,:),   allocatable :: normal_bdy
      real,    dimension(:),     allocatable :: detwei_bdy
      real,    dimension(:,:),   allocatable :: x_ele_bdy
      real,    dimension(:),     allocatable :: p_rhs_local_bdy
      real,    dimension(:,:),   allocatable :: p_matrix_local_bdy
      integer, dimension(:),     allocatable :: p_nodes_bdy
      integer, parameter :: PRESSURE_BC_TYPE_WEAKDIRICHLET   = 1
      integer, parameter :: SATURATION_BC_TYPE_WEAKDIRICHLET = 1
      integer, parameter :: VPHI_BC_TYPE_NORMAL_FLOW         = 1, VPHI_BC_TYPE_NO_NORMAL_FLOW = 2
            
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      allocate(p_dshape(ele_loc(di%pressure_mesh,1), di%x_cvshape%ngi, mesh_dim(di%pressure_mesh)))
      allocate(normal(di%ndim,di%x_cvshape%ngi))
      allocate(detwei(di%x_cvshape%ngi))      
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(p_mat_local(ele_loc(di%pressure_mesh,1), ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(old_saturation_ele(ele_loc(di%pressure_mesh,1)))
      allocate(old_relperm_ele(ele_loc(di%pressure_mesh,1)))
      allocate(vphi_face(di%ndim,di%p_cvshape%ngi))
      allocate(cfl_ele(ele_loc(di%pressure_mesh,1)))
      allocate(grad_pressure_ele(di%ndim,1))

      allocate(old_relperm_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(grad_pressure_ele_bdy(di%ndim,1))
      allocate(old_saturation_ele_bdy(face_loc(di%pressure_mesh,1)))      
      allocate(ghost_old_saturation_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(bc_sele_val(face_loc(di%pressure_mesh,1)))
      allocate(inv_char_len_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(detwei_bdy(di%x_cvbdyshape%ngi))
      allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      allocate(p_rhs_local_bdy(face_loc(di%pressure_mesh,1)))
      allocate(p_matrix_local_bdy(face_loc(di%pressure_mesh,1),face_loc(di%pressure_mesh,1)))
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(p_nodes_bdy(face_loc(di%pressure_mesh,1)))
      
      ewrite(1,*) 'Solve first phase Pressure'
      
      ! Initialise the matrix and rhs
      call zero(di%pressure_matrix)
      call zero(di%rhs)
      
      ! Decide if a face value needs to be determined
      determine_face_value = di%nonlinear_iter == 1
           
      ! Assemble a contribution from each phase to form a global continuity equation to solve for first phase pressure
      phase_loop: do p = 1, di%number_phase
 
         if (determine_face_value) then

            ! Initialise the cached phase face values
            di%cached_phase_face_value_domain(p)%value   = 0.0
            di%cached_phase_face_value_boundary(p)%value = 0.0
         
            ! Determine the upwind saturation values if required for higher order CV face value
            if(need_upwind_values(di%saturation_cv_options)) then

              call find_upwind_values(di%state, &
                                      di%positions_pressure_mesh, &
                                      di%old_saturation(p)%ptr, &
                                      di%old_sfield_upwind, &
                                      di%old_saturation(p)%ptr, &
                                      di%old_sfield_upwind, &
                                      option_path = trim(di%saturation(p)%ptr%option_path))

            else

              call zero(di%old_sfield_upwind)

            end if
         
         end if
                  
         ! Initialise optimisation flag used in finding upwind value in high resolution schemes
         upwind_pos = 0
            
         ! Loop volume elements assembling local contributions     
         vol_element_loop: do vele = 1,element_count(di%pressure_mesh)

            ! Initialise the face value counter for this ele
            f_ele_counter = 0
            
            ! Find the face value ele base
            f_ele_base = di%cached_phase_face_value_domain(p)%ele_base(vele)
            
            ! get the coordinate values for this element for each positions local node
            x_ele = ele_val(di%positions, vele)         

            ! get the coordinate values for this element for each quadrature point
            x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)

            ! The node indices of the pressure mesh
            p_nodes => ele_nodes(di%pressure_mesh, vele)

            ! The node indices of the positions projected to the pressure mesh
            x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)

            if (determine_face_value) then

               ! get the old_relperm ele values for this phase
               old_relperm_ele = ele_val(di%old_relative_permeability(p)%ptr, vele)

               ! get the old_viscosity value for this element for this phase
               old_visc_ele = ele_val(di%old_viscosity(p)%ptr, vele)

               ! get the old_absolute permeability value for this element
               old_absperm_ele = ele_val(di%old_absolute_permeability, vele)         

               ! get the gradient pressure value for this element
               grad_pressure_ele = ele_val(di%gradient_pressure(p)%ptr, vele)         

               ! get the old saturation ele values for this phase
               old_saturation_ele = ele_val(di%old_saturation(p)%ptr, vele)
            
               ! get the CFL values for this element
               cfl_ele = ele_val(di%cfl(p)%ptr, vele)
            
               ! Determine the node numbers to use to determine the
               ! Saturation and RelPerm upwind values
               if((di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
                  (di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

                  upwind_nodes => x_pmesh_nodes

               else

                  upwind_nodes => p_nodes

               end if
            
               ! the latest intersitial velocity_porosity at the quadrature points, used to determine upwind
               do dim = 1,di%ndim

                  vphi_face(dim,:) =  - ele_val_at_quad(di%old_relative_permeability(p)%ptr, vele, di%p_cvshape) * &
                                        old_absperm_ele(1) * grad_pressure_ele(dim,1) / old_visc_ele(1)     

               end do
            
            end if
            
            ! obtain the transformed determinant*weight and normals
            call transform_cvsurf_to_physical(x_ele, di%x_cvshape, detwei, normal, di%cvfaces)

            ! obtain the derivative of the pressure mesh shape function at the CV face quadrature points
            call transform_to_physical(di%positions, vele, x_shape = di%x_cvshape_full, &
                                       shape = di%p_cvshape_full, dshape = p_dshape)
            
            ! Initialise array for the quadrature points of this 
            ! element for whether it has already been visited
            notvisited = .true.

            ! Initialise the local p matrix to assemble for this element
            p_mat_local = 0.0

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
                      
                      f_ele_counter = f_ele_counter + 1
                      
                      ! correct the orientation of the normal so it points away from iloc
                      normgi = orientate_cvsurf_normgi(node_val(di%positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                      &x_face_quad(:,ggi), normal(:,ggi))
                      
                      if (determine_face_value) then

                         ! determine if the flow is in or out of the face at this quadrature
                         ! with respect to the normal orientation using the latest vphi
                         vphi_dot_n = dot_product(vphi_face(:,ggi), normgi(:))

                         inflow = (vphi_dot_n<=0.0)

                         income = merge(1.0,0.0,inflow)
                      
                         ! evaluate the nonlinear face value for saturation
                         call evaluate_face_val(old_saturation_face_value, &
                                                old_saturation_face_value, & 
                                                iloc, &
                                                oloc, &
                                                ggi, &
                                                upwind_nodes, &
                                                di%p_cvshape, &
                                                old_saturation_ele, &
                                                old_saturation_ele, &
                                                di%old_sfield_upwind, &
                                                di%old_sfield_upwind, &
                                                inflow, &
                                                cfl_ele, &
                                                di%saturation_cv_options, &
                                                save_pos = upwind_pos)
                      
                         ! Evaluate the face value for old_relperm (assuming upwind for now) 
                         old_relperm_face_value = income*old_relperm_ele(oloc) + (1.0-income)*old_relperm_ele(iloc)

                         ! face value = S*old_relperm*old_absperm/old_visc, where old_absperm is phase independent
                         ! (if absperm and viscosity are considered tensors this requires modifying below)
                         face_value = old_saturation_face_value * old_relperm_face_value * old_absperm_ele(1) / old_visc_ele(1)

                         ! cache the phase face value for consistent use later
                         di%cached_phase_face_value_domain(p)%value((f_ele_base - 1) + f_ele_counter) = face_value

                      else 
                      
                         face_value = di%cached_phase_face_value_domain(p)%value((f_ele_base - 1) + f_ele_counter)
                      
                      end if 
                      
                      face_value = face_value * detwei(ggi)
                      
                      ! Form the local matrix given by - n_i . sum_{phase} ( S*old_relperm*old_absperm/old_visc ) dP/dx_j
                      do jloc = 1,di%pressure_mesh%shape%loc

                         p_mat_local(iloc,jloc) = p_mat_local(iloc,jloc) - &
                                                  sum(p_dshape(jloc, ggi, :)*normgi, 1)*face_value

                         p_mat_local(oloc,jloc) = p_mat_local(oloc,jloc) - &
                                                  sum(p_dshape(jloc, ggi, :)*(-normgi), 1)*face_value

                      end do

                    end if check_visited

                  end do quadrature_loop

                end if is_neigh

              end do face_loop

            end do nodal_loop_i

            ! Add volume element contribution to global pressure matrix
            call addto(di%pressure_matrix, p_nodes, p_nodes, p_mat_local)

         end do vol_element_loop

      end do phase_loop
      
      ! Add rate of change of porosity to rhs    
      call set(di%rhs, di%cv_mass_pressure_mesh_with_old_porosity)
                  
      call addto(di%rhs, di%cv_mass_pressure_mesh_with_porosity, scale = -1.0)
      
      call scale(di%rhs, 1.0/di%dt)
                  
      ! Add normal vphi BC integrals for each phase, else include if required 
      ! a weak pressure BC (which if not given for phase 1 is assumed zero).
      ! This all requires a value  of the saturation of each phase on the boundary 
      ! which must take account of any weak BC it has.

      call get_entire_saturation_or_pressure_boundary_condition(di%pressure(1)%ptr, &
                                                                (/"weakdirichlet"/), &
                                                                di%bc_surface_mesh, &
                                                                di%pressure_bc_value, &
                                                                di%pressure_bc_flag)
      phase_loop_bc: do p = 1, di%number_phase
         
         ! Get this phase vphi BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
         call get_vphi_boundary_condition(di%inter_velocity_porosity(p)%ptr, &
                                          (/"normal_flow   ", &
                                            "no_normal_flow"/), &
                                          di%bc_surface_mesh, &
                                          di%vphi_bc_value, &
                                          di%vphi_bc_flag)
         
         ! Get this phase old_satutation BC info - only weak dirichlet
         call get_entire_saturation_or_pressure_boundary_condition(di%saturation(p)%ptr, &
                                                                   (/"weakdirichlet"/), &
                                                                   di%bc_surface_mesh, &
                                                                   di%saturation_bc_value, &
                                                                   di%saturation_bc_flag)
         
         vphi_sele_loop: do sele = 1, surface_element_count(di%pressure_mesh)
            
            ! A no_normal_flow BC adds nothing to the matrix and rhs so cycle sele loop
            if (di%vphi_bc_flag(sele) == VPHI_BC_TYPE_NO_NORMAL_FLOW) cycle vphi_sele_loop
                     
            ! Initialise the face value counter for this ele
            f_ele_counter = 0

            ! Find the face value ele base
            f_ele_base = di%cached_phase_face_value_boundary(p)%ele_base(sele)
            
            if (di%vphi_bc_flag(sele) == VPHI_BC_TYPE_NORMAL_FLOW) then
                              
               if (determine_face_value) then
               
                  bc_sele_val = ele_val(di%vphi_bc_value, sele)
                                                   
               end if
            
            else
            
               ! get the old_relperm sele values for this phase
               old_relperm_ele_bdy = face_val(di%old_relative_permeability(p)%ptr, sele)

               ! get the old_viscosity value for this sele for this phase
               old_visc_ele_bdy = face_val(di%old_viscosity(p)%ptr, sele)

               ! get the old_absolute permeability value for this sele
               old_absperm_ele_bdy = face_val(di%old_absolute_permeability, sele)         
               
               ! get the latest gradient pressure for this sele 
               ! - used to determine upwind direction and for the bc integral (really this should be implicit)
               grad_pressure_ele_bdy = face_val(di%gradient_pressure(p)%ptr, sele)
               
               ! the inverse characteristic length used for pressure weak bc
               inv_char_len_ele_bdy = ele_val(di%inverse_characteristic_length, sele)
                              
               ! Get the pressure BC value if required
               if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then
               
                  bc_sele_val = ele_val(di%pressure_bc_value, sele)
               
               end if
               
            end if 

            ! get the saturation domain and boundary (in case weak) on the boundary
            old_saturation_ele_bdy = face_val(di%old_saturation(p)%ptr, sele)

            if (di%saturation_bc_flag(sele) == SATURATION_BC_TYPE_WEAKDIRICHLET) then
               ghost_old_saturation_ele_bdy = ele_val(di%saturation_bc_value, sele)
            else
               ghost_old_saturation_ele_bdy = face_val(di%old_saturation(p)%ptr, sele)
            end if
               
            x_ele       = ele_val(di%positions, face_ele(di%positions, sele))
            x_ele_bdy   = face_val(di%positions, sele)
            p_nodes_bdy = face_global_nodes(di%pressure_mesh, sele)

            call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)

            p_rhs_local_bdy    = 0.0
            p_matrix_local_bdy = 0.0

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi

                        f_ele_counter = f_ele_counter + 1
                        
                        
                        if (di%vphi_bc_flag(sele) == VPHI_BC_TYPE_NORMAL_FLOW) then
                           
                           ! If have normal_flow BC for this phase then include CV integral of value * saturation
                           
                           if (determine_face_value) then

                              ! Determine the upwind saturation boundary values to use from the sign 
                              ! of the normal flow BC. If positive, hence outflow, always 
                              ! use domain value. If negative, hence inflow, use the boundary 
                              ! value which is the weak BC value if appropriate

                              if (bc_sele_val(iloc) > 0.0) then

                                 face_value = bc_sele_val(iloc) * old_saturation_ele_bdy(iloc)
                              
                              else

                                 face_value = bc_sele_val(iloc) * ghost_old_saturation_ele_bdy(iloc)
                              
                              end if 
                            
                              ! cache the phase face value for consistent use later
                              di%cached_phase_face_value_boundary(p)%value((f_ele_base - 1) + f_ele_counter) = face_value
                              
                           else
                              
                              face_value = di%cached_phase_face_value_boundary(p)%value((f_ele_base - 1) + f_ele_counter)
                           
                           end if
                           
                           p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - face_value * detwei_bdy(ggi)
                        
                        else 
                           
                           ! Else apply a weak pressure BC, which has a value of zero by default (if not given)
                           ! - NOTE for non first phases P = P_1 + P_C, and P_C = 0.0 on boundary so
                           !        all phases pressure boundary conditions related to phase 1
                           
                           old_relperm_absperm_over_visc = old_relperm_ele_bdy(iloc) * old_absperm_ele_bdy(1) / old_visc_ele_bdy(1)
                           
                           vphi_dot_n = - old_relperm_absperm_over_visc * dot_product(grad_pressure_ele_bdy(:,1), normal_bdy(:,ggi))

                           inflow = (vphi_dot_n<=0.0)

                           income = merge(1.0,0.0,inflow)
                           
                           ! Find upwind old saturation face value
                           old_saturation_face_value = income*ghost_old_saturation_ele_bdy(iloc) + &
                                                      &(1.0-income)*old_saturation_ele_bdy(iloc)
                           
                           ! If have phase pressure BC then include integral of BC value in rhs
                           ! - if not given then this implies adding a zero weak pressure BC
                           if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then
                                                            
                              do jloc = 1,di%pressure_mesh%faces%shape%loc 
                              
                                 p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - &
                                                         inv_char_len_ele_bdy(jloc) *&
                                                         di%p_cvbdyshape%n(jloc,ggi) * &
                                                         bc_sele_val(jloc) * &
                                                         old_relperm_absperm_over_visc * &
                                                         old_saturation_face_value * &
                                                         sum(normal_bdy(:,ggi)) * &
                                                         detwei_bdy(ggi) * &
                                                         di%weak_pressure_bc_coeff
                                                            
                              end do
                              
                           end if 
                           
                           ! Include weak pressure BC implicit term in matrix
                           ! - if no pressure BC adding this term implies a 
                           !   a default of weak zero pressure BC
                           do jloc = 1,di%pressure_mesh%faces%shape%loc 

                              p_matrix_local_bdy(iloc,jloc) = p_matrix_local_bdy(iloc,jloc) - &
                                                              inv_char_len_ele_bdy(jloc) *&
                                                              di%p_cvbdyshape%n(jloc,ggi) * &
                                                              old_relperm_absperm_over_visc * &
                                                              old_saturation_face_value * &
                                                              sum(normal_bdy(:,ggi)) * &
                                                              detwei_bdy(ggi) * &
                                                              di%weak_pressure_bc_coeff

                           end do
                           
                           ! Add coeff *gradient of pressure term to rhs via vphi_dot_n 
                           !  - which may include gravity and capilliary pressure term already                           
                           ! (to be fully implicit this should be added to matrix but 2 non linear iterations may save us ...)
                           p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - vphi_dot_n * detwei_bdy(ggi)
                                                      
                        end if
                        
                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop

            call addto(di%rhs, p_nodes_bdy, p_rhs_local_bdy)
            
            call addto(di%pressure_matrix, p_nodes_bdy, p_nodes_bdy, p_matrix_local_bdy)

         end do vphi_sele_loop      

      end do phase_loop_bc
            
      ! Apply any strong dirichlet BC's
      call apply_dirichlet_conditions(di%pressure_matrix, di%rhs, di%pressure(1)%ptr)
      
      ! Solve the pressure
      call petsc_solve(di%pressure(1)%ptr, di%pressure_matrix, di%rhs, di%state(1))

      ! deallocate local variables as required
      deallocate(x_ele)
      deallocate(p_dshape)
      deallocate(normal)
      deallocate(detwei)
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(p_mat_local)
      deallocate(x_face_quad)
      deallocate(old_saturation_ele)
      deallocate(old_relperm_ele)
      deallocate(vphi_face)
      deallocate(cfl_ele)
      deallocate(grad_pressure_ele)

      deallocate(old_relperm_ele_bdy)
      deallocate(grad_pressure_ele_bdy)
      deallocate(old_saturation_ele_bdy)
      deallocate(ghost_old_saturation_ele_bdy)
      deallocate(bc_sele_val)
      deallocate(inv_char_len_ele_bdy)
      deallocate(detwei_bdy)
      deallocate(normal_bdy)
      deallocate(p_rhs_local_bdy)
      deallocate(p_matrix_local_bdy)
      deallocate(x_ele_bdy)      
      deallocate(p_nodes_bdy)
      
      ewrite(1,*) 'Finished solve first phase Pressure'
      
      ewrite_minmax(di%pressure(1)%ptr)
      
   end subroutine darcy_impes_assemble_and_solve_first_phase_pressure

! ----------------------------------------------------------------------------

   subroutine get_vphi_boundary_condition(vphi, &
                                          types, &
                                          bc_surface_mesh, &
                                          vphi_bc_value, &
                                          vphi_bc_flag)

      !!< Form the data associated with any vphi BC by returning 
      !!< full surface mesh arrays of a flag indicating  either 
      !!< no_normal_flow or normal_flow and the value.

      type(vector_field),               intent(in),   target :: vphi
      character(len=*),   dimension(:), intent(in)           :: types
      type(mesh_type),                  intent(in)           :: bc_surface_mesh
      type(scalar_field),               intent(inout)        :: vphi_bc_value
      integer,            dimension(:), intent(inout)        :: vphi_bc_flag

      ! Local variables
      character(len=FIELD_NAME_LEN)                       :: bctype
      type(scalar_field),                         pointer :: scalar_surface_field
      integer,                      dimension(:), pointer :: surface_element_list
      integer                                             :: i, j, k, sele

      ewrite(1,*) 'Get InterstitialVelocityPorosity boundary condition data'

      ! Zero the normal_flow value surface field on whole boundary mesh  
      call zero(vphi_bc_value)

      ! Initialise flag for whether surface element has normal_flow BC
      vphi_bc_flag = 0

      ! Loop each BC object instance for the vphi
      ! (May have multiple normal_flow BC's applied to different surface id's)
      BC_loop: do i=1, get_boundary_condition_count(vphi)

         ! Get this BC info
         call get_boundary_condition(vphi, i, type = bctype, &
                                     surface_element_list = surface_element_list)

         ! check this is a normal_flow BC  or no_normal_flow (nothing else is permitted)
         if ((trim(bctype) /= 'normal_flow') .and. (trim(bctype) /= 'no_normal_flow')) then
            FLAbort('Have unknown BC type for a InterstitialvphiPorosity')
         end if

         ! Extract the scalar_surface_field for this BC for normal_flow
         if (trim(bctype) == 'normal_flow') then
            if (associated(vphi%bc%boundary_condition(i)%scalar_surface_fields)) then
               scalar_surface_field => vphi%bc%boundary_condition(i)%scalar_surface_fields(1)
            else
               FLAbort('Component scalar_surface_fields for InterstitialVelocityPorosity BC type not associated')
            end if
         end if 
         
         ! Loop the surface elements associated with this BC instance
         ! and place the required BC value in whole boundary field
         ! for the normal_flow BC type
         BC_sele_loop: do k = 1, size(surface_element_list)

            ! Find the whole domain surface element number
            sele = surface_element_list(k)

            ! Check that there is only 1 BC applied per surface element
            if (vphi_bc_flag(sele) /= 0) then             
               FLExit('Cannot apply more than 1 BC to a surface element for InterstitialVelocityPorosity')
            end if

            ! Set the sele flag to indicate the BC type
            do j = 1, size(types)
               if (trim(types(j)) == trim(bctype)) exit
            end do
            if (j > size(types)) then
               FLAbort('Cannot find no_normal_flow or normal_flow bctype for InterstitialVelocityPorosity')
            end if
            
            vphi_bc_flag(sele) = j

            ! Set the normal_flow field values from this BC for its sele
            if (trim(bctype) == 'normal_flow') then
               call set(vphi_bc_value, &
                        ele_nodes(bc_surface_mesh, sele), &
                        ele_val(scalar_surface_field, k))
            end if
            
         end do BC_sele_loop

      end do BC_loop

      ewrite_minmax(vphi_bc_value)

      ewrite(1,*) 'Finished get InterstitialVelocityPorosity boundary condition data'

   end subroutine get_vphi_boundary_condition
  
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

   end subroutine darcy_impes_calculate_average_pressure

! ----------------------------------------------------------------------------

   subroutine darcy_impes_calculate_gradient_pressures(di)
      
      !!< Calculate the gradient pressures for each phase

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
      
   end subroutine darcy_impes_calculate_gradient_pressures

! ----------------------------------------------------------------------------

   subroutine darcy_impes_calculate_divergence_total_darcy_velocity(di)
      
      !!< Calculate the divergence of the total divergence velocity, 
      !!< which may include terms summed over subcycles

      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: vele, p, iloc, oloc, face, gi, ggi, sele, f_ele_counter, f_ele_base
      real    :: darcy_vel_face_value_dot_n
      real,    dimension(:,:),   allocatable :: grad_pressure_ele
      real,    dimension(:,:),   allocatable :: x_ele
      real,    dimension(:,:),   allocatable :: normal
      real,    dimension(:),     allocatable :: detwei
      real,    dimension(:),     allocatable :: normgi
      logical, dimension(:),     allocatable :: notvisited
      real,    dimension(:),     allocatable :: div_tvphi_rhs_local
      real,    dimension(:,:),   allocatable :: x_face_quad
      integer, dimension(:),     pointer     :: x_pmesh_nodes
      integer, dimension(:),     pointer     :: p_nodes      
      real,    dimension(:,:),   allocatable :: normal_bdy
      real,    dimension(:),     allocatable :: detwei_bdy
      real,    dimension(:,:),   allocatable :: x_ele_bdy
      real,    dimension(:),     allocatable :: div_tvphi_rhs_local_bdy
      integer, dimension(:),     allocatable :: p_nodes_bdy
      integer, parameter :: VPHI_BC_TYPE_NORMAL_FLOW = 1, VPHI_BC_TYPE_NO_NORMAL_FLOW = 2
                        
      ewrite(1,*) 'Calculate the DivergenceTotalDarcyVelocity'
      
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      allocate(normal(di%ndim,di%x_cvshape%ngi))
      allocate(detwei(di%x_cvshape%ngi))      
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(div_tvphi_rhs_local(ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(grad_pressure_ele(di%ndim,1))
      
      allocate(detwei_bdy(di%x_cvbdyshape%ngi))
      allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      allocate(div_tvphi_rhs_local_bdy(face_loc(di%pressure_mesh,1)))      
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(p_nodes_bdy(face_loc(di%pressure_mesh,1)))

      call zero(di%div_total_darcy_velocity)

      phase_loop: do p = 1, di%number_phase
                       
         vele_loop: do vele = 1, element_count(di%pressure_mesh)

            ! Initialise the face value counter for this ele
            f_ele_counter = 0
            
            ! Find the face value ele base
            f_ele_base = di%cached_phase_face_value_domain(p)%ele_base(vele)

            ! get the coordinate values for this element for each positions local node
            x_ele = ele_val(di%positions, vele)         

            ! get the coordinate values for this element for each quadrature point
            x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)

            ! The node indices of the pressure field
            p_nodes => ele_nodes(di%pressure_mesh, vele)

            ! The node indices of the positions projected to the pressure mesh
            x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)

            ! get the gradient pressure value for this element
            grad_pressure_ele = ele_val(di%gradient_pressure(p)%ptr, vele)         

            ! obtain the transformed determinant*weight and normals
            call transform_cvsurf_to_physical(x_ele, di%x_cvshape, &
                                              detwei, normal, di%cvfaces)

            ! Initialise array for the quadrature points of this 
            ! element for whether it has already been visited
            notvisited = .true.

            ! Initialise the local rhs to assemble for this element
            div_tvphi_rhs_local = 0.0

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
                      
                      f_ele_counter = f_ele_counter + 1
                      
                      ! correct the orientation of the normal so it points away from iloc
                      normgi = orientate_cvsurf_normgi(node_val(di%positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                      &x_face_quad(:,ggi), normal(:,ggi))

                      darcy_vel_face_value_dot_n = dot_product(grad_pressure_ele(:,1), normgi(:)) * &
                     &di%cached_phase_face_value_domain(p)%value((f_ele_base - 1) + f_ele_counter)
                      
                      div_tvphi_rhs_local(iloc) = div_tvphi_rhs_local(iloc) - &
                                                  darcy_vel_face_value_dot_n * &
                                                  detwei(ggi)

                      div_tvphi_rhs_local(oloc) = div_tvphi_rhs_local(oloc) + &
                                                  darcy_vel_face_value_dot_n * &
                                                  detwei(ggi)

                    end if check_visited

                  end do quadrature_loop

                end if is_neigh

              end do face_loop

            end do nodal_loop_i

            call addto(di%div_total_darcy_velocity, p_nodes, div_tvphi_rhs_local)

         end do vele_loop

         ! Get this phase vphi BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
         call get_vphi_boundary_condition(di%inter_velocity_porosity(p)%ptr, &
                                          (/"normal_flow   ", &
                                            "no_normal_flow"/), &
                                          di%bc_surface_mesh, &
                                          di%vphi_bc_value, &
                                          di%vphi_bc_flag)
         
         sele_loop: do sele = 1, surface_element_count(di%pressure_mesh)
            
            if (di%vphi_bc_flag(sele) == VPHI_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop
         
            ! Initialise the face value counter for this ele
            f_ele_counter = 0

            ! Find the face value ele base
            f_ele_base = di%cached_phase_face_value_boundary(p)%ele_base(sele)
            
            x_ele       = ele_val(di%positions, face_ele(di%positions, sele))
            x_ele_bdy   = face_val(di%positions, sele)
            p_nodes_bdy = face_global_nodes(di%pressure_mesh, sele)

            call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)

            ! Initialise the local rhs to assemble for this element
            div_tvphi_rhs_local_bdy = 0.0

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi
                       
                        f_ele_counter = f_ele_counter + 1
                        
                        if (di%vphi_bc_flag(sele) == VPHI_BC_TYPE_NORMAL_FLOW) then

                           div_tvphi_rhs_local_bdy(iloc) = div_tvphi_rhs_local_bdy(iloc) + &
                          &di%cached_phase_face_value_boundary(p)%value((f_ele_base - 1) + f_ele_counter) * detwei_bdy(ggi)                   
                        
                        else
                        
                        
                        
                        end if
                        
                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop

            call addto(di%div_total_darcy_velocity, p_nodes_bdy, div_tvphi_rhs_local_bdy)

         end do sele_loop
         
      end do phase_loop
      
      call scale(di%div_total_darcy_velocity, di%inverse_cv_mass_pressure_mesh)
      
      ewrite_minmax(di%div_total_darcy_velocity)
      
      ! deallocate local variables as required
      deallocate(x_ele)
      deallocate(normal)
      deallocate(detwei)
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(div_tvphi_rhs_local)      
      deallocate(x_face_quad)
      deallocate(grad_pressure_ele)

      deallocate(detwei_bdy)
      deallocate(normal_bdy)
      deallocate(div_tvphi_rhs_local_bdy)
      deallocate(x_ele_bdy)
      deallocate(p_nodes_bdy)
            
   end subroutine darcy_impes_calculate_divergence_total_darcy_velocity

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve_phase_saturations(di)
      
      !!< Assemble and solve the phase saturations
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p, number_subcycle
      real    :: dt_subcycle, max_cfl
      
      ! solve the saturation for each phase but not the first
      s_phase_loop: do p = 2, di%number_phase
         
         ! Deduce the number of subcycles to do and the subcycle time step size
                  
         if (di%subcy_opt_sat%have_number) then
         
            number_subcycle = di%subcy_opt_sat%number_advection_subcycle
            
            dt_subcycle = di%dt/real(number_subcycle)

            call set(di%cfl_subcycle, di%cfl(p)%ptr)

            call scale(di%cfl_subcycle, 1.0/real(number_subcycle))            

         else if (di%subcy_opt_sat%have_max_cfl) then
         
            ! Find the max cfl number for this phase (accounting for parallel)
            max_cfl = maxval(di%cfl(p)%ptr)

            call allmax(max_cfl)

            number_subcycle = max(1,ceiling(max_cfl/di%subcy_opt_sat%max_courant_per_advection_subcycle))

            dt_subcycle = di%dt/real(number_subcycle)

            call set(di%cfl_subcycle, di%cfl(p)%ptr)

            call scale(di%cfl_subcycle, 1.0/real(number_subcycle))

                  
         else 
            
            number_subcycle = 1
            
            dt_subcycle = di%dt

            call set(di%cfl_subcycle, di%cfl(p)%ptr)
         
         end if 
         
         ewrite(1,*) 'Solve phase ',p,' Saturation with number_subcycle ',number_subcycle,' and dt_subcycle ',dt_subcycle
         
         call solve_phase_saturation(di, &
                                     number_subcycle, &
                                     dt_subcycle, &
                                     p)
         
         ewrite(1,*) 'Finished solve phase ',p,' Saturation'
         
         ewrite_minmax(di%saturation(p)%ptr)
                  
      end do s_phase_loop
      
      ! Either solve the first phase saturation or calculate if diagnostic
      s1_solve_if: if (di%phase_one_saturation_diagnostic) then
      
         call darcy_impes_calculate_phase_one_saturation_diagnostic(di)
      
      else s1_solve_if

         ! Deduce the number of subcycles to do and the subcycle time step size
         
         if (di%subcy_opt_sat%have_number) then
         
            number_subcycle = di%subcy_opt_sat%number_advection_subcycle
            
            dt_subcycle = di%dt/real(number_subcycle)

            call set(di%cfl_subcycle, di%cfl(1)%ptr)

            call scale(di%cfl_subcycle, 1.0/real(number_subcycle))            

         else if (di%subcy_opt_sat%have_max_cfl) then
         
            ! Find the max cfl number for this phase (accounting for parallel)
            max_cfl = maxval(di%cfl(1)%ptr)

            call allmax(max_cfl)

            number_subcycle = max(1,ceiling(max_cfl/di%subcy_opt_sat%max_courant_per_advection_subcycle))

            dt_subcycle = di%dt/real(number_subcycle)

            call set(di%cfl_subcycle, di%cfl(1)%ptr)

            call scale(di%cfl_subcycle, 1.0/real(number_subcycle))
         
         else 

            number_subcycle = 1
            
            dt_subcycle = di%dt

            call set(di%cfl_subcycle, di%cfl(1)%ptr)
         
         end if 
         
         ewrite(1,*) 'Solve phase 1 Saturation with number_subcycle ',number_subcycle,' and dt_subcycle ',dt_subcycle

         call solve_phase_saturation(di, &
                                     number_subcycle, &
                                     dt_subcycle, &
                                     p = 1)

         ewrite(1,*) 'Finished solve phase 1 Saturation'
         
         ewrite_minmax(di%saturation(1)%ptr)
      
      end if s1_solve_if
      
   end subroutine darcy_impes_assemble_and_solve_phase_saturations

! ----------------------------------------------------------------------------

   subroutine solve_phase_saturation(di, &
                                     number_subcycle, &
                                     dt_subcycle, &
                                     p)
      
      !!< Assemble and solve the saturation for the given phase number p

      type(darcy_impes_type), intent(inout) :: di
      integer,                intent(in)    :: number_subcycle
      real,                   intent(in)    :: dt_subcycle
      integer,                intent(in)    :: p

      call solve_scalar_field_using_cv_cts_mesh(di%saturation(p)%ptr, &
                                                di%old_saturation(p)%ptr, &
                                                di%gradient_pressure(p)%ptr, &
                                                di%old_viscosity(p)%ptr, &
                                                di%old_relative_permeability(p)%ptr, &
                                                di%old_absolute_permeability, &
                                                di%positions, &
                                                di%positions_pressure_mesh, &
                                                di%old_sfield_subcycle, &
                                                di%cv_mass_pressure_mesh_with_porosity, &
                                                di%cv_mass_pressure_mesh_with_old_porosity, &
                                                di%old_sfield_upwind, &
                                                di%cfl(p)%ptr, &
                                                di%inter_velocity_porosity(p)%ptr, &
                                                di%bc_surface_mesh, &
                                                di%vphi_bc_value, &
                                                di%vphi_bc_flag, &
                                                di%rhs_adv, &
                                                di%rhs_time, &
                                                di%rhs, &
                                                di%lhs, &
                                                di%x_cvshape, &
                                                di%p_cvshape, &
                                                di%x_cvbdyshape, &
                                                di%cvfaces, &
                                                di%saturation_cv_options, &
                                                di%state, &
                                                di%ndim, &
                                                number_subcycle, &
                                                dt_subcycle, &
                                                di%cached_phase_face_value_domain(p), &
                                                di%cached_phase_face_value_boundary(p))
      
   end subroutine solve_phase_saturation

! ----------------------------------------------------------------------------

   subroutine solve_scalar_field_using_cv_cts_mesh(sfield, &
                                                   old_sfield, &
                                                   gradient_pressure, &
                                                   old_viscosity, &
                                                   old_relative_permeability, &
                                                   old_absolute_permeability, &
                                                   positions, &
                                                   positions_sfield_mesh, &
                                                   old_sfield_subcycle, &
                                                   cv_mass_sfield_mesh_with_porosity, &
                                                   cv_mass_sfield_mesh_with_old_porosity, &
                                                   old_sfield_upwind, &
                                                   cfl, &
                                                   vphi, &
                                                   bc_surface_mesh, &
                                                   vphi_bc_value, &
                                                   vphi_bc_flag, &
                                                   rhs_adv, &
                                                   rhs_time, &
                                                   rhs, &
                                                   lhs, &
                                                   x_cvshape, &
                                                   s_cvshape, &
                                                   x_cvbdyshape, &
                                                   cvfaces, &
                                                   sfield_cv_options, &
                                                   state, &
                                                   ndim, &
                                                   number_subcycle, &
                                                   dt_subcycle, &
                                                   cached_face_value_domain, &
                                                   cached_face_value_boundary)
      
      !!< Assemble and solve a time+advection equation for the scalar field 
      !!< using CV on a continous mesh with advecting velocity calculated 
      !!< on the go from the relation vphi = - Kk/nu (grad P)
      
      type(scalar_field),                  intent(inout) :: sfield
      type(scalar_field),                  intent(inout) :: old_sfield
      type(vector_field),                  intent(in)    :: gradient_pressure
      type(scalar_field),                  intent(in)    :: old_viscosity
      type(scalar_field),                  intent(in)    :: old_relative_permeability
      type(scalar_field),                  intent(in)    :: old_absolute_permeability
      type(vector_field),                  intent(in)    :: positions
      type(vector_field),                  intent(inout) :: positions_sfield_mesh
      type(scalar_field),                  intent(inout) :: old_sfield_subcycle
      type(scalar_field),                  intent(in)    :: cv_mass_sfield_mesh_with_porosity
      type(scalar_field),                  intent(in)    :: cv_mass_sfield_mesh_with_old_porosity
      type(csr_matrix),                    intent(inout) :: old_sfield_upwind
      type(scalar_field),                  intent(in)    :: cfl
      type(vector_field),                  intent(in)    :: vphi
      type(mesh_type),                     intent(in)    :: bc_surface_mesh
      type(scalar_field),                  intent(inout) :: vphi_bc_value
      integer,               dimension(:), intent(inout) :: vphi_bc_flag
      type(scalar_field),                  intent(inout) :: rhs_adv
      type(scalar_field),                  intent(inout) :: rhs_time
      type(scalar_field),                  intent(inout) :: rhs
      type(scalar_field),                  intent(inout) :: lhs
      type(element_type),                  intent(in)    :: x_cvshape
      type(element_type),                  intent(in)    :: s_cvshape
      type(element_type),                  intent(in)    :: x_cvbdyshape      
      type(cv_faces_type),                 intent(in)    :: cvfaces
      type(cv_options_type),               intent(in)    :: sfield_cv_options      
      type(state_type),      dimension(:), intent(inout) :: state
      integer,                             intent(in)    :: ndim
      integer,                             intent(in)    :: number_subcycle
      real,                                intent(in)    :: dt_subcycle
      type(cached_face_value_type),        intent(inout), optional :: cached_face_value_domain
      type(cached_face_value_type),        intent(inout), optional :: cached_face_value_boundary

      ! local variables
      integer :: vele, iloc, oloc, face, gi, ggi, sele, isub, f_ele_counter, f_ele_base, upwind_pos, dim
      real    :: income, face_value, vphi_dot_n, alpha_start, alpha_end
      real    :: old_sfield_face_value, old_relperm_face_value
      logical :: inflow, cached_face_value_present, determine_face_value
      real,    dimension(1)                  :: old_absperm_ele, old_visc_ele
      real,    dimension(:,:),   allocatable :: grad_pressure_ele      
      real,    dimension(:,:),   allocatable :: vphi_face
      real,    dimension(:),     allocatable :: old_sfield_subcycle_ele
      real,    dimension(:),     allocatable :: old_relperm_ele
      real,    dimension(:),     allocatable :: cfl_ele
      real,    dimension(:,:),   allocatable :: x_ele
      real,    dimension(:,:),   allocatable :: normal
      real,    dimension(:),     allocatable :: detwei
      real,    dimension(:),     allocatable :: normgi
      logical, dimension(:),     allocatable :: notvisited
      real,    dimension(:),     allocatable :: s_rhs_local
      real,    dimension(:,:),   allocatable :: x_face_quad
      integer, dimension(:),     pointer     :: x_smesh_nodes
      integer, dimension(:),     pointer     :: s_nodes      
      integer, dimension(:),     pointer     :: upwind_nodes
      real,    dimension(1)                  :: old_absperm_ele_bdy, old_visc_ele_bdy
      real,    dimension(:),     allocatable :: old_sfield_subcycle_ele_bdy
      real,    dimension(:),     allocatable :: ghost_old_sfield_subcycle_ele_bdy      
      real,    dimension(:,:),   allocatable :: normal_bdy
      real,    dimension(:),     allocatable :: detwei_bdy
      real,    dimension(:,:),   allocatable :: x_ele_bdy
      real,    dimension(:),     allocatable :: s_rhs_local_bdy
      integer, dimension(:),     allocatable :: s_nodes_bdy
      real,    dimension(:,:),   allocatable :: grad_pressure_ele_bdy
      real,    dimension(:),     allocatable :: old_relperm_ele_bdy
      real,    dimension(:),     allocatable :: bc_sele_val
      integer, dimension(:),     allocatable :: sfield_bc_type
      type(scalar_field) :: sfield_bc
      integer, parameter :: S_BC_TYPE_WEAKDIRICHLET  = 1, S_BC_TYPE_DIRICHLET         = 2, S_BC_TYPE_ZERO_FLUX = 3
      integer, parameter :: VPHI_BC_TYPE_NORMAL_FLOW = 1, VPHI_BC_TYPE_NO_NORMAL_FLOW = 2
             
      ! set a local flag for whether cached face value is present
      cached_face_value_present = present(cached_face_value_domain)
      
      if (cached_face_value_present) then
         
         if (.not. present(cached_face_value_boundary)) then
         
            FLAbort('Issue in solve_scalar_field_using_cv_cts_mesh, must have both domain and boundary cached face value present if using them')
         
         end if 
         
      end if
      
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(ndim, ele_loc(positions,1)))      
      allocate(normal(ndim,x_cvshape%ngi))
      allocate(detwei(x_cvshape%ngi))      
      allocate(normgi(ndim))
      allocate(notvisited(x_cvshape%ngi))
      allocate(s_rhs_local(ele_loc(sfield,1)))
      allocate(x_face_quad(ndim, x_cvshape%ngi))
      allocate(old_sfield_subcycle_ele(ele_loc(sfield,1)))
      allocate(cfl_ele(ele_loc(cfl,1)))
      allocate(old_relperm_ele(ele_loc(sfield,1)))
      allocate(vphi_face(ndim,s_cvshape%ngi))
      allocate(grad_pressure_ele(ndim,1))

      allocate(old_sfield_subcycle_ele_bdy(face_loc(sfield,1)))      
      allocate(ghost_old_sfield_subcycle_ele_bdy(face_loc(sfield,1)))
      allocate(detwei_bdy(x_cvbdyshape%ngi))
      allocate(normal_bdy(ndim, x_cvbdyshape%ngi))
      allocate(s_rhs_local_bdy(face_loc(sfield,1)))
      allocate(x_ele_bdy(ndim, face_loc(positions,1)))      
      allocate(grad_pressure_ele_bdy(ndim,1))
      allocate(s_nodes_bdy(face_loc(sfield,1)))
      allocate(old_relperm_ele_bdy(face_loc(sfield,1)))
      allocate(bc_sele_val(face_loc(sfield,1)))
            
      ! Solve the sfield for each subcycle via:            
      !  - Form the lhs = cv_mass_sfield_mesh_with_porosity / dt
      !  - Assemble the rhs_adv contribution via discretising using CV
      !  - and rhs = rhs_adv + rhs_time, where rhs_time = cv_mass_sfield_mesh_with_old_porosity * old_sfield / dt
      !  - solve for latest sfield and copy to start subcycle sfield step
      
      ! Note: the porosity at the start and end of a subcycle time step
      ! are linearly interpolated values from the main time step start and end
      
      call set(old_sfield_subcycle, old_sfield)

      
      ! Allocate and get the BC data. If the BC is time dependent then
      ! an end of overall time step values are used for all subcycles. 
      allocate(sfield_bc_type(surface_element_count(sfield)))         

      call get_entire_boundary_condition(sfield, &
                                         (/"weakdirichlet", &
                                           "dirichlet    ", &
                                           "zero_flux    "/), &
                                         sfield_bc, &
                                         sfield_bc_type)

      ! Get this phase vphi BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
      call get_vphi_boundary_condition(vphi, &
                                       (/"normal_flow   ", &
                                         "no_normal_flow"/), &
                                       bc_surface_mesh, &
                                       vphi_bc_value, &
                                       vphi_bc_flag)

      sub_loop: do isub = 1,number_subcycle
         
         ! form the start and end of subcycle dt
         ! porosity linear interpolents
         alpha_start = (isub - 1) / number_subcycle
         alpha_end   = isub / number_subcycle
         
         call set(lhs, cv_mass_sfield_mesh_with_porosity)
         
         call scale(lhs, alpha_end)
         
         call addto(lhs, cv_mass_sfield_mesh_with_old_porosity, scale = (1.0 - alpha_end))

         call scale(lhs, 1.0/dt_subcycle)
         
         ! form the rhs_time contribution
         call set(rhs_time, cv_mass_sfield_mesh_with_porosity)
         
         call scale(rhs_time, alpha_start)
         
         call addto(rhs_time, cv_mass_sfield_mesh_with_old_porosity, scale = (1.0 - alpha_start))
         
         call scale(rhs_time, 1.0/dt_subcycle)

         ! add rhs time term
         call set(rhs, rhs_time)         
         call scale(rhs, old_sfield_subcycle)
         
         ! assemble the rhs_adv contribution
         call assemble_rhs_adv()

         ! add rhs advection term
         call addto(rhs, rhs_adv)

         ! apply strong dirichlet BC
         call apply_dirichlet_conditions(lhs, rhs, sfield)

         ! Solve for the sfield
         sfield%val = rhs%val / lhs%val
         
         ! Copy latest solution to old subcycle
         call set(old_sfield_subcycle, sfield)
         
      end do sub_loop
           
      ! deallocate local variables as required
      deallocate(x_ele)
      deallocate(normal)
      deallocate(detwei)
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(s_rhs_local)
      deallocate(x_face_quad)
      deallocate(old_sfield_subcycle_ele)
      deallocate(cfl_ele)
      deallocate(old_relperm_ele)
      deallocate(vphi_face)
      deallocate(grad_pressure_ele)

      deallocate(old_sfield_subcycle_ele_bdy)
      deallocate(ghost_old_sfield_subcycle_ele_bdy)
      deallocate(detwei_bdy)
      deallocate(normal_bdy)
      deallocate(s_rhs_local_bdy)
      deallocate(x_ele_bdy)
      deallocate(grad_pressure_ele_bdy)      
      deallocate(s_nodes_bdy)
      deallocate(old_relperm_ele_bdy)
      deallocate(bc_sele_val)

      deallocate(sfield_bc_type)
      call deallocate(sfield_bc)
      
   contains
         
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         
         subroutine assemble_rhs_adv()
            
            !!< Assemble the rhs advection contribtion for sfield

            ! Inititalise rhs advection field
            call zero(rhs_adv)
            
            ! Decide if a face value needs to be determined
            determine_face_value = (cached_face_value_present .and. isub > 1) .or. &
                                  &(.not. cached_face_value_present)
            
            if (determine_face_value) then
            
               ! Determine the upwind saturation values if required for higher order CV face value
               if(need_upwind_values(sfield_cv_options)) then

                 call find_upwind_values(state, &
                                         positions_sfield_mesh, &
                                         old_sfield, &
                                         old_sfield_upwind, &
                                         old_sfield, &
                                         old_sfield_upwind, &
                                         option_path = trim(sfield%option_path))

               else

                 call zero(old_sfield_upwind)

               end if
            
            end if
      
            ! Initialise optimisation flag used in finding upwind value in high resolution schemes
            upwind_pos = 0
            
            ! Loop volume elements assembling local contributions    
            vol_element_loop: do vele = 1,element_count(sfield)

               ! Initialise the face value counter for this ele
               f_ele_counter = 0

               ! Find the face value ele base
               f_ele_base = cached_face_value_domain%ele_base(vele)

               ! get the gradient pressure value for this element
               grad_pressure_ele = ele_val(gradient_pressure, vele)         

               ! get the coordinate values for this element for each positions local node
               x_ele = ele_val(positions, vele)         

               ! get the coordinate values for this element for each quadrature point
               x_face_quad = ele_val_at_quad(positions, vele, x_cvshape)

               ! The node indices of the sfield field
               s_nodes => ele_nodes(sfield, vele)

               ! The node indices of the positions projected to the sfield mesh
               x_smesh_nodes => ele_nodes(positions_sfield_mesh, vele)

               if (determine_face_value) then

                  ! get the old_relperm ele values
                  old_relperm_ele = ele_val(old_relative_permeability, vele)

                  ! get the old_viscosity value for this element
                  old_visc_ele = ele_val(old_viscosity, vele)

                  ! get the old_absolute permeability value for this element
                  old_absperm_ele = ele_val(old_absolute_permeability, vele)         
               
                  ! get the old sfield ele values from start of subcycle
                  old_sfield_subcycle_ele = ele_val(old_sfield_subcycle, vele)

                  ! get the CFL values for this element
                  cfl_ele = ele_val(cfl, vele)
               
                  ! Determine the node numbers to use to determine the upwind values
                  if((sfield_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
                     (sfield_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

                     upwind_nodes => x_smesh_nodes

                  else

                     upwind_nodes => s_nodes

                  end if

                  ! the latest intersitial velocity_porosity at the quadrature points, used to determine upwind
                  do dim = 1,ndim

                     vphi_face(dim,:) =  - ele_val_at_quad(old_relative_permeability, vele, s_cvshape) * &
                                           old_absperm_ele(1) * grad_pressure_ele(dim,1) / old_visc_ele(1)     

                  end do
               
               end if
               
               ! obtain the transformed determinant*weight and normals
               call transform_cvsurf_to_physical(x_ele, x_cvshape, &
                                                 detwei, normal, cvfaces)

               ! Initialise array for the quadrature points of this 
               ! element for whether it has already been visited
               notvisited = .true.

               ! Initialise the local rhs to assemble for this element
               s_rhs_local = 0.0

               ! loop over local nodes within this element
               nodal_loop_i: do iloc = 1, sfield%mesh%shape%loc

                 ! loop over CV faces internal to this element
                 face_loop: do face = 1, cvfaces%faces

                   ! is this a face neighbouring iloc?
                   is_neigh: if(cvfaces%neiloc(iloc, face) /= 0) then

                     ! find the opposing local node across the CV face
                     oloc = cvfaces%neiloc(iloc, face)

                     ! loop over gauss points on face
                     quadrature_loop: do gi = 1, cvfaces%shape%ngi

                       ! global gauss pt index for this element
                       ggi = (face-1)*cvfaces%shape%ngi + gi

                       ! check if this quadrature point has already been visited
                       check_visited: if(notvisited(ggi)) then

                         notvisited(ggi) = .false.
                         
                         f_ele_counter = f_ele_counter + 1
                         
                         ! correct the orientation of the normal so it points away from iloc
                         normgi = orientate_cvsurf_normgi(node_val(positions_sfield_mesh, x_smesh_nodes(iloc)), &
                                                         &x_face_quad(:,ggi), normal(:,ggi))

                         if (determine_face_value) then

                            ! determine if the flow is in or out of the face at this quadrature
                            ! with respect to the normal orientation using latest vphi
                            vphi_dot_n = dot_product(vphi_face(:,ggi), normgi(:))

                            inflow = (vphi_dot_n<=0.0)

                            income = merge(1.0,0.0,inflow)

                            ! evaluate the nonlinear face value for sfield
                            call evaluate_face_val(old_sfield_face_value, &
                                                   old_sfield_face_value, & 
                                                   iloc, &
                                                   oloc, &
                                                   ggi, &
                                                   upwind_nodes, &
                                                   s_cvshape, &
                                                   old_sfield_subcycle_ele, &
                                                   old_sfield_subcycle_ele, &
                                                   old_sfield_upwind, &
                                                   old_sfield_upwind, &
                                                   inflow, &
                                                   cfl_ele, &
                                                   sfield_cv_options, &
                                                   save_pos = upwind_pos)

                            ! Evaluate the face value for old_relperm (assuming upwind for now) 
                            old_relperm_face_value = income*old_relperm_ele(oloc) + (1.0-income)*old_relperm_ele(iloc)

                            ! face value = S*old_relperm*old_absperm/old_visc, where old_absperm is phase independent
                            ! (if old_absperm and old_viscosity are considered tensors this requires modifying below)
                            face_value = old_sfield_face_value * old_relperm_face_value * old_absperm_ele(1) / old_visc_ele(1)

                            if (cached_face_value_present) then
                            
                               cached_face_value_domain%value((f_ele_base - 1) + f_ele_counter) = &
                              &cached_face_value_domain%value((f_ele_base - 1) + f_ele_counter) + face_value
                            
                            end if
                         
                         else 
                         
                            face_value = cached_face_value_domain%value((f_ele_base - 1) + f_ele_counter)
                         
                         end if

                         face_value = face_value * detwei(ggi) * dot_product(grad_pressure_ele(:,1), normgi)

                         ! Form the local rhs for iloc and opposing oloc with normal vector sign change
                         s_rhs_local(iloc) = s_rhs_local(iloc) + face_value

                         s_rhs_local(oloc) = s_rhs_local(oloc) - face_value

                       end if check_visited

                     end do quadrature_loop

                   end if is_neigh

                 end do face_loop

               end do nodal_loop_i

               ! Add volume element contribution to global rhs advection field
               call addto(rhs_adv, s_nodes, s_rhs_local)

            end do vol_element_loop

            ! Add BC integrals
            
            sele_loop: do sele = 1, surface_element_count(sfield)

               if (sfield_bc_type(sele) == S_BC_TYPE_DIRICHLET) cycle sele_loop

               if (sfield_bc_type(sele) == S_BC_TYPE_ZERO_FLUX) cycle sele_loop
               
               if (vphi_bc_flag(sele) == VPHI_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop

               ! Initialise the face value counter for this ele
               f_ele_counter = 0

               ! Find the face value ele base
               f_ele_base = cached_face_value_boundary%ele_base(sele)
                              
               if (vphi_bc_flag(sele) == VPHI_BC_TYPE_NORMAL_FLOW) then
               
                  if (determine_face_value) then

                     bc_sele_val = ele_val(vphi_bc_value, sele)

                  end if                
               
               else

                  old_relperm_ele_bdy   = face_val(old_relative_permeability, sele)
                  old_visc_ele_bdy      = face_val(old_viscosity, sele)
                  old_absperm_ele_bdy   = face_val(old_absolute_permeability, sele)
                  grad_pressure_ele_bdy = face_val(gradient_pressure, sele)
               
               end if

               if (sfield_bc_type(sele) == S_BC_TYPE_WEAKDIRICHLET) then
                  ghost_old_sfield_subcycle_ele_bdy = ele_val(sfield_bc, sele)
               else
                  ghost_old_sfield_subcycle_ele_bdy = face_val(old_sfield_subcycle, sele)               
               end if

               old_sfield_subcycle_ele_bdy = face_val(old_sfield_subcycle, sele)
                             
               x_ele       = ele_val(positions, face_ele(positions, sele))
               x_ele_bdy   = face_val(positions, sele)
               s_nodes_bdy = face_global_nodes(sfield%mesh, sele)

               call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, x_cvbdyshape, normal_bdy, detwei_bdy)

               ! Initialise the local rhs to assemble for this element
               s_rhs_local_bdy = 0.0

               bc_iloc_loop: do iloc = 1, sfield%mesh%faces%shape%loc

                  bc_face_loop: do face = 1, cvfaces%sfaces

                     bc_neigh_if: if(cvfaces%sneiloc(iloc,face)/=0) then

                        bc_quad_loop: do gi = 1, cvfaces%shape%ngi

                           ggi = (face-1)*cvfaces%shape%ngi + gi
                           
                           f_ele_counter = f_ele_counter + 1
                           
                           if (vphi_bc_flag(sele) == VPHI_BC_TYPE_NORMAL_FLOW) then
                           
                              if (determine_face_value) then

                                 ! Determine the upwind saturation boundary values to use from the sign 
                                 ! of the normal flow BC. If positive, hence outflow, always 
                                 ! use domain value. If negative, hence inflow, use the boundary 
                                 ! value which is the weak BC value if appropriate                                 
                                 if (bc_sele_val(iloc) > 0.0) then
                                 
                                    face_value = bc_sele_val(iloc) * old_sfield_subcycle_ele_bdy(iloc)
                                 
                                 else 
                                 
                                    face_value = bc_sele_val(iloc) * ghost_old_sfield_subcycle_ele_bdy(iloc)
                                 
                                 end if 
                                 
                                 if (cached_face_value_present) then

                                    cached_face_value_boundary%value((f_ele_base - 1) + f_ele_counter) = &
                                   &cached_face_value_boundary%value((f_ele_base - 1) + f_ele_counter) + face_value 

                                 end if 

                              else

                                 face_value = cached_face_value_boundary%value((f_ele_base - 1) + f_ele_counter)

                              end if 

                              s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - face_value * detwei_bdy(ggi)
                           
                           else                            
                                
                              vphi_dot_n = - old_relperm_ele_bdy(iloc) * old_absperm_ele_bdy(1) * &
                                             dot_product(grad_pressure_ele_bdy(:,1), normal_bdy(:,ggi)) / old_visc_ele_bdy(1)

                              inflow = (vphi_dot_n<=0.0)

                              income = merge(1.0,0.0,inflow)

                              old_sfield_face_value = income*ghost_old_sfield_subcycle_ele_bdy(iloc) + &
                                                      (1.0 - income)*old_sfield_subcycle_ele_bdy(iloc)

                              face_value = old_sfield_face_value * vphi_dot_n

                              s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - face_value * detwei_bdy(ggi)                           

                           end if
                           
                        end do bc_quad_loop

                     end if bc_neigh_if

                  end do bc_face_loop

               end do bc_iloc_loop

               call addto(rhs_adv, s_nodes_bdy, s_rhs_local_bdy)

            end do sele_loop
         
         end subroutine assemble_rhs_adv
      
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
   end subroutine solve_scalar_field_using_cv_cts_mesh

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
   
   subroutine darcy_impes_calculate_velocity_and_cfl_fields(di)
      
      !!< Calculate the various velocities, including fractional flow, and CFL fields
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      logical :: inflow
      integer :: vele, p, dim, iloc, oloc, face, gi, ggi, sele
      real    :: income, vphi_dot_n, vphi_face_value_dot_n, relperm_face_value
      real,    dimension(1)                  :: visc_ele, absperm_ele
      real,    dimension(:,:),   allocatable :: grad_pressure_ele
      real,    dimension(:),     allocatable :: relperm_ele
      real,    dimension(:),     allocatable :: inter_velocity_porosity_local
      real,    dimension(:,:),   allocatable :: vphi_face
      real,    dimension(:,:),   allocatable :: x_ele
      real,    dimension(:,:),   allocatable :: normal
      real,    dimension(:),     allocatable :: detwei
      real,    dimension(:),     allocatable :: normgi
      logical, dimension(:),     allocatable :: notvisited
      real,    dimension(:),     allocatable :: cfl_rhs_local
      real,    dimension(:,:),   allocatable :: x_face_quad
      integer, dimension(:),     pointer     :: x_pmesh_nodes
      integer, dimension(:),     pointer     :: p_nodes      
      real,    dimension(1)                  :: visc_ele_bdy, absperm_ele_bdy
      real,    dimension(:,:),   allocatable :: normal_bdy
      real,    dimension(:),     allocatable :: detwei_bdy
      real,    dimension(:,:),   allocatable :: x_ele_bdy
      real,    dimension(:),     allocatable :: cfl_rhs_local_bdy
      integer, dimension(:),     allocatable :: p_nodes_bdy
      real,    dimension(:,:),   allocatable :: grad_pressure_ele_bdy
      real,    dimension(:),     allocatable :: relperm_ele_bdy
      real,    dimension(:),     allocatable :: bc_sele_val
      integer, parameter :: VPHI_BC_TYPE_NORMAL_FLOW = 1, VPHI_BC_TYPE_NO_NORMAL_FLOW = 2
            
      ewrite(1,*) 'Calculate Velocity, FractionalFlow and CFL fields'
      
      ewrite(1,*) 'Calculate InterstitialVelocityPorosity'
      
      ! Calculate the inter_velocity_porosity field = - (1.0/sigma) * grad_pressure, where sigma = visc / absperm * relperm
      allocate(relperm_ele(ele_loc(di%pressure_mesh,1)))
      allocate(grad_pressure_ele(di%ndim,1))
      
      allocate(inter_velocity_porosity_local(ele_loc(di%inter_velocity_porosity(1)%ptr,1)))      
      
      do vele = 1, element_count(di%pressure_mesh)
         
         ! Find the element wise abs perm
         
         absperm_ele = ele_val(di%absolute_permeability, vele)
                  
         ! Form the inter_velocity_porosity
         do p = 1,di%number_phase

            ! Find the element wise gradient pressure 
            grad_pressure_ele = ele_val(di%gradient_pressure(p)%ptr, vele)
         
            ! Find the element wise viscosity
            visc_ele = ele_val(di%viscosity(p)%ptr, vele)
            
            ! Find the element wise local values for relperm
            relperm_ele = ele_val(di%relative_permeability(p)%ptr, vele)
            
            ! find the local DG values for inter_velocity_porosity for each geometric dimension
            ! (if absperm and viscosity are considered tensors this requires modifying below)
            do dim = 1,di%ndim
               
               inter_velocity_porosity_local(:) = - relperm_ele(:) * absperm_ele(1) * grad_pressure_ele(dim,1) / visc_ele(1)
            
               call set(di%inter_velocity_porosity(p)%ptr, &
                        dim, &
                        ele_nodes(di%inter_velocity_porosity(p)%ptr, vele), &
                        inter_velocity_porosity_local)
                           
            end do
            
         end do
         
      end do
            
      do p = 1,di%number_phase
         ewrite_minmax(di%inter_velocity_porosity(p)%ptr)
      end do
            
      deallocate(relperm_ele)
      deallocate(grad_pressure_ele)
      deallocate(inter_velocity_porosity_local)
      
      ! calculate the darcy velocity = inter_velocity_porosity * ssaturation
      do p = 1, di%number_phase         
         
         ewrite(1,*) 'Calculate DarcyVelocity for phase ',p
         
         call zero(di%darcy_velocity(p)%ptr)
         
         do vele = 1, element_count(di%pressure_mesh)
            
            do dim = 1, di%ndim
            
               call addto(di%darcy_velocity(p)%ptr, &
                          dim, &
                          ele_nodes(di%darcy_velocity(p)%ptr, vele), &
                          ele_val(di%saturation(p)%ptr, vele))
            
            end do
            
         end do 

         call scale(di%darcy_velocity(p)%ptr, di%inter_velocity_porosity(p)%ptr)
         
         ewrite_minmax(di%darcy_velocity(p)%ptr)
         
      end do 
      
      ewrite(1,*) 'Calculate TotalDarcyVelocity'
      
      ! calculate the total darcy velocity
      call set(di%total_darcy_velocity, di%darcy_velocity(1)%ptr)
      
      do p = 2, di%number_phase
         
         call addto(di%total_darcy_velocity, di%darcy_velocity(p)%ptr)
         
      end do
      
      ewrite_minmax(di%total_darcy_velocity)
            
      ! calculate the fractional flow for each phase
      do p = 1, di%number_phase
         
         ewrite(1,*) 'Calculate FractionalFlow for phase ',p
         
         call set(di%fractional_flow(p)%ptr, di%total_darcy_velocity)
         
         call invert(di%fractional_flow(p)%ptr)
         
         call scale(di%fractional_flow(p)%ptr, di%darcy_velocity(p)%ptr)
         
         ewrite_minmax(di%fractional_flow(p)%ptr)
      
      end do      

      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      allocate(normal(di%ndim,di%x_cvshape%ngi))
      allocate(detwei(di%x_cvshape%ngi))      
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(cfl_rhs_local(ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(vphi_face(di%ndim,di%p_cvshape%ngi))
      allocate(grad_pressure_ele(di%ndim,1))
      allocate(relperm_ele(ele_loc(di%pressure_mesh,1)))
      
      allocate(detwei_bdy(di%x_cvbdyshape%ngi))
      allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      allocate(cfl_rhs_local_bdy(face_loc(di%pressure_mesh,1)))
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(p_nodes_bdy(face_loc(di%pressure_mesh,1)))
      allocate(grad_pressure_ele_bdy(di%ndim,1))
      allocate(relperm_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(bc_sele_val(face_loc(di%pressure_mesh,1)))

      phase_loop: do p = 1, di%number_phase
         
         call zero(di%cfl(p)%ptr)
              
         vele_loop: do vele = 1, element_count(di%pressure_mesh)

            ! get the gradient pressure value for this element for this phase
            grad_pressure_ele = ele_val(di%gradient_pressure(p)%ptr, vele)         
            
            ! get the relperm ele values for this phase
            relperm_ele = ele_val(di%relative_permeability(p)%ptr, vele)

            ! get the viscosity value for this element for this phase
            visc_ele = ele_val(di%viscosity(p)%ptr, vele)

            ! get the absolute permeability value for this element
            absperm_ele = ele_val(di%absolute_permeability, vele)         
            
            ! get the coordinate values for this element for each positions local node
            x_ele = ele_val(di%positions, vele)         

            ! get the coordinate values for this element for each quadrature point
            x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)

            ! The node indices of the pressure field
            p_nodes => ele_nodes(di%pressure_mesh, vele)

            ! The node indices of the positions projected to the pressure mesh
            x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)

            ! the latest intersitial velocity_porosity at the quadrature points, used to determine upwind
            do dim = 1,di%ndim

               vphi_face(dim,:) = - ele_val_at_quad(di%relative_permeability(p)%ptr, vele, di%p_cvshape) * &
                                    absperm_ele(1) * grad_pressure_ele(dim,1) / visc_ele(1)     

            end do

            ! obtain the transformed determinant*weight and normals
            call transform_cvsurf_to_physical(x_ele, di%x_cvshape, &
                                              detwei, normal, di%cvfaces)

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
                                            
                      ! correct the orientation of the normal so it points away from iloc
                      normgi = orientate_cvsurf_normgi(node_val(di%positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                      &x_face_quad(:,ggi), normal(:,ggi))

                      ! determine if the flow is in or out of the face at this quadrature
                      ! with respect to the normal orientation using the latest vphi
                      vphi_dot_n = dot_product(vphi_face(:,ggi), normgi(:))

                      inflow = (vphi_dot_n<=0.0)

                      income = merge(1.0,0.0,inflow)

                      ! Evaluate the face value for relperm (assuming upwind for now) 
                      relperm_face_value = income*relperm_ele(oloc) + (1.0-income)*relperm_ele(iloc)

                      vphi_face_value_dot_n = - relperm_face_value * absperm_ele(1) * &
                                                dot_product(grad_pressure_ele(:,1), normgi) / visc_ele(1) 

                      cfl_rhs_local(iloc) = cfl_rhs_local(iloc) + &
                                            abs(vphi_face_value_dot_n) * &
                                            detwei(ggi) * &
                                            (1.0 - income)

                      cfl_rhs_local(oloc) = cfl_rhs_local(oloc) + &
                                            abs(vphi_face_value_dot_n) * &
                                            detwei(ggi) * &
                                            income

                    end if check_visited

                  end do quadrature_loop

                end if is_neigh

              end do face_loop

            end do nodal_loop_i

            call addto(di%cfl(p)%ptr, p_nodes, cfl_rhs_local)

         end do vele_loop

         ! Get this phase vphi BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
         call get_vphi_boundary_condition(di%inter_velocity_porosity(p)%ptr, &
                                          (/"normal_flow   ", &
                                            "no_normal_flow"/), &
                                          di%bc_surface_mesh, &
                                          di%vphi_bc_value, &
                                          di%vphi_bc_flag)
         
         sele_loop: do sele = 1, surface_element_count(di%pressure_mesh)
            
            if (di%vphi_bc_flag(sele) == VPHI_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop
            
            if (di%vphi_bc_flag(sele) == VPHI_BC_TYPE_NORMAL_FLOW) then
            
               bc_sele_val = ele_val(di%vphi_bc_value, sele)
            
            else
            
               relperm_ele_bdy       = face_val(di%relative_permeability(p)%ptr, sele)
               visc_ele_bdy          = face_val(di%viscosity(p)%ptr, sele)
               absperm_ele_bdy       = face_val(di%absolute_permeability, sele)
               grad_pressure_ele_bdy = face_val(di%gradient_pressure(p)%ptr, sele)
            
            end if
            
            x_ele       = ele_val(di%positions, face_ele(di%positions, sele))
            x_ele_bdy   = face_val(di%positions, sele)
            p_nodes_bdy = face_global_nodes(di%pressure_mesh, sele)

            call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)

            ! Initialise the local rhs to assemble for this element
            cfl_rhs_local_bdy = 0.0

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi
                        
                        if (di%vphi_bc_flag(sele) == VPHI_BC_TYPE_NORMAL_FLOW) then
                        
                           vphi_face_value_dot_n = bc_sele_val(iloc)
                           
                           if (vphi_face_value_dot_n > 0.0) then
                              income = 0.0
                           else
                              income = 1.0
                           end if
                           
                        else
                        
                           vphi_face_value_dot_n = - relperm_ele_bdy(iloc) * absperm_ele_bdy(1) * &
                                                     dot_product(grad_pressure_ele_bdy(:,1), normal_bdy(:,ggi)) / visc_ele_bdy(1)

                           inflow = (vphi_face_value_dot_n<=0.0)

                           income = merge(1.0,0.0,inflow)
                        
                        end if
                        
                        cfl_rhs_local_bdy(iloc) = cfl_rhs_local_bdy(iloc) + &
                                                 &abs(vphi_face_value_dot_n) * detwei_bdy(ggi) * (1.0 - income)                     

                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop

            call addto(di%cfl(p)%ptr, p_nodes_bdy, cfl_rhs_local_bdy)

         end do sele_loop
                  
         di%cfl(p)%ptr%val = di%cfl(p)%ptr%val * di%dt / di%cv_mass_pressure_mesh_with_porosity%val

         ewrite_minmax(di%cfl(p)%ptr)

      end do phase_loop      
      
      ! deallocate local variables as required
      deallocate(x_ele)
      deallocate(normal)
      deallocate(detwei)
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(cfl_rhs_local)
      deallocate(x_face_quad)
      deallocate(vphi_face)
      deallocate(grad_pressure_ele)
      deallocate(relperm_ele)

      deallocate(detwei_bdy)
      deallocate(normal_bdy)
      deallocate(x_ele_bdy)
      deallocate(p_nodes_bdy)
      deallocate(cfl_rhs_local_bdy)
      deallocate(grad_pressure_ele_bdy)      
      deallocate(relperm_ele_bdy)      
      deallocate(bc_sele_val)      
             
   end subroutine darcy_impes_calculate_velocity_and_cfl_fields

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
       
   end subroutine darcy_impes_calculate_cflnumber_field_based_dt

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_initialise_cached_phase_face_value(di)
   
      !!< Initialise the cached_phase_face_value variable's

      type(darcy_impes_type), intent(inout) :: di
      
      ! local variable
      integer :: p, face_value_counter, vele, sele, iloc, face, gi, ggi
      logical, dimension(:), pointer :: notvisited

      allocate(notvisited(di%x_cvshape%ngi))
      
      allocate(di%cached_phase_face_value_domain(di%number_phase))
      allocate(di%cached_phase_face_value_boundary(di%number_phase))
      
      phase_loop: do p = 1, di%number_phase
      
         ! initialise face value counter
         face_value_counter = 0

         allocate(di%cached_phase_face_value_domain(p)%ele_base(element_count(di%pressure_mesh)))

         vele_loop: do vele = 1, element_count(di%pressure_mesh)

            ! Set the ele_base value
            di%cached_phase_face_value_domain(p)%ele_base(vele) = face_value_counter + 1

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

                      face_value_counter = face_value_counter + 1

                    end if check_visited

                  end do quadrature_loop

                end if is_neigh

              end do face_loop

            end do nodal_loop_i

         end do vele_loop

         allocate(di%cached_phase_face_value_domain(p)%value(face_value_counter))

         ! initialise face value counter
         face_value_counter = 0

         allocate(di%cached_phase_face_value_boundary(p)%ele_base(surface_element_count(di%pressure_mesh)))

         sele_loop: do sele = 1, surface_element_count(di%pressure_mesh)

            ! Set the ele_base value
            di%cached_phase_face_value_boundary(p)%ele_base(sele) = face_value_counter + 1

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        face_value_counter = face_value_counter + 1

                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop

         end do sele_loop

         allocate(di%cached_phase_face_value_boundary(p)%value(face_value_counter))
      
      end do phase_loop

      deallocate(notvisited)
      
   end subroutine darcy_impes_initialise_cached_phase_face_value

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_inverse_characteristic_length(di)
   
      !!< Calculate the inverse_characteristic_length variable
      
      ! This is only set up for linear mesh
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variable
      integer :: vele, sele
      real    :: subcv_area, subcv_volume
      real, dimension(:), pointer :: len_vals
      
      allocate(len_vals(face_loc(di%pressure_mesh,1)))
      
      sele_loop: do sele = 1, surface_element_count(di%pressure_mesh)
         
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
      
   end subroutine darcy_impes_calculate_inverse_characteristic_length

! ----------------------------------------------------------------------------

end module darcy_impes_assemble_module
