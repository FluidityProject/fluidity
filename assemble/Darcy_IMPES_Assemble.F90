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
             darcy_impes_copy_to_old, &
             darcy_impes_assemble_and_solve, &
             darcy_impes_calculate_gradient_pressures, &
             darcy_impes_calculate_non_first_phase_pressures, &
             darcy_impes_calculate_phase_one_saturation_diagnostic, &
             darcy_impes_calculate_velocity_and_cfl_fields, &
             darcy_impes_calculate_sum_saturation, &
             darcy_impes_calculate_densities, &
             darcy_impes_calculate_cflnumber_field_based_dt, &
             darcy_impes_initialise_cached_phase_face_value, &
             darcy_impes_calculate_inverse_characteristic_length, &
             darcy_impes_calculate_modified_relative_permeability
   
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
      type(vector_field_pointer), dimension(:), pointer :: darcy_velocity_over_saturation
      type(vector_field_pointer), dimension(:), pointer :: old_darcy_velocity_over_saturation
      type(vector_field_pointer), dimension(:), pointer :: darcy_velocity   
      type(vector_field_pointer), dimension(:), pointer :: fractional_flow
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
      type(scalar_field), pointer :: sum_saturation
      type(scalar_field), pointer :: old_sum_saturation
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
      type(scalar_field) :: cfl_subcycle
      type(scalar_field) :: old_sfield_subcycle
      ! *** The modified relative permeability fields for each phase ***
      type(scalar_field_pointer), dimension(:), pointer :: modified_relative_permeability
      type(scalar_field_pointer), dimension(:), pointer :: old_modified_relative_permeability
      ! *** Minimum value of saturation used in denominator of modified relperm from options ***
      real :: min_sat_value_for_mod_relperm
      ! *** Data associated with v_over_s, pressure and saturation BC allocated here ***
      type(mesh_type)                          :: bc_surface_mesh
      type(scalar_field)                       :: v_over_s_bc_value
      integer,           dimension(:), pointer :: v_over_s_bc_flag
      type(scalar_field)                       :: pressure_bc_value
      integer,           dimension(:), pointer :: pressure_bc_flag
      type(scalar_field)                       :: saturation_bc_value
      integer,           dimension(:), pointer :: saturation_bc_flag
      type(scalar_field)                       :: inverse_characteristic_length
      real                                     :: weak_pressure_bc_coeff
      ! *** The number of phase and the CV surface quadrature degree to use ***
      integer :: number_phase, cv_surface_quaddegree   
      ! *** Data specifically associated with the CV discretisation allocated here ***
      type(cv_options_type)                   :: saturation_cv_options
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
      logical                                 :: phase_one_saturation_diagnostic
      type(darcy_impes_subcycle_options_type) :: subcy_opt_sat
      type(csr_matrix)                        :: old_sfield_upwind
      ! *** Flag for whether the first phase pressure is prognostic, else it is prescribed *** 
      logical :: first_phase_pressure_prognostic
      ! *** Flag for whether the Porosity is diagnostic, else it is prescribed *** 
      logical :: porosity_is_diagnostic
      ! *** Flag for whether the AbsolutePermeability is diagnostic, else it is prescribed *** 
      logical :: absolute_permeability_is_diagnostic
      ! *** Flag for whether the phase saturation sources are diagnostic, else it is prescribed ***
      logical, dimension(:), pointer :: saturation_source_is_diagnostic
      ! *** The cached phase face value at each quadrature point summed over subcycles if necessary for each phase ***
      type(cached_face_value_type), dimension(:), pointer :: cached_phase_face_value_domain
      type(cached_face_value_type), dimension(:), pointer :: cached_phase_face_value_boundary
      ! *** The subcycle timestep size and number of them ***
      real    :: dt_subcycle
      integer :: number_subcycle
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
      
      phase_loop: do p = 1, di%number_phase
      
         call set(di%old_modified_relative_permeability(p)%ptr, di%modified_relative_permeability(p)%ptr)
      
      end do phase_loop
      
   end subroutine darcy_impes_copy_to_old

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve(di)
      
      !!< Assemble and solve the Darcy equations using an CIMPESS algorithm
      !!< which is a modification of the IMPES to include consistent subcycling.
      
      type(darcy_impes_type), intent(inout) :: di
      
      ewrite(1,*) 'Start Darcy IMPES assemble and solve'
                  
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
      !!< due to the capilliary pressures of non first phases as well as gravity 
      !!< and the rate of change of porosity and the individual phase sources.
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      logical :: inflow, determine_face_value
      integer :: p, vele, sele, iloc, oloc, jloc, face, gi, ggi, f_ele_counter, f_ele_base, upwind_pos, dim
      real    :: income, face_value, iter_v_over_s_dot_n, modrelperm_absperm_over_visc, grad_cap_p_dot_n
      real    :: old_saturation_face_value, modrelperm_face_value, g_dot_n, den_face_value
      real,    dimension(1)                  :: absperm_ele, visc_ele
      real,    dimension(:),     allocatable :: old_saturation_ele
      real,    dimension(:),     allocatable :: modrelperm_ele
      real,    dimension(:),     allocatable :: den_ele
      real,    dimension(:,:),   allocatable :: grav_ele
      real,    dimension(:),     allocatable :: cfl_ele
      real,    dimension(:,:),   allocatable :: iter_grad_pressure_face_quad
      real,    dimension(:,:),   allocatable :: grad_cap_pressure_face_quad
      real,    dimension(:,:),   allocatable :: iter_v_over_s_face_quad
      real,    dimension(:,:),   allocatable :: x_ele
      real,    dimension(:,:,:), allocatable :: p_dshape
      real,    dimension(:,:),   allocatable :: normal
      real,    dimension(:),     allocatable :: detwei
      real,    dimension(:),     allocatable :: normgi
      logical, dimension(:),     allocatable :: notvisited
      real,    dimension(:,:),   allocatable :: p_mat_local
      real,    dimension(:),     allocatable :: p_rhs_local
      real,    dimension(:,:),   allocatable :: x_face_quad
      integer, dimension(:),     pointer     :: x_pmesh_nodes
      integer, dimension(:),     pointer     :: p_nodes      
      integer, dimension(:),     pointer     :: upwind_nodes
      real,    dimension(1)                  :: absperm_ele_bdy, visc_ele_bdy
      real,    dimension(:),     allocatable :: modrelperm_ele_bdy
      real,    dimension(:,:),   allocatable :: grav_ele_bdy
      real,    dimension(:,:),   allocatable :: iter_grad_pressure_face_quad_bdy
      real,    dimension(:,:),   allocatable :: iter_v_over_s_face_quad_bdy
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
      integer, parameter :: V_OVER_S_BC_TYPE_NORMAL_FLOW     = 1, V_OVER_S_BC_TYPE_NO_NORMAL_FLOW = 2
            
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      allocate(p_dshape(ele_loc(di%pressure_mesh,1), di%x_cvshape%ngi, mesh_dim(di%pressure_mesh)))
      allocate(normal(di%ndim,di%x_cvshape%ngi))
      allocate(detwei(di%x_cvshape%ngi))      
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(p_mat_local(ele_loc(di%pressure_mesh,1), ele_loc(di%pressure_mesh,1)))
      allocate(p_rhs_local(ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(old_saturation_ele(ele_loc(di%pressure_mesh,1)))
      allocate(modrelperm_ele(ele_loc(di%pressure_mesh,1)))
      allocate(den_ele(ele_loc(di%pressure_mesh,1)))
      allocate(grav_ele(di%ndim,1))
      allocate(cfl_ele(ele_loc(di%pressure_mesh,1)))
      allocate(iter_grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(grad_cap_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(iter_v_over_s_face_quad(di%ndim,di%p_cvshape%ngi))

      allocate(modrelperm_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(grav_ele_bdy(di%ndim,1))
      allocate(iter_grad_pressure_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      allocate(iter_v_over_s_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
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
            
      ! Decide if a face value needs to be determined
      determine_face_value = di%nonlinear_iter == 1

      ! Assemble a contribution from each phase to form a global continuity equation to solve for first phase pressure
      phase_loop: do p = 1, di%number_phase
         
         ewrite(1,*) 'Assemble volume contribution to global continuity from phase ',p
         
         if (determine_face_value) then

            ! Initialise the cached phase face values
            di%cached_phase_face_value_domain(p)%value   = 0.0
            di%cached_phase_face_value_boundary(p)%value = 0.0
         
            ! Determine the upwind saturation values if required for higher order CV face value
            if(need_upwind_values(di%saturation_cv_options)) then
              
              ! NOTE only old upwind values are required but this procedure 
              !      expects latest and old. So pass in old twice - not optimal
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

               ! get the old saturation ele values for this phase
               old_saturation_ele = ele_val(di%old_saturation(p)%ptr, vele)
            
               ! get the CFL values for this element
               cfl_ele = ele_val(di%cfl(p)%ptr, vele)
            
               ! Determine the node numbers to use to determine the Saturation upwind values
               if((di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
                  (di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

                  upwind_nodes => x_pmesh_nodes

               else

                  upwind_nodes => p_nodes

               end if
                        
            end if

            ! get the modrelperm ele values for this phase
            modrelperm_ele = ele_val(di%modified_relative_permeability(p)%ptr, vele)

            ! get the density value for this element for this phase
            den_ele = ele_val(di%density(p)%ptr, vele)
            
            ! The gravity values for this element for each direction
            grav_ele = ele_val(di%gravity, vele) 

            ! get the viscosity value for this element for this phase
            visc_ele = ele_val(di%viscosity(p)%ptr, vele)

            ! get the absolute permeability value for this element
            absperm_ele = ele_val(di%absolute_permeability, vele)         

            ! get the iterated gradient pressure at the cv surface quadrature points for each direction for this phase
            iter_grad_pressure_face_quad = ele_val_at_quad(di%iterated_gradient_pressure(p)%ptr, vele, di%gradp_cvshape)
            
            ! obtain the transformed determinant*weight and normals
            call transform_cvsurf_to_physical(x_ele, di%x_cvshape, detwei, normal, di%cvfaces)

            ! obtain the derivative of the pressure mesh shape function at the CV face quadrature points
            call transform_to_physical(di%positions, vele, x_shape = di%x_cvshape_full, &
                                       shape = di%p_cvshape_full, dshape = p_dshape)
            
            ! get the old gradient capilliary pressure at the cv surface quadrature points for each direction 
            if (p > 1) then
               grad_cap_pressure_face_quad = &
              &darcy_impes_ele_grad_at_quad_scalar(di%capilliary_pressure(p)%ptr, vele, dn = p_dshape)
            end if
            
            ! the iterated DarcyVelocityOverSaturation at the quadrature points
            ! determined from FE interpolation of each component, only used to determine upwind.
            do dim = 1,di%ndim

               iter_v_over_s_face_quad(dim,:) =  &
- (ele_val_at_quad(di%modified_relative_permeability(p)%ptr, vele, di%p_cvshape) * absperm_ele(1) / visc_ele(1)) * &
  (iter_grad_pressure_face_quad(dim,:) - ele_val_at_quad(di%density(p)%ptr, vele, di%p_cvshape) * grav_ele(dim,1))

            end do
            
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

                       f_ele_counter = f_ele_counter + 1

                       ! correct the orientation of the normal so it points away from iloc
                       normgi = orientate_cvsurf_normgi(node_val(di%positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                       &x_face_quad(:,ggi), normal(:,ggi))

                       ! determine if the flow is in or out of the face at this quadrature
                       ! with respect to the normal orientation using the iterated v_over_s
                       iter_v_over_s_dot_n = dot_product(iter_v_over_s_face_quad(:,ggi), normgi)

                       inflow = (iter_v_over_s_dot_n<=0.0)

                       income = merge(1.0,0.0,inflow)

                       if (determine_face_value) then

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

                          ! cache the phase face value for consistent use later
                          di%cached_phase_face_value_domain(p)%value((f_ele_base - 1) + f_ele_counter) = old_saturation_face_value

                       else 

                          old_saturation_face_value = di%cached_phase_face_value_domain(p)%value((f_ele_base - 1) + f_ele_counter)

                       end if 

                       ! Evaluate the face value for modrelperm (taking upwind) 
                       modrelperm_face_value = income*modrelperm_ele(oloc) + (1.0-income)*modrelperm_ele(iloc)
                       
                       ! Form the face value = detwei * (S*modrelperm*absperm/visc)
                       face_value = detwei(ggi)*old_saturation_face_value*modrelperm_face_value*absperm_ele(1)/visc_ele(1) 

                       ! Form the local matrix given by - n_i . sum_{phase} ( S*modrelperm*absperm/visc ) dP_1/dx_j
                       do jloc = 1,di%pressure_mesh%shape%loc

                          p_mat_local(iloc,jloc) = p_mat_local(iloc,jloc) - &
                                                   sum(p_dshape(jloc, ggi, :)*normgi, 1)*face_value

                          p_mat_local(oloc,jloc) = p_mat_local(oloc,jloc) - &
                                                   sum(p_dshape(jloc, ggi, :)*(-normgi), 1)*face_value

                       end do

                       ! Add gravity term to rhs = - n_i . sum_{phase} ( S*modrelperm*absperm/visc ) * den*grav

                       ! Find g dot n
                       g_dot_n = dot_product(grav_ele(:,1), normgi)

                       ! Find the density face value (taking upwind) 
                       den_face_value = income*den_ele(oloc) + (1.0-income)*den_ele(iloc)

                       p_rhs_local(iloc)  = p_rhs_local(iloc) - &
                                            face_value * &
                                            den_face_value * &
                                            g_dot_n

                       ! Add capilliary pressure term to rhs = n_i . sum_{phase} ( S*modrelperm*absperm/visc ) dP_c/dx_j
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
                  
      ! Add normal v_over_s * S BC integrals for each phase, else include if required 
      ! a weak pressure BC (which if not given is assumed zero).

      phase_loop_bc: do p = 1, di%number_phase

         ewrite(1,*) 'Assemble boundary contribution to global continuity from phase ',p
         
         ! Get the phase pressure BC - if no v_over_s and weak pressure then extra integrals are added
         call get_entire_saturation_or_pressure_boundary_condition(di%pressure(p)%ptr, &
                                                                   (/"weakdirichlet"/), &
                                                                   di%bc_surface_mesh, &
                                                                   di%pressure_bc_value, &
                                                                   di%pressure_bc_flag)
         
         ! Get this phase v_over_s BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
         call get_v_over_s_boundary_condition(di%darcy_velocity_over_saturation(p)%ptr, &
                                              (/"normal_flow   ", &
                                                "no_normal_flow"/), &
                                              di%bc_surface_mesh, &
                                              di%v_over_s_bc_value, &
                                              di%v_over_s_bc_flag)
         
         ! Get this phase satutation BC info - only weak dirichlet
         call get_entire_saturation_or_pressure_boundary_condition(di%saturation(p)%ptr, &
                                                                   (/"weakdirichlet"/), &
                                                                   di%bc_surface_mesh, &
                                                                   di%saturation_bc_value, &
                                                                   di%saturation_bc_flag)
         
         sele_loop: do sele = 1, surface_element_count(di%pressure_mesh)
            
            ! A no_normal_flow BC adds nothing to the matrix and rhs so cycle sele loop
            if (di%v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop
                     
            ! Initialise the face value counter for this ele
            f_ele_counter = 0

            ! Find the face value ele base
            f_ele_base = di%cached_phase_face_value_boundary(p)%ele_base(sele)
            
            if (di%v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NORMAL_FLOW) then
                              
               bc_sele_val = ele_val(di%v_over_s_bc_value, sele)
                        
            else
            
               ! get the modrelperm sele values for this phase
               modrelperm_ele_bdy = face_val(di%modified_relative_permeability(p)%ptr, sele)

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

               ! get the iterated gradient pressure at the cv surface quadrature points for each direction
               iter_grad_pressure_face_quad_bdy = face_val_at_quad(di%iterated_gradient_pressure(p)%ptr, sele, di%gradp_cvbdyshape)         

               ! The gravity values for this element for each direction
               grav_ele_bdy = face_val(di%gravity, sele) 

               ! the iterated DarcyVelocityOverSaturation at the quadrature points
               ! determined from FE interpolation of each component, 
               ! used to determine upwind and also added in certain integrals
               ! associated with weak pressure BCs.
               do dim = 1,di%ndim

                  iter_v_over_s_face_quad_bdy(dim,:) =  &
- (face_val_at_quad(di%modified_relative_permeability(p)%ptr, sele, di%p_cvbdyshape) * absperm_ele_bdy(1) / visc_ele_bdy(1)) * &
  (iter_grad_pressure_face_quad_bdy(dim,:) - face_val_at_quad(di%density(p)%ptr, sele, di%p_cvbdyshape) * grav_ele_bdy(dim,1))

               end do
               
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

            p_matrix_local_bdy = 0.0
            p_rhs_local_bdy    = 0.0

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi

                        f_ele_counter = f_ele_counter + 1
                                                
                        if (di%v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NORMAL_FLOW) then
                           
                           ! If have normal_flow BC for this phase then include CV integral of value * saturation
                           
                           if (determine_face_value) then

                              ! Determine the upwind saturation boundary values to use from the sign 
                              ! of the normal flow BC. If positive, hence outflow, always 
                              ! use domain value. If negative, hence inflow, use the boundary 
                              ! value which is the weak BC value if appropriate

                              if (bc_sele_val(iloc) > 0.0) then

                                 old_saturation_face_value = old_saturation_ele_bdy(iloc)
                              
                              else

                                 old_saturation_face_value = ghost_old_saturation_ele_bdy(iloc)
                              
                              end if 
                            
                              ! cache the phase face value for consistent use later
                              di%cached_phase_face_value_boundary(p)%value((f_ele_base - 1) + f_ele_counter) = old_saturation_face_value
                              
                           else
                              
                              old_saturation_face_value = di%cached_phase_face_value_boundary(p)%value((f_ele_base - 1) + f_ele_counter)
                           
                           end if
                           
                           p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - old_saturation_face_value * bc_sele_val(iloc) * detwei_bdy(ggi)
                        
                        else 
                           
                           ! Else apply a weak pressure BC, which has a value of zero by default (if not given)
                           ! - NOTE for non first phases P = P_1 + P_C, and P_C = 0.0 on boundary so
                           !        all phases pressure boundary conditions related to phase 1
                           
                           ! determine if the flow is in or out of the face at this quadrature
                           ! with respect to the normal orientation using the iterated v_over_s
                           iter_v_over_s_dot_n = dot_product(iter_v_over_s_face_quad_bdy(:,ggi), normal_bdy(:,ggi))
                           
                           ! Find modrelperm * absperm / visc
                           modrelperm_absperm_over_visc = modrelperm_ele_bdy(iloc) * absperm_ele_bdy(1) / visc_ele_bdy(1)
                           
                           ! Find upwind old saturation face value
                           if (determine_face_value) then

                              inflow = (iter_v_over_s_dot_n<=0.0)

                              income = merge(1.0,0.0,inflow)
                           
                              old_saturation_face_value = income*ghost_old_saturation_ele_bdy(iloc) + &
                                                         &(1.0-income)*old_saturation_ele_bdy(iloc)
                                                            
                              di%cached_phase_face_value_boundary(p)%value((f_ele_base - 1) + f_ele_counter) = old_saturation_face_value
                           
                           else
                           
                              old_saturation_face_value = di%cached_phase_face_value_boundary(p)%value((f_ele_base - 1) + f_ele_counter)
                           
                           end if 

                           ! Add coeff * gradient of pressure term to rhs via iter_v_over_s_dot_n * old_saturation_face_value
                           !  - which will include gravity and capilliary pressure term as required                         
                           ! (to be fully implicit part of this should be added to matrix but 2 non linear iterations may save us ...)
                           p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - &
                                                 & iter_v_over_s_dot_n * old_saturation_face_value * detwei_bdy(ggi)
if (.false.) then
                           ! If have phase pressure BC then include integral of BC value in rhs
                           ! - if not given then this implies adding a zero weak pressure BC
                           if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then
                                                            
                              do jloc = 1,di%pressure_mesh%faces%shape%loc 
                              
                                 p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - &
                                                         inv_char_len_ele_bdy(jloc) * &
                                                         di%p_cvbdyshape%n(jloc,ggi) * &
                                                         bc_sele_val(jloc) * &
                                                         modrelperm_absperm_over_visc * &
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
                                                              inv_char_len_ele_bdy(jloc) * &
                                                              di%p_cvbdyshape%n(jloc,ggi) * &
                                                              modrelperm_absperm_over_visc * &
                                                              old_saturation_face_value * &
                                                              sum(normal_bdy(:,ggi)) * &
                                                              detwei_bdy(ggi) * &
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
            
      ! Apply any strong dirichlet BC's that can only be applied to the first phase pressure
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
      deallocate(p_rhs_local)
      deallocate(x_face_quad)
      deallocate(old_saturation_ele)
      deallocate(modrelperm_ele)
      deallocate(den_ele)
      deallocate(grav_ele)
      deallocate(cfl_ele)
      deallocate(iter_grad_pressure_face_quad)
      deallocate(grad_cap_pressure_face_quad)
      deallocate(iter_v_over_s_face_quad)

      deallocate(modrelperm_ele_bdy)
      deallocate(grav_ele_bdy)
      deallocate(iter_grad_pressure_face_quad_bdy)
      deallocate(iter_v_over_s_face_quad_bdy)
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

   subroutine get_v_over_s_boundary_condition(v_over_s, &
                                              types, &
                                              bc_surface_mesh, &
                                              v_over_s_bc_value, &
                                              v_over_s_bc_flag)

      !!< Form the data associated with any v_over_s BC by returning 
      !!< full surface mesh arrays of a flag indicating  either 
      !!< no_normal_flow or normal_flow and the value.

      type(vector_field),               intent(in),   target :: v_over_s
      character(len=*),   dimension(:), intent(in)           :: types
      type(mesh_type),                  intent(in)           :: bc_surface_mesh
      type(scalar_field),               intent(inout)        :: v_over_s_bc_value
      integer,            dimension(:), intent(inout)        :: v_over_s_bc_flag

      ! Local variables
      character(len=FIELD_NAME_LEN)                       :: bctype
      type(scalar_field),                         pointer :: scalar_surface_field
      integer,                      dimension(:), pointer :: surface_element_list
      integer                                             :: i, j, k, sele

      ewrite(1,*) 'Get DarcyVelocityOverSaturation boundary condition data'

      ! Zero the normal_flow value surface field on whole boundary mesh  
      call zero(v_over_s_bc_value)

      ! Initialise flag for whether surface element has normal_flow BC
      v_over_s_bc_flag = 0

      ! Loop each BC object instance for the v_over_s
      ! (May have multiple normal_flow BC's applied to different surface id's)
      BC_loop: do i=1, get_boundary_condition_count(v_over_s)

         ! Get this BC info
         call get_boundary_condition(v_over_s, i, type = bctype, &
                                     surface_element_list = surface_element_list)

         ! check this is a normal_flow BC  or no_normal_flow (nothing else is permitted)
         if ((trim(bctype) /= 'normal_flow') .and. (trim(bctype) /= 'no_normal_flow')) then
            FLAbort('Have unknown BC type for a DarcyVelocityOverSaturation')
         end if

         ! Extract the scalar_surface_field for this BC for normal_flow
         if (trim(bctype) == 'normal_flow') then
            if (associated(v_over_s%bc%boundary_condition(i)%scalar_surface_fields)) then
               scalar_surface_field => v_over_s%bc%boundary_condition(i)%scalar_surface_fields(1)
            else
               FLAbort('Component scalar_surface_fields for DarcyVelocityOverSaturation BC type not associated')
            end if
         end if 
         
         ! Loop the surface elements associated with this BC instance
         ! and place the required BC value in whole boundary field
         ! for the normal_flow BC type
         BC_sele_loop: do k = 1, size(surface_element_list)

            ! Find the whole domain surface element number
            sele = surface_element_list(k)

            ! Check that there is only 1 BC applied per surface element
            if (v_over_s_bc_flag(sele) /= 0) then             
               FLExit('Cannot apply more than 1 BC to a surface element for DarcyVelocityOverSaturation')
            end if

            ! Set the sele flag to indicate the BC type
            do j = 1, size(types)
               if (trim(types(j)) == trim(bctype)) exit
            end do
            if (j > size(types)) then
               FLAbort('Cannot find no_normal_flow or normal_flow bctype for DarcyVelocityOverSaturation')
            end if
            
            v_over_s_bc_flag(sele) = j

            ! Set the normal_flow field values from this BC for its sele
            if (trim(bctype) == 'normal_flow') then
               call set(v_over_s_bc_value, &
                        ele_nodes(bc_surface_mesh, sele), &
                        ele_val(scalar_surface_field, k))
            end if
            
         end do BC_sele_loop

      end do BC_loop

      ewrite_minmax(v_over_s_bc_value)

      ewrite(1,*) 'Finished get DarcyVelocityOverSaturation boundary condition data'

   end subroutine get_v_over_s_boundary_condition
  
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
      logical :: inflow
      integer :: vele, p, iloc, oloc, jloc, face, gi, ggi, sele, f_ele_counter, f_ele_base, dim
      real    :: income, darcy_vel_face_value_dot_n, v_over_s_dot_n, iter_v_over_s_dot_n, den_face_value
      real    :: old_saturation_face_value, modrelperm_absperm_over_visc, modrelperm_face_value, face_value
      real,    dimension(1)                :: absperm_ele, visc_ele
      real,    dimension(:,:), allocatable :: grad_pressure_face_quad
      real,    dimension(:,:), allocatable :: iter_grad_pressure_face_quad
      real,    dimension(:),   allocatable :: modrelperm_ele
      real,    dimension(:),   allocatable :: den_ele
      real,    dimension(:,:), allocatable :: grav_ele
      real,    dimension(:,:), allocatable :: iter_v_over_s_face_quad
      real,    dimension(:,:), allocatable :: x_ele
      real,    dimension(:,:), allocatable :: normal
      real,    dimension(:),   allocatable :: detwei
      real,    dimension(:),   allocatable :: normgi
      logical, dimension(:),   allocatable :: notvisited
      real,    dimension(:),   allocatable :: div_tvel_rhs_local
      real,    dimension(:,:), allocatable :: x_face_quad
      integer, dimension(:),   pointer     :: x_pmesh_nodes
      integer, dimension(:),   pointer     :: p_nodes      
      real,    dimension(1)                :: absperm_ele_bdy, visc_ele_bdy
      real,    dimension(:),   allocatable :: modrelperm_ele_bdy
      real,    dimension(:,:), allocatable :: grav_ele_bdy
      real,    dimension(:,:), allocatable :: grad_pressure_face_quad_bdy
      real,    dimension(:,:), allocatable :: v_over_s_face_quad_bdy
      real,    dimension(:),   allocatable :: bc_sele_val
      real,    dimension(:),   allocatable :: pressure_ele_bdy
      real,    dimension(:),   allocatable :: inv_char_len_ele_bdy
      real,    dimension(:,:), allocatable :: normal_bdy
      real,    dimension(:),   allocatable :: detwei_bdy
      real,    dimension(:,:), allocatable :: x_ele_bdy
      real,    dimension(:),   allocatable :: div_tvel_rhs_local_bdy
      integer, dimension(:),   allocatable :: p_nodes_bdy
      integer, parameter :: V_OVER_S_BC_TYPE_NORMAL_FLOW   = 1, V_OVER_S_BC_TYPE_NO_NORMAL_FLOW = 2
      integer, parameter :: PRESSURE_BC_TYPE_WEAKDIRICHLET = 1
                        
      ewrite(1,*) 'Calculate the DivergenceTotalDarcyVelocity'
      
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      allocate(normal(di%ndim,di%x_cvshape%ngi))
      allocate(detwei(di%x_cvshape%ngi))      
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(div_tvel_rhs_local(ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(iter_grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(modrelperm_ele(ele_loc(di%pressure_mesh,1)))
      allocate(den_ele(ele_loc(di%pressure_mesh,1)))
      allocate(grav_ele(di%ndim,1))
      allocate(iter_v_over_s_face_quad(di%ndim,di%p_cvshape%ngi))
      
      allocate(modrelperm_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(grav_ele_bdy(di%ndim,1))
      allocate(grad_pressure_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      allocate(v_over_s_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      allocate(bc_sele_val(face_loc(di%pressure_mesh,1)))
      allocate(pressure_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(inv_char_len_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(detwei_bdy(di%x_cvbdyshape%ngi))
      allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      allocate(div_tvel_rhs_local_bdy(face_loc(di%pressure_mesh,1)))      
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

            ! get the modrelperm ele values for this phase
            modrelperm_ele = ele_val(di%modified_relative_permeability(p)%ptr, vele)

            ! get the density value for this element for this phase
            den_ele = ele_val(di%density(p)%ptr, vele)
            
            ! The gravity values for this element for each direction
            grav_ele = ele_val(di%gravity, vele) 

            ! get the viscosity value for this element for this phase
            visc_ele = ele_val(di%viscosity(p)%ptr, vele)

            ! get the absolute permeability value for this element
            absperm_ele = ele_val(di%absolute_permeability, vele)         

            ! get the latest gradient pressure at the cv surface quadrature points for each direction for this phase
            grad_pressure_face_quad = ele_val_at_quad(di%gradient_pressure(p)%ptr, vele, di%gradp_cvshape)

            ! get the iterated gradient pressure at the cv surface quadrature points for each direction for this phase
            ! - only used to determine upwind values, same as pressure assemble
            iter_grad_pressure_face_quad = ele_val_at_quad(di%iterated_gradient_pressure(p)%ptr, vele, di%gradp_cvshape)

            ! obtain the transformed determinant*weight and normals
            call transform_cvsurf_to_physical(x_ele, di%x_cvshape, &
                                              detwei, normal, di%cvfaces)

            ! the Iterated DarcyVelocityOverSaturation at the quadrature points
            ! determined from FE interpolation of each component, only used to determine upwind, same as pressure assemble
            do dim = 1,di%ndim

               iter_v_over_s_face_quad(dim,:) =  &
- (ele_val_at_quad(di%modified_relative_permeability(p)%ptr, vele, di%p_cvshape) * absperm_ele(1) / visc_ele(1)) * &
  (iter_grad_pressure_face_quad(dim,:) - ele_val_at_quad(di%density(p)%ptr, vele, di%p_cvshape) * grav_ele(dim,1))   

            end do

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

                       f_ele_counter = f_ele_counter + 1

                       ! correct the orientation of the normal so it points away from iloc
                       normgi = orientate_cvsurf_normgi(node_val(di%positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                       &x_face_quad(:,ggi), normal(:,ggi))

                       ! determine if the flow is in or out of the face at this quadrature
                       ! with respect to the normal orientation using the iterated v_over_s
                       iter_v_over_s_dot_n = dot_product(iter_v_over_s_face_quad(:,ggi), normgi)

                       inflow = (iter_v_over_s_dot_n<=0.0)

                       income = merge(1.0,0.0,inflow)

                       ! Find the density face value (taking upwind) 
                       den_face_value = income*den_ele(oloc) + (1.0-income)*den_ele(iloc)

                       ! Evaluate the face value for modrelperm (taking upwind) 
                       modrelperm_face_value = income*modrelperm_ele(oloc) + (1.0-income)*modrelperm_ele(iloc)
                       
                       old_saturation_face_value = di%cached_phase_face_value_domain(p)%value((f_ele_base - 1) + f_ele_counter)
                       
                       ! Form the face value = detwei * (S*modrelperm*absperm/visc)
                       face_value = detwei(ggi)*old_saturation_face_value*modrelperm_face_value*absperm_ele(1)/visc_ele(1)                       
                       
                       darcy_vel_face_value_dot_n = - dot_product( normgi, &
                      &face_value * (grad_pressure_face_quad(:,ggi) - den_face_value * grav_ele(:,1)) )

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

         ! Get this phase v_over_s BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
         call get_v_over_s_boundary_condition(di%darcy_velocity_over_saturation(p)%ptr, &
                                              (/"normal_flow   ", &
                                                "no_normal_flow"/), &
                                              di%bc_surface_mesh, &
                                              di%v_over_s_bc_value, &
                                              di%v_over_s_bc_flag)
         
         ! Get the pressure BC - required if no v_over_s given and weak pressure dirichlet given for extra integrals
         call get_entire_saturation_or_pressure_boundary_condition(di%pressure(p)%ptr, &
                                                                   (/"weakdirichlet"/), &
                                                                   di%bc_surface_mesh, &
                                                                   di%pressure_bc_value, &
                                                                   di%pressure_bc_flag)
         
         sele_loop: do sele = 1, surface_element_count(di%pressure_mesh)
            
            if (di%v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop
         
            ! Initialise the face value counter for this ele
            f_ele_counter = 0

            ! Find the face value ele base
            f_ele_base = di%cached_phase_face_value_boundary(p)%ele_base(sele)
            
            if (di%v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NORMAL_FLOW) then
              
               bc_sele_val = ele_val(di%v_over_s_bc_value, sele)
            
            else
            
               modrelperm_ele_bdy   = face_val(di%modified_relative_permeability(p)%ptr, sele)
               visc_ele_bdy         = face_val(di%viscosity(p)%ptr, sele)
               absperm_ele_bdy      = face_val(di%absolute_permeability, sele)
               inv_char_len_ele_bdy = ele_val(di%inverse_characteristic_length, sele)
               if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then
                   bc_sele_val = ele_val(di%pressure_bc_value, sele)     
               end if
               grad_pressure_face_quad_bdy = face_val_at_quad(di%gradient_pressure(p)%ptr, sele, di%gradp_cvbdyshape)         
               grav_ele_bdy                = face_val(di%gravity, sele) 
               pressure_ele_bdy            = face_val(di%pressure(p)%ptr, sele) 

               ! the latest DarcyVelocityOverSaturation at the quadrature points
               do dim = 1,di%ndim

                  v_over_s_face_quad_bdy(dim,:) =  &
- (face_val_at_quad(di%modified_relative_permeability(p)%ptr, sele, di%p_cvbdyshape) * absperm_ele_bdy(1) / visc_ele_bdy(1)) * &
(grad_pressure_face_quad_bdy(dim,:) - face_val_at_quad(di%density(p)%ptr, sele, di%p_cvbdyshape) * grav_ele_bdy(dim,1))   

               end do

            end if
            
            x_ele       = ele_val(di%positions, face_ele(di%positions, sele))
            x_ele_bdy   = face_val(di%positions, sele)
            p_nodes_bdy = face_global_nodes(di%pressure_mesh, sele)

            call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)

            ! Initialise the local rhs to assemble for this element
            div_tvel_rhs_local_bdy = 0.0

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi
                       
                        f_ele_counter = f_ele_counter + 1
                        
                        old_saturation_face_value = di%cached_phase_face_value_boundary(p)%value((f_ele_base - 1) + f_ele_counter)
                        
                        if (di%v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NORMAL_FLOW) then

                           div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) +&
                          &old_saturation_face_value * bc_sele_val(iloc) * detwei_bdy(ggi)                   
                        
                        else
                         
                           v_over_s_dot_n = dot_product(v_over_s_face_quad_bdy(:,ggi), normal_bdy(:,ggi))
                        
                           div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) + &
                          &old_saturation_face_value * v_over_s_dot_n * detwei_bdy(ggi)
if (.false.) then                           
                           modrelperm_absperm_over_visc = modrelperm_ele_bdy(iloc) * absperm_ele_bdy(1) / visc_ele_bdy(1)

                           if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then

                              do jloc = 1, di%pressure_mesh%faces%shape%loc 

                                 div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) - &
                                                                inv_char_len_ele_bdy(jloc) * &
                                                                di%p_cvbdyshape%n(jloc,ggi) * &
                                                                bc_sele_val(jloc) * &
                                                                modrelperm_absperm_over_visc * &
                                                                old_saturation_face_value * &
                                                                sum(normal_bdy(:,ggi)) * &
                                                                detwei_bdy(ggi) * &
                                                                di%weak_pressure_bc_coeff

                              end do

                           end if 

                           do jloc = 1, di%pressure_mesh%faces%shape%loc 

                              div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) + &
                                                             inv_char_len_ele_bdy(jloc) * &
                                                             di%p_cvbdyshape%n(jloc,ggi) * &
                                                             pressure_ele_bdy(jloc) * &
                                                             modrelperm_absperm_over_visc * &
                                                             old_saturation_face_value * &
                                                             sum(normal_bdy(:,ggi)) * &
                                                             detwei_bdy(ggi) * &
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
      deallocate(normal)
      deallocate(detwei)
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(div_tvel_rhs_local)      
      deallocate(x_face_quad)
      deallocate(grad_pressure_face_quad)
      deallocate(iter_grad_pressure_face_quad)
      deallocate(modrelperm_ele)
      deallocate(den_ele)
      deallocate(grav_ele)
      deallocate(iter_v_over_s_face_quad)

      deallocate(modrelperm_ele_bdy)
      deallocate(grav_ele_bdy)
      deallocate(grad_pressure_face_quad_bdy)
      deallocate(v_over_s_face_quad_bdy)
      deallocate(bc_sele_val)
      deallocate(pressure_ele_bdy)
      deallocate(inv_char_len_ele_bdy)
      deallocate(detwei_bdy)
      deallocate(normal_bdy)
      deallocate(div_tvel_rhs_local_bdy)
      deallocate(x_ele_bdy)
      deallocate(p_nodes_bdy)

      ewrite(1,*) 'Finished Calculate the DivergenceTotalDarcyVelocity'
            
   end subroutine darcy_impes_calculate_divergence_total_darcy_velocity

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve_phase_saturations(di)
      
      !!< Assemble and solve the phase saturations. Phase 1 is special
      !!< as it may be either solved prognostically or calculated 
      !!< diagnostically from the other phase saturation values that 
      !!< are solved first.
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p

      ewrite(1,*) 'Assemble and solve saturations'
         
      ! Deduce the number of subcycles to do and the subcycle time step size
      call deduce_subcycle_number_dt_and_cfl_field(di)
      
      ! solve the saturation for each phase but not the first (this is done after)
      s_phase_loop: do p = 2, di%number_phase
                  
         ewrite(1,*) 'Solve phase ',p,' Saturation with number_subcycle ',di%number_subcycle,' and dt_subcycle ',di%dt_subcycle
         
         call solve_phase_saturation(di, p)
         
         ewrite(1,*) 'Finished solve phase ',p,' Saturation'
         
         ewrite_minmax(di%saturation(p)%ptr)
                  
      end do s_phase_loop
      
      ! Either solve the first phase saturation or calculate if diagnostic
      s1_solve_if: if (di%phase_one_saturation_diagnostic) then
      
         call darcy_impes_calculate_phase_one_saturation_diagnostic(di)
      
      else s1_solve_if
         
         ewrite(1,*) 'Solve phase 1 Saturation with number_subcycle ',di%number_subcycle,' and dt_subcycle ',di%dt_subcycle

         call solve_phase_saturation(di, p = 1)

         ewrite(1,*) 'Finished solve phase 1 Saturation'
         
         ewrite_minmax(di%saturation(1)%ptr)
      
      end if s1_solve_if

      ewrite(1,*) 'Finished Assemble and solve saturations'
      
   end subroutine darcy_impes_assemble_and_solve_phase_saturations

! ----------------------------------------------------------------------------
   
   subroutine deduce_subcycle_number_dt_and_cfl_field(di)
      
      !!< Deduce the number of subcycles to do and the subcycle time step size

      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p
      real    :: max_cfl

      ewrite(1,*) 'Deduce subcycle number and dt to use for all phase saturations'
      
      if (di%subcy_opt_sat%have_number) then

         di%number_subcycle = di%subcy_opt_sat%number_advection_subcycle

         di%dt_subcycle = di%dt/real(di%number_subcycle)

      else if (di%subcy_opt_sat%have_max_cfl) then
         
         ! Initialise - as will take the max number of all phases
         di%number_subcycle = 1
         
         phase_loop: do p = 1, di%number_phase
         
            ! Find the max cfl number for this phase (accounting for parallel)
            max_cfl = maxval(di%cfl(p)%ptr)

            call allmax(max_cfl)

            di%number_subcycle = max(di%number_subcycle, & 
                                    &max(1, ceiling(max_cfl/di%subcy_opt_sat%max_courant_per_advection_subcycle)))

         end do phase_loop
         
         ! Find min subcycle timestep of all phases
         di%dt_subcycle = di%dt/real(di%number_subcycle)

      else 

         di%number_subcycle = 1

         di%dt_subcycle     = di%dt

      end if 

      ewrite(1,*) 'Finished Deduce subcycle number and dt to use for all phase saturations'
      
   end subroutine deduce_subcycle_number_dt_and_cfl_field

! ----------------------------------------------------------------------------

   subroutine solve_phase_saturation(di, p)
      
      !!< Assemble and solve the saturation for the given phase number p

      type(darcy_impes_type), intent(inout) :: di
      integer,                intent(in)    :: p
      
      ! Set the suitable subcycled CFL number - may be used for CV face limiting
      call set(di%cfl_subcycle, di%cfl(p)%ptr)

      call scale(di%cfl_subcycle, 1.0/real(di%number_subcycle))            

      call solve_scalar_field_using_cv_cts_mesh(di%saturation(p)%ptr, &
                                                di%old_saturation(p)%ptr, &
                                                di%gradient_pressure(p)%ptr, &
                                                di%iterated_gradient_pressure(p)%ptr, &
                                                di%viscosity(p)%ptr, &
                                                di%modified_relative_permeability(p)%ptr, &
                                                di%absolute_permeability, &
                                                di%density(p)%ptr, &
                                                di%gravity, &
                                                di%positions, &
                                                di%positions_pressure_mesh, &
                                                di%old_sfield_subcycle, &
                                                di%cv_mass_pressure_mesh_with_porosity, &
                                                di%cv_mass_pressure_mesh_with_old_porosity, &
                                                di%old_sfield_upwind, &
                                                di%cfl_subcycle, &
                                                di%darcy_velocity_over_saturation(p)%ptr, &
                                                di%bc_surface_mesh, &
                                                di%v_over_s_bc_value, &
                                                di%v_over_s_bc_flag, &
                                                di%pressure(p)%ptr, &
                                                di%pressure_bc_value, &
                                                di%pressure_bc_flag, &
                                                di%inverse_characteristic_length, &
                                                di%weak_pressure_bc_coeff, &
                                                di%saturation_source(p)%ptr, &
                                                di%cv_mass_pressure_mesh, &
                                                di%rhs_adv, &
                                                di%rhs_time, &
                                                di%rhs, &
                                                di%lhs, &
                                                di%x_cvshape, &
                                                di%p_cvshape, &
                                                di%gradp_cvshape, &
                                                di%x_cvbdyshape, &
                                                di%p_cvbdyshape, &
                                                di%gradp_cvbdyshape, &
                                                di%cvfaces, &
                                                di%saturation_cv_options, &
                                                di%state, &
                                                di%ndim, &
                                                di%number_subcycle, &
                                                di%dt_subcycle, &
                                                di%max_nonlinear_iter_this_timestep, &
                                                di%cached_phase_face_value_domain(p), &
                                                di%cached_phase_face_value_boundary(p))
      
   end subroutine solve_phase_saturation

! ----------------------------------------------------------------------------

   subroutine solve_scalar_field_using_cv_cts_mesh(sfield, &
                                                   old_sfield, &
                                                   gradient_pressure, &
                                                   iterated_gradient_pressure, &
                                                   viscosity, &
                                                   modified_relative_permeability, &
                                                   absolute_permeability, &
                                                   density, &
                                                   gravity, &
                                                   positions, &
                                                   positions_sfield_mesh, &
                                                   old_sfield_subcycle, &
                                                   cv_mass_sfield_mesh_with_porosity, &
                                                   cv_mass_sfield_mesh_with_old_porosity, &
                                                   old_sfield_upwind, &
                                                   cfl_subcycle, &
                                                   darcy_velocity_over_saturation, &
                                                   bc_surface_mesh, &
                                                   v_over_s_bc_value, &
                                                   v_over_s_bc_flag, &
                                                   pressure, &
                                                   pressure_bc_value, &
                                                   pressure_bc_flag, &
                                                   inverse_characteristic_length, &
                                                   weak_pressure_bc_coeff, &
                                                   s_source, &
                                                   cv_mass_sfield_mesh, &
                                                   rhs_adv, &
                                                   rhs_time, &
                                                   rhs, &
                                                   lhs, &
                                                   x_cvshape, &
                                                   p_cvshape, &
                                                   gradp_cvshape, &
                                                   x_cvbdyshape, &
                                                   p_cvbdyshape, &
                                                   gradp_cvbdyshape, &
                                                   cvfaces, &
                                                   sfield_cv_options, &
                                                   state, &
                                                   ndim, &
                                                   number_subcycle, &
                                                   dt_subcycle, &
                                                   max_nonlinear_iter_this_timestep, &
                                                   cached_face_value_domain, &
                                                   cached_face_value_boundary)
      
      !!< Assemble and solve a time+advection equation for the scalar field 
      !!< using CV on a continous mesh with advecting velocity calculated 
      !!< on the go from the relation v_over_s = - K*mod_k/nu (grad P - den*g)
      !!< as for higher than linear meshes this will avoid complicated basis 
      !!< (ie overlapping) expansions for the v_over_s.
      
      type(scalar_field),                  intent(inout) :: sfield
      type(scalar_field),                  intent(inout) :: old_sfield
      type(vector_field),                  intent(in)    :: gradient_pressure
      type(vector_field),                  intent(in)    :: iterated_gradient_pressure
      type(scalar_field),                  intent(in)    :: viscosity
      type(scalar_field),                  intent(in)    :: modified_relative_permeability
      type(scalar_field),                  intent(in)    :: absolute_permeability
      type(scalar_field),                  intent(in)    :: density
      type(vector_field),                  intent(in)    :: gravity
      type(vector_field),                  intent(in)    :: positions
      type(vector_field),                  intent(inout) :: positions_sfield_mesh
      type(scalar_field),                  intent(inout) :: old_sfield_subcycle
      type(scalar_field),                  intent(in)    :: cv_mass_sfield_mesh_with_porosity
      type(scalar_field),                  intent(in)    :: cv_mass_sfield_mesh_with_old_porosity
      type(csr_matrix),                    intent(inout) :: old_sfield_upwind
      type(scalar_field),                  intent(in)    :: cfl_subcycle
      type(vector_field),                  intent(in)    :: darcy_velocity_over_saturation
      type(mesh_type),                     intent(in)    :: bc_surface_mesh
      type(scalar_field),                  intent(inout) :: v_over_s_bc_value
      integer,               dimension(:), intent(inout) :: v_over_s_bc_flag
      type(scalar_field),                  intent(in)    :: pressure
      type(scalar_field),                  intent(inout) :: pressure_bc_value
      integer,               dimension(:), intent(inout) :: pressure_bc_flag
      type(scalar_field),                  intent(in)    :: inverse_characteristic_length
      real,                                intent(in)    :: weak_pressure_bc_coeff
      type(scalar_field),                  intent(in)    :: s_source
      type(scalar_field),                  intent(in)    :: cv_mass_sfield_mesh
      type(scalar_field),                  intent(inout) :: rhs_adv
      type(scalar_field),                  intent(inout) :: rhs_time
      type(scalar_field),                  intent(inout) :: rhs
      type(scalar_field),                  intent(inout) :: lhs
      type(element_type),                  intent(in)    :: x_cvshape
      type(element_type),                  intent(in)    :: p_cvshape
      type(element_type),                  intent(in)    :: gradp_cvshape
      type(element_type),                  intent(in)    :: x_cvbdyshape
      type(element_type),                  intent(in)    :: p_cvbdyshape
      type(element_type),                  intent(in)    :: gradp_cvbdyshape
      type(cv_faces_type),                 intent(in)    :: cvfaces
      type(cv_options_type),               intent(in)    :: sfield_cv_options      
      type(state_type),      dimension(:), intent(inout) :: state
      integer,                             intent(in)    :: ndim
      integer,                             intent(in)    :: number_subcycle
      real,                                intent(in)    :: dt_subcycle
      integer,                             intent(in)    :: max_nonlinear_iter_this_timestep
      type(cached_face_value_type),        intent(inout), optional :: cached_face_value_domain
      type(cached_face_value_type),        intent(inout), optional :: cached_face_value_boundary

      ! local variables
      logical :: inflow, cached_face_value_present, determine_face_value
      integer :: vele, iloc, oloc, jloc, face, gi, ggi, sele, isub, f_ele_counter, f_ele_base, upwind_pos, dim
      real    :: income, face_value, v_over_s_dot_n, iter_v_over_s_dot_n
      real    :: alpha_start, alpha_end, modrelperm_absperm_over_visc
      real    :: old_sfield_face_value, modrelperm_face_value, den_face_value
      real,    dimension(1)                :: absperm_ele, visc_ele
      real,    dimension(:,:), allocatable :: grad_pressure_face_quad
      real,    dimension(:,:), allocatable :: iter_grad_pressure_face_quad
      real,    dimension(:),   allocatable :: den_ele
      real,    dimension(:,:), allocatable :: grav_ele
      real,    dimension(:,:), allocatable :: v_over_s_face_quad
      real,    dimension(:,:), allocatable :: iter_v_over_s_face_quad
      real,    dimension(:),   allocatable :: old_sfield_subcycle_ele
      real,    dimension(:),   allocatable :: modrelperm_ele
      real,    dimension(:),   allocatable :: cfl_subcycle_ele
      real,    dimension(:,:), allocatable :: x_ele
      real,    dimension(:,:), allocatable :: normal
      real,    dimension(:),   allocatable :: detwei
      real,    dimension(:),   allocatable :: normgi
      logical, dimension(:),   allocatable :: notvisited
      real,    dimension(:),   allocatable :: s_rhs_local
      real,    dimension(:,:), allocatable :: x_face_quad
      integer, dimension(:),   pointer     :: x_smesh_nodes
      integer, dimension(:),   pointer     :: s_nodes      
      integer, dimension(:),   pointer     :: upwind_nodes
      real,    dimension(1)                :: absperm_ele_bdy, visc_ele_bdy
      real,    dimension(:),   allocatable :: modrelperm_ele_bdy
      real,    dimension(:,:), allocatable :: grav_ele_bdy
      real,    dimension(:,:), allocatable :: grad_pressure_face_quad_bdy
      real,    dimension(:,:), allocatable :: iter_grad_pressure_face_quad_bdy
      real,    dimension(:,:), allocatable :: v_over_s_face_quad_bdy
      real,    dimension(:,:), allocatable :: iter_v_over_s_face_quad_bdy
      real,    dimension(:),   allocatable :: bc_sele_val
      real,    dimension(:),   allocatable :: pressure_ele_bdy
      real,    dimension(:),   allocatable :: inv_char_len_ele_bdy
      real,    dimension(:),   allocatable :: old_sfield_subcycle_ele_bdy
      real,    dimension(:),   allocatable :: ghost_old_sfield_subcycle_ele_bdy      
      real,    dimension(:,:), allocatable :: normal_bdy
      real,    dimension(:),   allocatable :: detwei_bdy
      real,    dimension(:,:), allocatable :: x_ele_bdy
      real,    dimension(:),   allocatable :: s_rhs_local_bdy
      integer, dimension(:),   allocatable :: s_nodes_bdy
      integer, dimension(:),   allocatable :: sfield_bc_flag
      type(scalar_field) :: sfield_bc
      integer, parameter :: S_BC_TYPE_WEAKDIRICHLET        = 1, S_BC_TYPE_DIRICHLET             = 2
      integer, parameter :: V_OVER_S_BC_TYPE_NORMAL_FLOW   = 1, V_OVER_S_BC_TYPE_NO_NORMAL_FLOW = 2
      integer, parameter :: PRESSURE_BC_TYPE_WEAKDIRICHLET = 1
             
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
      allocate(cfl_subcycle_ele(ele_loc(cfl_subcycle,1)))
      allocate(modrelperm_ele(ele_loc(sfield,1)))
      allocate(grad_pressure_face_quad(ndim, p_cvshape%ngi))
      allocate(iter_grad_pressure_face_quad(ndim, p_cvshape%ngi))
      allocate(den_ele(ele_loc(sfield,1)))
      allocate(grav_ele(ndim,1))
      allocate(v_over_s_face_quad(ndim,p_cvshape%ngi))
      allocate(iter_v_over_s_face_quad(ndim,p_cvshape%ngi))

      allocate(modrelperm_ele_bdy(face_loc(sfield,1)))
      allocate(grav_ele_bdy(ndim,1))
      allocate(grad_pressure_face_quad_bdy(ndim, p_cvbdyshape%ngi))
      allocate(iter_grad_pressure_face_quad_bdy(ndim, p_cvbdyshape%ngi))
      allocate(v_over_s_face_quad_bdy(ndim, p_cvbdyshape%ngi))
      allocate(iter_v_over_s_face_quad_bdy(ndim, p_cvbdyshape%ngi))
      allocate(bc_sele_val(face_loc(sfield,1)))
      allocate(pressure_ele_bdy(face_loc(sfield,1)))
      allocate(inv_char_len_ele_bdy(face_loc(sfield,1)))
      allocate(old_sfield_subcycle_ele_bdy(face_loc(sfield,1)))      
      allocate(ghost_old_sfield_subcycle_ele_bdy(face_loc(sfield,1)))
      allocate(detwei_bdy(x_cvbdyshape%ngi))
      allocate(normal_bdy(ndim, x_cvbdyshape%ngi))
      allocate(s_rhs_local_bdy(face_loc(sfield,1)))
      allocate(x_ele_bdy(ndim, face_loc(positions,1)))      
      allocate(s_nodes_bdy(face_loc(sfield,1)))
            
      ! Solve the sfield for each subcycle via:            
      !  - Form the lhs = cv_mass_sfield_mesh_with_porosity / dt
      !  - Add the s_source to rhs
      !  - Form the rhs_time = cv_mass_sfield_mesh_with_old_porosity * old_sfield / dt and add to rhs
      !  - Assemble the rhs_adv contribution and add to rhs 
      !  - Apply any strong dirichlet BCs
      !  - solve for latest sfield and copy to start subcycle sfield step
      
      ! Note: the porosity at the start and end of a subcycle time step
      ! are linearly interpolated values from the main time step start and end
      
      call set(old_sfield_subcycle, old_sfield)
      
      ! Allocate and get the BC data. If the BC is time dependent then
      ! the end of overall time step values are used for all subcycles. 
      
      ! Get the sfield BC info
      allocate(sfield_bc_flag(surface_element_count(sfield)))
      call get_entire_boundary_condition(sfield, &
                                         (/"weakdirichlet", &
                                           "dirichlet    "/), &
                                         sfield_bc, &
                                         sfield_bc_flag)

      ! Get this phase v_over_s BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
      call get_v_over_s_boundary_condition(darcy_velocity_over_saturation, &
                                           (/"normal_flow   ", &
                                             "no_normal_flow"/), &
                                           bc_surface_mesh, &
                                           v_over_s_bc_value, &
                                           v_over_s_bc_flag)
      
      ! Get the pressure BC - required if no v_over_s given and weak pressure dirichlet given for extra integrals
      call get_entire_saturation_or_pressure_boundary_condition(pressure, &
                                                                (/"weakdirichlet"/), &
                                                                bc_surface_mesh, &
                                                                pressure_bc_value, &
                                                                pressure_bc_flag)

      sub_loop: do isub = 1,number_subcycle
         
         ! form the start and end of subcycle dt
         ! porosity linear interpolents
         alpha_start = (isub - 1) / number_subcycle
         alpha_end   = isub / number_subcycle
         
         call set(lhs, cv_mass_sfield_mesh_with_porosity)
         
         call scale(lhs, alpha_end)
         
         call addto(lhs, cv_mass_sfield_mesh_with_old_porosity, scale = (1.0 - alpha_end))

         call scale(lhs, 1.0/dt_subcycle)
         
         ! Add the s_source contribution
         call set(rhs, s_source)
         
         call scale(rhs, cv_mass_sfield_mesh)
         
         ! form the rhs_time contribution and add
         call set(rhs_time, cv_mass_sfield_mesh_with_porosity)
         
         call scale(rhs_time, alpha_start)
         
         call addto(rhs_time, cv_mass_sfield_mesh_with_old_porosity, scale = (1.0 - alpha_start))
         
         call scale(rhs_time, 1.0/dt_subcycle)

         call scale(rhs_time, old_sfield_subcycle)

         call addto(rhs, rhs_time)
         
         ! assemble the rhs_adv contribution and add
         call assemble_rhs_adv()

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
      deallocate(cfl_subcycle_ele)
      deallocate(modrelperm_ele)
      deallocate(grad_pressure_face_quad)
      deallocate(iter_grad_pressure_face_quad)
      deallocate(den_ele)
      deallocate(grav_ele)
      deallocate(v_over_s_face_quad)
      deallocate(iter_v_over_s_face_quad)

      deallocate(modrelperm_ele_bdy)
      deallocate(grav_ele_bdy)
      deallocate(grad_pressure_face_quad_bdy)
      deallocate(iter_grad_pressure_face_quad_bdy)
      deallocate(v_over_s_face_quad_bdy)
      deallocate(iter_v_over_s_face_quad_bdy)
      deallocate(bc_sele_val)
      deallocate(pressure_ele_bdy)
      deallocate(inv_char_len_ele_bdy)
      deallocate(old_sfield_subcycle_ele_bdy)      
      deallocate(ghost_old_sfield_subcycle_ele_bdy)
      deallocate(detwei_bdy)
      deallocate(normal_bdy)
      deallocate(s_rhs_local_bdy)
      deallocate(x_ele_bdy)
      deallocate(s_nodes_bdy)

      deallocate(sfield_bc_flag)
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
                                         old_sfield_subcycle, &
                                         old_sfield_upwind, &
                                         old_sfield_subcycle, &
                                         old_sfield_upwind, &
                                         option_path = trim(sfield%option_path))

               else

                 call zero(old_sfield_upwind)

               end if
            
            end if
      
            ! Initialise optimisation flag used in finding upwind value in high resolution schemes
            upwind_pos = 0
            
            ! Loop volume elements assembling local contributions    
            vol_element_loop: do vele = 1, element_count(sfield)

               ! Initialise the face value counter for this ele
               f_ele_counter = 0

               ! Find the face value ele base
               f_ele_base = cached_face_value_domain%ele_base(vele)

               ! get the coordinate values for this element for each positions local node
               x_ele = ele_val(positions, vele)         

               ! get the coordinate values for this element for each quadrature point
               x_face_quad = ele_val_at_quad(positions, vele, x_cvshape)

               ! The node indices of the sfield field
               s_nodes => ele_nodes(sfield, vele)

               ! The node indices of the positions projected to the sfield mesh
               x_smesh_nodes => ele_nodes(positions_sfield_mesh, vele)

               ! get the modrelperm ele values
               modrelperm_ele = ele_val(modified_relative_permeability, vele)

               ! get the density value for this element
               den_ele = ele_val(density, vele)

               ! The gravity values for this element for each direction
               grav_ele = ele_val(gravity, vele) 

               ! get the latest gradient pressure at the cv surface quadrature points for each direction 
               grad_pressure_face_quad = ele_val_at_quad(gradient_pressure, vele, gradp_cvshape)

               ! get the iterated gradient pressure at the cv surface quadrature points for each direction 
               iter_grad_pressure_face_quad = ele_val_at_quad(iterated_gradient_pressure, vele, gradp_cvshape)

               ! get the viscosity value for this element
               visc_ele = ele_val(viscosity, vele)

               ! get the absolute permeability value for this element
               absperm_ele = ele_val(absolute_permeability, vele)         
               
               ! Get necessary element data for determining sfield CV face value
               if (determine_face_value) then
               
                  ! get the old sfield ele values from start of subcycle
                  old_sfield_subcycle_ele = ele_val(old_sfield_subcycle, vele)

                  ! get the CFL values for this element
                  cfl_subcycle_ele = ele_val(cfl_subcycle, vele)
               
                  ! Determine the node numbers to use to determine the upwind values
                  if((sfield_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
                     (sfield_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

                     upwind_nodes => x_smesh_nodes

                  else

                     upwind_nodes => s_nodes

                  end if

               end if

               ! the iterated DarcyVelocityOverSaturation at the quadrature points
               ! determined from FE interpolation of each component, only used to determine upwind.
               do dim = 1,ndim

                  iter_v_over_s_face_quad(dim,:) =  &
- (ele_val_at_quad(modified_relative_permeability, vele, p_cvshape) * absperm_ele(1) / visc_ele(1)) * &
  (iter_grad_pressure_face_quad(dim,:) - ele_val_at_quad(density, vele, p_cvshape) * grav_ele(dim,1))   

               end do
                              
               ! obtain the transformed determinant*weight and normals
               call transform_cvsurf_to_physical(x_ele, x_cvshape, detwei, normal, cvfaces)

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

                         ! determine if the flow is in or out of the face at this quadrature
                         ! with respect to the normal orientation using iterated v_over_s
                         iter_v_over_s_dot_n = dot_product(iter_v_over_s_face_quad(:,ggi), normgi(:))

                         inflow = (iter_v_over_s_dot_n<=0.0)

                         income = merge(1.0,0.0,inflow)
                         
                         ! Determine the sfield CV face value
                         if (determine_face_value) then

                            ! evaluate the nonlinear face value for sfield
                            call evaluate_face_val(old_sfield_face_value, &
                                                   old_sfield_face_value, & 
                                                   iloc, &
                                                   oloc, &
                                                   ggi, &
                                                   upwind_nodes, &
                                                   p_cvshape, &
                                                   old_sfield_subcycle_ele, &
                                                   old_sfield_subcycle_ele, &
                                                   old_sfield_upwind, &
                                                   old_sfield_upwind, &
                                                   inflow, &
                                                   cfl_subcycle_ele, &
                                                   sfield_cv_options, &
                                                   save_pos = upwind_pos)
                            
                            ! Cache the phase face value if present and non linear iter > 1, else no need
                            if (cached_face_value_present .and. (max_nonlinear_iter_this_timestep > 1)) then
                            
                               cached_face_value_domain%value((f_ele_base - 1) + f_ele_counter) = &
                              &cached_face_value_domain%value((f_ele_base - 1) + f_ele_counter) + old_sfield_face_value
                            
                            end if
                         
                         else 
                         
                            old_sfield_face_value = cached_face_value_domain%value((f_ele_base - 1) + f_ele_counter)
                         
                         end if

                         ! Evaluate the face value for modrelperm (taking upwind) 
                         modrelperm_face_value = income*modrelperm_ele(oloc) + (1.0-income)*modrelperm_ele(iloc)

                         ! Find the density face value (taking upwind) 
                         den_face_value = income*den_ele(oloc) + (1.0-income)*den_ele(iloc)

                         face_value = detwei(ggi) * old_sfield_face_value * modrelperm_face_value * absperm_ele(1) * &
                        &dot_product((grad_pressure_face_quad(:,ggi) - den_face_value * grav_ele(:,1)), normgi)/ visc_ele(1)

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
               
               if (v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop

               ! Initialise the face value counter for this ele
               f_ele_counter = 0

               ! Find the face value ele base
               f_ele_base = cached_face_value_boundary%ele_base(sele)
               
               if (v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NORMAL_FLOW) then
               
                  bc_sele_val = ele_val(v_over_s_bc_value, sele)
               
               else

                  modrelperm_ele_bdy   = face_val(modified_relative_permeability, sele)
                  visc_ele_bdy         = face_val(viscosity, sele)
                  absperm_ele_bdy      = face_val(absolute_permeability, sele)
                  inv_char_len_ele_bdy = ele_val(inverse_characteristic_length, sele)
                  if (pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then
                      bc_sele_val = ele_val(pressure_bc_value, sele)     
                  end if
                  grad_pressure_face_quad_bdy      = face_val_at_quad(gradient_pressure, sele, gradp_cvbdyshape)         
                  iter_grad_pressure_face_quad_bdy = face_val_at_quad(iterated_gradient_pressure, sele, gradp_cvbdyshape)
                  grav_ele_bdy                     = face_val(gravity, sele) 
                  pressure_ele_bdy                 = face_val(pressure, sele) 

                  ! the latest and iterated DarcyVelocityOverSaturation at the quadrature points
                  ! determined from FE interpolation of each component, 
                  ! used to determine upwind (iterated) and also added in certain integrals (latest).
                  do dim = 1,ndim

                     v_over_s_face_quad_bdy(dim,:) =  &
- (face_val_at_quad(modified_relative_permeability, sele, p_cvbdyshape) * absperm_ele_bdy(1) / visc_ele_bdy(1)) * &
  (grad_pressure_face_quad_bdy(dim,:) - face_val_at_quad(density, sele, p_cvbdyshape) * grav_ele_bdy(dim,1))   

                     iter_v_over_s_face_quad_bdy(dim,:) =  &
- (face_val_at_quad(modified_relative_permeability, sele, p_cvbdyshape) * absperm_ele_bdy(1) / visc_ele_bdy(1)) * &
  (iter_grad_pressure_face_quad_bdy(dim,:) - face_val_at_quad(density, sele, p_cvbdyshape) * grav_ele_bdy(dim,1))   

                  end do
               
               end if

               if (sfield_bc_flag(sele) == S_BC_TYPE_WEAKDIRICHLET) then
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
                           
                           if (v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NORMAL_FLOW) then
                           
                              if (determine_face_value) then

                                 ! Determine the upwind saturation boundary values to use from the sign 
                                 ! of the normal flow BC. If positive, hence outflow, always 
                                 ! use domain value. If negative, hence inflow, use the boundary 
                                 ! value which is the weak BC value if appropriate                                 
                                 if (bc_sele_val(iloc) > 0.0) then
                                 
                                    old_sfield_face_value = old_sfield_subcycle_ele_bdy(iloc)
                                 
                                 else 
                                 
                                    old_sfield_face_value = ghost_old_sfield_subcycle_ele_bdy(iloc)
                                 
                                 end if 
                                 
                                 ! Cache the phase face value if present and non linear iter > 1, else no need
                                 if (cached_face_value_present .and. (max_nonlinear_iter_this_timestep > 1)) then

                                    cached_face_value_boundary%value((f_ele_base - 1) + f_ele_counter) = &
                                   &cached_face_value_boundary%value((f_ele_base - 1) + f_ele_counter) + old_sfield_face_value 

                                 end if 

                              else

                                 old_sfield_face_value = cached_face_value_boundary%value((f_ele_base - 1) + f_ele_counter)

                              end if 

                              s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - old_sfield_face_value * bc_sele_val(iloc) * detwei_bdy(ggi)
                           
                           else                            
                              
                              ! determine if the flow is in or out of the face at this quadrature
                              ! with respect to the normal orientation using the iterated v_over_s
                              iter_v_over_s_dot_n = dot_product(iter_v_over_s_face_quad_bdy(:,ggi), normal_bdy(:,ggi))
                              
                              if (determine_face_value) then
                              
                                 inflow = (iter_v_over_s_dot_n<=0.0)

                                 income = merge(1.0,0.0,inflow)

                                 old_sfield_face_value = income*ghost_old_sfield_subcycle_ele_bdy(iloc) + &
                                                         (1.0 - income)*old_sfield_subcycle_ele_bdy(iloc)
                                                               
                                 ! Cache the phase face value if present and non linear iter > 1, else no need
                                 if (cached_face_value_present .and. (max_nonlinear_iter_this_timestep > 1)) then

                                    cached_face_value_boundary%value((f_ele_base - 1) + f_ele_counter) = &
                                   &cached_face_value_boundary%value((f_ele_base - 1) + f_ele_counter) + old_sfield_face_value 

                                 end if                                 
                                 
                              else 
                              
                                 old_sfield_face_value = cached_face_value_boundary%value((f_ele_base - 1) + f_ele_counter)                              
                              
                              end if

                              v_over_s_dot_n = dot_product(v_over_s_face_quad_bdy(:,ggi), normal_bdy(:,ggi))

                              s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - old_sfield_face_value * v_over_s_dot_n * detwei_bdy(ggi)
if (.false.) then
                              modrelperm_absperm_over_visc = modrelperm_ele_bdy(iloc) * absperm_ele_bdy(1) / visc_ele_bdy(1)

                              if (pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then

                                 do jloc = 1, sfield%mesh%faces%shape%loc 

                                    s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - &
                                                            inv_char_len_ele_bdy(jloc) * &
                                                            p_cvbdyshape%n(jloc,ggi) * &
                                                            bc_sele_val(jloc) * &
                                                            modrelperm_absperm_over_visc * &
                                                            old_sfield_face_value * &
                                                            sum(normal_bdy(:,ggi)) * &
                                                            detwei_bdy(ggi) * &
                                                            weak_pressure_bc_coeff

                                 end do

                              end if 

                              do jloc = 1,sfield%mesh%faces%shape%loc 

                                 s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) + &
                                                         inv_char_len_ele_bdy(jloc) * &
                                                         p_cvbdyshape%n(jloc,ggi) * &
                                                         pressure_ele_bdy(jloc) * &
                                                         modrelperm_absperm_over_visc * &
                                                         old_sfield_face_value * &
                                                         sum(normal_bdy(:,ggi)) * &
                                                         detwei_bdy(ggi) * &
                                                         weak_pressure_bc_coeff

                              end do
end if
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
      real    :: income, v_over_s_dot_n, v_over_s_face_value_dot_n
      real    :: modrelperm_face_value, den_face_value, v_over_s_dot_n_face_value
      real,    dimension(1)                  :: visc_ele, absperm_ele
      real,    dimension(:,:),   allocatable :: grad_pressure_face_quad
      real,    dimension(:,:),   allocatable :: grad_pressure_ele
      real,    dimension(:),     allocatable :: v_over_s_local
      real,    dimension(:,:),   allocatable :: grav_ele
      real,    dimension(:),     allocatable :: den_ele
      real,    dimension(:),     allocatable :: modrelperm_ele      
      real,    dimension(:,:),   allocatable :: v_over_s_face_quad
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
      real,    dimension(:,:),   allocatable :: grav_ele_bdy
      real,    dimension(:,:),   allocatable :: grad_pressure_face_quad_bdy
      real,    dimension(:,:),   allocatable :: v_over_s_face_quad_bdy
      real,    dimension(:,:),   allocatable :: normal_bdy
      real,    dimension(:),     allocatable :: detwei_bdy
      real,    dimension(:,:),   allocatable :: x_ele_bdy
      real,    dimension(:),     allocatable :: cfl_rhs_local_bdy
      integer, dimension(:),     allocatable :: p_nodes_bdy
      real,    dimension(:),     allocatable :: bc_sele_val
      integer, parameter :: V_OVER_S_BC_TYPE_NORMAL_FLOW = 1, V_OVER_S_BC_TYPE_NO_NORMAL_FLOW = 2
            
      ewrite(1,*) 'Calculate CFL, Velocity and FractionalFlow fields'

      ! allocate arrays used in CFL assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      allocate(normal(di%ndim,di%x_cvshape%ngi))
      allocate(detwei(di%x_cvshape%ngi))      
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(cfl_rhs_local(ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(grav_ele(di%ndim,1))
      allocate(v_over_s_face_quad(di%ndim,di%p_cvshape%ngi))
      
      allocate(grav_ele_bdy(di%ndim,1))
      allocate(grad_pressure_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      allocate(v_over_s_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      allocate(detwei_bdy(di%x_cvbdyshape%ngi))
      allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      allocate(cfl_rhs_local_bdy(face_loc(di%pressure_mesh,1)))
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(p_nodes_bdy(face_loc(di%pressure_mesh,1)))
      allocate(bc_sele_val(face_loc(di%pressure_mesh,1)))
      
      allocate(den_ele(ele_loc(di%pressure_mesh,1)))
      allocate(modrelperm_ele(ele_loc(di%pressure_mesh,1)))
      allocate(grad_pressure_ele(di%ndim,1))
      allocate(v_over_s_local(ele_loc(di%pressure_mesh,1)))

      phase_loop: do p = 1, di%number_phase
         
         call zero(di%cfl(p)%ptr)
              
         vele_loop: do vele = 1, element_count(di%pressure_mesh)
            
            ! get the coordinate values for this element for each positions local node
            x_ele = ele_val(di%positions, vele)         

            ! get the coordinate values for this element for each quadrature point
            x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)

            ! The node indices of the pressure field
            p_nodes => ele_nodes(di%pressure_mesh, vele)

            ! The node indices of the positions projected to the pressure mesh
            x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)

            ! get the modrelperm ele values
            modrelperm_ele = ele_val(di%modified_relative_permeability(p)%ptr, vele)

            ! get the density value for this element
            den_ele = ele_val(di%density(p)%ptr, vele)
            
            ! The gravity values for this element for each direction
            grav_ele = ele_val(di%gravity, vele) 

            ! get the viscosity value for this element for this phase
            visc_ele = ele_val(di%viscosity(p)%ptr, vele)

            ! get the absolute permeability value for this element
            absperm_ele = ele_val(di%absolute_permeability, vele)         

            ! get the latest gradient pressure at the cv surface quadrature points for each direction for this phase
            grad_pressure_face_quad = ele_val_at_quad(di%gradient_pressure(p)%ptr, vele, di%gradp_cvshape)

            ! the latest DarcyVelocityOverSaturation at the quadrature points
            ! determined from FE interpolation of each component
            do dim = 1,di%ndim

               v_over_s_face_quad(dim,:) =  &
- (ele_val_at_quad(di%modified_relative_permeability(p)%ptr, vele, di%p_cvshape) * absperm_ele(1) / visc_ele(1)) * &
  (grad_pressure_face_quad(dim,:) - ele_val_at_quad(di%density(p)%ptr, vele, di%p_cvshape) * grav_ele(dim,1))   

            end do

            ! obtain the transformed determinant*weight and normals
            call transform_cvsurf_to_physical(x_ele, di%x_cvshape, detwei, normal, di%cvfaces)

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
                      ! with respect to the normal orientation using the latest v_over_s
                      v_over_s_dot_n = dot_product(v_over_s_face_quad(:,ggi), normgi(:))

                      inflow = (v_over_s_dot_n<=0.0)

                      income = merge(1.0,0.0,inflow)

                      ! Evaluate the face value for modrelperm (taking upwind) 
                      modrelperm_face_value = income*modrelperm_ele(oloc) + (1.0-income)*modrelperm_ele(iloc)

                      ! Find the density face value (taking upwind) 
                      den_face_value = income*den_ele(oloc) + (1.0-income)*den_ele(iloc)

                      v_over_s_dot_n_face_value = detwei(ggi) * modrelperm_face_value * absperm_ele(1) * &
                     &dot_product((grad_pressure_face_quad(:,ggi) - den_face_value * grav_ele(:,1)), normgi)/ visc_ele(1)

                      cfl_rhs_local(iloc) = cfl_rhs_local(iloc) + &
                                            abs(v_over_s_dot_n_face_value) * &
                                            (1.0 - income)

                      cfl_rhs_local(oloc) = cfl_rhs_local(oloc) + &
                                            abs(v_over_s_dot_n_face_value) * &
                                            income

                    end if check_visited

                  end do quadrature_loop

                end if is_neigh

              end do face_loop

            end do nodal_loop_i

            call addto(di%cfl(p)%ptr, p_nodes, cfl_rhs_local)

         end do vele_loop

         ! Get this phase v_over_s BC info - only for no_normal_flow and normal_flow which is special as it is a scalar
         call get_v_over_s_boundary_condition(di%darcy_velocity_over_saturation(p)%ptr, &
                                              (/"normal_flow   ", &
                                                "no_normal_flow"/), &
                                              di%bc_surface_mesh, &
                                              di%v_over_s_bc_value, &
                                              di%v_over_s_bc_flag)
         
         sele_loop: do sele = 1, surface_element_count(di%pressure_mesh)
            
            if (di%v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NO_NORMAL_FLOW) cycle sele_loop
            
            if (di%v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NORMAL_FLOW) then
            
               bc_sele_val = ele_val(di%v_over_s_bc_value, sele)
            
            else
            
               visc_ele_bdy                = face_val(di%viscosity(p)%ptr, sele)
               absperm_ele_bdy             = face_val(di%absolute_permeability, sele)
               grad_pressure_face_quad_bdy = face_val_at_quad(di%gradient_pressure(p)%ptr, sele, di%gradp_cvbdyshape)         
               grav_ele_bdy                = face_val(di%gravity, sele) 

               ! the latest DarcyVelocityOverSaturation at the quadrature points
               ! determined from FE interpolation of each component
               do dim = 1,di%ndim

                  v_over_s_face_quad_bdy(dim,:) =  &
- (face_val_at_quad(di%modified_relative_permeability(p)%ptr, sele, di%p_cvbdyshape) * absperm_ele_bdy(1) / visc_ele_bdy(1)) * &
  (grad_pressure_face_quad_bdy(dim,:) - face_val_at_quad(di%density(p)%ptr, sele, di%p_cvbdyshape) * grav_ele_bdy(dim,1))   

               end do
            
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
                        
                        if (di%v_over_s_bc_flag(sele) == V_OVER_S_BC_TYPE_NORMAL_FLOW) then
                        
                           v_over_s_face_value_dot_n = bc_sele_val(iloc)
                           
                           if (v_over_s_dot_n > 0.0) then
                              income = 0.0
                           else
                              income = 1.0
                           end if
                           
                        else
                        
                           ! determine if the flow is in or out of the face at this quadrature
                           ! with respect to the normal orientation using the latest v_over_s
                           v_over_s_face_value_dot_n = dot_product(v_over_s_face_quad_bdy(:,ggi), normal_bdy(:,ggi))

                           inflow = (v_over_s_dot_n<=0.0)

                           income = merge(1.0,0.0,inflow)
                                                      
                        end if
                        
                        cfl_rhs_local_bdy(iloc) = cfl_rhs_local_bdy(iloc) + &
                                                 &abs(v_over_s_face_value_dot_n) * detwei_bdy(ggi) * (1.0 - income)                     

                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop

            call addto(di%cfl(p)%ptr, p_nodes_bdy, cfl_rhs_local_bdy)

         end do sele_loop
                  
         di%cfl(p)%ptr%val = di%cfl(p)%ptr%val * di%dt / di%cv_mass_pressure_mesh_with_porosity%val

         ewrite_minmax(di%cfl(p)%ptr)

      end do phase_loop      

      
      ewrite(1,*) 'Calculate DarcyVelocityOverSaturation'
      
      ! Calculate the darcy_velocity_over_saturation field = - (1.0/sigma) * ( grad_pressure - den * grav), 
      ! where sigma = visc / (absperm * modrelperm)
      
      do p = 1, di%number_phase
                  
         do vele = 1, element_count(di%velocity_mesh)

            ! Find the element wise absolute_permeability
            absperm_ele = ele_val(di%absolute_permeability, vele)

            ! Find the element wise viscosity
            visc_ele = ele_val(di%viscosity(p)%ptr, vele)

            ! Find the element wise local values for modrelperm
            modrelperm_ele = ele_val(di%modified_relative_permeability(p)%ptr, vele)

            ! Find the element wise local values for density
            den_ele = ele_val(di%density(p)%ptr, vele)
         
            ! The gravity values for this element for each direction
            grav_ele = ele_val(di%gravity, vele) 

            ! Find the element wise gradient pressure 
            grad_pressure_ele = ele_val(di%gradient_pressure(p)%ptr, vele)

            do dim = 1,di%ndim

               do iloc = 1, di%velocity_mesh%shape%loc
               
                  v_over_s_local(iloc) = - (modrelperm_ele(iloc) * absperm_ele(1) / visc_ele(1)) * &
                                         & ( grad_pressure_ele(dim,1) - den_ele(iloc) * grav_ele(dim,1))
                                 
               end do

               call set(di%darcy_velocity_over_saturation(p)%ptr, &
                        dim, &
                        ele_nodes(di%velocity_mesh, vele), &
                        v_over_s_local)               
               
            end do

         end do
         
      end do 
                  
      do p = 1,di%number_phase
         ewrite_minmax(di%darcy_velocity_over_saturation(p)%ptr)
      end do
      
      ! calculate the darcy velocity = darcy_velocity_over_saturation * old_saturation
      do p = 1, di%number_phase         
         
         ewrite(1,*) 'Calculate DarcyVelocity for phase ',p
                  
         call zero(di%darcy_velocity(p)%ptr)
         
         do vele = 1, element_count(di%velocity_mesh)
                        
            do dim = 1, di%ndim
                              
               call addto(di%darcy_velocity(p)%ptr, &
                          dim, &
                          ele_nodes(di%velocity_mesh, vele), &
                          ele_val(di%old_saturation(p)%ptr, vele))
                           
            end do
            
         end do 

         call scale(di%darcy_velocity(p)%ptr, di%darcy_velocity_over_saturation(p)%ptr)
                  
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
      
      ! deallocate local variables as required
      deallocate(x_ele)
      deallocate(normal)
      deallocate(detwei)
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(cfl_rhs_local)
      deallocate(x_face_quad)
      deallocate(grad_pressure_face_quad)
      deallocate(grav_ele)
      deallocate(v_over_s_face_quad)
      
      deallocate(grav_ele_bdy)
      deallocate(grad_pressure_face_quad_bdy)
      deallocate(v_over_s_face_quad_bdy)
      deallocate(detwei_bdy)
      deallocate(normal_bdy)
      deallocate(x_ele_bdy)
      deallocate(p_nodes_bdy)
      deallocate(cfl_rhs_local_bdy)
      deallocate(bc_sele_val)      

      deallocate(den_ele)
      deallocate(modrelperm_ele)
      deallocate(grad_pressure_ele)
      deallocate(v_over_s_local)
             
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
   
   subroutine darcy_impes_initialise_cached_phase_face_value(di)
   
      !!< Initialise the cached_phase_face_value variable's

      type(darcy_impes_type), intent(inout) :: di
      
      ! local variable
      integer :: p, face_value_counter, vele, sele, iloc, face, gi, ggi
      logical, dimension(:), pointer :: notvisited

      ewrite(1,*) 'Allocate cached phase face value types'

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

      ewrite(1,*) 'Finished Allocate cached phase face value types'
      
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

      ewrite(1,*) 'Calculate inverse characteristic length'
      
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
      
      ewrite_minmax(di%inverse_characteristic_length)
      
      ewrite(1,*) 'Finished Calculate inverse characteristic length'
      
   end subroutine darcy_impes_calculate_inverse_characteristic_length

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_modified_relative_permeability(di)
   
      !!< Calculate the modified relative permeability = relperm / max( S, S_min) for each phase
            
      type(darcy_impes_type), intent(inout) :: di 
      
      ! local variable
      integer :: p, i

      ewrite(1,*) 'Calculate modified relative permeabilities'
      
      phase_loop: do p = 1, di%number_phase
         
         node_loop: do i = 1, node_count(di%relative_permeability(p)%ptr)

            di%modified_relative_permeability(p)%ptr%val(i) = &
            di%relative_permeability(p)%ptr%val(i) / &
            max(di%saturation(p)%ptr%val(i), di%min_sat_value_for_mod_relperm)
            
         end do node_loop
         
         ewrite_minmax(di%modified_relative_permeability(p)%ptr)
         
      end do phase_loop

      ewrite(1,*) 'Finished Calculate modified relative permeabilities'
      
   end subroutine darcy_impes_calculate_modified_relative_permeability

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
