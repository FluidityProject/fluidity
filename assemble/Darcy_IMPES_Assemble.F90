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
   use sparse_tools_petsc
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
   use sparse_matrices_fields
   use transform_elements, only: transform_cvsurf_to_physical, &
                                 transform_cvsurf_facet_to_physical
   use state_module
   use fldebug
   use spud
   use porous_media
   use parallel_tools
   use adaptive_timestepping
   use halos
   use signal_vars, only : SIG_INT
   use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN

   !!*****03 July 2014 LCai ************use Leaching Chemical model****!
   use darcy_impes_leaching_chemical_model, only: add_leach_chemical_prog_src_to_rhs
   use darcy_impes_assemble_type
   !!************Finish********************

   implicit none

   private

   public :: darcy_impes_copy_to_old, &
             darcy_impes_copy_to_iterated, &
             darcy_impes_assemble_and_solve_part_one, &
             darcy_impes_assemble_and_solve_part_two, &
             darcy_impes_assemble_and_solve_part_three, &             
             darcy_impes_calculate_gradient_pressures, &
             darcy_impes_calculate_non_first_phase_pressures, &
             darcy_impes_calculate_phase_one_saturation_diagnostic, &
             darcy_impes_calculate_vel_mob_ff_and_cfl_fields, &
             darcy_impes_calculate_sum_saturation, &
             darcy_impes_calculate_relperm_fields, &
             darcy_impes_calculate_densities, &
             darcy_impes_calculate_cflnumber_field_based_dt, &
             darcy_impes_initialise_cached_face_value, &
             darcy_impes_calculate_inverse_characteristic_length, &
             darcy_impes_calculate_relperm_den_first_face_values, &
             darcy_impes_calculate_relperm_isub_face_values, &
             darcy_impes_get_entire_sfield_boundary_condition, &
             darcy_impes_get_v_boundary_condition, &
             darcy_impes_assemble_saturation_rhs_adv, &
             darcy_impes_assemble_and_solve_phase_pressures, &
             darcy_impes_calculate_divergence_total_darcy_velocity, &
             darcy_impes_calculate_inverse_cv_sa, &
             darcy_trans_MIM_assemble_and_solve_mobile_saturation, &  !**LCai 27 July 2013
             darcy_trans_assemble_galerkin_projection_elemesh_to_pmesh ! ** LCai 22 Aug 2013
   
   ! Parameters defining Darcy IMPES cached options
   integer, parameter, public :: RELPERM_CORRELATION_POWER                 = 1, &
                                 RELPERM_CORRELATION_COREY2PHASE           = 2, &
                                 RELPERM_CORRELATION_COREY2PHASEOPPOSITE   = 3, &
                                 RELPERM_CORRELATION_MINERAL               = 4, &
                                 RELPERM_CORRELATION_VANGENUCHTEN          = 5, &
                                 RELPERM_CORRELATION_JACKSON2PHASE         = 6, &
                                 RELPERM_CORRELATION_JACKSON2PHASEOPPOSITE = 7
   
   ! Parameters defining Darcy IMPES CV face value schemes
   integer, parameter, public :: DARCY_IMPES_CV_FACEVALUE_NONE                        = 0, &
                                 DARCY_IMPES_CV_FACEVALUE_FIRSTORDERUPWIND            = 1, &
                                 DARCY_IMPES_CV_FACEVALUE_FINITEELEMENT               = 2, &
                                 DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATUPWIND        = 3, &
                                 DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATFINITEELEMENT = 4, &
                                 DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATCORRELATION   = 5
   
   contains

! ----------------------------------------------------------------------------
  
   subroutine darcy_impes_copy_to_old(di)
            
      !!< Copy to old darcy impes data not associated with di%state

      type(darcy_impes_type), intent(inout) :: di
            
      call set(di%cv_mass_pressure_mesh_with_old_porosity, di%cv_mass_pressure_mesh_with_porosity)
      
      ! *****************22 Aug 2013 *** LCai ********************************** 
      ! copy the porosity used to calculate the Mobile_immobile_model to old
      if (size(di%MIM_options%immobile_prog_sfield) > 0) then
        
        if (di%prt_is_constant) then
          di%old_porosity_cnt = di%porosity_cnt
        else
          call set(di%old_porosity_pmesh, di%porosity_pmesh)
        end if

      end if
      ! *************Finish **** LCai ******************************************
   end subroutine darcy_impes_copy_to_old

! ----------------------------------------------------------------------------
  
   subroutine darcy_impes_copy_to_iterated(di)
            
      !!< Copy to iterated darcy impes data not associated with di%state

      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p
      
      do p = 1, di%number_phase
      
         call set(di%iterated_gradient_pressure(p)%ptr, di%gradient_pressure(p)%ptr)
      
      end do
      
   end subroutine darcy_impes_copy_to_iterated

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve_part_one(di, have_dual)
      
      !!< Assemble and solve the Darcy equations using an CIMPESS algorithm
      !!< which is a modification of the IMPES to include consistent subcycling.
      !!< Part one: form cv mass and initial saturation solve for subcycling
      
      type(darcy_impes_type), intent(inout) :: di
      logical ,               intent(in)    :: have_dual      
      
      ! local variables
      logical :: form_new_subcycle_relperm_face_values
      
      ewrite(1,*) 'Start Darcy IMPES assemble and solve part one'
                  
      ! Calculate the latest CV mass on the pressure mesh with porosity included
      call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_porosity, di%porosity)      
      
      ! If it is the first iteration and there is subcycling (>1)
      ! and a consistent global continuity then solve the saturations 
      if ((di%nonlinear_iter == 1) .and. &
          (di%subcy_opt_sat%number > 1) .and. &
          (di%subcy_opt_sat%consistent)) then
         
         form_new_subcycle_relperm_face_values  = .true.
         
         call darcy_impes_assemble_and_solve_phase_saturations(di, &
                                                               form_new_subcycle_relperm_face_values, &
                                                               have_dual)
         
      end if
                        
      ewrite(1,*) 'Finished Darcy IMPES assemble and solve part one'
           
   end subroutine darcy_impes_assemble_and_solve_part_one

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve_part_two(di, di_dual, have_dual, solve_dual_pressure)
      
      !!< Assemble and solve the Darcy equations using an CIMPESS algorithm
      !!< which is a modification of the IMPES to include consistent subcycling.
      !!< Part two: solve for the pressure
      
      type(darcy_impes_type), intent(inout) :: di
      type(darcy_impes_type), intent(inout) :: di_dual
      logical ,               intent(in)    :: have_dual      
      logical ,               intent(in)    :: solve_dual_pressure
            
      ewrite(1,*) 'Start Darcy IMPES assemble and solve part two'
                  
      ! Assemble and solve the phase pressures
      call darcy_impes_assemble_and_solve_phase_pressures(di, di_dual, have_dual, solve_dual_pressure)
                  
      ewrite(1,*) 'Finished Darcy IMPES assemble and solve part two'
           
   end subroutine darcy_impes_assemble_and_solve_part_two

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve_part_three(di, have_dual)
      
      !!< Assemble and solve the Darcy equations using an CIMPESS algorithm
      !!< which is a modification of the IMPES to include consistent subcycling.
      !!< Part three: find pressure gradients, calc div.u, solve for saturation
      !!<             and generic scalar fields, calc sum S, calc density, calc relperm.
      
      type(darcy_impes_type), intent(inout) :: di
      logical ,               intent(in)    :: have_dual      
      
      ! local variables
      logical :: form_new_subcycle_relperm_face_values
      
      ewrite(1,*) 'Start Darcy IMPES assemble and solve part three'
                  
      ! Calculate the gradient pressures 
      call darcy_impes_calculate_gradient_pressures(di)

      ! Calculate the divergence total darcy velocity
      call darcy_impes_calculate_divergence_total_darcy_velocity(di)
      
      ! Assemble and solve the phase saturations
            
      if ((di%nonlinear_iter == di%max_nonlinear_iter_this_timestep) .and. &
          (di%subcy_opt_sat%consistent)) then

         form_new_subcycle_relperm_face_values = .false.

      else 

         form_new_subcycle_relperm_face_values = .true.

      end if 
            
      call darcy_impes_assemble_and_solve_phase_saturations(di, &
                                                            form_new_subcycle_relperm_face_values, &
                                                            have_dual)
      
      call darcy_impes_assemble_and_solve_generic_prog_sfields(di)
      
      ! Calculate the sum of the saturations
      call darcy_impes_calculate_sum_saturation(di)
      
      ! Calculate the density field of each phase
      call darcy_impes_calculate_densities(di)

      ! Calculate the relative permeabilities of each phase
      call darcy_impes_calculate_relperm_fields(di)
            
      ewrite(1,*) 'Finished Darcy IMPES assemble and solve part three'
           
   end subroutine darcy_impes_assemble_and_solve_part_three

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve_phase_pressures(di, di_dual, have_dual, solve_dual_pressure)
      
      !!< Assemble and solve the phase pressures. The first phase pressure 
      !!< is solved via a matrix equation (if prognostic) and the other
      !!< phases pressures are calculated from the first phase pressure
      !!< and their own capilliary pressures
      
      type(darcy_impes_type), intent(inout) :: di
      type(darcy_impes_type), intent(inout) :: di_dual
      logical ,               intent(in)    :: have_dual  
      logical ,               intent(in)    :: solve_dual_pressure
      
      ! Assemble and solve the first phase pressure if it is prognostic
      if (di%first_phase_pressure_prognostic) call darcy_impes_assemble_and_solve_first_phase_pressure(di, di_dual, have_dual, solve_dual_pressure)
      
      ! Calculate the non first phase pressure's
      call darcy_impes_calculate_non_first_phase_pressures(di) 
      if (have_dual)  call darcy_impes_calculate_non_first_phase_pressures(di_dual)     
      
      ! Calculate the average pressure
      call darcy_impes_calculate_average_pressure(di)
      if (have_dual) call darcy_impes_calculate_average_pressure(di_dual)
      
   end subroutine darcy_impes_assemble_and_solve_phase_pressures

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve_first_phase_pressure(di, di_dual, have_dual, solve_dual_pressure)
      
      !!< Assemble and solve the first phase pressure. A source is included 
      !!< due to the capilliary pressures of non first phases as well as gravity 
      !!< and the rate of change of porosity and the individual phase sources.
      !!< Face values for relative permeability and density have already been evaluated.
      !!< If have_dual and solve_dual_pressure are both true then a block matrix using 
      !!< di%dual_block_pressure_matrix is solved for. If have_dual is true but 
      !!< solve_dual_pressure is false then a single matrix using di%pressure_matrix 
      !!< is solved for the prime pressure phase 1 then the dual pressure phase 1 
      !!< is set to this. If have_dual is false then di%pressure_matrix is used.
      !!< (This is because an error is occuring for petsc_csr_matrix in 2d)
      
      type(darcy_impes_type), intent(inout) :: di
      type(darcy_impes_type), intent(inout) :: di_dual
      logical ,               intent(in)    :: have_dual  
      logical ,               intent(in)    :: solve_dual_pressure
      
      ! local variables
      type(scalar_field_pointer), dimension(:), pointer :: all_pressure
      type(scalar_field_pointer), dimension(:), pointer :: all_rhs
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
      integer, parameter :: V_BC_TYPE_PRESCRIBED_NORMAL_FLOW = 1, V_BC_TYPE_NO_NORMAL_FLOW   = 2
            
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
      
      ! Initialise the matrix and rhs and solution pointer
      call zero(di%pressure_rhs)      
      
      if (have_dual .and. solve_dual_pressure) then
         call zero(di%dual_block_pressure_matrix)
         call zero(di_dual%pressure_rhs)
         
         allocate(all_pressure(2))
         allocate(all_rhs(2))
         all_pressure(1)%ptr => di%pressure(1)%ptr
         all_pressure(2)%ptr => di_dual%pressure(1)%ptr
         all_rhs(1)%ptr => di%pressure_rhs
         all_rhs(2)%ptr => di_dual%pressure_rhs         
      else
         call zero(di%pressure_matrix)
         allocate(all_pressure(1))
         allocate(all_rhs(1)) 
         all_pressure(1)%ptr => di%pressure(1)%ptr
         all_rhs(1)%ptr => di%pressure_rhs                 
      end if

      ewrite(1,*) 'Add phase sources to global continuity'
      
      src_phase_loop: do p = 1, di%number_phase
         
         ! Should this take account of subcycling?!?!
         call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_source, &
              di%saturation_source(p)%ptr)
         
         ! Should this include porosity ?!?!
         call addto(di%pressure_rhs, di%cv_mass_pressure_mesh_with_source)
         
         ewrite(2,*) '   Phase = ', p
         ewrite(2,*) '      Added saturation source, norm2 = ', &
              norm2(di%saturation_source(p)%ptr, di%positions)
         ewrite(2,*) '      Resulting cv_mass_pressure_mesh_with_source, norm2 = ', &
              norm2(di%cv_mass_pressure_mesh_with_source, di%positions)
         ewrite(2,*) '      Resulting pressure_rhs, norm2 = ', &
              norm2(di%pressure_rhs, di%positions)
         
      end do src_phase_loop
      
      if (have_dual .and. solve_dual_pressure) then
      
         dual_src_phase_loop: do p = 1, di_dual%number_phase

            ! Should this take account of subcycling?!?!
            call compute_cv_mass(di_dual%positions, di_dual%cv_mass_pressure_mesh_with_source, di_dual%saturation_source(p)%ptr)

            ! Should this include porosity ?!?!
            call addto(di_dual%pressure_rhs, di_dual%cv_mass_pressure_mesh_with_source)
         
         ewrite(2,*) '   Dual p-field.  Phase = ', p
         ewrite(2,*) '      Added saturation source, norm2 = ', &
              norm2(di%saturation_source(p)%ptr, di%positions)
         ewrite(2,*) '      Resulting cv_mass_pressure_mesh_with_source, norm2 = ', &
              norm2(di%cv_mass_pressure_mesh_with_source, di%positions)
         ewrite(2,*) '      Resulting pressure_rhs, norm2 = ', &
              norm2(di%pressure_rhs, di%positions)

         end do dual_src_phase_loop
      
      end if

      ! Add the transmissibility lambda dual to each block diagonal      
      if (have_dual .and. solve_dual_pressure) then         
         transmiss_phase_loop: do p = 1, di_dual%number_phase
            call compute_cv_mass(di_dual%positions, di_dual%cv_mass_pressure_mesh_with_lambda_dual, di_dual%transmissibility_lambda_dual(p)%ptr)

            call addto_diag(di%dual_block_pressure_matrix, 1, 1, di_dual%cv_mass_pressure_mesh_with_lambda_dual, scale =  1.0)
            call addto_diag(di%dual_block_pressure_matrix, 1, 2, di_dual%cv_mass_pressure_mesh_with_lambda_dual, scale = -1.0)
            call addto_diag(di%dual_block_pressure_matrix, 2, 1, di_dual%cv_mass_pressure_mesh_with_lambda_dual, scale = -1.0)
            call addto_diag(di%dual_block_pressure_matrix, 2, 2, di_dual%cv_mass_pressure_mesh_with_lambda_dual, scale =  1.0)
            
            if (p > 1) then
               ! Add lambda*(PC_dual - PC_prime) to prime rhs
               call set(di%rhs, di_dual%capilliary_pressure(p)%ptr)
               call addto(di%rhs, di%capilliary_pressure(p)%ptr, scale = -1.0)
               call scale(di%rhs, di_dual%cv_mass_pressure_mesh_with_lambda_dual)               
               call addto(di%pressure_rhs, di%rhs)

               ! Add lambda*(PC_prime- PC_dual) to dual rhs
               call set(di%rhs, di%capilliary_pressure(p)%ptr)
               call addto(di%rhs, di_dual%capilliary_pressure(p)%ptr, scale = -1.0)
               call scale(di%rhs, di_dual%cv_mass_pressure_mesh_with_lambda_dual)
               
               call addto(di_dual%pressure_rhs, di%rhs)               
            end if
         end do transmiss_phase_loop
      end if
                  
      ewrite(1,*) 'Add rate of change of porosity to global continuity equation'
      
      ! Add rate of change of porosity to rhs.  First, -dpor = por{n} - por{n+1}
      call set(di%work_array_of_size_pressure_mesh, &
           di%cv_mass_pressure_mesh_with_old_porosity)
      call addto(di%work_array_of_size_pressure_mesh, &
           di%cv_mass_pressure_mesh_with_porosity, scale = -1.0)
      ! then rhs += -dpor/dt
      call addto(di%pressure_rhs, &
           di%work_array_of_size_pressure_mesh, scale = 1.0/di%dt)
      ewrite(2,*) '      -dpor, norm2 = ', &
           norm2(di%work_array_of_size_pressure_mesh, di%positions)
      ewrite(2,*) '      1/dt = ', 1.0/di%dt
      ewrite(2,*) '      Resulting pressure_rhs, norm2 = ', &
           norm2(di%pressure_rhs, di%positions)
      
      if (have_dual .and. solve_dual_pressure) then
         ! repeat 
         call set(di_dual%work_array_of_size_pressure_mesh, &
              di_dual%cv_mass_pressure_mesh_with_old_porosity)
         call addto(di_dual%work_array_of_size_pressure_mesh, &
              di_dual%cv_mass_pressure_mesh_with_porosity, scale = -1.0)
         call addto(di_dual%pressure_rhs, di_dual%&
              work_array_of_size_pressure_mesh, scale = 1.0/di%dt)
         ewrite(2,*) '      Dual p-field.  -dpor, norm2 = ', &
              norm2(di_dual%work_array_of_size_pressure_mesh, di_dual%positions)
         ewrite(2,*) '      1/dt = ', 1.0/di%dt
         ewrite(2,*) '      Resulting pressure_rhs, norm2 = ', &
              norm2(di_dual%pressure_rhs, di_dual%positions)
      end if
      
      ! Assemble a contribution from each phase to form a global continuity equation to solve for first phase pressure
      phase_loop: do p = 1, di%number_phase
         
         ewrite(1,*) 'Assemble volume contribution to global continuity from phase ',p
            
         ! Loop volume elements assembling local contributions     
         vol_element_loop: do vele = 1, di%number_vele
                       
            call pressure_volume_element_local_assemble(di) 
            
            ! Add volume element contribution to global pressure matrix and rhs            
            if (have_dual .and. solve_dual_pressure) then 
               call addto(di%dual_block_pressure_matrix, 1, 1, p_nodes, p_nodes, p_mat_local)
            
               call addto(di%pressure_rhs, p_nodes, p_rhs_local)
 
               call pressure_volume_element_local_assemble(di_dual) 
               
               call addto(di%dual_block_pressure_matrix, 2, 2, p_nodes, p_nodes, p_mat_local)
               
               call addto(di_dual%pressure_rhs, p_nodes, p_rhs_local)
            
            else 

              call addto(di%pressure_matrix, p_nodes, p_nodes, p_mat_local)
            
              call addto(di%pressure_rhs, p_nodes, p_rhs_local)
            
            end if

         end do vol_element_loop
         
      end do phase_loop
                  
      ! Add normal v BC integrals for each phase, else include if required 
      ! a weak pressure BC (which if not given is assumed zero).

      ! Get the pressure BC - if no v and weak pressure then extra integrals are added
      ! - also get strong pressure BC as this implies no need to add any surface integrals
      call darcy_impes_get_entire_sfield_boundary_condition(di%pressure(1)%ptr, &
                                                            (/"weakdirichlet", &
                                                              "dirichlet    "/), &
                                                            di%bc_surface_mesh, &
                                                            di%pressure_bc_value, &
                                                            di%pressure_bc_flag)
      
      if (have_dual .and. solve_dual_pressure) then
         call darcy_impes_get_entire_sfield_boundary_condition(di_dual%pressure(1)%ptr, &
                                                               (/"weakdirichlet", &
                                                                 "dirichlet    "/), &
                                                               di_dual%bc_surface_mesh, &
                                                               di_dual%pressure_bc_value, &
                                                               di_dual%pressure_bc_flag)      
      end if 
      
      phase_loop_bc: do p = 1, di%number_phase

         ewrite(1,*) 'Assemble boundary contribution to global continuity from phase ',p
                           
         ! Get this phase v BC info - only for no_normal_flow and prescribed_normal_flow which is special as it is a scalar
         call darcy_impes_get_v_boundary_condition(di%darcy_velocity(p)%ptr, &
                                                   (/"prescribed_normal_flow", &
                                                     "no_normal_flow        "/), &
                                                   di%bc_surface_mesh, &
                                                   di%v_bc_value, &
                                                   di%v_bc_flag)
         
         if (have_dual .and. solve_dual_pressure) then         
            call darcy_impes_get_v_boundary_condition(di_dual%darcy_velocity(p)%ptr, &
                                                      (/"prescribed_normal_flow", &
                                                        "no_normal_flow        "/), &
                                                      di_dual%bc_surface_mesh, &
                                                      di_dual%v_bc_value, &
                                                      di_dual%v_bc_flag)
         end if

         sele_loop: do sele = 1, di%number_sele
            
            call pressure_surface_element_local_assemble(di)
            
            if (have_dual .and. solve_dual_pressure) then            

               call addto(di%dual_block_pressure_matrix, 1, 1, p_nodes_bdy, p_nodes_bdy, p_matrix_local_bdy)

               call addto(di%pressure_rhs, p_nodes_bdy, p_rhs_local_bdy)

               call pressure_surface_element_local_assemble(di_dual) 
                       
               call addto(di%dual_block_pressure_matrix, 2, 2, p_nodes_bdy, p_nodes_bdy, p_matrix_local_bdy)

               call addto(di_dual%pressure_rhs, p_nodes_bdy, p_rhs_local_bdy)            
            
            else

               call addto(di%pressure_matrix, p_nodes_bdy, p_nodes_bdy, p_matrix_local_bdy)

               call addto(di%pressure_rhs, p_nodes_bdy, p_rhs_local_bdy)
            
            end if

         end do sele_loop      

      end do phase_loop_bc
            
      ! Apply any strong dirichlet BC's - that can only be applied to the first phase pressure
      ! of both the prime and dual porous media
      if (have_dual .and. solve_dual_pressure) then
         call apply_dirichlet_conditions(di%dual_block_pressure_matrix, all_rhs, all_pressure)
      else
         call apply_dirichlet_conditions(di%pressure_matrix, all_rhs(1)%ptr, all_pressure(1)%ptr)      
      end if
      
      ! Solve the pressure(s)
      if (have_dual .and. solve_dual_pressure) then
         call petsc_solve(all_pressure, di%dual_block_pressure_matrix, all_rhs, di%state(1))
      else
         call petsc_solve(all_pressure(1)%ptr, di%pressure_matrix, all_rhs(1)%ptr, di%state(1))      
      end if
      
      ! Set the strong BC nodes to the values to be consistent
      call set_dirichlet_consistent(di%pressure(1)%ptr) 
      
      if (have_dual) then
         
         if (solve_dual_pressure) then
         
            call set_dirichlet_consistent(di_dual%pressure(1)%ptr) 
         
         else
            
            call set(di_dual%pressure(1)%ptr, di%pressure(1)%ptr)
         
         end if
         
      end if
          
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

      deallocate(all_pressure)
      deallocate(all_rhs)
      
      ewrite(1,*) 'Finished solve first phase Pressure'
      
      ewrite_minmax(di%pressure(1)%ptr)
      
      if (have_dual) then
         ewrite_minmax(di_dual%pressure(1)%ptr)
      end if
      
      contains
      
      ! --------------------------------------------------------------------------------------
      
         subroutine pressure_volume_element_local_assemble(di) 

            !!< Assemble the local volume element pressure matrix and rhs
         
            type(darcy_impes_type), intent(inout) :: di
                        
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
                           if (di%subcy_opt_sat%have .and. di%subcy_opt_sat%consistent) then

                              face_value = detwei(ggi) * &
                                           sum(di%cached_face_value%relperm(:,ggi,vele,p)) * &
                                           absperm_ele(1) / &
                                           visc_ele(1)
                           
                           else 
                           
                              face_value = detwei(ggi) * &
                                           di%cached_face_value%relperm(1,ggi,vele,p) * &
                                           absperm_ele(1) / &
                                           visc_ele(1)
                           
                           end if
                           
                           ! Form the local matrix given by - n_i . sum_{phase} ( relperm*absperm/visc ) dP_1/dx_j
                           do jloc = 1,di%pressure_mesh%shape%loc
                              
                              p_mat_local(iloc,jloc) = p_mat_local(iloc,jloc) - &
                                                       sum(p_dshape(jloc, ggi, :)*normgi, 1) * &
                                                       face_value

                              p_mat_local(oloc,jloc) = p_mat_local(oloc,jloc) - &
                                                       sum(p_dshape(jloc, ggi, :)*(-normgi), 1) * &
                                                       face_value

                           end do

                           ! Add gravity term to rhs = - n_i . sum_{phase} ( relperm*absperm/visc ) * den*grav

                           ! Find g dot n
                           g_dot_n = dot_product(grav_ele(:,1), normgi)

                           p_rhs_local(iloc)  = p_rhs_local(iloc) - &
                                                face_value * &
                                                di%cached_face_value%den(ggi,vele,p) * &
                                                g_dot_n
                           
                           ! Add opposite node value with change of sign due to opposite normal direction
                           p_rhs_local(oloc)  = p_rhs_local(oloc) + &
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

                              ! Add opposite node value with change of sign due to opposite normal direction
                              p_rhs_local(oloc) = p_rhs_local(oloc) - &
                                                  face_value * &
                                                  grad_cap_p_dot_n

                           end if

                        end if check_visited

                     end do quadrature_loop

                  end if is_neigh

               end do face_loop

            end do nodal_loop_i
            
         end subroutine pressure_volume_element_local_assemble 

      ! --------------------------------------------------------------------------------------
      
         subroutine pressure_surface_element_local_assemble(di) 

            !!< Assemble the local surface element pressure matrix and rhs
         
            type(darcy_impes_type), intent(inout) :: di
                        
            p_matrix_local_bdy = 0.0
            p_rhs_local_bdy    = 0.0

            p_nodes_bdy = face_global_nodes(di%pressure_mesh, sele)
            
            ! A no_normal_flow BC adds nothing to the matrix and rhs so return zeros
            if (di%v_bc_flag(sele) == V_BC_TYPE_NO_NORMAL_FLOW) return
            
            ! If phase 1 and strong pressure BC then no need to add integrals so return zeros
            if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_DIRICHLET) return
            
            ! obtain the transformed determinant*weight and normals
            if (di%cached_face_value%cached_detwei_normal) then

               detwei_bdy => di%cached_face_value%detwei_bdy(:,sele)
               
               normal_bdy => di%cached_face_value%normal_bdy(:,:,sele)
            
            else

               x_ele     = ele_val(di%positions, face_ele(di%positions, sele))
               x_ele_bdy = face_val(di%positions, sele)
            
               call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)
            
            end if
            
            if (di%v_bc_flag(sele) == V_BC_TYPE_PRESCRIBED_NORMAL_FLOW) then
                              
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

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi
                                                
                        prescribed_normal_flow: if (di%v_bc_flag(sele) == V_BC_TYPE_PRESCRIBED_NORMAL_FLOW) then
                           
                           ! If have prescribed_normal_flow BC for this phase then include CV integral of value
                           
                           if (di%subcy_opt_sat%have .and. di%subcy_opt_sat%consistent) then
                                                      
                              p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - &
                                                      di%subcy_opt_sat%number * &
                                                      bc_sele_val(iloc) * &
                                                      detwei_bdy(ggi)
                           
                           else 
                           
                              p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - &
                                                      bc_sele_val(iloc) * &
                                                      detwei_bdy(ggi)
                           
                           end if
                           
                        else prescribed_normal_flow
                                                                               
                           ! Form the face value = detwei * (relperm*absperm/visc)
                           if (di%subcy_opt_sat%have .and. di%subcy_opt_sat%consistent) then
                           
                              face_value = detwei_bdy(ggi) * &
                                           sum(di%cached_face_value%relperm_bdy(:,ggi,sele,p)) * &
                                           absperm_ele_bdy(1) / &
                                           visc_ele_bdy(1) 
                           
                           else 
                           
                              face_value = detwei_bdy(ggi) * &
                                           di%cached_face_value%relperm_bdy(1,ggi,sele,p) * &
                                           absperm_ele_bdy(1) / &
                                           visc_ele_bdy(1) 
                           
                           end if
                           
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
                           pressure_weak: if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then
                                                            
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

                           end if pressure_weak

                        end if prescribed_normal_flow
                        
                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop            
            
         end subroutine pressure_surface_element_local_assemble 

      ! --------------------------------------------------------------------------------------
      
   end subroutine darcy_impes_assemble_and_solve_first_phase_pressure

! ----------------------------------------------------------------------------

   subroutine darcy_impes_get_v_boundary_condition(v, &
                                                   types, &
                                                   bc_surface_mesh, &
                                                   v_bc_value, &
                                                   v_bc_flag)

      !!< Form the data associated with any v BC by returning 
      !!< full surface mesh arrays of a flag indicating  either 
      !!< no_normal_flow or prescribed_normal_flow and the value.

      type(vector_field),               intent(in),   target :: v
      character(len=*),   dimension(:), intent(in)           :: types
      type(mesh_type),                  intent(in)           :: bc_surface_mesh
      type(scalar_field),               intent(inout)        :: v_bc_value
      integer,            dimension(:), intent(inout)        :: v_bc_flag

      ! Local variables
      character(len=FIELD_NAME_LEN)                       :: bctype
      type(scalar_field)                                  :: scalar_surface_field
      integer,                      dimension(:), pointer :: surface_element_list
      integer                                             :: i, j, k, sele

      ewrite(1,*) 'Get DarcyVelocity boundary condition data'

      ! Zero the prescribed_normal_flow value surface field on whole boundary mesh  
      call zero(v_bc_value)

      ! Initialise flag for whether surface element has prescribed_normal_flow BC
      v_bc_flag = 0

      ! Loop each BC object instance for the v
      ! (May have multiple prescribed_normal_flow BC's applied to different surface id's)
      BC_loop: do i=1, get_boundary_condition_count(v)

         ! Get this BC info
         call get_boundary_condition(v, i, type = bctype, &
                                     surface_element_list = surface_element_list)

         ! check this is a prescribed_normal_flow BC  or no_normal_flow (nothing else is permitted)
         if ((trim(bctype) /= 'prescribed_normal_flow') .and. (trim(bctype) /= 'no_normal_flow')) then
            FLAbort('Have unknown BC type for a DarcyVelocity')
         end if

         ! Extract the surface_field for this BC for prescribed_normal_flow
         if (trim(bctype) == 'prescribed_normal_flow') then
            if (associated(v%bc%boundary_condition(i)%surface_fields)) then
               scalar_surface_field = extract_scalar_field(&
                    v%bc%boundary_condition(i)%surface_fields(1), 1)
            else
               FLAbort('Component surface_fields for DarcyVelocity BC type not associated')
            end if
         end if 
         
         ! Loop the surface elements associated with this BC instance
         ! and place the required BC value in whole boundary field
         ! for the prescribed_normal_flow BC type
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
               FLAbort('Cannot find no_normal_flow or prescribed_normal_flow bctype for DarcyVelocity')
            end if
            
            v_bc_flag(sele) = j

            ! Set the prescribed_normal_flow field values from this BC for its sele
            if (trim(bctype) == 'prescribed_normal_flow') then
               call set(v_bc_value, &
                        ele_nodes(bc_surface_mesh, sele), &
                        ele_val(scalar_surface_field, k))
            end if
            
         end do BC_sele_loop

      end do BC_loop

      ewrite_minmax(v_bc_value)

      ewrite(1,*) 'Finished get DarcyVelocity boundary condition data'

   end subroutine darcy_impes_get_v_boundary_condition
  
! ----------------------------------------------------------------------------

   subroutine darcy_impes_get_entire_sfield_boundary_condition(sfield, &
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
      
   end subroutine darcy_impes_get_entire_sfield_boundary_condition

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
      integer, parameter :: V_BC_TYPE_PRESCRIBED_NORMAL_FLOW = 1, V_BC_TYPE_NO_NORMAL_FLOW = 2
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

      ! Get the pressure BC - required if no v given and weak pressure dirichlet given for extra integrals
      call darcy_impes_get_entire_sfield_boundary_condition(di%pressure(1)%ptr, &
                                                            (/"weakdirichlet"/), &
                                                            di%bc_surface_mesh, &
                                                            di%pressure_bc_value, &
                                                            di%pressure_bc_flag)

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
                           if (di%subcy_opt_sat%have .and. di%subcy_opt_sat%consistent) then
                           
                              face_value = detwei(ggi) * &
                                           sum(di%cached_face_value%relperm(:,ggi,vele,p)) * &
                                           absperm_ele(1) / &
                                           visc_ele(1) 
                           
                           else
                           
                              face_value = detwei(ggi) * &
                                           di%cached_face_value%relperm(1,ggi,vele,p) * &
                                           absperm_ele(1) / &
                                           visc_ele(1) 
                           
                           end if
                           
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

         ! Get this phase v BC info - only for no_normal_flow and prescribed_normal_flow which is special as it is a scalar
         call darcy_impes_get_v_boundary_condition(di%darcy_velocity(p)%ptr, &
                                                   (/"prescribed_normal_flow", &
                                                     "no_normal_flow        "/), &
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
                     
            if (di%v_bc_flag(sele) == V_BC_TYPE_PRESCRIBED_NORMAL_FLOW) then
              
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
                                                
                        prescribed_normal_flow: if (di%v_bc_flag(sele) == V_BC_TYPE_PRESCRIBED_NORMAL_FLOW) then

                           if (di%subcy_opt_sat%have .and. di%subcy_opt_sat%consistent) then
                           
                              div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) + &
                                                             di%subcy_opt_sat%number * &
                                                             bc_sele_val(iloc) * &
                                                             detwei_bdy(ggi)                   
                           
                           else
                           
                              div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) + &
                                                             bc_sele_val(iloc) * &
                                                             detwei_bdy(ggi)                   
                           
                           end if
                           
                        else prescribed_normal_flow
                         
                           ! NOTE this v_over_relperm_dot_n does not contain the relperm terms
                           v_over_relperm_dot_n = dot_product(v_over_relperm_face_quad_bdy(:,ggi), normal_bdy(:,ggi))
                           
                           ! Adding this term includes the gradient of pressure, gradient of
                           ! capilliary pressure and the gravity term. 
                           if (di%subcy_opt_sat%have .and. di%subcy_opt_sat%consistent) then
                           
                              div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) + &
                                                             sum(di%cached_face_value%relperm_bdy(:,ggi,sele,p)) * &
                                                             v_over_relperm_dot_n * &
                                                             detwei_bdy(ggi)
                              
                           else
                              
                              div_tvel_rhs_local_bdy(iloc) = div_tvel_rhs_local_bdy(iloc) + &
                                                             di%cached_face_value%relperm_bdy(1,ggi,sele,p) * &
                                                             v_over_relperm_dot_n * &
                                                             detwei_bdy(ggi)
                           
                           end if
                           
                           ! Add two extra terms associated with weak pressure BC's
                                                      
                           ! Add two extra terms associated with weak pressure BC's
                           pressure_weak: if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then

                              ! Form the face value = detwei * (relperm*absperm/visc)
                              if (di%subcy_opt_sat%have .and. di%subcy_opt_sat%consistent) then
                              
                                 face_value = detwei_bdy(ggi) * &
                                              sum(di%cached_face_value%relperm_bdy(:,ggi,sele,p)) * &
                                              absperm_ele_bdy(1) / &
                                              visc_ele_bdy(1) 
                                 
                              else
                                 
                                 face_value = detwei_bdy(ggi) * &
                                              di%cached_face_value%relperm_bdy(1,ggi,sele,p) * &
                                              absperm_ele_bdy(1) / &
                                              visc_ele_bdy(1) 
                              
                              end if
                              
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

                           end if pressure_weak
                        
                        end if prescribed_normal_flow
                        
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

   subroutine darcy_impes_assemble_and_solve_phase_saturations(di, &
                                                               form_new_subcycle_relperm_face_values, &
                                                               have_dual)
      
      !!< Assemble and solve the saturations for each subcycle.
      !!< If form_new_subcycle_relperm_face_values is true
      !!< and there are subcycles then new relperm face values are 
      !!< formed. Phase 1 is special as it may be either solved 
      !!< prognostically or calculated diagnostically from the 
      !!< other phase saturation values that are solved first (for each subcycle).
      
      type(darcy_impes_type), intent(inout) :: di
      logical,                intent(in)    :: form_new_subcycle_relperm_face_values
      logical,                intent(in)    :: have_dual      

      ! local variables
      integer :: p, isub, isub_index
      real    :: alpha_start, alpha_end

      ewrite(1,*) 'Assemble and solve the saturations for each phase'
                  
      ! Solve the saturations for each subcycle via:            
      !  - Form the relperm face values for each phase if required
      !  - For each phase:
      !     - Form the lhs = cv_mass_sfield_mesh_with_porosity / dt
      !     - Form the rhs time = cv_mass_pressure_mesh_with_old_porosity * old_saturation_field / dt and add to rhs
      !     - Add the s_source to rhs, which may include the dual perm transmissibility term
      !     - Assemble and add the rhs adv to rhs 
      !     - Apply any strong dirichlet BCs
      !     - solve for latest saturations and copy to start subcycle saturations step
      
      ! Note: the porosity at the start and end of a subcycle time step
      ! are linearly interpolated values from the main time step start and end

      ! Deduce the subcycle time step size
      if (di%subcy_opt_sat%have) then

         di%dt_subcycle = di%dt/real(di%subcy_opt_sat%number)

      else 

         di%dt_subcycle = di%dt

      end if 

      ! Set the initial old saturation subcycle field
      do p = 1, di%number_phase
         call set(di%old_saturation_subcycle(p)%ptr, di%old_saturation(p)%ptr)
      end do

      ! Get the pressure BC - required if no v given and weak pressure dirichlet given for extra integrals
      call darcy_impes_get_entire_sfield_boundary_condition(di%pressure(1)%ptr, &
                                                            (/"weakdirichlet"/), &
                                                            di%bc_surface_mesh, &
                                                            di%pressure_bc_value, &
                                                            di%pressure_bc_flag)

      sub_loop: do isub = 1, di%subcy_opt_sat%number
         
         ewrite(1,*) 'Start subcycle ',isub
         
         ! form the start and end of subcycle dt
         ! porosity linear interpolents
         alpha_start = (isub - 1) / di%subcy_opt_sat%number
         alpha_end   = isub / di%subcy_opt_sat%number
         
         ! Determine the index for cached face values depending on subcycling
         if (di%subcy_opt_sat%have .and. di%subcy_opt_sat%consistent) then
            isub_index = isub
         else
            isub_index = 1
         end if            

         ! If this is not the first subcycle then form the relperms 
         ! for the start of this subcycle from the end of last subcycle sat
         ! and face values if required
         if (form_new_subcycle_relperm_face_values .and. (isub > 1)) then
            call darcy_impes_calculate_relperm_fields(di)
            
            call darcy_impes_calculate_relperm_isub_face_values(di, isub_index)
         end if
       
         phase_loop: do p = di%number_phase, 1, -1
            
            ewrite(1,*) 'Assemble and solve phase ',p
            
            ! If this is phase 1 and it is diagnostic then calculate it and exit loop
            if ((p == 1) .and. di%phase_one_saturation_diagnostic) then
               
               call darcy_impes_calculate_phase_one_saturation_diagnostic(di)

               ! Copy latest phase 1 solution to old subcycle field
               call set(di%old_saturation_subcycle(1)%ptr, di%saturation(1)%ptr)
               
               exit phase_loop
               
            end if
            
            ! Allocate and get the BC data. If the BC is time dependent then
            ! the end of overall time step values are used for all subcycles. 

            ! Get this phase v BC info - only for no_normal_flow and prescribed_normal_flow which is special as it is a scalar
            call darcy_impes_get_v_boundary_condition(di%darcy_velocity(p)%ptr, &
                                                      (/"prescribed_normal_flow", &
                                                        "no_normal_flow        "/), &
                                                      di%bc_surface_mesh, &
                                                      di%v_bc_value, &
                                                      di%v_bc_flag)
            
            call zero(di%lhs)
            call zero(di%rhs)
            
            call addto(di%lhs, di%cv_mass_pressure_mesh_with_porosity, scale = alpha_end / di%dt_subcycle)

            call addto(di%lhs, di%cv_mass_pressure_mesh_with_old_porosity, scale = (1.0 - alpha_end) / di%dt_subcycle)

            ! form the rhs time contribution and add
            call addto(di%rhs, di%cv_mass_pressure_mesh_with_porosity, scale = alpha_start / di%dt_subcycle)

            call addto(di%rhs, di%cv_mass_pressure_mesh_with_old_porosity, scale = (1.0 - alpha_start) / di%dt_subcycle)

            call scale(di%rhs, di%old_saturation_subcycle(p)%ptr)

            ! Add the saturation_source contribution
            call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_source, di%saturation_source(p)%ptr)
            
            ! Should this include porosity ?!?!
            call addto(di%rhs, di%cv_mass_pressure_mesh_with_source)
            
            ! Add the dual source term
            if (have_dual) then               
               call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_lambda_dual, di%transmissibility_lambda_dual(p)%ptr)

               call set(di%rhs_dual, di%pressure_other_porous_media(p)%ptr)
               call addto(di%rhs_dual, di%pressure(p)%ptr, scale = -1.0)            
               call scale(di%rhs_dual, di%cv_mass_pressure_mesh_with_lambda_dual)
               
               call addto(di%rhs, di%rhs_dual)
            end if
            
            ! assemble the rhs adv contribution and add
            call darcy_impes_assemble_saturation_rhs_adv(di, &
                                                         p, &
                                                         isub_index)

            ! apply strong dirichlet BC
            call apply_dirichlet_conditions(di%lhs, di%rhs, di%saturation(p)%ptr)

            ! Solve for the saturation
            di%saturation(p)%ptr%val = di%rhs%val / di%lhs%val
            
            ! Update the halos
            call halo_update(di%saturation(p)%ptr)
            
            ! Set the strong BC nodes to the values to be consistent
            call set_dirichlet_consistent(di%saturation(p)%ptr) 

            ! Copy latest solution to old subcycle
            call set(di%old_saturation_subcycle(p)%ptr, di%saturation(p)%ptr)
         
         end do phase_loop
         
      end do sub_loop
      
      !******************** 27 July 2013 LCai **********************!
      !Check wether there is MIM model, if true calculate the Mobile saturation
      if (di%MIM_options%have_MIM_phase) call darcy_trans_MIM_assemble_and_solve_mobile_saturation(di)
      	
      !******Finish********* 27 July 2013 LCai **********************!
      
      ewrite(1,*) 'Finished assemble and solve the saturations for each phase'
            
   end subroutine darcy_impes_assemble_and_solve_phase_saturations

! ----------------------------------------------------------------------------
                  
   subroutine darcy_impes_assemble_saturation_rhs_adv(di, &
                                                      p, &
                                                      isub_index)

      !!< Assemble the rhs advection contribtion for saturation phase p
      
      type(darcy_impes_type), intent(inout) :: di
      integer,                intent(in)    :: p      
      integer,                intent(in)    :: isub_index    

      ! local variables
      integer :: vele, iloc, oloc, jloc, face, gi, ggi, sele, dim
      real    :: face_value, v_over_relperm_dot_n
      real,    dimension(1)            :: absperm_ele, visc_ele
      real,    dimension(:,:), pointer :: grad_pressure_face_quad
      real,    dimension(:,:), pointer :: grav_ele
      real,    dimension(:,:), pointer :: x_ele
      real,    dimension(:,:), pointer :: normal
      real,    dimension(:),   pointer :: detwei
      real,    dimension(:),   pointer :: normgi
      logical, dimension(:),   pointer :: notvisited
      real,    dimension(:),   pointer :: s_rhs_local
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
      real,    dimension(:),   pointer :: s_rhs_local_bdy
      integer, dimension(:),   pointer :: p_nodes_bdy
      integer, parameter :: V_BC_TYPE_PRESCRIBED_NORMAL_FLOW = 1, V_BC_TYPE_NO_NORMAL_FLOW = 2
      integer, parameter :: PRESSURE_BC_TYPE_WEAKDIRICHLET = 1
      
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
      allocate(s_rhs_local_bdy(face_loc(di%pressure_mesh,1)))
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(p_nodes_bdy(face_loc(di%pressure_mesh,1)))      
      
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

                        face_value = detwei(ggi) * di%cached_face_value%relperm(isub_index,ggi,vele,p) * absperm_ele(1) * &
dot_product((grad_pressure_face_quad(:,ggi) - di%cached_face_value%den(ggi,vele,p) * grav_ele(:,1)), normgi)/ visc_ele(1)

                        ! Form the local rhs for iloc and opposing oloc with normal vector sign change
                        s_rhs_local(iloc) = s_rhs_local(iloc) + face_value

                        s_rhs_local(oloc) = s_rhs_local(oloc) - face_value

                     end if check_visited

                  end do quadrature_loop

               end if is_neigh

            end do face_loop

         end do nodal_loop_i

         ! Add volume element contribution to global rhs field
         call addto(di%rhs, p_nodes, s_rhs_local)

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

         if (di%v_bc_flag(sele) == V_BC_TYPE_PRESCRIBED_NORMAL_FLOW) then

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

                     prescribed_normal_flow: if (di%v_bc_flag(sele) == V_BC_TYPE_PRESCRIBED_NORMAL_FLOW) then

                        s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - &
                                                bc_sele_val(iloc) * &
                                                detwei_bdy(ggi)

                     else prescribed_normal_flow                           

                        ! NOTE this v_over_relperm_dot_n does not contain the relperm terms
                        v_over_relperm_dot_n = dot_product(v_over_relperm_face_quad_bdy(:,ggi), normal_bdy(:,ggi))

                        ! This term includes the gradient pressure, gradient capilliary pressure 
                        ! and the gravity term. 
                        s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - &
                                                di%cached_face_value%relperm_bdy(isub_index,ggi,sele,p) * &
                                                v_over_relperm_dot_n * &
                                                detwei_bdy(ggi)

                        ! Add two extra terms associated with weak pressure BC's                              
                        pressure_weak: if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then

                           ! Form the face value = detwei * (relperm*absperm/visc)
                           face_value = detwei_bdy(ggi) * &
                                        di%cached_face_value%relperm_bdy(isub_index,ggi,sele,p) * &
                                        absperm_ele_bdy(1) / &
                                        visc_ele_bdy(1) 

                           do jloc = 1, di%pressure_mesh%faces%shape%loc 

                              s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - &
                                                      inv_char_len_ele_bdy(jloc) * &
                                                      di%p_cvbdyshape%n(jloc,ggi) * &
                                                      bc_sele_val(jloc) * &
                                                      sum(normal_bdy(:,ggi)) * &
                                                      face_value * &
                                                      di%weak_pressure_bc_coeff

                           end do

                           do jloc = 1, di%pressure_mesh%faces%shape%loc 

                              s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) + &
                                                      inv_char_len_ele_bdy(jloc) * &
                                                      di%p_cvbdyshape%n(jloc,ggi) * &
                                                      pressure_ele_bdy(jloc) * &
                                                      sum(normal_bdy(:,ggi)) * &
                                                      face_value * &
                                                      di%weak_pressure_bc_coeff

                           end do

                        end if pressure_weak

                     end if prescribed_normal_flow

                  end do bc_quad_loop

               end if bc_neigh_if

            end do bc_face_loop

         end do bc_iloc_loop

         call addto(di%rhs, p_nodes_bdy, s_rhs_local_bdy)

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
      deallocate(s_rhs_local_bdy)
      deallocate(x_ele_bdy)
      deallocate(p_nodes_bdy)

   end subroutine darcy_impes_assemble_saturation_rhs_adv
   
! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve_generic_prog_sfields(di)
      
      !!< Assemble and solve generic prognostic scalar fields
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: f, p
      
      ewrite(1,*) 'Assemble and solve generic prognostic scalar fields'
      
      sfield_loop: do f = 1, size(di%generic_prog_sfield)
         
         p = di%generic_prog_sfield(f)%phase
         
         call darcy_impes_assemble_and_solve_generic_prog_sfield(di, f, p)
                  
      end do sfield_loop

      ewrite(1,*) 'Finished assemble and solve generic prognostic scalar fields'

   end subroutine darcy_impes_assemble_and_solve_generic_prog_sfields

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve_generic_prog_sfield(di, f, p)
      
      !!< Assemble and a solve generic prognostic scalar field given by index f and p
      
      type(darcy_impes_type), intent(inout) :: di
      integer,                intent(in)    :: f
      integer,                intent(in)    :: p
      
      ! local variables
      integer :: i, imf ! ***imf is the index of immobile prog sfield ** LCai **

      ! *** 16 Aug 2013 ***LCai*****************************************************
      character(len=FIELD_NAME_LEN) :: cp_imfield_name  !*** the name of the coupled immobile field with mobile field 
                                                   !*** 16 Aug 2013 *** LCai ***
      

      type(scalar_field) :: temp_MIM_src !*** the temperory term for source term which is '1/(theta_s+alpha*dt)'

      type(scalar_field) :: leach_src !**the source term for leaching chemical model**16 July 2014 Lcai

      ! ***Finish *** LCai *********************************************************


      ewrite(1,*) 'Assemble and solve sfield ',trim(di%generic_prog_sfield(f)%sfield%name),' of phase ',p

      ! Get this phase v BC info - only for no_normal_flow and prescribed_normal_flow which is special as it is a scalar
      call darcy_impes_get_v_boundary_condition(di%darcy_velocity(p)%ptr, &
                                                (/"prescribed_normal_flow", &
                                                  "no_normal_flow        "/), &
                                                di%bc_surface_mesh, &
                                                di%v_bc_value, &
                                                di%v_bc_flag)

      ! Get the pressure BC - required if no v given and weak pressure dirichlet given for extra integrals
      call darcy_impes_get_entire_sfield_boundary_condition(di%pressure(1)%ptr, &
                                                            (/"weakdirichlet"/), &
                                                            di%bc_surface_mesh, &
                                                            di%pressure_bc_value, &
                                                            di%pressure_bc_flag)

      ! Get the sfield BC
      call darcy_impes_get_entire_sfield_boundary_condition(di%generic_prog_sfield(f)%sfield, &
                                                            (/"weakdirichlet"/), &
                                                            di%bc_surface_mesh, &
                                                            di%sfield_bc_value, &
                                                            di%sfield_bc_flag)
            
      call zero(di%lhs)         
      call zero(di%matrix)
      call zero(di%rhs)
      call zero(di%rhs_full)
      
      ! **************** 27 July& 16 Aug 2013 LCai *******************!
      !check the saturation will be used to calculate the ADE, total saturation OR mobile saturation
      if (di%MIM_options%have_MIM(p)) then
          ewrite(1,*) 'Assemble and solve prog sfield, use the mobile saturation of phase', p
        di%sat_ADE => di%MIM_options%mobile_saturation(p)%ptr
        di%old_sat_ADE => di%MIM_options%old_mobile_saturation(p)%ptr
      else
        ewrite(1,*) 'Assemble and solve prog sfield, use the total saturation of phase', p
        di%sat_ADE => di%saturation(p)%ptr
        di%old_sat_ADE =>  di%old_saturation(p)%ptr
      end if
       ! *******Finish*** 27 July 2013 LCai *******************!

          
      
      ! Add porosity*saturation(absorption + 1/dt) to lhs 
      if (di%generic_prog_sfield(f)%have_abs) then
         
         call addto(di%lhs, di%generic_prog_sfield(f)%sfield_abs)         
         
      end if
      
      call addto(di%lhs, 1.0/di%dt)            
      call scale(di%lhs, di%cv_mass_pressure_mesh_with_porosity)
      call scale(di%lhs, di%sat_ADE)  ! *****27 July 2013 LCai***change the saturation into pointed one
            
      ! Add old_porosity*old_saturation*old_sfield/dt to rhs      
      call addto(di%rhs, 1.0/di%dt)
      call scale(di%rhs, di%cv_mass_pressure_mesh_with_old_porosity)
      call scale(di%rhs, di%old_sat_ADE) ! *****27 July 2013 LCai***change the old saturation into pointed one
      call scale(di%rhs, di%generic_prog_sfield(f)%old_sfield)
      
      ! Add source to rhs   
      if (di%generic_prog_sfield(f)%have_src) then
         
         call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_source, di%generic_prog_sfield(f)%sfield_src)
         
         call addto(di%rhs, di%cv_mass_pressure_mesh_with_source)
         
      end if

      !**************16 July 2014 lcai**********Leaching chemical model*************!
      !Add leaching chemical source term to rhs
      if (di%lc%have_leach_chem_model) then
         call add_leach_chemical_prog_src_to_rhs(di,f)
      end if
      
      !*************end 16 July 2014 lcai**********Leaching chemical model**********!
      

      ! ************** 15 & 22 Aug 2013 LCai ********************************!
      !If the immobile prognostic field exist, solve the source term of Mobile-immobile model implicitly
      if (di%generic_prog_sfield(f)%have_MIM_source) then
         
         call allocate(temp_MIM_src, di%pressure_mesh)
         call zero(temp_MIM_src)

         call zero(di%MIM_options%MIM_src)
         call zero(di%MIM_options%MIM_src_s)
         

         !Addto the lhs matrix
         !calculate 1/(theta_s+alpha*dt)
         call set(di%MIM_options%MIM_src, di%MIM_options%immobile_saturation(p)%ptr)
         
         !check wether to scale with the porosity as a constant of a scalar field
         if (di%prt_is_constant) then
          call scale(di%MIM_options%MIM_src, di%porosity_cnt(1))
         else
          call scale(di%MIM_options%MIM_src, di%porosity_pmesh)
         end if


         call addto(di%MIM_options%MIM_src, di%MIM_options%mass_trans_coef(p)%ptr, scale=di%dt)
         call invert(di%MIM_options%MIM_src, temp_MIM_src) !temp_MIM_src is '1/(theta_s+alpha*dt)'

         call set(di%MIM_options%MIM_src, temp_MIM_src)
         call scale(di%MIM_options%MIM_src, di%MIM_options%mass_trans_coef(p)%ptr) 
         call scale(di%MIM_options%MIM_src, di%MIM_options%mass_trans_coef(p)%ptr)! repeated to scale with alpha**2
         call scale(di%MIM_options%MIM_src, di%dt)
         call addto(di%MIM_options%MIM_src, di%MIM_options%mass_trans_coef(p)%ptr, scale=-1.0)
           
         call compute_cv_mass(di%positions, di%MIM_options%MIM_src_s, di%MIM_options%MIM_src)

         call addto(di%lhs, di%MIM_options%MIM_src_s, scale=-1.0)

         !Addto the rhs matrix
         !before start to compute rhs, zero the field to be used
         call zero(di%MIM_options%MIM_src)

         !extract the idex of coupled immobile prog sfield of this generic prog sfield
         cp_imfield_name = di%generic_prog_sfield(f)%source_name
         
         imf= 0

         do imf= 1, size(di%MIM_options%immobile_prog_sfield)
           if (di%MIM_options%immobile_prog_sfield(imf)%phase == p) then
             if (trim(di%MIM_options%immobile_prog_sfield(imf)%sfield%name) == trim(cp_imfield_name)) exit  
           end if
         end do

         !check the selected immobile sfield is correct
         if (.not.( trim(di%MIM_options%immobile_prog_sfield(imf)%source_name) == trim(di%generic_prog_sfield(f)%sfield%name))) then
         FLExit('The source terms of the mobile-immobile concentration should coincide')
         end if

         !Add the source term with old immobile concentration 
         call set(di%MIM_options%MIM_src,temp_MIM_src)
         call scale(di%MIM_options%MIM_src, di%MIM_options%mass_trans_coef(p)%ptr)
         call scale(di%MIM_options%MIM_src, di%MIM_options%old_immobile_saturation(p)%ptr)
         call scale(di%MIM_options%MIM_src, di%MIM_options%immobile_prog_sfield(imf)%old_sfield)
         call scale(di%MIM_options%MIM_src, di%cv_mass_pressure_mesh_with_old_porosity) ! this has already inlcuded the cv pmesh
         call addto (di%rhs, di%MIM_options%MIM_src)

      end if
      ! *************Finish *** LCai **********************************!
      
      
      ! Add diagonal lhs to matrix
      call addto_diag(di%matrix, di%lhs)
      
      ! Add implicit low order advection and diffusion terms to matrix and rhs
      if (di%generic_prog_sfield(f)%have_adv .or. di%generic_prog_sfield(f)%have_diff) then
         
         call darcy_impes_assemble_generic_prog_sfield_adv_diff(di, f, p)
      
      end if
      
      ! Solve for each face value iteration (default is 1)
      do i = 1, di%generic_prog_sfield(f)%sfield_cv_options%number_face_value_iteration
         
         call set(di%rhs_full, di%rhs)
         
         ! Assmemble the high resolution rhs and to to rhs_full
         if (di%generic_prog_sfield(f)%have_adv .and. &
             di%generic_prog_sfield(f)%sfield_cv_options%facevalue == DARCY_IMPES_CV_FACEVALUE_FINITEELEMENT) then
            
            call zero(di%rhs_high_resolution)
            
            call darcy_impes_assemble_generic_prog_sfield_rhs_high_resolution(di, f, p)
            
            call addto(di%rhs_full, di%rhs_high_resolution)
                        
         end if
         
         ! Apply any strong dirichlet BC's
         call apply_dirichlet_conditions(di%matrix, di%rhs_full, di%generic_prog_sfield(f)%sfield)

         ! Solve the sfield
         call petsc_solve(di%generic_prog_sfield(f)%sfield, di%matrix, di%rhs_full, di%state(p))

         ! Set the strong BC nodes to the values to be consistent
         call set_dirichlet_consistent(di%generic_prog_sfield(f)%sfield) 
      
      end do
      
      ! **************** 27 July 2013 LCai *****************!
      nullify(di%sat_ADE) 
      nullify(di%old_sat_ADE)
       ! ***** Finish*** 27 July 2013 LCai *****************!     
      
      ewrite(1,*) 'Finished assemble and solve sfield ',trim(di%generic_prog_sfield(f)%sfield%name),' of phase ',p

      
      ! ******* 16 Aug 2013 *** LCai **********************
      !Solve the immobile prog sfield if it exist
      if (di%generic_prog_sfield(f)%have_MIM_source) then 
          
          call darcy_trans_assemble_and_solve_immobile_prog_sfield(di, p, f, imf, temp_MIM_src) 

          call deallocate(temp_MIM_src)

      end if

      !****** Finished ******* LCai **********************

   end subroutine darcy_impes_assemble_and_solve_generic_prog_sfield

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_generic_prog_sfield_adv_diff(di, f, p)
      
      !!< Assemble the advection and diffusion terms for the 
      !!< generic prognostic scalar field given by index f and p
      
      type(darcy_impes_type), intent(inout) :: di
      integer,                intent(in)    :: f
      integer,                intent(in)    :: p
      
      ! local variables
      logical :: inflow
      integer :: vele, iloc, oloc, jloc, face, gi, ggi, sele
      real    :: face_value, v_dot_n, income, sat_face_value
      real,    dimension(1)            :: absperm_ele, visc_ele, porosity_ele
      real,    dimension(:,:),   pointer :: grad_pressure_face_quad
      real,    dimension(:,:),   pointer :: grav_ele
      real,    dimension(:,:,:), pointer :: diff_ele
      real,    dimension(:),     pointer :: sat_ele            
      real,    dimension(:,:),   pointer :: x_ele
      real,    dimension(:,:,:), pointer :: p_dshape
      real,    dimension(:,:),   pointer :: normal
      real,    dimension(:),     pointer :: detwei
      real,    dimension(:),     pointer :: normgi
      logical, dimension(:),     pointer :: notvisited
      real,    dimension(:,:),   pointer :: matrix_local
      real,    dimension(:,:),   pointer :: x_face_quad
      integer, dimension(:),     pointer :: x_pmesh_nodes
      integer, dimension(:),     pointer :: p_nodes      
      real,    dimension(1)              :: absperm_ele_bdy, visc_ele_bdy
      real,    dimension(:,:),   pointer :: grav_ele_bdy
      real,    dimension(:),     pointer :: ghost_sfield_ele_bdy      
      real,    dimension(:,:),   pointer :: grad_pressure_face_quad_bdy
      real,    dimension(:,:),   pointer :: v_over_relperm_face_quad_bdy
      real,    dimension(:),     pointer :: bc_sele_val
      real,    dimension(:),     pointer :: pressure_ele_bdy
      real,    dimension(:),     pointer :: inv_char_len_ele_bdy
      real,    dimension(:,:),   pointer :: normal_bdy
      real,    dimension(:),     pointer :: detwei_bdy
      real,    dimension(:,:),   pointer :: x_ele_bdy
      real,    dimension(:),     pointer :: rhs_local_bdy
      real,    dimension(:),     pointer :: matrix_local_bdy
      integer, dimension(:),     pointer :: p_nodes_bdy
      integer, parameter :: V_BC_TYPE_PRESCRIBED_NORMAL_FLOW = 1, V_BC_TYPE_NO_NORMAL_FLOW = 2
      integer, parameter :: PRESSURE_BC_TYPE_WEAKDIRICHLET = 1
      integer, parameter :: SFIELD_BC_TYPE_WEAKDIRICHLET   = 1
      
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
      allocate(matrix_local(ele_loc(di%pressure_mesh,1),ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(grav_ele(di%ndim,1))
      allocate(diff_ele(di%ndim,di%ndim,1))
      allocate(sat_ele(ele_loc(di%positions,1)))      

      allocate(grav_ele_bdy(di%ndim,1))
      allocate(ghost_sfield_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(grad_pressure_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      allocate(v_over_relperm_face_quad_bdy(di%ndim, di%p_cvbdyshape%ngi))
      allocate(bc_sele_val(face_loc(di%pressure_mesh,1)))
      allocate(pressure_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(inv_char_len_ele_bdy(face_loc(di%pressure_mesh,1)))
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(detwei_bdy(di%x_cvbdyshape%ngi))
         allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      end if
      allocate(rhs_local_bdy(face_loc(di%pressure_mesh,1)))
      allocate(matrix_local_bdy(face_loc(di%pressure_mesh,1)))
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(p_nodes_bdy(face_loc(di%pressure_mesh,1)))      
      
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
                  
         if (di%generic_prog_sfield(f)%have_diff) then
       
            ! get the sfield diffusivity 
            diff_ele = ele_val(di%generic_prog_sfield(f)%sfield_diff,vele)
            
            ! get the element values for porosity
            porosity_ele = ele_val(di%porosity, vele)
            
            ! get the element values for saturation
            sat_ele = ele_val(di%sat_ADE, vele)  ! *** 27 July 2013 LCai **change the saturation***!
            
            ! obtain the derivative of the pressure mesh shape function at the CV face quadrature points
            if (di%cached_face_value%cached_p_dshape) then

               p_dshape => di%cached_face_value%p_dshape(:,:,:,vele)

            else

               call transform_to_physical(di%positions, vele, x_shape = di%x_cvshape_full, &
                                          shape = di%p_cvshape_full, dshape = p_dshape)

            end if
         
         end if
         
         ! Initialise array for the quadrature points of this 
         ! element for whether it has already been visited
         notvisited = .true.

         ! Initialise the local matrix to assemble for this element
         matrix_local = 0.0
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

                        v_dot_n = - di%cached_face_value%relperm(1,ggi,vele,p) * absperm_ele(1) * &
dot_product((grad_pressure_face_quad(:,ggi) - di%cached_face_value%den(ggi,vele,p) * grav_ele(:,1)), normgi)/ visc_ele(1)

                        inflow = (v_dot_n<=0.0)

                        income = merge(1.0,0.0,inflow)                        
                        
                        if (di%generic_prog_sfield(f)%have_adv) then

                           matrix_local(iloc, oloc) = matrix_local(iloc, oloc) + &
                                                      detwei(ggi) * &
                                                      v_dot_n * &
                                                      income

                           matrix_local(oloc, iloc) = matrix_local(oloc, iloc) + &
                                                      detwei(ggi) * &
                                                      (-v_dot_n) * &
                                                      (1.0-income)

                           matrix_local(iloc, iloc) = matrix_local(iloc, iloc) + &
                                                      detwei(ggi) * &
                                                      v_dot_n * &
                                                      (1.0-income)

                           matrix_local(oloc, oloc) = matrix_local(oloc, oloc) + &
                                                      detwei(ggi) * &
                                                      (-v_dot_n) * &
                                                      income
                        
                        end if
                        
                        if (di%generic_prog_sfield(f)%have_diff) then
                           
                           ! Find the diffusion term saturation face value (taking upwind)
                           sat_face_value = income*sat_ele(oloc) + (1.0-income)*sat_ele(iloc)
                           
                           do jloc=1, di%pressure_mesh%shape%loc
                             
                             matrix_local(iloc,jloc) = matrix_local(iloc,jloc) - &
                                                       sum(matmul(diff_ele(:,:,1), p_dshape(jloc, ggi, :))*normgi, 1) * &
                                                       porosity_ele(1) * &
                                                       sat_face_value * &
                                                       detwei(ggi)

                             matrix_local(oloc, jloc) = matrix_local(oloc,jloc) - &
                                                        sum(matmul(diff_ele(:,:,1), p_dshape(jloc, ggi, :))*(-normgi), 1) * &
                                                        porosity_ele(1) * &
                                                        sat_face_value * &
                                                        detwei(ggi)
                           
                           end do                           
                           
                        end if
                        
                     end if check_visited

                  end do quadrature_loop

               end if is_neigh

            end do face_loop

         end do nodal_loop_i

         call addto(di%matrix, p_nodes, p_nodes, matrix_local)

      end do vol_element_loop

      ! Add BC integrals associated with an advection term (diffusive flux is assumed zero on all boundary)
      
      have_adv_so_add_bc: if (di%generic_prog_sfield(f)%have_adv) then
      
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

            if (di%v_bc_flag(sele) == V_BC_TYPE_PRESCRIBED_NORMAL_FLOW) then

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

            end if

            if (di%sfield_bc_flag(sele) == SFIELD_BC_TYPE_WEAKDIRICHLET) then
               ghost_sfield_ele_bdy = ele_val(di%sfield_bc_value, sele)
            else
               ghost_sfield_ele_bdy = face_val(di%generic_prog_sfield(f)%sfield, sele)
            end if

            ! Initialise the local arrays to assemble for this element
            rhs_local_bdy    = 0.0
            matrix_local_bdy = 0.0

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi

                        prescribed_normal_flow: if (di%v_bc_flag(sele) == V_BC_TYPE_PRESCRIBED_NORMAL_FLOW) then

                           if (bc_sele_val(iloc) <= 0) then

                              income = 1.0

                           else

                              income = 0.0

                           end if 

                           rhs_local_bdy(iloc) = rhs_local_bdy(iloc) - &
                                                 bc_sele_val(iloc) * &
                                                 ghost_sfield_ele_bdy(iloc) * &
                                                 detwei_bdy(ggi) * &
                                                 income

                           matrix_local_bdy(iloc) = matrix_local_bdy(iloc) + &
                                                    bc_sele_val(iloc) * &
                                                    detwei_bdy(ggi) * &
                                                    (1.0 - income)

                        else prescribed_normal_flow                          

                           v_dot_n = - di%cached_face_value%relperm_bdy(1,ggi,sele,p) * absperm_ele_bdy(1) * &
dot_product((grad_pressure_face_quad_bdy(:,ggi) - di%cached_face_value%den_bdy(ggi,sele,p) * grav_ele_bdy(:,1)), normal_bdy(:,ggi))/ visc_ele_bdy(1)

                           if (v_dot_n <= 0) then

                              income = 1.0

                           else

                              income = 0.0

                           end if                         

                           rhs_local_bdy(iloc) = rhs_local_bdy(iloc) - &
                                                 v_dot_n * &
                                                 ghost_sfield_ele_bdy(iloc) * &
                                                 detwei_bdy(ggi) * &
                                                 income

                           matrix_local_bdy(iloc) = matrix_local_bdy(iloc) + &
                                                    v_dot_n * &
                                                    detwei_bdy(ggi) * &
                                                    (1.0 - income)

                           ! Add two extra terms associated with weak pressure BC's                              
                           pressure_weak: if (di%pressure_bc_flag(sele) == PRESSURE_BC_TYPE_WEAKDIRICHLET) then

                              ! Form the face value = detwei * (relperm*absperm/visc)
                              face_value = detwei_bdy(ggi) * &
                                           di%cached_face_value%relperm_bdy(1,ggi,sele,p) * &
                                           absperm_ele_bdy(1) / &
                                           visc_ele_bdy(1) 

                              do jloc = 1, di%pressure_mesh%faces%shape%loc 

                                 rhs_local_bdy(iloc) = rhs_local_bdy(iloc) - &
                                                       inv_char_len_ele_bdy(jloc) * &
                                                       di%p_cvbdyshape%n(jloc,ggi) * &
                                                       bc_sele_val(jloc) * &
                                                       sum(normal_bdy(:,ggi)) * &
                                                       face_value * &
                                                       ghost_sfield_ele_bdy(iloc) * &
                                                       di%weak_pressure_bc_coeff * &
                                                       income

                                 matrix_local_bdy(iloc) = matrix_local_bdy(iloc) - &
                                                          inv_char_len_ele_bdy(jloc) * &
                                                          di%p_cvbdyshape%n(jloc,ggi) * &
                                                          bc_sele_val(jloc) * &
                                                          sum(normal_bdy(:,ggi)) * &
                                                          face_value * &
                                                          di%weak_pressure_bc_coeff * &
                                                          (1.0 - income)

                              end do

                              do jloc = 1, di%pressure_mesh%faces%shape%loc 

                                 rhs_local_bdy(iloc) = rhs_local_bdy(iloc) + &
                                                       inv_char_len_ele_bdy(jloc) * &
                                                       di%p_cvbdyshape%n(jloc,ggi) * &
                                                       pressure_ele_bdy(jloc) * &
                                                       sum(normal_bdy(:,ggi)) * &
                                                       face_value * &
                                                       ghost_sfield_ele_bdy(iloc) * &
                                                       di%weak_pressure_bc_coeff * &
                                                       income

                                 matrix_local_bdy(iloc) = matrix_local_bdy(iloc) + &
                                                          inv_char_len_ele_bdy(jloc) * &
                                                          di%p_cvbdyshape%n(jloc,ggi) * &
                                                          pressure_ele_bdy(jloc) * &
                                                          sum(normal_bdy(:,ggi)) * &
                                                          face_value * &
                                                          di%weak_pressure_bc_coeff * &
                                                          (1.0 - income)

                              end do

                           end if pressure_weak

                        end if prescribed_normal_flow

                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop

            call addto(di%rhs, p_nodes_bdy, rhs_local_bdy)

            call addto_diag(di%matrix, p_nodes_bdy, matrix_local_bdy)

         end do sele_loop
      
      end if have_adv_so_add_bc
      
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
      deallocate(matrix_local)
      deallocate(x_face_quad)
      deallocate(grad_pressure_face_quad)
      deallocate(grav_ele)
      deallocate(diff_ele)
      deallocate(sat_ele)      

      deallocate(grav_ele_bdy)
      deallocate(ghost_sfield_ele_bdy)
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
      deallocate(rhs_local_bdy)
      deallocate(matrix_local_bdy)
      deallocate(x_ele_bdy)
      deallocate(p_nodes_bdy)      
      
   end subroutine darcy_impes_assemble_generic_prog_sfield_adv_diff

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_generic_prog_sfield_rhs_high_resolution(di, f, p)
      
      !!< Assemble the advection high resolution terms for the 
      !!< generic prognostic scalar field given by index f and p
      
      type(darcy_impes_type), intent(inout) :: di
      integer,                intent(in)    :: f
      integer,                intent(in)    :: p
      
      ! local variables
      logical :: inflow
      integer :: vele, iloc, oloc, face, gi, ggi, s_upwind_pos
      real    :: v_dot_n, income, sfield_face_value, sfield_low_face_value, sfield_high_face_value
      real,    dimension(1)            :: absperm_ele, visc_ele
      real,    dimension(:,:),   pointer :: grad_pressure_face_quad
      real,    dimension(:,:),   pointer :: grav_ele
      real,    dimension(:),     pointer :: sfield_ele            
      real,    dimension(:,:),   pointer :: x_ele
      real,    dimension(:,:,:), pointer :: p_dshape
      real,    dimension(:,:),   pointer :: normal
      real,    dimension(:),     pointer :: detwei
      real,    dimension(:),     pointer :: normgi
      logical, dimension(:),     pointer :: notvisited
      real,    dimension(:),     pointer :: rhs_local
      real,    dimension(:,:),   pointer :: x_face_quad
      integer, dimension(:),     pointer :: x_pmesh_nodes
      integer, dimension(:),     pointer :: p_nodes      
      integer, dimension(:),     pointer :: s_upwind_nodes
      
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
      allocate(rhs_local(ele_loc(di%pressure_mesh,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(grav_ele(di%ndim,1))
      allocate(sfield_ele(ele_loc(di%positions,1)))      
      
      if(di%generic_prog_sfield(f)%sfield_cv_options%limit_facevalue) then

        ! NOTE only new upwind values are required but this procedure 
        !      expects latest and old. So pass in new twice - not optimal
        call find_upwind_values(di%state, &
                                di%positions_pressure_mesh, &
                                di%generic_prog_sfield(f)%sfield, &
                                di%sfield_upwind, &
                                di%generic_prog_sfield(f)%sfield, &
                                di%sfield_upwind, &
                                option_path = trim(di%generic_prog_sfield(f)%sfield%option_path))

      end if
      
      ! Initialise optimisation flag used in finding upwind value in high resolution schemes
      s_upwind_pos = 0
      
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

         ! get the sfield value for this element
         sfield_ele = ele_val(di%generic_prog_sfield(f)%sfield, vele)

         ! Determine the node numbers to use to determine the sfield upwind values
         if((di%generic_prog_sfield(f)%sfield_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
            (di%generic_prog_sfield(f)%sfield_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

            s_upwind_nodes => x_pmesh_nodes

         else

            s_upwind_nodes => p_nodes

         end if
                          
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
                  
         ! Initialise array for the quadrature points of this 
         ! element for whether it has already been visited
         notvisited = .true.

         ! Initialise the local rhs to assemble for this element
         rhs_local = 0.0

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

                        v_dot_n = - di%cached_face_value%relperm(1,ggi,vele,p) * absperm_ele(1) * &
dot_product((grad_pressure_face_quad(:,ggi) - di%cached_face_value%den(ggi,vele,p) * grav_ele(:,1)), normgi)/ visc_ele(1)

                        inflow = (v_dot_n<=0.0)

                        income = merge(1.0,0.0,inflow)                        

                        ! form the high resolution face value
                        call darcy_impes_evaluate_sfield_face_val(sfield_high_face_value, &
                                                                  iloc, &
                                                                  oloc, &
                                                                  ggi, &
                                                                  s_upwind_nodes, &
                                                                  di%p_cvshape, &
                                                                  sfield_ele, &
                                                                  di%sfield_upwind, &
                                                                  inflow, &
                                                                  income, &
                                                                  di%generic_prog_sfield(f)%sfield_cv_options, &
                                                                  save_pos = s_upwind_pos)                        
                                                
                        ! form and add the lower order iterated term to rhs
                        sfield_low_face_value = income*sfield_ele(oloc) + (1.0-income)*sfield_ele(iloc)
                        
                        ! form the rhs face value and add to the rhs field
                        sfield_face_value = sfield_low_face_value - sfield_high_face_value
                        
                        rhs_local(iloc) = rhs_local(iloc) + &
                                          sfield_face_value * &
                                          v_dot_n * &
                                          detwei(ggi)

                        rhs_local(oloc) = rhs_local(oloc)  - &
                                          sfield_face_value * &
                                          v_dot_n * &
                                          detwei(ggi)
                                                                       
                     end if check_visited

                  end do quadrature_loop

               end if is_neigh

            end do face_loop

         end do nodal_loop_i

         call addto(di%rhs_high_resolution, p_nodes, rhs_local)

      end do vol_element_loop
      
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
      deallocate(rhs_local)
      deallocate(x_face_quad)
      deallocate(grad_pressure_face_quad)
      deallocate(grav_ele)
      deallocate(sfield_ele)      
      
   end subroutine darcy_impes_assemble_generic_prog_sfield_rhs_high_resolution

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
      type(vector_field) :: aux_vfield
      logical :: inflow
      integer :: vele, p, dim, iloc, oloc, face, gi, ggi, sele
      real    :: income, v_dot_n
      real,    dimension(1)              :: visc_ele, absperm_ele
      real,    dimension(:,:),   pointer :: grad_pressure_face_quad
      real,    dimension(:,:),   pointer :: v_local
      real,    dimension(:),     pointer :: m_local
      real,    dimension(:,:),   pointer :: grav_ele
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
      real,    dimension(:,:),   pointer :: v_local_bdy
      real,    dimension(:),     pointer :: m_local_bdy
      integer, dimension(:),     pointer :: p_nodes_bdy
      real,    dimension(:),     pointer :: bc_sele_val
      integer, parameter :: V_BC_TYPE_PRESCRIBED_NORMAL_FLOW = 1, V_BC_TYPE_NO_NORMAL_FLOW = 2
            
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
      
      allocate(v_local(di%ndim,ele_loc(di%pressure_mesh,1)))
      allocate(m_local(ele_loc(di%pressure_mesh,1)))
      allocate(v_local_bdy(di%ndim,face_loc(di%pressure_mesh,1)))
      allocate(m_local_bdy(face_loc(di%pressure_mesh,1)))

      phase_loop: do p = 1, di%number_phase
         
         call zero(di%cfl(p)%ptr)

         call zero(di%darcy_velocity(p)%ptr)
         
         call zero(di%mobility(p)%ptr)
                               
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

            ! Initialise the local assemble arrays for this element
            cfl_rhs_local = 0.0
            v_local       = 0.0
            m_local       = 0.0
                        
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

                           v_dot_n = &
- di%cached_face_value%relperm(1,ggi,vele,p) * absperm_ele(1) * &
dot_product((grad_pressure_face_quad(:,ggi) - di%cached_face_value%den(ggi,vele,p) * grav_ele(:,1)), normgi) / &
visc_ele(1)

                           inflow = (v_dot_n<=0.0)

                           income = merge(1.0,0.0,inflow)

                           cfl_rhs_local(iloc) = cfl_rhs_local(iloc) + &
                                                 abs(v_dot_n) * &
                                                 detwei(ggi) * &
                                                 (1.0 - income)

                           cfl_rhs_local(oloc) = cfl_rhs_local(oloc) + &
                                                 abs(v_dot_n) * &
                                                 detwei(ggi) * &
                                                 income

                           do dim = 1, di%ndim

                              v_local(dim,iloc) = v_local(dim,iloc) - &
detwei(ggi) * &
di%cached_face_value%relperm(1,ggi,vele,p) * &
absperm_ele(1) * &
(grad_pressure_face_quad(dim,ggi) - di%cached_face_value%den(ggi,vele,p) * grav_ele(dim,1)) / &
visc_ele(1)

                              v_local(dim,oloc) = v_local(dim,oloc) - &
detwei(ggi) * &
di%cached_face_value%relperm(1,ggi,vele,p) * &
absperm_ele(1) * &
(grad_pressure_face_quad(dim,ggi) - di%cached_face_value%den(ggi,vele,p) * grav_ele(dim,1)) / &
visc_ele(1)
                           
                           end do
                           
                           m_local(iloc) = m_local(iloc) + &
                                           detwei(ggi) * &
                                           di%cached_face_value%relperm(1,ggi,vele,p) / &
                                           visc_ele(1)

                           m_local(oloc) = m_local(oloc) + &
                                           detwei(ggi) * &
                                           di%cached_face_value%relperm(1,ggi,vele,p) / &
                                           visc_ele(1)

                        end if check_visited

                     end do quadrature_loop

                  end if is_neigh

               end do face_loop

            end do nodal_loop_i
            
            call addto(di%cfl(p)%ptr, p_nodes, cfl_rhs_local)
            
            call addto(di%darcy_velocity(p)%ptr, p_nodes, v_local)
            
            call addto(di%mobility(p)%ptr, p_nodes, m_local)

         end do vele_loop

         ! Get this phase v BC info - only for no_normal_flow and
         ! prescribed_normal_flow which is special as it is a scalar
         call darcy_impes_get_v_boundary_condition(&
              di%darcy_velocity(p)%ptr, &
              (/"prescribed_normal_flow", "no_normal_flow        "/), &
              di%bc_surface_mesh, di%v_bc_value, di%v_bc_flag)
         
         sele_loop: do sele = 1, di%number_sele
            
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
            
            if (di%v_bc_flag(sele) == V_BC_TYPE_PRESCRIBED_NORMAL_FLOW) then
            
               bc_sele_val = ele_val(di%v_bc_value, sele)
            
            end if
            
            visc_ele_bdy                = face_val(di%viscosity(p)%ptr, sele)
            absperm_ele_bdy             = face_val(di%absolute_permeability, sele)
            grad_pressure_face_quad_bdy = face_val_at_quad(di%gradient_pressure(p)%ptr, sele, di%gradp_cvbdyshape)         
            grav_ele_bdy                = face_val(di%gravity, sele)                         
            
            ! Initialise the local assemble arrays for this element
            cfl_rhs_local_bdy = 0.0
            v_local_bdy       = 0.0
            m_local_bdy       = 0.0

            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi
                        
                        if (di%v_bc_flag(sele) == V_BC_TYPE_NO_NORMAL_FLOW) then
                        
                           v_dot_n = 0.0
                        
                        else if (di%v_bc_flag(sele) == V_BC_TYPE_PRESCRIBED_NORMAL_FLOW) then
                        
                           v_dot_n = bc_sele_val(iloc)
                                                      
                        else

                           v_dot_n = &
- di%cached_face_value%relperm_bdy(1,ggi,sele,p) * absperm_ele_bdy(1) * &
dot_product((grad_pressure_face_quad_bdy(:,ggi) - di%cached_face_value%den_bdy(ggi,sele,p) * grav_ele_bdy(:,1)), normal_bdy(:,ggi)) / &
visc_ele_bdy(1)
                                                                              
                        end if
                        
                        if (v_dot_n > 0.0) then
                           income = 0.0
                        else
                           income = 1.0
                        end if
                        
                        cfl_rhs_local_bdy(iloc) = cfl_rhs_local_bdy(iloc) + &
                                                  abs(v_dot_n) * &
                                                  detwei_bdy(ggi) * &
                                                  (1.0 - income)                     

                        do dim = 1, di%ndim

                           v_local_bdy(dim,iloc) = v_local_bdy(dim,iloc) - &
detwei_bdy(ggi) * &
di%cached_face_value%relperm_bdy(1,ggi,sele,p) * &
absperm_ele_bdy(1) * &
(grad_pressure_face_quad_bdy(dim,ggi) - di%cached_face_value%den_bdy(ggi,sele,p) * grav_ele_bdy(dim,1)) / &
visc_ele_bdy(1)
                           
                        end do

                        m_local_bdy(iloc) = m_local_bdy(iloc) + &
                                            detwei_bdy(ggi) * &
                                            di%cached_face_value%relperm_bdy(1,ggi,sele,p) / &
                                            visc_ele_bdy(1)

                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop

            call addto(di%cfl(p)%ptr, p_nodes_bdy, cfl_rhs_local_bdy)
            
            call addto(di%darcy_velocity(p)%ptr, p_nodes_bdy, v_local_bdy)
            
            call addto(di%mobility(p)%ptr, p_nodes_bdy, m_local_bdy)

         end do sele_loop
                  
         di%cfl(p)%ptr%val = di%cfl(p)%ptr%val * di%dt / di%cv_mass_pressure_mesh_with_porosity%val
         
         call scale(di%darcy_velocity(p)%ptr, di%inverse_cv_sa_pressure_mesh)
         
         call scale(di%mobility(p)%ptr, di%inverse_cv_sa_pressure_mesh)
         
         ewrite_minmax(di%cfl(p)%ptr)

         ewrite_minmax(di%darcy_velocity(p)%ptr)

         ewrite_minmax(di%mobility(p)%ptr)

      end do phase_loop      
            
      ewrite(1,*) 'Calculate TotalDarcyVelocity and TotalMobility'
      
      call set(di%total_darcy_velocity, di%darcy_velocity(1)%ptr)
      call set(di%total_mobility, di%mobility(1)%ptr)
      
      do p = 2, di%number_phase
         
         call addto(di%total_darcy_velocity, di%darcy_velocity(p)%ptr)
         call addto(di%total_mobility, di%mobility(p)%ptr)
         
      end do

      ewrite_minmax(di%total_darcy_velocity)
      ewrite_minmax(di%total_mobility)
      
      ewrite(1,*) 'Calculate BulkDarcyVelocity'
      
      call allocate(aux_vfield, di%ndim, di%pressure_mesh)
      
      call zero(di%bulk_darcy_velocity)
      
      do p = 1, di%number_phase
      
         call set(aux_vfield, di%darcy_velocity(p)%ptr)
         call scale(aux_vfield, di%saturation(p)%ptr)
         
         call addto(di%bulk_darcy_velocity, aux_vfield)
         
      end do
      
      call deallocate(aux_vfield)
      
      ewrite_minmax(di%bulk_darcy_velocity)
      
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

      deallocate(v_local)
      deallocate(m_local)
      deallocate(v_local_bdy)
      deallocate(m_local_bdy)
             
   end subroutine darcy_impes_calculate_vel_mob_ff_and_cfl_fields

! ----------------------------------------------------------------------------

   subroutine darcy_impes_calculate_inverse_cv_sa(di)
      
      !!< Calculate the inverse of the surface area of the CV's on the pressure mesh
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: vele, iloc, oloc, face, gi, ggi, sele
      real,    dimension(:),     pointer :: sa_local
      real,    dimension(:,:),   pointer :: x_ele
      real,    dimension(:,:),   pointer :: normal
      real,    dimension(:),     pointer :: detwei
      logical, dimension(:),     pointer :: notvisited
      integer, dimension(:),     pointer :: p_nodes      
      real,    dimension(:,:),   pointer :: normal_bdy
      real,    dimension(:),     pointer :: detwei_bdy
      real,    dimension(:,:),   pointer :: x_ele_bdy
      real,    dimension(:),     pointer :: sa_local_bdy
      integer, dimension(:),     pointer :: p_nodes_bdy
            
      ewrite(1,*) 'Calculate the inverse of the surface area of CV on the pressure mesh'

      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(normal(di%ndim,di%x_cvshape%ngi))
         allocate(detwei(di%x_cvshape%ngi))      
      end if
      allocate(notvisited(di%x_cvshape%ngi))
      
      if (.not. di%cached_face_value%cached_detwei_normal) then
         allocate(detwei_bdy(di%x_cvbdyshape%ngi))
         allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      end if
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(p_nodes_bdy(face_loc(di%pressure_mesh,1)))      
      allocate(sa_local(ele_loc(di%pressure_mesh,1)))
      allocate(sa_local_bdy(face_loc(di%pressure_mesh,1)))
         
      call zero(di%inverse_cv_sa_pressure_mesh)

      vele_loop: do vele = 1, di%number_vele

         ! The node indices of the pressure field
         p_nodes => ele_nodes(di%pressure_mesh, vele)

         ! obtain the transformed determinant*weight and normals
         if (di%cached_face_value%cached_detwei_normal) then

            detwei => di%cached_face_value%detwei(:,vele)

            normal => di%cached_face_value%normal(:,:,vele)

         else

            ! get the coordinate values for this element for each positions local node
            x_ele = ele_val(di%positions, vele)         

            call transform_cvsurf_to_physical(x_ele, di%x_cvshape, detwei, normal, di%cvfaces)

         end if

         ! Initialise array for the quadrature points of this 
         ! element for whether it has already been visited
         notvisited = .true.

         sa_local = 0.0

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

                        sa_local(iloc) = sa_local(iloc) + &
                                         detwei(ggi)

                        sa_local(oloc) = sa_local(oloc) + &
                                         detwei(ggi)

                     end if check_visited

                  end do quadrature_loop

               end if is_neigh

            end do face_loop

         end do nodal_loop_i

         call addto(di%inverse_cv_sa_pressure_mesh, p_nodes, sa_local)

      end do vele_loop
      
      sele_loop: do sele = 1, di%number_sele

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

         sa_local_bdy = 0.0

         bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

            bc_face_loop: do face = 1, di%cvfaces%sfaces

               bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                  bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                     ggi = (face-1)*di%cvfaces%shape%ngi + gi

                     sa_local_bdy(iloc) = sa_local_bdy(iloc) + &
                                          detwei_bdy(ggi)                    

                  end do bc_quad_loop

               end if bc_neigh_if

            end do bc_face_loop

         end do bc_iloc_loop

         call addto(di%inverse_cv_sa_pressure_mesh, p_nodes_bdy, sa_local_bdy)

      end do sele_loop
      
      call invert(di%inverse_cv_sa_pressure_mesh)
            
      ! deallocate local variables as required
      deallocate(x_ele)
      if (di%cached_face_value%cached_detwei_normal) then
         nullify(detwei)
         nullify(normal)      
      else
         deallocate(detwei)
         deallocate(normal)
      end if
      deallocate(notvisited)
      
      if (di%cached_face_value%cached_detwei_normal) then
         nullify(detwei_bdy)
         nullify(normal_bdy)      
      else
         deallocate(detwei_bdy)
         deallocate(normal_bdy)
      end if
      deallocate(x_ele_bdy)
      deallocate(p_nodes_bdy)

      deallocate(sa_local)
      deallocate(sa_local_bdy)
      
      ewrite_minmax(di%inverse_cv_sa_pressure_mesh)
      
      ewrite(1,*) 'Finished Calculate the inverse of the surface area of CV on the pressure mesh'
      
   end subroutine darcy_impes_calculate_inverse_cv_sa

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

   subroutine darcy_impes_calculate_relperm_fields(di)
      
      !!< Calculate the latest relperm fields from the latest saturations.
      !!< This is assumes that the relperm and saturation fields are 
      !!< on the same mesh, which is all the input permits
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer                     :: p, node
      real                        :: relperm_node_val
      real, dimension(:), pointer :: sat_node_val_all_phases

      ewrite(1,*) 'Calculate RelativePermeability field of each phase'
      
      allocate(sat_node_val_all_phases(di%number_phase))
      
      ! Calc relperm phase values for each node
      node_loop: do node = 1, di%number_pmesh_node
         
         ! Find all the phase saturation node values
         do p = 1, di%number_phase
            
            sat_node_val_all_phases(p) = node_val(di%saturation(p)%ptr, node)
            
         end do
         
         ! Calc and set relperm node value for each phase
         do p = 1, di%number_phase
            
            call darcy_impes_calculate_relperm_value(relperm_node_val, &
                                                     sat_node_val_all_phases, &
                                                     p, &
                                                     di%relperm_corr_options%type, &
                                                     di%relperm_corr_options%exponents, &
                                                     di%relperm_corr_options%residual_saturations, &
                                                     di%relperm_corr_options%cutoff_saturations, &
                                                     di%relperm_corr_options%scaling_coefficients)
            
            call set(di%relative_permeability(p)%ptr, &
                     node, &
                     relperm_node_val)
               
         end do
         
      end do node_loop

      deallocate(sat_node_val_all_phases)
      
      ewrite(1,*) 'Finished Calculate RelativePermeability field of each phase'
            
   end subroutine darcy_impes_calculate_relperm_fields

! ----------------------------------------------------------------------------

   subroutine darcy_impes_calculate_relperm_value(relperm_val, &
                                                  sat_val_all_phases, &
                                                  p, &
                                                  relperm_corr_type, &
                                                  relperm_corr_exponents, &
                                                  relperm_corr_residual_sats, &
                                                  relperm_corr_cutoff_sats, &
                                                  relperm_corr_scaling_coefficients)
      
      !!< Calculate the latest relperm value for phase p for the 
      !!< given saturation values of all phases using the given options.
      
      real,                  intent(out) :: relperm_val
      real,    dimension(:), intent(in)  :: sat_val_all_phases
      integer,               intent(in)  :: p
      integer,               intent(in)  :: relperm_corr_type
      real,    dimension(:), intent(in)  :: relperm_corr_exponents
      real,    dimension(:), intent(in)  :: relperm_corr_residual_sats
      real,    dimension(:), intent(in)  :: relperm_corr_cutoff_sats
      real,    dimension(:), intent(in)  :: relperm_corr_scaling_coefficients
    
      ! local variables
      real :: sat_minus_res_sat
      real :: sat_effective
         
      select case (relperm_corr_type)
           
      case (RELPERM_CORRELATION_POWER)

         relperm_val = (sat_val_all_phases(p) - relperm_corr_residual_sats(p)) ** relperm_corr_exponents(p)

      case (RELPERM_CORRELATION_COREY2PHASE)
         
         sat_minus_res_sat = sat_val_all_phases(1) - relperm_corr_residual_sats(1)
         
         if (p == 1) then

            relperm_val = sat_minus_res_sat ** 4

         else if (p == 2) then

            relperm_val = (1.0 - sat_minus_res_sat ** 2) * (1.0 - sat_minus_res_sat) ** 2

         else 

            ! This has already been option checked so should not happen
            FLAbort('Trying to use Corey2Phase relative permeabiltiy correlation for simulation with more than 2 phases')

         end if

      case (RELPERM_CORRELATION_COREY2PHASEOPPOSITE)

         sat_minus_res_sat = sat_val_all_phases(2) - relperm_corr_residual_sats(2)

         if (p == 1) then

            relperm_val = (1.0 - sat_minus_res_sat ** 2) * (1.0 - sat_minus_res_sat) ** 2

         else if (p == 2) then

            relperm_val = sat_minus_res_sat ** 4

         else 

            ! This has already been option checked so should not happen
            FLAbort('Trying to use Corey2PhaseOpposite relative permeabiltiy correlation for simulation with more than 2 phases')

         end if

      case (RELPERM_CORRELATION_MINERAL)        
         
         relperm_val = ((max(relperm_corr_cutoff_sats(p), sat_val_all_phases(p)) - relperm_corr_residual_sats(p)) **&
                      &relperm_corr_exponents(p)) /  ((1.0 - relperm_corr_residual_sats(p)) ** relperm_corr_exponents(p))

      case (RELPERM_CORRELATION_VANGENUCHTEN)     

         sat_effective = (sat_val_all_phases(2) - relperm_corr_residual_sats(2) ) / ( 1.0 - relperm_corr_residual_sats(1) - &
                         &relperm_corr_residual_sats(2))
         
         if (p == 1) then

            relperm_val = ( 1.0 - sat_effective )**(1.0/3.0) * (1.0 - sat_effective ** ( 1.0/ relperm_corr_exponents(p)) ) ** (2* relperm_corr_exponents(p))

         else if (p == 2) then

            relperm_val =  sat_effective ** (1.0/2.0) * ( 1.0 - ( 1.0 - sat_effective ** ( 1.0 / relperm_corr_exponents(p))) ** relperm_corr_exponents(p))**2.0

         else 

            ! This has already been option checked so should not happen
            FLAbort('Trying to use VanGenuchten relative permeabiltiy correlation for simulation with more than 2 phases')

         end if         

      case (RELPERM_CORRELATION_JACKSON2PHASE)
         
         sat_effective = (sat_val_all_phases(2) - relperm_corr_residual_sats(2)) / &
                         (1.0 - relperm_corr_residual_sats(1) - relperm_corr_residual_sats(2))
	 if ( sat_effective >= 1.0 ) then
	    sat_effective = 1.0
	 end if 
         
	 if ( sat_effective <= 0.0 ) then
	    sat_effective = 0.0
	 end if 

         if (p == 1) then

            relperm_val = relperm_corr_scaling_coefficients(1) * (1.0 - sat_effective) ** relperm_corr_exponents(1)

         else if (p == 2) then

            relperm_val = relperm_corr_scaling_coefficients(2) * sat_effective ** relperm_corr_exponents(2)

         else 

            ! This has already been option checked so should not happen
            FLAbort('Trying to use Jackson2Phase relative permeabiltiy correlation for simulation with more than 2 phases')

         end if

      case (RELPERM_CORRELATION_JACKSON2PHASEOPPOSITE)

         sat_effective = (sat_val_all_phases(1) - relperm_corr_residual_sats(1)) / &
                         (1.0 - relperm_corr_residual_sats(1) - relperm_corr_residual_sats(2))
         
         if (p == 1) then

            relperm_val = relperm_corr_scaling_coefficients(1) * sat_effective ** relperm_corr_exponents(1)

         else if (p == 2) then

            relperm_val = relperm_corr_scaling_coefficients(2) * (1.0 - sat_effective) ** relperm_corr_exponents(2)

         else 

            ! This has already been option checked so should not happen
            FLAbort('Trying to use Jackson2PhaseOpposite relative permeabiltiy correlation for simulation with more than 2 phases')

         end if
     
     end select
          
   end subroutine darcy_impes_calculate_relperm_value

! ----------------------------------------------------------------------------
    
   subroutine darcy_impes_calculate_cflnumber_field_based_dt(di)
      
      !!< Calculate the cfl number field based adaptive time step
      
      type(darcy_impes_type), intent(inout) :: di

     !! NEW NoT TESTED FROM HERE
      type(logical) :: diff_cfl


      real :: ds, dpds, diff, deltax
      type(scalar_field), pointer::sat, phi, perm, vis
     
      type(vector_field), pointer:: X

      real, dimension(:), allocatable :: detwei
      real, dimension(:), pointer :: phi_ele, perm_ele, sat_ele, vis_ele
      integer, dimension(:), pointer :: sat_ele_nodes
      integer :: i , ele
      real, dimension(:, :, :), allocatable :: invJ
      real :: sat_node
      type(state_type) :: state

     !! TILL HERE
  
      
      ! local variables
      integer :: p
      real, dimension(di%number_phase) :: phase_dt, dtdiff
      
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


     ! NEW NOT TESTED FROM HERE 
!     diff_cfl = .False.        	
!     if (diff_cfl) then 

!      do p = 1, di%number_phase
     
!        sat=>extract_scalar_field(state(p), "Saturation")
!        phi=>extract_scalar_field(state(p), "Porosity")
!        perm=>extract_scalar_field(state(p), "Permeability")
!        vis=>extract_scalar_field(state(p), "Viscosity")

!        X=>extract_vector_field(state(1), "Coordinate")
!        allocate (invj(mesh_dim(x), mesh_dim(x), ele_ngi(x, 1)))
!        allocate (detwei(ele_ngi(X,ele)))


!        ds = 1.e-10
!        dpds=0.
    
!        diff= 0.

!        do ele = 1, ele_count(sat)
        
!          sat_ele_nodes => ele_nodes(sat, ele)
 
!          phi_ele = ele_val(phi, ele)
!          perm_ele = ele_val(perm, ele)
!          vis_ele = ele_val(vis, ele)

 !         call compute_inverse_jacobian(ele_val(X,ele), ele_shape(X,ele), &
  !                                     detwei=detwei, invJ=invJ)


 
!          deltax = minval(invJ) 

!          do i  = 1, size(sat_ele_nodes)
 
!             sat_node = node_val(sat, sat_ele_nodes(i) )
!             dpds =  abs( ( (sat_node+ds)**(-0.5) - (sat_node-ds)**(-0.5) ) /(2.*ds) )

 !            diff = max( diff, ( 1./phi_ele(ele)) * (perm_ele(ele) / vis_ele(ele)) * dpds )

!          end do 
         
!          dtdiff(p) = deltax**2 / (2*diff)

 !       end do

!       deallocate (invj)
!       deallocate (detwei)

!     end do
     

      
!      di%dt = min(minval(dtdiff), minval(phase_dt))

!     end if 

! UNTIL HERE

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
            
      if (di%subcy_opt_sat%have .and. di%subcy_opt_sat%consistent) then
      
         allocate(di%cached_face_value%relperm(di%subcy_opt_sat%number, di%p_cvshape%ngi, di%number_vele, di%number_phase))
         allocate(di%cached_face_value%relperm_bdy(di%subcy_opt_sat%number, di%p_cvbdyshape%ngi, di%number_sele, di%number_phase))
         
      else
         
         allocate(di%cached_face_value%relperm(1, di%p_cvshape%ngi, di%number_vele, di%number_phase))
         allocate(di%cached_face_value%relperm_bdy(1, di%p_cvbdyshape%ngi, di%number_sele, di%number_phase))
      
      end if
      
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
      integer :: r_upwind_pos, d_upwind_pos, s_upwind_pos, dim, sat_p
      real    :: income, v_over_relperm_dot_n, old_relperm_face_value, den_face_value, old_sat_face_value
      real,    dimension(1)              :: absperm_ele, visc_ele
      real,    dimension(:),     pointer :: old_relperm_ele
      real,    dimension(:),     pointer :: den_ele
      real,    dimension(:,:),   pointer :: old_sat_ele
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
      integer, dimension(:),     pointer :: s_upwind_nodes
      real,    dimension(:),     pointer :: den_ele_bdy
      real,    dimension(:),     pointer :: old_relperm_ele_bdy
            
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
      allocate(old_sat_ele(ele_loc(di%pressure_mesh,1), di%number_phase))
      allocate(grav_ele(di%ndim,1))
      allocate(grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(v_over_relperm_face_quad(di%ndim,di%p_cvshape%ngi))

      allocate(den_ele_bdy(face_loc(di%pressure_mesh,1)))
      allocate(old_relperm_ele_bdy(face_loc(di%pressure_mesh,1)))      
      
      ewrite(1,*) 'Calculate relperm and density and first face values'      

      ! Initialise the cached face values
      di%cached_face_value%relperm(1,:,:,:)     = 0.0
      di%cached_face_value%relperm_bdy(1,:,:,:) = 0.0
      di%cached_face_value%den                  = 0.0
      di%cached_face_value%den_bdy              = 0.0

      phase_loop: do p = 1, di%number_phase
          
         ! Determine the upwind relperm values of all node pairs if required for higher order CV face value.
         ! If the relperm face value is a RELPERMOVERSAT... scheme then the upwind values 
         ! of modrelperm are required (= relperm / sat)
         if(di%relperm_cv_options%limit_facevalue) then
           
            if (di%determine_saturation_face_values) then
               
               call darcy_impes_calculate_old_modified_relative_permeability(di, p)
               
               ! NOTE only old values are required but this procedure 
               !      expects latest and old. So pass in old twice - not optimal
               call find_upwind_values(di%state, &
                                       di%positions_pressure_mesh, &
                                       di%modified_relative_permeability, &
                                       di%relperm_upwind, &
                                       di%modified_relative_permeability, &
                                       di%relperm_upwind, &
                                       option_path = trim(di%relative_permeability(p)%ptr%option_path))
           
            else 
           
               ! NOTE only old values are required but this procedure 
               !      expects latest and old. So pass in old twice - not optimal
               call find_upwind_values(di%state, &
                                       di%positions_pressure_mesh, &
                                       di%old_relative_permeability(p)%ptr, &
                                       di%relperm_upwind, &
                                       di%old_relative_permeability(p)%ptr, &
                                       di%relperm_upwind, &
                                       option_path = trim(di%relative_permeability(p)%ptr%option_path))

            end if
            
         end if

         ! If the relperm face value is a RELPERMOVERSAT... scheme then the saturation 
         ! face values need determining which may require the upwind values. 
         if (di%determine_saturation_face_values .and. di%saturation_cv_options%limit_facevalue) then

            ! NOTE only old values are required but this procedure       
            !      expects latest and old. So pass in old twice - not optimal
            call find_upwind_values(di%state, &
                                    di%positions_pressure_mesh, &
                                    di%old_saturation(p)%ptr, &
                                    di%sfield_upwind, &
                                    di%old_saturation(p)%ptr, &
                                    di%sfield_upwind, &
                                    option_path = trim(di%saturation(p)%ptr%option_path))            

         end if 
                          
         ! Initialise optimisation flag used in finding upwind value in high resolution schemes
         r_upwind_pos = 0
         s_upwind_pos = 0
            
         vol_element_loop_relperm: do vele = 1, di%number_vele
           
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

            if (di%determine_saturation_face_values) then
            
               ! Determine the node numbers to use to determine the saturation upwind values
               if((di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
                  (di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

                  s_upwind_nodes => x_pmesh_nodes

               else

                  s_upwind_nodes => p_nodes

               end if
               
               do sat_p = 1, di%number_phase
               
                  ! get the old saturation ele values for this phase
                  old_sat_ele(:,sat_p) = ele_val(di%old_saturation(sat_p)%ptr, vele)
               
               end do
               
            end if

            ! get the old_relperm ele values for this phase
            old_relperm_ele = ele_val(di%old_relative_permeability(p)%ptr, vele)
            
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
            nodal_loop_i_relperm: do iloc = 1, di%pressure_mesh%shape%loc

              ! loop over CV faces internal to this element
              face_loop_relperm: do face = 1, di%cvfaces%faces

                ! is this a face neighbouring iloc?
                is_neigh_relperm: if(di%cvfaces%neiloc(iloc, face) /= 0) then

                  ! find the opposing local node across the CV face
                  oloc = di%cvfaces%neiloc(iloc, face)

                  ! loop over gauss points on face
                  quadrature_loop_relperm: do gi = 1, di%cvfaces%shape%ngi

                    ! global gauss pt index for this element
                    ggi = (face-1)*di%cvfaces%shape%ngi + gi

                    ! check if this quadrature point has already been visited
                    check_visited_relperm: if(notvisited(ggi)) then

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

                       ! Evaluate the face value for relperm 
                       call darcy_impes_evaluate_relperm_face_val(old_relperm_face_value, &
                                                                  iloc, &
                                                                  oloc, &
                                                                  ggi, &
                                                                  r_upwind_nodes, &
                                                                  di%p_cvshape, &
                                                                  old_relperm_ele, &
                                                                  di%relperm_upwind, &
                                                                  inflow, &
                                                                  income, &
                                                                  old_sat_ele, &
                                                                  di%minimum_denominator_saturation_value, &
                                                                  p, &
                                                                  di%number_phase, &
                                                                  di%relperm_corr_options, &
                                                                  di%relperm_cv_options, &
                                                                  save_pos = r_upwind_pos)

                       ! Evaluate the face value for saturation if required, else set to 1.0
                       if (di%determine_saturation_face_values) then
                       
                          call darcy_impes_evaluate_sfield_face_val(old_sat_face_value, &
                                                                    iloc, &
                                                                    oloc, &
                                                                    ggi, &
                                                                    s_upwind_nodes, &
                                                                    di%p_cvshape, &
                                                                    old_sat_ele(:,p), &
                                                                    di%sfield_upwind, &
                                                                    inflow, &
                                                                    income, &
                                                                    di%saturation_cv_options, &
                                                                    save_pos = s_upwind_pos)                          
                       
                       else
                       
                          old_sat_face_value = 1.0
                                              
                       end if
                                              
                       ! Cache the old relperm (or modrelperm*sat) face value
                       di%cached_face_value%relperm(1,ggi,vele,p) = old_relperm_face_value * old_sat_face_value
                                           
                    end if check_visited_relperm

                  end do quadrature_loop_relperm

                end if is_neigh_relperm

              end do face_loop_relperm

            end do nodal_loop_i_relperm

         end do vol_element_loop_relperm
 
         ! Determine the upwind density values of all node pairs if required for higher order CV face value
         if(di%density_cv_options%limit_facevalue) then

           ! NOTE only new upwind values are required but this procedure 
           !      expects latest and old. So pass in new twice - not optimal
           call find_upwind_values(di%state, &
                                   di%positions_pressure_mesh, &
                                   di%density(p)%ptr, &
                                   di%sfield_upwind, &
                                   di%density(p)%ptr, &
                                   di%sfield_upwind, &
                                   option_path = trim(di%density(p)%ptr%option_path))

         end if
                          
         ! Initialise optimisation flag used in finding upwind value in high resolution schemes
         d_upwind_pos = 0
            
         vol_element_loop_den: do vele = 1, di%number_vele
           
            ! The node indices of the pressure mesh
            p_nodes => ele_nodes(di%pressure_mesh, vele)

            ! The node indices of the positions projected to the pressure mesh
            x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)
            
            ! Determine the node numbers to use to determine the density upwind values
            if((di%density_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
               (di%density_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

               d_upwind_nodes => x_pmesh_nodes

            else

               d_upwind_nodes => p_nodes

            end if

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
            nodal_loop_i_den: do iloc = 1, di%pressure_mesh%shape%loc

              ! loop over CV faces internal to this element
              face_loop_den: do face = 1, di%cvfaces%faces

                ! is this a face neighbouring iloc?
                is_neigh_den: if(di%cvfaces%neiloc(iloc, face) /= 0) then

                  ! find the opposing local node across the CV face
                  oloc = di%cvfaces%neiloc(iloc, face)

                  ! loop over gauss points on face
                  quadrature_loop_den: do gi = 1, di%cvfaces%shape%ngi

                    ! global gauss pt index for this element
                    ggi = (face-1)*di%cvfaces%shape%ngi + gi

                    ! check if this quadrature point has already been visited
                    check_visited_den: if(notvisited(ggi)) then

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
                                              
                       ! Evaluate the density face value 
                       call darcy_impes_evaluate_density_face_val(den_face_value, &
                                                                  iloc, &
                                                                  oloc, &
                                                                  ggi, &
                                                                  d_upwind_nodes, &
                                                                  di%p_cvshape, &
                                                                  den_ele, &
                                                                  di%sfield_upwind, &
                                                                  inflow, &
                                                                  income, &
                                                                  di%density_cv_options, &
                                                                  save_pos = d_upwind_pos)
                                             
                       ! Cache the density face value
                       di%cached_face_value%den(ggi,vele,p) = den_face_value
                                           
                    end if check_visited_den

                  end do quadrature_loop_den

                end if is_neigh_den

              end do face_loop_den

            end do nodal_loop_i_den

         end do vol_element_loop_den
                  
         sele_loop: do sele = 1, di%number_sele
            
            ! get the density sele values for this phase
            den_ele_bdy = face_val(di%density(p)%ptr, sele)
            
            ! get the old relperm domain value
            old_relperm_ele_bdy = face_val(di%old_relative_permeability(p)%ptr, sele)
                        
            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi
                                                       
                        ! Evaluate the relperm and density face value via taking the CV value
                        ! - no other choice currently ...

                        di%cached_face_value%relperm_bdy(1,ggi,sele,p) = old_relperm_ele_bdy(iloc)

                        di%cached_face_value%den_bdy(ggi,sele,p) = den_ele_bdy(iloc)
                                                   
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
      deallocate(old_sat_ele)
      deallocate(grav_ele)
      deallocate(grad_pressure_face_quad)
      deallocate(v_over_relperm_face_quad)

      deallocate(den_ele_bdy)
      deallocate(old_relperm_ele_bdy)
      
      ewrite(1,*) 'Finished calculate relperm and density first face values'      
      
   end subroutine darcy_impes_calculate_relperm_den_first_face_values

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_relperm_isub_face_values(di, isub_index)
   
      !!< Find the relperm  CV face values for all phase and cache for the 
      !!< start of a subcycle. isub_index flags where to cache the value.
            
      type(darcy_impes_type), intent(inout) :: di
      integer,                intent(in)    :: isub_index
      
      ! local variables
      logical :: inflow
      integer :: p, sat_p, vele, sele, iloc, oloc, face, gi, ggi, r_upwind_pos, s_upwind_pos, dim
      real    :: income, v_over_relperm_dot_n, relperm_face_value, sat_face_value
      real,    dimension(1)              :: absperm_ele, visc_ele
      real,    dimension(:),     pointer :: relperm_ele
      real,    dimension(:,:),   pointer :: sat_ele
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
      integer, dimension(:),     pointer :: s_upwind_nodes
      real,    dimension(:),     pointer :: relperm_ele_bdy
            
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
      allocate(relperm_ele(ele_loc(di%pressure_mesh,1)))
      allocate(sat_ele(ele_loc(di%pressure_mesh,1), di%number_phase))
      allocate(grav_ele(di%ndim,1))
      allocate(grad_pressure_face_quad(di%ndim, di%p_cvshape%ngi))
      allocate(v_over_relperm_face_quad(di%ndim,di%p_cvshape%ngi))

      allocate(relperm_ele_bdy(face_loc(di%pressure_mesh,1)))      
      
      ewrite(1,*) 'Calculate relperm and first face values for subcycle'      

      ! Initialise the cached face values
      di%cached_face_value%relperm(isub_index,:,:,:)     = 0.0
      di%cached_face_value%relperm_bdy(isub_index,:,:,:) = 0.0
            
      phase_loop: do p = 1, di%number_phase
          
         ! Determine the upwind relperm values of all node pairs if required for higher order CV face value.
         ! If the relperm face value is a RELPERMOVERSAT... scheme then the upwind values 
         ! of modrelperm are required (= relperm / sat)
         if(di%relperm_cv_options%limit_facevalue) then
           
            if (di%determine_saturation_face_values) then
               
               call darcy_impes_calculate_modified_relative_permeability(di, p)
               
               ! NOTE only new values are required but this procedure 
               !      expects latest and old. So pass in new twice - not optimal
               call find_upwind_values(di%state, &
                                       di%positions_pressure_mesh, &
                                       di%modified_relative_permeability, &
                                       di%relperm_upwind, &
                                       di%modified_relative_permeability, &
                                       di%relperm_upwind, &
                                       option_path = trim(di%relative_permeability(p)%ptr%option_path))
           
            else 
           
               ! NOTE only new values are required but this procedure 
               !      expects latest and old. So pass in new twice - not optimal
               call find_upwind_values(di%state, &
                                       di%positions_pressure_mesh, &
                                       di%relative_permeability(p)%ptr, &
                                       di%relperm_upwind, &
                                       di%relative_permeability(p)%ptr, &
                                       di%relperm_upwind, &
                                       option_path = trim(di%relative_permeability(p)%ptr%option_path))

            end if
            
         end if

         ! If the relperm face value is a RELPERMOVERSAT... scheme then the saturation 
         ! face values need determining which may require the upwind values. 
         if (di%determine_saturation_face_values .and. di%saturation_cv_options%limit_facevalue) then

            ! NOTE only new values are required but this procedure       
            !      expects latest and old. So pass in new twice - not optimal
            call find_upwind_values(di%state, &
                                    di%positions_pressure_mesh, &
                                    di%saturation(p)%ptr, &
                                    di%sfield_upwind, &
                                    di%saturation(p)%ptr, &
                                    di%sfield_upwind, &
                                    option_path = trim(di%saturation(p)%ptr%option_path))            

         end if 
                           
         ! Initialise optimisation flags used in finding upwind value in high resolution schemes
         r_upwind_pos = 0
         s_upwind_pos = 0
            
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
            
            if (di%determine_saturation_face_values) then
            
               ! Determine the node numbers to use to determine the saturation upwind values
               if((di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
                  (di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

                  s_upwind_nodes => x_pmesh_nodes

               else

                  s_upwind_nodes => p_nodes

               end if
               
               do sat_p = 1, di%number_phase
               
                  ! get the saturation ele values for this phase
                  sat_ele(:,sat_p) = ele_val(di%saturation(sat_p)%ptr, vele)
               
               end do
               
            end if

            ! get the relperm ele values for this phase
            relperm_ele = ele_val(di%relative_permeability(p)%ptr, vele)
            
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

                           ! Evaluate the face value for relperm 
                           call darcy_impes_evaluate_relperm_face_val(relperm_face_value, &
                                                                      iloc, &
                                                                      oloc, &
                                                                      ggi, &
                                                                      r_upwind_nodes, &
                                                                      di%p_cvshape, &
                                                                      relperm_ele, &
                                                                      di%relperm_upwind, &
                                                                      inflow, &
                                                                      income, &
                                                                      sat_ele, &
                                                                      di%minimum_denominator_saturation_value, &
                                                                      p, &
                                                                      di%number_phase, &
                                                                      di%relperm_corr_options, &
                                                                      di%relperm_cv_options, &
                                                                      save_pos = r_upwind_pos)

                           ! Evaluate the face value for saturation if required, else set to 1.0
                           if (di%determine_saturation_face_values) then

                              call darcy_impes_evaluate_sfield_face_val(sat_face_value, &
                                                                        iloc, &
                                                                        oloc, &
                                                                        ggi, &
                                                                        s_upwind_nodes, &
                                                                        di%p_cvshape, &
                                                                        sat_ele(:,p), &
                                                                        di%sfield_upwind, &
                                                                        inflow, &
                                                                        income, &
                                                                        di%saturation_cv_options, &
                                                                        save_pos = s_upwind_pos)                          

                           else

                              sat_face_value = 1.0

                           end if

                           ! Cache the relperm (or modrelperm*sat) face value
                           di%cached_face_value%relperm(isub_index,ggi,vele,p) = relperm_face_value * sat_face_value

                        end if check_visited

                     end do quadrature_loop

                  end if is_neigh

               end do face_loop

            end do nodal_loop_i

         end do vol_element_loop
         
         sele_loop: do sele = 1, di%number_sele
                        
            ! get the relperm domain value
            relperm_ele_bdy = face_val(di%relative_permeability(p)%ptr, sele)
                        
            bc_iloc_loop: do iloc = 1, di%pressure_mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi
                                                       
                        ! Evaluate the relperm face value via taking the CV value
                        ! - no other choice currently ...

                        di%cached_face_value%relperm_bdy(isub_index,ggi,sele,p) = relperm_ele_bdy(iloc)
                                                   
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
      deallocate(relperm_ele)
      deallocate(sat_ele)
      deallocate(grav_ele)
      deallocate(grad_pressure_face_quad)
      deallocate(v_over_relperm_face_quad)

      deallocate(relperm_ele_bdy)
      
      ewrite(1,*) 'Finished calculate relperm face values for subcycle'      
      
   end subroutine darcy_impes_calculate_relperm_isub_face_values

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
   
   subroutine darcy_impes_evaluate_sfield_face_val(sfield_face_val, &
                                                   iloc, &
                                                   oloc, &
                                                   ggi, &
                                                   upwind_nodes, &
                                                   cvshape, &
                                                   sfield_ele, &
                                                   sfield_upwind, &
                                                   inflow, &
                                                   income, &
                                                   sfield_cv_options, &
                                                   save_pos)
   
      !!< Evaluate the high resolution sfield CV face value deduced from the 
      !!< the cv_options and solution field
            
      real,                              intent(out) :: sfield_face_val
      integer,                           intent(in)  :: iloc
      integer,                           intent(in)  :: oloc
      integer,                           intent(in)  :: ggi
      integer,            dimension(:),  intent(in)  :: upwind_nodes
      type(element_type),                intent(in)  :: cvshape
      real,               dimension(:),  intent(in)  :: sfield_ele
      type(csr_matrix),                  intent(in)  :: sfield_upwind
      logical,                           intent(in)  :: inflow
      real,                              intent(in)  :: income
      type(darcy_impes_cv_options_type), intent(in)  :: sfield_cv_options 
      integer,                           intent(inout), optional :: save_pos
      
      ! local variables
      real    :: upwind_val, donor_val, downwind_val, cfl_donor
      integer :: l_save_pos
      
      if(present(save_pos)) then
        ! an attempt at optimising the val calls by saving the matrix position
        l_save_pos=save_pos 
      else
        l_save_pos = 0
      end if

      select case(sfield_cv_options%facevalue)

      case (DARCY_IMPES_CV_FACEVALUE_FIRSTORDERUPWIND)

         sfield_face_val = income*sfield_ele(oloc) + (1.-income)*sfield_ele(iloc)

      case (DARCY_IMPES_CV_FACEVALUE_FINITEELEMENT)

         sfield_face_val = dot_product(cvshape%n(:,ggi), sfield_ele)

         if(sfield_cv_options%limit_facevalue) then

           downwind_val = income*sfield_ele(iloc) + (1.-income)*sfield_ele(oloc)

           donor_val = income*sfield_ele(oloc) + (1.-income)*sfield_ele(iloc)
                      
           if(inflow) then
             upwind_val = val(sfield_upwind, upwind_nodes(oloc), upwind_nodes(iloc), save_pos=l_save_pos)
           else
             upwind_val = val(sfield_upwind, upwind_nodes(iloc), upwind_nodes(oloc), save_pos=l_save_pos)
           end if

           sfield_face_val = limit_val(upwind_val, &
                                       donor_val, &
                                       downwind_val, &
                                       sfield_face_val, &
                                       sfield_cv_options%limiter, &
                                       cfl_donor, &
                                       sfield_cv_options%limiter_slopes)

         end if

         if(present(save_pos)) then
           save_pos = l_save_pos
         end if

      case default

         FLAbort('Unknown CV facevalue scheme for sfield')
            
      end select
                  
   end subroutine darcy_impes_evaluate_sfield_face_val

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_evaluate_density_face_val(density_face_val, &
                                                    iloc, &
                                                    oloc, &
                                                    ggi, &
                                                    upwind_nodes, &
                                                    cvshape, &
                                                    density_ele, &
                                                    density_upwind, &
                                                    inflow, &
                                                    income, &
                                                    density_cv_options, &
                                                    save_pos)
   
      !!< Evaluate the high resolution density CV face value deduced from the 
      !!< the cv_options and solution field
            
      real,                              intent(out) :: density_face_val
      integer,                           intent(in)  :: iloc
      integer,                           intent(in)  :: oloc
      integer,                           intent(in)  :: ggi
      integer,            dimension(:),  intent(in)  :: upwind_nodes
      type(element_type),                intent(in)  :: cvshape
      real,               dimension(:),  intent(in)  :: density_ele
      type(csr_matrix),                  intent(in)  :: density_upwind
      logical,                           intent(in)  :: inflow
      real,                              intent(in)  :: income
      type(darcy_impes_cv_options_type), intent(in)  :: density_cv_options 
      integer,                           intent(inout), optional :: save_pos
      
      ! local variables
      real    :: upwind_val, donor_val, downwind_val, cfl_donor
      integer :: l_save_pos
      
      if(present(save_pos)) then
        ! an attempt at optimising the val calls by saving the matrix position
        l_save_pos=save_pos 
      else
        l_save_pos = 0
      end if

      select case(density_cv_options%facevalue)

      case (DARCY_IMPES_CV_FACEVALUE_FIRSTORDERUPWIND)

         density_face_val = income*density_ele(oloc) + (1.-income)*density_ele(iloc)

      case (DARCY_IMPES_CV_FACEVALUE_FINITEELEMENT)

         density_face_val = dot_product(cvshape%n(:,ggi), density_ele)

         if(density_cv_options%limit_facevalue) then

           downwind_val = income*density_ele(iloc) + (1.-income)*density_ele(oloc)

           donor_val = income*density_ele(oloc) + (1.-income)*density_ele(iloc)
                      
           if(inflow) then
             upwind_val = val(density_upwind, upwind_nodes(oloc), upwind_nodes(iloc), save_pos=l_save_pos)
           else
             upwind_val = val(density_upwind, upwind_nodes(iloc), upwind_nodes(oloc), save_pos=l_save_pos)
           end if

           density_face_val = limit_val(upwind_val, &
                                        donor_val, &
                                        downwind_val, &
                                        density_face_val, &
                                        density_cv_options%limiter, &
                                        cfl_donor, &
                                        density_cv_options%limiter_slopes)

         end if

         if(present(save_pos)) then
           save_pos = l_save_pos
         end if

      case default

         FLAbort('Unknown CV facevalue scheme for density')

      end select
                  
   end subroutine darcy_impes_evaluate_density_face_val

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_evaluate_relperm_face_val(relperm_face_val, &
                                                    iloc, &
                                                    oloc, &
                                                    ggi, &
                                                    upwind_nodes, &
                                                    cvshape, &
                                                    relperm_ele, &
                                                    relperm_upwind, &
                                                    inflow, &
                                                    income, &
                                                    sat_ele_all_phases, &
                                                    min_denom_sat, &
                                                    p, &
                                                    number_phase, &
                                                    relperm_corr_options, &
                                                    relperm_cv_options, &
                                                    save_pos)
   
      !!< Evaluate the high resolution relperm CV face value deduced from the 
      !!< the cv_options, solution field and perhaps the given saturation 
      !!< face and element values.
            
      real,                                                        intent(out) :: relperm_face_val
      integer,                                                     intent(in)  :: iloc
      integer,                                                     intent(in)  :: oloc
      integer,                                                     intent(in)  :: ggi
      integer,                                     dimension(:),   intent(in)  :: upwind_nodes
      type(element_type),                                          intent(in)  :: cvshape
      real,                                        dimension(:),   intent(in)  :: relperm_ele
      type(csr_matrix),                                            intent(in)  :: relperm_upwind
      logical,                                                     intent(in)  :: inflow
      real,                                                        intent(in)  :: income
      real,                                        dimension(:,:), intent(in)  :: sat_ele_all_phases
      real,                                                        intent(in)  :: min_denom_sat
      integer,                                                     intent(in)  :: p
      integer,                                                     intent(in)  :: number_phase
      type(darcy_impes_relperm_corr_options_type),                 intent(in)  :: relperm_corr_options
      type(darcy_impes_cv_options_type),                           intent(in)  :: relperm_cv_options 
      integer,                                                     intent(inout), optional :: save_pos
      
      ! local variables
      integer :: l_save_pos, sat_p
      real    :: upwind_val, donor_val, downwind_val, cfl_donor
      real, dimension(number_phase) :: sat_face_val_all_phases
      
      if(present(save_pos)) then
        ! an attempt at optimising the val calls by saving the matrix position
        l_save_pos=save_pos 
      else
        l_save_pos = 0
      end if

      select case(relperm_cv_options%facevalue)

      case (DARCY_IMPES_CV_FACEVALUE_FIRSTORDERUPWIND)

         relperm_face_val = income*relperm_ele(oloc) + (1.-income)*relperm_ele(iloc)
      
      case (DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATUPWIND)

         relperm_face_val = income*(relperm_ele(oloc) / max(sat_ele_all_phases(oloc,p),min_denom_sat)) + &
                           (1.-income)*(relperm_ele(iloc) / max(sat_ele_all_phases(iloc,p),min_denom_sat))

      case (DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATFINITEELEMENT)

         relperm_face_val = dot_product(cvshape%n(:,ggi), relperm_ele) / &
                            max(dot_product(cvshape%n(:,ggi), sat_ele_all_phases(:,p)),min_denom_sat)

      case (DARCY_IMPES_CV_FACEVALUE_RELPERMOVERSATCORRELATION)
         
         ! Use the FE interpolation of all phases saturations
         ! to calculate the modified relative permeability face value.
         do sat_p = 1, number_phase
         
            sat_face_val_all_phases(sat_p) = dot_product(cvshape%n(:,ggi), sat_ele_all_phases(:,sat_p))
         
         end do
         
         call darcy_impes_calculate_relperm_value(relperm_face_val, &
                                                  sat_face_val_all_phases, &
                                                  p, &
                                                  relperm_corr_options%type, &
                                                  relperm_corr_options%exponents, &
                                                  relperm_corr_options%residual_saturations, &
                                                  relperm_corr_options%cutoff_saturations, &
                                                  relperm_corr_options%scaling_coefficients)
         
         relperm_face_val = relperm_face_val / max(sat_face_val_all_phases(p),min_denom_sat)
         
      case default

         FLAbort('Unknown CV facevalue scheme for relative permeability')

      end select

      if(relperm_cv_options%limit_facevalue) then

        downwind_val = income*relperm_ele(iloc) + (1.-income)*relperm_ele(oloc)

        donor_val = income*relperm_ele(oloc) + (1.-income)*relperm_ele(iloc)

        if(inflow) then
          upwind_val = val(relperm_upwind, upwind_nodes(oloc), upwind_nodes(iloc), save_pos=l_save_pos)
        else
          upwind_val = val(relperm_upwind, upwind_nodes(iloc), upwind_nodes(oloc), save_pos=l_save_pos)
        end if

        relperm_face_val = limit_val(upwind_val, &
                                     donor_val, &
                                     downwind_val, &
                                     relperm_face_val, &
                                     relperm_cv_options%limiter, &
                                     cfl_donor, &
                                     relperm_cv_options%limiter_slopes)

      end if

      if(present(save_pos)) then
        save_pos = l_save_pos
      end if
                  
   end subroutine darcy_impes_evaluate_relperm_face_val

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_modified_relative_permeability(di, p)
   
      !!< Calculate the modified relative permeability = relperm / max( S, S_min) for phase p
            
      type(darcy_impes_type), intent(inout) :: di 
      integer,                intent(in)    :: p
      
      ! local variable
      integer :: i

      ewrite(1,*) 'Calculate modified relative permeability for phase ',p
               
      node_loop: do i = 1, node_count(di%pressure_mesh)

         di%modified_relative_permeability%val(i) = &
         di%relative_permeability(p)%ptr%val(i) / &
         max(di%saturation(p)%ptr%val(i), di%minimum_denominator_saturation_value)

      end do node_loop

      ewrite_minmax(di%modified_relative_permeability)
      
      ewrite(1,*) 'Finished Calculate modified relative permeability for phase ',p
      
   end subroutine darcy_impes_calculate_modified_relative_permeability

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_old_modified_relative_permeability(di, p)
   
      !!< Calculate the old modified relative permeability = old_relperm / max( old_S, S_min) for phase p
            
      type(darcy_impes_type), intent(inout) :: di 
      integer,                intent(in)    :: p
      
      ! local variable
      integer :: i

      ewrite(1,*) 'Calculate old modified relative permeability for phase ',p
               
      node_loop: do i = 1, node_count(di%pressure_mesh)

         di%modified_relative_permeability%val(i) = &
         di%old_relative_permeability(p)%ptr%val(i) / &
         max(di%old_saturation(p)%ptr%val(i), di%minimum_denominator_saturation_value)

      end do node_loop

      ewrite_minmax(di%modified_relative_permeability)
      
      ewrite(1,*) 'Finished Calculate old modified relative permeability for phase ',p
      
   end subroutine darcy_impes_calculate_old_modified_relative_permeability

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

! ******************26 July 2013 LCai ****************************************!
!Slove the mobile saturation if MIM exist
subroutine darcy_trans_MIM_assemble_and_solve_mobile_saturation(di)

        type(darcy_impes_type), intent(inout) :: di
        integer :: i
        type(scalar_field), pointer :: total_sat => null()  !total saturation 
        type(scalar_field), pointer :: immobile_sat  => null()  ! immobile saturation
        type(scalar_field)  :: mobile_sat   !mobile saturation

        call allocate(mobile_sat, di%pressure_mesh)
        
        do i= 2, di%number_phase

          total_sat      => di%saturation(i)%ptr
          immobile_sat   => di%MIM_options%immobile_saturation(i)%ptr

          if (di%MIM_options%have_MIM(i)) then
             ewrite(1, *) "calculate the mobile saturation of phase: ", i
             call set(mobile_sat, total_sat)
             call addto(mobile_sat, immobile_sat, scale=-1.0)
             call set(di%MIM_options%mobile_saturation(i)%ptr, mobile_sat)
          end if
          nullify(immobile_sat, total_sat)
          call zero(mobile_sat)
        end do
        
        call deallocate(mobile_sat)

end subroutine darcy_trans_MIM_assemble_and_solve_mobile_saturation


! ******** 16 July 2013 LCai *************************************************
!solve the immobile prog sfield 
subroutine darcy_trans_assemble_and_solve_immobile_prog_sfield(di, p, f, imf, temp_MIM_src)
       
       type(darcy_impes_type), intent(inout) :: di
       type(scalar_field), intent(in) :: temp_MIM_src
       integer, intent(in) :: p
       integer, intent(in) :: f
       integer, intent(in) :: imf

       ! local variable
       type(scalar_field) :: temp_rhs
       
       call allocate(temp_rhs, di%pressure_mesh)
       
       ! calculate '(old_theta_s*old_C_s)/(theta_s+alpha*dt)'
       call set(temp_rhs, temp_MIM_src)
       call scale(temp_rhs, di%MIM_options%old_immobile_saturation(p)%ptr)

       if (di%prt_is_constant) then 
         call scale(temp_rhs, di%old_porosity_cnt(1))
       else
         call scale(temp_rhs, di%old_porosity_pmesh)
       end if

       call scale(temp_rhs, di%MIM_options%immobile_prog_sfield(imf)%old_sfield)

       call set(di%MIM_options%immobile_prog_sfield(imf)%sfield, temp_rhs)
       
       !calculate '(alpha*dt*C_d)/(theta_s+alpha*dt)'
       call set(temp_rhs, temp_MIM_src)
       call scale(temp_rhs,di%MIM_options%mass_trans_coef(p)%ptr)
       call scale(temp_rhs, di%dt)
       call scale(temp_rhs, di%generic_prog_sfield(f)%sfield) ! this use the C_d value at most recent time step n+1

       call addto(di%MIM_options%immobile_prog_sfield(imf)%sfield, temp_rhs)

       call deallocate(temp_rhs)

       ewrite(1,*) 'Finished assemble and solve immobile prog sfield ',trim(di%MIM_options%immobile_prog_sfield(imf)%sfield%name),' of phase ',p

end subroutine darcy_trans_assemble_and_solve_immobile_prog_sfield


! ********* 22 Aug 2013 *******LCai ***********************************

subroutine darcy_trans_assemble_galerkin_projection_elemesh_to_pmesh(field, projected_field, positions, ele)

        type(scalar_field), intent(inout) :: field
        type(scalar_field), intent(in) :: projected_field
        type(vector_field), intent(in) :: positions
        integer, intent(in) :: ele
        type(element_type), pointer :: field_shape, proj_field_shape
        real, dimension(ele_loc(field, ele)) :: little_rhs
        real, dimension(ele_loc(field, ele), ele_loc(field, ele)) :: little_mass
        real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba
        real, dimension(ele_loc(field, ele), ele_loc(projected_field, ele)) :: little_mba_int
        real, dimension(ele_ngi(field, ele)) :: detwei
        real, dimension(ele_loc(projected_field, ele)) :: proj_field_val 
        
        integer :: i, j, k


          field_shape => ele_shape(field, ele)
          proj_field_shape => ele_shape(projected_field, ele)

          call transform_to_physical(positions, ele, detwei=detwei)

          little_mass = shape_shape(field_shape, field_shape, detwei)

          ! And compute the product of the basis functions
          little_mba = 0
          do i=1,ele_ngi(field, ele)
           forall(j=1:ele_loc(field, ele), k=1:ele_loc(projected_field, ele))
             little_mba_int(j, k) = field_shape%n(j, i) * proj_field_shape%n(k, i)
           end forall
           little_mba = little_mba + little_mba_int * detwei(i)
          end do

          proj_field_val = ele_val(projected_field, ele)
          little_rhs = matmul(little_mba, proj_field_val)

          call solve(little_mass, little_rhs)
          call set(field, ele_nodes(field, ele), little_rhs)
 

end subroutine darcy_trans_assemble_galerkin_projection_elemesh_to_pmesh

! ********************Finish*** LCai **********************************

! *****Finished **** LCai *********************************************
end module darcy_impes_assemble_module
