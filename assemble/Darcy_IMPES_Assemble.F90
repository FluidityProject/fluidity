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
   use global_parameters, only : FIELD_NAME_LEN, OPTION_PATH_LEN

   implicit none

   private

   public :: darcy_impes_type, &
             darcy_impes_assemble_and_solve, &
             darcy_impes_calculate_gradient_pressure_etc, &
             darcy_impes_calculate_phase_one_saturation_diagnostic, &
             darcy_impes_calculate_sum_saturation
   
   type darcy_impes_type
      type(vector_field_pointer), dimension(:), pointer :: inter_velocity_porosity, old_inter_velocity_porosity
      type(vector_field_pointer), dimension(:), pointer :: darcy_velocity, fractional_flow   
      type(scalar_field_pointer), dimension(:), pointer :: saturation, old_saturation
      type(scalar_field_pointer), dimension(:), pointer :: relative_permeability
      type(scalar_field_pointer), dimension(:), pointer :: viscosity
      type(scalar_field_pointer), dimension(:), pointer :: cfl, old_cfl
      type(scalar_field), pointer :: pressure, old_pressure
      type(vector_field), pointer :: gradient_pressure, old_gradient_pressure
      type(scalar_field), pointer :: porosity, old_porosity
      type(scalar_field), pointer :: absolute_permeability
      type(vector_field), pointer :: positions
      type(vector_field), pointer :: total_darcy_velocity   
      type(scalar_field), pointer :: sum_saturation, old_sum_saturation
      type(scalar_field), pointer :: div_total_darcy_velocity
      type(vector_field) :: positions_pressure_mesh
      type(vector_field) :: inter_velocity_porosity_tmp
      type(csr_matrix) :: pressure_matrix
      type(scalar_field) :: lhs, rhs, rhs_adv, rhs_time
      type(scalar_field) :: inverse_cv_mass_cfl_mesh
      type(scalar_field) :: inverse_cv_mass_pressure_mesh
      type(scalar_field) :: cv_mass_pressure_mesh_with_porosity   
      type(scalar_field) :: cv_mass_pressure_mesh_with_old_porosity 
      type(scalar_field) :: cv_mass_velocity_mesh
      type(scalar_field) :: cfl_subcycle
      type(scalar_field) :: old_saturation_subcycle   
      type(scalar_field), dimension(:), pointer :: darcy_velocity_normal_flow_bc_value
      type(scalar_field) :: total_darcy_velocity_normal_flow_bc_value
      type(mesh_type), pointer :: darcy_velocity_surface_mesh
      integer, dimension(:,:), pointer :: darcy_velocity_normal_flow_bc_flag
      integer, dimension(:), pointer :: total_darcy_velocity_normal_flow_bc_flag      
      type(csr_sparsity), pointer :: sparsity
      integer :: p, number_phase, quaddegree   
      type(cv_options_type) :: saturation_cv_options
      type(cv_faces_type) :: cvfaces
      type(element_type) :: x_cvshape_full, p_cvshape_full
      type(element_type) :: x_cvshape, p_cvshape
      type(element_type) :: x_cvbdyshape
      logical :: phase_one_saturation_diagnostic
      logical, dimension(:), pointer :: vphi_average_over_CV
      real :: saturation_max_courant_per_subcycle
      real :: dt
      integer :: ndim
      type(state_type), dimension(:), pointer :: state
   end type darcy_impes_type
   
   contains

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve(di)
      
      !!< Assemble and solve the Darcy equations using an IMPES algorithm
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p, number_subcycle
      real    :: dt_subcycle, max_cfl
      
      ewrite(1,*) 'Start Darcy IMPES assemble and solve'
      
      ! Copy to Old the CV mass on the pressure mesh with porosity 
      call set(di%cv_mass_pressure_mesh_with_old_porosity, di%cv_mass_pressure_mesh_with_porosity)
      
      ! Calculate the latest CV mass on the pressure mesh with porosity
      call compute_cv_mass(di%positions, di%cv_mass_pressure_mesh_with_porosity, di%porosity)      
            
      ! solve the pressure 
      call solve_pressure(di)
      
      ! Calculate the gradient pressure, inter_velocity_porosity ... etc
      call darcy_impes_calculate_gradient_pressure_etc(di)
            
      ! solve the saturation for each phase but not the first
      s_phase_loop: do p = 2, di%number_phase
         
         ! Deduce the number of subcycles to do and the subcycle time step size
         
         ! Find the max cfl number for this phase (accounting for parallel)
         max_cfl = maxval(di%cfl(p)%ptr)
         
         call allmax(max_cfl)
         
         ! Find the subcycle dt and cfl field (which may be used for face value limiting)
         number_subcycle = max(1,ceiling(max_cfl/di%saturation_max_courant_per_subcycle))
           
         dt_subcycle = di%dt/real(number_subcycle)

         call set(di%cfl_subcycle, di%cfl(p)%ptr)
         
         call scale(di%cfl_subcycle, 1.0/real(number_subcycle))
         
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
         
         ! Find the max cfl number for this phase (accounting for parallel)
         max_cfl = maxval(di%cfl(1)%ptr)
         
         call allmax(max_cfl)
         
         ! Find the subcycle dt and cfl field (which may be used for face value limiting)
         number_subcycle = max(1,ceiling(max_cfl/di%saturation_max_courant_per_subcycle))
           
         dt_subcycle = di%dt/real(number_subcycle)
         
         call set(di%cfl_subcycle, di%cfl(1)%ptr)
         
         call scale(di%cfl_subcycle, 1.0/real(number_subcycle))
         
         ewrite(1,*) 'Solve phase 1 Saturation with number_subcycle ',number_subcycle,' and dt_subcycle ',dt_subcycle

         call solve_phase_saturation(di, &
                                     number_subcycle, &
                                     dt_subcycle, &
                                     p = 1)

         ewrite(1,*) 'Finished solve phase 1 Saturation'
         
         ewrite_minmax(di%saturation(1)%ptr)
      
      end if s1_solve_if

      ! Calculate the sum of the saturations
      call darcy_impes_calculate_sum_saturation(di)
      
      ewrite(1,*) 'Finished Darcy IMPES assemble and solve'
           
   end subroutine darcy_impes_assemble_and_solve

! ----------------------------------------------------------------------------

   subroutine solve_pressure(di)
      
      !!< Assemble and solve the pressure
      
      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      integer :: p, vele, sele, iloc, oloc, jloc, face, gi, ggi
      real    :: absperm_vele, income, face_value, old_vphi_dot_n
      real    :: old_saturation_face_value, relperm_face_value
      logical :: inflow
      real,    dimension(1) :: dummy
      real,    dimension(:,:,:), allocatable :: old_vphi_face
      real,    dimension(:,:),   allocatable :: old_saturation_ele
      real,    dimension(:,:),   allocatable :: relperm_ele
      real,    dimension(:),     allocatable :: visc_vele
      real,    dimension(:,:),   allocatable :: x_ele
      real,    dimension(:,:,:), allocatable :: p_dshape
      real,    dimension(:,:),   allocatable :: normal
      real,    dimension(:),     allocatable :: detwei
      real,    dimension(:),     allocatable :: normgi
      logical, dimension(:),     allocatable :: notvisited
      real,    dimension(:,:),   allocatable :: p_mat_local
      real,    dimension(:,:),   allocatable :: x_face_quad
      real,    dimension(:,:),   allocatable :: minus_grad_old_p_quad
      integer, dimension(:),     pointer     :: x_pmesh_nodes
      integer, dimension(:),     pointer     :: p_nodes      
      integer, dimension(:),     pointer     :: upwind_nodes
      real,    dimension(:),     allocatable :: bc_sele_val
      real,    dimension(:,:),   allocatable :: normal_bdy
      real,    dimension(:),     allocatable :: detwei_bdy
      real,    dimension(:,:),   allocatable :: x_ele_bdy
      real,    dimension(:),     allocatable :: p_rhs_local_bdy
      integer, dimension(:),     allocatable :: p_nodes_bdy
      integer, parameter :: BC_TYPE_NORMAL_FLOW = 1
            
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      allocate(p_dshape(ele_loc(di%pressure,1), di%x_cvshape%ngi, mesh_dim(di%pressure)))
      allocate(normal(di%ndim,di%x_cvshape%ngi))
      allocate(detwei(di%x_cvshape%ngi))      
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(p_mat_local(ele_loc(di%pressure,1), ele_loc(di%pressure,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(minus_grad_old_p_quad(di%ndim, di%p_cvshape%ngi))
      allocate(visc_vele(di%number_phase))
      allocate(old_saturation_ele(di%number_phase,ele_loc(di%pressure,1)))
      allocate(relperm_ele(di%number_phase,ele_loc(di%pressure,1)))
      allocate(old_vphi_face(di%ndim,di%p_cvshape%ngi,di%number_phase))

      allocate(bc_sele_val(face_loc(di%pressure,1)))
      allocate(detwei_bdy(di%x_cvbdyshape%ngi))
      allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      allocate(p_rhs_local_bdy(face_loc(di%pressure,1)))
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(p_nodes_bdy(face_loc(di%pressure%mesh,1)))
      
      ewrite(1,*) 'Solve Pressure'
      
      ! Initialise the matrix and rhs
      call zero(di%pressure_matrix)
      call zero(di%rhs)
                 
      ! Loop volume elements assembling local contributions     
      vol_element_loop: do vele = 1,element_count(di%pressure)
         
         ! get the saturation ele values for each phase
         do p = 1,di%number_phase
            old_saturation_ele(p,:) = ele_val(di%old_saturation(p)%ptr, vele)
         end do
         
         ! get the relperm ele values for each phase
         do p = 1,di%number_phase
            relperm_ele(p,:) = ele_val(di%relative_permeability(p)%ptr, vele)
         end do
         
         ! get the viscosity value for this element for each phase
         do p = 1,di%number_phase
            dummy(1:1) = ele_val(di%viscosity(p)%ptr, vele)
            visc_vele(p) = dummy(1)
         end do
         
         ! get the absolute permeability value for this element
         dummy(1:1) = ele_val(di%absolute_permeability, vele)         
         absperm_vele = dummy(1)
         
         ! get the coordinate values for this element for each positions local node
         x_ele = ele_val(di%positions, vele)         
         
         ! get the coordinate values for this element for each quadrature point
         x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)
         
         ! The node indices of the pressure field
         p_nodes => ele_nodes(di%pressure, vele)
         
         ! The node indices of the positions projected to the pressure mesh
         x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)
         
         ! Determine the node numbers to use to determine the
         ! Saturation and RelPerm upwind values
         if((di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
            (di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then
            
            upwind_nodes => x_pmesh_nodes
        
         else
            
            upwind_nodes => p_nodes
         
         end if
         
         ! obtain the transformed determinant*weight and normals
         call transform_cvsurf_to_physical(x_ele, di%x_cvshape, &
                                           detwei, normal, di%cvfaces)
         
         ! obtain the derivative of the pressure shape function at 
         ! the CV face quadrature points
         call transform_to_physical(di%positions, vele, x_shape = di%x_cvshape_full, &
                                    shape = di%p_cvshape_full, dshape = p_dshape)
         
         ! the old intersitial velocity_porosity at the quadrature points, used to determine upwind
         ! via a FE interpolation using the pressure basis functions (which for linear is fine)
         do p = 1,di%number_phase            
            old_vphi_face(:,:,p) = ele_val_at_quad(di%old_inter_velocity_porosity(p)%ptr, vele, di%p_cvshape)         
         end do 
         
         ! Initialise array for the quadrature points of this 
         ! element for whether it has already been visited
         notvisited = .true.
         
         ! Initialise the local p matrix to assemble for this element
         p_mat_local = 0.0

         ! loop over local nodes within this element
         nodal_loop_i: do iloc = 1, di%pressure%mesh%shape%loc

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
                                      
                   ! Evaluate the face value for old saturation and old relperm for each phase
                   ! (assuming upwind for now) and evaluate the face value.

                   ! face value = sum_{phase} ( S*relperm*absperm/visc ), where absperm is phase independent
                   ! (if absperm and viscosity are considered tensors this requires modifying below)
                   
                   face_value = 0.0
                      
                   do p = 1,di%number_phase

                      ! determine if the flow is in or out of the face at this quadrature
                      ! with respect to the normal orientation using the old vphi - same as saturation assemble                     
                      old_vphi_dot_n = dot_product(old_vphi_face(:,ggi,p), normgi(:))

                      inflow = (old_vphi_dot_n<=0.0)

                      income = merge(1.0,0.0,inflow)
                                            
                      old_saturation_face_value = income*old_saturation_ele(p,oloc) + (1.0-income)*old_saturation_ele(p,iloc)

                      relperm_face_value = income*relperm_ele(p,oloc) + (1.0-income)*relperm_ele(p,iloc)

                      face_value = face_value + old_saturation_face_value * relperm_face_value / visc_vele(p)
                   
                   end do
                   
                   face_value = face_value * absperm_vele
                   
                   ! Form the local matrix given by - n_i . sum_{phase} ( S*relperm*absperm/visc ) dP/dx_j
                   do jloc = 1,di%pressure%mesh%shape%loc
                   
                      p_mat_local(iloc,jloc) = p_mat_local(iloc,jloc) - &
                                               sum(p_dshape(jloc, ggi, :)*normgi, 1)*detwei(ggi)*face_value

                      p_mat_local(oloc,jloc) = p_mat_local(oloc,jloc) - &
                                               sum(p_dshape(jloc, ggi, :)*(-normgi), 1)*detwei(ggi)*face_value
                   
                   end do
                   
                 end if check_visited
                 
               end do quadrature_loop

             end if is_neigh
             
           end do face_loop
         
         end do nodal_loop_i
         
         ! Add volume element contribution to global pressure matrix
         call addto(di%pressure_matrix, p_nodes, p_nodes, p_mat_local)
         
      end do vol_element_loop
      
      ! Add rate of change of porosity to rhs    
      call set(di%rhs, di%cv_mass_pressure_mesh_with_old_porosity)
                  
      call addto(di%rhs, di%cv_mass_pressure_mesh_with_porosity, scale = -1.0)
      
      call scale(di%rhs, 1.0/di%dt)
      
      ! Add normal darcy velocity or normal total darcy velocity BC integrals
      ! - NOTE that there cannot be individual phase darcy velocity BCs
      !        when there is a total darcy velocity BC on a surface element. 
      
      ! Get the total darcy velocity BC info
      call get_darcy_velocity_normal_flow_boundary_condition(di%total_darcy_velocity, &
                                                             di%darcy_velocity_surface_mesh, &
                                                             di%total_darcy_velocity_normal_flow_bc_value, &
                                                             di%total_darcy_velocity_normal_flow_bc_flag)
      
      ! Get the darcy velocity BC info for each phase
      do p = 1, di%number_phase
      
         call get_darcy_velocity_normal_flow_boundary_condition(di%darcy_velocity(p)%ptr, &
                                                                di%darcy_velocity_surface_mesh, &
                                                                di%darcy_velocity_normal_flow_bc_value(p), &
                                                                di%darcy_velocity_normal_flow_bc_flag(:,p))
      
      end do
      
      sele_loop: do sele = 1, surface_element_count(di%pressure)
         
         ! If total darcy BC then apply, else apply individual phase darcy BC
         ! and where there is no velocity BC include necessary pressure integral.
                 
         have_total_darcy_vel_bc: if (di%total_darcy_velocity_normal_flow_bc_flag(sele) == BC_TYPE_NORMAL_FLOW) then
            
            ! Check that there are no individual darcy velocity phase BC
            do p = 1, di%number_phase
               if (di%darcy_velocity_normal_flow_bc_flag(sele,p) == BC_TYPE_NORMAL_FLOW) then
                  ewrite(-1,*) 'Issue with DarcyVelocity normal_flow BC for phase ',p
                  FLExit('Cannot apply individual phase DarcyVelocity normal_flow BC when one is applied to the TotalDarcyVelocity')
               end if
            end do

            bc_sele_val = ele_val(di%total_darcy_velocity_normal_flow_bc_value, sele)

            x_ele       = ele_val(di%positions, face_ele(di%positions, sele))
            x_ele_bdy   = face_val(di%positions, sele)
            p_nodes_bdy = face_global_nodes(di%pressure%mesh, sele)

            call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)

            p_rhs_local_bdy = 0.0

            total_bc_iloc_loop: do iloc = 1, di%pressure%mesh%faces%shape%loc

               total_bc_face_loop: do face = 1, di%cvfaces%sfaces

                  total_bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     total_bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi

                        p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - bc_sele_val(iloc)* detwei_bdy(ggi)                    

                     end do total_bc_quad_loop

                  end if total_bc_neigh_if

               end do total_bc_face_loop

            end do total_bc_iloc_loop

            call addto(di%rhs, p_nodes_bdy, p_rhs_local_bdy)
             
         else have_total_darcy_vel_bc
         
            loop_phase_darcy_vel_bc: do p = 1, di%number_phase
               
               have_phase_darcy_vel_bc: if (di%darcy_velocity_normal_flow_bc_flag(sele,p) == BC_TYPE_NORMAL_FLOW) then

                  bc_sele_val = ele_val(di%darcy_velocity_normal_flow_bc_value(p), sele)

                  x_ele       = ele_val(di%positions, face_ele(di%positions, sele))
                  x_ele_bdy   = face_val(di%positions, sele)
                  p_nodes_bdy = face_global_nodes(di%pressure%mesh, sele)

                  call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)

                  p_rhs_local_bdy = 0.0
                  
                  bc_iloc_loop: do iloc = 1, di%pressure%mesh%faces%shape%loc

                     bc_face_loop: do face = 1, di%cvfaces%sfaces

                        bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                           bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                              ggi = (face-1)*di%cvfaces%shape%ngi + gi

                              p_rhs_local_bdy(iloc) = p_rhs_local_bdy(iloc) - bc_sele_val(iloc)* detwei_bdy(ggi)                    

                           end do bc_quad_loop

                        end if bc_neigh_if

                     end do bc_face_loop

                  end do bc_iloc_loop

                  call addto(di%rhs, p_nodes_bdy, p_rhs_local_bdy)                  
                  
               else have_phase_darcy_vel_bc
                  
                  ! ***** the integral of the phase coeff * grad_p over the boundary cv surface requires coding here *****
                  
                  
                  
               end if have_phase_darcy_vel_bc
            
            end do loop_phase_darcy_vel_bc
         
         end if have_total_darcy_vel_bc
                           
      end do sele_loop      
      
      ! Apply any strong dirichlet BC's
      call apply_dirichlet_conditions(di%pressure_matrix, di%rhs, di%pressure)
      
      ! Solve the pressure
      call petsc_solve(di%pressure, di%pressure_matrix, di%rhs, di%state(1))
      
      ! deallocate local variables as required
      deallocate(x_ele)
      deallocate(p_dshape)
      deallocate(normal)
      deallocate(detwei)
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(p_mat_local)
      deallocate(x_face_quad)
      deallocate(minus_grad_old_p_quad)
      deallocate(visc_vele)
      deallocate(old_saturation_ele)
      deallocate(relperm_ele)
      deallocate(old_vphi_face)
      deallocate(bc_sele_val)
      deallocate(detwei_bdy)
      deallocate(normal_bdy)
      deallocate(p_rhs_local_bdy)
      deallocate(x_ele_bdy)      
      deallocate(p_nodes_bdy)
      
      ewrite(1,*) 'Finished solve Pressure'
      
      ewrite_minmax(di%pressure)
      
   end subroutine solve_pressure

! ----------------------------------------------------------------------------

   subroutine darcy_impes_calculate_gradient_pressure_etc(di)
      
      !!< Calculate the gradient pressure, inter_velocity_porosity, CFL and darcy velocity fields

      type(darcy_impes_type), intent(inout) :: di
      
      ! local variables
      logical :: inflow
      integer :: vele, p, nloc, dim, iloc, oloc, face, gi, ggi, sele
      real :: visc_ele, absperm_ele, darcy_vel_face_value_dot_n
      real    :: income, old_vphi_dot_n, vphi_face_value_dot_n
      real, dimension(1) :: dummy
      real, dimension(di%ndim) :: grad_pressure_ele
      real, dimension(:), allocatable :: relperm_ele
      real, dimension(:), allocatable :: inter_velocity_porosity_local
      real,    dimension(:,:),   allocatable :: old_vphi_face
      real,    dimension(:,:),   allocatable :: x_ele
      real,    dimension(:,:),   allocatable :: normal
      real,    dimension(:),     allocatable :: detwei
      real,    dimension(:),     allocatable :: normgi
      logical, dimension(:),     allocatable :: notvisited
      real,    dimension(:),     allocatable :: cfl_rhs_local
      real,    dimension(:),     allocatable :: div_tvphi_rhs_local
      real,    dimension(:,:),   allocatable :: x_face_quad
      real,    dimension(:,:),   allocatable :: vphi_ele
      real,    dimension(:),     allocatable :: vphi_face_value
      real,    dimension(:,:),   allocatable :: darcy_vel_ele
      real,    dimension(:),     allocatable :: darcy_vel_face_value
      integer, dimension(:),     pointer     :: x_pmesh_nodes
      integer, dimension(:),     pointer     :: p_nodes      
      integer, dimension(:),     pointer     :: upwind_nodes
      real,    dimension(:,:),   allocatable :: normal_bdy
      real,    dimension(:),     allocatable :: detwei_bdy
      real,    dimension(:,:),   allocatable :: x_ele_bdy
      real,    dimension(:),     allocatable :: cfl_rhs_local_bdy
      real,    dimension(:),     allocatable :: div_tvphi_rhs_local_bdy
      integer, dimension(:),     allocatable :: p_nodes_bdy
      real,    dimension(:,:),   allocatable :: vphi_sele
      real,    dimension(:,:),   allocatable :: darcy_vel_sele
      
      ewrite(1,*) 'Calculate GradientPressure'
      
      ! Calculate the gradient pressure field assumed DG
      call grad(di%pressure, &
                di%positions, &
                di%gradient_pressure)
      
      ewrite_minmax(di%gradient_pressure)
      
      ewrite(1,*) 'Calculate InterstitialVelocityPorosity'
      
      ! Calculate the inter_velocity_porosity field = - (1.0/old_sigma) * grad_pressure, where sigma = visc / absperm * relperm
      
      allocate(relperm_ele(ele_loc(di%pressure,1)))
      
      ! deduce the number of local nodes (all elements assumed same) for inter_velocity_porosity and CFL
      nloc = di%inter_velocity_porosity(1)%ptr%mesh%shape%loc

      allocate(inter_velocity_porosity_local(nloc))      
      
      do vele = 1, element_count(di%pressure)
         
         ! Find the element wise abs perm, porosity and grad_pressure
         ! (a dummy is used such as to not assume that the node number
         !  is the same as the element number)
         
         dummy = ele_val(di%absolute_permeability, vele)
         absperm_ele = dummy(1)
                  
         do dim = 1,di%ndim         
            dummy = ele_val(di%gradient_pressure, dim, vele)
            grad_pressure_ele(dim) = dummy(1)         
         end do
         
         ! Form the inter_velocity_porosity and CFL for each phase
         do p = 1,di%number_phase
         
            ! Find the element wise viscosity
            dummy = ele_val(di%viscosity(p)%ptr, vele)
            visc_ele = dummy(1)
            
            ! Find the element wise local values for relperm
            relperm_ele = ele_val(di%relative_permeability(p)%ptr, vele)
            
            ! find the local DG values for inter_velocity_porosity for each geometric dimension
            ! (if absperm and viscosity are considered tensors this requires modifying below)
            do dim = 1,di%ndim
               
               inter_velocity_porosity_local(:) = - relperm_ele(:) * absperm_ele * grad_pressure_ele(dim) / visc_ele
            
               call set(di%inter_velocity_porosity(p)%ptr, &
                        dim, &
                        ele_nodes(di%inter_velocity_porosity(p)%ptr, vele), &
                        inter_velocity_porosity_local)
                           
            end do
            
         end do
         
      end do
      
      ! Average the interstitalvelocity_porosity over the CV if required
      do p = 1, di%number_phase
         
         if (di%vphi_average_over_CV(p)) then
         
            ewrite(1,*) 'Average InterstitialVelocityPorosity over CV for phase ',p
            
            call set(di%inter_velocity_porosity_tmp, &
                     di%inter_velocity_porosity(p)%ptr)
                        
            allocate(vphi_ele(di%ndim, ele_loc(di%inter_velocity_porosity(p)%ptr,1)))

            call zero(di%inter_velocity_porosity(p)%ptr)

            do vele = 1, element_count(di%pressure)

               vphi_ele = ele_val(di%inter_velocity_porosity_tmp,vele)

               do dim = 1, di%ndim

                  call addto(di%inter_velocity_porosity(p)%ptr, &
                             dim, &
                             ele_nodes(di%inter_velocity_porosity(p)%ptr,vele), &
                             vphi_ele(dim,:)*&
                            &ele_val(di%cv_mass_velocity_mesh,vele)*&
                            &ele_val(di%inverse_cv_mass_pressure_mesh,vele))

               end do

            end do 

            deallocate(vphi_ele)
            
         end if
         
      end do
      
      do p = 1,di%number_phase
         ewrite_minmax(di%inter_velocity_porosity(p)%ptr)
      end do
      
      ! calculate the darcy velocity = inter_velocity_porosity * old_saturation
      do p = 1, di%number_phase         
         
         ewrite(1,*) 'Calculate DarcyVelocity for phase ',p
         
         call zero(di%darcy_velocity(p)%ptr)
         
         do vele = 1, element_count(di%pressure)
            
            do dim = 1, di%ndim
            
               call addto(di%darcy_velocity(p)%ptr, &
                          dim, &
                          ele_nodes(di%darcy_velocity(p)%ptr, vele), &
                          ele_val(di%old_saturation(p)%ptr, vele))
            
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
      
      ewrite(1,*) 'Calculate the DivergenceTotalDarcyVelocity'
      
      ! calculate the divergence of the total darcy velocity
      ! and the CFL number of each phase using CV surface integrals

      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      allocate(normal(di%ndim,di%x_cvshape%ngi))
      allocate(detwei(di%x_cvshape%ngi))      
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(cfl_rhs_local(ele_loc(di%pressure,1)))
      allocate(div_tvphi_rhs_local(ele_loc(di%pressure,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(vphi_ele(di%ndim,ele_loc(di%inter_velocity_porosity(1)%ptr,1)))
      allocate(darcy_vel_ele(di%ndim,ele_loc(di%darcy_velocity(1)%ptr,1)))
      allocate(old_vphi_face(di%ndim,di%p_cvshape%ngi))
      allocate(vphi_face_value(di%ndim))
      allocate(darcy_vel_face_value(di%ndim))
      
      allocate(detwei_bdy(di%x_cvbdyshape%ngi))
      allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      allocate(cfl_rhs_local_bdy(face_loc(di%pressure,1)))
      allocate(div_tvphi_rhs_local_bdy(face_loc(di%pressure,1)))      
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(vphi_sele(di%ndim,face_loc(di%inter_velocity_porosity(1)%ptr,1)))
      allocate(darcy_vel_sele(di%ndim,face_loc(di%darcy_velocity(1)%ptr,1)))
      allocate(p_nodes_bdy(face_loc(di%pressure%mesh,1)))

      call zero(di%div_total_darcy_velocity)

      do p = 1, di%number_phase
         
         call zero(di%cfl(p)%ptr)
               
         do vele = 1, element_count(di%pressure)

            ! The interstitial velocity porosity ele value
            vphi_ele(:,:) = ele_val(di%inter_velocity_porosity(p)%ptr, vele)

            ! The darcy velocity ele value
            darcy_vel_ele(:,:) = ele_val(di%darcy_velocity(p)%ptr, vele)

            ! the old intersitial velocity_porosity at the quadrature points, used to determine upwind
            ! via a FE interpolation using the pressure basis functions (which for linear is fine)
            old_vphi_face(:,:) = ele_val_at_quad(di%old_inter_velocity_porosity(p)%ptr, vele, di%p_cvshape)         

            ! get the coordinate values for this element for each positions local node
            x_ele = ele_val(di%positions, vele)         

            ! get the coordinate values for this element for each quadrature point
            x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)

            ! The node indices of the pressure field
            p_nodes => ele_nodes(di%pressure, vele)

            ! The node indices of the positions projected to the pressure mesh
            x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)

            ! Determine the node numbers to use to determine the upwind values
            if((di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
               (di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then

               upwind_nodes => x_pmesh_nodes

            else

               upwind_nodes => p_nodes

            end if

            ! obtain the transformed determinant*weight and normals
            call transform_cvsurf_to_physical(x_ele, di%x_cvshape, &
                                              detwei, normal, di%cvfaces)

            ! Initialise array for the quadrature points of this 
            ! element for whether it has already been visited
            notvisited = .true.

            ! Initialise the local rhs's to assemble for this element
            cfl_rhs_local       = 0.0
            div_tvphi_rhs_local = 0.0

            ! loop over local nodes within this element
            nodal_loop_i: do iloc = 1, di%pressure%mesh%shape%loc

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
                      ! with respect to the normal orientation using the old vphi - same as pressure assemble
                      old_vphi_dot_n = dot_product(old_vphi_face(:,ggi), normgi(:))

                      inflow = (old_vphi_dot_n<=0.0)

                      income = merge(1.0,0.0,inflow)

                      ! Evaluate the face value for interstitial_velocity_porosity (assuming upwind for now)

                      vphi_face_value(:) = income*vphi_ele(:,oloc) + (1.0-income)*vphi_ele(:,iloc)

                      vphi_face_value_dot_n = dot_product(vphi_face_value, normgi)

                      darcy_vel_face_value(:) = income*darcy_vel_ele(:,oloc) + (1.0-income)*darcy_vel_ele(:,iloc)

                      darcy_vel_face_value_dot_n = dot_product(darcy_vel_face_value, normgi)

                      div_tvphi_rhs_local(iloc) = div_tvphi_rhs_local(iloc) + darcy_vel_face_value_dot_n * detwei(ggi)

                      div_tvphi_rhs_local(oloc) = div_tvphi_rhs_local(oloc) - darcy_vel_face_value_dot_n * detwei(ggi)

                      cfl_rhs_local(iloc) = cfl_rhs_local(iloc) + abs(vphi_face_value_dot_n) * detwei(ggi) * (1.0 - income)

                      cfl_rhs_local(oloc) = cfl_rhs_local(oloc) + abs(vphi_face_value_dot_n) * detwei(ggi) * income

                    end if check_visited

                  end do quadrature_loop

                end if is_neigh

              end do face_loop

            end do nodal_loop_i

            call addto(di%div_total_darcy_velocity, p_nodes, div_tvphi_rhs_local)

            call addto(di%cfl(p)%ptr, p_nodes, cfl_rhs_local)

         end do
         
         sele_loop: do sele = 1, surface_element_count(di%pressure)

            vphi_sele(:,:) = face_val(di%inter_velocity_porosity(p)%ptr, sele)

            darcy_vel_sele(:,:) = face_val(di%darcy_velocity(p)%ptr, sele)

            x_ele       = ele_val(di%positions, face_ele(di%positions, sele))
            x_ele_bdy   = face_val(di%positions, sele)
            p_nodes_bdy = face_global_nodes(di%pressure%mesh, sele)

            call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)

            ! Initialise the local rhs's to assemble for this element
            cfl_rhs_local_bdy       = 0.0
            div_tvphi_rhs_local_bdy = 0.0

            bc_iloc_loop: do iloc = 1, di%pressure%mesh%faces%shape%loc

               bc_face_loop: do face = 1, di%cvfaces%sfaces

                  bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                     bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                        ggi = (face-1)*di%cvfaces%shape%ngi + gi

                        vphi_face_value_dot_n = dot_product(vphi_sele(:,iloc), normal_bdy(:,ggi))

                        darcy_vel_face_value_dot_n = dot_product(darcy_vel_sele(:,iloc), normal_bdy(:,ggi))

                        div_tvphi_rhs_local_bdy(iloc) = div_tvphi_rhs_local_bdy(iloc) + &
                                                       &darcy_vel_face_value_dot_n * detwei_bdy(ggi)                   

                        cfl_rhs_local_bdy(iloc) = cfl_rhs_local_bdy(iloc) + &
                                                 &abs(vphi_face_value_dot_n) * detwei_bdy(ggi) * (1.0 - income)                     

                     end do bc_quad_loop

                  end if bc_neigh_if

               end do bc_face_loop

            end do bc_iloc_loop

            call addto(di%cfl(p)%ptr, p_nodes_bdy, cfl_rhs_local_bdy)

         end do sele_loop
                  
         di%cfl(p)%ptr%val = di%cfl(p)%ptr%val * di%dt / di%cv_mass_pressure_mesh_with_porosity%val
         
      end do
      
      call scale(di%div_total_darcy_velocity, di%inverse_cv_mass_pressure_mesh)
      
      ewrite_minmax(di%div_total_darcy_velocity)
      
      ! deallocate local variables as required
      deallocate(x_ele)
      deallocate(normal)
      deallocate(detwei)
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(cfl_rhs_local)
      deallocate(div_tvphi_rhs_local)      
      deallocate(x_face_quad)
      deallocate(vphi_ele)
      deallocate(darcy_vel_ele)
      deallocate(old_vphi_face)
      deallocate(vphi_face_value)
      deallocate(darcy_vel_face_value)

      deallocate(detwei_bdy)
      deallocate(normal_bdy)
      deallocate(x_ele_bdy)
      deallocate(p_nodes_bdy)
      deallocate(cfl_rhs_local_bdy)
      deallocate(div_tvphi_rhs_local_bdy)      
      deallocate(vphi_sele)      
      deallocate(darcy_vel_sele)      
            
      deallocate(relperm_ele)
      deallocate(inter_velocity_porosity_local)
      
   end subroutine darcy_impes_calculate_gradient_pressure_etc

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

      ! local variables
      integer :: vele, iloc, oloc, face, gi, ggi, sele, isub
      real    :: income, face_value, old_vphi_dot_n, alpha_start, alpha_end
      real    :: old_saturation_face_value, vphi_face_value_dot_n
      logical :: inflow
      real,    dimension(:,:),   allocatable :: old_vphi_face
      real,    dimension(:),     allocatable :: old_saturation_ele
      real,    dimension(:,:),   allocatable :: x_ele
      real,    dimension(:,:),   allocatable :: normal
      real,    dimension(:),     allocatable :: detwei
      real,    dimension(:),     allocatable :: normgi
      logical, dimension(:),     allocatable :: notvisited
      real,    dimension(:),     allocatable :: s_rhs_local
      real,    dimension(:,:),   allocatable :: x_face_quad
      real,    dimension(:,:),   allocatable :: vphi_ele
      real,    dimension(:),     allocatable :: vphi_face_value
      integer, dimension(:),     pointer     :: x_pmesh_nodes
      integer, dimension(:),     pointer     :: p_nodes      
      integer, dimension(:),     pointer     :: upwind_nodes
      real,    dimension(:),     allocatable :: old_saturation_ele_bdy
      real,    dimension(:,:),   allocatable :: normal_bdy
      real,    dimension(:),     allocatable :: detwei_bdy
      real,    dimension(:,:),   allocatable :: x_ele_bdy
      real,    dimension(:),     allocatable :: s_rhs_local_bdy
      integer, dimension(:),     allocatable :: p_nodes_bdy
      real,    dimension(:,:),   allocatable :: vphi_sele
      integer, dimension(:),     allocatable :: saturation_bc_type
      type(scalar_field) :: saturation_bc
      integer, parameter :: BC_TYPE_WEAKDIRICHLET = 1, BC_TYPE_ZERO_FLUX = 2, BC_TYPE_DIRICHLET = 3
         
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(di%ndim, ele_loc(di%positions,1)))      
      allocate(normal(di%ndim,di%x_cvshape%ngi))
      allocate(detwei(di%x_cvshape%ngi))      
      allocate(normgi(di%ndim))
      allocate(notvisited(di%x_cvshape%ngi))
      allocate(s_rhs_local(ele_loc(di%pressure,1)))
      allocate(x_face_quad(di%ndim, di%x_cvshape%ngi))
      allocate(old_saturation_ele(ele_loc(di%pressure,1)))
      allocate(vphi_ele(di%ndim,ele_loc(di%inter_velocity_porosity(p)%ptr,1)))
      allocate(old_vphi_face(di%ndim,di%p_cvshape%ngi))
      allocate(vphi_face_value(di%ndim))
      
      allocate(old_saturation_ele_bdy(face_loc(di%pressure,1)))
      allocate(detwei_bdy(di%x_cvbdyshape%ngi))
      allocate(normal_bdy(di%ndim, di%x_cvbdyshape%ngi))
      allocate(s_rhs_local_bdy(face_loc(di%pressure,1)))
      allocate(x_ele_bdy(di%ndim, face_loc(di%positions,1)))      
      allocate(vphi_sele(di%ndim,face_loc(di%inter_velocity_porosity(p)%ptr,1)))
      allocate(p_nodes_bdy(face_loc(di%pressure%mesh,1)))
      
      ! Inititalise rhs advection field
      call zero(di%rhs_adv)
                 
      ! Loop volume elements assembling local contributions    
      vol_element_loop: do vele = 1,element_count(di%pressure)
         
         ! get the saturation ele values
         old_saturation_ele(:) = ele_val(di%old_saturation(p)%ptr, vele)

         ! The interstitial velocity porosity ele value
         vphi_ele(:,:) = ele_val(di%inter_velocity_porosity(p)%ptr, vele)

         ! the old intersitial velocity_porosity at the quadrature points, used to determine upwind
         ! via a FE interpolation using the pressure basis functions (which for linear is fine)
         old_vphi_face(:,:) = ele_val_at_quad(di%old_inter_velocity_porosity(p)%ptr, vele, di%p_cvshape)         
                  
         ! get the coordinate values for this element for each positions local node
         x_ele = ele_val(di%positions, vele)         
         
         ! get the coordinate values for this element for each quadrature point
         x_face_quad = ele_val_at_quad(di%positions, vele, di%x_cvshape)
         
         ! The node indices of the pressure field
         p_nodes => ele_nodes(di%pressure, vele)
         
         ! The node indices of the positions projected to the pressure mesh
         x_pmesh_nodes => ele_nodes(di%positions_pressure_mesh, vele)
         
         ! Determine the node numbers to use to determine the upwind values
         if((di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
            (di%saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then
            
            upwind_nodes => x_pmesh_nodes
        
         else
            
            upwind_nodes => p_nodes
         
         end if
         
         ! obtain the transformed determinant*weight and normals
         call transform_cvsurf_to_physical(x_ele, di%x_cvshape, &
                                           detwei, normal, di%cvfaces)
                  
         ! Initialise array for the quadrature points of this 
         ! element for whether it has already been visited
         notvisited = .true.
         
         ! Initialise the local rhs to assemble for this element
         s_rhs_local = 0.0

         ! loop over local nodes within this element
         nodal_loop_i: do iloc = 1, di%pressure%mesh%shape%loc

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
                   ! with respect to the normal orientation using the old vphi - same as pressure assemble
                   old_vphi_dot_n = dot_product(old_vphi_face(:,ggi), normgi(:))
                                      
                   inflow = (old_vphi_dot_n<=0.0)
                   
                   income = merge(1.0,0.0,inflow)

                   ! Evaluate the face value for old saturation and interstitial_velocity_porosity
                   ! (assuming upwind for now) and evaluate the full face value.
                                                                  
                   old_saturation_face_value = income*old_saturation_ele(oloc) + (1.0-income)*old_saturation_ele(iloc)
                   
                   vphi_face_value(:) = income*vphi_ele(:,oloc) + (1.0-income)*vphi_ele(:,iloc)

                   vphi_face_value_dot_n = dot_product(vphi_face_value, normgi)

                   face_value = old_saturation_face_value * vphi_face_value_dot_n * detwei(ggi)
                                                         
                   ! Form the local rhs for iloc
                   s_rhs_local(iloc) = s_rhs_local(iloc) - face_value
                   
                   ! Also form contribution to opposing local node with sign change due to reverse normal vector
                   s_rhs_local(oloc) = s_rhs_local(oloc) + face_value
                   
                 end if check_visited
                 
               end do quadrature_loop

             end if is_neigh
             
           end do face_loop
         
         end do nodal_loop_i
         
         ! Add volume element contribution to global rhs advection field
         call addto(di%rhs_adv, p_nodes, s_rhs_local)
         
      end do vol_element_loop
      
      ! Add BC integrals

      allocate(saturation_bc_type(surface_element_count(di%saturation(p)%ptr)))         
      
      call get_entire_boundary_condition(di%saturation(p)%ptr, (/"zero_flux", &
                                                                 "dirichlet"/), saturation_bc, saturation_bc_type)
     
      sele_loop: do sele = 1, surface_element_count(di%saturation(p)%ptr)
         
         if (saturation_bc_type(sele) == BC_TYPE_ZERO_FLUX) cycle
         
         if (saturation_bc_type(sele) == BC_TYPE_DIRICHLET) cycle
         
         old_saturation_ele_bdy = face_val(di%old_saturation(p)%ptr, sele)
                           
         vphi_sele(:,:) = face_val(di%inter_velocity_porosity(p)%ptr, sele)
         
         x_ele       = ele_val(di%positions, face_ele(di%positions, sele))
         x_ele_bdy   = face_val(di%positions, sele)
         p_nodes_bdy = face_global_nodes(di%pressure%mesh, sele)
         
         call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, di%x_cvbdyshape, normal_bdy, detwei_bdy)

         ! Initialise the local rhs to assemble for this element
         s_rhs_local_bdy = 0.0

         bc_iloc_loop: do iloc = 1, di%pressure%mesh%faces%shape%loc

            bc_face_loop: do face = 1, di%cvfaces%sfaces

               bc_neigh_if: if(di%cvfaces%sneiloc(iloc,face)/=0) then

                  bc_quad_loop: do gi = 1, di%cvfaces%shape%ngi

                     ggi = (face-1)*di%cvfaces%shape%ngi + gi

                     vphi_face_value_dot_n = dot_product(vphi_sele(:,iloc), normal_bdy(:,ggi))
                   
                     face_value = old_saturation_ele_bdy(iloc) * vphi_face_value_dot_n * detwei_bdy(ggi)
                     
                     s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - face_value                     
                     
                  end do bc_quad_loop

               end if bc_neigh_if

            end do bc_face_loop

         end do bc_iloc_loop

         call addto(di%rhs_adv, p_nodes_bdy, s_rhs_local_bdy)
                  
      end do sele_loop

      ! Solve the saturation for each subcycle via:            
      !  - Form the lhs = cv_mass_pressure_mesh_with_porosity / dt
      !  - and rhs = rhs_adv + rhs_time, where rhs_time = cv_mass_pressure_mesh_with_old_porosity * old_saturation / dt
      !  - solve for latest saturation and copy to start subcycle saturation step
      
      ! Note: the porosity at the start and end of a subcycle time step
      ! are linearly interpolated values from the main time step start and end
      
      call set(di%old_saturation_subcycle, di%old_saturation(p)%ptr)

      sub_loop: do isub = 1,number_subcycle
         
         ! Form the start and end of subcycle dt
         ! porosity linear interpolents
         alpha_start = (isub - 1) / number_subcycle
         alpha_end   = isub / number_subcycle
         
         call set(di%lhs, di%cv_mass_pressure_mesh_with_porosity)
         
         call scale(di%lhs, alpha_end)
         
         call addto(di%lhs, di%cv_mass_pressure_mesh_with_old_porosity, scale = (1.0 - alpha_end))

         call scale(di%lhs, 1.0/dt_subcycle)

         call set(di%rhs_time, di%cv_mass_pressure_mesh_with_porosity)
         
         call scale(di%rhs_time, alpha_start)
         
         call addto(di%rhs_time, di%cv_mass_pressure_mesh_with_old_porosity, scale = (1.0 - alpha_start))
         
         call scale(di%rhs_time, 1.0/dt_subcycle)

         ! add rhs time term
         call set(di%rhs, di%rhs_time)         
         call scale(di%rhs, di%old_saturation_subcycle)

         ! add rhs advection term
         call addto(di%rhs, di%rhs_adv)

         ! apply strong dirichlet BC
         call apply_dirichlet_conditions(di%lhs, di%rhs, di%saturation(p)%ptr)

         ! Solve for the saturation
         di%saturation(p)%ptr%val = di%rhs%val / di%lhs%val
         
         ! Copy latest solution to old subcycle
         call set(di%old_saturation_subcycle, di%saturation(p)%ptr)
         
      end do sub_loop
           
      ! deallocate local variables as required
      deallocate(x_ele)
      deallocate(normal)
      deallocate(detwei)
      deallocate(normgi)
      deallocate(notvisited)
      deallocate(s_rhs_local)
      deallocate(x_face_quad)
      deallocate(old_saturation_ele)
      deallocate(vphi_ele)
      deallocate(old_vphi_face)
      deallocate(vphi_face_value)

      deallocate(old_saturation_ele_bdy)
      deallocate(detwei_bdy)
      deallocate(normal_bdy)
      deallocate(x_ele_bdy)
      deallocate(p_nodes_bdy)
      deallocate(s_rhs_local_bdy)
      deallocate(vphi_sele)      
      deallocate(saturation_bc_type)
      call deallocate(saturation_bc)
      
   end subroutine solve_phase_saturation

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

  subroutine get_darcy_velocity_normal_flow_boundary_condition(darcy_velocity, &
                                                               darcy_velocity_surface_mesh, &
                                                               darcy_velocity_normal_flow_bc_value, &
                                                               darcy_velocity_normal_flow_bc_flag)
    
    !!< Form the data associated with darcy velocity normal_flow BC by returning 
    !!< full surface mesh arrays of a flag indicating normal_flow and the value.
    !!< This can be used for the phase DarcyVelocity and the TotalDarcyVelocity.
    
    type(vector_field),               intent(in),   target :: darcy_velocity
    type(mesh_type),                  intent(in)           :: darcy_velocity_surface_mesh
    type(scalar_field),               intent(inout)        :: darcy_velocity_normal_flow_bc_value
    integer,            dimension(:), intent(inout)        :: darcy_velocity_normal_flow_bc_flag
    
    ! Local variables
    type(scalar_field),                         pointer :: scalar_surface_field
    integer,                      dimension(:), pointer :: surface_element_list
    integer                                             :: i, k, sele
    character(len=FIELD_NAME_LEN)                       :: bctype
    
    ewrite(1,*) 'Get ',trim(darcy_velocity%name),' normal_flow boundary condition data'
    
    ! Zero the normal_flow value surface field on whole boundary mesh  
    call zero(darcy_velocity_normal_flow_bc_value)
    
    ! Initialise flag for whether surface element has normal_flow BC
    darcy_velocity_normal_flow_bc_flag = 0

    ! Loop each BC object instance for the darcy_velocity
    ! (May have multiple normal_flow BC's applied to different surface id's)
    BC_loop: do i=1, get_boundary_condition_count(darcy_velocity)
       
       ! Get this BC info
       call get_boundary_condition(darcy_velocity, i, type = bctype, &
                                   surface_element_list = surface_element_list)
          
       ! check this is a normal_flow BC (nothing else is permitted)
       if (trim(bctype) /= 'normal_flow') then
          FLAbort('Have unknown BC type for either DarcyVelocity or TotalDarcyVelocity')
       end if
       
       ! Extract the scalar_surface_field for this BC 
       if (associated(darcy_velocity%bc%boundary_condition(i)%scalar_surface_fields)) then
          scalar_surface_field => darcy_velocity%bc%boundary_condition(i)%scalar_surface_fields(1)
       else
          FLAbort('Component scalar_surface_fields for DarcyVelocity or TotalDarcyVelocity BC type not associated')
       end if
       
       ! Loop the surface elements associated with this BC instance
       ! and place the required BC value in whole boundary field
       BC_sele_loop: do k=1, size(surface_element_list)
          
          ! Find the whole domain surface element number
          sele = surface_element_list(k)
          
          ! Check that there is only 1 BC applied per surface element
          if (darcy_velocity_normal_flow_bc_flag(sele) /= 0) then             
             FLExit('Cannot apply more than 1 BC to a surface element for DarcyVelocity and TotalDarcyVelocity')
          end if
          
          ! Set the sele flag to indicate it has a normal_flow BC
          darcy_velocity_normal_flow_bc_flag(sele) = 1
          
          ! Set the normal_flow field values from this BC for its sele
          call set(darcy_velocity_normal_flow_bc_value, &
                   ele_nodes(darcy_velocity_surface_mesh, sele), &
                   ele_val(scalar_surface_field, k))
          
       end do BC_sele_loop
    
    end do BC_loop
    
    ewrite_minmax(darcy_velocity_normal_flow_bc_value)
    
    ewrite(1,*) 'Finished get ',trim(darcy_velocity%name),' normal_flow boundary condition data'

  end subroutine get_darcy_velocity_normal_flow_boundary_condition
  
! ----------------------------------------------------------------------------

end module darcy_impes_assemble_module
