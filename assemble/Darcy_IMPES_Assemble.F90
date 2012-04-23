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

   public :: darcy_impes_assemble_and_solve, &
             darcy_impes_calculate_gradient_pressure_etc, &
             darcy_impes_calculate_phase_one_saturation_diagnostic, &
             darcy_impes_calculate_sum_saturation

   contains

! ----------------------------------------------------------------------------

   subroutine darcy_impes_assemble_and_solve(pressure, &
                                             old_pressure, &
                                             pressure_matrix, &
                                             lhs, &
                                             rhs, &
                                             rhs_adv, &
                                             rhs_time, &
                                             positions, &
                                             positions_pressure_mesh, &
                                             porosity, &
                                             old_porosity, &
                                             absolute_permeability, &
                                             saturation, &
                                             old_saturation, &
                                             old_saturation_subcycle, &
                                             sum_saturation, &
                                             old_sum_saturation, &
                                             relative_permeability, &
                                             viscosity, &
                                             inter_velocity_porosity, &
                                             old_inter_velocity_porosity, &
                                             gradient_pressure, &
                                             old_gradient_pressure, &
                                             darcy_velocity, &
                                             total_darcy_velocity, &
                                             div_total_darcy_velocity, &
                                             inverse_cv_mass_cfl_mesh, &
                                             inverse_cv_mass_pressure_mesh, &
                                             cv_mass_pressure_mesh_with_porosity, &
                                             cv_mass_pressure_mesh_with_old_porosity, &
                                             cfl, &
                                             old_cfl, &
                                             cfl_subcycle, &
                                             saturation_max_courant_per_subcycle, &
                                             saturation_cv_options, &
                                             cvfaces, &
                                             x_cvshape_full, &
                                             p_cvshape_full, &
                                             x_cvshape, &
                                             p_cvshape, &
                                             x_cvbdyshape, &                                 
                                             dt, &
                                             number_phase, &
                                             phase_one_saturation_diagnostic, &
                                             state_one)
      
      !!< Assemble and solve the Darcy equations using an IMPES algorithm
      
      type(scalar_field),                       intent(inout) :: pressure
      type(scalar_field),                       intent(in)    :: old_pressure
      type(csr_matrix),                         intent(inout) :: pressure_matrix
      type(scalar_field),                       intent(inout) :: lhs
      type(scalar_field),                       intent(inout) :: rhs
      type(scalar_field),                       intent(inout) :: rhs_adv
      type(scalar_field),                       intent(inout) :: rhs_time 
      type(vector_field),                       intent(in)    :: positions
      type(vector_field),                       intent(in)    :: positions_pressure_mesh
      type(scalar_field),                       intent(in)    :: porosity
      type(scalar_field),                       intent(in)    :: old_porosity
      type(scalar_field),                       intent(in)    :: absolute_permeability
      type(scalar_field_pointer), dimension(:), intent(inout) :: saturation
      type(scalar_field_pointer), dimension(:), intent(in)    :: old_saturation
      type(scalar_field),                       intent(inout) :: old_saturation_subcycle 
      type(scalar_field),                       intent(inout) :: sum_saturation
      type(scalar_field),                       intent(in)    :: old_sum_saturation      
      type(scalar_field_pointer), dimension(:), intent(in)    :: relative_permeability
      type(scalar_field_pointer), dimension(:), intent(in)    :: viscosity
      type(vector_field_pointer), dimension(:), intent(inout) :: inter_velocity_porosity
      type(vector_field_pointer), dimension(:), intent(in)    :: old_inter_velocity_porosity
      type(vector_field),                       intent(inout) :: gradient_pressure
      type(vector_field),                       intent(in)    :: old_gradient_pressure
      type(vector_field_pointer), dimension(:), intent(inout) :: darcy_velocity
      type(vector_field),                       intent(inout) :: total_darcy_velocity
      type(scalar_field),                       intent(inout) :: div_total_darcy_velocity
      type(scalar_field),                       intent(in)    :: inverse_cv_mass_cfl_mesh
      type(scalar_field),                       intent(in)    :: inverse_cv_mass_pressure_mesh
      type(scalar_field),                       intent(inout) :: cv_mass_pressure_mesh_with_porosity
      type(scalar_field),                       intent(inout) :: cv_mass_pressure_mesh_with_old_porosity
      type(scalar_field_pointer), dimension(:), intent(inout) :: cfl
      type(scalar_field_pointer), dimension(:), intent(in)    :: old_cfl
      type(scalar_field),                       intent(inout) :: cfl_subcycle
      real,                                     intent(in)    :: saturation_max_courant_per_subcycle
      type(cv_options_type),                    intent(in)    :: saturation_cv_options
      type(cv_faces_type),                      intent(in)    :: cvfaces
      type(element_type),                       intent(in)    :: x_cvshape_full
      type(element_type),                       intent(in)    :: p_cvshape_full
      type(element_type),                       intent(in)    :: x_cvshape
      type(element_type),                       intent(in)    :: p_cvshape
      type(element_type),                       intent(in)    :: x_cvbdyshape
      real,                                     intent(in)    :: dt
      integer,                                  intent(in)    :: number_phase
      logical,                                  intent(in)    :: phase_one_saturation_diagnostic
      type(state_type),                         intent(inout) :: state_one ! this is only here for the petsc interface
      
      ! local variables
      integer :: p, number_subcycle
      real    :: dt_subcycle, max_cfl

      ! Copy to Old the CV mass on the pressure mesh with porosity 
      call set(cv_mass_pressure_mesh_with_old_porosity, cv_mass_pressure_mesh_with_porosity)
      
      ! Calculate the latest CV mass on the pressure mesh with porosity
      call compute_cv_mass(positions, cv_mass_pressure_mesh_with_porosity, porosity)      
            
      ! solve the pressure 
      call solve_pressure(pressure, &
                          pressure_matrix, &
                          rhs, &
                          old_inter_velocity_porosity, &
                          positions, &
                          positions_pressure_mesh, &
                          absolute_permeability, &
                          old_saturation, &
                          relative_permeability, &
                          viscosity, &
                          cv_mass_pressure_mesh_with_porosity, &
                          cv_mass_pressure_mesh_with_old_porosity, &
                          old_sum_saturation, &
                          saturation_cv_options, &
                          cvfaces, &
                          x_cvshape_full, &
                          p_cvshape_full, &
                          x_cvshape, &
                          p_cvshape, &
                          number_phase, &
                          state_one, &
                          dt)
      
      ! Calculate the gradient pressure, inter_velocity_porosity ... etc
      call darcy_impes_calculate_gradient_pressure_etc(pressure, & 
                                                       positions,&
                                                       porosity, &
                                                       absolute_permeability,&
                                                       relative_permeability,&
                                                       viscosity,&
                                                       inter_velocity_porosity, &
                                                       gradient_pressure, &
                                                       darcy_velocity, &
                                                       total_darcy_velocity, &
                                                       div_total_darcy_velocity, &
                                                       old_saturation, &
                                                       inverse_cv_mass_cfl_mesh, &
                                                       inverse_cv_mass_pressure_mesh, &
                                                       cfl, &
                                                       saturation_cv_options, &
                                                       cvfaces, &
                                                       x_cvshape, &
                                                       p_cvshape, &
                                                       number_phase, &
                                                       dt)
            
      ! solve the saturation for each phase but not the first
      s_phase_loop: do p = 2, number_phase
         
         ! Deduce the number of subcycles to do and the subcycle time step size
         
         ! Find the max cfl number for this phase (accounting for parallel)
         max_cfl = maxval(cfl(p)%ptr)
         
         call allmax(max_cfl)
         
         ! Find the subcycle dt and cfl field (which may be used for face value limiting)
         number_subcycle = max(1,ceiling(max_cfl/saturation_max_courant_per_subcycle))
           
         dt_subcycle = dt/real(number_subcycle)

         call set(cfl_subcycle, cfl(p)%ptr)
         
         call scale(cfl_subcycle, 1.0/real(number_subcycle))
         
         call solve_phase_saturation(saturation(p)%ptr, &
                                     pressure, &
                                     inter_velocity_porosity(p)%ptr, &
                                     old_inter_velocity_porosity(p)%ptr, &
                                     lhs, &
                                     rhs, &
                                     rhs_adv, &
                                     rhs_time, &
                                     old_saturation_subcycle, &
                                     positions, &
                                     positions_pressure_mesh, &
                                     cv_mass_pressure_mesh_with_porosity, &
                                     cv_mass_pressure_mesh_with_old_porosity, &
                                     old_saturation(p)%ptr, &
                                     saturation_cv_options, &
                                     cvfaces, &
                                     x_cvshape, &
                                     p_cvshape, &
                                     x_cvbdyshape, &
                                     cfl_subcycle, &
                                     dt_subcycle, &
                                     number_subcycle)
                  
      end do s_phase_loop
      
      ! Either solve the first phase saturation or calculate if diagnostic
      s1_solve_if: if (phase_one_saturation_diagnostic) then
      
         call darcy_impes_calculate_phase_one_saturation_diagnostic(saturation, &
                                                                    number_phase)
      
      else s1_solve_if

         ! Deduce the number of subcycles to do and the subcycle time step size
         
         ! Find the max cfl number for this phase (accounting for parallel)
         max_cfl = maxval(cfl(1)%ptr)
         
         call allmax(max_cfl)
         
         ! Find the subcycle dt and cfl field (which may be used for face value limiting)
         number_subcycle = max(1,ceiling(max_cfl/saturation_max_courant_per_subcycle))
           
         dt_subcycle = dt/real(number_subcycle)
         
         call set(cfl_subcycle, cfl(1)%ptr)
         
         call scale(cfl_subcycle, 1.0/real(number_subcycle))
         
         call solve_phase_saturation(saturation(1)%ptr, &
                                     pressure, &
                                     inter_velocity_porosity(1)%ptr, &
                                     old_inter_velocity_porosity(1)%ptr, &
                                     lhs, &
                                     rhs, &
                                     rhs_adv, &
                                     rhs_time, &
                                     old_saturation_subcycle, &
                                     positions, &
                                     positions_pressure_mesh, &
                                     cv_mass_pressure_mesh_with_porosity, &
                                     cv_mass_pressure_mesh_with_old_porosity, &
                                     old_saturation(1)%ptr, &
                                     saturation_cv_options, &
                                     cvfaces, &
                                     x_cvshape, &
                                     p_cvshape, &
                                     x_cvbdyshape, &
                                     cfl_subcycle, &
                                     dt_subcycle, &
                                     number_subcycle)
      
      end if s1_solve_if

      ! Calculate the sum of the saturations
      call darcy_impes_calculate_sum_saturation(sum_saturation, &
                                                saturation, &
                                                number_phase)
            
   end subroutine darcy_impes_assemble_and_solve

! ----------------------------------------------------------------------------

   subroutine solve_pressure(pressure, &
                             pressure_matrix, &
                             rhs, &
                             old_inter_velocity_porosity, &
                             positions, &
                             positions_pressure_mesh, &
                             absolute_permeability, &
                             old_saturation, &
                             relative_permeability, &
                             viscosity, &
                             cv_mass_pressure_mesh_with_porosity, &
                             cv_mass_pressure_mesh_with_old_porosity, &
                             old_sum_saturation, &
                             saturation_cv_options, &
                             cvfaces, &
                             x_cvshape_full, &
                             p_cvshape_full, &
                             x_cvshape, &
                             p_cvshape, &
                             number_phase, &
                             state_one, &
                             dt)
      
      !!< Assemble and solve the pressure
      
      type(scalar_field),                       intent(inout) :: pressure
      type(csr_matrix),                         intent(inout) :: pressure_matrix
      type(scalar_field),                       intent(inout) :: rhs
      type(vector_field_pointer), dimension(:), intent(in)    :: old_inter_velocity_porosity
      type(vector_field),                       intent(in)    :: positions
      type(vector_field),                       intent(in)    :: positions_pressure_mesh
      type(scalar_field),                       intent(in)    :: absolute_permeability
      type(scalar_field_pointer), dimension(:), intent(in)    :: old_saturation
      type(scalar_field_pointer), dimension(:), intent(in)    :: relative_permeability
      type(scalar_field_pointer), dimension(:), intent(in)    :: viscosity
      type(scalar_field),                       intent(in)    :: cv_mass_pressure_mesh_with_porosity
      type(scalar_field),                       intent(in)    :: cv_mass_pressure_mesh_with_old_porosity
      type(scalar_field),                       intent(in)    :: old_sum_saturation
      type(cv_options_type),                    intent(in)    :: saturation_cv_options
      type(cv_faces_type),                      intent(in)    :: cvfaces
      type(element_type),                       intent(in)    :: x_cvshape_full
      type(element_type),                       intent(in)    :: p_cvshape_full
      type(element_type),                       intent(in)    :: x_cvshape
      type(element_type),                       intent(in)    :: p_cvshape
      integer,                                  intent(in)    :: number_phase
      type(state_type),                         intent(inout) :: state_one ! this is only here for the petsc interface
      real,                                     intent(in)    :: dt
      
      ! local variables
      integer :: p, vele, iloc, oloc, jloc, face, gi, ggi
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
      
      ! allocate arrays used in assemble process - many assume
      ! that all elements are the same type.
      allocate(x_ele(positions%dim, ele_loc(positions,1)))      
      allocate(p_dshape(ele_loc(pressure,1), x_cvshape%ngi, mesh_dim(pressure)))
      allocate(normal(positions%dim,x_cvshape%ngi))
      allocate(detwei(x_cvshape%ngi))      
      allocate(normgi(positions%dim))
      allocate(notvisited(x_cvshape%ngi))
      allocate(p_mat_local(ele_loc(pressure,1), ele_loc(pressure,1)))
      allocate(x_face_quad(positions%dim, x_cvshape%ngi))
      allocate(minus_grad_old_p_quad(positions%dim, p_cvshape%ngi))
      allocate(visc_vele(number_phase))
      allocate(old_saturation_ele(number_phase,ele_loc(pressure,1)))
      allocate(relperm_ele(number_phase,ele_loc(pressure,1)))
      allocate(old_vphi_face(positions%dim,p_cvshape%ngi,number_phase))
      
      ! Initialise the matrix and rhs
      call zero(pressure_matrix)
      call zero(rhs)
                 
      ! Loop volume elements assembling local contributions     
      vol_element_loop: do vele = 1,element_count(pressure)
         
         ! get the saturation ele values for each phase
         do p = 1,number_phase
            old_saturation_ele(p,:) = ele_val(old_saturation(p)%ptr, vele)
         end do
         
         ! get the relperm ele values for each phase
         do p = 1,number_phase
            relperm_ele(p,:) = ele_val(relative_permeability(p)%ptr, vele)
         end do
         
         ! get the viscosity value for this element for each phase
         do p = 1,number_phase
            dummy(1:1) = ele_val(viscosity(p)%ptr, vele)
            visc_vele(p) = dummy(1)
         end do
         
         ! get the absolute permeability value for this element
         dummy(1:1) = ele_val(absolute_permeability, vele)         
         absperm_vele = dummy(1)
         
         ! get the coordinate values for this element for each positions local node
         x_ele = ele_val(positions, vele)         
         
         ! get the coordinate values for this element for each quadrature point
         x_face_quad = ele_val_at_quad(positions, vele, x_cvshape)
         
         ! The node indices of the pressure field
         p_nodes => ele_nodes(pressure, vele)
         
         ! The node indices of the positions projected to the pressure mesh
         x_pmesh_nodes => ele_nodes(positions_pressure_mesh, vele)
         
         ! Determine the node numbers to use to determine the
         ! Saturation and RelPerm upwind values
         if((saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
            (saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then
            
            upwind_nodes => x_pmesh_nodes
        
         else
            
            upwind_nodes => p_nodes
         
         end if
         
         ! obtain the transformed determinant*weight and normals
         call transform_cvsurf_to_physical(x_ele, x_cvshape, &
                                           detwei, normal, cvfaces)
         
         ! obtain the derivative of the pressure shape function at 
         ! the CV face quadrature points
         call transform_to_physical(positions, vele, x_shape = x_cvshape_full, &
                                    shape = p_cvshape_full, dshape = p_dshape)
         
         ! the old intersitial velocity_porosity at the quadrature points, used to determine upwind
         ! via a FE interpolation using the pressure basis functions (which for linear is fine)
         do p = 1,number_phase            
            old_vphi_face(:,:,p) = ele_val_at_quad(old_inter_velocity_porosity(p)%ptr, vele, p_cvshape)         
         end do 
         
         ! Initialise array for the quadrature points of this 
         ! element for whether it has already been visited
         notvisited = .true.
         
         ! Initialise the local p matrix to assemble for this element
         p_mat_local = 0.0

         ! loop over local nodes within this element
         nodal_loop_i: do iloc = 1, pressure%mesh%shape%loc

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

                   ! correct the orientation of the normal so it points away from iloc
                   normgi = orientate_cvsurf_normgi(node_val(positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                   &x_face_quad(:,ggi), normal(:,ggi))
                                      
                   ! Evaluate the face value for old saturation and old relperm for each phase
                   ! (assuming upwind for now) and evaluate the face value.

                   ! face value = sum_{phase} ( S*relperm*absperm/visc ), where absperm is phase independent
                   ! (if absperm and viscosity are considered tensors this requires modifying below)
                   
                   face_value = 0.0
                      
                   do p = 1,number_phase

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
                   do jloc = 1,pressure%mesh%shape%loc
                   
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
         call addto(pressure_matrix, p_nodes, p_nodes, p_mat_local)
         
      end do vol_element_loop
      
      ! Add rate of change of porosity to rhs    
      call set(rhs, cv_mass_pressure_mesh_with_old_porosity)
                  
      call addto(rhs, cv_mass_pressure_mesh_with_porosity, scale = -1.0)
      
      call scale(rhs, 1.0/dt)
      
      ! Apply any strong dirichlet BC's
      call apply_dirichlet_conditions(pressure_matrix, rhs, pressure)
      
      ! Solve the pressure
      call petsc_solve(pressure, pressure_matrix, rhs, state_one)
      
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
   
   end subroutine solve_pressure

! ----------------------------------------------------------------------------

   subroutine darcy_impes_calculate_gradient_pressure_etc(pressure, & 
                                                          positions,&
                                                          porosity, &
                                                          absolute_permeability,&
                                                          relative_permeability,&
                                                          viscosity,&
                                                          inter_velocity_porosity, &
                                                          gradient_pressure, &
                                                          darcy_velocity, &
                                                          total_darcy_velocity, &
                                                          div_total_darcy_velocity, &
                                                          old_saturation, &
                                                          inverse_cv_mass_cfl_mesh, &
                                                          inverse_cv_mass_pressure_mesh, &
                                                          cfl, &
                                                          saturation_cv_options, &
                                                          cvfaces, &
                                                          x_cvshape, &
                                                          p_cvshape, &
                                                          number_phase, &
                                                          dt)
      
      !!< Calculate the gradient pressure, inter_velocity_porosity, CFL and darcy velocity fields

      type(scalar_field),                       intent(in)    :: pressure
      type(vector_field),                       intent(in)    :: positions
      type(scalar_field),                       intent(in)    :: porosity      
      type(scalar_field),                       intent(in)    :: absolute_permeability
      type(scalar_field_pointer), dimension(:), intent(in)    :: relative_permeability
      type(scalar_field_pointer), dimension(:), intent(in)    :: viscosity
      type(vector_field_pointer), dimension(:), intent(inout) :: inter_velocity_porosity
      type(vector_field),                       intent(inout) :: gradient_pressure
      type(vector_field_pointer), dimension(:), intent(inout) :: darcy_velocity
      type(vector_field),                       intent(inout) :: total_darcy_velocity
      type(scalar_field),                       intent(inout) :: div_total_darcy_velocity
      type(scalar_field_pointer), dimension(:), intent(in)    :: old_saturation      
      type(scalar_field),                       intent(in)    :: inverse_cv_mass_cfl_mesh
      type(scalar_field),                       intent(in)    :: inverse_cv_mass_pressure_mesh
      type(scalar_field_pointer), dimension(:), intent(inout) :: cfl      
      type(cv_options_type),                    intent(in)    :: saturation_cv_options
      type(cv_faces_type),                      intent(in)    :: cvfaces
      type(element_type),                       intent(in)    :: x_cvshape
      type(element_type),                       intent(in)    :: p_cvshape
      integer,                                  intent(in)    :: number_phase
      real,                                     intent(in)    :: dt
      
      ! local variables
      integer :: vele, p, nloc, dim
      real :: visc_ele, absperm_ele
      real, dimension(1) :: dummy
      real, dimension(positions%dim) :: grad_pressure_ele
      real, dimension(:), allocatable :: relperm_ele
      real, dimension(:), allocatable :: inter_velocity_porosity_local
      
      ! Calculate the gradient pressure field assumed DG
      call grad(pressure, &
                positions, &
                gradient_pressure)
      
      ! Calculate the inter_velocity_porosity field = - (1.0/old_sigma) * grad_pressure, where sigma = visc / absperm * relperm
      
      allocate(relperm_ele(ele_loc(pressure,1)))
      
      ! deduce the number of local nodes (all elements assumed same) for inter_velocity_porosity and CFL
      nloc = inter_velocity_porosity(1)%ptr%mesh%shape%loc

      allocate(inter_velocity_porosity_local(nloc))      
      
      vol_element_loop: do vele = 1, element_count(pressure)
         
         ! Find the element wise abs perm, porosity and grad_pressure
         ! (a dummy is used such as to not assume that the node number
         !  is the same as the element number)
         
         dummy = ele_val(absolute_permeability, vele)
         absperm_ele = dummy(1)
                  
         do dim = 1,positions%dim         
            dummy = ele_val(gradient_pressure, dim, vele)
            grad_pressure_ele(dim) = dummy(1)         
         end do
         
         ! Form the inter_velocity_porosity and CFL for each phase
         phase_loop: do p = 1,number_phase
         
            ! Find the element wise viscosity
            dummy = ele_val(viscosity(p)%ptr, vele)
            visc_ele = dummy(1)
            
            ! Find the element wise local values for relperm
            relperm_ele = ele_val(relative_permeability(p)%ptr, vele)
            
            ! find the local DG values for inter_velocity_porosity for each geometric dimension
            ! (if absperm and viscosity are considered tensors this requires modifying below)
            dim_loop: do dim = 1,positions%dim
               
               inter_velocity_porosity_local(:) = - relperm_ele(:) * absperm_ele * grad_pressure_ele(dim) / visc_ele
            
               call set(inter_velocity_porosity(p)%ptr, &
                        dim, &
                        ele_nodes(inter_velocity_porosity(p)%ptr, vele), &
                        inter_velocity_porosity_local)
                           
            end do dim_loop
            
         end do phase_loop
         
      end do vol_element_loop
      
      ! calculate the darcy velocity = inter_velocity_porosity * old_saturation
      do p = 1, number_phase         
         
         call zero(darcy_velocity(p)%ptr)
         
         do vele = 1, element_count(pressure)
            
            do dim = 1,positions%dim
            
               call addto(darcy_velocity(p)%ptr, &
                          dim, &
                          ele_nodes(darcy_velocity(p)%ptr, vele), &
                          ele_val(old_saturation(p)%ptr, vele))
            
            end do
            
         end do 

         call scale(darcy_velocity(p)%ptr, inter_velocity_porosity(p)%ptr)
         
      end do 
      
      ! calculate the total darcy velocity
      call set(total_darcy_velocity, darcy_velocity(1)%ptr)
      
      do p = 2, number_phase
         
         call addto(total_darcy_velocity, darcy_velocity(p)%ptr)
         
      end do
      
      ! calculate the divergence of the total darcy velocity
      call zero(div_total_darcy_velocity)
      
      do vele = 1, element_count(pressure)
         
         
         
      end do
      
      call scale(div_total_darcy_velocity, inverse_cv_mass_pressure_mesh)
      
      ! calculate the CFL number using CV surface integrals rather than volume ones
      do p = 1, number_phase
         call zero(cfl(p)%ptr)
      end do
      
      deallocate(relperm_ele)
      deallocate(inter_velocity_porosity_local)
      
   end subroutine darcy_impes_calculate_gradient_pressure_etc

! ----------------------------------------------------------------------------

   subroutine solve_phase_saturation(saturation, &
                                     pressure, &
                                     inter_velocity_porosity, &
                                     old_inter_velocity_porosity, &
                                     lhs, &
                                     rhs, &
                                     rhs_adv, &
                                     rhs_time, &
                                     old_saturation_subcycle, &
                                     positions, &
                                     positions_pressure_mesh, &
                                     cv_mass_pressure_mesh_with_porosity, &
                                     cv_mass_pressure_mesh_with_old_porosity, &
                                     old_saturation, &
                                     saturation_cv_options, &
                                     cvfaces, &
                                     x_cvshape, &
                                     p_cvshape, & 
                                     x_cvbdyshape, &                                    
                                     cfl_subcycle, &
                                     dt_subcycle, &
                                     number_subcycle)
      
      !!< Assemble and solve the saturation

      type(scalar_field),    intent(inout) :: saturation      
      type(scalar_field),    intent(in)    :: pressure
      type(vector_field),    intent(in)    :: inter_velocity_porosity
      type(vector_field),    intent(in)    :: old_inter_velocity_porosity
      type(scalar_field),    intent(inout) :: lhs   
      type(scalar_field),    intent(inout) :: rhs
      type(scalar_field),    intent(inout) :: rhs_adv
      type(scalar_field),    intent(inout) :: rhs_time
      type(scalar_field),    intent(inout) :: old_saturation_subcycle
      type(vector_field),    intent(in)    :: positions
      type(vector_field),    intent(in)    :: positions_pressure_mesh
      type(scalar_field),    intent(in)    :: cv_mass_pressure_mesh_with_porosity
      type(scalar_field),    intent(in)    :: cv_mass_pressure_mesh_with_old_porosity
      type(scalar_field),    intent(in)    :: old_saturation
      type(cv_options_type), intent(in)    :: saturation_cv_options
      type(cv_faces_type),   intent(in)    :: cvfaces
      type(element_type),    intent(in)    :: x_cvshape
      type(element_type),    intent(in)    :: p_cvshape
      type(element_type),    intent(in)    :: x_cvbdyshape
      type(scalar_field),    intent(in)    :: cfl_subcycle
      real,                  intent(in)    :: dt_subcycle
      integer,               intent(in)    :: number_subcycle

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
      real,    dimension(:),     allocatable :: ghost_old_saturation_ele_bdy
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
      allocate(x_ele(positions%dim, ele_loc(positions,1)))      
      allocate(normal(positions%dim,x_cvshape%ngi))
      allocate(detwei(x_cvshape%ngi))      
      allocate(normgi(positions%dim))
      allocate(notvisited(x_cvshape%ngi))
      allocate(s_rhs_local(ele_loc(pressure,1)))
      allocate(x_face_quad(positions%dim, x_cvshape%ngi))
      allocate(old_saturation_ele(ele_loc(pressure,1)))
      allocate(vphi_ele(positions%dim,ele_loc(inter_velocity_porosity,1)))
      allocate(old_vphi_face(positions%dim,p_cvshape%ngi))
      allocate(vphi_face_value(positions%dim))
      
      allocate(ghost_old_saturation_ele_bdy(face_loc(pressure,1)))
      allocate(old_saturation_ele_bdy(face_loc(pressure,1)))
      allocate(detwei_bdy(x_cvbdyshape%ngi))
      allocate(normal_bdy(positions%dim, x_cvbdyshape%ngi))
      allocate(s_rhs_local_bdy(face_loc(pressure,1)))
      allocate(x_ele_bdy(positions%dim, face_loc(positions,1)))      
      allocate(vphi_sele(positions%dim,face_loc(inter_velocity_porosity,1)))
      allocate(p_nodes_bdy(face_loc(pressure%mesh,1)))
      
      ! Inititalise rhs advection field
      call zero(rhs_adv)
                 
      ! Loop volume elements assembling local contributions    
      vol_element_loop: do vele = 1,element_count(pressure)
         
         ! get the saturation ele values
         old_saturation_ele(:) = ele_val(old_saturation, vele)

         ! The interstitial velocity porosity ele value
         vphi_ele(:,:) = ele_val(inter_velocity_porosity, vele)

         ! the old intersitial velocity_porosity at the quadrature points, used to determine upwind
         ! via a FE interpolation using the pressure basis functions (which for linear is fine)
         old_vphi_face(:,:) = ele_val_at_quad(old_inter_velocity_porosity, vele, p_cvshape)         
                  
         ! get the coordinate values for this element for each positions local node
         x_ele = ele_val(positions, vele)         
         
         ! get the coordinate values for this element for each quadrature point
         x_face_quad = ele_val_at_quad(positions, vele, x_cvshape)
         
         ! The node indices of the pressure field
         p_nodes => ele_nodes(pressure, vele)
         
         ! The node indices of the positions projected to the pressure mesh
         x_pmesh_nodes => ele_nodes(positions_pressure_mesh, vele)
         
         ! Determine the node numbers to use to determine the
         ! Saturation and RelPerm upwind values
         if((saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_POINT).or.&
            (saturation_cv_options%upwind_scheme == CV_UPWINDVALUE_PROJECT_GRAD)) then
            
            upwind_nodes => x_pmesh_nodes
        
         else
            
            upwind_nodes => p_nodes
         
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
         nodal_loop_i: do iloc = 1, pressure%mesh%shape%loc

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

                   ! correct the orientation of the normal so it points away from iloc
                   normgi = orientate_cvsurf_normgi(node_val(positions_pressure_mesh, x_pmesh_nodes(iloc)), &
                                                   &x_face_quad(:,ggi), normal(:,ggi))
                   
                   ! determine if the flow is in or out of the face at this quadrature
                   ! with respect to the normal orientation using the old vphi - same as pressure assemble
                   old_vphi_dot_n = dot_product(old_vphi_face(:,ggi), normgi(:))
                                      
                   inflow = (old_vphi_dot_n<=0.0)
                   
                   income = merge(1.0,0.0,inflow)

                   ! Evaluate the face value for old saturation and interstitial velocity
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
         call addto(rhs_adv, p_nodes, s_rhs_local)
         
      end do vol_element_loop
      
      ! Add BC surface integrals

      allocate(saturation_bc_type(surface_element_count(saturation)))         
      
      call get_entire_boundary_condition(saturation, (/"zero_flux    ", &
                                                       "dirichlet    "/), saturation_bc, saturation_bc_type)
     
      sele_loop: do sele = 1, surface_element_count(saturation)
         
         if (saturation_bc_type(sele) == BC_TYPE_ZERO_FLUX) cycle
         
         if (saturation_bc_type(sele) == BC_TYPE_DIRICHLET) cycle
         
         old_saturation_ele_bdy = face_val(old_saturation, sele)
                           
         vphi_sele(:,:) = face_val(inter_velocity_porosity, sele)
         
         x_ele       = ele_val(positions, face_ele(positions, sele))
         x_ele_bdy   = face_val(positions, sele)
         p_nodes_bdy = face_global_nodes(pressure%mesh, sele)
         
         call transform_cvsurf_facet_to_physical(x_ele, x_ele_bdy, x_cvbdyshape, normal_bdy, detwei_bdy)

         ! Initialise the local rhs to assemble for this element
         s_rhs_local_bdy = 0.0

         bc_iloc_loop: do iloc = 1, pressure%mesh%faces%shape%loc

            bc_face_loop: do face = 1, cvfaces%sfaces

               bc_neigh_if: if(cvfaces%sneiloc(iloc,face)/=0) then

                  bc_quad_loop: do gi = 1, cvfaces%shape%ngi

                     ggi = (face-1)*cvfaces%shape%ngi + gi

                     vphi_face_value_dot_n = dot_product(vphi_sele(:,iloc), normal_bdy(:,ggi))
                   
                     face_value = old_saturation_ele_bdy(iloc) * vphi_face_value_dot_n * detwei_bdy(ggi)
                     
                     s_rhs_local_bdy(iloc) = s_rhs_local_bdy(iloc) - face_value                     
                     
                  end do bc_quad_loop

               end if bc_neigh_if

            end do bc_face_loop

         end do bc_iloc_loop

         call addto(rhs_adv, p_nodes_bdy, s_rhs_local_bdy)
                  
      end do sele_loop

      ! Solve the saturation for each subcycle via:            
      !  - Form the lhs = cv_mass_pressure_mesh_with_porosity / dt
      !  - and rhs = rhs_adv + rhs_time, where rhs_time = cv_mass_pressure_mesh_with_old_porosity * old_saturation / dt
      !  - solve for latest saturation and copy to start subcycle saturation step
      
      ! Note: the porosity at the start and end of a subcycle time step
      ! are linearly interpolated values from the main time step start and end
      
      call set(old_saturation_subcycle, old_saturation)

      sub_loop: do isub = 1,number_subcycle
         
         ! Form the start and end of subcycle dt
         ! porosity linear interpolents
         alpha_start = (isub - 1) / number_subcycle
         alpha_end   = isub / number_subcycle
         
         call set(lhs, cv_mass_pressure_mesh_with_porosity)
         
         call scale(lhs, alpha_end)
         
         call addto(lhs, cv_mass_pressure_mesh_with_old_porosity, scale = (1.0 - alpha_end))

         call scale(lhs, 1.0/dt_subcycle)

         call set(rhs_time, cv_mass_pressure_mesh_with_porosity)
         
         call scale(rhs_time, alpha_start)
         
         call addto(rhs_time, cv_mass_pressure_mesh_with_old_porosity, scale = (1.0 - alpha_start))
         
         call scale(rhs_time, 1.0/dt_subcycle)

         ! add rhs time term
         call set(rhs, rhs_time)         
         call scale(rhs, old_saturation_subcycle)

         ! add rhs advection term
         call addto(rhs, rhs_adv)

         ! apply strong dirichlet BC
         call apply_dirichlet_conditions(lhs, rhs, saturation)

         ! Solve for the saturation
         saturation%val = rhs%val / lhs%val
         
         ! Copy latest solution to old subcycle
         call set(old_saturation_subcycle, saturation)
         
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

      deallocate(ghost_old_saturation_ele_bdy)
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
   
   subroutine darcy_impes_calculate_phase_one_saturation_diagnostic(saturation, &
                                                                    number_phase)
      
      !!< Calculate the first phase diagnostic saturation as 1.0 - the rest
      
      type(scalar_field_pointer), dimension(:), intent(inout) :: saturation
      integer,                                  intent(in)    :: number_phase
      
      ! local variables
      integer :: p
      
      call set(saturation(1)%ptr, 1.0)
      
      do p = 2, number_phase
         
         call addto(saturation(1)%ptr, saturation(p)%ptr, scale = -1.0)
         
      end do
       
   end subroutine darcy_impes_calculate_phase_one_saturation_diagnostic

! ----------------------------------------------------------------------------
   
   subroutine darcy_impes_calculate_sum_saturation(sum_saturation, &
                                                   saturation, &                 
                                                   number_phase)
      
      !!< Calculate the sum of the saturation fields
      
      type(scalar_field),                       intent(inout) :: sum_saturation
      type(scalar_field_pointer), dimension(:), intent(in)    :: saturation      
      integer,                                  intent(in)    :: number_phase
      
      ! local variables
      integer :: p
      
      call zero(sum_saturation)
      
      do p = 1, number_phase
         
         call addto(sum_saturation, saturation(p)%ptr)
         
      end do
       
   end subroutine darcy_impes_calculate_sum_saturation

! ----------------------------------------------------------------------------

end module darcy_impes_assemble_module
