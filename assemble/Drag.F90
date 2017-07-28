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
module drag_module
use fldebug
use global_parameters, only : OPTION_PATH_LEN
use spud
use futils, only: int2str, vmean
use parallel_tools
use sparse_tools
use fetools
use parallel_fields
use fields
use sparse_tools_petsc
use state_module
use boundary_conditions
use k_epsilon, only: get_friction_velocity 

implicit none

private
public drag_surface

contains

subroutine drag_surface(bigm, rhs, state, density)
!!< Applies quadratic or linear drag at boundaries with bc of "drag" type
!!< This version applies the proper surface integral.
   type(petsc_csr_matrix), intent(inout):: bigm
   type(vector_field), intent(inout):: rhs
   type(state_type), intent(in):: state
   type(scalar_field), intent(in) :: density
   
   type(vector_field), pointer:: velocity, nl_velocity, position, old_velocity
   type(tensor_field), pointer:: bg_visc
   type(scalar_field), pointer:: drag_coefficient, distance_top, distance_bottom
   character(len=OPTION_PATH_LEN) bctype
   type(vector_field) :: bc_value
   integer, dimension(:,:), allocatable :: bc_type
   real, dimension(:), allocatable:: face_detwei, coefficient, density_face_gi
   real, dimension(:,:), allocatable:: drag_mat
   real dt, theta, gravity_magnitude
   integer, dimension(:), allocatable:: faceglobalnodes
   integer, dimension(:), pointer:: surface_element_list
   integer i, j, k, nobcs, stat
   integer snloc, sele, sngi, log_bc_count
   logical:: parallel_dg, have_distance_bottom, have_distance_top, have_gravity, manning_strickler

   real                            :: yPlus, y, C_mu, kappa, beta
   type(scalar_field), pointer     :: TKE
   real, dimension(:), allocatable :: friction_velocity

   ewrite(1,*) 'Inside drag_surface'
   
   velocity => extract_vector_field(state, "Velocity")
   ! velocity at the beginning of the time step
   old_velocity => extract_vector_field(state, "OldVelocity")
   ! velocity weighted between old and new with theta
   nl_velocity => extract_vector_field(state, "NonlinearVelocity")
   position => extract_vector_field(state, "Coordinate")
   distance_bottom => extract_scalar_field(state, "DistanceToBottom", stat)
   have_distance_bottom = stat == 0
   distance_top => extract_scalar_field(state, "DistanceToTop", stat)
   have_distance_top = stat == 0
   bg_visc => extract_tensor_field(state, "BackgroundViscosity", stat)
   
   call get_option("/timestepping/timestep", dt)
   call get_option(trim(velocity%option_path)//"/prognostic/temporal_discretisation/theta", &
                      theta)
   call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude, &
          stat=stat)
   have_gravity = stat == 0
   parallel_dg=continuity(velocity)<0 .and. IsParallel()
                      
   sngi=face_ngi(velocity, 1)
   snloc=face_loc(velocity,1)
   
   allocate(faceglobalnodes(1:snloc), &
     face_detwei(1:sngi), coefficient(1:sngi), &
     drag_mat(1:snloc,1:snloc), density_face_gi(1:sngi))

!   if(have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon")) then
!     TKE => extract_scalar_field(state, "TurbulentKineticEnergy")
!   end if
   allocate(friction_velocity(1:sngi))

!   log_bc_count = option_count(trim(velocity%option_path)//'/prognostic/boundary_conditions/type/linear_drag')
   
   nobcs=option_count(trim(velocity%option_path)//'/prognostic/boundary_conditions')

   allocate(bc_type(velocity%dim, surface_element_count(velocity)))
   call get_entire_boundary_condition(velocity, ['drag'], bc_value, bc_type)

!   if (log_bc_count>0.0) then
!      call drag_friction_bc(velocity, u_tau)
!   end if

   do i=1, nobcs
      call get_boundary_condition(velocity, i, type=bctype, &
         surface_element_list=surface_element_list)
      if (bctype=='drag') then
         manning_strickler=have_option(trim(velocity%option_path)//&
               '/prognostic/boundary_conditions['//int2str(i-1)//']/type[0]/quadratic_drag/manning-strickler')
         if (manning_strickler) then
            if (.not. have_distance_bottom .or. .not. have_distance_top .or. .not. have_gravity) then
               ewrite(-1,*) "Manning-strickler drag needs DistanceToTop and DistanceToBottom fields and gravity."
               FLExit("Turn on ocean_boundaries underneath geometry.")
            end if
         end if
         drag_coefficient => extract_scalar_surface_field(velocity, i, "DragCoefficient")
         call get_option(trim(velocity%option_path)//&
              '/prognostic/boundary_conditions['//int2str(i-1)//']/type[0]/log_law_friction_velocity/yPlus', yPlus, default = 11.06) !A! grab yPlus from diamond
         call get_option(trim(velocity%option_path)//&
              '/prognostic/boundary_conditions['//int2str(i-1)//']/type[0]/log_law_friction_velocity/y', y, default = 0.0)
         call get_option(trim(velocity%option_path)//&
              '/prognostic/boundary_conditions['//int2str(i-1)//']/type[0]/log_law_friction_velocity/kappa', kappa, default = 0.41)
         call get_option(trim(velocity%option_path)//&
              '/prognostic/boundary_conditions['//int2str(i-1)//']/type[0]/log_law_friction_velocity/C_mu', C_mu, default = 0.09)
         call get_option(trim(velocity%option_path)//&
              '/prognostic/boundary_conditions['//int2str(i-1)//']/type[0]/log_law_friction_velocity/beta', beta, default = 5.2)

         do j=1, size(surface_element_list)
           
            sele=surface_element_list(j)
            if (parallel_dg) then
              if (.not. element_owned(velocity, face_ele(velocity, sele))) cycle
            end if

            call transform_facet_to_physical(position, sele, face_detwei)
            
            faceglobalnodes=face_global_nodes(nl_velocity, sele)
            
            if(have_option(trim(velocity%option_path)//&
               '/prognostic/boundary_conditions['//int2str(i-1)//']/type[0]/linear_drag')) then
              ! drag coefficient: C_D
              coefficient=ele_val_at_quad(drag_coefficient, j)

            elseif(have_option(trim(velocity%option_path)//&
               '/prognostic/boundary_conditions['//int2str(i-1)//']/type[0]/log_law_friction_velocity')) then

              if(have_option("/material_phase[0]/subgridscale_parameterisations/k-epsilon")) then
                TKE => extract_scalar_field(state, "TurbulentKineticEnergy")

                 ! calc friction velocity: u_tau_1 = |u_wall|/yPlus, u_tau_2 = C_mu^0.25*sqrt(k), u_tau = max ( u_tau_1 , u_tau_2 )
                 friction_velocity = get_friction_velocity(face_val(nl_velocity, sele), &
                      vmean(pack(ele_val(bg_visc, sele), .true.)), C_mu, y, kappa, beta, face_val(TKE,sele))
! max( sqrt(sum(face_val_at_quad(nl_velocity, sele)**2, dim=1)) / yPlus , &
!                                          sqrt(face_val_at_quad(TKE, sele)) * 0.09**0.25 )
              else
                 ! calc friction velocity: u_tau = u_tau_1
                 friction_velocity = get_friction_velocity(face_val(nl_velocity, sele), &
                      vmean(pack(ele_val(bg_visc, sele), .true.)), C_mu, y, kappa, beta)
!                 friction_velocity = sqrt(sum(face_val_at_quad(nl_velocity, sele)**2, dim=1)) / yPlus
              end if

              ! calc wall shear stress: tau_wall = - (u_tau/yPlus)*|u_wall|
              !coefficient = (friction_velocity/yPlus)*sqrt(sum(face_val_at_quad(nl_velocity, sele)**2, dim=1))
!              coefficient = (friction_velocity/yPlus)
              coefficient = friction_velocity(1)**2*sqrt(face_loc(nl_velocity, sele)/max(1.0e-8,sum(face_val(nl_velocity, sele)**2)))
              !ewrite(1,*) 'AMIN: Are we here yet?', yPlus, coefficient, coefficient*sqrt(sum(face_val_at_quad(nl_velocity, sele)**2, dim=1))
            else ! default to quadratic_drag
              ! drag coefficient: C_D * |u|
              coefficient=ele_val_at_quad(drag_coefficient, j)* &
                sqrt(sum((face_val_at_quad(nl_velocity, sele)&
                -ele_val_at_quad(bc_value,sele))**2, dim=1))
              if (manning_strickler) then
                 ! The manning-strickler formulation takes the form n**2g|u|u/(H**0.3333), where H is the water level, g is gravity and n is the Manning coefficient
                 ! Note that distance_bottom+distance_top is the current water level H
                 coefficient=ele_val_at_quad(drag_coefficient, j)*gravity_magnitude*coefficient/((face_val_at_quad(distance_bottom, sele)+face_val_at_quad(distance_top, sele))**(1./3.))
               end if
            end if
               
            ! density to turn this into a momentum absorption term
            ! (of course this will just be 1 with boussinesq)
            density_face_gi = face_val_at_quad(density, sele)
               
            drag_mat=shape_shape(face_shape(velocity, sele), &
               face_shape(velocity, sele), coefficient*face_detwei*density_face_gi)
               
            do k=1, velocity%dim
               Call addto(bigm, k, k, faceglobalnodes, faceglobalnodes, &
                  dt*theta*drag_mat)
               call addto(rhs, k, faceglobalnodes, &
                  -matmul(drag_mat, face_val(old_velocity, k, sele)))
               if (bc_type(k,sele) == 1) then
                  call addto(rhs, k, faceglobalnodes, &
                       +matmul(drag_mat, ele_val(bc_value,k,sele)))
               end if
            end do
            
         end do
            
      end if
   end do

   call deallocate(bc_value)
   deallocate(bc_type)
     
   deallocate(faceglobalnodes, face_detwei, coefficient, drag_mat)
   
end subroutine drag_surface

subroutine drag_friction_bc(velocity,u_tau, bc_value)
  type(vector_field) :: velocity, bc_value
  type(scalar_field) :: u_tau
 
  type(vector_field) :: surface_velocity
  type(scalar_field) :: y, yplus, nu

  integer :: sele
  integer, dimension(:), pointer :: snodes

!  call allocate(surface_velocity, bc_value%mesh, dim=velocity%dim)
  call allocate(y, bc_value%mesh)
  call allocate(yplus, bc_value%mesh)
  call allocate(nu, bc_value%mesh)

  do sele=1, surface_element_count(velocity)

     snodes => ele_nodes(bc_value,sele)
     call set(surface_velocity, snodes, &
          face_val(velocity, sele) - ele_val(bc_value, sele))
     call set(y, snodes, 1.0e-4)

  end do

!  call find_friction_velocity(0.41, 5.5, 

end subroutine drag_friction_bc

end module drag_module
