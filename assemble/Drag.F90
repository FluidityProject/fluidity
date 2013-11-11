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
use fields_data_types
use fields
use state_module
use boundary_conditions
use fetools
use global_parameters, only : OPTION_PATH_LEN
use spud
use sparse_tools
use sparse_tools_petsc
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
   type(scalar_field), pointer:: drag_coefficient, distance_top, distance_bottom
   character(len=OPTION_PATH_LEN) bctype
   real, dimension(:), allocatable:: face_detwei, coefficient, density_face_gi
   real, dimension(:,:), allocatable:: drag_mat
   real dt, theta, gravity_magnitude
   integer, dimension(:), allocatable:: faceglobalnodes
   integer, dimension(:), pointer:: surface_element_list
   integer i, j, k, nobcs, stat
   integer snloc, sele, sngi
   logical:: parallel_dg, have_distance_bottom, have_distance_top, have_gravity, manning_strickler
   !wetting and drying, extra_drag represents a very large drag force enforced on the bottom surface
   integer::node
   logical :: have_wd, have_extra_drag  
   real :: d0, extra_drag
   type(scalar_field), pointer :: dtt,dtb
   real, dimension(:), allocatable:: depth_at_quads,extra_drag_mat
   type(scalar_field) :: depth
   have_wd=have_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying")
   have_extra_drag=have_option(trim(velocity%option_path)//&
    '/prognostic/boundary_conditions['//int2str(i-1)//']/type[0]/extra_drag')
   if (have_wd) then
     call get_option("/mesh_adaptivity/mesh_movement/free_surface/wetting_and_drying/d0", d0)
     if (have_extra_drag) then
       call get_option(trim(velocity%option_path)//&
        '/prognostic/boundary_conditions['//int2str(i-1)//']/type[0]/extra_drag',extra_drag)
     end if
   end if
   dtt => extract_scalar_field(state, "DistanceToTop")
   dtb => extract_scalar_field(state, "DistanceToBottom")
   call allocate(depth, dtt%mesh, "Depth")
   do node=1,node_count(dtt)
      call set(depth, node, node_val(dtt, node)+node_val(dtb, node))
   end do
   
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
   allocate(depth_at_quads(1:sngi),extra_drag_mat(1:sngi))
   nobcs=option_count(trim(velocity%option_path)//'/prognostic/boundary_conditions')
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
            else ! default to quadratic_drag
              ! drag coefficient: C_D * |u|
              coefficient=ele_val_at_quad(drag_coefficient, j)* &
                sqrt(sum(face_val_at_quad(nl_velocity, sele)**2, dim=1))
              if (manning_strickler) then
                 ! The manning-strickler formulation takes the form n**2g|u|u/(H**0.3333), where H is the water level, g is gravity and n is the Manning coefficient      
                 ! Note that distance_bottom+distance_top is the current water level H 
                 if(have_wd .and. have_extra_drag) then
                   depth_at_quads=face_val_at_quad(depth, sele)
                   do k=1, sngi
                     !if (depth_at_quads(k)<d0+eps) then
                       !coefficient(k)=1.0e+8   
                     extra_drag_mat(k)= max((2*d0-depth_at_quads(k))/d0,real(0))*extra_drag
                   end do            
                   coefficient =(ele_val_at_quad(drag_coefficient, j)+extra_drag_mat)*gravity_magnitude*coefficient/((face_val_at_quad(distance_bottom, sele)+face_val_at_quad(distance_top, sele))**(1./3.))  
                    print *, 'extra_drag_mat=',extra_drag_mat
                     !end if
                 else
                   coefficient=ele_val_at_quad(drag_coefficient, j)*gravity_magnitude*coefficient/((face_val_at_quad(distance_bottom, sele)+face_val_at_quad(distance_top, sele))**(1./3.))  
                    
                  end if
               end if
            end if
               
            ! density to turn this into a momentum absorption term
            ! (of course this will just be 1 with boussinesq)
            density_face_gi = face_val_at_quad(density, sele)
               
            drag_mat=shape_shape(face_shape(velocity, sele), &
               face_shape(velocity, sele), coefficient*face_detwei*density_face_gi)
               
            do k=1, velocity%dim
               call addto(bigm, k, k, faceglobalnodes, faceglobalnodes, &
                  dt*theta*drag_mat)
               call addto(rhs, k, faceglobalnodes, &
                  -matmul(drag_mat, face_val(old_velocity, k, sele)) )
            end do
            
         end do
            
      end if
   end do
     
   deallocate(faceglobalnodes, face_detwei, coefficient, drag_mat,density_face_gi)
   deallocate(depth_at_quads,extra_drag_mat)
end subroutine drag_surface

end module drag_module
