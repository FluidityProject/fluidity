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

module multiphase_diagnostics
   !!< Module containing all generic multiphase-related diagnostic algorithms.

   use field_options
   use fields
   use fldebug
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_fields_module
   use state_module

   implicit none
   
   private

   public :: calculate_particle_reynolds_number
  
   contains

      subroutine calculate_particle_reynolds_number(states, state_index, s_field)
         !!< Calculates the particle Reynolds number, (density_c * |u_c - u_p| * d)/viscosity_c
         !!< where _c and _p denote the continuous and particle phases respectively,
         !!< and d is the particle diameter

         type(state_type), dimension(:), intent(inout) :: states
         integer, intent(in) :: state_index ! Index of the particle phase in states
         type(scalar_field), intent(inout) :: s_field ! Particle Reynolds number field

         !! Sub-options of the diagnostic field
         real :: d ! Particle diameter
         character(len = OPTION_PATH_LEN) :: continuous_phase_name

         !! Other local variables
         integer :: i, continuous_state_index

         ! Velocities of the continuous and particle phases
         type(vector_field), pointer :: u_continuous, u_particle
         type(vector_field), pointer :: x
         type(tensor_field), pointer :: viscosity
         type(scalar_field), pointer :: density

         ! Counters over the elements and Gauss points
         integer :: ele, gi
         ! Transformed quadrature weights.
         real, dimension(ele_ngi(s_field, 1)) :: detwei
         ! Inverse of the local coordinate change matrix.
         real, dimension(mesh_dim(s_field), mesh_dim(s_field), ele_ngi(s_field, 1)) :: J

         ! Field values at each quadrature point.
         real, dimension(mesh_dim(s_field), ele_ngi(s_field, 1)) :: particle_re_gi
         real, dimension(mesh_dim(s_field), ele_ngi(s_field, 1)) :: u_continuous_gi, u_particle_gi
         real, dimension(mesh_dim(s_field), mesh_dim(s_field), ele_ngi(s_field, 1)) :: viscosity_gi
         real, dimension(ele_ngi(s_field, 1)) :: density_gi

         ! Current element global node numbers.
         integer, dimension(:), pointer :: particle_re_nodes
         ! Current particle_re element shape
         type(element_type), pointer :: particle_re_shape

         ! Local particle_re matrix on the current element.
         real, dimension(ele_loc(s_field, 1),ele_loc(s_field, 1)) :: particle_re_mat


         ewrite(1,*) 'Entering calculate_particle_reynolds_number'

         ! Get sub-options from under the diagnostic field in the options tree
         call get_option(trim(s_field%option_path)//'/diagnostic/algorithm::particle_reynolds_number/particle_diameter', d)
         call get_option(trim(s_field%option_path)//'/diagnostic/algorithm::particle_reynolds_number/continuous_phase_name', continuous_phase_name)

         ! Find the index of the continuous phase in states
         do i = 1, size(states)
            if(trim(states(i)%name) == continuous_phase_name) then
               continuous_state_index = i
            end if
         end do

         ! Zero particle Reynolds number field
         call zero(s_field)

         u_continuous => extract_vector_field(states(continuous_state_index), "Velocity")
         u_particle => extract_vector_field(states(state_index), "Velocity")
         x => extract_vector_field(states(state_index), "Coordinate")
         viscosity => extract_tensor_field(states(continuous_state_index),"Viscosity")
         density => extract_scalar_field(states(continuous_state_index),"Density")

         ! Loop through and integrate over each element
         do ele = 1, element_count(s_field)

            particle_re_nodes => ele_nodes(s_field, ele)
            particle_re_shape => ele_shape(s_field, ele)

            call compute_jacobian(ele_val(x,ele), ele_shape(x,ele), J=J, detwei=detwei)

            ! Calculate the particle_re number at each quadrature point.
            u_particle_gi = ele_val_at_quad(u_particle, ele)
            u_continuous_gi = ele_val_at_quad(u_continuous, ele)
            viscosity_gi = ele_val_at_quad(viscosity, ele)
            density_gi = ele_val_at_quad(density, ele)
            do gi = 1, size(detwei)
               particle_re_gi(:,gi) = matmul(u_continuous_gi(:,gi) - u_particle_gi(:,gi), J(:,:,gi))

               ! Divide through by the viscosity of the particle
               particle_re_gi(:,gi) = matmul(inverse(viscosity_gi(:,:,gi)), particle_re_gi(:,gi)) 
            end do

            ! Invert the mass matrix to get the particle_re value at each node
            particle_re_mat = matmul(inverse(shape_shape(particle_re_shape, particle_re_shape, detwei)), &
                  shape_shape(particle_re_shape, particle_re_shape, detwei*d*density_gi*maxval(abs(particle_re_gi),1)))

            ! (Taken from the GridReynoldsNumber field above)
            ! particle_re is inherently discontinuous. In the case where a continuous
            ! mesh is provided for particle_re, the following takes the safest option
            ! of taking the maximum value at a node.
            s_field%val(particle_re_nodes) = max(s_field%val(particle_re_nodes), sum(particle_re_mat,2))

         end do

         ewrite(1,*) 'Exiting calculate_particle_reynolds_number'

      end subroutine calculate_particle_reynolds_number

end module multiphase_diagnostics