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

   module multiphase_module
      !! This module contains various subroutines and functions for 
      !! multiphase flow simulations
      use fldebug
      use state_module
      use fields
      use spud
      use global_parameters, only: OPTION_PATH_LEN
      use field_priority_lists
      use field_options
      use diagnostic_fields_matrices
      use sparse_tools
      use fields_manipulation
      use sparse_matrices_fields
      use state_matrices_module
      use field_derivatives
      use sparsity_patterns, only: make_sparsity
      use divergence_matrix_cg, only: assemble_divergence_matrix_cg
      use solvers
      implicit none

      private
      public :: get_phase_submaterials, get_nonlinear_volume_fraction, &
                calculate_diagnostic_phase_volume_fraction, &
                calculate_sum_velocity_divergence

   contains

      subroutine get_phase_submaterials(state, istate, submaterials, phase_istate)
         !!< Sets up an array of the submaterials of a phase.
         !!< NB: This includes the current state itself (i.e. state(istate)).

         type(state_type), dimension(:), target, intent(inout) :: state
         integer, intent(in) :: istate
         type(state_type), dimension(:), pointer :: submaterials
         integer, intent(inout), optional :: phase_istate

         !! Local variables
         integer :: i, next, stat, material_count, material_istate
         type(vector_field), pointer :: u
         character(len=OPTION_PATH_LEN) :: phase_name, target_name
         logical, dimension(:), pointer :: is_submaterial
         
         ewrite(1,*) 'Entering get_phase_submaterials'

         !! Store whether state(i) is a material or not in an array of logicals
         !! to save on computations in the second loop
         allocate(is_submaterial(size(state)))

         !! Get the number of submaterials so we can make submaterials the correct size

         material_count = 1 ! We will include state(istate) as one of the materials

         call get_option(trim(state(istate)%option_path)//"/name", phase_name)

         do i = 1, size(state)
        
            if(i == istate) then
               is_submaterial(i) = .true.
               cycle ! Already counted the current state
            end if
            
            u => extract_vector_field(state(i), "Velocity", stat)
            is_submaterial(i) = .false.

            ! If velocity field exists and is aliased to a phase's velocity field...
            if(stat == 0) then
               if(aliased(u)) then
                  ! ...then find out which phase it is aliased to. If it's the current phase,
                  ! then we have found one more submaterial.

                  ! Save the name of the phase that the current state is aliased to
                  call get_option(trim(state(i)%option_path)//"/vector_field::Velocity/aliased/material_phase_name", target_name)

                  if(target_name == phase_name) then
                     ! Found one more submaterial!
                     material_count = material_count + 1
                     is_submaterial(i) = .true.
                  end if
               end if
            end if

         end do

         ewrite(1,*) 'Number of sub-materials = ', material_count

         !! Allocate submaterials array
         allocate(submaterials(material_count))
         
         !! Assign the states to the submaterials array
         next = 1 ! Keep track of where we are in the submaterials array
         do i = 1, size(state)
            if(is_submaterial(i)) then
               submaterials(next) = state(i)

               ! Keep track of the phase's index in the new submaterials array
               if(present(phase_istate) .and. (i == istate)) then
                  phase_istate = next
               end if

               next = next + 1
            end if
         end do

         deallocate(is_submaterial)

         ewrite(1,*) 'Exiting get_phase_submaterials'

      end subroutine get_phase_submaterials


      subroutine get_nonlinear_volume_fraction(state, nvfrac)
         !!< Computes the nonlinear approximation to the phase volume fraction
         !!< and stores it in a locally allocated field before assembling the momentum equation.

         type(state_type), intent(inout) :: state
         type(scalar_field), intent(inout) :: nvfrac

         type(scalar_field), pointer :: volumefraction, oldvolumefraction
         type(vector_field), pointer :: velocity
         type(scalar_field) :: remapvfrac

         integer :: stat
         real :: theta

         ewrite(1,*) 'Entering get_nonlinear_volume_fraction'


         volumefraction => extract_scalar_field(state, 'PhaseVolumeFraction')
         
         call remap_field(volumefraction, nvfrac)
              
         velocity => extract_vector_field(state, 'Velocity', stat=stat)
         if(stat==0) then
            call get_option(trim(velocity%option_path)//'/prognostic/temporal_discretisation/theta', &
                              theta, stat)
            if(stat==0) then
               call allocate(remapvfrac, nvfrac%mesh, "RemppedPhaseVolumeFraction")
               
               oldvolumefraction => extract_scalar_field(state, 'OldPhaseVolumeFraction')

               ewrite_minmax(oldvolumefraction)
               ewrite_minmax(volumefraction)

               call remap_field(oldvolumefraction, remapvfrac)
               
               call scale(nvfrac, theta)
               call addto(nvfrac, remapvfrac, (1.0-theta))
               
               call deallocate(remapvfrac)
            end if
         end if
        
         ewrite_minmax(nvfrac%val)

         ewrite(1,*) 'Exiting get_nonlinear_volume_fraction'
      
      end subroutine get_nonlinear_volume_fraction


      subroutine calculate_diagnostic_phase_volume_fraction(state)
         !!< Searches for the state with the diagnostic PhaseVolumeFraction field,
         !!< and then computes it using the formula:
         !!< diagnostic volume fraction = 1.0 - (sum of all other volume fractions)

         type(state_type), dimension(:), intent(inout) :: state
         
         ! Local variables
         type(scalar_field), pointer :: phasevolumefraction
         integer :: i, stat, diagnostic_count
         type(scalar_field) :: sumvolumefractions
         type(scalar_field), pointer :: sfield
         logical :: diagnostic

         ewrite(1,*) 'Entering calculate_diagnostic_phase_volume_fraction'

         diagnostic_count = option_count("/material_phase/scalar_field::PhaseVolumeFraction/diagnostic")
         if(diagnostic_count>1) then
            ewrite(0,*) diagnostic_count, ' diagnostic PhaseVolumeFractions'
            FLExit("Only one diagnostic PhaseVolumeFraction permitted.")
         end if

         if(diagnostic_count==1) then
            ! Find the diagnostic volume fraction
            state_loop: do i = 1, size(state)
               phasevolumefraction=>extract_scalar_field(state(i), 'PhaseVolumeFraction', stat)
               if (stat==0) then
                  diagnostic = (have_option(trim(phasevolumefraction%option_path)//'/diagnostic'))
                  if((.not. aliased(phasevolumefraction)).and. diagnostic) then
                     exit state_loop
                  end if
               end if
            end do state_loop
            
            call allocate(sumvolumefractions, phasevolumefraction%mesh, 'Sum of volume fractions')
            call zero(sumvolumefractions)
            
            do i = 1, size(state)
               sfield=>extract_scalar_field(state(i),'PhaseVolumeFraction',stat)
               diagnostic=(have_option(trim(sfield%option_path)//'/diagnostic'))
               if ( (stat==0).and.(.not. aliased(sfield)).and.(.not.diagnostic)) then
                  call addto(sumvolumefractions, sfield)
               end if
            end do
            
            call set(phasevolumefraction, 1.0)
            call addto(phasevolumefraction, sumvolumefractions, -1.0)
            call deallocate(sumvolumefractions)
         end if

         ewrite(1,*) 'Exiting calculate_diagnostic_phase_volume_fraction'

      end subroutine calculate_diagnostic_phase_volume_fraction

      subroutine calculate_sum_velocity_divergence(state, sum_velocity_divergence)
         !!< Calculates \sum{div(vfrac*u)}, where we sum over each prognostic velocity 
         !!< field (i.e. each phase)

         type(state_type), dimension(:), intent(inout) :: state
         type(scalar_field), pointer :: sum_velocity_divergence

         ! Local variables
         type(vector_field), pointer :: u, x
         type(scalar_field), pointer :: vfrac
         integer :: i, stat

         type(csr_sparsity) :: divergence_sparsity
         type(block_csr_matrix) :: ct_m

         type(csr_sparsity) :: mass_sparsity
         type(csr_matrix) :: mass
         type(scalar_field) :: ctfield, ct_rhs, temp

         ewrite(1,*) 'Entering calculate_sum_velocity_divergence'

         ! Allocate memory for matrices and sparsity patterns
         call allocate(ctfield, sum_velocity_divergence%mesh, name="CTField")
         call zero (ctfield)
         call allocate(temp, sum_velocity_divergence%mesh, name="Temp")

         mass_sparsity=make_sparsity(sum_velocity_divergence%mesh, sum_velocity_divergence%mesh, "MassSparsity")
         call allocate(mass, mass_sparsity, name="MassMatrix")
         call zero(mass)

         ! Sum up over the div's
         do i = 1, size(state)
            u => extract_vector_field(state(i), "Velocity", stat)

            ! If there's no velocity then cycle
            if(stat/=0) cycle
            ! If this is an aliased velocity then cycle
            if(aliased(u)) cycle
            ! If the velocity isn't prognostic then cycle
            if(.not.have_option(trim(u%option_path)//"/prognostic")) cycle

            ! If velocity field is prognostic, begin calculations below
            x => extract_vector_field(state(i), "Coordinate")

            ! Allocate sparsity patterns, C^T matrix and C^T RHS for current state
            divergence_sparsity=make_sparsity(sum_velocity_divergence%mesh, u%mesh, "DivergenceSparsity")
            call allocate(ct_m, divergence_sparsity, (/1, u%dim/), name="DivergenceMatrix" )
            call allocate(ct_rhs, sum_velocity_divergence%mesh, name="CTRHS")

            ! Reassemble C^T matrix here
            if(i==1) then ! Construct the mass matrix (just do this once)          
               call assemble_divergence_matrix_cg(ct_m, state(i), ct_rhs=ct_rhs, &
                                 test_mesh=sum_velocity_divergence%mesh, field=u, &
                                 option_path=sum_velocity_divergence%option_path, div_mass=mass)
            else
               call assemble_divergence_matrix_cg(ct_m, state(i), ct_rhs=ct_rhs, &
                                 test_mesh=sum_velocity_divergence%mesh, field=u, &
                                 option_path=sum_velocity_divergence%option_path)
            end if

            ! Construct the linear system of equations to solve for \sum{div(vfrac*u)}
            call zero(temp)
            call mult(temp, ct_m, u)
            call addto(temp, ct_rhs, -1.0)

            ! Now add it to the sum
            call addto(ctfield, temp)

            call deallocate(ct_m)
            call deallocate(ct_rhs)
            call deallocate(divergence_sparsity)

         end do

         ! Solve for sum_velocity_divergence ( = \sum{div(vfrac*u)} )
         call zero(sum_velocity_divergence)
         call petsc_solve(sum_velocity_divergence, mass, ctfield)

         ! Deallocate memory
         call deallocate(mass_sparsity)
         call deallocate(mass)
         call deallocate(ctfield)
         call deallocate(temp)

         ewrite(1,*) 'Exiting calculate_sum_velocity_divergence'
            
      end subroutine calculate_sum_velocity_divergence

   end module multiphase_module
