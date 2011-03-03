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
!    but WITHOUT ANY WARRANTY; without even the implied arranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module radiation_solve_module
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   
   use radiation_materials
   use radiation_materials_interpolation
   use radiation_solve_power_iteration
   use radiation_normalise_flux

   implicit none
   
   private 

   public :: radiation_solve_eigenvalue, &
             radiation_solve_time
        
contains

   ! --------------------------------------------------------------------------

   subroutine radiation_solve_eigenvalue(state, &
                                         np_radmat_ii, &
                                         np_radmat) 
      
      !!< Solve the radiation eigenvalue problem for this np object
      
      type(state_type), intent(inout) :: state
      type(np_radmat_ii_type), intent(inout) :: np_radmat_ii
      type(np_radmat_type), intent(in) :: np_radmat

      ! local variables
      real :: keff
      
      ! form the material data interpolation/mixing instructions
      call form(np_radmat_ii, &
                np_radmat, &
                state)      
      
      ! solve the problem via a particular algorithm
      eig_solver: if (have_option(trim(np_radmat%option_path)//'/eigenvalue_run/power_iteration')) then
      
         call eigenvalue_power_iteration(state, &
                                         np_radmat, &
                                         np_radmat_ii, &
                                         keff)
             
      else eig_solver
      
         FLAbort('Unknown radiation eigenvalue solver')
      
      end if eig_solver
            
      ! normalise the neutral particle flux solution
      call normalise_neutral_particle_flux(state, &
                                           np_radmat)
                        
      ! calculate output diagnostics
            
      ! output the solution
            
   end subroutine radiation_solve_eigenvalue

   ! --------------------------------------------------------------------------

   subroutine radiation_solve_time(state, &
                                   np_radmat_ii, &
                                   np_radmat) 
      
      !!< Solve the radiation timestep problem for for this np object
      
      type(state_type), intent(inout) :: state
      type(np_radmat_ii_type), intent(inout) :: np_radmat_ii
      type(np_radmat_type), intent(in) :: np_radmat
      
            
      
      
   end subroutine radiation_solve_time

   ! --------------------------------------------------------------------------

end module radiation_solve_module
