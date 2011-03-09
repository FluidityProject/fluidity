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

module radiation
   
   !!< Radiation model initialisation, cleanup and solve procedures
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module   
   
   use radiation_materials
   use radiation_materials_interpolation
   use radiation_check_options_module
   use radiation_diagnostics
   use radiation_solve_module
   
   implicit none
   
   private 

   public :: radiation_initialise, &
             radiation_cleanup, &
             radiation_solve, &
             np_radmat_type, &
             np_radmat_ii_type
         
contains

   ! --------------------------------------------------------------------------

   subroutine radiation_initialise(state, &
                                   np_radmats, &
                                   np_radmats_ii) 
      
      !!< Initialise the radiation model. First check the options and then 
      !!< create and set the radiation material databases as well is interpolation 
      !!< instructions for each particle type.
      
      type(state_type), intent(in) :: state
      type(np_radmat_type), dimension(:), allocatable, intent(inout) :: np_radmats
      type (np_radmat_ii_type), dimension(:), allocatable, intent(inout) :: np_radmats_ii
      
      ! local variable
      integer :: p 
      character(len=OPTION_PATH_LEN) :: np_radmat_option_path   
      
      ewrite(1,*) 'Initialise radiation model'
                     
      allocate(np_radmats(option_count('/radiation/particle_type')))
      
      allocate(np_radmats_ii(option_count('/radiation/particle_type')))
      
      call radiation_check_options()

      particle_type_loop: do p = 1,option_count('/radiation/particle_type')
    
         ! - 1 needed as options count from 0
         np_radmat_option_path = '/radiation/particle_type['//int2str(p - 1)//']'
      
         call create(np_radmats(p), &
                     trim(np_radmat_option_path))
         
         call np_radmat_read(np_radmats(p))
         
         call create(np_radmats_ii(p), &
                     np_radmats(p), &
                     state)
      
         ! register the radiation diagnostics for stat file for this p
         call radiation_register_diagnostics(trim(np_radmat_option_path))
                  
      end do particle_type_loop
            
      ewrite(1,*) 'Finished initialise radiation model'
      
   end subroutine radiation_initialise

   ! --------------------------------------------------------------------------

   subroutine radiation_cleanup(np_radmats, &
                                np_radmats_ii) 
   
      !! Finalise the radiation model

      type(np_radmat_type), dimension(:), allocatable, intent(inout) :: np_radmats
      type (np_radmat_ii_type), dimension(:), allocatable, intent(inout) :: np_radmats_ii
      
      ! local variable
      integer :: p 
      
      ewrite(1,*) 'Cleanup radiation model'
      
      particle_type_loop: do p = 1,size(np_radmats)
      
         call destroy(np_radmats(p))
         
         call destroy(np_radmats_ii(p))
      
      end do particle_type_loop 
      
      if (allocated(np_radmats)) deallocate(np_radmats)

      if (allocated(np_radmats_ii)) deallocate(np_radmats_ii)
          
   end subroutine radiation_cleanup

   ! --------------------------------------------------------------------------

   subroutine radiation_solve(state, &
                              np_radmats, &
                              np_radmats_ii, &
                              invoke_eigenvalue_solve) 
                              
      !!< Invoke the relevant radiation solver for each particle type
      
      type(state_type), intent(inout) :: state
      type(np_radmat_type), dimension(:), allocatable, intent(in) :: np_radmats
      type (np_radmat_ii_type), dimension(:), allocatable, intent(inout) :: np_radmats_ii      
      logical, intent(in) :: invoke_eigenvalue_solve

      ! local variable
      integer :: p 
      
      particle_type_loop: do p = 1,size(np_radmats)
      
         which_solve: if (invoke_eigenvalue_solve) then
            
            solve_eig: if (have_option(trim(np_radmats(p)%option_path)//'/eigenvalue_run')) then
               
               ewrite(1,*) 'Solve radiation model eigenvalue for particle type ',trim(np_radmats(p)%name)
               
               call radiation_solve_eigenvalue(state, &
                                               np_radmats_ii(p), &
                                               np_radmats(p))
            
            end if solve_eig
            
         else which_solve

            solve_time: if (have_option(trim(np_radmats(p)%option_path)//'/time_run')) then
               
               ewrite(1,*) 'Solve radiation model time for particle type ',trim(np_radmats(p)%name)
               
               call radiation_solve_time(state, &
                                         np_radmats_ii(p), &
                                         np_radmats(p))

            end if solve_time
      
         end if which_solve

      end do particle_type_loop 
            
   end subroutine radiation_solve

   ! --------------------------------------------------------------------------

end module radiation
