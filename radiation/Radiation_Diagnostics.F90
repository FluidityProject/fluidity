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

module radiation_diagnostics
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use Fields
   use diagnostic_variables   
   
   implicit none
   
   private 

   public :: radiation_diagnostics_register_diagnostic, &
             radiation_eigenvalue_set_diagnostics
        
contains

   ! --------------------------------------------------------------------------

   subroutine radiation_diagnostics_register_diagnostic  
      
      !!< Register the radiation diagnostics for stat file 
      !!< for all particle types
      
      ! local variables
      integer :: p
      character(len=OPTION_PATH_LEN) :: particle_option_path         
      character(len=OPTION_PATH_LEN) :: particle_name   
      character(len=OPTION_PATH_LEN) :: equation_type

      have_radiation: if (have_option('/embedded_models/radiation')) then
      
         particle_loop: do p = 1,option_count('/embedded_models/radiation/particle_type')
         
            ! set the particle path
            particle_option_path = '/embedded_models/radiation/particle_type['//int2str(p-1)//']'
         
            call get_option(trim(particle_option_path)//'/equation/name',equation_type)  
      
            eig_register: if (trim(equation_type) == 'Eigenvalue') then
         
               ! get the particle name
               call get_option(trim(particle_option_path)//'/name',particle_name)
         
               call radiation_eigenvalue_register_diagnostics(trim(particle_name))
      
            else eig_register
      
               ! fill in ...
            
            end if eig_register
      
         end do particle_loop

      end if have_radiation
           
   end subroutine radiation_diagnostics_register_diagnostic

   ! --------------------------------------------------------------------------
   
   subroutine radiation_eigenvalue_register_diagnostics(particle_name)
      
      !!< Register the diagnostics variables for the stat file
      !!< associated with the eigenvalue run for this particle type given by particle_name 
      
      character(len=*), intent(in) :: particle_name
      
      ewrite(1,*) 'Radiation eigenvalue run register diagnostics for particle type ',trim(particle_name)
      
      call register_diagnostic(dim       = 1, &
                               name      = 'ParticleKeff'//trim(particle_name), &
                               statistic = 'Value')
      
   end subroutine radiation_eigenvalue_register_diagnostics
   
   ! --------------------------------------------------------------------------
   
   subroutine radiation_eigenvalue_set_diagnostics(state, &
                                                   particle_name)
         
      !!< Set the diagnostics variables for the stat file
      !!< associated with the eigenvalue run for this particle given by particle_name 

      type(state_type), intent(in) :: state      
      character(len=*), intent(in) :: particle_name
      
      ! local variables
      type(scalar_field), pointer:: keff_field      
      real :: keff

      ewrite(1,*) 'Radiation eigenvalue set registered diagnostics for particle type ',trim(particle_name)
      
      keff_field => extract_scalar_field(state, &
                                         'ParticleKeff'//trim(particle_name))
      
      keff = node_val(keff_field,1)
 
      call set_diagnostic(name      = 'ParticleKeff'//trim(particle_name),  &
                          statistic = 'Value', &
                          value     = (/keff/))
      
   end subroutine radiation_eigenvalue_set_diagnostics
   
   ! --------------------------------------------------------------------------

end module radiation_diagnostics
