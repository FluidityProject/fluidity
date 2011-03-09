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

module radiation_normalise_flux

   !!< This module contains procedures associated with normalising the radiation particle flux 
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   
   use radiation_materials

   implicit none
   
   private 

   public :: normalise_particle_flux

contains

   ! --------------------------------------------------------------------------

   subroutine normalise_particle_flux(state,np_radmat) 
   
      !!< Normalise the particle flux

      type(state_type), intent(inout) :: state
      type(np_radmat_type), intent(in) :: np_radmat
      
      
      
   end subroutine normalise_particle_flux

   ! --------------------------------------------------------------------------

end module radiation_normalise_flux
