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

module porous_media
   
   !!< Module containing porous media related procedures.

   use fldebug
   use global_parameters, only: OPTION_PATH_LEN
   use futils
   use spud
   use fields
   use state_module

   implicit none
   
   private

   public :: form_porosity_theta
  
contains

! ----------------------------------------------------------------------------------

   subroutine form_porosity_theta(porosity_theta, state, option_path)
      
      !!< Form the porosity theta averaged field determined by the input option path
      
      type(scalar_field), intent(inout) :: porosity_theta
      type(state_type), intent(in) :: state
      character(len=*) :: option_path

      ! local variables
      type(scalar_field), pointer :: porosity_old, porosity_new
      character(len=OPTION_PATH_LEN) :: porosity_name
      real :: porosity_theta_value
      integer :: stat
      
      call get_option(trim(option_path)//'/porosity_field_name', &
                      porosity_name, &
                      default = 'Porosity')
         
      ! get the porosity theta value
      call get_option(trim(option_path)//'/porosity_temporal_theta', &
                      porosity_theta_value, &
                      default = 0.0)
         
      porosity_new => extract_scalar_field(state, trim(porosity_name), stat = stat)                  
  
      if (stat /=0) then 
         FLExit('Failed to extract Porosity from state to be used for forming the theta averaged porosity field')
      end if
         
      porosity_old => extract_scalar_field(state, "Old"//trim(porosity_name), stat = stat)

      if (stat /=0) then 
         FLExit('Failed to extract OldPorosity from state to be used for forming the theta averaged porosity field')
      end if
         
      call allocate(porosity_theta, porosity_new%mesh)
         
      call set(porosity_theta, porosity_new, porosity_old, porosity_theta_value)
         
      ewrite_minmax(porosity_theta)
      
   end subroutine form_porosity_theta
      
! ----------------------------------------------------------------------------------

end module porous_media

