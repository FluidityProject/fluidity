!     Copyright (C) 2006 Imperial College London and others.
!     
!     Please see the AUTHORS file in the main source directory for a full list
!     of copyright holders.
!     
!     Prof. C Pain
!     Applied Modelling and Computation Group
!     Department of Earth Science and Engineering
!     Imperial College London
!     
!     C.Pain@Imperial.ac.uk
!     
!     This library is free software; you can redistribute it and/or
!     modify it under the terms of the GNU Lesser General Public
!     License as published by the Free Software Foundation,
!     version 2.1 of the License.
!     
!     This library is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!     Lesser General Public License for more details.
!     
!     You should have received a copy of the GNU Lesser General Public
!     License along with this library; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!     USA

#include "fdebug.h"

module AdvectionDiffusion

  use fields
  use global_parameters, only : OPTION_PATH_LEN
  use spud
  use state_module

  implicit none
  
  private
  
  public :: get_copied_field
  
contains

  subroutine get_copied_field(fieldname, state)

    type(state_type), intent(in) :: state
    character(len=*), intent(in) :: fieldname

    type(scalar_field), pointer :: copiedfield
    type(scalar_field), pointer :: tmpfield
    character(len=OPTION_PATH_LEN) :: tmpstring

    if(trim(fieldname)=="CopiedField") then
      copiedfield=>extract_scalar_field(state, "CopiedField")
      call get_option(trim(copiedfield%option_path)//"/prognostic/copy_from_field", &
                              tmpstring)
      tmpfield=>extract_scalar_field(state, "Old"//trim(tmpstring))
      call set(copiedfield, tmpfield)
    end if

  end subroutine get_copied_field

end module AdvectionDiffusion
