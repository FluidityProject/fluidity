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
!    C.Pain@Imperial.ac.uk
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

module legacy_tensor_fields
  use fields
  implicit none

  private

  public copy_tensor_from_legacy, copy_tensor_to_legacy

contains
  
  subroutine copy_tensor_from_legacy(field, xx, xy, xz, yy, yz, zz)
    !!< Copy the contents of the symmetric tensor components xx..zz into the
    !!< tensor field provided.
    type(tensor_field), intent(inout) :: field
    real, dimension(:), intent(in) :: xx, xy, xz, yy, yz, zz
    
    field%val(1,1,:) = xx
    
    if (field%dim>1) then
       field%val(1,2,:) = xy
       field%val(2,1,:) = xy
       field%val(2,2,:) = yy
       
       if (field%dim>2) then
          field%val(1,3,:) = xz
          field%val(3,1,:) = xz
          field%val(2,3,:) = yz
          field%val(3,2,:) = yz
          field%val(3,3,:) = zz
       end if
       
    end if
    
  end subroutine copy_tensor_from_legacy

  subroutine copy_tensor_to_legacy(field, xx, xy, xz, yy, yz, zz)
    !!< Copy the contents of the symmetric tensor field provided into xx..zz
    type(tensor_field), intent(in) :: field
    real, dimension(:), intent(inout) :: xx, xy, xz, yy, yz, zz
    
    xx=field%val(1,1,:)
    
    if (field%dim>1) then
       xy = field%val(1,2,:) 
       yy = field%val(2,2,:)
       
       if (field%dim>2) then
          xz=field%val(1,3,:)
          yz=field%val(2,3,:)
          zz=field%val(3,3,:)
       end if

    end if

  end subroutine copy_tensor_to_legacy

end module legacy_tensor_fields
