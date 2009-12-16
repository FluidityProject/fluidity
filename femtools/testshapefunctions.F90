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
program testshapefunctions
  use shape_functions
  implicit none

  type(element_type) :: element
  type(quadrature_type) :: quad
  integer :: dim, deg
  character(len=42) :: messg

  do dim=1,3
     do deg=0,7
        quad=make_quadrature(loc=dim+1, dimension=dim, degree=max(deg,1))
     
        write (messg, '(a,i0,a,i0)') "Testing ",dim,"D P",deg
        
        ewrite(3,*)  trim(messg)
        element=make_element_shape&
             (loc=dim+1, dimension=dim, degree=deg, quad=quad)       

     end do
  end do

  do deg=1,1
     
     quad=make_quadrature(loc=3, dimension=2, degree=max(deg,1))

     ewrite(3,*)  'Testing P1_NC'
     element=make_element_shape(loc=3, dimension=2, degree=deg, quad=quad,&
          & type=ELEMENT_NONCONFORMING)  

  end do

end program testshapefunctions
