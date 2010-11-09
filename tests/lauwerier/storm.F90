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
#include "confdefs.h"
#include "fdebug.h"

program create_stress
  use NetCDFWriter
  implicit none
  
  ! integer times(1)
  real latitude(181), longitude(361), variable1(361,181,1), variable2(361,181,1)
  character(len=256)::name1="nsss", long_name1="North-South surface stress", units1="N m**-2 s"
  character(len=256)::name2="ewss", long_name2="East-West surface stress", units2="N m**-2 s"
  integer i, j
  integer::t=1
  real::r
  
  do j=0, 180
     latitude(j+1) = j - 90.0
     do i=0, 360
        longitude(i+1) = i - 180
        
        variable1(i+1, j+1, t) = sqrt(float(i*i+j*j))
        variable2(i+1, j+1, t) = sqrt(float((i-45)*(i-45)+j*j))
     end do
  end do
  
  call NetCDFWriter_init("WindStress.nc", longitude, latitude)
  call NetCDFWriter_write_variable(name1, long_name1, variable1, units1)
  call NetCDFWriter_write_variable(name2, long_name2, variable2, units2)

  return
end program create_stress
