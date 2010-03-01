!    Copyright (C) 2010 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Dr Gerard Gorman
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    g.gorman@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
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

module Profiler
  use global_parameters, only : real_8
  implicit none
contains
  
  real(kind = real_8) function profiler_get(key)
    character(len=*), intent(in)::key
    call cprofiler_get(key, len_trim(key), profiler_get)
  end function profiler_get

  subroutine profiler_tic(key)
    character(len=*), intent(in)::key
    call cprofiler_tic(key, len_trim(key))
  end subroutine profiler_tic

  subroutine profiler_toc(key)
    character(len=*), intent(in)::key
    call cprofiler_toc(key, len_trim(key))
  end subroutine profiler_toc

  subroutine profiler_zero()
    call cprofiler_zero()
  end subroutine profiler_zero

end module Profiler
