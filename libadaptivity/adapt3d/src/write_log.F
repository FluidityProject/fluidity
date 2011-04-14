C     Copyright (C) 2006 Imperial College London and others.
C     
C     Please see the AUTHORS file in the main source directory for a full list
C     of copyright holders.
C     
C     Adrian Umpleby
C     Applied Modelling and Computation Group
C     Department of Earth Science and Engineering
C     Imperial College London
C     
C     adrian@Imperial.ac.uk
C     
C     This library is free software; you can redistribute it and/or
C     modify it under the terms of the GNU Lesser General Public
C     License as published by the Free Software Foundation; either
C     version 2.1 of the License.
C     
C     This library is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C     Lesser General Public License for more details.
C     
C     You should have received a copy of the GNU Lesser General Public
C     License along with this library; if not, write to the Free Software
C     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C     USA

      module write_log
      implicit none
      integer::cdl=3
C     These defaults are an assumption which may not be true on all platforms.
      integer, save :: debug_error_unit=0, debug_log_unit=6

      contains

      function dunit(priority)
C     Decide where to send output based on the level of the error.
      implicit none
      integer :: dunit
      integer, intent(in) :: priority
      
      if (priority<1) then
         dunit=debug_error_unit
      else
         dunit=debug_log_unit
      end if
      end function dunit

      subroutine get_debug_level(level)
      integer, intent(out)::level
      level = cdl
      end subroutine get_debug_level
      
      subroutine set_debug_level(level)
      integer, intent(in)::level
      cdl = level
      end subroutine set_debug_level
      end module write_log
      
