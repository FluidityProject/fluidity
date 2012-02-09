! Copyright (C) 2006 Imperial College London and others.
! 
! Please see the AUTHORS file in the main source directory for a full list
! of copyright holders.
! 
! Adrian Umpleby
! Applied Modelling and Computation Group
! Department of Earth Science and Engineering
! Imperial College London
! 
! adrian@Imperial.ac.uk
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA

#include "fdebug_adapt.h"

      module addebug
  
        use write_log
  
        implicit none
  
        private
  
        public :: adabort_pinpoint
  
      contains

        SUBROUTINE adabort_pinpoint(ErrorStr, FromFile, LineNumber)

          CHARACTER*(*) ErrorStr, FromFile
          INTEGER LineNumber
    
          ewrite(-1,FMT='(A)') "*** LIBADAPT ERROR ***"
          ewrite(-1,FMT='(3A,I5,A)') "Source location: (",FromFile,",",LineNumber,")"
          ewrite(-1,FMT='(2A)') "Error message: ",ErrorStr
          ewrite(-1,FMT='(A)') "Error is terminal."
    
          STOP
        END SUBROUTINE adabort_pinpoint
  
      end module addebug
