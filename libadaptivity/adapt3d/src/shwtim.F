C Copyright (C) 2006 Imperial College London and others.
C 
C Please see the AUTHORS file in the main source directory for a full list
C of copyright holders.
C 
C Adrian Umpleby
C Applied Modelling and Computation Group
C Department of Earth Science and Engineering
C Imperial College London
C 
C adrian@Imperial.ac.uk
C 
C This library is free software; you can redistribute it and/or
C modify it under the terms of the GNU Lesser General Public
C License as published by the Free Software Foundation; either
C version 2.1 of the License.
C 
C This library is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Lesser General Public License for more details.
C 
C You should have received a copy of the GNU Lesser General Public
C License along with this library; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C USA
#include "ewrite.h"
      SUBROUTINE SHWTIM
C-----------------------------------------------------------------------
C
C
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      INCLUDE 'blktim.i'
      INCLUDE 'blkerr.i'
C
      ewrite(3,*) 
     :     ' '
      ewrite(3,*) 
     :     '-------------------------------------------------'
      ewrite(3,*) 
     :     ' '
      ewrite(3,*) 
     :     'TIME FOR THIS ADAPT: ',TIMTOT
      ewrite(3,*) 
     :     ' '
      IF( TOTEDE .GT. 0 ) then
         ewrite(3,*) 
     :        ' EDGELS: ',TOTEDE,TIMEDE
      end if
      IF( TOTNCN .GT. 0 ) then
         ewrite(3,*) 
     :        ' NDCNCT: ',TOTNCN,TIMNCN
      end if
      IF( TOTSPL .GT. 0 ) then
         ewrite(3,*) 
     :     ' TSTSPL: ',TOTSPL,TIMSPL
      end if
      IF( TOTSHR .GT. 0 ) then
         ewrite(3,*) 
     :     ' TSTSHR: ',TOTSHR,TIMSHR
      end if
      IF( TOTE2F .GT. 0 ) then
         ewrite(3,*) 
     :        ' TSTE2F: ',TOTE2F,TIME2F
      end if
      IF( TOTF2E .GT. 0 ) then 
         ewrite(3,*) ' TSTF2E: ',TOTF2E,TIMF2E
      end if
      IF( TOTE2E .GT. 0 ) then
         ewrite(3,*) ' TSTE2E: ',TOTE2E,TIME2E
      end if
      IF( TOTNOD .GT. 0 ) then
         ewrite(3,*) ' TSTNOD: ',TOTNOD,TIMNOD
      end if

      ewrite(3,*) 
     :     ' CHANGE: ',TOTCHG,TIMCHG
      ewrite(3,*) 
     :     '  ADDEL: ',TOTADD,TIMADD
      ewrite(3,*) 
     :     '  SUBEL: ',TOTSUB,TIMSUB
      IF( TOTFND .GT. 0 ) then
         ewrite(3,*) ' FNDELM: ',TOTFND,TIMFND
      end if

      ewrite(3,*) ' '
      ewrite(3,*) 
     :     '-------------------------------------------------'
      ewrite(3,*) ' '
C
      RETURN
      END
C
