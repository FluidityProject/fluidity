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
      SUBROUTINE INITIM
C-----------------------------------------------------------------------
C
C
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blktim.i'
C
      TOTADD = 0
      TOTSUB = 0
      TOTCHG = 0
      TOTSPL = 0
      TOTSHR = 0
      TOTE2F = 0
      TOTF2E = 0
      TOTE2E = 0
      TOTNOD = 0
      TOTEDE = 0
      TOTNCN = 0
      TOTFND = 0
C
      TIMTOT = 0.0
      TIMFND = 0.0
      TIMADD = 0.0
      TIMSUB = 0.0
      TIMCHG = 0.0
      TIMSPL = 0.0
      TIMSHR = 0.0
      TIME2F = 0.0
      TIMF2E = 0.0
      TIME2E = 0.0
      TIMNOD = 0.0
      TIMEDE = 0.0
      TIMNCN = 0.0
C
      RETURN
      END
C
