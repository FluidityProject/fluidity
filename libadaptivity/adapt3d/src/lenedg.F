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
      REAL FUNCTION LENEDG( BIGLST, NODLST, IEDG )
C-----------------------------------------------------------------------
C
C - This function 
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER IEDG
C
      INTEGER NODS(2), INEW, IFLAG, NXT, LST
C
      REAL DX, DY, DZ, DOT
C
      LENEDG = 0.0
C
      CALL EDGINF( BIGLST, NODLST, IEDG, NODS,
     :             INEW, IFLAG, NXT, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( NODS(1) .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** LENEDG: INVALID NODE POINTER 1'
         RETURN
      END IF
C
      IF( NODS(2) .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** LENEDG: INVALID NODE POINTER 2'
         RETURN
      END IF
C
      DX = NODLST(1,NODS(1)) - NODLST(1,NODS(2))
      DY = NODLST(2,NODS(1)) - NODLST(2,NODS(2))
      DZ = NODLST(3,NODS(1)) - NODLST(3,NODS(2))
C
      DOT = DX*DX + DY*DY + DZ*DZ
C
      LENEDG = SQRT(DOT)
C
      RETURN
      END
C
