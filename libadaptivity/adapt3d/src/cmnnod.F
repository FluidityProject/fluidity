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
      SUBROUTINE CMNNOD( BIGLST, NODLST, IEDG1, IEDG2,
     :                   INOD, INOD1, INOD2 )
C-----------------------------------------------------------------------
C
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
      INTEGER IEDG1, IEDG2, INOD, INOD1, INOD2
C
      INTEGER NOD1, NOD2
C
      INOD  = 0
      INOD1 = 0
      INOD2 = 0
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IEDG1 .LE. 0 ) THEN
         IERR = -151
         PRINT*,'*** CMNNOD: INVALID POINTER FOR EDGE 1'
         RETURN
      END IF
C
      IF( IEDG2 .LE. 0 ) THEN
         IERR = -151
         PRINT*,'*** CMNNOD: INVALID POINTER FOR EDGE 2'
         RETURN
      END IF
C
      NOD1 = BIGLST(1,IEDG1)
      NOD2 = BIGLST(2,IEDG1)
C
      IF( BIGLST(1,IEDG2) .EQ. NOD1 ) THEN
         INOD  = NOD1
         INOD1 = NOD2
         INOD2 = BIGLST(2,IEDG2)
      ELSE IF( BIGLST(2,IEDG2) .EQ. NOD1 ) THEN
         INOD  = NOD1
         INOD1 = NOD2
         INOD2 = BIGLST(1,IEDG2)
      ELSE IF( BIGLST(1,IEDG2) .EQ. NOD2 ) THEN
         INOD  = NOD2
         INOD1 = NOD1
         INOD2 = BIGLST(2,IEDG2)
      ELSE IF( BIGLST(2,IEDG2) .EQ. NOD2 ) THEN
         INOD  = NOD2
         INOD1 = NOD1
         INOD2 = BIGLST(1,IEDG2)
      END IF
C
      RETURN
      END
C
