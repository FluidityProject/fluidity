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
      SUBROUTINE INIRLS( EMTSTT, THELST, SIZBLK, NXT, LST, MAXBLK )
C-----------------------------------------------------------------------
C
C - This subroutine initialises a linked block-list, ready to receive
C - new objects.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blkerr.i'
C
      INTEGER SIZBLK, MAXBLK, EMTSTT, NXT, LST
C
      REAL THELST(SIZBLK,MAXBLK)
C
      INTEGER I, J
c
c      print*,'Into INIRLS ',ierr
c      print*,'   SIZ,MAX: ',SIZBLK,MAXBLK
c      PRINT*,'   NXT,LST: ',NXT,LST
C
      IF( IERR .NE. 0 ) RETURN
C
      EMTSTT = 1
C
      DO I = 1, MAXBLK
C
         DO J = 1, SIZBLK
            THELST(J,I) = 0
         END DO
C
         THELST(LST,I) = FLOAT(I-1)
         THELST(NXT,I) = FLOAT(I+1)
C
      END DO
C
      THELST(LST,1)      = -1.0
      THELST(NXT,MAXBLK) = -1.0
C
      RETURN
      END
C
