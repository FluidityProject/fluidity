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
      SUBROUTINE CHGEDS( BIGLST, NODLST, IELM, ETR1, ETR2, ETA1, ETA2 )
C-----------------------------------------------------------------------
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
      INTEGER IELM, ETR1, ETR2, ETA1, ETA2
C
      INTEGER NXT, NXT2, I
C
      IF( IELM .LE. 0 ) THEN
         PRINT*,'*** CHGEDS: INVALID ELEMENT'
         print*,ielm,ETR1,ETR2,ETA1,ETA2
         IERR = -1
         RETURN
      END IF
C
      IF( ETR1 .LE. 0 ) THEN
         PRINT*,'*** CHGEDS: INVALID EDGE'
         print*,ielm,ETR1,ETR2,ETA1,ETA2
         IERR = -2
         RETURN
      END IF
C
      IF( ETR2 .LE. 0 ) THEN
         PRINT*,'*** CHGEDS: INVALID EDGE'
         print*,ielm,ETR1,ETR2,ETA1,ETA2
         IERR = -2
         RETURN
      END IF
C
      IF( ETA1 .LE. 0 ) THEN
         PRINT*,'*** CHGEDS: INVALID EDGE'
         print*,ielm,ETR1,ETR2,ETA1,ETA2
         IERR = -2
         RETURN
      END IF
C
      IF( ETA2 .LE. 0 ) THEN
         PRINT*,'*** CHGEDS: INVALID EDGE'
         print*,ielm,ETR1,ETR2,ETA1,ETA2
         IERR = -2
         RETURN
      END IF
C
      NXT = BIGLST(NXTBIG,IELM)
      NXT2 = BIGLST(NXTBIG,NXT)
C
      print*,'changing eds for elm ',ielm
      print*,etr1,etr2,eta1,eta2
      print*,(biglst(i,nxt),i=1,4),(biglst(i,nxt2),i=1,2)
c
      DO I = 1, 4
         IF(BIGLST(I,NXT) .EQ. ETR1) THEN
            BIGLST(I,NXT) = ETA1
         ELSE IF(BIGLST(I,NXT) .EQ. ETR2 ) THEN
            BIGLST(I,NXT) = ETA2
         END IF
      END DO
C
      DO I = 1, 2
         IF( BIGLST(I,NXT2) .EQ. ETR1 ) THEN
            BIGLST(I,NXT2) = ETA1
         ELSE IF( BIGLST(I,NXT2) .EQ. ETR2 ) THEN
            BIGLST(I,NXT2) = ETA2
         END IF
      END DO
C
      RETURN
      END
C
