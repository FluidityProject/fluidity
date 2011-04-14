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
      INTEGER FUNCTION WCHFAC( BIGLST, NODLST, IELM, IED1, IED2, IED3 )
C-----------------------------------------------------------------------
C
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
      INTEGER IELM, IFAC, IED1, IED2, IED3
C
      INTEGER I, NXT, NXT2, LED1, LED2, LED3, LSUM
C
      LOGICAL GOTONE
C
      WCHFAC = 0
      IF( IERR .NE. 0 ) RETURN
C
      IF( IELM .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** WCHFAC: INVALID ELEMENT POINTER GIVEN'
         print*,ielm, ied1, ied2, ied3
         RETURN
      END IF
C
      IF( IED1 .LE. 0 .OR. IED2 .LE. 0 .OR. IED3 .LE. 0 ) THEN
         IERR = -2
         PRINT*,'*** WCHFAC: INVALID EDGE POINTER GIVEN'
         print*,ielm, ied1, ied2, ied3
         RETURN
      END IF
C
      IF( BIGLST(4,IELM) .LT. 0 ) THEN
         IERR = -3
         PRINT*,'*** WCHFAC: FIRST BLOCK IS NOT AN ELEMENT'
         print*,ielm,(biglst(nxt,ielm),nxt=1,4)
         RETURN
      END IF
C
      NXT = BIGLST(NXTBIG,IELM)
C
      IF( NXT .LE. 0 ) THEN
         IERR = -4
         PRINT*,'*** WCHFAC: INVALID POINTER TO SECOND BLOCK'
         print*,ielm, ied1, ied2, ied3, nxt
         RETURN
      END IF
C
      IF( BIGLST(4,NXT) .LT. 0 ) THEN
         IERR = -5
         PRINT*,'*** WCHFAC: SECOND BLOCK IS NOT AN ELEMENT'
         print*,ielm,nxt,biglst(4,nxt)
         RETURN
      END IF
C
      LED1 = 0
      LED2 = 0
      LED3 = 0
C
      DO I = 1, 4
         IF( IED1 .EQ. BIGLST(I,NXT) ) LED1 = I
         IF( IED2 .EQ. BIGLST(I,NXT) ) LED2 = I
         IF( IED3 .EQ. BIGLST(I,NXT) ) LED3 = I
      END DO
C
      IF( LED1 .EQ. 0 .OR. LED2 .EQ. 0 .OR. LED3 .EQ. 0 ) THEN
C
         NXT2 = BIGLST(NXTBIG,NXT)
C
         IF( NXT2 .LE. 0 ) THEN
            IERR = -6
            PRINT*,'*** WCHFAC: INVALID POINTER TO THIRD BLOCK'
            print*,ielm,ied1,ied2,ied3,nxt,nxt2
            RETURN
         END IF
C
         IF( BIGLST(4,NXT2) .LT. 0 ) THEN
            IERR = -7
            PRINT*,'*** WCHFAC: THIRD BLOCK IS NOT AN ELEMENT'
            print*,ielm,ied1,ied2,ied3,nxt,nxt2,biglst(4,nxt2)
            RETURN
         END IF
C
         IF( LED1 .EQ. 0 ) THEN
            IF( IED1 .EQ. BIGLST(1,NXT2) ) THEN
               LED1 = 5
            ELSE IF( IED1 .EQ. BIGLST(2,NXT2) ) THEN
               LED1 = 6
            END IF
         END IF
C
         IF( LED2 .EQ. 0 ) THEN
            IF( IED2 .EQ. BIGLST(1,NXT2) ) THEN
               LED2 = 5
            ELSE IF( IED2 .EQ. BIGLST(2,NXT2) ) THEN
               LED2 = 6
            END IF
         END IF
C
         IF( LED3 .EQ. 0 ) THEN
            IF( IED3 .EQ. BIGLST(1,NXT2) ) THEN
               LED3 = 5
            ELSE IF( IED3 .EQ. BIGLST(2,NXT2) ) THEN
               LED3 = 6
            END IF
         END IF
C
      ELSE
C
         NXT2 = 0
C
      END IF
C
      IF( LED1 .EQ. 0 .OR. LED2 .EQ. 0 .OR. LED3 .EQ. 0 ) THEN
         IERR = -5431
         PRINT*,'*** WCHFAC: Failed to find all edges in element!'
         print*,'elem & edges: ',ielm,ied1,ied2,ied3
         print*,'local positions found: ',led1,led2,led3
         NXT2 = BIGLST(NXTBIG,NXT)
         print*,'conelm: ',(biglst(i,ielm),i=1,3),biglst(4,ielm)/2
         print*,'edglst: ',(biglst(i,nxt),i=1,4),
     :                      biglst(1,nxt2),biglst(2,nxt2)
         RETURN
      ENDIF
C
      LSUM = LED1 + LED2 + LED3
C
      GOTONE = (LED1 .EQ. 1) .OR. (LED2 .EQ. 1) .OR. (LED3 .EQ. 1)
C
      IF( GOTONE ) THEN
C
         IF( LSUM .EQ. 6 ) THEN
C
            WCHFAC = 1
C
         ELSE IF( LSUM .EQ. 11 ) THEN
C
            WCHFAC = 2
C
         ELSE
C
            IERR = -5432
            PRINT*,'*** WCHFAC: Failed to find face for edges (1)'
            print*,ielm,ied1,ied2,ied3
            NXT2 = BIGLST(NXTBIG,NXT)
            print*,(biglst(i,nxt),i=1,4),biglst(1,nxt2),biglst(2,nxt2)
C
         END IF
C
      ELSE IF( LSUM .EQ. 11 ) THEN
C
         WCHFAC = 4
C
      ELSE IF( LSUM .EQ. 14 ) THEN
C
         WCHFAC = 3
C
      ELSE
C
         IERR = -5433
         PRINT*,'*** WCHFAC: Failed to find face for edges (2)'
         print*,ielm,ied1,ied2,ied3
         NXT2 = BIGLST(NXTBIG,NXT)
         print*,(biglst(i,nxt),i=1,4),biglst(1,nxt2),biglst(2,nxt2)
C
      END IF
C
      RETURN
      END
C
