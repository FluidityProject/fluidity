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
      SUBROUTINE SHWHST( BIGLST, NODLST, OUTFL1, OUTFL2, OUTFL3 )
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      INTEGER DIVS
      PARAMETER( DIVS = 81 )
C
      INTEGER OUTFL1, OUTFL2, OUTFL3
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER IPTR, NXT, LST, NODS(4), IFLAG, INEW, I, J, CNT(3,DIVS)
C
      REAL V, F, ELMFNC, EDGFNC, X(DIVS)
C
      LOGICAL LSTEDG, GETFLG
C
      LSTEDG = EDGON
C
      DO I = 1, DIVS
         X(I) = 0.05*(I-1)
         CNT(1,I) = 0
         CNT(2,I) = 0
         CNT(3,I) = 0
      END DO
C
      IPTR = STTBIG
C
  10  IF( BIGLST(4,IPTR) .LT. 0 ) THEN
C
         NXT = BIGLST(NXTBIG,IPTR)
C
         EDGON = .TRUE.
C
         F = EDGFNC( BIGLST, NODLST, BIGLST(1,IPTR), BIGLST(2,IPTR), V )
C
         I = DIVS
  20     IF( I .GT. 1 ) THEN
            I = I - 1
            IF( V .LT. X(I) ) GOTO 20
            CNT(1,I+1) = CNT(1,I+1) + 1
         ELSE
            print*,'1: AAAAAAAARRRGH! ',V,x
            stop
         END IF
C
      ELSE IF( BIGLST(4,IPTR) .GT. 0 ) THEN
C
         CALL ELNODS( BIGLST, NODLST, IPTR, NODS, NXT, LST, .false. )
C
         EDGON = .FALSE.
C
         F = ELMFNC( BIGLST, NODLST, iptr,
     :               NODS(1), NODS(2), NODS(3), NODS(4), V )
         V = ABS(V)*SQRT(24.0)
C
         I = DIVS
  30     IF( I .GT. 1 ) THEN
            I = I - 1
            IF( V .LT. X(I) ) GOTO 30
            CNT(2,I+1) = CNT(2,I+1) + 1
         ELSE
            print*,'2: AAAAAAAARRRGH! ',V,x
            stop
         END IF
C
         I = DIVS
  40     IF( I .GT. 1 ) THEN
            I = I - 1
            IF( F .LT. X(I) ) GOTO 40
            CNT(3,I+1) = CNT(3,I+1) + 1
         ELSE
            print*,'3: AAAAAAAARRRGH! ',F,x
            stop
         END IF
C
      END IF
C
      IPTR = NXT
C
      IF( GETFLG(BIGLST,NODLST,IPTR,1) ) GOTO 10
C
      EDGON = LSTEDG
C
C - now dump the results to the files
C
      DO I = 1, DIVS-1
         ewrite(3,900) X(I),CNT(1,I+1)
         ewrite(3,900) X(I),CNT(2,I+1)
         ewrite(3,900) X(I),CNT(3,I+1)
         ewrite(3,900) X(I+1),CNT(1,I+1)
         ewrite(3,900) X(I+1),CNT(2,I+1)
         ewrite(3,900) X(I+1),CNT(3,I+1)
      END DO
c
  900 format( 1pe15.7, i7 )
C
      RETURN
      END
C
