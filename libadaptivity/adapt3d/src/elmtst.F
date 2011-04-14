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
      SUBROUTINE ELMTST( BIGLST, NODLST,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   IELM, NODS, CURDIF, CURORG, AVEDIF, AVEORG )
C-----------------------------------------------------------------------
C
C - This subroutine works out the element functionals produced when
C - splitting an element (IELM) into four elements, with a new node
C - at the average of the original element's four nodes.
C
C *** NOT CURRENTLY USED BY ADAPTIVITY *** - see lines 257-261 adptvy.F
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
      INCLUDE 'blkorg.i'
C
      INTEGER IELM, NODS(4), EDGLST(6)
C
      REAL CURDIF, CURORG, AVEDIF, AVEORG
C
      INTEGER INOD, I, INEW, INELM
C
      REAL EDGFNC, ELMFNC, U, U1, U2, U3, U4, D1, D2, D3, D4, X, Y, Z,
     :     F0, F1, F2, F3, F4, F5
C
      CURORG = -1E+30
      CURDIF = 1E+31
      AVEORG = -1E+30
      AVEDIF = 1E+31
C
      IF( IERR .NE. 0 ) RETURN
C
C - current functional contribution of this element
C
      F0 = ELMFNC( BIGLST, NODLST, IELM,
     :             NODS(1), NODS(2), NODS(3), NODS(4), U )
C
C - don't divide if element is already very small
C
c      IF( U .LT. 0.04 ) GOTO 900
C
      X = 0.0
      Y = 0.0
      Z = 0.0
C
      DO I = 1, 4
         X = X + NODLST(1,NODS(I))
         Y = Y + NODLST(2,NODS(I))
         Z = Z + NODLST(3,NODS(I))
      END DO
C
      X = X / 4
      Y = Y / 4
      Z = Z / 4
C
      INEW = -1
      INELM = INT(NODLST(16,NODS(1)))  ! start search from this elem
      CALL ADDND( BIGLST, NODLST,
     :            ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :            SZENLS, SZNELS, NNOD,   NELM,
     :            NODX,   NODY,   NODZ,   ORGMTX,
     :            X, Y, Z, INOD,  INEW,   INELM,  .TRUE. )
C
      IF( IERR .NE. 0 ) RETURN
C
C - functional contribution of new elements
C
      F2 = ELMFNC( BIGLST, NODLST, 0,
     :             NODS(1), NODS(2), NODS(3), INOD, U1 )
      F2 = MAX( F2, ELMFNC( BIGLST, NODLST, 0,
     :                      NODS(1), NODS(2), NODS(4), INOD, U2 ) )
      F2 = MAX( F2, ELMFNC( BIGLST, NODLST, 0,
     :                      NODS(1), NODS(3), NODS(4), INOD, U3 ) )
      F2 = MAX( F2, ELMFNC( BIGLST, NODLST, 0,
     :                      NODS(2), NODS(3), NODS(4), INOD, U4 ) )
C
      CALL RMVND( BIGLST, NODLST, INOD )
c      goto 900
C
C - don't want to create very small elements
C
c      IF( ABS(U1) .LT. 0.05 ) GOTO 900
c      IF( ABS(U2) .LT. 0.05 ) GOTO 900
c      IF( ABS(U3) .LT. 0.05 ) GOTO 900
c      IF( ABS(U4) .LT. 0.05 ) GOTO 900
C
      CURORG = F0
      CURDIF = F2 - CURORG
c
c      print*,'Element split functional: ',curdif,curorg
C
C - force division if element is far too large
C
c      IF( U .GT. 3.0 ) THEN
c         print*,'Forcing split of element ',ielm
c         CURDIF = -1E+25
c         CURORG = 0.0
c      END IF
C
 900  CONTINUE
C
      RETURN
      END
C
