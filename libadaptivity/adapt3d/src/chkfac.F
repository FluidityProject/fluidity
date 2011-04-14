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
      SUBROUTINE CHKFAC( BIGLST, NODLST, N12, N1, N2,
     :                   F1, F2, F3, F4, F5, F6,
     :                   MXRAD, MXRND, MXRDS )
C-----------------------------------------------------------------------
C
C - This subroutine checks to see if a new face, defined by nodes N1,N2,N12)
C - is allowable.
C - This is done by looping through all front edges (via the faces) and
C - calling EDGDST to find the distance to an edge.
C - F1 to F6 are faces that should be excluded from checking: typically
C - they are faces that are joined to the current face, and any others that
C - may be close by that should be left unchecked.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkfrt.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER N12, N1, N2, MXRND, F1, F2, F3, F4, F5, F6
C
      REAL MXRAD, MXRDS
C
      INTEGER IFLAG, INEWED, NXT, LST, NOD1(2), NOD2(2), IELM,
     :        IPT, DSTCNT, EDGLST(3), FACLST(3), NXTFAC, INEWFC, I,
     :        IE1, IE2, IE3
C
      REAL DX1, DY1, DZ1, DX2, DY2, DZ2, DOT, ONX, ONY, ONZ,
     :     vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3,
     :     LN1, LN2, LN3, AREA, DST, MNDST, MNDST2, VAL1, VAL2, VAL3,
     :     DIRX, DIRY, DIRZ, LNGTH, X(3), Y(3), Z(3), EDGDST
C
      LOGICAL GETFLG
C
c      PRINT*,' '
c      print*,'Into CHKFAC... ',ierr
c      print*,'  N12,N1,N2:',N12,N1,N2
c      print*,'   F1,F2,F3:',F1,F2,F3
c      print*,'   F4,F5,F6:',F4,F5,F6
c      print*,'      MXRDS:',MXRDS
C
      IF( IERR .NE. 0 ) RETURN
C
      MNDST = MXRDS
      MNDST2 = MNDST*MNDST
C
      CALL FACNRM( BIGLST, NODLST, N12, N1, N2, ONX, ONY, ONZ )
C
      IF( IERR .NE. 0 ) RETURN
C
      DIRX = NODLST(1,N2) - NODLST(1,N1)
      DIRY = NODLST(2,N2) - NODLST(2,N1)
      DIRZ = NODLST(3,N2) - NODLST(3,N1)
C
      LNGTH = DIRX*DIRX + DIRY*DIRY + DIRZ*DIRZ
C
      LNGTH = SQRT(LNGTH)
C
      DIRX = DIRX/LNGTH
      DIRY = DIRY/LNGTH
      DIRZ = DIRZ/LNGTH
C
      X(1) = NODLST(1,N12)
      Y(1) = NODLST(2,N12)
      Z(1) = NODLST(3,N12)
C
      X(2) = NODLST(1,N1)
      Y(2) = NODLST(2,N1)
      Z(2) = NODLST(3,N1)
C
      X(3) = NODLST(1,N2)
      Y(3) = NODLST(2,N2)
      Z(3) = NODLST(3,N2)
C
      GOTO 200
C
      IPT = STFRND
C
      MXRAD = -1E+30
C
      DO I = 1, NFRTND
C
C - position of this front node
C
         print*,'=== Checking front node ',ipt
C
         DX2 = NODLST(1,IPT)
         DY2 = NODLST(2,IPT)
         DZ2 = NODLST(3,IPT)
C
C - don't want nodes that belong to this face
C
         IF( IPT .EQ. N1 ) GOTO 100
         IF( IPT .EQ. N2 ) GOTO 100
         IF( IPT .EQ. N12 ) GOTO 100
C
C - vector (3) from front node to node N12
C
         VX3 = NODLST(1,N12) - DX2
         VY3 = NODLST(2,N12) - DY2
         VZ3 = NODLST(3,N12) - DZ2
C
C - distance from front node to plane of face
C
         DST = -VX3*ONX - VY3*ONY - VZ3*ONZ
C
C - if front node is far enough from plane of face then let it pass
C
         IF( ABS(DST) .GT. MNDST ) GOTO 100
C
         DSTCNT = 0
C
C - normalise the vector (3) from front node to common node
C
         DOT = VX3*VX3 + VY3*VY3 + VZ3*VZ3
C
         IF( DOT .GT. MNDST2 ) DSTCNT = 1
C
         DOT = 1.0/SQRT(DOT)
C
         VX3 = VX3*DOT
         VY3 = VY3*DOT
         VZ3 = VZ3*DOT
C
C  - vector (1) from front node to other node of EDG1
C
         VX1 = NODLST(1,N1) - DX2
         VY1 = NODLST(2,N1) - DY2
         VZ1 = NODLST(3,N1) - DZ2
C
C - normalise vector (1)
C
         DOT = VX1*VX1 + VY1*VY1 + VZ1*VZ1
C
         IF( DOT .GT. MNDST2 ) DSTCNT = DSTCNT + 1
C
         DOT = 1.0/SQRT(DOT)
C
C - vector (1) becomes difference of vector (3) and vector (1)
C
         VX1 = VX1*DOT - VX3
         VY1 = VY1*DOT - VY3
         VZ1 = VZ1*DOT - VZ3
C
C  - vector (2) from front node to other node of EDG2
C
         VX2 = NODLST(1,N2) - DX2
         VY2 = NODLST(2,N2) - DY2
         VZ2 = NODLST(3,N2) - DZ2
C
C - normalise vector (2)
C
         DOT = VX2*VX2 + VY2*VY2 + VZ2*VZ2
C
         IF( DOT .GT. MNDST2 ) DSTCNT = DSTCNT + 1
C
C - if the front node is far enough away from all three face nodes
C - then let it pass
C
         IF( DSTCNT .EQ. 3 ) GOTO 100
C
         DOT = 1.0/SQRT(DOT)
C
C - vector (2) becomes difference of vector (3) and vector (2)
C
         VX2 = VX2*DOT - VX3
         VY2 = VY2*DOT - VY3
         VZ2 = VZ2*DOT - VZ3
C
C - 'length-squared's of vectors (1&2) and the length^2 of their difference
C
         LN1 = VX1*VX1 + VY1*VY1 + VZ1*VZ1
         LN2 = VX2*VX2 + VY2*VY2 + VZ2*VZ2
         LN3 = (VX2-VX1)*(VX2-VX1) + (VY2-VY1)*(VY2-VY1)
     :       + (VZ2-VZ1)*(VZ2-VZ1)
C
C - cross product of vector (1) and vector (2) to find area
C
         DX1 = VY1*VZ2 - VZ1*VY2
         DY1 = VZ1*VX2 - VX1*VZ2
         DZ1 = VX1*VY2 - VY1*VX2
C
C - This is 4 times square of area of the triangle of points where the
C - vectors from the front node to the face corners hit the unit sphere
C - centred on the front node.
C
         AREA = DX1*DX1 + DY1*DY1 + DZ1*DZ1
C
C - This should (according to Del.Tri.& Mesh.) be the radius of the
C - circumcircle of the above triangle. This is supposed to give a
C - rough estimate of the acceptability of the front node being checked.
C - The radius should range from zero (the node is very distant) to pi
C - (the node is very close to the plane of the face).
C
         DOT = SQRT(LN1*LN2*LN3/AREA)/2
C
         print*,'      distance: ',dst
         print*,'        radius: ',dot
c
         IF( DOT .GT. MXRAD ) THEN
C
            MXRND = IPT
            MXRAD = DOT
            MXRDS = DST
C
         END IF
C
         GOTO 101
C
 100     CONTINUE
C
         PRINT*,'      No need to check it'
C
 101     CONTINUE
C
         LST = INT(NODLST(LSTNOD,IPT))
         IPT = INT(NODLST(NXTNOD,IPT))
C
      END DO
C
 200  CONTINUE
C
C - loop through faces making sure all edges have bit 4 set to zero
C
      IPT = STTFRT
C
      DO I = 1, NFACES
C
         CALL FACINF( IPT, FACLST, EDGLST, IELM,
     :                INEWFC, IFLAG, NXTFAC, LST )
C
         CALL SETFLG( BIGLST, NODLST, EDGLST(1), 4, .FALSE. )
         CALL SETFLG( BIGLST, NODLST, EDGLST(2), 4, .FALSE. )
         CALL SETFLG( BIGLST, NODLST, EDGLST(3), 4, .FALSE. )
C
         IPT = NXTFAC
C
      END DO
C
C - loop through faces looking at each edge that has not already been done
C
      IPT = STTFRT
C
      DO I = 1, NFACES
C
         CALL FACINF( IPT, FACLST, EDGLST, IELM,
     :                INEWFC, IFLAG, NXTFAC, LST )
C
         IF( IERR .NE. 0 ) RETURN
c         print*,'>>>Checking face ',IPT,edglst
C
         IF( IPT .EQ. F1  ) GOTO 500
         IF( IPT .EQ. F2  ) GOTO 500
         IF( IPT .EQ. F3  ) GOTO 500
         IF( IPT .EQ. F4  ) GOTO 500
         IF( IPT .EQ. F5  ) GOTO 500
         IF( IPT .EQ. F6  ) GOTO 500
C
         IE1 = EDGLST(1)
         IE2 = EDGLST(2)
         IE3 = EDGLST(3)
C
         VAL1 = 1E+30
         VAL2 = 1E+30
         VAL3 = 1E+30
C
         IF( .NOT. GETFLG( BIGLST, NODLST, IE1, 4 ) ) THEN
            CALL SETFLG( BIGLST, NODLST, IE1, 4, .TRUE. )
            IF( BIGLST(1,IE1) .EQ. N1  ) GOTO 300
            IF( BIGLST(2,IE1) .EQ. N1  ) GOTO 302
            IF( BIGLST(1,IE1) .EQ. N2  ) GOTO 304
            IF( BIGLST(2,IE1) .EQ. N2  ) GOTO 306
            IF( BIGLST(1,IE1) .EQ. N12 ) GOTO 310
            IF( BIGLST(2,IE1) .EQ. N12 ) GOTO 310
c            print*,'   Checking edge ',IE1
            VAL1 = EDGDST( BIGLST, NODLST, IE1, X, Y, Z, ONX, ONY, ONZ,
     :                     DIRX, DIRY, DIRZ, LNGTH, MNDST, MXRAD )
c            print*,'      Got val1: ',val1
         END IF
C
         GOTO 310
C
 300     IF( BIGLST(2,IE1) .EQ. N2 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
         GOTO 310
C
 302     IF( BIGLST(1,IE1) .EQ. N2 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
         GOTO 310
C
 304     IF( BIGLST(2,IE1) .EQ. N1 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
         GOTO 310
C
 306     IF( BIGLST(1,IE1) .EQ. N1 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
C
 310     IF( .NOT. GETFLG( BIGLST, NODLST, IE2, 4 ) ) THEN
            CALL SETFLG( BIGLST, NODLST, IE2, 4, .TRUE. )
            IF( BIGLST(1,IE2) .EQ. N1  ) GOTO 350
            IF( BIGLST(2,IE2) .EQ. N1  ) GOTO 352
            IF( BIGLST(1,IE2) .EQ. N2  ) GOTO 354
            IF( BIGLST(2,IE2) .EQ. N2  ) GOTO 356
            IF( BIGLST(1,IE2) .EQ. N12 ) GOTO 360
            IF( BIGLST(2,IE2) .EQ. N12 ) GOTO 360
c            print*,'   Checking edge ',IE2
            VAL2 = EDGDST( BIGLST, NODLST, IE2, X, Y, Z, ONX, ONY, ONZ,
     :                     DIRX, DIRY, DIRZ, LNGTH, MNDST, MXRAD )
c            print*,'      Got val2: ',val2
         END IF
C
         GOTO 360
C
 350     IF( BIGLST(2,IE1) .EQ. N2 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
         GOTO 360
C
 352     IF( BIGLST(1,IE1) .EQ. N2 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
         GOTO 360
C
 354     IF( BIGLST(2,IE1) .EQ. N1 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
         GOTO 360
C
 356     IF( BIGLST(1,IE1) .EQ. N1 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
C
 360     IF( .NOT. GETFLG( BIGLST, NODLST, IE3, 4 ) ) THEN
            CALL SETFLG( BIGLST, NODLST, IE3, 4, .TRUE. )
            IF( BIGLST(1,IE3) .EQ. N1  ) GOTO 400
            IF( BIGLST(2,IE3) .EQ. N1  ) GOTO 402
            IF( BIGLST(1,IE3) .EQ. N2  ) GOTO 404
            IF( BIGLST(2,IE3) .EQ. N2  ) GOTO 406
            IF( BIGLST(1,IE3) .EQ. N12 ) GOTO 410
            IF( BIGLST(2,IE3) .EQ. N12 ) GOTO 410
c            print*,'   Checking edge ',IE3
            VAL3 = EDGDST( BIGLST, NODLST, IE3, X, Y, Z, ONX, ONY, ONZ,
     :                     DIRX, DIRY, DIRZ, LNGTH, MNDST, MXRAD )
c            print*,'      Got val3: ',val3
         END IF
C
         GOTO 410
C
 400     IF( BIGLST(2,IE1) .EQ. N2 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
         GOTO 410
C
 402     IF( BIGLST(1,IE1) .EQ. N2 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
         GOTO 410
C
 404     IF( BIGLST(2,IE1) .EQ. N1 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
         GOTO 410
C
 406     IF( BIGLST(1,IE1) .EQ. N1 ) THEN
            MXRDS = 0.0
            GOTO 600
         END IF
C
 410     MXRDS = MIN(MXRDS,VAL1)
         MXRDS = MIN(MXRDS,VAL2)
         MXRDS = MIN(MXRDS,VAL3)
C
         IF( MXRDS .EQ. 0.0 ) GOTO 600
C
 500     IPT = NXTFAC
C
      END DO
C
 600  CONTINUE
c      print*,'Leaving CHKFAC... '
c      print*,'  MXRAD:',MXRAD
c      print*,'  MXRND:',MXRND
c      print*,'  MXRDS:',MXRDS
c
      RETURN
      END
C
