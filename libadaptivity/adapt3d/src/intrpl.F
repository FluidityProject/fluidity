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
      SUBROUTINE INTRPL( XN, YN, ZN, ELMGES, X, Y, Z,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD, NELM, SEARCH, BASIS )
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      use write_log
      IMPLICIT NONE
C
      INTEGER SZENLS, INOD, ELMGES, NNOD, NELM, SZNELS, ENLBAS(NELM+1),
     :        ENLIST(SZENLS), NELBAS(NNOD+1), NELIST(SZNELS),
     :        EELIST(4*NELM)
C
      REAL XN, YN, ZN, X(NNOD), Y(NNOD), Z(NNOD), BASIS(5)
C
      LOGICAL SEARCH
C
      INCLUDE 'blkerr.i'
C
      INCLUDE 'blktim.i'
C
      CHARACTER*100 COORDS2STR
C
      REAL XSHIFT, YSHIFT, ZSHIFT, XSCALE, YSCALE, ZSCALE
      COMMON / BLKSCL / XSHIFT, YSHIFT, ZSHIFT, XSCALE, YSCALE, ZSCALE
C
      INTEGER, PARAMETER :: MXSRCH = 200000, MXPRV = 20000
C      INTEGER SRCHLS(MXSRCH), PRVELM(MXPRV)
      integer, dimension(:), allocatable :: srchls, prvelm
C
      INTEGER CURELM, ILN, IBAS, IELM, BSTELM, I1, I2, I3, I4, I,
     :        ITRACK, NXTELM, ISRCH, IPOS, MINTRK, CNTR
C
      REAL XX(10), YY(10), ZZ(10), XB(10), YB(10), ZB(10), MINBSS,
     :     BSTBSS(4), SUM, TRKBSS, TRYBSS
C
      REAL, PARAMETER :: TOLBSS = -2E-4
C
      LOGICAL SEQNTL
C
      allocate(srchls(mxsrch), prvelm(mxprv))
C
      XX(5) = XN
      YY(5) = YN
      ZZ(5) = ZN
C
      CURELM = ELMGES
      BSTELM = ELMGES
      MINBSS = -1E+30
      BSTBSS(1:4) = -1E+30
C
      ITRACK = 1
      ISRCH  = 0
      CNTR   = 0
C
      SEQNTL = .FALSE.
C
      DO
C
        IBAS = ENLBAS(CURELM)
C
        IF( SEQNTL ) THEN
          IBAS = ABS(IBAS)
        ELSE IF( IBAS .GE. 0 ) THEN
          IF( ISRCH .LT. MXSRCH ) THEN
            ISRCH = ISRCH + 1
            SRCHLS(ISRCH)  = CURELM
            ENLBAS(CURELM) = -IBAS
            IF( MOD(ISRCH,5000) .EQ. 0 ) THEN
              IF( MINBSS .GT. TOLBSS*2 ) EXIT  ! be a little more tolerant
              ewrite(1,*)
     :            '--- INTRPL: Doing lots of searching... '
c              ewrite(1,*) '  Operation: ',trim(opernam)
              ewrite(1,*)
     :            '  Best so far: ',BSTELM,' (out of ',ISRCH,')'
              ewrite(1,*)
     :            '  Bases: ',BSTBSS(1:4)
              ewrite(1,*)
     :            '  Point: ',XN,YN,ZN
            END IF
          ELSE
            IF( MINBSS .GT. TOLBSS*3 ) EXIT  ! be more tolerant
            SEQNTL = .TRUE.
            CURELM = 1
            ewrite(1,*)
     :            '+++ INTRPL: Reached end of search list '//
     :            ' - starting sequential search'
c            ewrite(1,*) '  Operation: ',trim(opernam)
            ewrite(1,*)
     :            '  Best so far: ',BSTELM,' (out of ',ISRCH,')'
            ewrite(1,*)
     :            '  Bases: ',BSTBSS(1:4)
            ewrite(1,*)
     :            '  Point: ',XN,YN,ZN
            CYCLE
          END IF
        END IF
C
C - get the node positions for the element, scaled to the cubified domain
C
        I1 = ENLIST(ABS(IBAS)+1)
        I2 = ENLIST(ABS(IBAS)+2)
        I3 = ENLIST(ABS(IBAS)+3)
        I4 = ENLIST(ABS(IBAS)+4)
C
        XX(1) = (X(I1)-XSHIFT)/XSCALE
        XX(2) = (X(I2)-XSHIFT)/XSCALE
        XX(3) = (X(I3)-XSHIFT)/XSCALE
        XX(4) = (X(I4)-XSHIFT)/XSCALE
C
        YY(1) = (Y(I1)-YSHIFT)/YSCALE
        YY(2) = (Y(I2)-YSHIFT)/YSCALE
        YY(3) = (Y(I3)-YSHIFT)/YSCALE
        YY(4) = (Y(I4)-YSHIFT)/YSCALE
C
        ZZ(1) = (Z(I1)-ZSHIFT)/ZSCALE
        ZZ(2) = (Z(I2)-ZSHIFT)/ZSCALE
        ZZ(3) = (Z(I3)-ZSHIFT)/ZSCALE
        ZZ(4) = (Z(I4)-ZSHIFT)/ZSCALE
C
C - find element's basis functions for searched point
C
        CALL VALS3D( XX, YY, ZZ, 5, XX(5), YY(5), ZZ(5), BASIS )
C
        IPOS   = CURELM*4 + 1
        TRYBSS = BASIS(1)
        SUM    = 0.0
        TRKBSS = 1E+30
        MINTRK = 0
C
        DO I = 1, 4
          SUM = SUM + BASIS(I)
          TRYBSS = MIN(TRYBSS,BASIS(I))
          IF( BASIS(I) .LT. TRKBSS ) THEN
            NXTELM = EELIST(IPOS-I)
            IF( NXTELM .GT. 0 ) THEN
              IF( ENLBAS(NXTELM) .GE. 0 ) THEN
                TRKBSS = BASIS(I)
                MINTRK = I
              END IF
            END IF
          END IF
        END DO
C
        SUM = SUM - 1.0   ! error in basis sum (maybe some use in tol test?)
C
        IF( TRYBSS .GT. MINBSS ) THEN  ! best so far
          BSTBSS(1:4) = BASIS(1:4)
          MINBSS = TRYBSS
          BSTELM = CURELM
          CNTR   = 0  ! reset the counter now we have a new best
          IF( TRYBSS .GE. 0.0 ) EXIT  ! point really is inside element
        ELSE IF( MINBSS .GT. TOLBSS ) THEN  ! best so far had pt nearly inside
          CNTR = CNTR + 1  ! count elems since pt was found nearly inside
          IF( CNTR .GT. 15 ) EXIT  ! give up and take what we got
        END IF
C
        IF( .NOT. SEARCH ) THEN  ! only check elem. provided in ELMGES
          IF( TRYBSS .LT. TOLBSS*5 ) THEN
            ewrite(1,*)
     :            '+++ INTRPL: WARNING: node outside element!'
c            ewrite(1,*) '  Operation: ',trim(opernam)
            ewrite(1,*)
     :            '  Elem: ',BSTELM
            ewrite(1,*)
     :            '  Bases: ',BSTBSS(1:4)
            ewrite(1,*)
     :            '  Point: ',XN,YN,ZN
          END IF
          EXIT
        END IF
C
        IF( SEQNTL ) THEN
          IF( TRYBSS .GT. TOLBSS*5 ) EXIT  ! be even more tolerant
          CURELM = CURELM + 1
          IF( CURELM .GT. NELM ) THEN
            ewrite(1,*)
     :            '+++ INTRPL: Reached end of sequential search!'
c            ewrite(1,*) '  Operation: ',trim(opernam)
            ewrite(1,*)
     :            '  Best: ',BSTELM
            ewrite(1,*)
     :            '  Bases: ',BSTBSS(1:4)
            ewrite(1,*)
     :            '  Point: ',XN,YN,ZN
            EXIT
          END IF
        ELSE
          IF( MINTRK .LT. 1 ) THEN
            ITRACK = ITRACK - 1
            IF( ITRACK .GT. 0 ) THEN
              CURELM = PRVELM(ITRACK)
            ELSE
              IF( MINBSS .GT. TOLBSS*3 ) EXIT  ! be more tolerant
              CURELM = 1
              SEQNTL = .TRUE.
              ewrite(1,*)
     :            '+++ INTRPL: Back to start of previous list '//
     :            ' - starting sequential search'
c              ewrite(1,*) '  Operation: ',trim(opernam)
              ewrite(1,*)
     :            '  Best so far: ',BSTELM,' (out of ',ISRCH,')'
              ewrite(1,*)
     :            '  Bases: ',BSTBSS(1:4)
              ewrite(1,*)
     :            '  Point: ',XN,YN,ZN
            END IF
          ELSE
            PRVELM(ITRACK) = CURELM
            ITRACK = ITRACK + 1
            IF( ITRACK .GT. MXPRV ) THEN
              IF( MINBSS .GT. TOLBSS*3 ) EXIT  ! be more tolerant
              CURELM = 1
              SEQNTL = .TRUE.
              ewrite(1,*)
     :            '+++ INTRPL: Reached end of previous list '//
     :            ' - starting sequential search'
c              ewrite(1,*) '  Operation: ',trim(opernam)
              ewrite(1,*)
     :            '  Best so far: ',BSTELM,' (out of ',ISRCH,')'
              ewrite(1,*)
     :            '  Bases: ',BSTBSS(1:4)
              ewrite(1,*)
     :            '  Point: ',XN,YN,ZN
            ELSE
              CURELM = EELIST(IPOS-MINTRK)
            END IF
          END IF
        END IF
C
      END DO
C
      ELMGES = BSTELM
      BASIS(1:4) = BSTBSS(1:4)
c      opernam = trim(opernam)//' INTRPL'
C
      IF( MINBSS .LT. 0.0 ) THEN
C
        SUM = 0.0
        DO I = 1, 4
          IF( BASIS(I) .LT. 0.0 ) BASIS(I) = 0.0
          SUM = SUM + BASIS(I)
        END DO
C
        SUM = 1.0/SUM
        BASIS(1:4) = BASIS(1:4)*SUM
C
        IF( BSTELM .NE. CURELM ) THEN
C
          IBAS = ABS(ENLBAS(ELMGES))
C
          I1 = ENLIST(ABS(IBAS)+1)
          I2 = ENLIST(ABS(IBAS)+2)
          I3 = ENLIST(ABS(IBAS)+3)
          I4 = ENLIST(ABS(IBAS)+4)
C
          XX(1) = (X(I1)-XSHIFT)/XSCALE
          XX(2) = (X(I2)-XSHIFT)/XSCALE
          XX(3) = (X(I3)-XSHIFT)/XSCALE
          XX(4) = (X(I4)-XSHIFT)/XSCALE
C
          YY(1) = (Y(I1)-YSHIFT)/YSCALE
          YY(2) = (Y(I2)-YSHIFT)/YSCALE
          YY(3) = (Y(I3)-YSHIFT)/YSCALE
          YY(4) = (Y(I4)-YSHIFT)/YSCALE
C
          ZZ(1) = (Z(I1)-ZSHIFT)/ZSCALE
          ZZ(2) = (Z(I2)-ZSHIFT)/ZSCALE
          ZZ(3) = (Z(I3)-ZSHIFT)/ZSCALE
          ZZ(4) = (Z(I4)-ZSHIFT)/ZSCALE
C
        END IF
C
        XN = XX(1)*BASIS(1) + XX(2)*BASIS(2)
     :     + XX(3)*BASIS(3) + XX(4)*BASIS(4)
        YN = YY(1)*BASIS(1) + YY(2)*BASIS(2)
     :     + YY(3)*BASIS(3) + YY(4)*BASIS(4)
        ZN = ZZ(1)*BASIS(1) + ZZ(2)*BASIS(2)
     :     + ZZ(3)*BASIS(3) + ZZ(4)*BASIS(4)
C
        XX(5) = XX(5) - XN
        YY(5) = YY(5) - YN
        ZZ(5) = ZZ(5) - ZN
C
        IF( ABS(XX(5)) .GT. 3E-3 .OR. ABS(YY(5)) .GT. 3E-3 .OR.
     :      ABS(ZZ(5)) .GT. 3E-3 ) THEN
c           ewrite(1,100)  -bug in gfortran for multiple Fx.y formatting?
           ewrite(1,'(A)')
     :              '+++ INTRPL: Shifted node by '//
     :              TRIM(COORDS2STR(XX(5),YY(5),ZZ(5),'(1PE12.4)'))//
     :              ' (to '//TRIM(COORDS2STR(XN,YN,ZN,'(F11.5)'))//')'
        END IF
C
      END IF
C
C - fix up all the flagging before we return
C
      DO I = 1, ISRCH
        CURELM = SRCHLS(I)
        ENLBAS(CURELM) = ABS(ENLBAS(CURELM))
      END DO
C
c 100  format (a,3(1pe12.4),a,3f11.5,a)
c
      deallocate(srchls, prvelm)
C
      RETURN
      END SUBROUTINE INTRPL
C
C

C
C - get around possible gfortran formatting bug by writing coords to string
C
      CHARACTER*100 FUNCTION COORDS2STR( X,Y,Z, FMT )
      IMPLICIT NONE
      REAL X, Y, Z
      CHARACTER*(*) FMT
      CHARACTER*20 STR1, STR2, STR3
      IF( FMT(1:1) .NE. '(' ) THEN
        WRITE( STR1, * ) X
        WRITE( STR2, * ) Y
        WRITE( STR3, * ) Z
      ELSE
        WRITE( STR1, FMT ) X
        WRITE( STR2, FMT ) Y
        WRITE( STR3, FMT ) Z
      END IF
      COORDS2STR = TRIM(ADJUSTL(STR1))//' '//TRIM(ADJUSTL(STR2))//' '//
     :             ADJUSTL(STR3)
      RETURN
      END FUNCTION COORDS2STR
C
