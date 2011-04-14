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
      SUBROUTINE FNDELM( ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD,   NELM,
     :                   NODX,   NODY,   NODZ,   ORGMTX,
     :                   XN, YN, ZN,     ELMGES, SEARCH, NEWMTX )
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER ELMGES
C
      REAL XN, YN, ZN, NEWMTX(9)
C
      LOGICAL SEARCH
C
      INCLUDE 'blkorg.i'
C
      INCLUDE 'blkerr.i'
C
      REAL XSHIFT, YSHIFT, ZSHIFT, XSCALE, YSCALE, ZSCALE
      COMMON / BLKSCL / XSHIFT, YSHIFT, ZSHIFT, XSCALE, YSCALE, ZSCALE
C
      INTEGER I, I1, I2, I3, I4, j
C
      REAL BASIS(5), SUM, v1, v2, v3, v4, det, getdtm
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** FNDELM: GOT ERROR ON ENTRY'
         RETURN
      END IF
C
      IF( ELMGES .LT. 0 ) THEN
         PRINT*,'+++ Warning: FNDELM GIVEN NEGATIVE ELEMENT:',ELMGES
         ELMGES = -ELMGES
      ELSE IF( ELMGES .EQ. 0 ) THEN
         PRINT*,'*** FNDELM: GIVEN ZERO ELEMENT GUESS'
         IERR = -2000
         RETURN
      END IF
C
C - scaling now done within INTRPL (so it scales the original node
C - positions to the cubified domain, rather than scaling the node
C - positions in the cubified domain back to the original domain)
C
c      XN = XN*XSCALE + XSHIFT
c      YN = YN*YSCALE + YSHIFT
c      ZN = ZN*ZSCALE + ZSHIFT
C
      CALL INTRPL( XN, YN, ZN, ELMGES, NODX, NODY, NODZ,
     :             ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :             SZENLS, SZNELS, NNOD, NELM, SEARCH, BASIS )
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** FNDELM: GOT ERROR FROM INTRPL'
         RETURN
      END IF
C
c      XN = (XN - XSHIFT)/XSCALE
c      YN = (YN - YSHIFT)/YSCALE
c      ZN = (ZN - ZSHIFT)/ZSCALE
C
      IF( ELMGES .LE. 0 ) THEN
         PRINT*,'*** FNDELM: INVALID ELEMENT RETURNED BY INTRPL'
         IERR = -2010
         RETURN
      END IF
C
      I1 = ENLIST(ENLBAS(ELMGES)+1)
      I2 = ENLIST(ENLBAS(ELMGES)+2)
      I3 = ENLIST(ENLBAS(ELMGES)+3)
      I4 = ENLIST(ENLBAS(ELMGES)+4)
C
      v1 = basis(1)
      v2 = basis(2)
      v3 = basis(3)
      v4 = basis(4)
      IF( BASIS(1) .LT. 0.0 ) BASIS(1) = 0.0
      IF( BASIS(2) .LT. 0.0 ) BASIS(2) = 0.0
      IF( BASIS(3) .LT. 0.0 ) BASIS(3) = 0.0
      IF( BASIS(4) .LT. 0.0 ) BASIS(4) = 0.0
      SUM = BASIS(1) + BASIS(2) + BASIS(3) + BASIS(4)
c
      SUM = 1.0/SUM
      BASIS(1) = BASIS(1)*SUM
      BASIS(2) = BASIS(2)*SUM
      BASIS(3) = BASIS(3)*SUM
      BASIS(4) = BASIS(4)*SUM
c
      if( abs(sum-1.0) .gt. 4e-4 ) then
         print*,'+++ FNDELM: Warning, sum not one...'
         print*,'   point: ',xn,yn,zn
         print*,'   elmnt: ',elmges,1.0/sum
         print*,'   bases: ',v1,v2,v3,v4
         print*,'   using: ',(basis(j),j=1,4)
         elmges = -elmges
      end if
C
      DO I = 1, 9
         NEWMTX(I) = ORGMTX(I,I1)*BASIS(1) + ORGMTX(I,I2)*BASIS(2)
     :             + ORGMTX(I,I3)*BASIS(3) + ORGMTX(I,I4)*BASIS(4)
      END DO
      NEWMTX(1) = NEWMTX(1)*XSCALE*XSCALE
      NEWMTX(2) = NEWMTX(2)*XSCALE*YSCALE
      NEWMTX(3) = NEWMTX(3)*XSCALE*ZSCALE
      NEWMTX(4) = NEWMTX(4)*YSCALE*XSCALE
      NEWMTX(5) = NEWMTX(5)*YSCALE*YSCALE
      NEWMTX(6) = NEWMTX(6)*YSCALE*ZSCALE
      NEWMTX(7) = NEWMTX(7)*ZSCALE*XSCALE
      NEWMTX(8) = NEWMTX(8)*ZSCALE*YSCALE
      NEWMTX(9) = NEWMTX(9)*ZSCALE*ZSCALE
c
      det = getdtm( newmtx )
c
      if( det .lt. 1e-20 ) then
         print*,'+++ FNDELM: Got very small determinant: ',det
         print*,newmtx(1),newmtx(4),newmtx(7)
         print*,newmtx(2),newmtx(5),newmtx(8)
         print*,newmtx(3),newmtx(6),newmtx(9)
      end if
c
C      if( debug ) then
C         print*,'     Metric: '
C         print*,'     ',newmtx(1),newmtx(4),newmtx(7)
C         print*,'     ',newmtx(2),newmtx(5),newmtx(8)
C         print*,'     ',newmtx(3),newmtx(6),newmtx(9)
C      end if
C
      RETURN
      END
C
