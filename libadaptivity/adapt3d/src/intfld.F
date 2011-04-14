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
      SUBROUTINE INTFLD( BIGLST, NODLST,
     :                   X, Y, Z, OLDFLD, M, NM, TOTM,
     :                   ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :                   SZENLS, SZNELS, NNOD, NELM, NHALO, NEWFLD )
C-----------------------------------------------------------------------
C
C - This interpolates, from the old mesh into NEWFLD, the NM fields in
C - OLDFLD with M(I) components in the Ith field.
C - TOTM is the sum of M(1) to M(NM). This allows a long vector with
C - moments, nodes and groups to be interpolated in one call by sending
C - TOTM as moments*groups, NM as no. of groups and each M(I) as the
C - number of moments.
C - OLDFLD has length TOTM*NNOD and NEWFLD should have length TOTM*NEWNND.
C - This can be called as many times as necessary for appropriate fields.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER NNOD, NM, M(NM), TOTM, SZENLS, SZNELS, NELM, NHALO,
     :        ENLBAS(NELM+1), ENLIST(SZENLS), NELBAS(NNOD+1),
     :        NELIST(SZNELS), EELIST(NELM*4)
C
      REAL X(NNOD), Y(NNOD), Z(NNOD), OLDFLD(TOTM*NNOD)
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      REAL NEWFLD(NUMNDS*TOTM)
C
      INCLUDE 'blkerr.i'
C
      INTEGER I1, I2, I3, I4, IPTR, IM1, IM2, IM3, IM4, ELMGES, INOD,
     :        NXT, II, IJ, IM, I, J, JNOD, KNOD, NODCNT
C
      REAL BASIS(5), XN, YN, ZN
C
      LOGICAL GTNDFL
C
      IF( IERR .NE. 0 ) RETURN
C
      if( debug ) then
         print*,'Into INTFLD...'
         print*,'NM,TOTM: ',NM,TOTM
C         print*,'M: ',M
c      pause
      end if
c
      INOD = -1
      JNOD = NUMNDS - NHALO - 1
C
      IPTR = STTNOD
      NODCNT = 0
C
  10  IF( GTNDFL(BIGLST,NODLST,IPTR,1) ) THEN
C
        NODCNT = NODCNT + 1
        IF( NODCNT .GT. NUMNDS ) THEN
           WRITE(0,*) '*** INTFLD: Too many nodes!'
           WRITE(0,*) nodcnt,numnds,iptr,topnod
           ierr = -9494
           return
        end if
C
        IF( GTNDFL(BIGLST,NODLST,IPTR,8) ) THEN
           JNOD = JNOD + 1
           KNOD = JNOD
C           print*,'Numbering halo node (',iptr,') to ',knod+1,
C     :             nodlst(6,iptr),nodlst(1,iptr),nodlst(2,iptr),
C     :             nodlst(3,iptr)
        ELSE
           INOD = INOD + 1
           KNOD = INOD
C           print*,'Numbering non-halo node (',iptr,') to ',knod+1,
C     :             nodlst(6,iptr),nodlst(1,iptr),nodlst(2,iptr),
C     :             nodlst(3,iptr)
        END IF
        NODLST(17,IPTR) = FLOAT(KNOD+1)
C
        XN     = NODLST(1,IPTR)
        YN     = NODLST(2,IPTR)
        ZN     = NODLST(3,IPTR)
C
        ELMGES = INT(NODLST(16,IPTR))
        NXT    = INT(NODLST(NXTNOD,IPTR))
C
        IF( ELMGES .LE. 0 ) THEN
           WRITE(0,*) '*** INTFLD: INVALID ELEMENT GUESS'
           WRITE(0,*) '  IPTR,KNOD,ELMGES: ',IPTR,KNOD+1,ELMGES
           IERR = -2010
           RETURN
        END IF
C
        CALL INTRPL( XN, YN, ZN, ELMGES, X, Y, Z,
     :               ENLBAS, ENLIST, NELBAS, NELIST, EELIST,
     :               SZENLS, SZNELS, NNOD, NELM, .FALSE., BASIS )
C
        IF( IERR .NE. 0 ) RETURN
C
        I1 = ENLIST(ENLBAS(ELMGES)+1) - 1
        I2 = ENLIST(ENLBAS(ELMGES)+2) - 1
        I3 = ENLIST(ENLBAS(ELMGES)+3) - 1
        I4 = ENLIST(ENLBAS(ELMGES)+4) - 1
C
        II = 0
        IJ = 0
        DO I = 1, NM
           IM  = M(I)*KNOD
           IM1 = M(I)*I1
           IM2 = M(I)*I2
           IM3 = M(I)*I3
           IM4 = M(I)*I4
           DO J = 1, M(I)
              NEWFLD(J+IM+IJ) = OLDFLD(J+IM1+II)*BASIS(1)
     :                        + OLDFLD(J+IM2+II)*BASIS(2)
     :                        + OLDFLD(J+IM3+II)*BASIS(3)
     :                        + OLDFLD(J+IM4+II)*BASIS(4)
           END DO
           IJ = IJ + M(I)*NUMNDS
           II = II + M(I)*NNOD
        END DO
C
        IPTR = NXT
        IF( IPTR .LE. TOPNOD ) GOTO 10
C
        IF( NODCNT .NE. NUMNDS ) THEN
           WRITE(0,*) '*** INTFLD: Reached end of node list too early'
           WRITE(0,*) nodcnt,numnds,iptr,topnod
           ierr = -9393
           return
        end if
C
      END IF
C
      RETURN
      END
C
