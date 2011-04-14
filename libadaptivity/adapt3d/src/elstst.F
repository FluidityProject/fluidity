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
      SUBROUTINE ELSTST( BIGLST, NODLST, IELM1, IREG1, CONEL1, EDGLS1,
     :                   FNCDIF, FNCORG, AVEDIF, AVEORG, ELM2,
     :                   DOTOP, MINCHG )
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
      INCLUDE 'blktim.i'
C
      INTEGER IELM1, CONEL1(4), EDGLS1(6), ELM2, IREG1
C
      REAL FNCDIF, FNCORG, AVEDIF, AVEORG, DOTOP, MINCHG
C
      INTEGER IELM2, INEW, IREG2, IFLAG, NXT, LST, ED1, ED2, ED3, ED4,
     :        ED5, ED6, ED7, ED8, ED9, OPED1, OPED2, OPED3, OPED4, I, J,
     :        OPED5, OPED6, N1, N2, N3, NE1, NE2, FAC1, FAC2, FAC3,
     :        FAC4, FAC5, NODS1(4), NODS2(4), CONEL2(4), EDGLS2(6),
     :        INFO(4)
C      IELM2, CONEL2(4), EDGLS2(6), INEW, IFLAG, NXT, LST, IREG2,
C     :        I, J, ED1, ED2, ED3, CONELS(6), N1, N2, N3, NE1, NE2,
C     :        OPED1, OPED2, OPED3, OPED4, OPED5, OPED6, SURFC1, SURFC2,
C     :        ND(4), tmpnxt, tmplst, hasnod
C
      REAL CURDIF, F0, F1, F2, F3, F4, F5, DX, DY, DZ, DOT1, DOT2, DOT3,
     :     ELMFNC, EDGFNC, VOL1, VOL2, UNDLEN,
     :     a,b,c,d,e,xx(4),yy(4),zz(4),tetvol
C
      INTEGER HASEDG
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** ELSTST: GOT ERROR ON ENTRY'
         RETURN
      END IF
C
      CALL ADAPT_SECOND( TIMSTT )
C
      FNCDIF = 1E+30
      FNCORG = 1E+31
      AVEDIF = 1E+30
      AVEORG = 1E+31
C
      IELM2 = CONEL1(ELM2)
C
      CALL ELMINF( BIGLST, NODLST, IELM2, CONEL2, EDGLS2,
     :             INEW, IFLAG, IREG2, NXT, LST )
C
C - can't do swap if elements have different 'materials'
C
      IF( IREG1 .NE. IREG2 ) RETURN
C
C - get the list of nodes for both elements
C
      CALL ELNODS( BIGLST, NODLST, IELM1, NODS1, NXT, LST, .TRUE. )
C
      CALL ELNODS( BIGLST, NODLST, IELM2, NODS2, NXT, LST, .TRUE. )
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** ELSTST: GOT ERROR FROM ELMINF/ELNODS'
         print*,'  elems: ',ielm1,ielm2
         print*,'  conel1: ',conel1
         print*,'  edgls1: ',edgls1
         print*,'  nods1 : ',nods1
         print*,'  conel2: ',conel2
         print*,'  edgls2: ',edgls2
         print*,'  nods2 : ',nods2
         RETURN
      END IF
C
c      debug = .false.
C      if( ielm1 .eq. 1164 .and. ielm2 .eq. 1674 ) debug = .true.
c
c      if( debug ) then
c         print*,'Face for edge test: ',ielm1,ielm2
c         print*,'   conelm for 1164: ',conel1
c         print*,'   edglst for 1164: ',edgls1
c         call elnods( ielm1, nd, nxt, lst )
c         print*,'   nodes  for 1164: ',nd
c         print*,'   conelm for 1674: ',conel2
c         print*,'   edglst for 1674: ',edgls2
c         call elnods( ielm2, nd, nxt, lst )
c         print*,'   nodes  for 1674: ',nd
c      end if
C
c      IF( IREG1 .NE. IREG2 ) RETURN
C
C
      INFO(1) = 5318
      INFO(2) = 6535
      INFO(3) = 15243
      INFO(4) = 10599
C
      FAC1 = ELM2
c      IF( CONEL1(1) .EQ. IELM2 ) THEN
c         FAC1 = 1
c      ELSE IF( CONEL1(2) .EQ. IELM2 ) THEN
c         FAC1 = 2
c      ELSE IF( CONEL1(3) .EQ. IELM2 ) THEN
c         FAC1 = 3
c      ELSE IF( CONEL1(4) .EQ. IELM2 ) THEN
c         FAC1 = 4
c      ELSE
c         PRINT*,'*** ELSTST: FAILED TO FIND CONNECTED ELEMENT (1)'
c         print*,ielm1,ielm2
c         print*,conel1
c         print*,conel2
c         print*,edgls1
c         print*,edgls2
c         IERR = -160
c         RETURN
c      END IF
C
C      FAC2 = 0
      IF( CONEL2(1) .EQ. IELM1 ) THEN
         FAC2 = 1
      ELSE IF( CONEL2(2) .EQ. IELM1 ) THEN
         FAC2 = 2
      ELSE IF( CONEL2(3) .EQ. IELM1 ) THEN
         FAC2 = 3
      ELSE IF( CONEL2(4) .EQ. IELM1 ) THEN
         FAC2 = 4
      ELSE
         PRINT*,'*** ELSTST: FAILED TO FIND CONNECTED ELEMENT'
         print*,ielm1,ielm2
         print*,conel1
         print*,conel2
         print*,edgls1
         print*,edgls2
         IERR = -160
         RETURN
      END IF
C
      NE1 = NODS1(5 - FAC1)
      NE2 = NODS2(5 - FAC2)
      I = INFO(FAC1)
C      CONELS(1) = CONEL1(MOD(I,4) + 1)
C      I = I / 4
C      CONELS(2) = CONEL1(MOD(I,4) + 1)
C      I = I / 4
C      CONELS(3) = CONEL1(MOD(I,4) + 1)
C      I = I / 4
      N1 = NODS1(MOD(I,4) + 1)
      I = I / 4
      N2 = NODS1(MOD(I,4) + 1)
      I = I / 4
      N3 = NODS1(MOD(I,4) + 1)
C
      F1 = ELMFNC( BIGLST, NODLST, IELM1, N1, N2, N3, NE1, VOL1 )
      F2 = ELMFNC( BIGLST, NODLST, IELM2, N1, N2, N3, NE2, VOL2 )
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** ELSTST: GOT ERROR FROM ELMFNC'
         RETURN
      END IF
C
      AVEORG = ( F1 + F2 ) / 2
      FNCORG = MAX( F1, F2 )
C
      IF( FNCORG .LE. DOTOP ) RETURN
C
      if( vol1*vol2 .ge. 0.0 ) then
         print*,'*** ELSTST: GOT BAD CONNECTED ELEMENTS'
         PRINT*,'   elements: ',ielm1,ielm2
         print*,'   nodes: ',n1,n2,n3,ne1,ne2
         print*,'   metric in-spheres: ',vol1,vol2
         print*,'   functionals: ',f1,f2
         debug = .true.
         biglst(7,ielm1) = 0
         biglst(7,ielm2) = 0
         f1 = ELMFNC( BIGLST, NODLST, IELM1, N1, N2, N3, NE1, VOL1 )
         f2 = ELMFNC( BIGLST, NODLST, IELM2, N1, N2, N3, NE2, VOL2 )
         print*,'   check 1: ',vol1,vol2
         f1 = ELMFNC( BIGLST, NODLST, 0, N1, N2, N3, NE1, VOL1 )
         f2 = ELMFNC( BIGLST, NODLST, 0, N1, N2, N3, NE2, VOL2 )
         print*,'   check 2: ',vol1,vol2
         print*,'   fnls: ',f1,f2
         xx(1) = nodlst(1,n1)
         yy(1) = nodlst(2,n1)
         zz(1) = nodlst(3,n1)
         xx(2) = nodlst(1,n2)
         yy(2) = nodlst(2,n2)
         zz(2) = nodlst(3,n2)
         xx(3) = nodlst(1,n3)
         yy(3) = nodlst(2,n3)
         zz(3) = nodlst(3,n3)
         xx(4) = nodlst(1,ne1)
         yy(4) = nodlst(2,ne1)
         zz(4) = nodlst(3,ne1)
         print*,'    common node 1: ',xx(1),yy(1),zz(1)
         print*,(nodlst(i,n1),i=7,9)
         print*,(nodlst(i,n1),i=10,12)
         print*,(nodlst(i,n1),i=13,15)
         print*,'    common node 2: ',xx(2),yy(2),zz(2)
         print*,(nodlst(i,n2),i=7,9)
         print*,(nodlst(i,n2),i=10,12)
         print*,(nodlst(i,n2),i=13,15)
         print*,'    common node 3: ',xx(3),yy(3),zz(3)
         print*,(nodlst(i,n3),i=7,9)
         print*,(nodlst(i,n3),i=10,12)
         print*,(nodlst(i,n3),i=13,15)
         print*,'     extra node 1: ',xx(4),yy(4),zz(4)
         print*,(nodlst(i,ne1),i=7,9)
         print*,(nodlst(i,ne1),i=10,12)
         print*,(nodlst(i,ne1),i=13,15)
         vol1 = tetvol( xx, yy, zz )
         xx(4) = nodlst(1,ne2)
         yy(4) = nodlst(2,ne2)
         zz(4) = nodlst(3,ne2)
         print*,'     extra node 2: ',xx(4),yy(4),zz(4)
         print*,(nodlst(i,ne2),i=7,9)
         print*,(nodlst(i,ne2),i=10,12)
         print*,(nodlst(i,ne2),i=13,15)
         vol2 = tetvol( xx, yy, zz )
         print*,'   flat space vols: ',vol1,vol2
         if( vol1*vol2 .ge. 0.0 ) then
          print*,'   Inside-out element detected - signs should differ!'
         else
          print*,'   Vols ok!! -CHECK YOUR METRICS ARE +VE DEF...'
         end if
         print*,'    conel1: ',conel1
         IERR = -402
         return
      end if
C
      IF( ABS(VOL1) .LT. ABS(VOL2) ) VOL2 = -VOL1
C
      A = ELMFNC( BIGLST, NODLST, 0, N2, N3, NE1, NE2, VOL1 )
      if( VOL1/VOL2 .LT. 1E-4 ) GOTO 950
C
      B = ELMFNC( BIGLST, NODLST, 0, N3, N1, NE1, NE2, VOL1 )
      if( VOL1/VOL2 .LT. 1E-4 ) GOTO 950
C
      C = ELMFNC( BIGLST, NODLST, 0, N1, N2, NE1, NE2, VOL1 )
      if( VOL1/VOL2 .LT. 1E-4 ) GOTO 950
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** ELSTST: GOT ERROR FROM ELMFNC (2)'
         RETURN
      END IF
C
      FNCDIF = MAX( MAX( A, B ), C ) - FNCORG
      IF( FNCDIF .LT. 0.0 ) AVEDIF = ( A + B + C ) / 3 - AVEORG
C
c      if( ielm1 .eq. 38431 .and. ielm2 .eq. 34490 ) then
c      if( debug ) then
c         print*,'   fc2ed fnl: ',ielm1,FNCDIF,F1,F3
c         print*,ne1,ne2,n1,n2,n3
c         print*,a,b,c,d,e
c      end if
C
  950 CONTINUE
C
      CALL ADAPT_SECOND( TIMEND )
C
      TIMF2E = TIMF2E + TIMEND - TIMSTT
      TOTF2E = TOTF2E + 1
C
      RETURN
      END
C
