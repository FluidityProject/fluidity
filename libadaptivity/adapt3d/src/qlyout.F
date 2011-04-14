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
      SUBROUTINE QLYOUT( BIGLST, NODLST, MSHNAM, LENFIL, COUNT )
C-----------------------------------------------------------------------
C
C
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER LENFIL, COUNT
C
      CHARACTER*120 MSHNAM
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INTEGER I, J, K, IBAS, QLYFL1, QLYFL2, QLYFL3
C
      LOGICAL DMPQLY, RSPNSE
C
      CHARACTER CNTCHR*3, NUMLST*10
C
      QLYFL1 = 19
      QLYFL2 = 20
      QLYFL3 = 21
C
      NUMLST = '0123456789'
C
      IBAS = COUNT
      I = MOD(IBAS,10)
      IBAS = ( IBAS - I )/10
      J = MOD(IBAS,10)
      IBAS = ( IBAS - J )/10
      K = MOD(IBAS,10)
      IBAS = ( IBAS - K )/10
      I = I + 1
      J = J + 1
      K = K + 1
      CNTCHR = NUMLST(K:K)//NUMLST(J:J)//NUMLST(I:I)
C
      DMPQLY = RSPNSE( 'Do you want to dump the histograms?' ,35 )
C
      IF( DMPQLY ) THEN
C
         MSHNAM = MSHNAM(1:LENFIL)//'.'//CNTCHR//'.hsted'
         OPEN( UNIT = QLYFL1, FILE = MSHNAM(1:LENFIL+10),
     :         STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
C
         MSHNAM = MSHNAM(1:LENFIL)//'.'//CNTCHR//'.hstel'
         OPEN( UNIT = QLYFL2, FILE = MSHNAM(1:LENFIL+10),
     :         STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
C
         MSHNAM = MSHNAM(1:LENFIL)//'.'//CNTCHR//'.hstfn'
         OPEN( UNIT = QLYFL3, FILE = MSHNAM(1:LENFIL+10),
     :         STATUS = 'UNKNOWN', FORM = 'FORMATTED' )
C
         CALL SHWHST( BIGLST, NODLST, QLYFL1, QLYFL2, QLYFL3 )
C
         CLOSE( QLYFL1 )
         CLOSE( QLYFL2 )
         CLOSE( QLYFL3 )
C
         PRINT*,'EDGE SIZES WRITTEN TO FILE ',
     :          MSHNAM(1:LENFIL+8),'ed'
C
         PRINT*,'IN-SPHERE RADII WRITTEN TO FILE ',
     :          MSHNAM(1:LENFIL+8),'el'
C
         PRINT*,'ELEMENT FUNCTIONALS WRITTEN TO FILE ',
     :          MSHNAM(1:LENFIL+8),'fn'
C
      END IF
C
      MSHNAM = MSHNAM(1:LENFIL)
C
      RETURN
      END
C
