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
      SUBROUTINE MK1ELE( BIGLST, NODLST,
     :                   IFAC1, IFAC2, IFAC3, IFAC4, IREG )
C-----------------------------------------------------------------------
C
C - This subroutine 
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
      INTEGER IFAC1, IFAC2, IFAC3, IFAC4, IREG
C
      INTEGER CONFC1(3), EDGLS1(3), IELM1, NXT1, CONELM(4),
     :        CONFC2(3), EDGLS2(3), IELM2, NXT2, EDGLST(6),
     :        CONFC3(3), EDGLS3(3), IELM3, NXT3, IFLAG, INEW,
     :        CONFC4(3), EDGLS4(3), IELM4, NXT4, LST, IELM,
     :        EDG12, EDG13, EDG23, EDG14, EDG24, EDG34
C
      PRINT*,'Into MK1ELE...',ierr
      print*,'   face1,face2:',ifac1,ifac2
      print*,'   face3,face4:',ifac3,ifac4
      print*,'        nfaces:',nfaces
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IFAC1 .EQ. IFAC2 ) THEN
         IERR = -60
      ELSE IF( IFAC1 .EQ. IFAC3 ) THEN
         IERR = -60
      ELSE IF( IFAC1 .EQ. IFAC4 ) THEN
         IERR = -60
      ELSE IF( IFAC2 .EQ. IFAC3 ) THEN
         IERR = -60
      ELSE IF( IFAC3 .EQ. IFAC4 ) THEN
         IERR = -60
      END IF
C
      IF( IERR .NE. 0 ) THEN
         PRINT*,'*** MK1ELE: FACES ARE NOT DISTINCT'
         RETURN
      END IF
c
      CALL FACINF( IFAC1, CONFC1, EDGLS1, IELM1, INEW, IFLAG,
     :             NXT1, LST )
c
      CALL FACINF( IFAC2, CONFC2, EDGLS2, IELM2, INEW, IFLAG,
     :             NXT2, LST )
c
      CALL FACINF( IFAC3, CONFC3, EDGLS3, IELM3, INEW, IFLAG,
     :             NXT3, LST )
c
      CALL FACINF( IFAC4, CONFC4, EDGLS4, IELM4, INEW, IFLAG,
     :             NXT4, LST )
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( CONFC1(1) .EQ. IFAC2 ) THEN
         EDG12 = 1
      ELSE IF( CONFC1(2) .EQ. IFAC2 ) THEN
         EDG12 = 2
      ELSE IF( CONFC1(3) .EQ. IFAC2 ) THEN
         EDG12 = 3
      ELSE
         IERR = -61
         PRINT*,'*** MK1ELE: FACE 1 AND FACE 2 ARE NOT CONNECTED'
         RETURN
      END IF
C
      IF( CONFC1(1) .EQ. IFAC3 ) THEN
         EDG13 = 1
      ELSE IF( CONFC1(2) .EQ. IFAC3 ) THEN
         EDG13 = 2
      ELSE IF( CONFC1(3) .EQ. IFAC3 ) THEN
         EDG13 = 3
      ELSE
         IERR = -61
         PRINT*,'*** MK1ELE: FACE 1 AND FACE 3 HAVE NO COMMON EDGE'
         RETURN
      END IF
C
      IF( EDG12 .EQ. EDG13 ) THEN
         IERR = -62
         PRINT*,'*** MK1ELE: INCONSISTENCY BETWEEN EDGES AND FACES'
         PRINT*,'edges of face1:',edgls1
         print*,'faces of face1:',confc1
         RETURN
      END IF
C
      IF( CONFC2(1) .EQ. IFAC3 ) THEN
         EDG23 = 1
      ELSE IF( CONFC2(2) .EQ. IFAC3 ) THEN
         EDG23 = 2
      ELSE IF( CONFC2(3) .EQ. IFAC3 ) THEN
         EDG23 = 3
      ELSE
         IERR = -61
         PRINT*,'*** MK1ELE: FACE 2 AND FACE 3 ARE NOT CONNECTED'
         RETURN
      END IF
C
      EDG14 = MOD(EDG12,3) + 1
      IF( EDG14 .EQ. EDG13 ) EDG14 = MOD(EDG14,3) + 1
C
      IF( CONFC1(EDG14) .NE. IFAC4 ) THEN
         IERR = -62
         PRINT*,'*** MK1ELE: INCONSISTENCY BETWWEN EDGES AND FACES'
         PRINT*,'edges of face1:',edgls1
         print*,'edges of face4:',edgls4
         print*,'faces of face1:',confc1
         print*,'faces of face4:',confc4
         RETURN
      END IF
C
      EDG24 = MOD(EDG23,3) + 1
      IF( EDGLS2(EDG24) .EQ. EDGLS1(EDG12) ) EDG24 = MOD(EDG24,3) + 1
C
      IF( CONFC3(1) .EQ. IFAC4 ) THEN
         EDG34 = 1
      ELSE IF( CONFC3(2) .EQ. IFAC4 ) THEN
         EDG34 = 2
      ELSE IF( CONFC3(3) .EQ. IFAC4 ) THEN
         EDG34 = 3
      ELSE
         IERR = -61
         PRINT*,'*** MK1ELE: FACE 3 AND FACE 4 ARE NOT CONNECTED'
         RETURN
      END IF
C
      EDGLST(1) = EDGLS1(EDG14)
      EDGLST(2) = EDGLS2(EDG24)
      EDGLST(3) = EDGLS3(EDG34)
      EDGLST(6) = EDGLS1(EDG13)
      EDGLST(5) = EDGLS2(EDG23)
      EDGLST(4) = EDGLS1(EDG12)
C
      CONELM(1) = IELM4
      CONELM(2) = IELM1
      CONELM(3) = IELM3
      CONELM(4) = IELM2
C
      CALL RMVFAC( IFAC1 )
      CALL RMVFAC( IFAC2 )
      CALL RMVFAC( IFAC3 )
      CALL RMVFAC( IFAC4 )
C
      CALL ADDELE( BIGLST, NODLST, CONELM, EDGLST, IELM, INEW, IREG )
      if( ierr .eq. 12345 ) then
         print*,'--- got bad oppedg in mk1ele'
         ierr = 0
      endif
C
      RETURN
      END
C
