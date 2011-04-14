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
      SUBROUTINE ADFRND( BIGLST, NODLST, X, Y, Z, IPTR, INEWND )
C-----------------------------------------------------------------------
C
C - This subroutine adds a node into nodlst.
C
C-----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER IPTR, INEWND, LST, NXT
C
      REAL X, Y, Z
C
      INCLUDE 'blknew.i'
C
      INCLUDE 'blkbig.i'
C
      INCLUDE 'blkerr.i'
C
      INTEGER I
c
c      print*,' '
c      print*,'Into ADFRND...',ierr
c      print*,'        X:',X
c      print*,'        Y:',Y
c      print*,'        Z:',Z
C
      IF( IERR .NE. 0 ) RETURN
C
      IPTR = EMTNOD
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** ADFRND: FAILED TO ADD NODE TO NODLST - FULL'
         RETURN
      END IF
C
      IF( NODLST(6,IPTR) .NE. 0.0 ) THEN
         IERR = -2
         PRINT*,'*** ADFRND: NODE BLOCK NOT EMPTY'
         RETURN
      END IF
C
      LST = INT(NODLST(LSTNOD,IPTR))
      NXT = INT(NODLST(NXTNOD,IPTR))
C
      IF( ENDNOD .NE. LST ) THEN
         IERR = -110
         PRINT*,'*** ADFRND: INCONSISTENCY AT END OF NODLST'
         RETURN
      END IF
C
      NODLST(1,IPTR) = X
      NODLST(2,IPTR) = Y
      NODLST(3,IPTR) = Z
C
      NODLST(6,IPTR) = 1.0
C
      CALL FIXEIG( NODLST(1,IPTR), NODLST(7,IPTR), ANISOT )
C      DO I = 1, 9
C         NODLST(I+6,IPTR) = 0.0
C      END DO
C      NODLST( 7,IPTR) = 44.44444
C      NODLST(11,IPTR) = 44.44444
C      NODLST(15,IPTR) = 400.0
C
      IF( LST .GT. 0 ) NODLST(NXTNOD,LST) = FLOAT(NXT)
      IF( NXT .GT. 0 ) NODLST(LSTNOD,NXT) = FLOAT(LST)
C
      EMTNOD = NXT
C
      NODLST(LSTNOD,IPTR) = -1
      NODLST(NXTNOD,IPTR) = STFRND
      STFRND = IPTR
C
      NFRTND = NFRTND + 1
      NUMNDS = NUMNDS + 1
c
c      print*,'Leaving ADFRND...'
c      print*,'     IPTR:',iptr
c      print*,'   INEWND:',inewnd
C
      RETURN
      END
C
