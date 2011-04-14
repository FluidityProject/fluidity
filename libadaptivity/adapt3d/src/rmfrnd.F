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
      SUBROUTINE RMFRND( BIGLST, NODLST, IPTR )
C-----------------------------------------------------------------------
C
C - This subroutine takes a node in nodlst from the 'front' to the back.
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
      INTEGER IPTR, INEWND, LST, NXT
C
      REAL X, Y, Z
C
c      print*,' '
c      print*,'Into RMFRND...',ierr
c      print*,'     IPTR:',iptr
C
      IF( IERR .NE. 0 ) RETURN
C
      IF( IPTR .LE. 0 ) THEN
         IERR = -1
         PRINT*,'*** RMFRND: INVALID NODE POINTER'
         RETURN
      END IF
C
      IF( IPTR .NE. ENDNOD ) THEN
C
         LST = INT(NODLST(LSTNOD,IPTR))
         NXT = INT(NODLST(NXTNOD,IPTR))
C
         IF( IPTR .EQ. STFRND ) STFRND = NXT
         IF( IPTR .EQ. STTNOD ) STTNOD = NXT
C
         IF( LST .GT. 0 ) NODLST(NXTNOD,LST) = FLOAT(NXT)
         IF( NXT .GT. 0 ) NODLST(LSTNOD,NXT) = FLOAT(LST)
C
         NODLST(LSTNOD,IPTR) = FLOAT(ENDNOD)
         NODLST(NXTNOD,IPTR) = FLOAT(EMTNOD)
C
         IF( ENDNOD .GT. 0 ) NODLST(NXTNOD,ENDNOD) = FLOAT(IPTR)
         IF( EMTNOD .GT. 0 ) NODLST(LSTNOD,EMTNOD) = FLOAT(IPTR)
C
         ENDNOD = IPTR
C
      END IF
C
      NFRTND = NFRTND - 1
C
c      print*,'Leaving RMFRND'
c
      RETURN
      END
C
